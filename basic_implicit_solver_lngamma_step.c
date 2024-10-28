#include <stdio.h>
#include <string.h>
#include <stdint.h>
#include <math.h>
#include <stdlib.h>
#include <stdbool.h>

#define BUFFER_SIZE 33554432
// CGS unit constants
#define m_e 9.109e-28
#define c 2.998e10
#define sigma_t 6.652e-25

typedef struct LeptonParams
{
    double *current_n;
    double *next_n;
    double *gamma;
    double delta_gamma;
} LeptonParams;

typedef struct SimulationParams
{
    // array parameters
    int64_t array_len;
    double max_gamma;
    double min_gamma;
    int64_t samples_per_decade;
    
    // time values
    double t;
    double end_t;

    // free parameters
    double dt;      // fixed time step
    double R;       // radius of system (spherical geo)
    double inject_power;// power law power
    double init_power;  // power to generate initial distribution
    double inject_min;  // injection gamma range min
    double inject_max;  // injection gamma range max
    double rho;     // background density for initial population calc
    double B;       // background magnetic field
    double L;       // external luminosity
    double end_tol; // tolerance for change in n to count as equilibrium

    double Q_e0;    // starting population
    double S;       // sync vs inv. compton
    double tau_esc; // free escape time
    double norm;    // normalize the prob dist
    double avg_gamma;   // average gamma of injected dist for Q_e0 calc
    double V;           // volume of system (spherical based on R)   

    // code specific values
    int32_t n_species;
    LeptonParams **Species;
    double change;
    bool end_sim;
    double final_time;
    char *buffer;
    size_t *buffer_index;

} SimulationParams;

void malloc_and_fill_gamma_array(SimulationParams *Sim, LeptonParams *Lepton)
{
    int64_t decades;
    // calculate number of decades in the gamma range
    decades = (int64_t)log10(Sim->max_gamma) - log10(Sim->min_gamma);
    // set array length based on decades and input samples per decade
    Sim->array_len = decades * Sim->samples_per_decade;

    // malloc gamma and delta gamma arrays
    Lepton->gamma = calloc(Sim->array_len + 1, sizeof(double));
    
    // fill gamma array with equal log step data
    for (int64_t i = 0; i < Sim->array_len + 1; i++)
    {
        Lepton->gamma[i] = pow(10, log10(Sim->min_gamma) + i / (double)Sim->samples_per_decade);
    }
    // calculate delta gamma for each point with the point ahead of it 
    // (as we are using FDM and stepping backwards)
    Lepton->delta_gamma = log(Lepton->gamma[2]) - log(Lepton->gamma[1]);
}

void malloc_Sim_arrays(SimulationParams *Sim)
{
    Sim->Species = malloc(Sim->n_species * sizeof(LeptonParams *));
    Sim->buffer = malloc(BUFFER_SIZE * sizeof(char));
    Sim->buffer_index = malloc(sizeof(*(Sim->buffer_index)));
    *Sim->buffer_index = 0;

    for (int32_t i = 0; i < Sim->n_species; i++)
    {
        Sim->Species[i] = malloc(sizeof(LeptonParams));
        malloc_and_fill_gamma_array(Sim, Sim->Species[i]);
        Sim->Species[i]->next_n = malloc((Sim->array_len + 1) * sizeof(double));
        Sim->Species[i]->current_n = malloc((Sim->array_len + 1) * sizeof(double));
    }
}

void free_Sim_arrays(SimulationParams *Sim)
{
    for (int32_t i = 0; i < Sim->n_species; i++)
    {
        free(Sim->Species[i]->next_n);
        free(Sim->Species[i]->current_n);
        free(Sim->Species[i]->gamma);
        free(Sim->Species[i]);
    }
    free(Sim->buffer);
    free(Sim->buffer_index);
    free(Sim->Species);
}

double I(double gamma, double min, double max, double power, SimulationParams *Sim)
{
    // inject power law within a specified range
    if (gamma < min || gamma > max)
    {
        return 0.;
    }
    else
    {
        return Sim->Q_e0 * Sim->norm * pow(gamma, (-1. * power));
    }
}

void normalize_power_law_dist(double power, SimulationParams *Sim)
{
        // normalise power law dist based on a given power
        Sim->norm =
        1. / 
        ((pow(Sim->max_gamma, 1. - power) / (1. - power)) 
        - (pow(Sim->min_gamma, 1. - power) / (1. - power)));
}

void set_initial_state(SimulationParams *Sim)
{

    for (int32_t lepton = 0; lepton < Sim->n_species; lepton++)
    {
        Sim->Species[lepton]->next_n[Sim->array_len] = 0.;
        for (int64_t i = 0; i < Sim->array_len; i++)
        {
            normalize_power_law_dist(Sim->init_power, Sim);
            
            // set initial population on a selected power law
            // number of photons included based on background density
            Sim->Species[lepton]->current_n[i] =
            (Sim->rho / (2.*m_e)) * 
            (I(Sim->Species[lepton]->gamma[i], Sim->min_gamma, Sim->max_gamma, Sim->init_power, Sim) 
            / Sim->Q_e0);
            
            normalize_power_law_dist(Sim->inject_power, Sim);
            //Sim->Species[lepton]->current_n[i] = 1.;
        }
    }
}

// file writing code
void flush_buffer(FILE *file, char *buffer, size_t *buffer_index)
{
    if (*buffer_index > 0)
    {
        fwrite(buffer, 1, *buffer_index, file);
        fflush(file);
        *buffer_index = 0;
    }
}

void save_data(FILE *file, SimulationParams *Sim)
{
    int buffer_space = BUFFER_SIZE - (*Sim->buffer_index);

    // Prepare the data to be saved
    int required_space = snprintf(Sim->buffer + (*Sim->buffer_index), buffer_space,
                                  "%lf,",
                                  Sim->t);

    for (int i = 0; i < Sim->array_len; i++)
    {
        required_space += snprintf(Sim->buffer + (*Sim->buffer_index) + required_space, buffer_space - required_space,
                                   "%e,",
                                   Sim->Species[0]->next_n[i]);
    }

    // After the loop, add a newline character
    required_space += snprintf(Sim->buffer + (*Sim->buffer_index) + required_space, buffer_space - required_space,
                               "\n");

    // If the data exceeds the buffer, flush the buffer to the file
    if (required_space >= buffer_space)
    {
        fwrite(Sim->buffer, 1, *Sim->buffer_index, file);
        fflush(file);
        *Sim->buffer_index = 0;
    }

    // Append the data to the buffer
    *Sim->buffer_index += required_space;

    // If the buffer is full, flush it
    if (*Sim->buffer_index >= BUFFER_SIZE)
    {
        fwrite(Sim->buffer, 1, *Sim->buffer_index, file);
        fflush(file);
        *Sim->buffer_index = 0;
    }
}

void calc_S(SimulationParams *Sim)
{
    // calculate synchrotron ratio based on physics
    Sim->S = -(4. / 3.) * c * sigma_t * (Sim->B * Sim->B / (8 * M_PI * m_e * c * c));
}

void calc_Q_e0(SimulationParams *Sim)
{
    // calculate a pre-solved integral of gamma * gamma^-p
    Sim->avg_gamma =
        Sim->norm * (pow(Sim->inject_max, 2. - Sim->inject_power) / (2. - Sim->inject_power)) 
        - Sim->norm * (pow(Sim->inject_min, 2. - Sim->inject_power) / (2. - Sim->inject_power));
    // correct for divide by zero errors when p = 2
    if (Sim->inject_power == 2.)
        Sim->avg_gamma = Sim->norm * log(Sim->inject_max) - Sim->norm * log(Sim->inject_min);
    
    // calculate volume
    Sim->V = (4./3.) * M_PI * Sim->R * Sim->R * Sim->R;
    // calculate injected Qe0
    Sim->Q_e0 = Sim->L / (Sim->V * Sim->avg_gamma * m_e * c * c);
}

void calc_tau_esc(SimulationParams *Sim)
{
    // calculate escape time for a spherical plasma based on the radius
    Sim->tau_esc = (3. / 4.) * (Sim->R / c);
}

void implicit_step(SimulationParams *Sim, LeptonParams *Lepton)
{
    for (int64_t i = Sim->array_len - 1; i >= 0; i--)
    {
        // implicit stepping regime
        Lepton->next_n[i] =
        Sim->tau_esc * 
        (Sim->S * Sim->dt * Lepton->gamma[i] * Lepton->next_n[i+1]
        - Lepton->delta_gamma * Sim->dt * I(Lepton->gamma[i], Sim->inject_min, Sim->inject_max, Sim->inject_power, Sim)
        - Lepton->delta_gamma * Lepton->current_n[i])
        /
        (-2 * Sim->S * Lepton->delta_gamma * Sim->tau_esc * Sim->dt
        +2 * Sim->S * Sim->tau_esc * Sim->dt * Lepton->gamma[i]
        - Lepton->delta_gamma * Sim->tau_esc
        - Lepton->delta_gamma * Sim->dt);
    }
}

void save_step_to_prev_n(SimulationParams *Sim, LeptonParams *Lepton)
{
    for (int64_t i = Sim->array_len - 1; i >= 0; i--)
    {
        // implicit stepping regime
        Lepton->current_n[i] = Lepton->next_n[i];
    }
}

bool equilibrium_check(SimulationParams *Sim, LeptonParams *Lepton)
{
    Sim->change = 0.;
    double dn;
    // calculate percentage change in n
    for (int64_t i = Sim->array_len - 1; i >= 0; i--)
    {
        dn = (Lepton->next_n[i]-Lepton->current_n[i]) / Lepton->current_n[i];
        Sim->change += pow(dn, 2.);
    }
    Sim->change /= (Sim->dt * Sim->dt);
    Sim->change = sqrt(Sim->change);
    
    // check change in population against specified end tolerance
    if (Sim->change < Sim->end_tol)
    {
        printf("equilibrium reached at t = %e, last change %e\n", Sim->t, Sim->change);
        Sim->end_sim=true;
        return true;
    }
    else
    {
        return false;
    }
}

void simulate(FILE *file, SimulationParams *Sim)
{
    // take input params and calculate coefficients
    calc_S(Sim);
    normalize_power_law_dist(Sim->inject_power, Sim);
    calc_Q_e0(Sim);
    calc_tau_esc(Sim);
    
    // ensure dt is not larger than tau_esc for stability
    if (Sim->dt > Sim->tau_esc)
    {
        Sim->dt = Sim->tau_esc;
        printf("dt larger than tau_esc. dt set to tau_esc to ensure stability");
    }

    set_initial_state(Sim);
    Sim->t = 0.;
    Sim->end_sim = false;
    // start simulation
    printf("Start Sim with C %e, tau %e B %.2lf, S %e g_array %lld\n", 
    Sim->Q_e0 * Sim->norm, Sim->tau_esc, Sim->B, Sim->S, Sim->array_len);
    save_data(file, Sim);
    while (Sim->t < Sim->end_t)
    {
        for (int32_t i = 0; i < Sim->n_species; i++)
        {
        
            implicit_step(Sim, Sim->Species[i]);

            equilibrium_check(Sim, Sim->Species[i]);

            if (Sim->end_sim)
            {
                Sim->t += Sim->dt;
                break;
            }
            save_step_to_prev_n(Sim, Sim->Species[i]);
            Sim->t += Sim->dt;
        }
        if (Sim->end_sim)
        {
            break;
        }
    }

    Sim->final_time=Sim->t;
    save_data(file, Sim);
    // flush remaining data in buffer to be written
    flush_buffer(file, Sim->buffer, Sim->buffer_index);
}

void write_gammas_to_file(FILE *file, SimulationParams *Sim)
{
    // print headers in csv file
    fprintf(file, "gamma,");
    for (int64_t i = 0; i < Sim->array_len; i++)
    {
        fprintf(file, "%lf,", Sim->Species[0]->gamma[i]);
    }
    fprintf(file, "\n");
    fflush(file);
}

void write_run_file(FILE *run_file, SimulationParams *Sim)
{
    if (run_file == NULL) {
        printf("Error: File pointer is NULL\n");
        return;
    }
    if (Sim == NULL) {
        printf("Error: SimulationParams pointer is NULL\n");
        return;
    }
    
    fprintf(run_file, 
    "dt,R,inject_p,inject_min,inject_max,rho,B,L,end_tol,Q_e0,S,tau_esc,norm,avg_gamma,V,array_len,max_gamma,min_gamma,init_p,samples_per_decade,final_time,change\n");
    fprintf(run_file,
    "%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%lld,%e,%e,%e,%lld,%e,%e\n",
    Sim->dt,Sim->R,Sim->inject_power,Sim->inject_min,Sim->inject_max,Sim->rho,Sim->B, 
    Sim->L,Sim->end_tol, Sim->Q_e0,Sim->S,Sim->tau_esc,Sim->norm,Sim->avg_gamma,Sim->V,
    Sim->array_len,Sim->max_gamma,Sim->min_gamma,Sim->init_power,Sim->samples_per_decade,
    Sim->final_time,Sim->change);
}

int main()
{
    SimulationParams *Sim = malloc(sizeof(SimulationParams));
    // simulation setup
    Sim->n_species = 1.;
    Sim->min_gamma = 1e1;
    Sim->max_gamma = 1e8;
    Sim->init_power = 2.;
    Sim->samples_per_decade = 10;
    Sim->dt = 10000.;
    Sim->end_t = 1e10;
    // free params
    Sim->inject_min = 1e4;
    Sim->inject_max = 1e8;
    Sim->inject_power = 2.3;
    Sim->B = 0.25;
    Sim->R = 1e16;
    Sim->L = 1e30;
    Sim->rho = 1e-38;
    Sim->end_tol = 1e-8;
    
    malloc_Sim_arrays(Sim);

    // generate the cooling test data
    double hold_B = Sim->B;
    double B[6] = {0.1,0.25,0.5,1.,1.5,2.};
    for (int i =0; i < 6; i++)
    {
        Sim->B = B[i];
        // generate file name based on B
        char filename[100], filepath[100], run_filepath[100];
        sprintf(filename, "B%4.0lf.csv", Sim->B*1000.);
        sprintf(filepath, "csv_data/%s", filename);
        sprintf(run_filepath, "csv_data/runs/run_%s", filename);
            
        FILE *file = fopen(filepath, "w");
        FILE *run_file = fopen(run_filepath, "w");

        // print gamma array in csv file as header
        write_gammas_to_file(file, Sim);

        simulate(file, Sim);

        write_run_file(run_file, Sim);
        fclose(file);
        fclose(run_file);
    }
    /*
    Sim->B = hold_B;
    FILE *file = fopen("csv_data/simulation_data.csv", "w");
    write_gammas_to_file(file, Sim);

    simulate(file, Sim);

    // flush remaining data in buffer to be written
    flush_buffer(file, Sim->buffer, Sim->buffer_index);
    fclose(file);
    */
    
    double param[6] = {5,10,20,40,80,120};
    for (int i =0; i < 6; i++)
    {
        Sim->samples_per_decade = param[i];
        malloc_and_fill_gamma_array(Sim, Sim->Species[0]);
        // generate file name based on B
        char filename[150], filepath[150], run_filepath[150];
        sprintf(filename, "samples_pd%.0lld.csv", Sim->samples_per_decade);
        sprintf(filepath, "csv_data/%s", filename);
        sprintf(run_filepath, "csv_data/runs/run_%s", filename);
            
        FILE *file = fopen(filepath, "w");
        FILE *run_file = fopen(run_filepath, "w");

        // print gamma array in csv file as header
        write_gammas_to_file(file, Sim);

        simulate(file, Sim);

        write_run_file(run_file, Sim);
        fclose(file);
        fclose(run_file);
    }
    
    // end program
    free_Sim_arrays(Sim);
    free(Sim);
    return 0;
}