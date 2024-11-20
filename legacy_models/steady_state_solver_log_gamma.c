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
    double *n;
    double *prev_n;
    double *gamma;
    double *ln_gamma;
    double delta_ln_gamma;
} LeptonParams;

typedef struct SimulationParams
{
    // array parameters
    int64_t array_len;
    double max_gamma;
    double min_gamma;
    int64_t samples_per_decade;

    // free parameters
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
    int64_t iter;
    int64_t max_iter;
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
    Lepton->gamma = calloc(Sim->array_len + 2, sizeof(double));
    // malloc gamma and delta gamma arrays
    Lepton->ln_gamma = calloc(Sim->array_len + 2, sizeof(double));

    // fill gamma array with equal log step data
    for (int64_t i = 1; i <= Sim->array_len + 1; i++)
    {
        Lepton->gamma[i] = pow(10, log10(Sim->min_gamma) + (i-1) / (double)Sim->samples_per_decade);
        Lepton->ln_gamma[i] = log(Lepton->gamma[i]);
    }

    // Calculate the step size in logs
    double log_step = 1.0 / (double)Sim->samples_per_decade;
    // Extrapolate gamma[0] based on gamma[1]
    Lepton->gamma[0] = pow(10, log10(Sim->min_gamma) - log_step);
    Lepton->ln_gamma[0] = log(Lepton->gamma[0]);

    Lepton->delta_ln_gamma = Lepton->ln_gamma[1] - Lepton->ln_gamma[0];
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
        Sim->Species[i]->n = malloc((Sim->array_len + 2) * sizeof(double));
        Sim->Species[i]->prev_n = malloc((Sim->array_len + 2) * sizeof(double));
    }
}

void free_Sim_arrays(SimulationParams *Sim)
{
    for (int32_t i = 0; i < Sim->n_species; i++)
    {
        if (Sim->Species[i] != NULL)
        {
            if (Sim->Species[i]->n != NULL)
                free(Sim->Species[i]->n);
            if (Sim->Species[i]->prev_n != NULL)
                free(Sim->Species[i]->prev_n);
            if (Sim->Species[i]->gamma != NULL)
                free(Sim->Species[i]->gamma);
            if (Sim->Species[i]->ln_gamma != NULL)
                free(Sim->Species[i]->ln_gamma);
            
            free(Sim->Species[i]);
        }
    }

    if (Sim->buffer != NULL)
        free(Sim->buffer);

    if (Sim->buffer_index != NULL)
        free(Sim->buffer_index);

    if (Sim->Species != NULL)
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

void normalize_inject_dist(double power, SimulationParams *Sim)
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
        for (int64_t i = 1; i <= Sim->array_len; i++)
        {
            // re-normalize dist for init power
            normalize_inject_dist(Sim->init_power, Sim);
            
            // set initial population on a selected power law
            // number of photons included based on background density
            Sim->Species[lepton]->n[i] =
            (Sim->rho / (2.*m_e)) * 
            (I(Sim->Species[lepton]->gamma[i], Sim->min_gamma, Sim->max_gamma, Sim->init_power, Sim) 
            / Sim->Q_e0);
            //printf("%e %e\n", Sim->Species[lepton]->gamma[i], Sim->Species[lepton]->n[i]);
            
            // get ready for injection
            normalize_inject_dist(Sim->inject_power, Sim);
            //Sim->Species[lepton]->n[i] = 0.;
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
                                  "");

    for (int i = 0; i < Sim->array_len + 2; i++)
    {
        required_space += snprintf(Sim->buffer + (*Sim->buffer_index) + required_space, buffer_space - required_space,
                                   "%e,",
                                   Sim->Species[0]->n[i]);
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

void cda_step(SimulationParams *Sim, LeptonParams *Lepton)
{
    for (int64_t i = 1; i <= Sim->array_len; i++)
    {
        // explicit stepping regime
        Lepton->n[i] = 
        (Sim->tau_esc * 
        (Sim->S * Lepton->gamma[i] * Lepton->n[i-1]
        - Sim->S * Lepton->gamma[i] * Lepton->n[i+1]
        + 2. * Lepton->delta_ln_gamma * I(Lepton->gamma[i], Sim->inject_min, Sim->inject_max, Sim->inject_power, Sim)))
         /
        (2. * Lepton->delta_ln_gamma * (1. + 2. * Sim->S * Sim->tau_esc));
    }
}

void fda_step(SimulationParams *Sim, LeptonParams *Lepton)
{
    for (int64_t i =  Sim->array_len; i >= 1; i--)
    {
        // explicit stepping regime
        Lepton->n[i] = 
        (Sim->tau_esc * 
        (Sim->S * Lepton->gamma[i+1] * Lepton->gamma[i+1] * Lepton->n[i+1]
        - 2. * Lepton->gamma[i] * Lepton->delta_ln_gamma 
        * I(Lepton->gamma[i], Sim->inject_min, Sim->inject_max, Sim->inject_power, Sim)))
         /
        (Lepton->gamma[i] * (Sim->S * Sim->tau_esc * Lepton->gamma[i] - Lepton->delta_ln_gamma));
    }
}

void save_step_to_prev_n(SimulationParams *Sim, LeptonParams *Lepton)
{
    for (int64_t i = 0; i <= Sim->array_len + 1; i++)
    {
        Lepton->prev_n[i] = Lepton->n[i];
    }
}

void equilibrium_check(SimulationParams *Sim, LeptonParams *Lepton)
{
    Sim->change = 0.;
    // calculate percentage change in n
    for (int64_t i = 1; i <= Sim->array_len; i++)
    {
        Sim->change += pow(
            (1. / Lepton->n[i]) * ((Lepton->n[i] - Lepton->prev_n[i])),
             2.);
    }
    Sim->change = sqrt(Sim->change);
    
    // check change in population against specified end tolerance
    if (Sim->change < Sim->end_tol)
    {
        printf("equilibrium reached at step:%lld, last change: %e\n", Sim->iter, Sim->change);
        Sim->end_sim=true;
    }
}

void impose_BCs(SimulationParams *Sim,LeptonParams *Lepton)
{
    Lepton->n[0] = Lepton->n[1];
    Lepton->n[Sim->array_len+1] = 0.;
}

void simulate(FILE *file, SimulationParams *Sim)
{
    // take input params and calculate coefficients
    calc_S(Sim);
    normalize_inject_dist(Sim->inject_power, Sim);
    calc_Q_e0(Sim);
    calc_tau_esc(Sim);
    
    set_initial_state(Sim);
    impose_BCs(Sim, Sim->Species[0]);

    // start simulation
    printf("Start Sim with C %e, tau %e B %.2lf, S %e g_array %lld\n", 
    Sim->Q_e0 * Sim->norm, Sim->tau_esc, Sim->B, Sim->S, Sim->array_len);
    save_data(file, Sim);
    Sim->end_sim = false;
    Sim->iter = 0;
    while (Sim->end_sim == false && Sim->iter < Sim->max_iter)
    {
        for (int32_t i = 0; i < Sim->n_species; i++)
        {
            save_step_to_prev_n(Sim, Sim->Species[i]);
            fda_step(Sim, Sim->Species[i]);
            equilibrium_check(Sim, Sim->Species[i]);
            impose_BCs(Sim, Sim->Species[i]);

            if (Sim->end_sim)
            {
                break;
            }
        }
        Sim->iter ++;
    }

    save_data(file, Sim);
    // flush remaining data in buffer to be written
    flush_buffer(file, Sim->buffer, Sim->buffer_index);
}

void write_gammas_to_file(FILE *file, SimulationParams *Sim)
{
    // print headers in csv file
    fprintf(file, "gamma,");
    for (int64_t i = 0; i < Sim->array_len+2; i++)
    {
        fprintf(file, "%lf,", Sim->Species[0]->gamma[i]);
    }
    fprintf(file, "\n");
    fflush(file);
}

void write_run_file(SimulationParams *Sim, char *filename)
{
    char filepath[100];
    sprintf(filepath, "csv_data/steady_state/runs/run_%s", filename);

    FILE *file = fopen(filepath, "w");
    fprintf(file, 
    "delta_ln_gamma,R,inject_p,inject_min,inject_max,rho,B,L,end_tol,Q_e0,S,tau_esc,norm,avg_gamma,V,array_len,max_gamma,min_gamma,init_p,samples_per_decade,change\n");
    fprintf(file,
    "%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%lld,%e,%e,%e,%lld,%e\n",
    Sim->Species[0]->delta_ln_gamma, Sim->R,Sim->inject_power,Sim->inject_min,Sim->inject_max,Sim->rho,Sim->B,Sim->L,
    Sim->end_tol, Sim->Q_e0,Sim->S,Sim->tau_esc,Sim->norm,Sim->avg_gamma,Sim->V,Sim->array_len,
    Sim->max_gamma,Sim->min_gamma,Sim->init_power,Sim->samples_per_decade, Sim->change);
    fclose(file);
}

int main()
{
    SimulationParams *Sim = malloc(sizeof(SimulationParams));
    // simulation setup
    Sim->n_species = 1.;
    Sim->min_gamma = 1e1;
    Sim->max_gamma = 1e8;
    Sim->init_power = 2.;
    Sim->samples_per_decade = 100;
    // free params
    Sim->inject_min = 1e4;
    Sim->inject_max = 1e8;
    Sim->inject_power = 2.3;
    Sim->B = 2.;
    Sim->R = 1e16;
    Sim->L = 1e30;
    Sim->rho = 1e-38;
    Sim->end_tol = 1e-8;
    Sim->max_iter = 1000;
    
    malloc_Sim_arrays(Sim);

    // generate the cooling test data
    
    double hold_B = Sim->B;
    double B[6] = {0.1,0.25,0.5,1.,1.5,2.};
    for (int i =0; i < 6; i++)
    {
        Sim->B = B[i];
        // generate file name based on B
        char filename[100], filepath[100];
        sprintf(filename, "B%4.0lf.csv", Sim->B*1000.);
        sprintf(filepath, "csv_data/steady_state/%s", filename);

        FILE *file = fopen(filepath, "w");
        // print gamma array in csv file as header
        write_gammas_to_file(file, Sim);

        simulate(file, Sim);

        write_run_file(Sim, filename);
        fclose(file);
    }
    
    Sim->B = hold_B;
    /*
    FILE *file2 = fopen("csv_data/steady_state/simulation_data.csv", "w");
    write_gammas_to_file(file2, Sim);

    simulate(file2, Sim);

    // flush remaining data in buffer to be written
    flush_buffer(file2, Sim->buffer, Sim->buffer_index);
    write_run_file(Sim);
    fclose(file2);
    */
    // end program
    free_Sim_arrays(Sim);
    free(Sim);
    return 0;
}