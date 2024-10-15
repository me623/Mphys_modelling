#include <stdio.h>
#include <string.h>
#include <stdint.h>
#include <math.h>
#include <stdlib.h>

#define BUFFER_SIZE 33554432
// for lambert W
#define MAX_ITER 100
#define EPSILON 1e-10

// CGS unit constants
#define m_e 9.109e-28
#define c 3e10
#define sigma_t 6.65e-25

typedef struct LeptonParams
{
    double *current_lnn;
    double *next_lnn;
    double *g;
    double *ln_g;
    double dln_g;
} LeptonParams;

typedef struct SimulationParams
{
    int64_t array_len;
    double max_gamma;
    double min_gamma;
    int64_t samples_per_decade;

    double t;
    double end_t;
    double dt;

    // free parameters
    double tau_esc;    // free escape time
    double R;          // radius of system (spherical geo)
    double inject_power;      // power law power
    double init_power; //
    double inject_min; // 
    double inject_max; //
    double rho;        // background density
    double B;          // background magnetic field
    double L;          // external luminosity
    double boundary_val;

    double Q_e0; // starting population
    double S;    // sync vs inv. compton
    double norm; // normalize the prob dist

    // temp values
    double lambertW;

    int32_t n_species;
    LeptonParams **Species;

    char *buffer;
    size_t *buffer_index;

} SimulationParams;

// Lambert W function for the principal branch W_0(x)
double lambertW(double x)
{
    if (x == 0.0)
    {
        return 0.0; // W(0) = 0
    }

    // Initial guess: log(x) for large x, or just x for smaller values
    double w = (x < 1.0) ? x : log(x);

    // Iterate using the Newton-Raphson method
    for (int i = 0; i < MAX_ITER; i++)
    {
        double ew = exp(w);
        double w_next = w - (w * ew - x) / (ew * (w + 1) - (w + 2) * (w * ew - x) / (2 * (w + 1))); // Halley's method for faster convergence

        // Check for convergence
        if (fabs(w_next - w) < EPSILON)
        {
            return w_next;
        }

        w = w_next;
    }

    // If it doesn't converge, return the last value computed
    return w;
}

void malloc_and_fill_gamma_array(SimulationParams *Sim, LeptonParams *Lepton)
{
    int64_t decades;

    // calculate number of decades in the gamma range
    decades = (int64_t)log10(Sim->max_gamma) - log10(Sim->min_gamma);
    // set array length based on decades and input samples per decade
    Sim->array_len = decades * Sim->samples_per_decade;

    printf("%lld\n", Sim->array_len);
    // malloc the gamma and ln_gamma arrays
    Lepton->g = calloc(Sim->array_len + 1, sizeof(double));
    Lepton->ln_g = calloc(Sim->array_len + 1, sizeof(double));

    // set the gamma values
    for (int64_t i = 0; i < Sim->array_len + 1; i++)
    {
        Lepton->g[i] = pow(10., log10(Sim->min_gamma) + (double)i / (double)Sim->samples_per_decade);
        Lepton->ln_g[i] = log(Lepton->g[i]);
    }
    // calculate delta ln gamma
    Lepton->dln_g = Lepton->ln_g[1] - Lepton->ln_g[0];
}

void malloc_Sim_Species_arrays(SimulationParams *Sim)
{
    Sim->Species = malloc(Sim->n_species * sizeof(LeptonParams *));
    Sim->buffer = malloc(BUFFER_SIZE * sizeof(char));
    Sim->buffer_index = malloc(sizeof(*(Sim->buffer_index)));
    *Sim->buffer_index = 0;

    for (int32_t i = 0; i < Sim->n_species; i++)
    {
        Sim->Species[i] = malloc(sizeof(LeptonParams));
        malloc_and_fill_gamma_array(Sim, Sim->Species[i]);
        Sim->Species[i]->next_lnn = calloc(Sim->array_len + 1, sizeof(double));
        Sim->Species[i]->current_lnn = calloc(Sim->array_len + 1, sizeof(double));
    }
}

void free_Sim_Species_array(SimulationParams *Sim)
{
    for (int32_t i = 0; i < Sim->n_species; i++)
    {
        free(Sim->Species[i]->next_lnn);
        free(Sim->Species[i]->current_lnn);
        free(Sim->Species[i]->g);
        free(Sim->Species[i]->ln_g);
        free(Sim->Species[i]);
    }
    free(Sim->buffer);
    free(Sim->buffer_index);
    free(Sim->Species);
}

double I(double gamma, double min, double max, double power, SimulationParams *Sim)
{
    // inject power law
    if (gamma < min || gamma > max)
    {
        return 0.;
    }
    else
    {
        return Sim->Q_e0 * Sim->norm * pow(gamma, (-1. * power));
    }
}

void iteration(SimulationParams *Sim, LeptonParams *Lepton)
{
    for (int64_t i = Sim->array_len - 1; i >= 0; i--)
    {
        // pre calculate lambertW
        Sim->lambertW = lambertW(
            -1. * ((I(Lepton->g[i], Sim->inject_min, Sim->inject_max, Sim->inject_power, Sim) * Sim->dt * Lepton->dln_g * exp(
                -1. * ((2. * Sim->S * Sim->tau_esc * Sim->dt * Lepton->g[i] * Lepton->dln_g 
                + Sim->S * Sim->tau_esc * Sim->dt * Lepton->g[i] * Lepton->next_lnn[i + 1] 
                - Sim->tau_esc * Lepton->current_lnn[i] * Lepton->dln_g 
                + Lepton->dln_g * Sim->dt) 
                / 
                (Sim->tau_esc * (Sim->S * Sim->dt * Lepton->g[i] - Lepton->dln_g))))) 
            /
            (Sim->S * Sim->dt * Lepton->g[i] - Lepton->dln_g)));

        // stepping regime
        Lepton->next_lnn[i] =
            (1. / (Sim->tau_esc * (Sim->S * Sim->dt * Lepton->g[i] - Lepton->dln_g))) 
            * (Sim->lambertW * Sim->S * Sim->tau_esc * Sim->dt * Lepton->g[i] 
            + Sim->S * Sim->tau_esc * Sim->dt * Lepton->g[i] * Lepton->next_lnn[i + 1] 
            + 2. * Sim->S * Sim->tau_esc * Sim->dt * Lepton->g[i] * Lepton->dln_g 
            - Sim->lambertW * Sim->tau_esc * Lepton->dln_g 
            - Sim->tau_esc * Lepton->current_lnn[i] * Lepton->dln_g 
            + Lepton->dln_g * Sim->dt);
        
        // copy new value to previous array
        Lepton->current_lnn[i] = Lepton->next_lnn[i];
    }
}

void set_initial_state(SimulationParams *Sim)
{

    for (int32_t lepton = 0; lepton < Sim->n_species; lepton++)
    {
        Sim->Species[lepton]->next_lnn[Sim->array_len] = Sim->boundary_val;
        for (int64_t i = 0; i < Sim->array_len; i++)
        {
            // set initial power law dist based on the background density
            Sim->Species[lepton]->current_lnn[i] = log((Sim->rho / m_e) * (I(Sim->Species[lepton]->g[i], Sim->min_gamma, Sim->max_gamma, Sim->init_power, Sim) / Sim->Q_e0));
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
                                   Sim->Species[0]->current_lnn[i]);
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

void simulate(FILE *file, SimulationParams *Sim)
{
    set_initial_state(Sim);
    save_data(file, Sim);
    Sim->t = 0.;
    while (Sim->t < Sim->end_t)
    {
        for (int32_t i = 0; i < Sim->n_species; i++)
        {
            iteration(Sim, Sim->Species[i]);
            Sim->t += Sim->dt;
            save_data(file, Sim);
        }
    }
}

void calc_S(SimulationParams *Sim)
{
    // calculate S
    Sim->S = -(4. / 3.) * c * sigma_t * (Sim->B * Sim->B / (2 * m_e * c * c));
}

int main()
{
    SimulationParams *Sim = malloc(sizeof(SimulationParams));
    // base sim params
    Sim->n_species = 1.;
    Sim->min_gamma = 1e1;
    Sim->max_gamma = 1e8;
    Sim->init_power=2.;
    Sim->samples_per_decade = 20;
    Sim->dt = 10000.;
    Sim->end_t = 10000000.;
    // free parameters
    Sim->inject_min = 1e4;
    Sim->inject_max = 1e8;
    Sim->boundary_val = 1e12;
    Sim->inject_power = 2.3;
    Sim->R = 1e16;
    Sim->L = 1e30;
    Sim->B = 0.01;
    Sim->rho = 1e-38;
    malloc_Sim_Species_arrays(Sim);

    calc_S(Sim);

    Sim->norm =
        1. /
        ((pow(Sim->max_gamma, 1. - Sim->inject_power) / (1. - Sim->inject_power)) - (pow(Sim->min_gamma, 1. - Sim->inject_power) / (1. - Sim->inject_power)));


    // calculate inputs
    double average_gamma =
        Sim->norm * (pow(Sim->inject_max, 2. - Sim->inject_power) / (2. - Sim->inject_power)) - Sim->norm * (pow(Sim->inject_min, 2. - Sim->inject_power) / (2. - Sim->inject_power));
    if (Sim->inject_power == 2.)
        average_gamma = Sim->norm * log(Sim->inject_max) - Sim->norm * log(Sim->inject_min);

    // find volume
    double volume = (4. / 3.) * M_PI * Sim->R * Sim->R * Sim->R;
    // calculate injected Qe0
    Sim->Q_e0 = (Sim->L / (volume * average_gamma * m_e * c * c));
    // calculate free escape time
    Sim->tau_esc = (3. / 4.) * (Sim->R / c);
    
    printf("norm_f%e, Q_e0%e, C%e, V%e, tau_esc%e \n", Sim->norm, Sim->Q_e0, Sim->Q_e0 * Sim->norm, volume, Sim->tau_esc);
    /*
    double B[6] = {0.1,0.25,0.5,1.,1.5,2.};
    for (int i =0; i < 6; i++)
    {
        Sim->B = B[i];
        calc_S(Sim);
        printf("%lf\n", Sim->B);
        char filename[20];
        sprintf(filename, "B%3.0lf.csv", Sim->B*100.);
        
        FILE *file = fopen(filename, "w");
        // print headers in csv file
        fprintf(file, "gamma,");
        for (int64_t i = 0; i < Sim->array_len; i++)
        {
            fprintf(file, "%lf,", Sim->Species[0]->g[i]);
        }
        fprintf(file, "\n");
        fflush(file);

        simulate(file, Sim);

        // flush remaining data in buffer to be written
        flush_buffer(file, Sim->buffer, Sim->buffer_index);
        fclose(file);
    }
    */

    FILE *file = fopen("log_step_simulation_data.csv", "w");
    // print headers in csv file
    fprintf(file, "gamma,");
    for (int64_t i = 0; i < Sim->array_len; i++)
    {
        fprintf(file, "%lf,", Sim->Species[0]->g[i]);
    }
    fprintf(file, "\n");
    fflush(file);

    simulate(file, Sim);

    // flush remaining data in buffer to be written
    flush_buffer(file, Sim->buffer, Sim->buffer_index);
    fclose(file);
    free_Sim_Species_array(Sim);
    free(Sim);

    return 0;
}