#include <stdio.h>
#include <string.h>
#include <stdint.h>
#include <math.h>
#include <stdlib.h>

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
    double *delta_gamma;
} LeptonParams;

typedef struct SimulationParams
{
    int64_t array_len;
    double max_gamma;
    double min_gamma;
    double init_power;
    int64_t samples_per_decade;

    double t;
    double end_t;
    double dt;

    // free parameters
    double R;           // radius of system (spherical geo)
    double inject_power;// power law power
    double inject_min;  // turning point energy
    double inject_max;
    double rho;     // background density
    double B;       // background magnetic field
    double L;       // external luminosity

    double Q_e0;    // starting population
    double S;       // sync vs inv. compton
    double tau_esc; // free escape time
    double norm;    // normalize the prob dist
    int32_t n_species;
    LeptonParams **Species;

    char *buffer;
    size_t *buffer_index;

} SimulationParams;

void malloc_and_fill_gamma_array(SimulationParams *Sim, LeptonParams *Lepton)
{
    int64_t decades;

    decades = (int64_t)log10(Sim->max_gamma) - log10(Sim->min_gamma);
    Sim->array_len = decades * Sim->samples_per_decade;

    printf("%lld\n", Sim->array_len);
    Lepton->gamma = calloc(Sim->array_len + 1, sizeof(double));
    Lepton->delta_gamma = calloc(Sim->array_len + 1, sizeof(double));
    for (int64_t i = 0; i < Sim->array_len + 1; i++)
    {
        Lepton->gamma[i] = pow(10, log10(Sim->min_gamma) + i / (double)Sim->samples_per_decade);
    }

    for (int64_t i = Sim->array_len; i >= 0; i--)
    {
        Lepton->delta_gamma[i] = Lepton->gamma[i + 1] - Lepton->gamma[i];
    }
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
        Sim->Species[i]->next_n = calloc(Sim->array_len + 1, sizeof(double));
        Sim->Species[i]->current_n = calloc(Sim->array_len + 1, sizeof(double));
    }
}

void free_Sim_Species_array(SimulationParams *Sim)
{
    for (int32_t i = 0; i < Sim->n_species; i++)
    {
        free(Sim->Species[i]->next_n);
        free(Sim->Species[i]->current_n);
        free(Sim->Species[i]->gamma);
        free(Sim->Species[i]->delta_gamma);
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
        Lepton->next_n[i] =
        Sim->tau_esc * 
        (Sim->S * Sim->dt * Lepton->gamma[i+1] * Lepton->gamma[i+1] * Lepton->next_n[i+1]
        - Lepton->delta_gamma[i] * Sim->dt * I(Lepton->gamma[i], Sim->inject_min, Sim->inject_max, Sim->inject_power, Sim)
        - Lepton->delta_gamma[i] * Lepton->current_n[i])
        /
        (Sim->S * Sim->tau_esc * Sim->dt * Lepton->gamma[i] * Lepton->gamma[i]
        - Lepton->delta_gamma[i] * Sim->tau_esc
        - Lepton->delta_gamma[i] * Sim->dt);

        // copy new value to previous array
        Lepton->current_n[i] = Lepton->next_n[i];
    }
}

void set_initial_state(SimulationParams *Sim)
{

    for (int32_t lepton = 0; lepton < Sim->n_species; lepton++)
    {
        Sim->Species[lepton]->next_n[Sim->array_len] = 0.;
        for (int64_t i = 0; i < Sim->array_len; i++)
        {
            //Sim->Species[lepton]->current_n[i] = 0.;
            Sim->Species[lepton]->current_n[i] =(Sim->rho / (2.*m_e)) * (I(Sim->Species[lepton]->gamma[i],Sim->min_gamma, Sim->max_gamma, Sim->init_power, Sim) / Sim->Q_e0);
        }
    }
}

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

    for (int i = 0; i < Sim->array_len+1; i++)
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

void simulate(FILE *file, SimulationParams *Sim)
{
    set_initial_state(Sim);
    Sim->t = 0.;
    save_data(file, Sim);
    while (Sim->t < Sim->end_t)
    {
        for (int32_t i = 0; i < Sim->n_species; i++)
        {
            iteration(Sim, Sim->Species[i]);
            Sim->t += Sim->dt;
        }
    }
    save_data(file, Sim);
}

void calc_S(SimulationParams *Sim)
{
    // calculate S
    Sim->S = - 2.9e-9 * (Sim->B * Sim->B);
}

int main()
{
    SimulationParams *Sim = malloc(sizeof(SimulationParams));
    Sim->n_species = 1.;
    Sim->min_gamma = 1e1;
    Sim->max_gamma = 1e8;
    Sim->init_power = 2.;
    Sim->samples_per_decade = 40;
    Sim->dt = 2500.;
    Sim->end_t = 1e7;

    Sim->inject_min = 1e4;
    Sim->inject_max = 1e8;
    Sim->inject_power = 2.3;
    Sim->R = 1e16;
    Sim->L = 1e30;
    Sim->rho = 1e-38;
    malloc_Sim_Species_arrays(Sim);

    // calculate S
    calc_S(Sim);

    Sim->norm =
        1. /
        ((pow(Sim->max_gamma, 1. - Sim->inject_power) / (1. - Sim->inject_power)) - (pow(Sim->min_gamma, 1. - Sim->inject_power) / (1. - Sim->inject_power)));
    
    // calculate inputs
    double average_gamma =
        Sim->norm * (pow(Sim->inject_max, 2. - Sim->inject_power) / (2. - Sim->inject_power)) - Sim->norm * (pow(Sim->inject_min, 2. - Sim->inject_power) / (2. - Sim->inject_power));
    if (Sim->inject_power == 2.)
        average_gamma = Sim->norm * log(Sim->inject_max) - Sim->norm * log(Sim->inject_min);
    
    double volume = (4./3.) * M_PI * Sim->R * Sim->R * Sim->R;
    // calculate injected Qe0
    Sim->Q_e0 = Sim->L / (volume * average_gamma * m_e * c * c);
    // calculate free escape time 
    Sim->tau_esc = (3. / 4.) * (Sim->R / c);

    printf("C %e, tau %e\n", Sim->Q_e0 * Sim->norm, Sim->tau_esc);
    
    double B[6] = {0.1,0.25,0.5,1.,1.5,2.};
    for (int i =0; i < 6; i++)
    {
        Sim->B = B[i];
        calc_S(Sim);
        printf("B %.2lf, S %e\n", Sim->B, Sim->S);
        char filename[20];
        sprintf(filename, "csv_data/B%4.0lf.csv", Sim->B*1000.);
        
        FILE *file = fopen(filename, "w");
        // print headers in csv file
        fprintf(file, "gamma,");
        for (int64_t i = 0; i < Sim->array_len+1; i++)
        {
            fprintf(file, "%lf,", Sim->Species[0]->gamma[i]);
        }
        fprintf(file, "\n");
        fflush(file);

        simulate(file, Sim);

        // flush remaining data in buffer to be written
        flush_buffer(file, Sim->buffer, Sim->buffer_index);
        fclose(file);
    }
    
    Sim->B = 0.1;
    FILE *file = fopen("csv_data/simulation_data.csv", "w");
    // print headers in csv file
    fprintf(file, "gamma,");
    for (int64_t i = 0; i < Sim->array_len; i++)
    {
        fprintf(file, "%lf,", Sim->Species[0]->gamma[i]);
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