#include <stdio.h>
#include <string.h>
#include <stdint.h>
#include <math.h>
#include <stdlib.h>

#define BUFFER_SIZE 33554432

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
    int64_t samples_per_decade;

    double t;
    double end_t;
    double dt;

    // free parameters
    double tau_esc; // free escape time
    double S;       // sync vs invers compton
    double C;       // power law constant
    double power;   // power law power
    double gamma_m; // turning point energy

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

double I(double gamma, SimulationParams *Sim)
{
    if (gamma < Sim->gamma_m)
    {
        return 0.;
    }
    else
    {
        return Sim->C * pow(gamma, (-1. * Sim->power));
    }
}

void iteration(SimulationParams *Sim, LeptonParams *Lepton)
{
    for (int64_t i = Sim->array_len; i >= 0; i--)
    {
        Lepton->next_n[i] =
            (1 / (1 - Sim->dt * ((Sim->S * Lepton->gamma[i] * Lepton->gamma[i] / Lepton->delta_gamma[i]) - (1 / Sim->tau_esc)))) *
            (Lepton->current_n[i] +
             Sim->dt * I(Lepton->gamma[i], Sim) -
             (Sim->S * Sim->dt * Lepton->gamma[i + 1] * Lepton->gamma[i + 1] / Lepton->delta_gamma[i]) * Lepton->next_n[i + 1]);

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
            Sim->Species[lepton]->current_n[i] = 1.;
            Sim->Species[lepton]->next_n[i] = 1.;
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

    for (int i = 0; i < Sim->array_len; i++)
    {
        required_space += snprintf(Sim->buffer + (*Sim->buffer_index) + required_space, buffer_space - required_space,
                                   "%lf,",
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

int main()
{
    SimulationParams *Sim = malloc(sizeof(SimulationParams));
    Sim->n_species = 1.;
    Sim->min_gamma = 0.1;
    Sim->max_gamma = 1e9;
    Sim->samples_per_decade = 10;
    Sim->dt = 1.;
    Sim->end_t = 1000.;

    Sim->gamma_m = 1e4;
    Sim->C = -1e30;
    Sim->power = 2;
    Sim->S = 5.45;
    Sim->tau_esc = 10.;
    malloc_Sim_Species_arrays(Sim);

    FILE *file = fopen("simulation_data.csv", "w");
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