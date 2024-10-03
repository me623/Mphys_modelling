#include <stdio.h>
#include <string.h>
#include <stdint.h>
#include <math.h>
#include <stdlib.h>

typedef struct LeptonParams
{
    int64_t *prev_n;
    int64_t *next_n;
    double *gamma;
    double *delta_gamma;
} LeptonParams;

typedef struct SimulationParams
{
    int64_t array_len;

    double t;
    double end_t;
    double dt;
    
    double tau_esc; // free escape time
    double S;   // sync vs invers compton
    double C;   // power law constant
    double power;   // power law power
    double gamma_m; // turning point energy

    int32_t n_species;
    LeptonParams **Species;

} SimulationParams;

void malloc_Sim_Species_arrays(SimulationParams *Sim)
{
    Sim->Species = malloc(Sim->n_species * sizeof(LeptonParams *));
    
    for (int32_t i = 0; i < Sim->n_species; i++)
    {
        Sim->Species[i] = malloc(sizeof(LeptonParams));
        Sim->Species[i]->next_n = calloc(Sim->array_len + 1, sizeof(int64_t));
        Sim->Species[i]->prev_n = calloc(Sim->array_len + 1, sizeof(int64_t));
        Sim->Species[i]->gamma = calloc(Sim->array_len, sizeof(double));

    }
}

double I(double gamma, SimulationParams *Sim)
{
    if (gamma < Sim->gamma_m)
    {
        return 0.;
    }
    else
    {
        return Sim->C * gamma * exp(-1. * Sim->power);
    }
}

void iteration(SimulationParams *Sim, LeptonParams *Lepton)
{
    for (int64_t i = Sim->array_len; i == 0; i--)
    {
        Lepton->next_n[i] = (1. / (2. * Sim->S * Lepton->gamma[i] + (1. / Sim->tau_esc) - ((Sim->S * Lepton->gamma[i] * Lepton->gamma[i]) / Lepton->delta_gamma[i]))) *
        (I(Lepton->gamma[i], Sim) - ((Sim->S * Lepton->gamma[i] * Lepton->gamma[i]) / Lepton->delta_gamma[i]) * Lepton->next_n[i+1]);

        // copy new value to previous array
        Lepton->prev_n[i] = Lepton->next_n[i];
    }
}

void set_initial_state(SimulationParams *Sim)
{
    
    for (int32_t i = 0; i < Sim->n_species; i++)
    {
        for (int64_t j = 0; j < Sim->array_len; j++)
        {
            Sim->Species[i]->next_n[Sim->array_len] = 0;
        }
    }
}

void simulate(SimulationParams *Sim)
{
    set_initial_state(Sim);

    for (Sim->t = 0.; Sim->t < Sim->end_t; Sim->t += Sim->dt)
    {
        for (int32_t i = 0; i < Sim->n_species; i++)
        {
            iteration(Sim, Sim->Species[i]);
        }
    }
}

int main()
{

}