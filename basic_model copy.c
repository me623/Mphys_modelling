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

    // physical parameters
    double t_esc;
    double t_dec;
    double *t_acc;
    double S; // sync vs invers compton
} LeptonParams;

typedef struct SimulationParams
{
    int64_t array_len;
    double delta_gamma;
    double delta_n;
    double t;
    double end_t;
    double dt;
    double free_escape_time;

    double background_density;
    double magnetic_field;
    double luminosity;
    double radius;

    double t_acc;
    double L;

    int32_t n_species;
    LeptonParams **Species;

} SimulationParams;

void malloc_Sim_Species_arrays(SimulationParams *Sim)
{
    Sim->Species = malloc(Sim->n_species * sizeof(LeptonParams *));
    
    for (int32_t i = 0; i < Sim->n_species; i++)
    {
        Sim->Species[i] = malloc(sizeof(LeptonParams));
        Sim->Species[i]->next_n = calloc(Sim->array_len, sizeof(int64_t));
        Sim->Species[i]->prev_n = calloc(Sim->array_len, sizeof(int64_t));
        Sim->Species[i]->gamma = calloc(Sim->array_len, sizeof(double));
        Sim->Species[i]->t_acc = malloc(sizeof(double *));
        Sim->Species[i]->t_acc = &Sim->t_acc;
    }
}

double power_law_generator(double uniform_deviate, double start, double finish, double power)
{
    return pow(((pow(finish, (power + 1)) - pow(start, (power + 1))) * uniform_deviate + pow(start, (power + 1))), (1 / (power + 1)));
}

double fill_distribution_power_law(int64_t *array, SimulationParams *Sim, int64_t N_particles, double start, double finish, double power)
{
    double uniform_deviate; 
    int64_t index, i = 0; 
    // clear array for binning the data too
    memset(array, 0 , Sim->array_len * sizeof(int64_t));
    
    // step over all required 
    while (i < N_particles)
    {
        uniform_deviate = (double)rand() / (double)RAND_MAX;
        index = (int64_t)power_law_generator(uniform_deviate, start, finish, power);
        if (Sim->array_len > index)
        {
            array[index]++;
            i++;
        }
    }

}

double A_prime(double gamma, SimulationParams *Sim, LeptonParams *Lepton)
{
    return Sim->L - 1. / Lepton->t_esc - 1. / Lepton->t_dec - 1. / Lepton->t_acc - 2. * Lepton->S * gamma;
}

double B_prime(double gamma, LeptonParams *Lepton)
{
    return -gamma / *Lepton->t_acc - Lepton->S * gamma * gamma;
}

void iteration(SimulationParams *Sim, LeptonParams *Lepton)
{
    for (int64_t i = Sim->array_len; i == 0; i--)
    {
        Lepton->next_n[i] = pow(exp(1.),
                            (log(Sim->delta_gamma) / B_prime(Lepton->gamma[i], Lepton)) *
                            ((1. / Lepton->prev_n[i]) + A_prime(Lepton->gamma[i], Sim, Lepton)) + // removed injection term Q for now
                            log(Lepton->next_n[i + 1]));

        // copy new value to previous array
        Lepton->prev_n[i] = Lepton->next_n[i];
    }
}

void set_initial_state(SimulationParams *Sim)
{
    
    for (int32_t i = 0; i < Sim->n_species; i++)
    {

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