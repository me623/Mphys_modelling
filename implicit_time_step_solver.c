
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <stdbool.h>
#include <sys/time.h>
#include <gsl/gsl_sf_hyperg.h>
#include <gsl/gsl_sf_bessel.h> // GNU Scientific Library for Bessel functions

// code settings
#define WRITE_STEPS 0
#define TOLERANCE 1e-6  // Tolerance for numerical integration
#define INFINITY_CUTOFF 100.0  // Approximate "infinity" for inner integral
#define APPROXIMATE_XCS 1   // use approximation avoiding long Whitaker function calculation

// CGS unit constants
#define m_e 9.109e-28
#define m_p 1.6726e-24
#define c 2.998e10
#define sigma_t 6.652e-25
#define q 4.803e-10
#define h 6.626e-27
#define H0 2.27e-18 // hubble constant in CGS
#define h_bar 1.055e-27
#define FINE_STRUCTURE_CONSTANT (1/137.035999)
#define CFE_RATIO 1
#define ELECTRON_ENERGY (m_e * c * c)
#define ELECTRON_RADIUS (q * q / ELECTRON_ENERGY)

typedef struct ElectronParams
{
    double *current_n;
    double *next_n;
    double *gamma;
    double *dgamma_fwd;
    double dln_gamma;
    double tau_esc;
} ElectronParams;

typedef struct PhotonParams
{
    double *n;
    double *eps;
    double *emission;
    double *absorption;
    double *xCS;
    double tau_esc;
} PhotonParams;

typedef struct SimulationParams
{
    /*** free parameters ***/

    // array parameters
    double max_gamma;           // max lorentz factor for electron population
    double min_gamma;           // min ___
    double max_eps;             // max normalised photon energy
    double min_eps;             // min ___
    int64_t samples_per_decade; // number of gamma samples to generate per decade

    // physical parameters
    double R;              // radius of system (spherical geo)
    double B;              // background magnetic field
    double L;              // external luminosity  
    double z;              // redshift
    double tau_acc;        // acceleration time scale
    double doppler_factor; // observational doppler factor
    double LHS_BC;         // left boundary condition

    // electron injection parameters
    double (*I)();  // injection function (gamma, Sim)
    double inject_power;    // power law power
    double inject_min;      // injection gamma range min
    double inject_max;      // injection gamma range max
    double inject_break;    // injection gamma for power law break
    double inject_power_2;  // power after the break

    // implicit unique params
    double dt;      // fixed time step
    double init_power;  // power to generate initial distribution
    double rho;     // background density for initial population calc
    double end_tol; // tolerance for change in n to count as equilibrium

    /*** behind the scenes parameters ***/
    // data storage structs
    ElectronParams *Electrons;
    PhotonParams *Photons;
    double *nu_flux;
    double *flux_eps;

    // physical parameters 
    double Q_e0;         // electron injection density
    double S;            // sync vs inv. compton proportionality factor
    double tau_esc_free; // free escape time
    double norm;         // normalization factor for the prob dist
    double avg_gamma;    // average gamma of injected dist
    double V;            // volume of system (spherical based on R)
    double gamma_eq;     // point of particle flow equilibrium
    double dL;           // luminosity distance
    double nu;           // background density proton ratio

    // time values
    double t;           // time to be stepped
    double end_t;       // hard simulated time limit
    double final_time;  // the last value of t

    char *filename;
    char *filepath;
    double change;
    bool end_sim;    // find the break point for the upwind and downwind solving
    int64_t array_len;      // every data array is this long (excluding ghost cells)
    int64_t gamma_eq_break; // the index of the equilibrium point
    int64_t iter;
    double sim_time;      // time for whole simulation
    double solve_time;   // time for just relaxation of electrons

} SimulationParams;

void fill_gamma_array(SimulationParams *Sim, ElectronParams *Electrons)
{
    // fill gamma array with equal log step data
    for (int64_t i = 0; i <= Sim->array_len + 1; i++)
    {
        Electrons->gamma[i] = pow(10, log10(Sim->min_gamma) + (i - 1) / (double)Sim->samples_per_decade);
    }

    // calculate gamma steps
    Electrons->dln_gamma = log(Electrons->gamma[1]) - log(Electrons->gamma[0]);
    for (int64_t i = 0; i <= Sim->array_len; i++)
    {
        Electrons->dgamma_fwd[i] = Electrons->gamma[i+1] - Electrons->gamma[i];
    }

}

void fill_eps_array(SimulationParams *Sim, PhotonParams *Photons)
{
    // calculate number of decades in the gamma range
    int64_t decades = (int64_t)log10(Sim->max_eps) - log10(Sim->min_eps);
    int64_t samples_per_decade = Sim->array_len / decades;

    // fill gamma array with equal log step data
    for (int64_t i = 0; i <= Sim->array_len + 1; i++)
    {
        Photons->eps[i] = pow(10, log10(Sim->min_eps) + (i - 1) / (double)samples_per_decade);
    }
}


void malloc_Sim_arrays(SimulationParams *Sim)
{
    Sim->Electrons = malloc(sizeof(ElectronParams));
    Sim->Photons = malloc(sizeof(PhotonParams));
    ElectronParams *Electrons = Sim->Electrons;
    PhotonParams *Photons = Sim->Photons;
    
    // calculate number of decades in the gamma range
    double decades = log10(Sim->max_gamma) - log10(Sim->min_gamma);
    // set array length based on decades and input samples per decade
    Sim->array_len = (int64_t)decades * Sim->samples_per_decade;

    // malloc gamma and delta gamma array
    Electrons->gamma = malloc((Sim->array_len + 2) * sizeof(double));
    Electrons->dgamma_fwd = malloc((Sim->array_len + 2) * sizeof(double));
    Electrons->current_n = malloc((Sim->array_len + 2) * sizeof(double));
    Electrons->next_n = malloc((Sim->array_len + 2) * sizeof(double));

    // malloc gamma and delta gamma arrays
    Photons->eps = calloc(Sim->array_len + 2, sizeof(double));
    Photons->n = calloc((Sim->array_len + 2), sizeof(double));
    Photons->emission = calloc(Sim->array_len + 2, sizeof(double));
    Photons->absorption = calloc(Sim->array_len + 2, sizeof(double));
    Photons->xCS = calloc(Sim->array_len + 2, sizeof(double));
    
    // malloc flux calculation arrays
    Sim->nu_flux = calloc(Sim->array_len + 2, sizeof(double));
    Sim->flux_eps = calloc(Sim->array_len + 2, sizeof(double));
    
    // fill in the energy values of the arrays
    fill_gamma_array(Sim, Sim->Electrons);
    fill_eps_array(Sim, Sim->Photons);
}

void free_Sim_arrays(SimulationParams *Sim)
{
    free(Sim->Electrons->current_n);
    free(Sim->Electrons->next_n);
    free(Sim->Electrons->gamma);
    free(Sim->Electrons);
    free(Sim->Photons->n);
    free(Sim->Photons->eps);
    free(Sim->Photons->emission);
    free(Sim->Photons->absorption);
    free(Sim->Photons->xCS);
    free(Sim->Photons);
    free(Sim->nu_flux);
}

int find_closest_index(double *array, int64_t array_len, double target)
{
    int64_t closest_index = 0;
    double minDiff = fabs(array[0] - target);
    for (int i = 1; i < array_len; i++)
    {
        double diff = fabs(array[i] - target);
        if (diff < minDiff)
        {
            minDiff = diff;
            closest_index = i;
        }
    }

    return closest_index;
}

double freq_from_eps(double eps)
{
    return (m_e * c * c * eps) / h;
}

double eps_from_freq(double freq)
{
    return (h * freq) / (m_e * c * c);
}

double power_law(double gamma, SimulationParams *Sim)
{
    // inject power law within a specified range
    if (gamma < Sim->inject_min || gamma > Sim->inject_max)
    {
        return 0.;
    }
    else
    {
        return Sim->Q_e0 * Sim->norm * pow(gamma, (-1. * Sim->init_power));
    }
}

double broken_power_law(double gamma, SimulationParams *Sim)
{
        // inject power law within a specified range
    if (gamma < Sim->inject_min || gamma > Sim->inject_max)
    {
        return 0.;
    }
    else
    {
      if(gamma > Sim->inject_break)
      {
        return Sim->Q_e0 * Sim->norm 
        * pow(Sim->inject_break / Sim->inject_min, -Sim->inject_power) 
        * pow(gamma / Sim->inject_break, -Sim->inject_power_2);
      }
      else
      {
        return Sim->Q_e0 * Sim->norm * pow(gamma / Sim->inject_min, -Sim->inject_power);
      }
    }
}

void normalize_inject_dist(SimulationParams *Sim)
{
    // normalize power law dist based on a given power
    double integral = 0;
    Sim->Q_e0 = 1.;
    Sim->norm = 1.;
    for (int64_t i = 0; i <= Sim->array_len; i++)
    {
        integral += 0.5 * (Sim->I(Sim->Electrons->gamma[i + 1], Sim) + Sim->I(Sim->Electrons->gamma[i], Sim)) * (Sim->Electrons->dgamma_fwd[i]);
    }
    Sim->norm /= integral;
}

void set_initial_state(SimulationParams *Sim)
{
    // set boundary conditions in ghost cells
    Sim->Electrons->next_n[Sim->array_len+1] = 0.;
    Sim->Electrons->next_n[0] = Sim->LHS_BC;

    for (int64_t i = 1; i <= Sim->array_len; i++)
    {
        // set initial population to be zero
        Sim->Electrons->current_n[i] = 0.;
    }
}

// file writing code
void write_column_to_csv(FILE *file, double *data, int rows, const char *header)
{
    if (!file || !data || !header || rows <= 0)
    {
        fprintf(stderr, "Error: Invalid input parameters\n");
        return;
    }

    // Temporary buffer for reading lines
    char line[4096];
    
    // Temporary file for writing updated content
    FILE *temp = tmpfile();
    if (!temp)
    {
        fprintf(stderr, "Error: Could not create temporary file\n");
        return;
    }

    // Check if original file is empty
    fseek(file, 0, SEEK_END);
    long file_size = ftell(file);
    rewind(file);

    int current_row = 0;
    int is_first_line = 1;

    if (file_size == 0)
    {
        // If file is empty, write header and data directly
        fprintf(temp, "%s\n", header);
        for (int i = 0; i < rows; i++)
        {
            fprintf(temp, "%e\n", data[i]);
        }
    }
    else
    {
        // Read and process existing file
        while (fgets(line, sizeof(line), file))
        {
            // Remove trailing newline
            size_t len = strlen(line);
            if (len > 0 && line[len - 1] == '\n')
            {
                line[len - 1] = '\0';
            }

            if (is_first_line)
            {
                // Add new header to first line
                fprintf(temp, "%s,%s\n", line, header);
                is_first_line = 0;
            }
            else
            {
                if (current_row < rows)
                {
                    // Add data to subsequent lines
                    fprintf(temp, "%s,%e\n", line, data[current_row]);
                    current_row++;
                }
                else
                {
                    // No more data to add
                    fprintf(temp, "%s,\n", line);
                }
            }
        }

        // Add remaining data if any
        while (current_row < rows)
        {
            fprintf(temp, ",%e\n", data[current_row]);
            current_row++;
        }
    }

    // Copy temporary file back to original file
    rewind(temp);
    rewind(file);
    
    // Truncate the original file
    if (ftruncate(fileno(file), 0) != 0)
    {
        fprintf(stderr, "Error: Could not truncate file\n");
        fclose(temp);
        return;
    }

    // Copy contents
    int ch;
    while ((ch = fgetc(temp)) != EOF)
    {
        fputc(ch, file);
    }

    fflush(file);
    fclose(temp);
}

void calc_S(SimulationParams *Sim)
{
    // calculate synchrotron ratio based on physics
    Sim->S = -(4. / 3.) * c * sigma_t * (Sim->B * Sim->B / (8 * M_PI * m_e * c * c));
}

double calc_tau_esc(SimulationParams *Sim)
{
    // calculate escape time for a spherical plasma based on the radius
    Sim->tau_esc_free = (3. * Sim->R) / (4. * c);

    return Sim->tau_esc_free;
}

void calc_V(SimulationParams *Sim)
{
    Sim->V = (4. / 3.) * M_PI * Sim->R * Sim->R * Sim->R;
}

void calc_Q_e0(SimulationParams *Sim)
{
    double hold_power, hold_g_max, hold_g_min, avg_gamma_protons = 0.;
    double (*hold_I)();
    hold_I = Sim->I;
    Sim->avg_gamma = 0.;
    Sim->Q_e0 = 1.;

    // calculate average gamma
    for (int64_t i = 0; i <= Sim->array_len + 1; i++)
    {
        Sim->avg_gamma += 0.5 * (Sim->Electrons->gamma[i + 1] * Sim->I(Sim->Electrons->gamma[i + 1], Sim) + Sim->Electrons->gamma[i] * Sim->I(Sim->Electrons->gamma[i], Sim)) * (Sim->Electrons->dgamma_fwd[i]);
    }

    // set up temp proton distribution based on KATU specifications
    hold_power = Sim->inject_power;
    hold_g_max = Sim->inject_max;
    hold_g_min = Sim->inject_min;

    Sim->inject_power = 2.;
    Sim->inject_min = 1e1;
    Sim->inject_max = 1e6;
    Sim->I = &power_law;
    // use these specified power law parameters to renormalise
    normalize_inject_dist(Sim);

    // calculate average proton energy based
    for (int64_t i = 0; i <= Sim->array_len + 1; i++)
    {
        avg_gamma_protons += 0.5 * (Sim->Electrons->gamma[i + 1] * power_law(Sim->Electrons->gamma[i + 1], Sim) + Sim->Electrons->gamma[i] * power_law(Sim->Electrons->gamma[i], Sim)) * (Sim->Electrons->dgamma_fwd[i]);
    }
    
    // reset injection distributions
    Sim->inject_power = hold_power;
    Sim->inject_max = hold_g_max;
    Sim->inject_min = hold_g_min;
    Sim->I = hold_I;
    normalize_inject_dist(Sim);

    // calculate injected Qe0
    Sim->Q_e0 = (Sim->L / Sim->V) / ((Sim->nu * avg_gamma_protons * m_p + Sim->avg_gamma * m_e)* c * c);
}

void implicit_step(SimulationParams *Sim, ElectronParams *Electrons)
{
    
    double tau_esc = Electrons->tau_esc;
    double dln_gamma = Electrons->dln_gamma;
    double *gamma = Electrons->gamma;
    double *next_n = Electrons->next_n;
    double *current_n = Electrons->current_n;

    for (int64_t i = 1; i <= Sim->gamma_eq_break; i++)
    {
        // implicit stepping regime
        next_n[i] =
            tau_esc * 
            (Sim->S * Sim->dt * pow(gamma[i-1], 2.) * next_n[i-1] * Sim->tau_acc
            + Sim->dt * dln_gamma * gamma[i] * Sim->tau_acc * Sim->I(gamma[i], Sim)
            + dln_gamma * gamma[i] * current_n[i] * Sim->tau_acc
            + Sim->dt * gamma[i-1] * next_n[i-1])
            /
            (gamma[i] * 
            (Sim->S * tau_esc * Sim->tau_acc * Sim->dt * gamma[i]
            + Sim->dt * dln_gamma * Sim->tau_acc
            + dln_gamma * Sim->tau_acc * tau_esc
            + tau_esc * Sim->dt));
    }
    for (int64_t i = Sim->array_len; i >= Sim->gamma_eq_break+1; i--)
    {            
        next_n[i] =
            tau_esc * 
            (Sim->S * Sim->dt * pow(gamma[i+1], 2.) * next_n[i+1] * Sim->tau_acc
            - Sim->dt * dln_gamma * gamma[i] * Sim->tau_acc* Sim->I(gamma[i], Sim)
            - dln_gamma * gamma[i] * current_n[i] * Sim->tau_acc
            + Sim->dt * gamma[i+1] * next_n[i+1])
            /
            (gamma[i] * 
            (Sim->S * tau_esc * Sim->tau_acc * Sim->dt * gamma[i]
            - Sim->dt * dln_gamma * Sim->tau_acc
            - dln_gamma * Sim->tau_acc * tau_esc
            + tau_esc * Sim->dt)); 
    }
}

void save_step_to_prev_n(SimulationParams *Sim, ElectronParams *Electrons)
{
    for (int64_t i = Sim->array_len + 1; i >= 0; i--)
    {
        // implicit stepping regime
        Electrons->current_n[i] = Electrons->next_n[i];
    }
}

bool equilibrium_check(SimulationParams *Sim, ElectronParams *Electrons)
{
    Sim->change = 0.;
    // calculate percentage change in n
    for (int64_t i = Sim->array_len; i >= 1; i--)
    {
        if (Electrons->next_n[i] < 1e-256 || Electrons->next_n[i+1] <  1e-256)
        {
            continue;
        }
        
        Sim->change += 0.5 * 
        (pow((1. / Electrons->next_n[i]) * ((Electrons->next_n[i] - Electrons->current_n[i])), 2.)
        + pow((1. / Electrons->next_n[i+1]) * ((Electrons->next_n[i+1] - Electrons->current_n[i+1])), 2.))
        * (Sim->Electrons->dgamma_fwd[i]);
    }
    Sim->change /= pow(Sim->dt, 2.);
    Sim->change = sqrt(Sim->change);
    // check change in population against specified end tolerance
    if (Sim->change < Sim->end_tol)
    {
        printf("equilibrium reached at t = %e, iter = %lld, last change %e\n", Sim->t, Sim->iter, Sim->change);
        Sim->end_sim=true;
    }
}

// critical frequency for photon calculation
double nu_crit(double gamma, SimulationParams *Sim)
{
    return (gamma * gamma * Sim->B * 3. * q) / (4. * M_PI * m_e * c);
}

#if APPROXIMATE_XCS == 0
// Function to compute the Whittaker function W
double W(double kappa, double mu, double z)
{
    // calculate a and b
    double a = 0.5 + mu - kappa;
    double b = 1.0 + 2.0 * mu;

    // Compute the confluent hypergeometric function
    double U_val = gsl_sf_hyperg_U(a, b, z);

    double scaling_factor = exp(-0.5 * z) * pow(z, 0.5 + mu);

    return scaling_factor * U_val;
}

double CS(double x)
{
    return W(0., 4. / 3., x) * W(0., 1. / 3., x) - W(1. / 2., 5. / 6., x) * W(-1. / 2., 5. / 6., x);
}

double fill_xCS_array(double freq, SimulationParams *Sim)
{
    for (int64_t i = 1; i <= Sim->array_len - 1; i++)
    {
        double *gamma = Sim->Electrons->gamma;
        double *xCS = Sim->Photons->xCS;
        // use trapezium rule
        double x = (freq / (nu_crit(gamma[i], Sim)));
        xCS[i] = x * CS(x);
    }
}
#else
void fill_xCS_array(double freq, SimulationParams *Sim)
{
    // use approximate 
    double *xCS = Sim->Photons->xCS;;
    for (int64_t i = 1; i <= Sim->array_len - 1; i++)
    {
    double x = (freq / (nu_crit(Sim->Electrons->gamma[i], Sim)));
    if(x < 0.291)
        xCS[i] = (29 * cbrt(x) - 27 * x) / 25;
    else if(x < 2.7675)
        xCS[i] = (135 - 77 * x + 12 * x * x) / 250;
    else
        xCS[i] = exp(-x) * (1 - 1 / (3 * x));
    }
}
#endif

double calc_absorption(double freq, SimulationParams *Sim)
{
    double integral = 0.0;
    double temp_val;
    double *n = Sim->Electrons->next_n;
    double *xCS = Sim->Photons->xCS;

    // following the optimized KATU calculation 
    for (int64_t i = 1; i <= Sim->array_len - 1; i++)
    {
        if (n[i+1] < 1e-256 || n[i] < 1e-256)
        {
            continue;
        }
        temp_val = (log(n[i+1]) - log(n[i])) / Sim->Electrons->dln_gamma;
        temp_val = n[i] * (temp_val - 2.);
        // use trapezium rule
        integral += 2. * xCS[i] * temp_val;
    }
    
    return integral * Sim->Electrons->dln_gamma / 2.;
}

double calc_emission(double freq, SimulationParams *Sim)
{
    double integral = 0.0;
    // following the optimized KATU calculation 
    for (int64_t i = 1; i <= Sim->array_len - 1; i++)
    {
        double *gamma = Sim->Electrons->gamma;
        double *n = Sim->Electrons->next_n;
        double *xCS = Sim->Photons->xCS;
        // use trapezium rule
        integral += 2. * gamma[i] * n[i] * xCS[i];
    }

    return integral * Sim->Electrons->dln_gamma / 2.;
}

void photon_calc(SimulationParams *Sim)
{
    // use KATU's optimized photon calculation
    double photon_gains_factor, freq_Hz, photon_losses_factor;
    double *emission = Sim->Photons->emission;
    double *absorption = Sim->Photons->absorption;

    // step through all frequencies of photon population and find the emission
    // or absorption from the electron population
    for (int64_t i = 1; i <= Sim->array_len; i++)
    {
        photon_gains_factor = sqrt(3.) * M_PI /  3. * FINE_STRUCTURE_CONSTANT * nu_crit(1., Sim);
        photon_losses_factor = sqrt(3.) * M_PI / 12. * ELECTRON_ENERGY * ELECTRON_RADIUS * nu_crit(1., Sim) / m_e;
        
        // find frequency for given normalized energy bin
        freq_Hz = freq_from_eps(Sim->Photons->eps[i]);
        // precalculate xCS to save computational time
        fill_xCS_array(freq_Hz, Sim);

        emission[i] = calc_emission(freq_Hz, Sim);
        //absorption[i] = calc_absorption(freq_Hz, Sim);

        Sim->Photons->n[i] = photon_gains_factor * emission[i] / Sim->Photons->eps[i];
        Sim->Photons->n[i] *= Sim->Photons->tau_esc;
        // unsure why this didnt work
        //printf("%.3e, %.3e, %.3e\n", freq_Hz, Sim->Photons->n[i],photon_losses_factor * absorption[i] / (freq_Hz * freq_Hz) * Sim->Photons->n[i]);
        //Sim->Photons->n[i] += Sim->Photons->tau_esc * photon_losses_factor * absorption[i] / (freq_Hz * freq_Hz) * Sim->Photons->n[i];
    }
}

void calc_flux(SimulationParams *Sim)
{
    // Calculate the luminosity distance from z<<1 approx
    Sim->dL = c * Sim->z / H0;

    double pre_factor, eps_obs;
    pre_factor = m_e * pow(c, 2.) * pow(Sim->doppler_factor, 4.) * Sim->V;
    pre_factor /= 4. * M_PI * (1. + Sim->z) * pow(Sim->dL, 2.) * Sim->Photons->tau_esc;

    // Iterate over the array of intrinsic frequencies
    for (int64_t i = 1; i <= Sim->array_len; i++)
    {
        // Calculate observed frequency for each intrinsic frequency
        eps_obs = Sim->doppler_factor * Sim->Photons->eps[i] / (1. + Sim->z);
        // Calculate the flux in the observer frame
        Sim->nu_flux[i] =
            pre_factor * (pow(Sim->Photons->eps[i], 2.) * Sim->Photons->n[i]);

        Sim->flux_eps[i] = eps_obs;
    }
}

void write_run_file(char *filename, SimulationParams *Sim)
{
    char filepath[100];
    sprintf(filepath, "csv_data/implicit_solve/runs/run_%s", filename);

    FILE *run_file = fopen(filepath, "w");
    
    fprintf(run_file, 
    "dt,R,inject_p,inject_min,inject_max,rho,B,L,end_tol,Q_e0,S,tau_esc,norm,avg_gamma,V,array_len,max_gamma,min_gamma,init_p,samples_per_decade,final_time,change,tau_acc,sim_time,solve_time\n");
    fprintf(run_file,
    "%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%lld,%e,%e,%e,%lld,%e,%e,%e,%e,%e\n",
    Sim->dt,Sim->R,Sim->inject_power,Sim->inject_min,Sim->inject_max,Sim->rho,Sim->B, 
    Sim->L,Sim->end_tol, Sim->Q_e0,Sim->S,Sim->tau_esc_free,Sim->norm,Sim->avg_gamma,Sim->V,
    Sim->array_len,Sim->max_gamma,Sim->min_gamma,Sim->init_power,Sim->samples_per_decade,
    Sim->final_time,Sim->change, Sim->tau_acc,Sim->sim_time,Sim->solve_time);
    fclose(run_file);
}

void simulate(char *filename, SimulationParams *Sim)
{
    struct timeval start1, end1, start2, end2;
    // start simulation timer
    gettimeofday(&start1, NULL);

    // allocate simulation memory
    malloc_Sim_arrays(Sim);
    char header[100];

    // write gammas to csv
    // create file
    Sim->filename = filename;
    Sim->filepath = malloc(sizeof(char) * 100);
    sprintf(Sim->filepath, "csv_data/implicit_solve/%s", filename);
    FILE *file = fopen(Sim->filepath, "w");
    // write gammas to csv
    write_column_to_csv(file, Sim->Electrons->gamma, Sim->array_len + 2, "gamma");
    fclose(file);
    
    // calculate physical parameters
    calc_S(Sim);
    normalize_inject_dist(Sim);
    calc_V(Sim);
    calc_Q_e0(Sim);
    calc_tau_esc(Sim);
    Sim->Electrons->tau_esc = Sim->tau_esc_free * CFE_RATIO;
    Sim->Photons->tau_esc = Sim->tau_esc_free;
    
    // find the break point for the upwind and downwind solving
    Sim->gamma_eq = -1. / (Sim->tau_acc * Sim->S);
    Sim->gamma_eq_break = find_closest_index(Sim->Electrons->gamma, Sim->array_len, Sim->gamma_eq);
    
    // ensure dt is not larger than tau_esc for stability
    if (Sim->dt > Sim->tau_esc_free)
    {
        Sim->dt = Sim->tau_esc_free;
        printf("dt larger than tau_esc. dt set to tau_esc to ensure stability\n");
    }

    // setup for simulation
    set_initial_state(Sim);
    Sim->t = 0.;
    Sim->iter = 0;
    Sim->end_sim = false;

    // start simulation
    printf("samples:%lld\n", Sim->array_len);
    gettimeofday(&start2, NULL);
    while (Sim->t < Sim->end_t && Sim->end_sim == false)
    {
        implicit_step(Sim, Sim->Electrons);
        equilibrium_check(Sim, Sim->Electrons);
        save_step_to_prev_n(Sim, Sim->Electrons);
        Sim->t += Sim->dt;
        Sim->iter ++;
        
        //if ((Sim->t > 1e6 && Sim->iter == 100) || (Sim->t > 1e5 && Sim->iter == 10 && Sim->t < 1e6) || Sim->t < 1e5)
        //{
        //    sprintf(header, "n_e t=%e", Sim->t);
        //    write_column_to_csv(filepath, Sim->Electrons->current_n, Sim->array_len+2, header);
        //    Sim->iter = 0;
        //}
        
    }
    gettimeofday(&end2, NULL);
    Sim->solve_time = (end2.tv_sec - start2.tv_sec) + (end2.tv_usec - start2.tv_usec) / 1000000.0;
    
    Sim->final_time=Sim->t;
    // generate photon population
    photon_calc(Sim);
    // generate the flux array based on photons
    calc_flux(Sim);

    FILE *file2 = fopen(Sim->filepath, "r+");
    write_column_to_csv(file2, Sim->Electrons->next_n, Sim->array_len + 2, "electron_n");
    write_column_to_csv(file2, Sim->Photons->eps, Sim->array_len + 2, "photon_eps");
    write_column_to_csv(file2, Sim->Photons->n, Sim->array_len + 2, "photon_n");
    write_column_to_csv(file2, Sim->Photons->emission, Sim->array_len + 2, "emission");
    write_column_to_csv(file2, Sim->Photons->absorption, Sim->array_len + 2, "absorption");
    write_column_to_csv(file2, Sim->flux_eps, Sim->array_len + 2, "flux_eps");
    write_column_to_csv(file2, Sim->nu_flux, Sim->array_len + 2, "nu_flux");
    fclose(file2);

    gettimeofday(&end1, NULL);
    Sim->sim_time = (end1.tv_sec - start1.tv_sec) + (end1.tv_usec - start1.tv_usec) / 1000000.0;
    printf("solving time: %e\n", Sim->solve_time);
    printf("sim time: %e\n", Sim->sim_time);
    write_run_file(Sim->filename, Sim);
}

int main()
{
    SimulationParams *Sim = malloc(sizeof(SimulationParams));
    
    // Cooling test
    Sim->nu = 1e-10;
    Sim->inject_min = 1e4;
    Sim->inject_max = 1e8;
    Sim->inject_power = 2.3;
    Sim->B = 0.1;
    Sim->R = 1e16;
    Sim->L = 1e30;
    Sim->doppler_factor = pow(10, 1.83);
    Sim->z = 0.33;
    Sim->tau_acc = 1e100;
    Sim->I = &power_law;
    Sim->LHS_BC = 0;

    Sim->nu = 1;
    Sim->inject_min = pow(10, 1.33);
    Sim->inject_break = pow(10, 4.03);
    Sim->inject_max = pow(10, 6.82);
    Sim->inject_power = 1.69;
    Sim->inject_power_2 = 4.29;
    Sim->B = pow(10, -1.01);
    Sim->R = pow(10, 16.45);
    //Sim->L = pow(10, 45.647);
    Sim->L = pow(10, 45.37);
    Sim->doppler_factor = pow(10, 1.44);
    Sim->tau_acc = 1e256;
    Sim->z = 0.34;
    Sim->I = &broken_power_law;
    Sim->LHS_BC = 0.;

    // simulation setup  
    Sim->min_gamma = 1e1;
    Sim->max_gamma = 1e8;
    Sim->min_eps = 1e-12;
    Sim->max_eps = 1e8;
    Sim->samples_per_decade = 80;
    Sim->init_power = 2.;
    Sim->dt = 1e100;    // max out the time step
    Sim->end_t = 1e10;
    Sim->end_tol = 1e-8;
    
    /*
    // Cooling test
    Sim->nu = 1e-10;
    Sim->inject_min = 1e4;
    Sim->inject_max = 1e8;
    Sim->inject_power = 2.3;
    Sim->B = 0.01;
    Sim->R = 1e16;
    Sim->L = 1e30;
    Sim->doppler_factor = pow(10, 1.83);
    Sim->z = 0.33;
    Sim->tau_acc = 1e100;
    Sim->I = &power_law;
    Sim->LHS_BC = 0;

    // generate the cooling test data
    hold = Sim->B;
    double B[6] = {0.1,0.25,0.5,1.,1.5,2.};
    for (int i =0; i < 6; i++)
    {
        Sim->B = B[i];
        // generate file name based on B
        char filename[100];
        sprintf(filename, "B%.2lf.csv", Sim->B);
        printf("\nSimulate: %s\n", filename);
        simulate(filename, Sim);
    }
    Sim->B = hold;
    */

    /*
    //Acceleration test
    Sim->nu = 1e-10;
    Sim->inject_min = 1e1;
    Sim->inject_max = 1e2;
    Sim->inject_power = 2.3;
    Sim->B = 0.0001;
    Sim->R = 1e16;
    Sim->L = 1e30;
    Sim->doppler_factor = pow(10, 1.83);
    Sim->z = 0.33;
    Sim->tau_acc = calc_tau_esc(Sim) * .5;
    Sim->I = &power_law;
    Sim->LHS_BC = 5e-10;

    // generate the acceleration test data
    hold = Sim->tau_acc;
    Sim->tau_esc_free = calc_tau_esc(Sim);
    double t_acc_multi[5] = {0.5,1.0,1.5,2.0,4.0};
    for (int i =0; i < 5; i++)
    {
        Sim->tau_acc = t_acc_multi[i] * Sim->tau_esc_free;
        // generate file name based on B
        char filename[100];
        sprintf(filename, "t_acc=%1.1lf_t_esc.csv", t_acc_multi[i]);
        printf("\nSimulate: %s\n",filename);
        simulate(filename, Sim);
    }
    Sim->tau_acc = hold;
    */

    simulate("simulation_data.csv", Sim);

    /* 
    hold = (double)Sim->samples_per_decade;
    double param[5] = {5,10,20,30,40};
    for (int i =0; i < 5; i++)
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
    Sim->samples_per_decade = (int64_t)hold;
    */

    /*
    hold = Sim->dt;
    double param[5] = {1e4,2.5e4,5e4,1e5,2.5e5};
    for (int i =0; i < 5; i++)
    {
        Sim->dt = param[i];
        malloc_and_fill_gamma_array(Sim, Sim->Species[0]);
        // generate file name based on B
        char filename[150], filepath[150], run_filepath[150];
        sprintf(filename, "dt%.0lf.csv", Sim->dt);
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
    Sim->dt = hold;
    */

    // end program
    free_Sim_arrays(Sim);
    free(Sim);
    return 0;
}