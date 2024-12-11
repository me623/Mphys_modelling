
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <stdbool.h>
#include <sys/time.h>
#include <gsl/gsl_sf_hyperg.h>

#define BUFFER_SIZE 33554432
// CGS unit constants
#define m_e 9.109e-28
#define c 2.998e10
#define sigma_t 6.652e-25
#define q 4.803e-10
#define h 6.626e-27
#define H0 2.27e-18 // hubble constant in CGS
#define h_bar 1.055e-27

typedef struct ElectronParams
{
    double *current_n;
    double *next_n;
    double *gamma;
    double *dgamma_fwd;
    double dln_gamma;
} ElectronParams;

typedef struct PhotonParams
{
    double *n;
    double *eps;
    double *delta_freq;
    double *j_nu;
} PhotonParams;

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
    double inject_power_2;
    double init_power;  // power to generate initial distribution
    double inject_min;  // injection gamma range min
    double inject_break;
    double inject_max;  // injection gamma range max
    double max_freq;
    double min_freq;
    double rho;     // background density for initial population calc
    double B;       // background magnetic field
    double L;       // external luminosity
    double tau_acc; // acceleration timescale
    double end_tol; // tolerance for change in n to count as equilibrium
    double doppler_factor;
    double z;

    double Q_e0;    // starting population
    double S;       // sync vs inv. compton
    double tau_esc; // free escape time
    double norm;    // normalize the prob dist
    double avg_gamma;   // average gamma of injected dist for Q_e0 calc
    double V;           // volume of system (spherical based on R)   
    double LHS_BC;
    double dL;
    double final_time;

    // code specific values
    ElectronParams *Electrons;
    PhotonParams *Photons;
    double (*I)();
    double *j_nu;
    double *nu_flux;
    double *flux_eps;
    double change;
    bool end_sim;    // find the break point for the upwind and downwind solving
    double gamma_eq;
    int64_t gamma_eq_break;
    int64_t iter;
    double sim_time;
    double solve_time;

} SimulationParams;

void malloc_and_fill_gamma_array(SimulationParams *Sim, ElectronParams *Electrons)
{
    int64_t decades;
    // calculate number of decades in the gamma range
    decades = (int64_t)log10(Sim->max_gamma) - log10(Sim->min_gamma);
    // set array length based on decades and input samples per decade
    Sim->array_len = decades * Sim->samples_per_decade;

    // malloc gamma and delta gamma array
    Electrons->gamma = malloc((Sim->array_len + 2) * sizeof(double));
    Electrons->dgamma_fwd = malloc((Sim->array_len + 2) * sizeof(double));
    // fill gamma array with equal log step data
    for (int64_t i = 1; i <= Sim->array_len + 1; i++)
    {
        Electrons->gamma[i] = pow(10, log10(Sim->min_gamma) + (i-1) / (double)Sim->samples_per_decade);
        if (i > 0)
        {
            Electrons->dgamma_fwd[i-1] = Electrons->gamma[i] - Electrons->gamma[i-1];
        }
    }
    // Calculate the step size in logs
    double log_step = 1.0 / (double)Sim->samples_per_decade;
    // Extrapolate gamma[0] based on gamma[1]
    Electrons->gamma[0] = pow(10, log10(Sim->min_gamma) - log_step);

    Electrons->dln_gamma = log(Electrons->gamma[1]) - log(Electrons->gamma[0]);
}
void malloc_and_fill_frequency_array(SimulationParams *Sim, PhotonParams *Photons)
{
    int64_t decades, samples_per_decade;
    // calculate number of decades in the gamma range
    decades = (int64_t)log10(Sim->max_freq) - log10(Sim->min_freq);
    samples_per_decade = Sim->array_len / decades;
    
    // malloc gamma and delta gamma arrays
    Photons->eps = calloc(Sim->array_len + 2, sizeof(double));
    Photons->delta_freq  = calloc(Sim->array_len + 2, sizeof(double));

    // fill gamma array with equal log step data
    for (int64_t i = 1; i <= Sim->array_len + 1; i++)
    {
        Photons->eps[i] = pow(10, log10(Sim->min_freq) + (i-1) / (double)samples_per_decade);
    }

    // Calculate the step size in logs
    double log_step = 1.0 / (double)Sim->samples_per_decade;
    // Extrapolate gamma[0] based on gamma[1]
    Photons->eps[0] = pow(10, log10(Sim->min_freq) - log_step);

   for (int64_t i = 1; i < Sim->array_len; i++)
    {
        Photons->delta_freq[i] = Photons->eps[i + 1] - Photons->eps[i];
    }
}

void malloc_Sim_arrays(SimulationParams *Sim)
{
    Sim->Electrons = malloc(sizeof(ElectronParams));
    malloc_and_fill_gamma_array(Sim, Sim->Electrons);
    Sim->Electrons->current_n = malloc((Sim->array_len + 2) * sizeof(double));
    Sim->Electrons->next_n = malloc((Sim->array_len + 2) * sizeof(double));

    Sim->Photons = malloc(sizeof(PhotonParams));
    malloc_and_fill_frequency_array(Sim, Sim->Photons);
    Sim->Photons->n = calloc((Sim->array_len + 2), sizeof(double));
    Sim->Photons->j_nu = calloc((Sim->array_len + 2), sizeof(double));
    Sim->nu_flux = calloc((Sim->array_len + 2), sizeof(double));
    Sim->flux_eps = calloc((Sim->array_len + 2), sizeof(double));
}

void free_Sim_arrays(SimulationParams *Sim)
{
    free(Sim->Electrons->current_n);
    free(Sim->Electrons->next_n);
    free(Sim->Electrons->gamma);
    free(Sim->Electrons->dgamma_fwd);
    free(Sim->Electrons);
    free(Sim->Photons->n);
    free(Sim->Photons->eps);
    free(Sim->Photons->delta_freq);
    free(Sim->Photons->j_nu);
    free(Sim->Photons);
    free(Sim->nu_flux);
    free(Sim->flux_eps);
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

void normalize_inject_dist(double power, SimulationParams *Sim)
{
        // normalize power law dist based on a given power
       double integral = 0;
       Sim->Q_e0 = 1.;
       Sim->norm = 1.;
       for (int64_t i = 0; i <= Sim->array_len + 2; i++)
       {
        integral += 0.5 
        * (Sim->I(Sim->Electrons->gamma[i + 1], Sim) 
            + Sim->I(Sim->Electrons->gamma[i], Sim))
        * (Sim->Electrons->gamma[i+1] - Sim->Electrons->gamma[i]);
       }
       Sim->norm /= integral;
       
}

void set_initial_state(SimulationParams *Sim)
{
    Sim->Electrons->next_n[Sim->array_len+1] = 0.;
    Sim->Electrons->next_n[0] = Sim->LHS_BC;
    for (int64_t i = 1; i <= Sim->array_len; i++)
    {
        // set initial population on a selected power law
        // number of photons included based on background density
        Sim->Electrons->current_n[i] =
        (Sim->rho / (2. * m_e)) * 
        (Sim->I(Sim->Electrons->gamma[i], Sim) 
        / Sim->Q_e0);
        Sim->Electrons->current_n[i] = 0.;
    }
}

// file writing code
void write_column_to_csv(const char *filename, double *data, int rows, const char *header) {
    // Try to open the file to check if it exists and has content
    FILE *original = fopen(filename, "r");
    char temp_filename[256];
    snprintf(temp_filename, sizeof(temp_filename), "%s.tmp", filename);
    FILE *temp = fopen(temp_filename, "w");
    
    if (!temp) {
        printf("Error: Cannot create temporary file\n");
        if (original) fclose(original);
        return;
    }

    // Check if file exists and has content
    int is_empty = 1;
    if (original) {
        // Check if file has any content
        fseek(original, 0, SEEK_END);
        is_empty = (ftell(original) == 0);
        rewind(original);
    }

    // If file doesn't exist or is empty, create new file with header and data
    if (!original || is_empty) {
        fprintf(temp, "%s\n", header);  // Write header first
        for (int i = 0; i < rows; i++) {
            fprintf(temp, "%e\n", data[i]);
        }
    } else {
        char line[65536];
        int current_row = 0;
        int is_first_line = 1;

        // Read each line from the original file
        while (fgets(line, sizeof(line), original)) {
            // Remove newline if present
            size_t len = strlen(line);
            if (len > 0 && line[len-1] == '\n') {
                line[len-1] = '\0';
                len--;
            }

            if (is_first_line) {
                // Add header to existing headers
                fprintf(temp, "%s,%s\n", line, header);
                is_first_line = 0;
            } else {
                // Add data
                if (current_row < rows) {
                    fprintf(temp, "%s,%e\n", line, data[current_row]);
                    current_row++;
                } else {
                    fprintf(temp, "%s,\n", line);
                }
            }
        }

        // If there are more data rows than existing CSV rows
        while (current_row < rows) {
            // Add commas for empty columns in original file
            int comma_count = 0;
            char *ptr = line;
            while (*ptr) {
                if (*ptr == ',') comma_count++;
                ptr++;
            }
            
            for (int i = 0; i < comma_count; i++) {
                fprintf(temp, ",");
            }
            fprintf(temp, "%e\n", data[current_row]);
            current_row++;
        }
    }

    if (original) fclose(original);
    fclose(temp);

    // Replace original file with temporary file
    remove(filename);
    rename(temp_filename, filename);
}

void calc_S(SimulationParams *Sim)
{
    // calculate synchrotron ratio based on physics
    Sim->S = -(4. / 3.) * c * sigma_t * (Sim->B * Sim->B / (8 * M_PI * m_e * c * c));
}

void calc_Q_e0(SimulationParams *Sim)
{
    // calculate a pre-solved    integral of gamma * gamma^-p
    Sim->avg_gamma = 0.;
    Sim->Q_e0 = 1.;
    for (int64_t i = 0; i <= Sim->array_len + 1; i++)
    {
        Sim->avg_gamma += 0.5 
        * (Sim->Electrons->gamma[i+1] * Sim->I(Sim->Electrons->gamma[i+1], Sim)
            + Sim->Electrons->gamma[i] * Sim->I(Sim->Electrons->gamma[i], Sim)) 
        * (Sim->Electrons->gamma[i+1] - Sim->Electrons->gamma[i]);
    }
    // calculate volume
    Sim->V = (4./3.) * M_PI * Sim->R * Sim->R * Sim->R;
    // calculate injected Qe0
    Sim->Q_e0 = Sim->L / (Sim->V * Sim->avg_gamma * m_e * c * c);
}

double calc_tau_esc(SimulationParams *Sim)
{
    // calculate escape time for a spherical plasma based on the radius
    Sim->tau_esc = (3. / 4.) * (Sim->R / c);
    return Sim->tau_esc;
}

void implicit_step(SimulationParams *Sim, ElectronParams *Electron)
{
    
    for (int64_t i = 1; i <= Sim->gamma_eq_break; i++)
    {
        // implicit stepping regime
        Electron->next_n[i] =
            Sim->tau_esc * 
            (Sim->S * Sim->dt * pow(Electron->gamma[i-1], 2.) * Electron->next_n[i-1] * Sim->tau_acc
            + Sim->dt * Electron->dln_gamma * Electron->gamma[i] * Sim->tau_acc * Sim->I(Electron->gamma[i], Sim)
            + Electron->dln_gamma * Electron->gamma[i] * Electron->current_n[i] * Sim->tau_acc
            + Sim->dt * Electron->gamma[i-1] * Electron->next_n[i-1])
            /
            (Electron->gamma[i] * 
            (Sim->S * Sim->tau_esc * Sim->tau_acc * Sim->dt * Electron->gamma[i]
            + Sim->dt * Electron->dln_gamma * Sim->tau_acc
            + Electron->dln_gamma * Sim->tau_acc * Sim->tau_esc
            + Sim->tau_esc * Sim->dt));
    }
    for (int64_t i = Sim->array_len; i >= Sim->gamma_eq_break+1; i--)
    {            
        Electron->next_n[i] =
            Sim->tau_esc * 
            (Sim->S * Sim->dt * pow(Electron->gamma[i+1], 2.) * Electron->next_n[i+1] * Sim->tau_acc
            - Sim->dt * Electron->dln_gamma * Electron->gamma[i] * Sim->tau_acc* Sim->I(Electron->gamma[i], Sim)
            - Electron->dln_gamma * Electron->gamma[i] * Electron->current_n[i] * Sim->tau_acc
            + Sim->dt * Electron->gamma[i+1] * Electron->next_n[i+1])
            /
            (Electron->gamma[i] * 
            (Sim->S * Sim->tau_esc * Sim->tau_acc * Sim->dt * Electron->gamma[i]
            - Sim->dt * Electron->dln_gamma * Sim->tau_acc
            - Electron->dln_gamma * Sim->tau_acc * Sim->tau_esc
            + Sim->tau_esc * Sim->dt)); 
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
        if (Electrons->next_n[i] == 0. || Electrons->next_n[i+1] == 0.)
        {
            continue;
        }
        
        Sim->change += 0.5 * 
        (pow((1. / Electrons->next_n[i]) * ((Electrons->next_n[i] - Electrons->current_n[i])), 2.)
        + pow((1. / Electrons->next_n[i+1]) * ((Electrons->next_n[i+1] - Electrons->current_n[i+1])), 2.))
        * (Sim->Electrons->dgamma_fwd[i]);
    }
    Sim->change /= (Sim->dt * Sim->dt);
    Sim->change = sqrt(Sim->change);
    // check change in population against specified end tolerance
    if (Sim->change < Sim->end_tol)
    {
        printf("equilibrium reached at t = %e, iter = %lld, last change %e\n", Sim->t, Sim->iter, Sim->change);
        Sim->end_sim=true;
        return true;
    }
    else
    {
        return false;
    }
}

// Function to compute the Whittaker function W
double W(double kappa, double mu, double z) {
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
    return W(0., 4./3., x) * W(0., 1./3., x) - W(1./2., 5./6., x) * W(-1./2., 5./6., x);
}   

// critical eps
double nu_crit(double gamma, SimulationParams *Sim)
{
    return (gamma * gamma * Sim->B * 3. * q) / (4. * M_PI * m_e * c);
}

// energy emitted by a particle with energy gamma at eps mu
double P_nu(double gamma, double freq, SimulationParams *Sim)
{
    double x = (freq / nu_crit(gamma, Sim));
    return ((sqrt(3.) * pow(q, 3.) * Sim->B) / (4. * M_PI * m_e * pow(c, 2.))) * x * CS(x);
    
}

double j_nu(double freq, SimulationParams *Sim)
{
    double integral = 0.0;

    for (int64_t i = 1; i < Sim->array_len - 1; i++)
    {
        double *gamma = Sim->Electrons->gamma;
        double *n = Sim->Electrons->next_n;
        double *dgamma = Sim->Electrons->dgamma_fwd;
        // use trapezium rule
        integral += 0.5 * (n[i] * P_nu(gamma[i], freq, Sim) + n[i+1] * P_nu(gamma[i+1], freq, Sim)) * dgamma[i];
    }

    return integral / (4.0 * M_PI);
}

double alpha_nu(double freq, SimulationParams *Sim) {
    double integrand, integrand_next, integral = 0.0;
    double *gamma = Sim->Electrons->gamma;
    double *dgamma = Sim->Electrons->dgamma_fwd;
    double *n = Sim->Electrons->next_n;
    for (int64_t i = 1; i < Sim->array_len; i++) {
        integrand = 
        pow(gamma[i], 2) 
        * P_nu(gamma[i], freq, Sim) 
        * ((n[i+1] / pow(gamma[i+1], 2)) - (n[i] / pow(gamma[i], 2))) / dgamma[i];
        integrand_next = 
        pow(gamma[i+1], 2) 
        * P_nu(gamma[i+1], freq, Sim) 
        * ((n[i+2] / pow(gamma[i+2], 2)) - (n[i+1] / pow(gamma[i+1], 2))) / dgamma[i+1];
        
        integral += 0.5 * (integrand + integrand_next) * dgamma[i];
    }
    
    integral /= (-8.0 * M_PI * m_e * pow(freq, 2));
    return integral;
}

double freq_from_eps(double eps)
{
    return (m_e * c * c * eps) / h;
}

double eps_from_freq(double freq)
{
    // Calculate intrinsic epsilon
    return (h * freq) / (m_e * c * c);
}

void photon_calc(SimulationParams *Sim)
{
    for (int64_t i = 1; i <= Sim->array_len; i++)
    {
        Sim->Photons->j_nu[i] = j_nu(freq_from_eps(Sim->Photons->eps[i]), Sim);
        Sim->Photons->n[i] = (2 * Sim->Photons->j_nu[i]) / (h_bar * Sim->Photons->eps[i]);
        //Sim->Photons->n[i] -= Sim->Photons->n[i] / Sim->tau_esc;
        //Sim->Photons->n[i] += c * alpha_nu(freq_from_eps(Sim->Photons->eps[i]), Sim) * Sim->Photons->n[i];
    }

}

void calc_flux(SimulationParams *Sim)
{
    // Calculate the luminosity distance from z<<1 approx
    Sim->dL = c * Sim->z / H0;
    //printf("dl %.2eGpc\n", Sim->dL / 3.086e27);
    
    double pre_fac, eps_obs;
    pre_fac = m_e * pow(c, 2.) * pow(Sim->doppler_factor, 4.) * Sim->V;
    pre_fac /= 4. * M_PI * (1. + Sim->z) * pow(Sim->dL, 2.) * Sim->tau_esc;
    
    // Iterate over the array of intrinsic frequencies
    for (int64_t i = 1; i <= Sim->array_len; i++)
    {
        // Calculate observed eps for each intrinsic eps
        eps_obs = Sim->doppler_factor * Sim->Photons->eps[i] / (1. + Sim->z);
        // Calculate the flux in the observer frame
        Sim->nu_flux[i] = 
            pre_fac * (pow(Sim->Photons->eps[i], 2.) * Sim->Photons->n[i]);
        
        Sim->flux_eps[i] = eps_obs;
    }
}


int find_closest_index(double *array,int64_t array_len, double target) {
    int64_t closest_index = 0;
    double minDiff = fabs(array[0] - target);
    for (int i = 1; i < array_len; i++) {
        double diff = fabs(array[i] - target);
        if (diff < minDiff) {
            minDiff = diff;
            closest_index = i;
        }
    }

    return closest_index;
}

void simulate(char *filepath, SimulationParams *Sim)
{
    malloc_Sim_arrays(Sim);
    char header[100];
    struct timeval start1, end1, start2, end2;
    gettimeofday(&start1, NULL);

    // write gammas to csv
    write_column_to_csv(filepath, Sim->Electrons->gamma, Sim->array_len+2, "gamma");
    
    // take input params and calculate coefficients
    calc_S(Sim);
    normalize_inject_dist(Sim->inject_power, Sim);
    calc_Q_e0(Sim);
    calc_tau_esc(Sim);
    
    // find the break point for the upwind and downwind solving
    Sim->gamma_eq = -1. / (Sim->tau_acc * Sim->S);
    Sim->gamma_eq_break = find_closest_index(Sim->Electrons->gamma, Sim->array_len, Sim->gamma_eq);
    
    // ensure dt is not larger than tau_esc for stability
    if (Sim->dt > Sim->tau_esc)
    {
        Sim->dt = Sim->tau_esc;
        printf("dt larger than tau_esc. dt set to tau_esc to ensure stability\n");
    }

    set_initial_state(Sim);
    Sim->t = 0.;
    Sim->iter = 0;
    Sim->end_sim = false;
    // start simulation
    printf("Start Sim with C %e, tau %e B %.2lf, S %e g_array %lld, dt %.2e,\n", 
    Sim->Q_e0 * Sim->norm, Sim->tau_esc, Sim->B, Sim->S, Sim->array_len, Sim->dt);
    // write initial state to file
    //sprintf(header, "n_e t=%e", Sim->t);
    //write_column_to_csv(filepath, Sim->Electrons->current_n, Sim->array_len+2, header);
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
    gettimeofday(&end1, NULL);
    Sim->solve_time = (end1.tv_sec - start1.tv_sec) + (end1.tv_usec - start1.tv_usec) / 1000000.0;
    printf("stepping time: %e\n", Sim->solve_time);

    Sim->final_time=Sim->t;
    sprintf(header, "electron_n", Sim->t);
    write_column_to_csv(filepath, Sim->Electrons->current_n, Sim->array_len+2, header);

    // generate photon population
    photon_calc(Sim);
    sprintf(header, "j_nu");
    write_column_to_csv(filepath, Sim->Photons->j_nu, Sim->array_len+2, header);
    sprintf(header, "photon_eps");
    write_column_to_csv(filepath, Sim->Photons->eps, Sim->array_len+2, header);
    sprintf(header, "photon_n");
    write_column_to_csv(filepath, Sim->Photons->n, Sim->array_len+2, header);

    calc_flux(Sim);
    sprintf(header, "flux_eps");
    write_column_to_csv(filepath, Sim->flux_eps, Sim->array_len+2, header);
    sprintf(header, "nu_flux");
    write_column_to_csv(filepath, Sim->nu_flux, Sim->array_len+2, header);

    gettimeofday(&end2, NULL);
    Sim->sim_time = (end2.tv_sec - start2.tv_sec) + (end2.tv_usec - start2.tv_usec) / 1000000.0;
    printf("whole time: %e\n", Sim->sim_time);
}

void write_gammas_to_file(FILE *file, SimulationParams *Sim)
{
    // print headers in csv file
    fprintf(file, "gamma,");
    for (int64_t i = 0; i <= Sim->array_len+1; i++)
    {
        fprintf(file, "%lf,", Sim->Electrons->gamma[i]);
    }
    fprintf(file, "\n");
    fflush(file);
}

void write_run_file(char *filename, SimulationParams *Sim)
{
    char filepath[100];
    sprintf(filepath, "csv_data/implicit_solve/runs/run_%s", filename);

    FILE *run_file = fopen(filepath, "w");
    
    fprintf(run_file, 
    "dt,R,inject_p,inject_min,inject_max,rho,B,L,end_tol,Q_e0,S,tau_esc,norm,avg_gamma,V,array_len,max_gamma,min_gamma,init_p,samples_per_decade,final_time,change,tau_acc,whole_t,stepping_t\n");
    fprintf(run_file,
    "%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%lld,%e,%e,%e,%lld,%e,%e,%e,%e,%e\n",
    Sim->dt,Sim->R,Sim->inject_power,Sim->inject_min,Sim->inject_max,Sim->rho,Sim->B, 
    Sim->L,Sim->end_tol, Sim->Q_e0,Sim->S,Sim->tau_esc,Sim->norm,Sim->avg_gamma,Sim->V,
    Sim->array_len,Sim->max_gamma,Sim->min_gamma,Sim->init_power,Sim->samples_per_decade,
    Sim->final_time,Sim->change, Sim->tau_acc,Sim->sim_time,Sim->solve_time);
    fclose(run_file);
}

int main()
{
    double hold;
    SimulationParams *Sim = malloc(sizeof(SimulationParams));
    // simulation setup  
    
    // acceleration test params
    Sim->inject_min = 1e1;
    Sim->inject_max = 1e2;
    Sim->inject_power = 2.3;
    Sim->B = 0.0001;
    Sim->R = 1e16;
    Sim->L = 1e30;
    Sim->tau_acc = calc_tau_esc(Sim) * 4.;
    Sim->I = &power_law;
    Sim->LHS_BC = 5e-10;    
    

    /*
    // cooling test params
    Sim->inject_min = 1e4;
    Sim->inject_max = 1e8;
    Sim->inject_power = 2.3;
    Sim->B = 0.1;
    Sim->R = 1e16;
    Sim->L = 1e30;
    Sim->rho = 1e-38;
    Sim->end_tol = 1e-8;
    Sim->tau_acc = 1e99;
    Sim->I=&power_law;
    */

    /*
    Sim->inject_min = pow(10, 3.95);
    Sim->inject_break = pow(10, 4.26);
    Sim->inject_max = pow(10, 7.68);
    Sim->inject_power = 1.22;
    Sim->inject_power_2 = 3.37;
    Sim->B = pow(10, -1.36);
    Sim->R = pow(10, 16.12);
    Sim->L = pow(10, 41.54);
    Sim->doppler_factor = pow(10, 1.83);
    Sim->tau_acc = 1e99;
    Sim->z = 0.33;
    Sim->I = &broken_power_law;
    Sim->LHS_BC = 0.;
    */

    Sim->min_gamma = 1e1;
    Sim->max_gamma = 1e8;
    Sim->min_freq = 1e6;
    Sim->max_freq = 1e17;
    Sim->init_power = 2.;
    Sim->samples_per_decade = 80;
    Sim->dt = 1e100; //calc_tau_esc(Sim);
    Sim->end_t = 1e10;
    Sim->end_tol = 1e-8;

    // generate the cooling test data
    /*
    hold = Sim->B;
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
    Sim->B = hold;
    */
    
    /*
    // generate the acceleration test data
    hold = Sim->tau_acc;
    double t_acc_multi[5] = {0.5,1.0,1.5,2.0,4.0};
    for (int i =0; i < 5; i++)
    {
        Sim->tau_acc = t_acc_multi[i] * Sim->tau_esc;
        // generate file name based on B
        char filename[100], filepath[100];
        sprintf(filename, "t_acc=%1.1lf_t_esc.csv", t_acc_multi[i]);
        printf("generating %s\n",filename);
        sprintf(filepath, "csv_data/implicit_solve/%s", filename);
        FILE *file = fopen(filepath, "w");
        fclose(file);
        simulate(filepath, Sim);

        write_run_file(filename, Sim);
    }
    Sim->tau_acc = hold;
    */
   
    /*
    // test time taken vs samples per decade
    hold = Sim->samples_per_decade;
    for (int i =1; i < 34; i++)
    {
        Sim->samples_per_decade = 16 * i;
        // generate file name based on B
        char filename[100], filepath[100];
        sprintf(filename, "samples_per_dec_%lld.csv", 16 * i);
        printf("generating %s\n",filename);
        sprintf(filepath, "csv_data/implicit_solve/%s", filename);
        FILE *file = fopen(filepath, "w");
        fclose(file);
        simulate(filepath, Sim);

        write_run_file(filename, Sim);
    }
    Sim->samples_per_decade = hold;
    */

    char filename[100], filepath[100];
    sprintf(filename, "simulation_data.csv");
    sprintf(filepath, "csv_data/implicit_solve/%s", filename);
    FILE *file = fopen(filepath, "w");
    fclose(file);

    simulate(filepath, Sim);
    write_run_file(filename, Sim);

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