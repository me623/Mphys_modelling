
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
    double *n;
    double *prev_n;
    double *gamma;
    double *dgamma_forward;
    double delta_ln_gamma;
} ElectronParams;

typedef struct PhotonParams
{
    double *n;
    double *eps;
    double *j_nu;
    double *alpha_nu;
} PhotonParams;

typedef struct SimulationParams
{
    // array parameters
    int64_t array_len;
    double max_gamma;
    double min_gamma;
    double max_eps;
    double min_eps;
    int64_t samples_per_decade;

    // free parameters
    double R;       // radius of system (spherical geo)
    double inject_power;// power law power
    double init_power;  // power to generate initial distribution
    double inject_min;  // injection gamma range min
    double inject_max;  // injection gamma range max
    double B;       // background magnetic field
    double L;       // external luminosity
    double end_tol; // tolerance for change in n to count as equilibrium
    double dL;      // luminosity distance
    double z;       // redshift
    double tau_acc;
    double doppler_factor;
    double (*I)();  // injection term!
    double LHS_BC;

    // broken power law parameters
    double inject_break;
    double inject_power_2;

    // behind the scenes parameters
    double Q_e0;    // starting population
    double S;       // sync vs inv. compton
    double tau_esc; // free escape time
    double norm;    // normalize the prob dist
    double avg_gamma;   // average gamma of injected dist for Q_e0 calc
    double V;           // volume of system (spherical based on R)   
    double gamma_eq;

    // code specific values
    ElectronParams *Electrons;
    PhotonParams *Photons;
    double *nu_flux;
    double *flux_eps;
    char *filepath;
    double change;
    bool end_sim;
    int64_t iter;
    int64_t max_iter;
    int64_t gamma_eq_break;

    double whole_time;
    double stepping_time;

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
    Electrons->dgamma_forward = malloc((Sim->array_len + 2) * sizeof(double));

    // fill gamma array with equal log step data
    for (int64_t i = 0; i <= Sim->array_len + 1; i++)
    {
        Electrons->gamma[i] = pow(10, log10(Sim->min_gamma) + (i-1) / (double)Sim->samples_per_decade);
        
        if (i > 0)
        {
            Electrons->dgamma_forward[i-1] = Electrons->gamma[i] - Electrons->gamma[i-1];
        }
    }

    Electrons->delta_ln_gamma = log(Electrons->gamma[1]) - log(Electrons->gamma[0]);
}
void fill_eps_array(SimulationParams *Sim, PhotonParams *Photons)
{
    int64_t decades, samples_per_decade;
    // calculate number of decades in the gamma range
    decades = (int64_t)log10(Sim->max_eps) - log10(Sim->min_eps);
    samples_per_decade = Sim->array_len / decades;
    
    // malloc gamma and delta gamma arrays
    Photons->eps = calloc(Sim->array_len + 2, sizeof(double));
    Sim->nu_flux = calloc(Sim->array_len + 2, sizeof(double));
    Sim->flux_eps = calloc(Sim->array_len + 2, sizeof(double));
    Photons->j_nu  = calloc(Sim->array_len + 2, sizeof(double));
    Photons->alpha_nu = calloc(Sim->array_len + 2, sizeof(double));
    // fill gamma array with equal log step data
    for (int64_t i = 0; i <= Sim->array_len + 1; i++)
    {
        Photons->eps[i] = pow(10, log10(Sim->min_eps) + (i-1) / (double)samples_per_decade);
    }

    // Calculate the step size in logs
    double log_step = 1.0 / (double)Sim->samples_per_decade;
    // Extrapolate gamma[0] based on gamma[1]
    Photons->eps[0] = pow(10, log10(Sim->min_eps) - log_step);
    printf("min eps %e, max eps %e\n", Photons->eps[0], Photons->eps[Sim->array_len + 1]);
}

void malloc_Sim_arrays(SimulationParams *Sim)
{
    Sim->Electrons = malloc(sizeof(ElectronParams));
    malloc_and_fill_gamma_array(Sim, Sim->Electrons);
    Sim->Electrons->n = malloc((Sim->array_len + 2) * sizeof(double));
    Sim->Electrons->prev_n = malloc((Sim->array_len + 2) * sizeof(double));

    Sim->Photons = malloc(sizeof(PhotonParams));
    fill_eps_array(Sim, Sim->Photons);
    Sim->Photons->n = calloc((Sim->array_len + 2), sizeof(double));

}

void free_Sim_arrays(SimulationParams *Sim)
{
    free(Sim->Electrons->n);
    free(Sim->Electrons->prev_n);
    free(Sim->Electrons->gamma);
    free(Sim->Electrons);
    free(Sim->Photons->n);
    free(Sim->Photons->eps);
    free(Sim->Photons->j_nu);
    free(Sim->Photons->alpha_nu);
    free(Sim->Photons);
    free(Sim->nu_flux);
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
        return Sim->Q_e0 * Sim->norm * pow(gamma, (-1. * Sim->inject_power));
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
       for (int64_t i = 0; i <= Sim->array_len + 1; i++)
       {
        integral += 0.5 
        * (Sim->I(Sim->Electrons->gamma[i + 1], Sim) 
            + Sim->I(Sim->Electrons->gamma[i], Sim))
        * (Sim->Electrons->dgamma_forward[i]);
       }
       Sim->norm /= integral;
       
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
        * (Sim->Electrons->dgamma_forward[i]);
    }
    // calculate volume
    Sim->V = (4./3.) * M_PI * Sim->R * Sim->R * Sim->R;
    // calculate injected Qe0
    Sim->Q_e0 = Sim->L / (Sim->V * Sim->avg_gamma * m_e * c * c);
}

double calc_tau_esc(SimulationParams *Sim)
{
    // calculate escape time for a spherical plasma based on the radius
    Sim->tau_esc = (3.* Sim->R) / (4. * c);
    return Sim->tau_esc;
}

int find_closest_index(double *array,int64_t array_len, double target) {
    int64_t closest_index = 0;
    double minDiff = fabs(array[0] - target);
    printf("target %e\n", target);
    for (int i = 1; i < array_len; i++) {
        double diff = fabs(array[i] - target);
        if (diff < minDiff) {
            minDiff = diff;
            closest_index = i;
        }
    }

    return closest_index;
}


void stepping_regime(SimulationParams *Sim, ElectronParams *Electron)
{
    char header[100];
    for (int64_t i = 1; i <= Sim->gamma_eq_break; i++)
    {
        //sprintf(header, "ne_iter_%lld", i);
        //write_column_to_csv(Sim->filepath, Sim->Electrons->n, Sim->array_len+2, header);
        // explicit stepping regime
        Electron->n[i] = 
        (Sim->tau_esc * 
        (Sim->S * pow(Electron->gamma[i-1], 2) * Electron->n[i-1] * Sim->tau_acc
        + Electron->delta_ln_gamma * Sim->I(Electron->gamma[i], Sim) * Electron->gamma[i] * Sim->tau_acc
        + Electron->gamma[i-1] * Electron->n[i-1]))
        /
        (Electron->gamma[i] 
        * (Sim->S *  Electron->gamma[i] * Sim->tau_esc * Sim->tau_acc 
        + Electron->delta_ln_gamma * Sim->tau_acc 
        + Sim->tau_esc));
    }
    
    for (int64_t i = Sim->array_len; i >= Sim->gamma_eq_break+1; i--)
    {      
        //sprintf(header, "ne_iter_%lld", i);
        //write_column_to_csv(Sim->filepath, Sim->Electrons->n, Sim->array_len+2, header);      
        // explicit stepping regime
        Electron->n[i] = 
        (Sim->tau_esc * 
        (Sim->S * pow(Electron->gamma[i+1], 2) * Electron->n[i+1] * Sim->tau_acc
        - Electron->delta_ln_gamma * Sim->I(Electron->gamma[i], Sim) * Electron->gamma[i] * Sim->tau_acc
        + Electron->gamma[i+1] * Electron->n[i+1]))
        /
            (Electron->gamma[i] 
            * (Sim->S *  Electron->gamma[i] * Sim->tau_esc * Sim->tau_acc 
            - Electron->delta_ln_gamma * Sim->tau_acc 
            + Sim->tau_esc));
    }
}

void save_step_to_prev_n(SimulationParams *Sim, ElectronParams *Electrons)
{
    for (int64_t i = 0; i <= Sim->array_len + 1; i++)
    {
        Electrons->prev_n[i] = Electrons->n[i];
    }
}

void equilibrium_check(SimulationParams *Sim, ElectronParams *Electrons)
{
    Sim->change = 0.;
    // calculate percentage change in n
    for (int64_t i = 1; i < Sim->array_len - 1; i++)
    {
        if (Electrons->n[i] == 0. || Electrons->n[i+1] == 0.)
        {
            continue;
        }
        
        Sim->change += 0.5 * 
        (pow((1. / Electrons->n[i]) * ((Electrons->n[i] - Electrons->prev_n[i])), 2.)
        + pow((1. / Electrons->n[i+1]) * ((Electrons->n[i+1] - Electrons->prev_n[i+1])), 2.))
        * (Electrons->dgamma_forward[i]);
    }
    Sim->change = sqrt(Sim->change);
    // check change in population against specified end tolerance
    if (Sim->change < Sim->end_tol)
    {
        printf("equilibrium reached at step:%lld, last change: %e\n", Sim->iter, Sim->change);
        Sim->end_sim=true;
    }
}

void impose_BCs(SimulationParams *Sim, ElectronParams *Electrons)
{
    Electrons->n[0] = Sim->LHS_BC;
    Electrons->n[Sim->array_len+1] = 0.;
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

// critical frequency
double nu_crit(double gamma, SimulationParams *Sim)
{
    return (gamma * gamma * Sim->B * 3. * q) / (4. * M_PI * m_e * c);
}

// energy emitted by a particle with energy gamma at frequency mu
double P_nu(double gamma, double freq, SimulationParams *Sim)
{
    double x = (freq / nu_crit(gamma, Sim));
    return ((sqrt(3.) * pow(q, 3.) * Sim->B) / (4. * M_PI * m_e * pow(c, 2.))) * x * CS(x);
    
}

double j_nu(double freq, SimulationParams *Sim)
{
    double integral = 0.0;

    for (int64_t i = 1; i <= Sim->array_len - 1; i++)
    {
        double *gamma = Sim->Electrons->gamma;
        double *n = Sim->Electrons->n;
        double *dgamma = Sim->Electrons->dgamma_forward;
        // use trapezium rule
        integral += 0.5 * (n[i] * P_nu(gamma[i], freq, Sim) + n[i+1] * P_nu(gamma[i+1], freq, Sim)) * dgamma[i];
    }

    return integral / (4.0 * M_PI);
}

double alpha_nu_integrand(SimulationParams *Sim, double freq, int64_t i)
{
    double *gamma = Sim->Electrons->gamma;
    double *dgamma = Sim->Electrons->dgamma_forward;
    double *n = Sim->Electrons->n;
    return 
    (pow(gamma[i], 2)
    * ((n[i+1] / pow(gamma[i+1], 2) - n[i] / pow(gamma[i], 2)) / dgamma[i])
    * P_nu(gamma[i], freq, Sim));
}

double alpha_nu(double freq, SimulationParams *Sim) {
    double integral = 0.0;
    double *dgamma = Sim->Electrons->dgamma_forward;
    for (int64_t i = 1; i < Sim->array_len; i++) {
        // integrate via trapesium rule
        integral += 0.5 * (alpha_nu_integrand(Sim, freq, i+1) + alpha_nu_integrand(Sim, freq, i)) * dgamma[i];
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
        Sim->Photons->alpha_nu[i] = alpha_nu(freq_from_eps(Sim->Photons->eps[i]), Sim);

        Sim->Photons->n[i] = (2 * Sim->Photons->j_nu[i]) / (h_bar * Sim->Photons->eps[i]);
        //Sim->Photons->n[i] -= Sim->Photons->n[i] / Sim->tau_esc;
        //Sim->Photons->n[i] += c * Sim->Photons->alpha_nu[i] * Sim->Photons->n[i];
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
        // Calculate observed frequency for each intrinsic frequency
        eps_obs = Sim->doppler_factor * Sim->Photons->eps[i] / (1. + Sim->z);
        // Calculate the flux in the observer frame
        Sim->nu_flux[i] = 
            pre_fac * (pow(Sim->Photons->eps[i], 2.) * Sim->Photons->n[i]);
        
        Sim->flux_eps[i] = eps_obs;
    }
}

void clear_n_e_array(SimulationParams *Sim)
{
    for (int64_t i =0; i <= Sim->array_len+1; i++)
    {
        Sim->Electrons->n[i] = 0.;
    }
}

void simulate(char *filepath, SimulationParams *Sim)
{
    malloc_Sim_arrays(Sim);
    Sim->filepath = filepath;
    Sim->end_sim = false;
    Sim->iter = 0;
    char header[100];

    struct timeval start1, end1, start2, end2;
    gettimeofday(&start1, NULL);
    // write gammas to csv
    write_column_to_csv(filepath, Sim->Electrons->gamma, Sim->array_len+2, "gamma");
    clear_n_e_array(Sim);
    
    // take input params and calculate coefficients
    calc_S(Sim);
    normalize_inject_dist(Sim->inject_power, Sim);
    calc_Q_e0(Sim);
    calc_tau_esc(Sim);

    // set BCs
    impose_BCs(Sim, Sim->Electrons);

    // find the break point for the upwind and downwind solving
    Sim->gamma_eq = -1. / (Sim->tau_acc * Sim->S);
    Sim->gamma_eq_break = find_closest_index(Sim->Electrons->gamma, Sim->array_len, Sim->gamma_eq);
    
    // start simulation
    printf("Start Sim with C %.3e, tau %.3e B %.3lf, S %.3e g_array %lld, t_acc %.3e, gamma_eq_i %lld\n", 
    Sim->Q_e0 * Sim->norm, Sim->tau_esc, Sim->B, Sim->S, Sim->array_len, Sim->tau_acc, Sim->gamma_eq_break);
    // write first csv column
    //sprintf(header, "ne_iter_%lld", Sim->iter);
    //write_column_to_csv(filepath, Sim->Electrons->n, Sim->array_len+2, header);
    gettimeofday(&start2, NULL);
    while (Sim->end_sim == false && Sim->iter < Sim->max_iter)
    {
        save_step_to_prev_n(Sim, Sim->Electrons);
        stepping_regime(Sim, Sim->Electrons);
        equilibrium_check(Sim, Sim->Electrons);
        impose_BCs(Sim, Sim->Electrons);
        //sprintf(header, "ne_iter_%e", Sim->iter);
        //write_column_to_csv(filepath, Sim->Electrons->n, Sim->array_len+2, header);
        if (Sim->end_sim)
        {
            break;
        }
        
        Sim->iter ++;
    }
    gettimeofday(&end1, NULL);
    Sim->stepping_time = (end1.tv_sec - start1.tv_sec) + (end1.tv_usec - start1.tv_usec) / 1000000.0;
    printf("stepping time: %e\n", Sim->stepping_time);
    
    sprintf(header, "electron_n");
    write_column_to_csv(filepath, Sim->Electrons->n, Sim->array_len+2, header);
    
    // generate photon population
    photon_calc(Sim);
    sprintf(header, "photon_eps");
    write_column_to_csv(filepath, Sim->Photons->eps, Sim->array_len+2, header);
    sprintf(header, "photon_n");
    write_column_to_csv(filepath, Sim->Photons->n, Sim->array_len+2, header);
    sprintf(header, "j_nu");
    write_column_to_csv(filepath, Sim->Photons->j_nu, Sim->array_len+2, header);
    sprintf(header, "alpha_nu");
    write_column_to_csv(filepath, Sim->Photons->alpha_nu, Sim->array_len+2, header);

    calc_flux(Sim);
    sprintf(header, "flux_eps");
    write_column_to_csv(filepath, Sim->flux_eps, Sim->array_len+2, header);
    sprintf(header, "nu_flux");
    write_column_to_csv(filepath, Sim->nu_flux, Sim->array_len+2, header);
    
    gettimeofday(&end2, NULL);
    Sim->whole_time = (end2.tv_sec - start2.tv_sec) + (end2.tv_usec - start2.tv_usec) / 1000000.0;
    printf("whole time: %e\n", Sim->whole_time);
}

void write_run_file(SimulationParams *Sim, char *filename)
{
    char filepath[100];
    sprintf(filepath, "csv_data/steady_state/runs/run_%s", filename);

    FILE *file = fopen(filepath, "w");
    fprintf(file, 
    "delta_ln_gamma,R,inject_p,inject_min,inject_max,B,L,end_tol,Q_e0,S,tau_esc,norm,avg_gamma,V,array_len,max_gamma,min_gamma,init_p,samples_per_decade,change,crit_freq,tau_acc,iter,gamma_eq,gamma_eq_break,whole_t,stepping_t\n");
    fprintf(file,
    "%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%lld,%e,%e,%e,%lld,%e,%e,%e,%lld,%e,%lld,%e,%e\n",
    Sim->Electrons->delta_ln_gamma, Sim->R,Sim->inject_power,Sim->inject_min,Sim->inject_max,Sim->B,Sim->L,
    Sim->end_tol, Sim->Q_e0,Sim->S,Sim->tau_esc,Sim->norm,Sim->avg_gamma,Sim->V,Sim->array_len,
    Sim->max_gamma,Sim->min_gamma,Sim->init_power,Sim->samples_per_decade, Sim->change,nu_crit(Sim->inject_min, Sim),
    Sim->tau_acc, Sim->iter, Sim->gamma_eq, Sim->gamma_eq_break, Sim->whole_time, Sim->stepping_time);
    fclose(file);
}

int main()
{
    double hold;
    SimulationParams *Sim = malloc(sizeof(SimulationParams));
    // simulation setup
    // free params
    /*
    // j_nu confirmation
    Sim->inject_min = 1e1;
    Sim->inject_max = 1e4;
    Sim->inject_power = 2.5;
    Sim->B = 1.;
    Sim->R = 1e16;
    Sim->L = 1e30;
    Sim->doppler_factor = pow(10, 1.83);
    Sim->z = 0.33;
    Sim->tau_acc = 1e99;
    Sim->I = &power_law;
    Sim->LHS_BC = 0.;
    */
    /*
    // KATU comparison Injection values
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
    */
    
    // KATU comparison Steady State values
    Sim->inject_min = pow(10, 1.33);
    Sim->inject_break = pow(10, 4.03);
    Sim->inject_max = pow(10, 6.82);
    Sim->inject_power = 1.69;
    Sim->inject_power_2 = 4.29;
    Sim->B = pow(10, -1.01);
    Sim->R = pow(10, 16.46);
    Sim->L = pow(10, 49.5);
    Sim->doppler_factor = pow(10, 1.44);
    Sim->tau_acc = 1e256;
    Sim->z = 0.33;
    Sim->I = &broken_power_law;
    Sim->LHS_BC = 0.;
    
    /*
    //Acceleration test
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
    */
    // array params
    Sim->min_gamma = 1e1;
    Sim->max_gamma = 1e8;
    Sim->min_eps = 1e-12;
    Sim->max_eps = 1e8;
    Sim->samples_per_decade = 80;
    Sim->end_tol = 1e-8;
    Sim->max_iter = 2;
    
    /*
    // generate the cooling test data
    hold = Sim->B;
    double B[6] = {0.1,0.25,0.5,1.,1.5,2.};
    for (int i =0; i < 6; i++)
    {
        Sim->B = B[i];
        // generate file name based on B
        char filename[100], filepath[100];
        sprintf(filename, "B%4.0lf.csv", Sim->B*1000.);
        sprintf(filepath, "csv_data/steady_state/%s", filename);
        FILE *file = fopen(filepath, "w");
        fclose(file);
        simulate(filepath, Sim);

        write_run_file(Sim, filename);
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
        sprintf(filepath, "csv_data/steady_state/%s", filename);
        FILE *file = fopen(filepath, "w");
        fclose(file);
        simulate(filepath, Sim);

        write_run_file(Sim, filename);
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
        sprintf(filepath, "csv_data/steady_state/%s", filename);
        FILE *file = fopen(filepath, "w");
        fclose(file);
        simulate(filepath, Sim);

        write_run_file(Sim, filename);
    }
    Sim->samples_per_decade = hold;
    */

    char filename[100], filepath[100];
    sprintf(filename, "simulation_data.csv");
    sprintf(filepath, "csv_data/steady_state/%s", filename);
    FILE *file = fopen(filepath, "w");
    fclose(file);

    simulate(filepath, Sim);

    write_run_file(Sim, filename);
    
   /*
    hold = Sim->inject_min;
    double param[6] = {1e3,2.5e3,5e3,1e4,2.5e4,5e4};
    for (int i =0; i < 6; i++)
    {
        Sim->inject_min = param[i];
        // generate file name based on B
        char filename[100], filepath[100];
        sprintf(filename, "inject_min%.0lf.csv", Sim->inject_min);
        sprintf(filepath, "csv_data/steady_state/%s", filename);
        FILE *file = fopen(filepath, "w");
        fclose(file);
        simulate(filepath, Sim);

        write_run_file(Sim, filename);
    }
    Sim->inject_min = hold;
    */
    // end program
    free_Sim_arrays(Sim);
    free(Sim);
    return 0;
}