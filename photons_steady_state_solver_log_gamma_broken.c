
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <stdbool.h>
#include <float.h>
#include <gsl/gsl_sf_hyperg.h>

#define BUFFER_SIZE 33554432
// CGS unit constants
#define m_e 9.109e-28
#define c 2.998e10
#define sigma_t 6.652e-25
#define q 4.803e-10 // charge of an electron in cgs
#define h 6.626e-27
#define h_bar 1.055e-27

#define H0 2.27e-18 // Hubble constant in cgs

typedef struct LeptonParams
{
    double *n;
    double *prev_n;
    double *gamma;
    double delta_ln_gamma;
} LeptonParams;

typedef struct PhotonParams
{
    double *n;
    double *frequency;
    double *j_v;
} PhotonParams;

typedef struct SimulationParams
{
    // array parameters
    int64_t array_len;
    double max_gamma;
    double min_gamma;
    double max_freq;
    double min_freq;
    int64_t samples_per_decade;

    // free parameters
    double R;       // radius of system (spherical geo)
    double inject_power;// power law pre break
    double inject_power_2;// power law post break
    double inject_min;  // injection gamma range min
    double inject_break;
    double inject_max;  // injection gamma range max
    double rho;     // background density for initial population calc
    double B;       // background magnetic field
    double L;       // external luminosity
    double end_tol; // tolerance for change in n to count as equilibrium
    double bulk_gamma;
    double z;       // redshift

    double Q_e0;    // starting population
    double S;       // sync vs inv. compton
    double tau_esc; // free escape time
    double norm;    // normalize the prob dist
    double avg_gamma;   // average gamma of injected dist for Q_e0 calc
    double V;           // volume of system (spherical based on R)   
    double dL;          // distance to plasma system
    double doppler_factor;


    // code specific values
    LeptonParams *Electrons;
    PhotonParams *Photons;
    double *nu_flux;
    double *flux_freq;
    double change;
    bool end_sim;
    int64_t iter;
    int64_t max_iter;

} SimulationParams;

void malloc_and_fill_gamma_array(SimulationParams *Sim, LeptonParams *Lepton)
{
    int64_t decades;
    // calculate number of decades in the gamma range
    decades = (int64_t)log10(Sim->max_gamma) - log10(Sim->min_gamma);
    // set array length based on decades and input samples per decade
    Sim->array_len = decades * Sim->samples_per_decade;

    // malloc gamma and delta gamma array
    Lepton->gamma = malloc((Sim->array_len + 2) * sizeof(double));
  
    // fill gamma array with equal log step data
    for (int64_t i = 0; i <= Sim->array_len + 1; i++)
    {
        Lepton->gamma[i] = pow(10, log10(Sim->min_gamma) + (i-1) / (double)Sim->samples_per_decade);
    }

    Lepton->delta_ln_gamma = log(Lepton->gamma[1]) - log(Lepton->gamma[0]);
}
void malloc_and_fill_frequency_array(SimulationParams *Sim, PhotonParams *Photons)
{
    int64_t decades, samples_per_decade;
    // calculate number of decades in the gamma range
    decades = (int64_t)log10(Sim->max_freq) - log10(Sim->min_freq);
    samples_per_decade = Sim->array_len / decades;
    
    // malloc gamma and delta gamma arrays
    Photons->frequency = calloc(Sim->array_len + 2, sizeof(double));
    Photons->j_v  = calloc(Sim->array_len + 2, sizeof(double));

    // fill gamma array with equal log step data
    for (int64_t i = 1; i <= Sim->array_len + 1; i++)
    {
        Photons->frequency[i] = pow(10, log10(Sim->min_freq) + (i-1) / (double)samples_per_decade);
    }

    // Calculate the step size in logs
    double log_step = 1.0 / (double)Sim->samples_per_decade;
    // Extrapolate gamma[0] based on gamma[1]
    Photons->frequency[0] = pow(10, log10(Sim->min_freq) - log_step);

   for (int64_t i = 1; i <= Sim->array_len; i++)
    {
        Photons->j_v[i] = Photons->frequency[i + 1] - Photons->frequency[i];
    }
}

void malloc_Sim_arrays(SimulationParams *Sim)
{
    Sim->Electrons = malloc(sizeof(LeptonParams));
    malloc_and_fill_gamma_array(Sim, Sim->Electrons);
    Sim->Electrons->n = calloc((Sim->array_len + 2), sizeof(double));
    Sim->Electrons->prev_n = malloc((Sim->array_len + 2) * sizeof(double));

    Sim->Photons = malloc(sizeof(PhotonParams));
    malloc_and_fill_frequency_array(Sim, Sim->Photons);
    Sim->Photons->n = calloc((Sim->array_len + 2), sizeof(double));
    Sim->nu_flux = calloc((Sim->array_len + 2), sizeof(double));
}

void free_Sim_arrays(SimulationParams *Sim)
{
    free(Sim->Electrons->n);
    free(Sim->Electrons->prev_n);
    free(Sim->Electrons->gamma);
    free(Sim->Electrons);
    free(Sim->Photons->n);
    free(Sim->Photons->frequency);
    free(Sim->Photons->j_v);
    free(Sim->Photons);
    free(Sim->nu_flux);
}


double I(double gamma, SimulationParams *Sim)
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
        return Sim->Q_e0 * Sim->norm * pow(gamma, -Sim->inject_power);
      }
    }
}

void normalize_power_law_dist(SimulationParams *Sim)
{
        Sim->norm = 1.;
        Sim->Q_e0 = 1.;

        double hold_min = Sim->inject_min;
        double hold_max = Sim->inject_max;

        Sim->inject_min = Sim->min_gamma;
        Sim->inject_max = Sim->max_gamma;
        double integral = 0.;
        // normalise power law dist based on a given power
        for (int64_t i = 0; i <= Sim->array_len + 1; i++)
        {
            integral += I(Sim->Electrons->gamma[i], Sim);
        }
        Sim->norm /= integral;
        Sim->inject_min = hold_min;
        Sim->inject_max = hold_max;
}

void set_initial_state(SimulationParams *Sim)
{
    Sim->Electrons->n[Sim->array_len + 1] = 0.;
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
    // calculate average gamma
    Sim->avg_gamma = 0.;
    
    for (int64_t i = 0; i <= Sim->array_len + 1; i++)
    {
        Sim->avg_gamma += Sim->Electrons->gamma[i] * I(Sim->Electrons->gamma[i], Sim);
    }

    // calculate volume
    Sim->V = (4./3.) * M_PI * Sim->R * Sim->R * Sim->R;
    // calculate injected Qe0
    Sim->Q_e0 = Sim->L / (Sim->V * Sim->avg_gamma * m_e * c * c);
}

void calc_tau_esc(SimulationParams *Sim)
{
    // calculate escape time for a spherical plasma based on the radius
    Sim->tau_esc = (3.* Sim->R) / (4. * c);
}

void calc_doppler_factor(SimulationParams *Sim)
{
    double beta = sqrt( 1 - 1 / pow(Sim->bulk_gamma , 2));
    Sim->doppler_factor = sqrt((1 + beta) / (1 - beta));
}

void fda_step(SimulationParams *Sim, LeptonParams *Lepton)
{
    for (int64_t i =  Sim->array_len; i >= 1; i--)
    {
        // explicit stepping regime
        Lepton->n[i] = 
        (Sim->tau_esc * 
        (Sim->S * Lepton->gamma[i+1] * Lepton->gamma[i+1] * Lepton->n[i+1]
        - Lepton->gamma[i] * Lepton->delta_ln_gamma 
        * I(Lepton->gamma[i], Sim)))
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
        printf("equilibrium reached at step: %lld, last change: %e\n", Sim->iter, Sim->change);
        Sim->end_sim=true;
    }
}

void impose_BCs(SimulationParams *Sim, LeptonParams *Lepton)
{
    Lepton->n[Sim->array_len+1] = 0.;
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
double vc(double gamma, SimulationParams *Sim)
{
    return ((3. * q * Sim->B) / (4. * M_PI * m_e * c)) * pow(gamma, 2.);
}

// energy emitted by a particle with energy gamma at frequency mu
double Pv(double gamma, double freq, SimulationParams *Sim)
{
    double x = (freq / vc(gamma, Sim));
    return ((sqrt(3.) * M_PI * pow(q, 3.) * Sim->B) / (2 * m_e * pow(c, 2.))) * x * CS(x);
    
}

double jv(double freq, SimulationParams *Sim)
{
    double integral = 0.;
    
    for (int64_t i = 1; i <= Sim->array_len; i++)
    {
        integral += Sim->Electrons->n[i] * Pv(Sim->Electrons->gamma[i], freq, Sim);
    }
    integral /= (4. * M_PI);
    
    return integral;
}
/*
double alphav(double freq, SimulationParams *Sim)
{
    double integral = 0.;
    
    for (int64_t i = 1; i <= Sim->array_len; i++)
    {
        integral += pow(Sim->Electrons->gamma[i], 2) * Pv(Sim->Electrons->gamma[i], freq, Sim) 
          * (Sim->Electrons->n[i+1] / Sim->Electrons->gamma[i + 1] 
            - Sim->Electrons->n[i] / Sim->Electrons->gamma[i]) 
            / 
            ( Sim->Electrons->gamma[i + 1] -  Sim->Electrons->gamma[i]);
    }
    integral /= (-1. * 8 * M_PI * m_e * freq * freq);
    
    return integral;
}
*/
double eps(double freq)
{
    return (h * freq) / (m_e * c * c);
}

void photon_calc(LeptonParams *Electrons, PhotonParams *Photons, SimulationParams *Sim)
{
    double max_val = 0.;
    char header[100];
    double j_v;
    for (int64_t i = 1; i <= Sim->array_len; i++)
    {
        j_v = jv(Photons->frequency[i], Sim);

        Photons->n[i] = (2 * j_v) / (h_bar * eps(Photons->frequency[i]));

        Photons->j_v[i] = j_v;
    }  
}

// Function to calculate the luminosity distance
double redshift_distance(double z) {
    // approx for z < 1
    return 1.32e28 * z;
}

void calc_flux(SimulationParams *Sim)
{
    Sim->dL = redshift_distance(Sim->z);

    for (int64_t i = 1; i <= Sim->array_len; i++)
    {
        Sim->nu_flux[i] = 
        (pow(eps(Sim->Photons->frequency[i]), 2) * m_e * pow(c, 2) * Sim->Photons->n[i] * pow(Sim->doppler_factor, 4) * Sim->V)
        /
        (4 * M_PI * pow(Sim->dL, 2) * (1 + Sim->z) * Sim->tau_esc);
    } 
}

void simulate(char *filepath, SimulationParams *Sim)
{
    char header[100];
    // write gammas to csv
    write_column_to_csv(filepath, Sim->Electrons->gamma, Sim->array_len+2, "gamma");
    // take input params and calculate coefficients
    calc_S(Sim);
    normalize_power_law_dist(Sim);
    calc_Q_e0(Sim);
    calc_tau_esc(Sim);

    set_initial_state(Sim);
    impose_BCs(Sim, Sim->Electrons);

    // start simulation
    printf("Start Sim with C %e, tau %e B %.2lf, S %e g_array %lld, minject%e\n", 
    Sim->Q_e0 * Sim->norm, Sim->tau_esc, Sim->B, Sim->S, Sim->array_len, Sim->inject_min);
    
    Sim->end_sim = false;
    Sim->iter = 0;
    // write first csv column
    sprintf(header, "iter %lld", Sim->iter);
    write_column_to_csv(filepath, Sim->Electrons->n, Sim->array_len+2, header);
    while (Sim->end_sim == false && Sim->iter < Sim->max_iter)
    {
        save_step_to_prev_n(Sim, Sim->Electrons);
        fda_step(Sim, Sim->Electrons);
        equilibrium_check(Sim, Sim->Electrons);
        impose_BCs(Sim, Sim->Electrons);

        if (Sim->end_sim)
        {
            break;
        }
        
        Sim->iter ++;
    }
    sprintf(header, "iter final", Sim->iter);
    write_column_to_csv(filepath, Sim->Electrons->n, Sim->array_len+2, header);
    // generate photon population
    photon_calc(Sim->Electrons, Sim->Photons, Sim);
    sprintf(header, "photon_freq");
    write_column_to_csv(filepath, Sim->Photons->frequency, Sim->array_len+2, header);
    sprintf(header, "photon_n");
    write_column_to_csv(filepath, Sim->Photons->n, Sim->array_len+2, header);
    sprintf(header, "j_v");
    write_column_to_csv(filepath, Sim->Photons->j_v, Sim->array_len+2, header);
    calc_flux(Sim);
    sprintf(header, "nu_flux");
    write_column_to_csv(filepath, Sim->nu_flux, Sim->array_len+2, header);

    printf("%e, distance pc\n", Sim->dL);
}

void write_run_file(SimulationParams *Sim, char *filename)
{
    char filepath[100];
    sprintf(filepath, "csv_data/steady_state/runs/run_%s", filename);

    FILE *file = fopen(filepath, "w");
    fprintf(file, 
    "delta_ln_gamma,R,inject_p,inject_p2,inject_min,inject_break,inject_max,rho,B,L,end_tol,Q_e0,S,tau_esc,norm,avg_gamma,V,array_len,max_gamma,min_gamma,samples_per_decade,change,bulk_gamma,doppler_factor,z,\n");
    fprintf(file,
    "%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%lld,%e,%e,%e,%lld,%e,%e,%e\n",
    Sim->Electrons->delta_ln_gamma, Sim->R,Sim->inject_power, 
    Sim->inject_power_2,Sim->inject_min,Sim->inject_break,Sim->inject_max,
    Sim->B,Sim->L,Sim->end_tol,Sim->Q_e0,Sim->S,Sim->tau_esc,
    Sim->norm,Sim->avg_gamma,Sim->V,Sim->array_len,Sim->max_gamma,
    Sim->min_gamma,Sim->samples_per_decade, Sim->change, Sim->bulk_gamma, 
    Sim->doppler_factor, Sim->z);
    fclose(file);
}

int main()
{
    double hold;
    SimulationParams *Sim = malloc(sizeof(SimulationParams));
    // simulation setup
    // free params
    Sim->inject_min = 1e4;
    Sim->inject_break = 1e99;
    Sim->inject_max = 1e8;
    Sim->inject_power = 2.3;
    Sim->inject_power_2 = 3.37;
    Sim->B = 0.1;
    Sim->R = 1e16;
    Sim->L = 1e30;
    Sim->doppler_factor = pow(10., 1.83);
    Sim->z = 0.15;
    // array params
    Sim->min_gamma = 1e1;
    Sim->max_gamma = 1e8;
    Sim->min_freq = 1e6;
    Sim->max_freq = 1e27;
    Sim->samples_per_decade = 120;
    Sim->end_tol = 1e-8;
    Sim->max_iter = 1000;

    malloc_Sim_arrays(Sim);

    // generate the cooling test data
    /*
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

    char filename[100], filepath[100];
    sprintf(filename, "simulation_data.csv");
    sprintf(filepath, "csv_data/steady_state/%s", filename);
    FILE *file = fopen(filepath, "w");
    fclose(file);
    printf("here\n");
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