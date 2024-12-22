
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <stdbool.h>
#include <windows.h>
#include <sys/time.h>
#include <gsl/gsl_sf_hyperg.h>
#include <gsl/gsl_sf_bessel.h> // GNU Scientific Library for Bessel functions

// code settings
#define WRITE_STEPS 0
#define TOLERANCE 1e-6  // Tolerance for numerical integration
#define INFINITY_CUTOFF 100.0  // Approximate "infinity" for inner integral
#define APPROXIMATE_XCS 1   // use approximation avoiding long Whitaker function calculation
#define THEORY_LEN 917504

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
    double *theory_n;
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

// Trim whitespace from both ends of a string
void trim(char *str) 
{
    char *start = str;
    char *end = str + strlen(str) - 1;

    // Trim leading whitespace
    while (isspace(*start)) start++;

    // Trim trailing whitespace
    while (end > start && isspace(*end)) end--;

    // Shift the trimmed string to the start
    memmove(str, start, end - start + 1);
    str[end - start + 1] = '\0';
}

int read_column_from_csv(const char *filename, double *data, int rows, const char *header)
{
    FILE *file = fopen(filename, "r");
    if (file == NULL) {
        fprintf(stderr, "Error: Unable to open file %s\n", filename);
        return -1;
    }

    char line[1024];  // Buffer to read each line
    char *token;
    int header_index = -1;
    int current_row = 0;

    // Read the header line to find the column index
    if (fgets(line, sizeof(line), file) != NULL) {
        // Copy the line to avoid modifying the original
        char *header_line = strdup(line);
        
        // Tokenize the header line
        token = strtok(header_line, ",");
        int col_index = 0;

        while (token != NULL) {
            // Trim whitespace from the token
            trim(token);

            // Debug print
            //printf("Found header: '%s'\n", token);

            if (strcmp(token, header) == 0) {
                header_index = col_index;
                break;
            }
            token = strtok(NULL, ",");
            col_index++;
        }

        free(header_line);

        if (header_index == -1) {
            fprintf(stderr, "Error: Header '%s' not found\n", header);
            fprintf(stderr, "CSV headers appear to be: %s", line);
            fclose(file);
            return -1;
        }
    }

    // Read data rows
    while (fgets(line, sizeof(line), file) && current_row < rows) {
        // Copy the line to avoid modifying the original
        char *data_line = strdup(line);
        
        token = strtok(data_line, ",");
        
        // Navigate to the correct column
        for (int i = 0; i < header_index; i++) {
            token = strtok(NULL, ",");
            if (token == NULL) break;
        }

        // Convert and store the value
        if (token != NULL) {
            trim(token);
            data[current_row] = atof(token);
            current_row++;
        }

        free(data_line);
    }

    fclose(file);

    // Check if we read the expected number of rows
    if (current_row < rows) {
        fprintf(stderr, "Warning: Only read %d rows instead of %d\n", current_row, rows);
        return current_row;
    }

    return rows;
}

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

void accuracy_check(SimulationParams *Sim, ElectronParams *Electrons)
{
    int64_t N = 0;
    Sim->change = 0.;
    // calculate RMS difference in n
    for (int64_t i = 1; i <= Sim->array_len; i++)
    {
        N++;
        if (Sim->Electrons->next_n[i] < 1e-256 || Sim->Electrons->theory_n[i] < 1e-256)
        {
            continue;
        }
        Sim->change += pow(Electrons->next_n[i] - Electrons->theory_n[i], 2.);
    }
    Sim->change /= N;
    Sim->change = sqrt(Sim->change);
    // check err against specified end tolerance
    if (Sim->change <= Sim->end_tol && Sim->end_sim == false)
    {
        printf("RMS err reached at step:%lld, last change: %e\n", Sim->iter, Sim->change);
        Sim->end_sim=true;
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
    Electrons->theory_n = malloc((Sim->array_len + 2) * sizeof(double));

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
    free(Sim->Electrons->theory_n);
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

// Function to find the closest index using binary search
int find_closest_index(double *arr, int len, double value) {
    int low = 0, high = len - 1;

    while (low < high) {
        int mid = (low + high) / 2;

        // Adjust binary search bounds
        if (arr[mid] < value) {
            low = mid + 1;
        } else {
            high = mid;
        }
    }

    // Edge case: closest can be between low-1 and low
    if (low > 0 && fabs(arr[low - 1] - value) < fabs(arr[low] - value)) {
        return low - 1;
    }
    return low;
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
        printf("%s\n",header);
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
    sprintf(filepath, "time_testing/implicit_solve/runs/run_%s", filename);

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

void set_theory_value(SimulationParams *Sim)
{  
    int64_t index;
    double *theory_n = malloc(sizeof(double) * THEORY_LEN);
    double *theory_gamma = malloc(sizeof(double) * THEORY_LEN);

    read_column_from_csv("csv_data/steady_state/best_value.csv", theory_n, THEORY_LEN, "electron_n");
    read_column_from_csv("csv_data/steady_state/best_value.csv", theory_gamma, THEORY_LEN, "gamma");

    for (int64_t i = 1; i < Sim->array_len+1; i++)
    {
        index = find_closest_index(theory_gamma, THEORY_LEN, Sim->Electrons->gamma[i]);
        Sim->Electrons->theory_n[i] = theory_n[index];
    }
}

void simulate(char *filename, SimulationParams *Sim)
{
    struct timeval start1, end1, start2, end2,starterr,enderr;
    // start simulation timer
    gettimeofday(&start1, NULL);

    // allocate simulation memory
    malloc_Sim_arrays(Sim);
    char header[100];

    // write gammas to csv
    // create file
    Sim->filename = filename;
    Sim->filepath = malloc(sizeof(char) * 100);
    sprintf(Sim->filepath, "time_testing/implicit_solve/%s", filename);
    /*
    FILE *file = fopen(Sim->filepath, "w");
    // write gammas to csv
    write_column_to_csv(file, Sim->Electrons->gamma, Sim->array_len + 2, "gamma");
    fclose(file);
    */
    gettimeofday(&starterr, NULL);
    set_theory_value(Sim);
    gettimeofday(&enderr, NULL);

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
    gettimeofday(&start2, NULL);
    while (Sim->t < Sim->end_t && Sim->end_sim == false)
    {
        implicit_step(Sim, Sim->Electrons);
        accuracy_check(Sim, Sim->Electrons);
        save_step_to_prev_n(Sim, Sim->Electrons);
        Sim->t += Sim->dt;
        Sim->iter ++;        
    }
    gettimeofday(&end2, NULL);
    Sim->solve_time = (end2.tv_sec - start2.tv_sec) + (end2.tv_usec - start2.tv_usec) / 1000000.0;
    
    Sim->final_time=Sim->t;
    // generate photon population
    //photon_calc(Sim);
    // generate the flux array based on photons
    //calc_flux(Sim);

    gettimeofday(&end1, NULL);
    /*
    FILE *file2 = fopen(Sim->filepath, "r+");
    write_column_to_csv(file2, Sim->Electrons->next_n, Sim->array_len + 2, "electron_n");
    write_column_to_csv(file2, Sim->Electrons->theory_n, Sim->array_len + 2, "theory");
    //write_column_to_csv(file2, Sim->Photons->eps, Sim->array_len + 2, "photon_eps");
    //write_column_to_csv(file2, Sim->Photons->n, Sim->array_len + 2, "photon_n");
    //write_column_to_csv(file2, Sim->Photons->emission, Sim->array_len + 2, "emission");
    //write_column_to_csv(file2, Sim->Photons->absorption, Sim->array_len + 2, "absorption");
    //write_column_to_csv(file2, Sim->flux_eps, Sim->array_len + 2, "flux_eps");
    //write_column_to_csv(file2, Sim->nu_flux, Sim->array_len + 2, "nu_flux");
    fclose(file2);
    */
    Sim->sim_time = ((end1.tv_sec - start1.tv_sec) + (end1.tv_usec - start1.tv_usec) / 1000000.0)
                    - ((enderr.tv_sec - starterr.tv_sec) + (enderr.tv_usec - starterr.tv_usec) / 1000000.0);
    printf("solving time: %e\n", Sim->solve_time);
    printf("sim time: %e\n", Sim->sim_time);
    write_run_file(Sim->filename, Sim);
}

int main()
{
    double hold;
    SimulationParams *Sim = malloc(sizeof(SimulationParams));
    
    // Cooling test
    Sim->nu = 1e-10;
    Sim->inject_min = 1e4;
    //Sim->inject_break = 1e99;
    Sim->inject_max = 1e8;
    Sim->inject_power = 2.3;
    //Sim->inject_power_2 = 4.29;
    Sim->B = 0.1;
    Sim->R = 1e16;
    Sim->L = 1e30;
    Sim->doppler_factor = pow(10, 1.44);
    Sim->tau_acc = 1e256;
    Sim->z = 0.33;
    Sim->I = &power_law;
    Sim->LHS_BC = 0.;

    // simulation setup  
    Sim->min_gamma = 1e1;
    Sim->max_gamma = 1e8;
    Sim->min_eps = 1e-12;
    Sim->max_eps = 1e8;
    Sim->samples_per_decade = 40;
    Sim->init_power = 2.;
    Sim->dt = 1e100;    // max out the time step
    Sim->end_t = 1e8;
    Sim->end_tol = 1e-8;


    // test time taken vs samples per decade
    hold = Sim->samples_per_decade;
    double errors[24] = {4.543821e-18, 4.362840e-18, 3.791769e-18, 3.096851e-18,
       2.392069e-18, 1.811216e-18, 1.335318e-18, 9.770250e-19,
       7.089302e-19, 5.104965e-19, 3.659197e-19, 2.608463e-19,
       1.855553e-19, 1.315693e-19, 9.291046e-20, 6.546602e-20,
       4.593581e-20, 3.210811e-20, 2.228042e-20, 1.536430e-20,
       1.041005e-20, 7.021828e-21, 4.464172e-21, 2.959807e-21};

    for (int iter=1; iter <= 50; iter++)
    {
    for (int i =8; i < 28; i++)
    {
        Sim->samples_per_decade = (int) pow(2., (float) i / 2.);
        Sim->end_tol = errors[i-8];
        // generate file name based on B
        char filename[100], filepath[100];
        sprintf(filename, "%lld_samples_per_dec_%lld.csv", iter, Sim->samples_per_decade);
        printf("generating %s\n",filename);
        simulate(filename, Sim);
        Sleep(5);
    }
    }
    Sim->samples_per_decade = hold;

    // end program
    free_Sim_arrays(Sim);
    free(Sim);
    return 0;
}