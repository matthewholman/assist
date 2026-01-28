/**
 * A simple benchmark to measure the speed of the integration (no interpolation is done).
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>
#include <math.h>
#include <time.h>
#include <unistd.h>
#include "rebound.h"
#include "assist.h"

static void usage(const char* prog){
    fprintf(stderr,
            "Usage: %s [--n_particles N] [--years Y] [--samples S] [--t0 T]\n"
            "\n"
            "Defaults match the historical benchmark:\n"
            "  N=10 particles, Y=10 years, S=5 samples, T=8416.5 (simulation time)\n",
            prog);
}

static int parse_int_arg(const char* prog, const char* flag, const char* s, int* out){
    if (s == NULL){
        fprintf(stderr, "Missing value for %s\n", flag);
        usage(prog);
        return 0;
    }
    errno = 0;
    char* end = NULL;
    long v = strtol(s, &end, 10);
    if (errno != 0 || end == s || *end != '\0' || v < 0 || v > 100000000){
        fprintf(stderr, "Invalid value for %s: '%s'\n", flag, s);
        usage(prog);
        return 0;
    }
    *out = (int)v;
    return 1;
}

static int parse_double_arg(const char* prog, const char* flag, const char* s, double* out){
    if (s == NULL){
        fprintf(stderr, "Missing value for %s\n", flag);
        usage(prog);
        return 0;
    }
    errno = 0;
    char* end = NULL;
    double v = strtod(s, &end);
    if (errno != 0 || end == s || *end != '\0' || !isfinite(v)){
        fprintf(stderr, "Invalid value for %s: '%s'\n", flag, s);
        usage(prog);
        return 0;
    }
    *out = v;
    return 1;
}

static double runtime(int n_particles, double years, double t0){
    struct assist_ephem* ephem = assist_ephem_create(
	    "../../data/de440.bsp",
            "../../data/sb441-n16.bsp");
    if (ephem == NULL){
        fprintf(stderr,"Error initializing assist_ephem.\n");
        exit(1);
    }
    struct reb_simulation* r = reb_simulation_create();
    struct assist_extras* ax = assist_attach(r, ephem);
    r->t = t0;

    // Asteroid Holman with slightly different initial conditions
    for (int i=0; i<n_particles; i++){
        reb_simulation_add_fmt(r, "x y z vx vy vz",
                -2.724183384883979E+00+i/1e-10, -3.523994546329214E-02, 9.036596202793466E-02, 
                -1.374545432301129E-04, -1.027075301472321E-02, -4.195690627695180E-03); 
    }
   
    struct timeval time_beginning;
    struct timeval time_end;
    
    gettimeofday(&time_beginning,NULL);

    reb_simulation_integrate(r, r->t + years*365.25);
    
    gettimeofday(&time_end,NULL);
    
    assist_free(ax);
    assist_ephem_free(ephem);
    reb_simulation_free(r);

    return time_end.tv_sec-time_beginning.tv_sec+(time_end.tv_usec-time_beginning.tv_usec)/1e6;
}

int main(int argc, char* argv[]){

    // Defaults (match historical benchmark)
    int n_particles = 10;
    double years = 10.0;
    int N_samples = 5;
    double t0 = 8416.5;

    // Parse CLI args
    for (int i = 1; i < argc; i++){
        const char* a = argv[i];
        if (strcmp(a, "--help") == 0 || strcmp(a, "-h") == 0){
            usage(argv[0]);
            return 0;
        } else if (strcmp(a, "--n_particles") == 0 || strcmp(a, "-n") == 0){
            if (!parse_int_arg(argv[0], a, (i+1 < argc) ? argv[i+1] : NULL, &n_particles)) return 2;
            i++;
        } else if (strcmp(a, "--years") == 0 || strcmp(a, "-y") == 0){
            if (!parse_double_arg(argv[0], a, (i+1 < argc) ? argv[i+1] : NULL, &years)) return 2;
            i++;
        } else if (strcmp(a, "--samples") == 0 || strcmp(a, "-s") == 0){
            if (!parse_int_arg(argv[0], a, (i+1 < argc) ? argv[i+1] : NULL, &N_samples)) return 2;
            i++;
        } else if (strcmp(a, "--t0") == 0){
            if (!parse_double_arg(argv[0], a, (i+1 < argc) ? argv[i+1] : NULL, &t0)) return 2;
            i++;
        } else {
            fprintf(stderr, "Unknown argument: %s\n", a);
            usage(argv[0]);
            return 2;
        }
    }

    if (N_samples < 2){
        fprintf(stderr, "--samples must be >= 2 (needed for stddev)\n");
        return 2;
    }
    if (n_particles < 0){
        fprintf(stderr, "--n_particles must be >= 0\n");
        return 2;
    }
    if (!(years >= 0.0)){
        fprintf(stderr, "--years must be >= 0\n");
        return 2;
    }

    double mean = 0;
    double M2 = 0;

    for (int i=0; i<N_samples; i++){
        double r = runtime(n_particles, years, t0);
        double delta = r - mean;
        mean += delta/(i+1);
        M2 += delta*(r-mean);
    }
    double std = sqrt(M2 / (N_samples -1));


    char hostname[1024];
    gethostname(hostname, 1024);
    char filename[2048];
    sprintf(filename, "benchmark_%s.txt", hostname); 

    FILE* of = fopen(filename, "a");
    printf("%.6fs %.6fs %s %s %s\n", mean, std, assist_version_str, assist_build_str, assist_githash_str);
    fprintf(of, "%.6fs %.6fs %s %s %s\n", mean, std, assist_version_str, assist_build_str, assist_githash_str);
    fclose(of);
    
}

