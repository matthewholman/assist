/**
 * A simple benchmark to measure the speed of the integration (no interpolation is done).
 */
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <unistd.h>
#include "rebound.h"
#include "assist.h"

double runtime(){
    struct assist_ephem* ephem = assist_ephem_create(
	    "../../data/linux_p1550p2650.440",
            "../../data/sb441-n16.bsp");
    if (ephem == NULL){
        fprintf(stderr,"Error initializing assist_ephem.\n");
        exit(1);
    }
    struct reb_simulation* r = reb_create_simulation();
    struct assist_extras* ax = assist_attach(r, ephem);
    r->t = 8416.5;

    // Asteroid Holman with slightly different initial conditions
    for (int i=0; i<10; i++){
        reb_add_fmt(r, "x y z vx vy vz",
                -2.724183384883979E+00+i/1e-10, -3.523994546329214E-02, 9.036596202793466E-02, 
                -1.374545432301129E-04, -1.027075301472321E-02, -4.195690627695180E-03); 
    }
   
    struct timeval time_beginning;
    struct timeval time_end;
    
    gettimeofday(&time_beginning,NULL);

    reb_integrate(r, r->t + 10*365.25); // 10 years
    
    gettimeofday(&time_end,NULL);
    
    assist_free(ax);
    assist_ephem_free(ephem);
    reb_free_simulation(r);

    return time_end.tv_sec-time_beginning.tv_sec+(time_end.tv_usec-time_beginning.tv_usec)/1e6;
}

int main(int argc, char* argv[]){

    int N_samples = 5;
    double mean = 0;
    double M2 = 0;

    for (int i=0; i<N_samples; i++){
        double r = runtime();
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

