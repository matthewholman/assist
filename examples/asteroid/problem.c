/**
 * This example demonstrates how to use ASSIST to integrate a test particle in 
 * the field of the Sun, planets, moon, and a set of massive asteroids, whose
 * positions come from JPL's DE440/441 ephemeris.
 */
#include <stdio.h>
#include <stdlib.h>
#include "rebound.h"
#include "assist.h"

int main(int argc, char* argv[]){
    // Initial conditions of asteroid
    int n_particles = 1;    
    double* state = malloc(n_particles*6*sizeof(double));
    state[0] = 3.3388753502614090e+00;  // x in AU
    state[1] = -9.1765182678903168e-01; // y
    state[2] = -5.0385906775843303e-01; // z
    state[3] = 2.8056633153049852e-03;  // vx in AU/day
    state[4] = 7.5504086883996860e-03;  // vy
    state[5] = 2.9800282074358684e-03;  // vz
   
    double jd_ref = 2451545.0;  // Reference date in JD
    double tstart = 7304.5;     // Initial time relative to jd_ref
    double tend = 7404.5;       // Final time relative to jd_ref
    double tstep = 20;          // Integration step size. 
    
    // Interpolate between integration steps on these substeps 
    int n_substeps = 4;
    double hg[n_substeps+1];
    for(int i=0; i<=n_substeps; i++){
        hg[i]=(1.0/n_substeps)*i;
    }

    // Allocate memory for output.
    // Put in a factor to ensure enough memory allocated
    int mem_fac = 2;
    int n_alloc = ((tend-tstart)/tstep + 1)*mem_fac*n_substeps;
    double* outstate = (double *) malloc((n_alloc)*6*sizeof(double));
    double* outtime  = (double *) malloc((n_alloc+1)*sizeof(double));

    struct assist_ephem* ephem = assist_ephem_create("../../data/linux_m13000p17000.441",
		      "../../data/sb441-n16.bsp");

    int n_steps_done;
    int status = assist_integrate(
            ephem,                  // assist_ephem. NULL will try to construct it using defaults.
            tstart, tend, tstep,    // Time range of integration 
            0,                      // 1=geocentric, 0=barycentric
            1e-9,                   // epsilon
            n_particles,            // number of particles
            state,                  // initial conditions
            NULL,                   // additional particle parameters for non-gravitational effects (not used in this example)
            0, NULL, NULL, NULL,    // No variational particles in this example
            n_alloc,                // Allocated space in output array
            &n_steps_done,          // Number of steps actually done
            n_substeps,             // Number of substeps 
            hg,                     // Location of substeps (in fractions of an integration time step)
            outtime,                // Output times
            outstate,               // Output states
            0.0                     // Minimum dt in integration
            ); 
    if (status != REB_EXIT_SUCCESS) {
       printf("Warning! Simulation did *not* run successfully.\n");
    }

    // Number of outputs generated
    int n_outs = n_steps_done*n_substeps + 1;

    for (int j=0; j<n_outs; j++){
        for (int k=0; k<n_particles; k++){
            int offset = j*(k+1)*6;
            printf("%d %lf %.12lf %.12lf %.12lf %.12lf %.12lf %.12lf\n",
                    j, outtime[j],
                    outstate[offset+0], outstate[offset+1], outstate[offset+2],
                    outstate[offset+3], outstate[offset+4], outstate[offset+5]);
        }
    }
}

