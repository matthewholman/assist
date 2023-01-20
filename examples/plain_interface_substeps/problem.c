/**
 * This example demonstrates how to use ASSIST to integrate two test particles in 
 * the field of the Sun, planets, moon, and a set of massive asteroids, whose
 * positions come from JPL's DE440/441 ephemeris.
 */
#include <stdio.h>
#include <stdlib.h>
#include "rebound.h"
#include "assist.h"

int main(int argc, char* argv[]){
    // Create a REBOUND simulation
    struct reb_simulation* r = reb_create_simulation();
    
    // Attach the ASSIST framework
    struct assist_extras* ax = assist_attach(r);

    // Initial time, relative to ax->jd_ref
    r->t = 7304.5;

    r->exact_finish_time = 0; // Let IAS overshoot. Use interpolation to get outputs.

    // Add a particle to the REBOUND simulation.
    reb_add_fmt(r, "x y z vx vy vz",
        3.3388753502614090e+00,   // x in AU
        -9.1765182678903168e-01,  // y
        -5.0385906775843303e-01,  // z
        2.8056633153049852e-03,   // vx in AU/day
        7.5504086883996860e-03,   // vy
        2.9800282074358684e-03);  // vz

    // Final integration time (relative to ax->jd_ref)
    double tend = 7404.5;  

    ///////////////////////////////////////////
    // Generate outputs with a fixed interval
    double next_output = r->t; 
    double output_interval = 8; 
    double* output_state = malloc(r->N*6);

    // Integrate until we reach tend or an error occurs
    while (r->t < tend && r->status <= 0){
        reb_integrate(r, next_output); // Do the actual integration. This might overshoot.

        while (r->t >= next_output) { // do interpolation to get outputs
            printf("t = %.3f     steps = %d\n", next_output,r->steps_done);
            double h = (next_output - (r->t - r->dt_last_done))/ r->dt_last_done;
            assist_interpolate(r, h, output_state);
            for(int i=0; i<r->N; i++){
                printf("\tx = %.12f \ty = %.12f \tz = %.12f\n", output_state[i*6+0], output_state[i*6+1], output_state[i*6+2]);
            }
            next_output += output_interval;
        }
    }

    /////////////////////////////////////////////////////////////
    // Generate outputs at nodes for Chebyshev interpolation 
    double interval = 30;
    tend += 200;
    const int N = 10;
    double nodes[N];
    for (int k=0; k<N; k++){
        nodes[k] = cos((2.0*k+1.0)*M_PI/(2.0*N))/2.0 + 0.5;
    } 
    int counter = 0; // steps through 0 to N-1
    

    // Integrate until we reach tend or an error occurs
    while (r->t < tend && r->status <= 0){
        reb_integrate(r, next_output); // Do the actual integration. This might overshoot.

        while (r->t >= next_output) { // do interpolation to get outputs
            printf("t = %.3f     steps = %d       counter = %d\n", next_output, r->steps_done, counter);
            double h = (next_output - (r->t - r->dt_last_done))/ r->dt_last_done;
            assist_interpolate(r, h, output_state);
            for(int i=0; i<r->N; i++){
                printf("\tx = %.12f \ty = %.12f \tz = %.12f\n", output_state[i*6+0], output_state[i*6+1], output_state[i*6+2]);
            }
            next_output += interval*nodes[counter];
            counter = (counter+1)%N;
        }
    }


    // Clean up memory
    assist_free(ax);
    reb_free_simulation(r);
}

