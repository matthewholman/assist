/**
 * A unit test integrating a variational particle.
 */
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include "rebound.h"
#include "assist.h"
    
const double au2meter = 149597870700;

int main(int argc, char* argv[]){

    struct assist_ephem* ephem = assist_ephem_create(
            "../../data/de440.bsp",
            "../../data/sb441-n16.bsp");
    if (ephem == NULL){
        fprintf(stderr,"Error initializing assist_ephem.\n");
        exit(1);
    }
    
    struct reb_simulation* r = reb_simulation_create();
    struct assist_extras* ax = assist_attach(r, ephem);
    r->t = 8416.5;
    
    // Add main particle
    struct reb_particle p0 = {
        .x = -2.724183384883979E+00, 
        .y = -3.523994546329214E-02, 
        .z = 9.036596202793466E-02, 
        .vx = -1.374545432301129E-04, 
        .vy = -1.027075301472321E-02, 
        .vz = -4.195690627695180E-03}; 
    reb_simulation_add(r, p0);

    // Add perturbed particle
    double dx = 10.0/au2meter; // 10 meter
    p0.x += dx;
    reb_simulation_add(r, p0);

    // Add variational particle
    int var = reb_simulation_add_variation_1st_order(r, 0); // create a variational particle corresponding to particle with index 0
    r->particles[var].x = 1.;
   
    // Integratr for 1 day
    reb_simulation_integrate(r, r->t + 1.0);

    // Compare outputs
    {
        double diff_1 = r->particles[1].x-r->particles[0].x;
        double diff_2 = r->particles[var].x*dx;
        double rel = fabs( (diff_1-diff_2)/diff_1);
        printf("Relative difference (1 day): %e\n",rel);
        assert(rel < 1e-5); // require that relative difference less than 1e-5
    }
    
    // Integratr for 100 days
    reb_simulation_integrate(r, r->t + 100.0);
    // Compare outputs
    {
        double diff_1 = r->particles[1].x-r->particles[0].x;
        double diff_2 = r->particles[var].x*dx;
        double rel = fabs( (diff_1-diff_2)/diff_1);
        printf("Relative difference (100 days): %e\n",rel);
        assert(rel < 1e-5); // require that relative difference less than 1e-5
    }

        
    assist_free(ax);
    assist_ephem_free(ephem);
    reb_simulation_free(r);
}

