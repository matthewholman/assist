#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include "rebound.h"
#include "assist.h"

int main(int argc, char* argv[]){
    struct assist_ephem* ephem = assist_ephem_create(
            "../../data/linux_p1550p2650.440",
            "../../data/sb441-n16.bsp");
    if (ephem == NULL){
        fprintf(stderr,"Error initializing assist_ephem.\n");
        exit(1);
    }
   
    // Create two copies of the initial simulation
    struct reb_simulation* r1 = reb_create_simulation();
    r1->t = 8416.5;
    reb_add_fmt(r1, "x y z vx vy vz",
        -2.724183384883979E+00, -3.523994546329214E-02, 9.036596202793466E-02, 
        -1.374545432301129E-04, -1.027075301472321E-02, -4.195690627695180E-03); 
    
    struct reb_simulation* r2 = reb_copy_simulation(r1);
   
    // Attach assist 
    struct assist_extras* ax1 = assist_attach(r1, ephem);
    struct assist_extras* ax2 = assist_attach(r2, ephem);

    // Create some random output times. 
    // These need to be sorted but don't need 
    // to be equidistantly spaced.
    double Noutput = 20;
    double* ts = malloc(sizeof(double)*Noutput);
    for (int i=0; i<Noutput; i++){
        ts[i] = r1->t + 40.*(i+1);
    }
    
    // Integrate to each time directly.
    double* x_direct = malloc(sizeof(double)*Noutput);
    for (int i=0; i<Noutput; i++){
        reb_integrate(r1, ts[i]);
        x_direct[i] = r1->particles[0].x; // store output 
    }   
    
    // Integrate or interpolate to each time.
    double* x_interpolate = malloc(sizeof(double)*Noutput);
    for (int i=0; i<Noutput; i++){
        assist_integrate_or_interpolate(ax2, ts[i]);
        x_interpolate[i] = r2->particles[0].x; 
    }   
    
    // Compare results
    for (int i=0; i<Noutput; i++){
        assert(fabs(x_direct[i]-x_interpolate[i]) < 1e-13); // Check that results agree
    }
    

    free(ts);
    free(x_direct);
    // Cleanup rebound and assist structures
    assist_free(ax1);
    assist_free(ax2);
    reb_free_simulation(r1);
    reb_free_simulation(r2);
    assist_ephem_free(ephem);

        
}

