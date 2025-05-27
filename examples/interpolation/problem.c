/**
 * Interpolation
 *
 * This example demonstrates how to interpolate between timesteps using 
 * ASSIST's built in interpolation feature.
 */
#include <stdio.h>
#include <stdlib.h>
#include "rebound.h"
#include "assist.h"

int main(int argc, char* argv[]){
    // Create a REBOUND simulation.
    struct reb_simulation* r = reb_simulation_create();

    // Load the ephemeris data.
    struct assist_ephem* ephem = assist_ephem_create(
            "../../data/de440.bsp",
            "../../data/sb441-n16.bsp");
    if (!ephem){
        printf("Cannot create ephemeris structure.\n");
        exit(-1);
    }
    
    // Attach the ASSIST framework.
    // This will set the additional force routine in the REBOUND simulation,
    // the IAS15 integrator, and the gravity module.
    struct assist_extras* ax = assist_attach(r, ephem);

    // Initial time, relative to ephem->jd_ref
    r->t = 7304.5;

    // Add a particle to the REBOUND simulation (this is asteroid Holman).
    reb_simulation_add_fmt(r, "x y z vx vy vz", 3.3388753502614090e+00, -9.1765182678903168e-01, -5.0385906775843303e-01, 2.8056633153049852e-03, 7.5504086883996860e-03, 2.9800282074358684e-03);


    // To integrate forward in time, we can use the reb_simulation_integrate()
    // function. By default REBOUND has the exact_finish_time parameter
    // set to 1:
    r->exact_finish_time = 1;
    reb_simulation_integrate(r, r->t + 200.0);
    
    // So by default, REBOUND integrates to the exact finish time that
    // we've requested, here 200.0 days into the future. The IAS15
    // integrator uses as many adaptive timestep as it needs to get 
    // to the end while maintaining the required accuracy. The last
    // timestep will be short to make sure the final time is exactly 
    // 7504.5 Julian Days.
    printf("Time: %f JD\nTimestep: %f days\n", r->t, r->dt_last_done);

    // Next, let us continue the integration but turn off 
    // exact_finish_time. As a result, we will slightly overshoot the 
    // target of 7704.5 Julian Days. However, the final timestep
    // can remain larger.
    r->exact_finish_time = 0;
    reb_simulation_integrate(r, r->t + 200.0);
    printf("Time: %f JD\nTimestep: %f days\n", r->t, r->dt_last_done);

    // If you need many outputs with a short cadence, you might not
    // want to use exact_finish_time = 1 because it results in 
    // small timestep, increasing runtime, and potentially decreasing
    // accuracy. The solution is to use the built-in function 
    // assist_integrate_or_interpolate(). This function is either 
    // integrating the the time required, or if possible, uses a high
    // order polynomial to interpolate between the last two completed
    // timesteps. It is used the same polynomial that IAS15 is using
    // internally, so there is minimal overhead and you are guaranteed
    // both high speed and accuracy. 
    //
    // The following example shows how to use this function to generate
    // a serious of closely spaced outputs. The are equidistantly 
    // spaced here, but can be arbitrary, as long as they are sorted 
    // in inreasing order.
    
    double t = r->t;
    for (int i=0; i<10; i++){
        assist_integrate_or_interpolate(ax, t);
        printf("Time: %f JD\tLast timestep: %f days\n", t, r->dt_last_done);
        // Generate output every 1.5 days
        t += 1.5;  
    }

    // Clean up memory
    assist_free(ax);
    assist_ephem_free(ephem);
    reb_simulation_free(r);
}

