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


    // Load the ephemeris data
    struct assist_ephem* ephem = assist_ephem_create(
            "../../data/linux_m13000p17000.441",
            "../../data/sb441-n16.bsp");
    
    // One has the option to change other setting.
    ephem->jd_ref = 2451545.0; // Reference JD. This line can be commented out as this is the default.
    
    // Attach the ASSIST framework
    // This will set the additional force routine in the REBOUND simulation,
    // the IAS15 integrator, the gravity module, etc.
    struct assist_extras* ax = assist_attach(r, ephem);

    
    // Initial time, relative to ephem->jd_ref
    r->t = 7304.5;

    // Set any other integrator settings, for example a minimum timestep
    r->ri_ias15.min_dt = 0.5;


    // Add a particle to the REBOUND simulation.
    // These are the barycentric coordinates of asteroid Holman.
    reb_add_fmt(r, "x y z vx vy vz",
        3.3388753502614090e+00,   // x in AU
        -9.1765182678903168e-01,  // y
        -5.0385906775843303e-01,  // z
        2.8056633153049852e-03,   // vx in AU/day
        7.5504086883996860e-03,   // vy
        2.9800282074358684e-03);  // vz


    // Query the position of the sun at current simulation time.
    struct reb_particle sun = assist_get_particle(ephem, 0, r->t); // 0 stans dfor sun
    
    // Add another test particle using orbital parameters relative to the sun
    reb_add_fmt(r, "a e omega primary",
        2.4,   // semi-major axis in AU
        0.1,   // eccentricity 
        0.4,   // periastron in radians
        sun);
    
    
    // Query the position of the earth at current simulation time.
    struct reb_particle earth = assist_get_particle(ephem, 3, r->t); // 3 stands for earth
    printf("%f %f %f \n", earth.x, earth.y, earth.z);
    
    // Add another test particle in orbit around the Earth
    reb_add_fmt(r, "a primary",
        0.0002,   // in AU. This roughly corresponds to a geostationary orbit.
        earth);


    // Final integration time (relative to ephem->jd_ref)
    double tend = 7404.5;  

    // Integrate until we reach tend or an error occurs
    while (r->t < tend && r->status <= 0){
       // Generate output every 20 days
       reb_integrate(r, r->t + 20.0);

       // Output test particle positions
       printf("t = %.1f\n", r->t + ephem->jd_ref);
       for(int i=0; i<r->N; i++){
           struct reb_particle p = r->particles[i];
           printf("particles[%d]:  \tx = %.12f \ty = %.12f \tz = %.12f\n", i, p.x, p.y, p.z);
       }
    }

    // Clean up memory
    assist_free(ax);
    reb_free_simulation(r);
}

