/**
 * A unit test checking if an ASSIST supported simulation can be converted to a pure REBOUND simuation. 
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
    r->t = 2458849.5-ephem->jd_ref;

    // Initial conditions of asteroid Holman
    reb_simulation_add_fmt(r, "x y z vx vy vz",
        3.3388753502614090e+00, -9.1765182678903168e-01, -5.0385906775843303e-01,
        2.8056633153049852e-03,  7.5504086883996860e-03,  2.9800282074358684e-03);
   
    // Integrate for some time.
    reb_simulation_integrate(r, r->t +100);

    // Create pure rebound simulation.
    struct reb_simulation* r2 = assist_simulation_convert_to_rebound(r, ephem, 0);
    
    // Integrate both simulations for some more time.
    reb_simulation_integrate(r, r->t +100);
    reb_simulation_integrate(r2, r2->t +100);

    double dx = r->particles[0].x - r2->particles[r2->N-1].x;
    double dy = r->particles[0].y - r2->particles[r2->N-1].y;
    double dz = r->particles[0].z - r2->particles[r2->N-1].z;

    double d = sqrt(dx*dx + dy*dy + dz*dz)*au2meter;
    assert(d < 200); // required accuracy in m
   
    assist_free(ax);
    reb_simulation_free(r);
    reb_simulation_free(r2);
    assist_ephem_free(ephem);
}

