/**
 * Tests the ephem_cache
 */
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include "rebound.h"
#include "assist.h"

struct reb_particle integrate(struct assist_ephem* ephem, int cache_on, int direction){
    struct reb_simulation* r = reb_simulation_create();
    struct assist_extras* ax = assist_attach(r, ephem);
    if (cache_on == 0){
        ax->ephem_cache = NULL;
    }

    double t0 = 2458849.5-ephem->jd_ref;
    r->t = t0; 

    // Initial conditions of asteroid Holman
    reb_simulation_add_fmt(r, "x y z vx vy vz",
        3.3388753502614090e+00, -9.1765182678903168e-01, -5.0385906775843303e-01,
        2.8056633153049852e-03,  7.5504086883996860e-03,  2.9800282074358684e-03);
   
    reb_simulation_integrate(r,  t0 + direction*1000);
    assert(r->t == t0+direction*1000);

    struct reb_particle p = r->particles[0];
    assist_free(ax);
    reb_simulation_free(r);
    return p;
}

int main(int argc, char* argv[]){

    struct assist_ephem* ephem = assist_ephem_create(
            "../../data/de440.bsp",
            "../../data/sb441-n16.bsp");
    if (ephem == NULL){
        fprintf(stderr,"Error initializing assist_ephem.\n");
        exit(1);
    }
   
    // Forward
    {
        struct reb_particle p0 = integrate(ephem, 0, 1); // Run without cache
        struct reb_particle p1 = integrate(ephem, 1, 1); // Run with cache
        // Check that ephem_cache does not affect result
        assert(p0.x == p1.x);
        assert(p0.y == p1.y);
        assert(p0.z == p1.z);
        assert(p0.vx == p1.vx);
        assert(p0.vy == p1.vy);
        assert(p0.vz == p1.vz);
        printf("Check passed.\n");
    }
    // Backward
    {
        struct reb_particle p0 = integrate(ephem, 0, -1); // Run without cache
        struct reb_particle p1 = integrate(ephem, 1, -1); // Run with cache
        // Check that ephem_cache does not affect result
        assert(p0.x == p1.x);
        assert(p0.y == p1.y);
        assert(p0.z == p1.z);
        assert(p0.vx == p1.vx);
        assert(p0.vy == p1.vy);
        assert(p0.vz == p1.vz);
        printf("Check passed.\n");
    }
     
    assist_ephem_free(ephem);
}

