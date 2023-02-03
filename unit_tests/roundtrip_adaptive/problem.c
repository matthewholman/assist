/**
 * A unit test integrating an asteroid forward, then backwards in time.
 */
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include "rebound.h"
#include "assist.h"

const double au2meter = 149597870700;

double roundtrip(struct assist_ephem* ephem, double trange){
    struct reb_simulation* r = reb_create_simulation();
    struct assist_extras* ax = assist_attach(r, ephem);
    double t0 = 2458849.5-ephem->jd_ref;
    double x0 = 3.3388753502614090e+00;
    double y0 = -9.1765182678903168e-01;
    double z0 = -5.0385906775843303e-01;
    double vx0 = 5.22000300435103972568e-03;
    double vy0 = 6.56760310074497024591e-03;
    double vz0 = 2.44038701006633581073e-03;
    r->t = t0; 
    // 3 copies with slightly different initial coonditions
    // Using 3 particles to achieve a smoother timestep. 
    const int N = 3;
    for (int i=0; i<N; i++){
        reb_add_fmt(r, "x y z vx vy vz", x0, y0, z0, vx0 +(double)i/1e4, vy0, vz0);
    }
   
    reb_integrate(r,  t0 + trange);
    reb_integrate(r,  t0);

    assert(r->t == t0);

    
    double d = 0;
    for (int i=0; i<N; i++){
        double dx = r->particles[i].x - x0;
        double dy = r->particles[i].y - y0;
        double dz = r->particles[i].z - z0;
        d += sqrt(dx*dx + dy*dy + dz*dz)*au2meter;
    }
   
    assist_free(ax);
    reb_free_simulation(r);
    return d/N;
}

int main(int argc, char* argv[]){

    struct assist_ephem* ephem = assist_ephem_create(
            "../../data/linux_p1550p2650.440",
            "../../data/sb441-n16.bsp");
    if (ephem == NULL){
        fprintf(stderr,"Error initializing assist_ephem.\n");
        exit(1);
    }
    {
        double r0 = roundtrip(ephem, 100);
        printf("distance: %em\n",r0);
        assert(r0 < 1e-2); // required accuracy in m
    }
    {
        double r0 = roundtrip(ephem, 1000);
        printf("distance: %em\n",r0);
        assert(r0 < 5e-2); // required accuracy in m
    }
    {
        double r0 = roundtrip(ephem, 10000);
        printf("distance: %em\n",r0);
        assert(r0 < 1e-1); // required accuracy in m
    }
    assist_ephem_free(ephem);
}

