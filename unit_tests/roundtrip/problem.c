/**
 * A unit test integrating the asteroid Holman for 30 days, then comparing the 
 * position and coordinates to JPL Horizon data.
 */
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include "rebound.h"
#include "assist.h"

double roundtrip(struct assist_ephem* ephem, double trange){
    struct reb_simulation* r = reb_create_simulation();
    struct assist_extras* ax = assist_attach(r, ephem);

    double t0 = 2458849.5-ephem->jd_ref;
    double x0 = 3.3388753502614090e+00;
    double y0 = -9.1765182678903168e-01;
    double z0 = -5.0385906775843303e-01;

    r->t = t0; 

    // Initial conditions of asteroid Holman
    reb_add_fmt(r, "x y z vx vy vz",
        x0, y0, z0,
        2.8056633153049852e-03,  7.5504086883996860e-03,  2.9800282074358684e-03);
   
    reb_integrate(r,  t0 + trange);
    assert(r->t == t0+trange);
    r->dt = 0.01; 
    reb_integrate(r,  t0);
    assert(r->t == t0);

    double dx = r->particles[0].x - x0;
    double dy = r->particles[0].y - y0;
    double dz = r->particles[0].z - z0;

    double au2meter = 149597870700;
    double d = sqrt(dx*dx + dy*dy + dz*dz)*au2meter;
   
    assist_free(ax);
    reb_free_simulation(r);
    return d;
}

int main(int argc, char* argv[]){

    struct assist_ephem* ephem = assist_ephem_init(
            "../../data/linux_p1550p2650.440",
            "../../data/sb441-n16.bsp");
    if (ephem == NULL){
        fprintf(stderr,"Error initializing assist_ephem.\n");
        exit(1);
    }
   
  //  {
  //      double r0 = roundtrip(ephem, 100);
  //      printf("distance: %fm\n",r0);
  //      assert(r0 < 0.002); // required accuracy in m
  //  }
  //  {
  //      double r0 = roundtrip(ephem, 1000);
  //      printf("distance: %fm\n",r0);
  //      assert(r0 < 0.01); // required accuracy in m
  //  }
  //  {
  //      double r0 = roundtrip(ephem, 10000);
  //      printf("distance: %fm\n",r0);
  //      assert(r0 < 0.05); // required accuracy in m
  //  }
    for (double trange = 10000; trange<1e5; trange*=1.05)
    {
        double r0 = roundtrip(ephem, trange);
        printf("trange = %f   error: %fm\n",trange, r0);
        //assert(r0 < 0.2); // required accuracy in m
    }
        
     
    assist_ephem_free(ephem);
}

