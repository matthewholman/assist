/**
 * A unit test fir tge 5303-Ceres encounter.
 */
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
    
    struct reb_simulation* r = reb_create_simulation();
    struct assist_extras* ax = assist_attach(r, ephem);
    r->t = 2449718.50 - ephem->jd_ref;
    
    double au2meter = 149597870700;

    // Decide whether to include Harmonics in this calculation
    // ax->forces ^= ASSIST_FORCE_EARTH_HARMONICS;
    // ax->forces ^= ASSIST_FORCE_SUN_HARMONICS;

    // Initial conditions of asteroid (5303) Parijskij
    reb_add_fmt(r, "x y z vx vy vz",
        -2.232847879711731E+00, 1.574146331186095E+00, 8.329414259670296E-01,
        -6.247432571575564E-03, -7.431073424167182E-03, -3.231725223736132E-03);
   
    reb_integrate(r, 2453371.5 - ephem->jd_ref);
   
    // Final data from NASA Horizons
    reb_add_fmt(r, "x y z vx vy vz",
            -2.656522667009432E+00, 8.168437454347069E-01, 4.947270505430544E-01,
            -3.333217972311964E-03, -8.880226801086633E-03, -4.036456444328579E-03);

    // Final data from JPL Small body code
    struct reb_particle sun = assist_get_particle(ephem, 0, r->t);
    reb_add_fmt(r, "x y z",
        -398052983.0882521/au2meter*1e3+sun.x, 
        122233037.4686653/au2meter*1e3+sun.y, 
        74042102.86177272/au2meter*1e3+sun.z);
    
    {
        double diff_x = fabs(r->particles[0].x-r->particles[1].x)*au2meter;
        double diff_y = fabs(r->particles[0].y-r->particles[1].y)*au2meter;
        double diff_z = fabs(r->particles[0].z-r->particles[1].z)*au2meter;
        double diff = sqrt(diff_x*diff_x + diff_y*diff_y + diff_z*diff_z);
        double tolerance = 501;// meter
        printf("diff to Horizons %fm\n",diff);
        assert(diff < tolerance); 
    }
    {
        double diff_x = fabs(r->particles[0].x-r->particles[2].x)*au2meter;
        double diff_y = fabs(r->particles[0].y-r->particles[2].y)*au2meter;
        double diff_z = fabs(r->particles[0].z-r->particles[2].z)*au2meter;
        double diff = sqrt(diff_x*diff_x + diff_y*diff_y + diff_z*diff_z);
        double tolerance = 500;// meter
        printf("diff to JPL Small Body code %fm\n",diff);
        assert(diff < tolerance); 
    }
        
    assist_free(ax);
    assist_ephem_free(ephem);
    reb_free_simulation(r);
}

