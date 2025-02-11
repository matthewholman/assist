/**
 * A unit test integrating the asteroid Apophis.
 * The integration includes non-gravitational forces.
 */
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include "rebound.h"
#include "assist.h"

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

    // Turn on GR from all planets
    ax->gr_eih_sources = ASSIST_BODY_NPLANETS;

    double t_initial = 2.4621385359989386E+06 - ephem->jd_ref;
    r->t = t_initial;

    // Initial data from JPL small body code
    struct reb_particle p_initial = {
                    .x = -5.5946538550488512E-01, 
                    .y =  8.5647564757574512E-01, 
                    .z =  3.0415066217102493E-01,
                    .vx= -1.3818324735921638E-02, 
                    .vy= -6.0088275597939191E-03, 
                    .vz= -2.5805044631309632E-03};
    
    // Move from heliocentric to barycentric frame
    struct reb_particle sun_initial = assist_get_particle(ephem, 0, t_initial);
    reb_particle_iadd(&p_initial, &sun_initial);

    // Add particle to simulation
    reb_simulation_add(r, p_initial);

    // Define non-gravitational parameters A1, A2, A3
    double params[] = {4.999999873689E-13, -2.901085508711E-14, 0.0}; 
    ax->particle_params = params;
    
    // Do the actual integration   
    r->ri_ias15.min_dt = 0.001; // in days (prevent timestep getting very small during encounter)
    double t_final = 2.4625030372426095E+06 -ephem->jd_ref; // 1 year later
    reb_simulation_integrate(r, t_final);

    // Final data from JPL small body code
    struct reb_particle p_final = {
                    .x =  1.7028330901729331E-02,
                    .y =  1.2193934090901304E+00,
                    .z =  4.7823589236374386E-01,
                    .vx= -1.3536187639388663E-02,
                    .vy=  5.3200999989786943E-04,
                    .vz= -1.6648346717629861E-05};

    // Move from heliocentric to barycentric frame
    struct reb_particle sun_final = assist_get_particle(ephem, 0, t_final);
    reb_particle_iadd(&p_final, &sun_final);
    
    // Calculate difference
    double au2meter = 149597870700;
    double diff = reb_particle_distance(&p_final, &r->particles[0]) * au2meter; // in meter

    // Ensure accuracy is better than 250m
    fprintf(stderr,"diff to JPL Small Body code %fm\n",diff);
    assert(diff < 200); 

    assist_free(ax);
    assist_ephem_free(ephem);
    reb_simulation_free(r);
}

