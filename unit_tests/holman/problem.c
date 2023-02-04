/**
 * A unit test integrating the asteroid Holman for 30 days, then comparing the 
 * position and coordinates to JPL Horizon data.
    
 * Batch file to generate the JPL data (last queried on Jan 16th 2023)
 *
 *    !$$SOF
 *    MAKE_EPHEM=YES
 *    COMMAND='DES=2003666;'
 *    EPHEM_TYPE=VECTORS
 *    CENTER='500@0'
 *    START_TIME='2023-01-17'
 *    STOP_TIME='2023-02-16'
 *    STEP_SIZE='1 DAYS'
 *    VEC_TABLE='3'
 *    REF_SYSTEM='ICRF'
 *    REF_PLANE='FRAME'
 *    VEC_CORR='NONE'
 *    OUT_UNITS='AU-D'
 *    VEC_LABELS='YES'
 *    VEC_DELTA_T='NO'
 *    CSV_FORMAT='NO'
 *    OBJ_DATA='YES'
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
    r->t = 8416.5;

    // Initial conditions of asteroid Holman
    reb_add_fmt(r, "x y z vx vy vz",
        -2.724183384883979E+00, -3.523994546329214E-02, 9.036596202793466E-02, 
        -1.374545432301129E-04, -1.027075301472321E-02, -4.195690627695180E-03); 
   
    reb_integrate(r, 8446.5);
   
    // Final data from NASA Horizons
    reb_add_fmt(r, "x y z vx vy vz",
        -2.710320457933958E+00, -3.424507930535848E-01, -3.582442972611413E-02,
        1.059255302926290E-03, -1.018748422976772E-02, -4.207712906489264E-03);
    
    double au2meter = 149597870700;

    double diff_x = fabs(r->particles[0].x-r->particles[1].x)*au2meter;
    assert(diff_x < 0.1); // require that integration error is less than 10 cm after 30 day integration
    double diff_y = fabs(r->particles[0].y-r->particles[1].y)*au2meter;
    assert(diff_y < 0.1); // require that integration error is less than 10 cm after 30 day integration
    double diff_z = fabs(r->particles[0].z-r->particles[1].z)*au2meter;
    assert(diff_z < 0.1); // require that integration error is less than 10 cm after 30 day integration
        
    double diff_vx = fabs(r->particles[0].vx-r->particles[1].vx)*au2meter;
    assert(diff_vx < 5e-3); // require that integration error is less than 5mm/day after 30 day integration
    double diff_vy = fabs(r->particles[0].vy-r->particles[1].vy)*au2meter;
    assert(diff_vy < 5e-3); // require that integration error is less than 5mm/day after 30 day integration
    double diff_vz = fabs(r->particles[0].vz-r->particles[1].vz)*au2meter;
    assert(diff_vz < 5e-3); // require that integration error is less than 5mm/day after 30 day integration
        
    assist_free(ax);
    assist_ephem_free(ephem);
    reb_free_simulation(r);
}

