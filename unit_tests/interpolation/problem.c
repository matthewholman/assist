/**
 * A unit test integrating the asteroid Holman for 30 days, then comparing the 
 * position and coordinates to JPL Horizon data.
 */
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include "rebound.h"
#include "assist.h"

int main(int argc, char* argv[]){
    // Cleanup previous file
    remove("out.bin"); 

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

    reb_simulationarchive_automate_step(r, "out.bin", 1);
   
    // Integrate past required output time
    reb_integrate(r, 9000.0);
    
    // Cleanup rebound and assist structures
    assist_free(ax);
    assist_ephem_free(ephem);
    reb_free_simulation(r);


    // Create an interpolated simulation from the archive 
    struct reb_simulationarchive* sa = reb_open_simulationarchive("out.bin");
    double tend = 8446.5;
    struct reb_simulation* r_interpolated = assist_create_interpolated_simulation(sa, tend);
    reb_close_simulationarchive(sa);
    
    // Get particle position from simulation
    struct reb_particle p = r_interpolated->particles[0];
    reb_free_simulation(r_interpolated);
   
    // Final data from NASA Horizons
    struct reb_particle ph = {.x= -2.710320457933958E+00, 
                              .y= -3.424507930535848E-01,
                              .z= -3.582442972611413E-02,
                              .vx= 1.059255302926290E-03, 
                              .vy= -1.018748422976772E-02, 
                              .vz= -4.207712906489264E-03};
    
    double au2meter = 149597870700;

    double diff_x = fabs(ph.x-p.x)*au2meter;
    printf("diff %.20f\n",diff_x);
    assert(diff_x < 0.1); // require that integration error is less than 10 cm after 30 day integration
    double diff_y = fabs(ph.y-p.y)*au2meter;
    printf("diff %.20f\n",diff_y);
    assert(diff_y < 0.1); // require that integration error is less than 10 cm after 30 day integration
    double diff_z = fabs(ph.z-p.z)*au2meter;
    printf("diff %.20f\n",diff_z);
    assert(diff_z < 0.1); // require that integration error is less than 10 cm after 30 day integration
        
    double diff_vx = fabs(ph.vx-p.vx)*au2meter;
    assert(diff_vx < 5e-3); // require that integration error is less than 5mm/day after 30 day integration
    double diff_vy = fabs(ph.vy-p.vy)*au2meter;
    assert(diff_vy < 5e-3); // require that integration error is less than 5mm/day after 30 day integration
    double diff_vz = fabs(ph.vz-p.vz)*au2meter;
    assert(diff_vz < 5e-3); // require that integration error is less than 5mm/day after 30 day integration
        
}

