/**
 * Post-Newtonian correction from general relativity
 * 
 * This example shows how to add post-newtonian corrections to REBOUND simulations with reboundx.
 * If you have GLUT installed for the visualization, press 'w' and/or 'c' for a clearer view of
 * the whole orbit.
 */
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include "rebound.h"
#include "assist.h"

int main(int argc, char* argv[]){
    struct reb_simulation* sim = reb_create_simulation();
    
    // Central object at origin with mass=1
    reb_add_fmt(sim, "m", 1.); 
    
    // Massless particle on circular orbit with a=1
    reb_add_fmt(sim, "a", 1.); 

    struct assist_extras* assist = assist_attach(sim);
    
    reb_integrate(sim, 100); 

    assist_free(assist);    // this explicitly frees all the memory allocated by ASSIST

    reb_free_simulation(sim);
}
