/**
 * Querying ephemeris data
 *
 * This example demonstrates how to use ASSIST to query the JPL DE440/441
 * ephemeris data. 
 */
#include <stdio.h>
#include <stdlib.h>
#include "rebound.h"
#include "assist.h"

int main(int argc, char* argv[]){
    // We first load the ephemeris structure
    // You might need to adjust the paths for the datafiles
    // if they are stored in a different location.
    // See the installation instruction if you still need to 
    // download the files.

    struct assist_ephem* ephem = assist_ephem_create(
            "../../data/de440.bsp",
            "../../data/sb441-n16.bsp");
    if (!ephem){
        printf("Cannot create ephemeris structure.\n");
        exit(-1);
    }
    
    // All the ephemeris queries are relative to this reference Julian Date.
    // The value below is the default. So this line has not effect.
    // You might want to change the reference Julian Date if you are looking
    // for high relative accuracy for different query times.
    ephem->jd_ref = 2451545.0; 
    
    // We can now query the position of planets and the 16 most important asteroids. 
    // Here are a few examples.   
    double t = 0; // time in Julian Days relative to ephem->jd_ref
    struct reb_particle sun = assist_get_particle(ephem, ASSIST_BODY_SUN, t);
    struct reb_particle earth = assist_get_particle(ephem, ASSIST_BODY_EARTH, t);
    struct reb_particle moon = assist_get_particle(ephem, ASSIST_BODY_MOON, t);
    
    // The return values are a reb_particle structure. 
    // The following prints out the position of the earth at 2451545.0 Julian Days
    printf("EARTH:  x=%f AU   y=%f AU   z=%f AU\n", earth.x, earth.y, earth.z); // Units are AU
    
    // Similarly, the velocity for the moon at the same time is:
    printf("MOON:   vx=%f AU/day   vy=%f AU/day   vz=%f AU/day\n", moon.vx, moon.vy, moon.vz); // Units are AU/day
                                                                                       
    // The m attribute of the partiles is not the mass itself, 
    // but the mass times the gravitational constant G.
    // We work in units of solarmasses, astronomical units, and days.
    // Thus G has the value of 0.00029591 AU^3/(day^2 sunmass).
    printf("SUN:    m=%f \n", sun.m); // Units are AU/day
    
    // Let's clean up the ephemeris structure to avoid memory leaks.
    assist_ephem_free(ephem);
}

