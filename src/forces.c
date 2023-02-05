/**
 * @file    forces.c
 * @brief   Functions implementing the various forces used in the integration
 * @author  Matthew Holman, Arya Akmal, Hanno Rein
 * 
 * @section     LICENSE
 * Copyright (c) 2022 Matthew Holman, Arya Akmal, Hanno Rein
 *
 * This file is part of ASSIST.
 *
 * ASSIST is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * ASSIST is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with ASSIST. If not, see <http://www.gnu.org/licenses/>.
 *
 */

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <stdint.h>
#include <math.h>
#include <limits.h>
#include <float.h>
#include "assist.h"
#include "rebound.h"
#include "const.h"
#include "spk.h"
#include "planets.h"
#include "forces.h"


// Forward function declarations
static void assist_additional_force_direct(struct reb_simulation* sim, double xo, double yo, double zo, FILE *outfile);
static void assist_additional_force_solar_J2(struct reb_simulation* sim, double xo, double yo, double zo, FILE *outfile);
static void assist_additional_force_earth_J2J4(struct reb_simulation* sim, double xo, double yo, double zo, FILE *outfile);
static void assist_additional_force_non_gravitational(struct reb_simulation* sim, double xo, double yo, double zo, double vxo, double vyo, double vzo, FILE *outfile);
static void assist_additional_force_potential_GR(struct reb_simulation* sim, double xo, double yo, double zo, double vxo, double vyo, double vzo, FILE *outfile);
static void assist_additional_force_simple_GR(struct reb_simulation* sim, double xo, double yo, double zo, double vxo, double vyo, double vzo, FILE *outfile);
static void assist_additional_force_eih_GR(struct reb_simulation* sim, int eih_loop_limit, double xo, double yo, double zo, double vxo, double vyo, double vzo, double axo, double ayo, double azo,	FILE *outfile, FILE *eih_file);

void assist_additional_forces(struct reb_simulation* sim){

    // implement additional_forces here

    struct assist_extras* assist = (struct assist_extras*) sim->extras;
    struct assist_ephem* ephem = assist->ephem;
    int geo = assist->geocentric;

    const unsigned int N = sim->N;  // N includes real+variational particles
    const unsigned int N_real = N - sim->N_var;    

    // The limit of the EIH GR limit should be a free
    // parameter
    int eih_loop_limit = ASSIST_BODY_NPLANETS; // 1;

    const double t = sim->t;

    double GM;

    double xo, yo, zo, vxo, vyo, vzo, axo, ayo, azo;

    // Check which center is used.
    // The current options are the barycenter (default) and geocenter.
    // We might consider adding the heliocenter.

    if(geo == 1){
	// geocentric
	// Get mass, position, velocity, and acceleration of the Earth for later use.
	// The offset position is used to adjust the particle positions.
	int flag = assist_all_ephem(ephem, assist->ephem_cache, ASSIST_BODY_EARTH, t, &GM, &xo, &yo, &zo, &vxo, &vyo, &vzo, &axo, &ayo, &azo);
	if(flag != ASSIST_SUCCESS){
	    char outstring[50];
	    sprintf(outstring, "%s %d %d\n", "Ephemeris error a ", 3, flag);
	    reb_error(sim, outstring);
	}
    }else{
	// barycentric
	xo = 0.0;  yo = 0.0;  zo = 0.0;
	vxo = 0.0; vyo = 0.0; vzo = 0.0;
	axo = 0.0; ayo = 0.0; azo = 0.0;	
    }

    // TODO: eliminate the output files after testing
    // or make this more flexible
    FILE *outfile = NULL;
    // Uncomment these lines and recompile for testing.    
    //outfile = fopen("acc.out", "a+");
    //FILE *vfile = NULL;    
    //vfile = fopen("vary_acc.out", "a+");

    // These should be executed in order from smallest
    // to largest

    // Pick one of the three GR routines
    //assist_additional_force_potential_GR(sim, xo, yo, zo, vxo, vyo, vzo, outfile);    
    //sim->force_is_velocity_dependent = 1;    
    //assist_additional_force_simple_GR(sim, xo, yo, zo, vxo, vyo, vzo, outfile);

    /*
    FILE *eih_file = NULL;
    // Uncomment this line and recompile for testing.
    //eih_file = fopen("eih_acc.out", "w");

    assist_additional_force_direct(sim, xo, yo, zo, outfile);
    assist_additional_force_earth_J2J4(sim, xo, yo, zo, outfile);
    assist_additional_force_solar_J2(sim, xo, yo, zo, outfile);        
    assist_additional_force_non_gravitational(sim, xo, yo, zo, vxo, vyo, vzo, outfile);    
    //assist_additional_force_simple_GR(sim, xo, yo, zo, vxo, vyo, vzo, outfile);    
    assist_additional_force_eih_GR(sim, eih_loop_limit,
	   xo, yo, zo, vxo, vyo, vzo, axo, ayo, azo,	   
	   outfile, eih_file);
    */

    if (assist->forces & ASSIST_FORCE_NON_GRAVITATIONAL){
        assist_additional_force_non_gravitational(sim, xo, yo, zo, vxo, vyo, vzo, outfile);
    }
    if (assist->forces & ASSIST_FORCE_EARTH_HARMONICS){
        assist_additional_force_earth_J2J4(sim, xo, yo, zo, outfile);
    }
    if (assist->forces & ASSIST_FORCE_SUN_HARMONICS){
        assist_additional_force_solar_J2(sim, xo, yo, zo, outfile);
    }
    
    FILE *eih_file = NULL;
    // Uncomment this line and recompile for testing.
    //eih_file = fopen("eih_acc.out", "w");

    if (assist->forces & ASSIST_FORCE_GR_EIH){
        assist_additional_force_eih_GR(sim, eih_loop_limit,
           xo, yo, zo, vxo, vyo, vzo, axo, ayo, azo,
           outfile, eih_file);
    }

    if (assist->forces & ASSIST_FORCE_GR_POTENTIAL){
        assist_additional_force_potential_GR(sim, xo, yo, zo, vxo, vyo, vzo, outfile);
    }
    if (assist->forces & ASSIST_FORCE_GR_SIMPLE){
        assist_additional_force_simple_GR(sim, xo, yo, zo, vxo, vyo, vzo, outfile);
    }
    if (assist->forces & (ASSIST_FORCE_SUN | ASSIST_FORCE_PLANETS | ASSIST_FORCE_ASTEROIDS)){
        assist_additional_force_direct(sim, xo, yo, zo, outfile);
    }
    
    // Uncomment one of these lines and recompile for testing.
    //test_vary(sim, vfile);
    //test_vary_2nd(sim, vfile);    

    //fclose(eih_file);
    //fflush(outfile);
    //fclose(outfile);

    if(geo == 1){
	// geocentric
	// TODO: This part will need work for the variational equations
	// to work properly.
	assist_all_ephem(ephem, assist->ephem_cache, ASSIST_BODY_EARTH, t, &GM, &xo, &yo, &zo, &vxo, &vyo, &vzo, &axo, &ayo, &azo);

	// This is the indirect term for geocentric equations
	// of motion.
	for (int j=0; j<N_real; j++){    

	    sim->particles[j].ax -= axo;
	    sim->particles[j].ay -= ayo;
	    sim->particles[j].az -= azo;

	}
    }
}

int assist_all_ephem(struct assist_ephem* ephem, struct assist_ephem_cache* ephem_cache, const int i, const double t, double* const GM,
		      double* const x, double* const y, double* const z,
		      double* const vx, double* const vy, double* const vz,
		      double* const ax, double* const ay, double* const az
        ){
    if (ephem_cache){
        double* ct = ephem_cache->t+7*i;
        for (int s=0; s<7; s+=1){
            if (t==ct[s]){
                struct assist_cache_item* items = ephem_cache->items+i*7+s;
                *GM = items->GM;
                *x  = items->x;
                *y  = items->y;
                *z  = items->z;
                *vx = items->vx;
                *vy = items->vy;
                *vz = items->vz;
                *ax = items->ax;
                *ay = items->ay;
                *az = items->az;
                return ASSIST_SUCCESS;
            }
        }
        // No match
    }

    const double jd_ref = ephem->jd_ref;

    // Get position and mass of massive body i.
    if(i < ASSIST_BODY_NPLANETS){
        int flag = assist_jpl_calc(ephem->jpl, jd_ref, t, i,  GM, x, y, z, vx, vy, vz, ax, ay, az);
        if(flag != ASSIST_SUCCESS) return(flag);
    }else{
        // Get position and mass of asteroid i-ASSIST_BODY_NPLANETS.
        int flag = assist_spk_calc(ephem->spl, jd_ref, t, i-ASSIST_BODY_NPLANETS, GM, x, y, z);
        if(flag != ASSIST_SUCCESS) return(flag);

        double GMs, xs, ys, zs;
        double vxs, vys, vzs, axs, ays, azs; // Not needed
        flag = assist_all_ephem(ephem, ephem_cache, ASSIST_BODY_SUN, t, &GMs, &xs, &ys, &zs, &vxs, &vys, &vzs, &axs, &ays, &azs);
        if(flag != ASSIST_SUCCESS) return(flag);		    

        // Translate massive asteroids from heliocentric to barycentric.
        *x += xs; *y += ys; *z += zs;
        // velocities and accelerations are not needed for these
        // bodies
        *vx = NAN; *vy = NAN; *vz = NAN;
        *ax = NAN; *ay = NAN; *az = NAN;

    }

    if (ephem_cache){
        double* ct = ephem_cache->t+7*i;
        // Find oldest
        int os = 0;
        double ot = ct[0];
        for (int s=1; s<7; s+=1){
            if (ct[s]<ot){
                ot = ct[s];
                os = s;
            }
        }
        // repolace oldest
        ct[os] = t;
        struct assist_cache_item* items = ephem_cache->items + i*7 +os;
        items->GM = *GM;
        items->x = *x;
        items->y = *y;
        items->z = *z;
        items->vx = *vx;
        items->vy = *vy;
        items->vz = *vz;
        items->ax = *ax;
        items->ay = *ay;
        items->az = *az;
    }
    return(ASSIST_SUCCESS);
}


static void assist_additional_force_direct(struct reb_simulation* sim, double xo, double yo, double zo, FILE *outfile){
    const unsigned int N = sim->N;  // N includes real+variational particles
    const unsigned int N_real = N - sim->N_var;

    struct assist_extras* assist = (struct assist_extras*) sim->extras;
    struct assist_ephem* ephem = assist->ephem;
    const double jd_ref = ephem->jd_ref;
    
    const double t = sim->t;    

    struct reb_particle* const particles = sim->particles;

    double GM;
    double x, y, z, vx, vy, vz, ax, ay, az;

    static const int order[ASSIST_BODY_NPLANETS + ASSIST_BODY_NASTEROIDS] = { 
        ASSIST_BODY_CYBELE,
        ASSIST_BODY_EUPHROSYNE,
        ASSIST_BODY_IRIS,
        ASSIST_BODY_THISBE,
        ASSIST_BODY_CAMILLA,
        ASSIST_BODY_PSYCHE,
        ASSIST_BODY_JUNO,
        ASSIST_BODY_EUNOMIA,
        ASSIST_BODY_SYLVIA,
        ASSIST_BODY_EUROPA,
        ASSIST_BODY_INTERAMNIA,
        ASSIST_BODY_DAVIDA,
        ASSIST_BODY_HYGIEA,
        ASSIST_BODY_PALLAS,
        ASSIST_BODY_VESTA,
        ASSIST_BODY_CERES,
        ASSIST_BODY_PLUTO,
        ASSIST_BODY_MOON,
        ASSIST_BODY_MARS,
        ASSIST_BODY_MERCURY,
        ASSIST_BODY_NEPTUNE,
        ASSIST_BODY_URANUS,
        ASSIST_BODY_EARTH,
        ASSIST_BODY_VENUS,
        ASSIST_BODY_SATURN,
        ASSIST_BODY_JUPITER,
        ASSIST_BODY_SUN
    };

    // Direct forces from massives bodies
    for (int k=0; k < ASSIST_BODY_NPLANETS + ASSIST_BODY_NASTEROIDS; k++){
        int i = order[k];
        if (i==ASSIST_BODY_SUN && !(assist->forces & ASSIST_FORCE_SUN)) continue;
        if (i>ASSIST_BODY_SUN && i<ASSIST_BODY_NPLANETS && !(assist->forces & ASSIST_FORCE_PLANETS)) continue;
        if (i>=ASSIST_BODY_NPLANETS && !(assist->forces & ASSIST_FORCE_ASTEROIDS)) continue;

        // Get position and mass of massive body i.
        // TOOD: make a version that returns the positions, velocities,
        // and accelerations for all the bodies at a given time.

        int flag = assist_all_ephem(ephem, assist->ephem_cache, i, t, &GM, &x, &y, &z, &vx, &vy, &vz, &ax, &ay, &az);

        if(flag != ASSIST_SUCCESS){
            char outstring[50];
            sprintf(outstring, "%s %d %d\n", "Ephemeris error b ", i, flag);	    
            reb_error(sim, outstring);
        }

        // Loop over test particles
        for (int j=0; j<N_real; j++){

            // Compute position vector of test particle j relative to massive body i.
            const double dx = particles[j].x + (xo - x); 
            const double dy = particles[j].y + (yo - y);
            const double dz = particles[j].z + (zo - z);
            const double r2 = dx*dx + dy*dy + dz*dz;
            const double _r  = sqrt(r2);
            const double prefac = GM/(_r*_r*_r);

            if(outfile){
                fprintf(outfile, "%3d %25.16le %25.16le %25.16le %25.16le %25.16le %25.16le %25.16le %25.16le\n", i, jd_ref+t, GM, dx, dy, dz, -prefac*dx, -prefac*dy, -prefac*dz);
		fflush(outfile);		
            }

            particles[j].ax -= prefac*dx;
            particles[j].ay -= prefac*dy;
            particles[j].az -= prefac*dz;

        }
    }

    // Acceleration of variational particles due to direct forces from massive bodies 
    // Loop over the perturbers
    // We should put a check at the top to see if there are any variational
    // particles.

    for (int k=0; k < ASSIST_BODY_NPLANETS + ASSIST_BODY_NASTEROIDS; k++){
	int i = order[k];
	//int i = k;

        // Get position and mass of massive body i.	
	int flag = assist_all_ephem(ephem, assist->ephem_cache, i, t, &GM, &x, &y, &z, &vx, &vy, &vz, &ax, &ay, &az);

	if(flag != ASSIST_SUCCESS){
	    char outstring[50];
	    sprintf(outstring, "%s %d %d\n", "Ephemeris error c ", i, flag);	    	    
	    reb_error(sim, outstring);
	}

    // Skip remainder of calculation if variational particles are not used
    if (sim->var_config_N == 0) return;

	// Loop over test particles
    for (int j=0; j<N_real; j++){

	    // This stuff was already computed above.
	    // We can recycle it.
	    //
	    // We could also skip real particles that
	    // have no corresponding variational particles. 
	    // For that we need an array, indexed on the
	    // real particles, that stores, for example, the
	    // number of associated variational particles.
	    //
	    const double dx = particles[j].x + (xo - x);
	    const double dy = particles[j].y + (yo - y);
	    const double dz = particles[j].z + (zo - z);
	    const double r2 = dx*dx + dy*dy + dz*dz;
	    const double _r  = sqrt(r2);
	    const double r3inv = 1./(r2*_r);
	    const double r5inv = 3.*r3inv/r2;

	    // Coefficients for variational equations
	    const double dxdx = dx*dx*r5inv - r3inv;
	    const double dydy = dy*dy*r5inv - r3inv;
	    const double dzdz = dz*dz*r5inv - r3inv;
	    const double dxdy = dx*dy*r5inv;
	    const double dxdz = dx*dz*r5inv;
	    const double dydz = dy*dz*r5inv;

	    // Loop over variational particles
	    // Update the accelerations for the variational
	    // particles that are associated with current
	    // real particle.  

	    for (int v=0; v < sim->var_config_N; v++){
		struct reb_variational_configuration const vc = sim->var_config[v];
		int tp = vc.testparticle;
		struct reb_particle* const particles_var1 = particles + vc.index;		
		if(tp == j){
	    
		    // Variational particle coords
		    const double ddx = particles_var1[0].x;
		    const double ddy = particles_var1[0].y;
		    const double ddz = particles_var1[0].z;

		    // Matrix multiplication
		    const double dax =   ddx * dxdx + ddy * dxdy + ddz * dxdz;
		    const double day =   ddx * dxdy + ddy * dydy + ddz * dydz;
		    const double daz =   ddx * dxdz + ddy * dydz + ddz * dzdz;

		    // No variational mass contributions for test particles!

		    // Accumulate acceleration terms
		    particles_var1[0].ax += GM * dax; 
		    particles_var1[0].ay += GM * day; 
		    particles_var1[0].az += GM * daz; 

		}
	    }
        }
    }
}

static void assist_additional_force_earth_J2J4(struct reb_simulation* sim, double xo, double yo, double zo, FILE *outfile){

    struct assist_extras* assist = (struct assist_extras*) sim->extras;
    struct assist_ephem* ephem = assist->ephem;
    const double jd_ref = ephem->jd_ref;
    
    //const double G = sim->G;
    const unsigned int N = sim->N;  // N includes real+variational particles
    const unsigned int N_real = N - sim->N_var;

    const double t = sim->t;    

    struct reb_particle* const particles = sim->particles;

    double GM;

    // 
    // Here is the treatment of the Earth's J2, J3, and J4.
    // J2 and J4 are based in part on gravitational_harmonics
    // example in reboundx.
    // Assumes the coordinates are geocentric.
    // Also assuming that Earth's pole is along the z
    // axis.  This is only precisely true at the J2000
    // epoch.

    // The geocenter is the reference for the Earth J2/J4 calculations.
    double xe, ye, ze, vxe, vye, vze, axe, aye, aze;    
    assist_all_ephem(ephem, assist->ephem_cache, ASSIST_BODY_EARTH, t, &GM, &xe, &ye, &ze, &vxe, &vye, &vze, &axe, &aye, &aze);
    const double GMearth = GM;

    double xr, yr, zr; //, vxr, vyr, vzr, axr, ayr, azr;
    xr = xe;  yr = ye;  zr = ze;

    const double J2e = JPL_EPHEM_J2E;
    const double J3e = JPL_EPHEM_J3E;
    const double J4e = JPL_EPHEM_J4E;
    const double au = JPL_EPHEM_CAU;
    const double Re_eq = JPL_EPHEM_RE/au;

    // Unit vector to equatorial pole at the epoch
    // Note also that the pole orientation is not changing during
    // the integration.

    double RAe =  359.87123273*M_PI/180.;
    double Dece =  89.88809752*M_PI/180.;

    // Reverting to J2000 equatorial plane
    // as is current done in JPL Horizons.
    RAe =  0.0*M_PI/180.;
    Dece =  90.0*M_PI/180.;

    double cosa = cos(RAe);
    double sina = sin(RAe);
    double cosd = cos(Dece);
    double sind = sin(Dece);
    
    // TODO: Rearrange this loop for efficiency
    // Loop over test particles        
    for (int j=0; j<N_real; j++){

        const struct reb_particle p = particles[j];
        double dx = p.x + (xo - xr);
        double dy = p.y + (yo - yr);
        double dz = p.z + (zo - zr);

        const double r2 = dx*dx + dy*dy + dz*dz;
        const double r = sqrt(r2);
	
	// Rotate to Earth equatorial frame
	double dxp =  - dx*sina      + dy*cosa;
	double dyp =  - dx*cosa*sind - dy*sina*sind + dz*cosd;
	double dzp =    dx*cosa*cosd + dy*sina*cosd + dz*sind;

	dx =  dxp;
	dy =  dyp;
	dz =  dzp;
	
	// Calculate acceleration in
	// Earth equatorial frame	

	// J2 terms
        const double costheta2 = dz*dz/r2;
        const double J2e_prefac = 3.*J2e*Re_eq*Re_eq/r2/r2/r/2.;
        const double J2e_fac = 5.*costheta2-1.;

	double resx = GMearth*J2e_prefac*J2e_fac*dx;
	double resy = GMearth*J2e_prefac*J2e_fac*dy;
	double resz = GMearth*J2e_prefac*(J2e_fac-2.)*dz;	

	// J3 terms
        const double J3e_prefac = 5.*J3e*Re_eq*Re_eq*Re_eq/r2/r2/r/2.;
        const double J3e_fac = 3.-7.*costheta2;
	
	resx += -GMearth*J3e_prefac*(1./r2)*J3e_fac*dx*dz;
        resy += -GMearth*J3e_prefac*(1./r2)*J3e_fac*dy*dz;
	resz += -GMearth*J3e_prefac*(6.*costheta2 - 7.*costheta2*costheta2-0.6);
	
	// J4 terms
        const double J4e_prefac = 5.*J4e*Re_eq*Re_eq*Re_eq*Re_eq/r2/r2/r2/r/8.;
        const double J4e_fac = 63.*costheta2*costheta2-42.*costheta2 + 3.;

        resx += GMearth*J4e_prefac*J4e_fac*dx;
        resy += GMearth*J4e_prefac*J4e_fac*dy;
        resz += GMearth*J4e_prefac*(J4e_fac+12.-28.*costheta2)*dz;

	// Rotate back to original frame
	double resxp = - resx*sina      - resy*cosa*sind + resz*cosa*cosd;
	double resyp =   resx*cosa      - resy*sina*sind + resz*sina*cosd;
	double reszp =                  + resy*cosd      + resz*sind;

	resx =  resxp;
	resy =  resyp;
	resz =  reszp;

	if(outfile){	
	    fprintf(outfile, "%3s %25.16le %25.16le %25.16le %25.16le\n", "J24", jd_ref+t, resx, resy, resz);
	    fflush(outfile);
	}
	
	// Accumulate final acceleration terms
  	particles[j].ax += resx;
        particles[j].ay += resy; 
        particles[j].az += resz;

	// Constants for variational equations
	// J2 terms
        const double J2e_fac2 = 7.*costheta2-1.;
        const double J2e_fac3 = 35.*costheta2*costheta2-30.*costheta2+3.;

	const double dxdx = GMearth*J2e_prefac*(J2e_fac-5.*J2e_fac2*dx*dx/r2);
	const double dydy = GMearth*J2e_prefac*(J2e_fac-5.*J2e_fac2*dy*dy/r2);
	const double dzdz = GMearth*J2e_prefac*(-1.)*J2e_fac3;
	const double dxdy = GMearth*J2e_prefac*(-5.)*J2e_fac2*dx*dy/r2;
	const double dydz = GMearth*J2e_prefac*(-5.)*(J2e_fac2-2.)*dy*dz/r2;
	const double dxdz = GMearth*J2e_prefac*(-5.)*(J2e_fac2-2.)*dx*dz/r2;

	// J3 terms
        const double costheta = dz/r;
        const double J3e_fac2 = 21*(-3.*costheta2+1.)/r2;
        const double J3e_fac3 = 3*(-21.*costheta2*costheta2+14.*costheta2-1.)/r2;
        const double J3e_fac4 = (-63.*costheta2*costheta2+70.*costheta2-15.)*costheta/r;

	const double dxdxJ3 = GMearth*J3e_prefac*costheta*(J3e_fac2*dx*dx-J3e_fac)/r;
	const double dydyJ3 = GMearth*J3e_prefac*costheta*(J3e_fac2*dy*dy-J3e_fac)/r;
	const double dzdzJ3 = GMearth*J3e_prefac*J3e_fac4;
	const double dxdyJ3 = GMearth*J3e_prefac*J3e_fac2*costheta*dx*dy/r;
	const double dydzJ3 = GMearth*J3e_prefac*J3e_fac3*dy;
	const double dxdzJ3 = GMearth*J3e_prefac*J3e_fac3*dx;

	// J4 terms
        const double J4e_fac2= 33.*costheta2*costheta2-18.*costheta2 + 1.;
        const double J4e_fac3= 33.*costheta2*costheta2-30.*costheta2 + 5.;
        const double J4e_fac4= 231.*costheta2*costheta2*costheta2-315.*costheta2*costheta2+105.*costheta2 - 5.;
	
	const double dxdxJ4 = GMearth*J4e_prefac*(J4e_fac-21.*J4e_fac2*dx*dx/r2);
	const double dydyJ4 = GMearth*J4e_prefac*(J4e_fac-21.*J4e_fac2*dy*dy/r2);
	const double dzdzJ4 = GMearth*J4e_prefac*(-3.)*J4e_fac4;
	const double dxdyJ4 = GMearth*J4e_prefac*(-21.)*J4e_fac2*dx*dy/r2;
	const double dydzJ4 = GMearth*J4e_prefac*(-21.)*J4e_fac3*dy*dz/r2;
	const double dxdzJ4 = GMearth*J4e_prefac*(-21.)*J4e_fac3*dx*dz/r2;

	for (int v=0; v < sim->var_config_N; v++){
	    struct reb_variational_configuration const vc = sim->var_config[v];
	    int tp = vc.testparticle;
	    struct reb_particle* const particles_var1 = particles + vc.index;		
	    if(tp == j){
	    
		// Variational particle coords
		const double ddx = particles_var1[0].x;
		const double ddy = particles_var1[0].y;
		const double ddz = particles_var1[0].z;

		// Rotate to Earth equatorial frame
		double ddxp =  - ddx*sina      + ddy*cosa;
		double ddyp =  - ddx*cosa*sind - ddy*sina*sind + ddz*cosd;
		double ddzp =    ddx*cosa*cosd + ddy*sina*cosd + ddz*sind;

		// Matrix multiplication
		// J2 part
		double dax =   ddxp * dxdx + ddyp * dxdy + ddzp * dxdz;
		double day =   ddxp * dxdy + ddyp * dydy + ddzp * dydz;
		double daz =   ddxp * dxdz + ddyp * dydz + ddzp * dzdz;

		// J3 part		
		dax +=   ddxp * dxdxJ3 + ddyp * dxdyJ3 + ddzp * dxdzJ3;
		day +=   ddxp * dxdyJ3 + ddyp * dydyJ3 + ddzp * dydzJ3;
		daz +=   ddxp * dxdzJ3 + ddyp * dydzJ3 + ddzp * dzdzJ3;

		// J4 part		
		dax +=   ddxp * dxdxJ4 + ddyp * dxdyJ4 + ddzp * dxdzJ4;
		day +=   ddxp * dxdyJ4 + ddyp * dydyJ4 + ddzp * dydzJ4;
		daz +=   ddxp * dxdzJ4 + ddyp * dydzJ4 + ddzp * dzdzJ4;

		// Rotate back to original frame
		double daxp = - dax*sina      - day*cosa*sind + daz*cosa*cosd;
		double dayp =   dax*cosa      - day*sina*sind + daz*sina*cosd;
		double dazp =                  + day*cosd      + daz*sind;

		// Accumulate acceleration terms
		particles_var1[0].ax += daxp; 
		particles_var1[0].ay += dayp; 
		particles_var1[0].az += dazp;
	    
	    }
        }
    }
}

static void assist_additional_force_solar_J2(struct reb_simulation* sim, double xo, double yo, double zo, FILE *outfile){

    struct assist_extras* assist = (struct assist_extras*) sim->extras;
    struct assist_ephem* ephem = assist->ephem;
    const double jd_ref = ephem->jd_ref;
    
    const unsigned int N = sim->N;  // N includes real+variational particles
    const unsigned int N_real = N - sim->N_var;

    const double t = sim->t;    

    struct reb_particle* const particles = sim->particles;

    double GM;

    // The Sun center is reference for these calculations.

    double xr, yr, zr, vxr, vyr, vzr, axr, ayr, azr;

    assist_all_ephem(ephem, assist->ephem_cache, ASSIST_BODY_SUN, t, &GM, &xr, &yr, &zr, &vxr, &vyr, &vzr, &axr, &ayr, &azr);
    const double GMsun = GM;    

    const double au = JPL_EPHEM_CAU;
    const double Rs_eq = JPL_EPHEM_ASUN/au;
    const double J2s = JPL_EPHEM_J2SUN;

    // Hard-coded constants.  BEWARE!
    double RAs = 286.13*M_PI/180.;
    double Decs = 63.87*M_PI/180.;

    double cosa = cos(RAs);
    double sina = sin(RAs);
    double cosd = cos(Decs);
    double sind = sin(Decs);

    //loop over test particles            
    for (int j=0; j<N_real; j++){

        const struct reb_particle p = particles[j];
        double dx = p.x + (xo - xr);
        double dy = p.y + (yo - yr);
        double dz = p.z + (zo - zr);

        const double r2 = dx*dx + dy*dy + dz*dz;
        const double r = sqrt(r2);

	// Rotate to solar equatorial frame
	double dxp =  - dx*sina      + dy*cosa;
	double dyp =  - dx*cosa*sind - dy*sina*sind + dz*cosd;
	double dzp =    dx*cosa*cosd + dy*sina*cosd + dz*sind;

	dx =  dxp;
	dy =  dyp;
	dz =  dzp;

	const double costheta2 = dz*dz/r2;
        const double J2s_prefac = 3.*J2s*Rs_eq*Rs_eq/r2/r2/r/2.;
        const double J2s_fac = 5.*costheta2-1.;
        const double J2s_fac2 = 7.*costheta2-1.;
        const double J2s_fac3 = 35.*costheta2*costheta2-30.*costheta2+3.;

	// Calculate acceleration
	double resx = GMsun*J2s_prefac*J2s_fac*dx;
	double resy = GMsun*J2s_prefac*J2s_fac*dy;
	double resz = GMsun*J2s_prefac*(J2s_fac-2.)*dz;

	// Rotate back to original frame
	double resxp = - resx*sina      - resy*cosa*sind + resz*cosa*cosd;
	double resyp =   resx*cosa      - resy*sina*sind + resz*sina*cosd;
	double reszp =                  + resy*cosd      + resz*sind;

	resx =  resxp;
	resy =  resyp;
	resz =  reszp;

	if(outfile){
	    fprintf(outfile, "%3s %25.16le %25.16le %25.16le %25.16le\n", "J2", jd_ref+t, resx, resy, resz);
	    fflush(outfile);
	}

        particles[j].ax += resx;
        particles[j].ay += resy;
        particles[j].az += resz;

	// Constants for variational equations
	// Only evaluate if there are variational particles
	const double dxdx = GMsun*J2s_prefac*(J2s_fac-5.*J2s_fac2*dx*dx/r2);
	const double dydy = GMsun*J2s_prefac*(J2s_fac-5.*J2s_fac2*dy*dy/r2);
	const double dzdz = GMsun*J2s_prefac*(-1.)*J2s_fac3;
	const double dxdy = GMsun*J2s_prefac*(-5.)*J2s_fac2*dx*dy/r2;
	const double dydz = GMsun*J2s_prefac*(-5.)*(J2s_fac2-2.)*dy*dz/r2;
	const double dxdz = GMsun*J2s_prefac*(-5.)*(J2s_fac2-2.)*dx*dz/r2;

	for (int v=0; v < sim->var_config_N; v++){
	    struct reb_variational_configuration const vc = sim->var_config[v];
	    int tp = vc.testparticle;
	    struct reb_particle* const particles_var1 = particles + vc.index;		
	    if(tp == j){
	    
		// Variational particle coords
		double ddx = particles_var1[0].x;
		double ddy = particles_var1[0].y;
		double ddz = particles_var1[0].z;
		
		// Rotate to solar equatorial frame
		double ddxp =  - ddx*sina      + ddy*cosa;
		double ddyp =  - ddx*cosa*sind - ddy*sina*sind + ddz*cosd;
		double ddzp =    ddx*cosa*cosd + ddy*sina*cosd + ddz*sind;

		ddx =  ddxp;
		ddy =  ddyp;
		ddz =  ddzp;
	    
		double daxp =   ddx * dxdx + ddy * dxdy + ddz * dxdz;
		double dayp =   ddx * dxdy + ddy * dydy + ddz * dydz;
		double dazp =   ddx * dxdz + ddy * dydz + ddz * dzdz;

		// Rotate back to original frame
		double dax = - daxp*sina      - dayp*cosa*sind + dazp*cosa*cosd;
		double day =   daxp*cosa      - dayp*sina*sind + dazp*sina*cosd;
		double daz =                  + dayp*cosd      + dazp*sind;

		// Accumulate acceleration terms
		particles_var1[0].ax += dax; 
		particles_var1[0].ay += day; 
		particles_var1[0].az += daz; 

	    }

        } 
       
    }

}

static void assist_additional_force_non_gravitational(struct reb_simulation* sim,
	       double xo, double yo, double zo,
	       double vxo, double vyo, double vzo,	       
	       FILE *outfile){

    const unsigned int N = sim->N;  // N includes real+variational particles
    const unsigned int N_real = N - sim->N_var;

    const double t = sim->t;    

    struct reb_particle* const particles = sim->particles;

    struct assist_extras* assist = (struct assist_extras*) sim->extras;
    struct assist_ephem* ephem = assist->ephem;
    const double jd_ref = ephem->jd_ref;

    //particle_params* part_params = NULL;
    double* part_params = NULL;        
    
    if(assist->particle_params == NULL)
	return;
    
    part_params = assist->particle_params;

    double GMsun;
    //double x, y, z, vx, vy, vz, ax, ay, az;

    // Here is the treatment of non-gravitational forces.

    double xr, yr, zr, vxr, vyr, vzr; //, axr, ayr, azr;
    
    // The Sun center is reference for these calculations.
    double xs, ys, zs, vxs, vys, vzs, axs, ays, azs;
    assist_all_ephem(ephem, assist->ephem_cache, ASSIST_BODY_SUN, t, &GMsun, &xs, &ys, &zs, &vxs, &vys, &vzs, &axs, &ays, &azs);

    xr = xs;  yr = ys;  zr = zs;
    vxr = vxs; vyr = vys; vzr = vzs;    

    // The non-grav parameters are specific to each object being
    // integrated.

    // 2020 CD3
    //double A1= 1.903810165823E-10;
    //double A2 = 0.0;
    //double A3 = 0.0;

    // Apophis
    //double A1 = 0.0;
    //double A2 = -5.592839897872E-14;
    //double A3 = 0.0;

    // 2020 SO
    //double A1 = 2.840852439404E-9;
    //double A2 = -2.521527931094E-10;
    //double A3= 2.317289821804E-10;

    // if no particles have non-zero non-grav
    // constants, skip the whole thing.
    
    // Loop over test particles
    for (int j=0; j<N_real; j++){

	//double A1 = part_params[j].A1;
	//double A2 = part_params[j].A2;
	//double A3 = part_params[j].A3;
	double A1 = part_params[3*j+0];
	double A2 = part_params[3*j+1];
	double A3 = part_params[3*j+2];

	//printf(" A123: %le %le %le\n", A1, A2, A3);
	//fflush(stdout);
	

	// If A1, A2, and A3 are zero, skip.
	if(A1==0. && A2==0. && A3==0.)
	    continue;
	
        const struct reb_particle p = particles[j];
        double dx = p.x + (xo - xr);
        double dy = p.y + (yo - yr);
        double dz = p.z + (zo - zr);

        const double r2 = dx*dx + dy*dy + dz*dz;
        const double r = sqrt(r2);

	// We may need to make this more general.
	//const double g = 1.0/r2;

	// Need to make these more flexible

	// 'Oumuamua
	//double alpha = 0.04083733261;
	//double nk = 2.6;
	//double nm = 2.0;
	//double nn = 3.0;
	//double r0 = 5.0;

	// Standard --> r^-2
	double alpha = 1.0;
	double nk = 0.0;
	double nm = 2.0;
	double nn = 5.093;
	double r0 = 1.0;
	
	const double g = alpha*pow(r/r0, -nm)*pow(1.0+pow(r/r0, nn), -nk);
	//const double g = 1.0/r2;	

	double dvx = p.vx + (vxo - vxr);
	double dvy = p.vy + (vyo - vyr);
	double dvz = p.vz + (vzo - vzr);

	double hx = dy*dvz - dz*dvy;
	double hy = dz*dvx - dx*dvz;
	double hz = dx*dvy - dy*dvx;	

	double h2 = hx*hx + hy*hy + hz*hz;
	double h = sqrt(h2);

	double tx = hy*dz - hz*dy;
	double ty = hz*dx - hx*dz;
	double tz = hx*dy - hy*dx;

        const double t2 = tx*tx + ty*ty + tz*tz;
        const double _t = sqrt(t2);

	if(outfile){	
	    fprintf(outfile, "%3s %25.16le %25.16le %25.16le %25.16le\n", "NG", jd_ref+t,
		    A1*g*dx/r + A2*g*tx/_t + A3*g*hx/h,
		    A1*g*dy/r + A2*g*ty/_t + A3*g*hy/h,
		    A1*g*dz/r + A2*g*tz/_t + A3*g*hz/h);
	    fflush(outfile);
	}

	particles[j].ax += A1*g*dx/r + A2*g*tx/_t + A3*g*hx/h;
        particles[j].ay += A1*g*dy/r + A2*g*ty/_t + A3*g*hy/h;
        particles[j].az += A1*g*dz/r + A2*g*tz/_t + A3*g*hz/h;

	// variational matrix elements
	// Only evaluate the constants if there are variational particles

        const double r3    = r*r*r;
        const double v2    = dvx*dvx + dvy*dvy + dvz*dvz;
        const double rdotv = dx*dvx  + dy*dvy  + dz*dvz;
        const double vdott = dvx*tx  + dvy*ty  + dvz*tz;

	// Need to update this for the new g(r) function. Done
	//const double dgdr = -2.*g/r;
	const double dgdr = (alpha/r0)*(-nm*pow(r/r0, -nm-1)*pow(1.0+pow(r/r0, nn), -nk)
                                      +pow(r/r0, -nm)*(-nk*nn)*pow(r/r0, nn-1)*pow(1.0+pow(r/r0, nn), -nk-1));
        const double dgx  = dgdr*dx/r;
        const double dgy  = dgdr*dy/r;
        const double dgz  = dgdr*dz/r;

        const double hxh3 = hx/(h*h*h);
        const double hyh3 = hy/(h*h*h);
        const double hzh3 = hz/(h*h*h);

        const double txt3 = tx/(_t*_t*_t);
        const double tyt3 = ty/(_t*_t*_t);
        const double tzt3 = tz/(_t*_t*_t);

	const double dxdA1 = g*dx/r;
	const double dydA1 = g*dy/r;
	const double dzdA1 = g*dz/r;

	const double dxdA2 = g*tx/_t;
	const double dydA2 = g*ty/_t;
	const double dzdA2 = g*tz/_t;

	const double dxdA3 = g*hx/h;
	const double dydA3 = g*hy/h;
	const double dzdA3 = g*hz/h;

	const double dxdx = A1*(dgx*dx/r + g*(1./r - dx*dx/r3)) 
	    + A2*(dgx*tx/_t + g*((dx*dvx - rdotv)/_t - txt3*(2.*dx*vdott - rdotv*tx)))
	    + A3*(dgx*hx/h + g*(-hxh3)*(v2*dx - rdotv*dvx));

	const double dydy = A1*(dgy*dy/r + g*(1./r - dy*dy/r3)) 
	    + A2*(dgy*ty/_t + g*((dy*dvy - rdotv)/_t - tyt3*(2.*dy*vdott - rdotv*ty)))
	    + A3*(dgy*hy/h + g*(-hyh3)*(v2*dy - rdotv*dvy));

	const double dzdz = A1*(dgz*dz/r + g*(1./r - dz*dz/r3)) 
	    + A2*(dgz*tz/_t + g*((dz*dvz - rdotv)/_t - tzt3*(2.*dz*vdott - rdotv*tz)))
	    + A3*(dgz*hz/h + g*(-hzh3)*(v2*dz - rdotv*dvz));

	const double dxdy = A1*(dgy*dx/r + g*(-dx*dy/r3))
	    + A2*(dgy*tx/_t + g*((2*dy*dvx - dx*dvy)/_t - txt3*(2*dy*vdott - rdotv*ty)))
	    + A3*(dgy*hx/h + g*(dvz/h -hxh3*(v2*dy - rdotv*dvy)));

	const double dydx = A1*(dgx*dy/r + g*(-dy*dx/r3))
	    + A2*(dgx*ty/_t + g*((2*dx*dvy - dy*dvx)/_t - tyt3*(2*dx*vdott - rdotv*tx)))
	    + A3*(dgx*hy/h + g*(-dvz/h -hyh3*(v2*dx - rdotv*dvx)));

	const double dxdz = A1*(dgz*dx/r + g*(-dx*dz/r3))
	    + A2*(dgz*tx/_t + g*((2*dz*dvx - dx*dvz)/_t - txt3*(2*dz*vdott - rdotv*tz)))
	    + A3*(dgz*hx/h + g*(-dvy/h -hxh3*(v2*dz - rdotv*dvz)));

	const double dzdx = A1*(dgx*dz/r + g*(-dz*dx/r3))
	    + A2*(dgx*tz/_t + g*((2*dx*dvz - dz*dvx)/_t - tzt3*(2*dx*vdott - rdotv*tx)))
	    + A3*(dgx*hz/h + g*(dvy/h -hzh3*(v2*dx - rdotv*dvx)));

	const double dydz = A1*(dgz*dy/r + g*(-dy*dz/r3))
	    + A2*(dgz*ty/_t + g*((2*dz*dvy - dy*dvz)/_t - tyt3*(2*dz*vdott - rdotv*tz)))
	    + A3*(dgz*hy/h + g*(dvx/h -hyh3*(v2*dz - rdotv*dvz)));

	const double dzdy = A1*(dgy*dz/r + g*(-dz*dy/r3))
	    + A2*(dgy*tz/_t + g*((2*dy*dvz - dz*dvy)/_t - tzt3*(2*dy*vdott - rdotv*ty)))
	    + A3*(dgy*hz/h + g*(-dvx/h -hzh3*(v2*dy - rdotv*dvy)));

	const double dxdvx = A1*(0.)
	    + A2*g*((dy*dy + dz*dz)/_t - txt3*r2*tx)
	    + A3*g*(-hxh3*(r2*dvx - dx*rdotv));

 	const double dydvy = A1*(0.)
	    + A2*g*((dx*dx + dz*dz)/_t - tyt3*r2*ty)
	    + A3*g*(-hyh3*(r2*dvy - dy*rdotv));

	const double dzdvz = A1*(0.)
	    + A2*g*((dx*dx + dy*dy)/_t - tzt3*r2*tz)
	    + A3*g*(-hzh3*(r2*dvz - dz*rdotv));

	const double dxdvy = A1*(0.)
	    + A2*g*(-dy*dx/_t - tyt3*r2*tx)
	    + A3*g*(-dz/h - hxh3*(r2*dvy - dy*rdotv));

	const double dydvx = A1*(0.)
	    + A2*g*(-dx*dy/_t - txt3*r2*ty)
	    + A3*g*(dz/h - hyh3*(r2*dvx - dx*rdotv));

	const double dxdvz = A1*(0.)
	    + A2*g*(-dz*dx/_t - tzt3*r2*tx)
	    + A3*g*(dy/h - hxh3*(r2*dvz - dz*rdotv));

	const double dzdvx = A1*(0.)
	    + A2*g*(-dx*dz/_t - txt3*r2*tz)
	    + A3*g*(-dy/h - hzh3*(r2*dvx - dx*rdotv));

	const double dydvz = A1*(0.)
	    + A2*g*(-dz*dy/_t - tzt3*r2*ty)
	    + A3*g*(-dx/h - hyh3*(r2*dvz - dz*rdotv));

	const double dzdvy = A1*(0.)
	    + A2*g*(-dy*dz/_t - tyt3*r2*tz)
	    + A3*g*(dx/h - hzh3*(r2*dvy - dy*rdotv));

	for (int v=0; v < sim->var_config_N; v++){
	    struct reb_variational_configuration const vc = sim->var_config[v];
	    int tp = vc.testparticle;
	    struct reb_particle* const particles_var1 = particles + vc.index;		
	    if(tp == j){
	
		// variational particle coords -- transformed to appropriate coord system.
		double ddx = particles_var1[0].x;
		double ddy = particles_var1[0].y;
		double ddz = particles_var1[0].z;
		double ddvx = particles_var1[0].vx;
		double ddvy = particles_var1[0].vy;
		double ddvz = particles_var1[0].vz;

		// Getting the variations in the non-grav params
		// There might be cleaner way to do this indexing.
		double dA1 = part_params[3*(N_real+v)+0];
		double dA2 = part_params[3*(N_real+v)+1];
		double dA3 = part_params[3*(N_real+v)+2];

		//printf("dA123: %lf %lf %lf\n", dA1, dA2, dA3);

		// Get the dA1, dA2, dA3 values.  These would normally be
		// 0 or 1, and they don't change with time

		// Matrix multiplication
		const double dax =   ddx  * dxdx  + ddy  * dxdy  + ddz  * dxdz
		    +   ddvx * dxdvx + ddvy * dxdvy + ddvz * dxdvz + dA1*dxdA1 + dA2*dxdA2 + dA3*dxdA3;
		const double day =   ddx  * dydx  + ddy  * dydy  + ddz  * dydz
		    +   ddvx * dydvx + ddvy * dydvy + ddvz * dydvz + dA1*dydA1 + dA2*dydA2 + dA3*dydA3;
		const double daz =   ddx  * dzdx  + ddy  * dzdy  + ddz  * dzdz
		    +   ddvx * dzdvx + ddvy * dzdvy + ddvz * dzdvz + dA1*dzdA1 + dA2*dzdA2 + dA3*dzdA3;

		// Accumulate acceleration terms
		particles_var1[0].ax += dax;
		particles_var1[0].ay += day;
		particles_var1[0].az += daz;

	    }
	}
	//  variational end
    }

}

static void assist_additional_force_potential_GR(struct reb_simulation* sim,
	       double xo, double yo, double zo,
	       double vxo, double vyo, double vzo,	       
	       FILE *outfile){
    
    struct assist_extras* assist = (struct assist_extras*) sim->extras;
    struct assist_ephem* ephem = assist->ephem;
    const double jd_ref = ephem->jd_ref;
    
    // Nobili and Roxburgh GR treatment

    const double au = JPL_EPHEM_CAU;    
    const double c = (JPL_EPHEM_CLIGHT/au)*86400;
    const double C2 = c*c;  
    
    const unsigned int N = sim->N;  // N includes real+variational particles
    const unsigned int N_real = N - sim->N_var;

    const double t = sim->t;    

    struct reb_particle* const particles = sim->particles;

    double GM;

    double xr, yr, zr;
    
    // The Sun center is reference for these calculations.
    double xs, ys, zs, vxs, vys, vzs, axs, ays, azs;

    assist_all_ephem(ephem, assist->ephem_cache, ASSIST_BODY_SUN, t, &GM, &xs, &ys, &zs, &vxs, &vys, &vzs, &axs, &ays, &azs);
    const double GMsun = GM;

    //xs = ys = zs = 0.0;

    xr  = xs;  yr  = ys;  zr = zs;

    for (int j=0; j<N_real; j++){

        struct reb_particle p = particles[j];

	p.x += (xo - xr);
	p.y += (yo - yr);
	p.z += (zo - zr);

        const double r2 = p.x*p.x + p.y*p.y + p.z*p.z;
        const double r = sqrt(r2);

	const double prefac = -6.0*GMsun*GMsun/(C2*r2*r2);

	particles[j].ax += prefac*p.x;
	particles[j].ay += prefac*p.y;
	particles[j].az += prefac*p.z;

	if(outfile){
	    fprintf(outfile, "%s %25.16le %25.16le %25.16le %25.16le\n", "potential GR", jd_ref+t,
		    prefac*p.x,
		    prefac*p.y,
		    prefac*p.z);
	    fflush(outfile);
	}

	// Constants for variational equations
	// Only evaluate if there are variational particles
	//const double dpdr = -prefac/r;

	// This section can be optimized.
	const double dxdx = prefac + -4.0*prefac*(p.x/r)*(p.x/r);
	const double dxdy =          -4.0*prefac*(p.x/r)*(p.y/r);
	const double dxdz =          -4.0*prefac*(p.x/r)*(p.z/r);

	// This section can be optimized.	
	const double dydx =          -4.0*prefac*(p.y/r)*(p.x/r);
	const double dydy = prefac + -4.0*prefac*(p.y/r)*(p.y/r);
	const double dydz =          -4.0*prefac*(p.y/r)*(p.z/r);

	// This section can be optimized.		
	const double dzdx =          -4.0*prefac*(p.z/r)*(p.x/r);
	const double dzdy =          -4.0*prefac*(p.z/r)*(p.y/r);
	const double dzdz = prefac + -4.0*prefac*(p.z/r)*(p.z/r);


	for (int v=0; v < sim->var_config_N; v++){
	    struct reb_variational_configuration const vc = sim->var_config[v];
	    int tp = vc.testparticle;
	    struct reb_particle* const particles_var1 = particles + vc.index;		
	    if(tp == j){
	    
		// variational particle coords
		const double ddx = particles_var1[0].x;
		const double ddy = particles_var1[0].y;
		const double ddz = particles_var1[0].z;

		// Matrix multiplication
		const double dax =   ddx  * dxdx  + ddy  * dxdy  + ddz  * dxdz;
		const double day =   ddx  * dydx  + ddy  * dydy  + ddz  * dydz;
		const double daz =   ddx  * dzdx  + ddy  * dzdy  + ddz  * dzdz;

		// Accumulate acceleration terms
		particles_var1[0].ax += dax;
		particles_var1[0].ay += day;
		particles_var1[0].az += daz;
		
	    }
	}
    }
}

static void assist_additional_force_simple_GR(struct reb_simulation* sim,
	       double xo, double yo, double zo,
	       double vxo, double vyo, double vzo,	       
	       FILE *outfile){
    
    struct assist_extras* assist = (struct assist_extras*) sim->extras;
    struct assist_ephem* ephem = assist->ephem;
    const double jd_ref = ephem->jd_ref;
    
    // Damour and Deruelle solar GR treatment

    const double au = JPL_EPHEM_CAU;    
    const double c = (JPL_EPHEM_CLIGHT/au)*86400;

    const double C2 = c*c; 
    
    const unsigned int N = sim->N;  // N includes real+variational particles
    const unsigned int N_real = N - sim->N_var;

    const double t = sim->t;    

    struct reb_particle* const particles = sim->particles;

    double GM;

    double xr, yr, zr, vxr, vyr, vzr;
    
    // The Sun center is reference for these calculations.
    double xs, ys, zs, vxs, vys, vzs, axs, ays, azs;

    assist_all_ephem(ephem, assist->ephem_cache, ASSIST_BODY_SUN, t, &GM, &xs, &ys, &zs, &vxs, &vys, &vzs, &axs, &ays, &azs);
    const double GMsun = GM;

    //xs = ys = zs = 0.0;

    xr  = xs;  yr  = ys;  zr = zs;
    vxr = vxs; vyr = vys; vzr = vzs;

    for (int j=0; j<N_real; j++){

        struct reb_particle p = particles[j];

	p.x += (xo - xr);
	p.y += (yo - yr);
	p.z += (zo - zr);
	p.vx += (vxo - vxr);
	p.vy += (vyo - vyr);
	p.vz += (vzo - vzr);
	
        const double v2 = p.vx*p.vx + p.vy*p.vy + p.vz*p.vz;
        const double r = sqrt(p.x*p.x + p.y*p.y + p.z*p.z);

	const double A = 4.0*GMsun/r - v2;
	const double B = 4.0*(p.x*p.vx + p.y*p.vy + p.z*p.vz);

	const double prefac = GMsun/(r*r*r*C2);

	particles[j].ax += prefac*(A*p.x + B*p.vx);
	particles[j].ay += prefac*(A*p.y + B*p.vy);
	particles[j].az += prefac*(A*p.z + B*p.vz);

	if(outfile){
	    fprintf(outfile, "%3s %25.16le %25.16le %25.16le %25.16le\n", "GR", jd_ref+t,
		    prefac*(A*p.x + B*p.vx),
		    prefac*(A*p.y + B*p.vy),
		    prefac*(A*p.z + B*p.vz));
	    fflush(outfile);
	}

	// Constants for variational equations
	// Only evaluate if there are variational particles
	const double dpdr = -3.0*prefac/r;

	// This section can be optimized.
	const double dxdx = dpdr*p.x/r * (A*p.x + B*p.vx) + prefac*(A - p.x*(p.x/r)*4.0*GMsun/(r*r) + 4.0*p.vx*p.vx);
	const double dxdy = dpdr*p.y/r * (A*p.x + B*p.vx) + prefac*(  - p.x*(p.y/r)*4.0*GMsun/(r*r) + 4.0*p.vy*p.vx);
	const double dxdz = dpdr*p.z/r * (A*p.x + B*p.vx) + prefac*(  - p.x*(p.z/r)*4.0*GMsun/(r*r) + 4.0*p.vz*p.vx);
	const double dxdvx =                                prefac*(  - 2.0*p.vx*p.x                + 4.0*p.x*p.vx + B);
	const double dxdvy =                                prefac*(  - 2.0*p.vy*p.x                + 4.0*p.y*p.vx    );
	const double dxdvz =                                prefac*(  - 2.0*p.vz*p.x                + 4.0*p.z*p.vx    );

	// This section can be optimized.	
	const double dydx = dpdr*p.x/r * (A*p.y + B*p.vy) + prefac*(  - p.y*(p.x/r)*4.0*GMsun/(r*r) + 4.0*p.vx*p.vy);
	const double dydy = dpdr*p.y/r * (A*p.y + B*p.vy) + prefac*(A - p.y*(p.y/r)*4.0*GMsun/(r*r) + 4.0*p.vy*p.vy);
	const double dydz = dpdr*p.z/r * (A*p.y + B*p.vy) + prefac*(  - p.y*(p.z/r)*4.0*GMsun/(r*r) + 4.0*p.vz*p.vy);
	const double dydvx =                                prefac*(  - 2.0*p.vx*p.y                + 4.0*p.x*p.vy    );
	const double dydvy =                                prefac*(  - 2.0*p.vy*p.y                + 4.0*p.y*p.vy + B);
	const double dydvz =                                prefac*(  - 2.0*p.vz*p.y                + 4.0*p.z*p.vy    );

	// This section can be optimized.		
	const double dzdx = dpdr*p.x/r * (A*p.z + B*p.vz) + prefac*(  - p.z*(p.x/r)*4.0*GMsun/(r*r) + 4.0*p.vx*p.vz);
	const double dzdy = dpdr*p.y/r * (A*p.z + B*p.vz) + prefac*(  - p.z*(p.y/r)*4.0*GMsun/(r*r) + 4.0*p.vy*p.vz);
	const double dzdz = dpdr*p.z/r * (A*p.z + B*p.vz) + prefac*(A - p.z*(p.z/r)*4.0*GMsun/(r*r) + 4.0*p.vz*p.vz);
	const double dzdvx =                                prefac*(  - 2.0*p.vx*p.z                + 4.0*p.x*p.vz    );
	const double dzdvy =                                prefac*(  - 2.0*p.vy*p.z                + 4.0*p.y*p.vz    );
	const double dzdvz =                                prefac*(  - 2.0*p.vz*p.z                + 4.0*p.z*p.vz + B);

	for (int v=0; v < sim->var_config_N; v++){
	    struct reb_variational_configuration const vc = sim->var_config[v];
	    int tp = vc.testparticle;
	    struct reb_particle* const particles_var1 = particles + vc.index;		
	    if(tp == j){
	    
		// variational particle coords
		const double ddx = particles_var1[0].x;
		const double ddy = particles_var1[0].y;
		const double ddz = particles_var1[0].z;
		const double ddvx = particles_var1[0].vx;
		const double ddvy = particles_var1[0].vy;
		const double ddvz = particles_var1[0].vz;

		// Matrix multiplication
		const double dax =   ddx  * dxdx  + ddy  * dxdy  + ddz  * dxdz
		    +   ddvx * dxdvx + ddvy * dxdvy + ddvz * dxdvz;
		const double day =   ddx  * dydx  + ddy  * dydy  + ddz  * dydz
		    +   ddvx * dydvx + ddvy * dydvy + ddvz * dydvz;
		const double daz =   ddx  * dzdx  + ddy  * dzdy  + ddz  * dzdz
		    +   ddvx * dzdvx + ddvy * dzdvy + ddvz * dzdvz;

		// Accumulate acceleration terms
		particles_var1[0].ax += dax;
		particles_var1[0].ay += day;
		particles_var1[0].az += daz;
		
	    }
	}
    }
}

static void assist_additional_force_eih_GR(struct reb_simulation* sim,
	    int eih_loop_limit,
	    double xo, double yo, double zo,
	    double vxo, double vyo, double vzo,
	    double axo, double ayo, double azo,	       	    
	    FILE *outfile,
	    FILE *eih_file){

    struct assist_extras* assist = (struct assist_extras*) sim->extras;
    struct assist_ephem* ephem = assist->ephem;
    const double jd_ref = ephem->jd_ref;

    // Einstein-Infeld-Hoffman PPN GR treatment
    // This is one of three options for GR.
    // This version is only rarely needed.

    const double au = JPL_EPHEM_CAU;    
    const double c = (JPL_EPHEM_CLIGHT/au)*86400;
    const double C2 = c*c;
    const double over_C2 = 1./(c*c);

    const unsigned int N = sim->N;  // N includes real+variational particles
    const unsigned int N_real = N - sim->N_var;

    const double t = sim->t;    

    struct reb_particle* const particles = sim->particles;

    const double beta = 1.0;
    const double gamma = 1.0;

    // First do the real particles
    // Loop over test particles        
    for (int i=0; i<N_real; i++){

	double GMj, xj, yj, zj, vxj, vyj, vzj, axj, ayj, azj;    
	double GMk, xk, yk, zk, vxk, vyk, vzk, axk, ayk, azk;

	double term7x_sum = 0.0;
	double term7y_sum = 0.0;
	double term7z_sum = 0.0;
	double term8x_sum = 0.0;
	double term8y_sum = 0.0;
	double term8z_sum = 0.0;

	double grx = 0.0;
	double gry = 0.0;
	double grz = 0.0;		

	for (int j=0; j<eih_loop_limit; j++){ // This is either 1 or ASSIST_BODY_NPLANETS

	    // Get position and mass of massive body j.
	    assist_all_ephem(ephem, assist->ephem_cache, j, t, &GMj,
		      &xj, &yj, &zj,
		      &vxj, &vyj, &vzj,
		      &axj, &ayj, &azj);

	    // Compute position vector of test particle i relative to massive body j.
	    const double dxij = particles[i].x + (xo - xj); 
	    const double dyij = particles[i].y + (yo - yj);
	    const double dzij = particles[i].z + (zo - zj);
	    const double rij2 = dxij*dxij + dyij*dyij + dzij*dzij;
	    const double _rij  = sqrt(rij2);
	    const double prefacij = GMj/(rij2*_rij);

	    // This is the place to do all the various i-j dot products
	    
	    const double vi2 = particles[i].vx*particles[i].vx +
		particles[i].vy*particles[i].vy +
		particles[i].vz*particles[i].vz;

	    const double term2 = gamma*over_C2*vi2;

	    const double vj2 = (vxj-vxo)*(vxj-vxo) + (vyj-vyo)*(vyj-vyo) + (vzj-vzo)*(vzj-vzo);

	    const double term3 = (1+gamma)*over_C2*vj2;
	    // Variational equations do not depend on term3

	    const double vidotvj = particles[i].vx*(vxj-vxo) +
		particles[i].vy*(vyj-vyo) +
		particles[i].vz*(vzj-vzo);

	    const double term4 = -2*(1+gamma)*over_C2*vidotvj;

	    const double rijdotvj = dxij*(vxj-vxo) + dyij*(vyj-vyo) + dzij*(vzj-vzo);

	    if(eih_file){
		fprintf(eih_file, " EIH_J%12d\n", j);	    
		fprintf(eih_file, "%25.16lE ", rijdotvj/_rij);
	    }

	    const double term5 = -1.5*over_C2*(rijdotvj*rijdotvj)/(_rij*_rij);

	    const double fx = (2+2*gamma)*particles[i].vx - (1+2*gamma)*(vxj-vxo);
	    const double fy = (2+2*gamma)*particles[i].vy - (1+2*gamma)*(vyj-vyo);
	    const double fz = (2+2*gamma)*particles[i].vz - (1+2*gamma)*(vzj-vzo);
	    const double f = dxij*fx + dyij*fy + dzij*fz;

	    const double prefacij_f = prefacij*f;
	    const double term7x = prefacij_f*(particles[i].vx-(vxj-vxo));
	    const double term7y = prefacij_f*(particles[i].vy-(vyj-vyo));
	    const double term7z = prefacij_f*(particles[i].vz-(vzj-vzo));
	    
	    term7x_sum += term7x;
	    term7y_sum += term7y;
	    term7z_sum += term7z;

	    double term0 = 0.0;
	    double term1 = 0.0;

	    axj = 0.0;
	    ayj = 0.0;
	    azj = 0.0;	    
	    
	    for (int k=0; k<ASSIST_BODY_NPLANETS; k++){

		// Get position and mass of massive body k.
		assist_all_ephem(ephem, assist->ephem_cache, k, t, &GMk,
			  &xk, &yk, &zk,
			  &vxk, &vyk, &vzk,
			  &axk, &ayk, &azk);

		// Compute position vector of test particle i relative to massive body k.
		const double dxik = particles[i].x + (xo - xk); 
		const double dyik = particles[i].y + (yo - yk);
		const double dzik = particles[i].z + (zo - zk);
		const double rik2 = dxik*dxik + dyik*dyik + dzik*dzik;
		const double _rik  = sqrt(rik2);

		// keep track of GM/rik sum
		term0 += GMk/_rik;

		if(k != j){
		    // Compute position vector of massive body j relative to massive body k.
		    const double dxjk = xj - xk;
		    const double dyjk = yj - yk;
		    const double dzjk = zj - zk;
		    const double rjk2 = dxjk*dxjk + dyjk*dyjk + dzjk*dzjk;
		    const double _rjk  = sqrt(rjk2);

		    // keep track of GM/rjk sum
		    term1 += GMk/_rjk;

		    const double fac = GMk/(rjk2*_rjk);
		    axj -= fac*dxjk;
		    ayj -= fac*dyjk;
		    azj -= fac*dzjk;

		}

	    }

	    term0 *= -2*(beta+gamma)*over_C2;
	    
	    term1 *= -(2*beta-1)*over_C2;

	    const double rijdotaj = dxij*(axj-axo) + dyij*(ayj-ayo) + dzij*(azj-azo);
	    const double term6 = -0.5*over_C2*rijdotaj;

	    const double term8_fac = GMj/_rij*(3+4*gamma)/2;
	    const double term8x = term8_fac*axj;
	    const double term8y = term8_fac*ayj;
	    const double term8z = term8_fac*azj;

	    term8x_sum += term8x;
	    term8y_sum += term8y;
	    term8z_sum += term8z;

	    const double factor = term0 + term1 + term2 + term3 + term4 + term5 + term6;

	    if(eih_file){
		fprintf(eih_file, "%24.16lE ", -factor*C2);
		fprintf(eih_file, "%24.16lE %24.16lE %24.16lE %24.16lE ",
			-factor*C2*prefacij*dxij,
			-factor*C2*prefacij*dyij,
			-factor*C2*prefacij*dzij,
			f);	    
		fprintf(eih_file, "%24.16lE %24.16lE %24.16lE ",
			prefacij*f*(particles[i].vx-(vxj-vxo)),
			prefacij*f*(particles[i].vy-(vyj-vyo)),
			prefacij*f*(particles[i].vz-(vzj-vzo)));	    
		fprintf(eih_file, "%24.16lE %24.16lE %24.16lE ",
			term8x,
			term8y,
			term8z);
		fprintf(eih_file, "%24.16lE %24.16lE %24.16lE\n",
			axj, ayj, azj);

		fflush(eih_file);
	    }

	    grx += -prefacij*dxij*factor;
	    gry += -prefacij*dyij*factor;
	    grz += -prefacij*dzij*factor;
	    
	    particles[i].ax += -prefacij*dxij*factor;
	    particles[i].ay += -prefacij*dyij*factor;
	    particles[i].az += -prefacij*dzij*factor;

        }

	grx += term7x_sum*over_C2 + term8x_sum*over_C2;
	gry += term7y_sum*over_C2 + term8y_sum*over_C2;
	grz += term7z_sum*over_C2 + term8z_sum*over_C2;

	if(outfile){
	    fprintf(outfile, "%3s %25.16le %25.16le %25.16le %25.16le\n", "GR", jd_ref+t,
		    grx, gry, grz);
	    fflush(outfile);
	}

	particles[i].ax += term7x_sum*over_C2 + term8x_sum*over_C2;
	particles[i].ay += term7y_sum*over_C2 + term8y_sum*over_C2;
	particles[i].az += term7z_sum*over_C2 + term8z_sum*over_C2;

    }

    if(sim->var_config_N==0)
	return;
    
    // Now do the variational particles
    // Loop over test particles        
    for (int i=0; i<N_real; i++){

	double GMj, xj, yj, zj, vxj, vyj, vzj, axj, ayj, azj;    
	double GMk, xk, yk, zk, vxk, vyk, vzk, axk, ayk, azk;

	// Declare and initialize variational terms
	// Only do this if the variational terms are needed.
	double dxdx = 0.0;
	double dxdy = 0.0;
	double dxdz = 0.0;    
	double dxdvx = 0.0;
	double dxdvy = 0.0;
	double dxdvz = 0.0;    
	double dydx = 0.0;
	double dydy = 0.0;
	double dydz = 0.0;    
	double dydvx = 0.0;
	double dydvy = 0.0;
	double dydvz = 0.0;    
	double dzdx = 0.0;
	double dzdy = 0.0;
	double dzdz = 0.0;    
	double dzdvx = 0.0;
	double dzdvy = 0.0;
	double dzdvz = 0.0;    

	double term7x_sum = 0.0;
	double dterm7x_sumdx = 0.0;
	double dterm7x_sumdy = 0.0;
	double dterm7x_sumdz = 0.0;		
	double dterm7x_sumdvx = 0.0;
	double dterm7x_sumdvy = 0.0;
	double dterm7x_sumdvz = 0.0;		

	double term7y_sum = 0.0;
	double dterm7y_sumdx = 0.0;
	double dterm7y_sumdy = 0.0;
	double dterm7y_sumdz = 0.0;		
	double dterm7y_sumdvx = 0.0;
	double dterm7y_sumdvy = 0.0;
	double dterm7y_sumdvz = 0.0;		
	
	double term7z_sum = 0.0;
	double dterm7z_sumdx = 0.0;
	double dterm7z_sumdy = 0.0;
	double dterm7z_sumdz = 0.0;		
	double dterm7z_sumdvx = 0.0;
	double dterm7z_sumdvy = 0.0;
	double dterm7z_sumdvz = 0.0;		
	
	double term8x_sum = 0.0;
	double dterm8x_sumdx = 0.0;
	double dterm8x_sumdy = 0.0;
	double dterm8x_sumdz = 0.0;		

	double term8y_sum = 0.0;
	double dterm8y_sumdx = 0.0;
	double dterm8y_sumdy = 0.0;
	double dterm8y_sumdz = 0.0;		
	
	double term8z_sum = 0.0;
	double dterm8z_sumdx = 0.0;
	double dterm8z_sumdy = 0.0;
	double dterm8z_sumdz = 0.0;		

	double grx = 0.0;
	double gry = 0.0;
	double grz = 0.0;		

	for (int j=0; j<eih_loop_limit; j++){ // This is either 1 or ASSIST_BODY_NPLANETS

	    // Get position and mass of massive body j.
	    assist_all_ephem(ephem, assist->ephem_cache, j, t, &GMj,
		      &xj, &yj, &zj,
		      &vxj, &vyj, &vzj,
		      &axj, &ayj, &azj);

	    // Compute position vector of test particle i relative to massive body j.
	    const double dxij = particles[i].x + (xo - xj); 
	    const double dyij = particles[i].y + (yo - yj);
	    const double dzij = particles[i].z + (zo - zj);
	    const double rij2 = dxij*dxij + dyij*dyij + dzij*dzij;
	    const double _rij  = sqrt(rij2);
	    const double prefacij = GMj/(_rij*_rij*_rij);

	    const double dprefacijdx = -3.0*GMj/(_rij*_rij*_rij*_rij*_rij)*dxij;
	    const double dprefacijdy = -3.0*GMj/(_rij*_rij*_rij*_rij*_rij)*dyij;
	    const double dprefacijdz = -3.0*GMj/(_rij*_rij*_rij*_rij*_rij)*dzij;

	    // This is the place to do all the various i-j dot products
	    
	    const double vi2 = particles[i].vx*particles[i].vx +
		particles[i].vy*particles[i].vy +
		particles[i].vz*particles[i].vz;

	    const double term2 = gamma*over_C2*vi2;
	    const double dterm2dvx = 2.0*gamma*over_C2*particles[i].vx;
	    const double dterm2dvy = 2.0*gamma*over_C2*particles[i].vy;
	    const double dterm2dvz = 2.0*gamma*over_C2*particles[i].vz;	    

	    const double vj2 = (vxj-vxo)*(vxj-vxo) + (vyj-vyo)*(vyj-vyo) + (vzj-vzo)*(vzj-vzo);

	    const double term3 = (1+gamma)*over_C2*vj2;
	    // Variational equations do not depend on term3

	    const double vidotvj = particles[i].vx*(vxj-vxo) +
		particles[i].vy*(vyj-vyo) +
		particles[i].vz*(vzj-vzo);

	    const double term4 = -2*(1+gamma)*over_C2*vidotvj;
	    const double dterm4dvx = -2*(1+gamma)*over_C2*(vxj-vxo);
	    const double dterm4dvy = -2*(1+gamma)*over_C2*(vyj-vyo);
	    const double dterm4dvz = -2*(1+gamma)*over_C2*(vzj-vzo);	    	    
	    

	    const double rijdotvj = dxij*(vxj-vxo) + dyij*(vyj-vyo) + dzij*(vzj-vzo);

	    if(eih_file){
		fprintf(eih_file, " EIH_J%12d\n", j);	    
		fprintf(eih_file, "%25.16lE ", rijdotvj/_rij);
	    }

	    const double term5 = -1.5*over_C2*(rijdotvj*rijdotvj)/(_rij*_rij);
	    const double term5_fac = 3.0*over_C2*rijdotvj/_rij;
	    const double dterm5dx = -term5_fac*((vxj-vxo)/_rij - rijdotvj*dxij/(_rij*_rij*_rij));
	    const double dterm5dy = -term5_fac*((vyj-vyo)/_rij - rijdotvj*dyij/(_rij*_rij*_rij));
	    const double dterm5dz = -term5_fac*((vzj-vzo)/_rij - rijdotvj*dzij/(_rij*_rij*_rij));	    	    

	    double fx = (2+2*gamma)*particles[i].vx - (1+2*gamma)*(vxj-vxo);
	    double fy = (2+2*gamma)*particles[i].vy - (1+2*gamma)*(vyj-vyo);
	    double fz = (2+2*gamma)*particles[i].vz - (1+2*gamma)*(vzj-vzo);
	    double f = dxij*fx + dyij*fy + dzij*fz;

	    double dfdx = fx;
	    double dfdy = fy;
	    double dfdz = fz;	    	    
	    double dfdvx = dxij*(2+2*gamma);
	    double dfdvy = dyij*(2+2*gamma);
	    double dfdvz = dzij*(2+2*gamma);

	    double term7x = prefacij*f*(particles[i].vx-(vxj-vxo));
	    double term7y = prefacij*f*(particles[i].vy-(vyj-vyo));
	    double term7z = prefacij*f*(particles[i].vz-(vzj-vzo));

	    double dterm7xdx = dprefacijdx * f * (particles[i].vx-(vxj-vxo))
		+ prefacij * dfdx * (particles[i].vx-(vxj-vxo));
	    double dterm7xdy = dprefacijdy * f * (particles[i].vx-(vxj-vxo))
		+ prefacij * dfdy * (particles[i].vx-(vxj-vxo));
	    double dterm7xdz = dprefacijdz * f * (particles[i].vx-(vxj-vxo))
		+ prefacij * dfdz * (particles[i].vx-(vxj-vxo));
	    double dterm7xdvx = prefacij * dfdvx * (particles[i].vx-(vxj-vxo))
		+ prefacij * f;
	    double dterm7xdvy = prefacij * dfdvy * (particles[i].vx-(vxj-vxo));

	    double dterm7xdvz = prefacij * dfdvz * (particles[i].vx-(vxj-vxo));	    

	    double dterm7ydx = dprefacijdx * f * (particles[i].vy-(vyj-vyo))
		+ prefacij * dfdx * (particles[i].vy-(vyj-vyo));
	    double dterm7ydy = dprefacijdy * f * (particles[i].vy-(vyj-vyo))
		+ prefacij * dfdy * (particles[i].vy-(vyj-vyo));
	    double dterm7ydz = dprefacijdz * f * (particles[i].vy-(vyj-vyo))
		+ prefacij * dfdz * (particles[i].vy-(vyj-vyo));
	    double dterm7ydvx = prefacij * dfdvx * (particles[i].vy-(vyj-vyo));

	    double dterm7ydvy = prefacij * dfdvy * (particles[i].vy-(vyj-vyo))
		+ prefacij * f;		
	    double dterm7ydvz = prefacij * dfdvz * (particles[i].vy-(vyj-vyo));	    

	    double dterm7zdx = dprefacijdx * f * (particles[i].vz-(vzj-vzo))
		+ prefacij * dfdx * (particles[i].vz-(vzj-vzo));
	    double dterm7zdy = dprefacijdy * f * (particles[i].vz-(vzj-vzo))
		+ prefacij * dfdy * (particles[i].vz-(vzj-vzo));
	    double dterm7zdz = dprefacijdz * f * (particles[i].vz-(vzj-vzo))
		+ prefacij * dfdz * (particles[i].vz-(vzj-vzo));

	    double dterm7zdvx = prefacij * dfdvx * (particles[i].vz-(vzj-vzo));

	    double dterm7zdvy = prefacij * dfdvy * (particles[i].vz-(vzj-vzo));

	    double dterm7zdvz = prefacij * dfdvz * (particles[i].vz-(vzj-vzo))
		+ prefacij * f;
	    
	    term7x_sum += term7x;
	    term7y_sum += term7y;
	    term7z_sum += term7z;

	    dterm7x_sumdx += dterm7xdx;
	    dterm7x_sumdy += dterm7xdy;
	    dterm7x_sumdz += dterm7xdz;	    
	    dterm7x_sumdvx += dterm7xdvx;
	    dterm7x_sumdvy += dterm7xdvy;
	    dterm7x_sumdvz += dterm7xdvz;	    

	    dterm7y_sumdx += dterm7ydx;
	    dterm7y_sumdy += dterm7ydy;
	    dterm7y_sumdz += dterm7ydz;
	    dterm7y_sumdvx += dterm7ydvx;
	    dterm7y_sumdvy += dterm7ydvy;
	    dterm7y_sumdvz += dterm7ydvz;

	    dterm7z_sumdx += dterm7zdx;
	    dterm7z_sumdy += dterm7zdy;
	    dterm7z_sumdz += dterm7zdz;	    
	    dterm7z_sumdvx += dterm7zdvx;
	    dterm7z_sumdvy += dterm7zdvy;
	    dterm7z_sumdvz += dterm7zdvz;	    
	    

	    double term0 = 0.0;
	    double dterm0dx = 0.0;
	    double dterm0dy = 0.0;
	    double dterm0dz = 0.0;	    

	    double term1 = 0.0;
	    double dterm1dx = 0.0;
	    double dterm1dy = 0.0;
	    double dterm1dz = 0.0;	    
	    double dterm1dvx = 0.0;
	    double dterm1dvy = 0.0;
	    double dterm1dvz = 0.0;	    

	    axj = 0.0;
	    ayj = 0.0;
	    azj = 0.0;	    
	    
	    for (int k=0; k<ASSIST_BODY_NPLANETS; k++){

		// Get position and mass of massive body k.
		assist_all_ephem(ephem, assist->ephem_cache, k, t, &GMk,
			  &xk, &yk, &zk,
			  &vxk, &vyk, &vzk,
			  &axk, &ayk, &azk);

		// Compute position vector of test particle i relative to massive body k.
		const double dxik = particles[i].x + (xo - xk); 
		const double dyik = particles[i].y + (yo - yk);
		const double dzik = particles[i].z + (zo - zk);
		const double rik2 = dxik*dxik + dyik*dyik + dzik*dzik;
		const double _rik  = sqrt(rik2);

		// keep track of GM/rik sum
		term0 += GMk/_rik;

		dterm0dx -= GMk/(_rik*_rik*_rik) * dxik;
		dterm0dy -= GMk/(_rik*_rik*_rik) * dyik;
		dterm0dz -= GMk/(_rik*_rik*_rik) * dzik;				

		if(k != j){
		    // Compute position vector of massive body j relative to massive body k.
		    const double dxjk = xj - xk;
		    const double dyjk = yj - yk;
		    const double dzjk = zj - zk;
		    const double rjk2 = dxjk*dxjk + dyjk*dyjk + dzjk*dzjk;
		    const double _rjk  = sqrt(rjk2);

		    // keep track of GM/rjk sum
		    term1 += GMk/_rjk;

		    axj -= GMk*dxjk/(_rjk*_rjk*_rjk);
		    ayj -= GMk*dyjk/(_rjk*_rjk*_rjk);
		    azj -= GMk*dzjk/(_rjk*_rjk*_rjk);		    		    

		}

	    }

	    term0 *= -2*(beta+gamma)*over_C2;
	    dterm0dx *= -2*(beta+gamma)*over_C2;
	    dterm0dy *= -2*(beta+gamma)*over_C2;
	    dterm0dz *= -2*(beta+gamma)*over_C2;	    	    
	    
	    term1 *= -(2*beta-1)*over_C2;

	    const double rijdotaj = dxij*(axj-axo) + dyij*(ayj-ayo) + dzij*(azj-azo);
	    const double term6 = -0.5*over_C2*rijdotaj;
	    const double dterm6dx = -0.5*over_C2*(axj-axo);
	    const double dterm6dy = -0.5*over_C2*(ayj-ayo);	    
	    const double dterm6dz = -0.5*over_C2*(azj-azo);
	    
	    double term8x = GMj*axj/_rij*(3+4*gamma)/2;
	    double dterm8xdx = -GMj*axj/(_rij*_rij*_rij)*dxij*(3+4*gamma)/2;
	    double dterm8xdy = -GMj*axj/(_rij*_rij*_rij)*dyij*(3+4*gamma)/2;
	    double dterm8xdz = -GMj*axj/(_rij*_rij*_rij)*dzij*(3+4*gamma)/2;	    	    

	    double term8y = GMj*ayj/_rij*(3+4*gamma)/2;
	    double dterm8ydx = -GMj*ayj/(_rij*_rij*_rij)*dxij*(3+4*gamma)/2;
	    double dterm8ydy = -GMj*ayj/(_rij*_rij*_rij)*dyij*(3+4*gamma)/2;
	    double dterm8ydz = -GMj*ayj/(_rij*_rij*_rij)*dzij*(3+4*gamma)/2;	    	    

	    double term8z = GMj*azj/_rij*(3+4*gamma)/2;
	    double dterm8zdx = -GMj*azj/(_rij*_rij*_rij)*dxij*(3+4*gamma)/2;
	    double dterm8zdy = -GMj*azj/(_rij*_rij*_rij)*dyij*(3+4*gamma)/2;
	    double dterm8zdz = -GMj*azj/(_rij*_rij*_rij)*dzij*(3+4*gamma)/2;	    	    

	    term8x_sum += term8x;
	    term8y_sum += term8y;
	    term8z_sum += term8z;

	    dterm8x_sumdx += dterm8xdx;
	    dterm8x_sumdy += dterm8xdy;
	    dterm8x_sumdz += dterm8xdz;	    

	    dterm8y_sumdx += dterm8ydx;
	    dterm8y_sumdy += dterm8ydy;
	    dterm8y_sumdz += dterm8ydz;

	    dterm8z_sumdx += dterm8zdx;
	    dterm8z_sumdy += dterm8zdy;
	    dterm8z_sumdz += dterm8zdz;	    
	    
	    double factor = term0 + term1 + term2 + term3 + term4 + term5 + term6;

	    double dfactordx = dterm0dx + dterm1dx + dterm5dx + dterm6dx;
	    double dfactordy = dterm0dy + dterm1dy + dterm5dy + dterm6dy;
	    double dfactordz = dterm0dz + dterm1dz + dterm5dz + dterm6dz;	    
	    double dfactordvx = dterm1dvx + dterm2dvx + dterm4dvx;
	    double dfactordvy = dterm1dvy + dterm2dvy + dterm4dvy;
	    double dfactordvz = dterm1dvz + dterm2dvz + dterm4dvz;

	    if(eih_file){
		fprintf(eih_file, "%24.16lE ", -factor*C2);
		fprintf(eih_file, "%24.16lE %24.16lE %24.16lE %24.16lE ",
			-factor*C2*prefacij*dxij,
			-factor*C2*prefacij*dyij,
			-factor*C2*prefacij*dzij,
			f);	    
		fprintf(eih_file, "%24.16lE %24.16lE %24.16lE ",
			prefacij*f*(particles[i].vx-(vxj-vxo)),
			prefacij*f*(particles[i].vy-(vyj-vyo)),
			prefacij*f*(particles[i].vz-(vzj-vzo)));	    
		fprintf(eih_file, "%24.16lE %24.16lE %24.16lE ",
			term8x,
			term8y,
			term8z);
		fprintf(eih_file, "%24.16lE %24.16lE %24.16lE\n",
			axj, ayj, azj);

		fflush(eih_file);
	    }

	    grx += -prefacij*dxij*factor;
	    gry += -prefacij*dyij*factor;
	    grz += -prefacij*dzij*factor;
	    
	    //particles[i].ax += -prefacij*dxij*factor;
	    //particles[i].ay += -prefacij*dyij*factor;
	    //particles[i].az += -prefacij*dzij*factor;

	    // Variational equation terms go here.

	    dxdx += -dprefacijdx*dxij*factor
		-prefacij*factor
		-prefacij*dxij*dfactordx;
	    
	    dxdy += -dprefacijdy*dxij*factor
		-prefacij*dxij*dfactordy;
	    
	    dxdz += -dprefacijdz*dxij*factor
		-prefacij*dxij*dfactordz;

	    dxdvx += 
		-prefacij*dxij*dfactordvx;

	    dxdvy += 
		-prefacij*dxij*dfactordvy;

	    dxdvz += 
		-prefacij*dxij*dfactordvz;

	    dydx += -dprefacijdx*dyij*factor
		-prefacij*dyij*dfactordx;
	    
	    dydy += -dprefacijdy*dyij*factor
		-prefacij*factor
		-prefacij*dyij*dfactordy;
	    
	    dydz += -dprefacijdz*dyij*factor
		-prefacij*dyij*dfactordz;

	    dydvx += 
		-prefacij*dyij*dfactordvx;

	    dydvy += 
		-prefacij*dyij*dfactordvy;

	    dydvz += 
		-prefacij*dyij*dfactordvz;
	    
	    dzdx += -dprefacijdx*dzij*factor
		-prefacij*dzij*dfactordx;

	    dzdy += -dprefacijdy*dzij*factor
		-prefacij*dzij*dfactordy;
	    
	    dzdz += -dprefacijdz*dzij*factor
		-prefacij*factor
		-prefacij*dzij*dfactordz;

	    dzdvx += 
		-prefacij*dzij*dfactordvx;

	    dzdvy += 
		-prefacij*dzij*dfactordvy;

	    dzdvz += 
		-prefacij*dzij*dfactordvz;

        }

	grx += term7x_sum*over_C2 + term8x_sum*over_C2;
	gry += term7y_sum*over_C2 + term8y_sum*over_C2;
	grz += term7z_sum*over_C2 + term8z_sum*over_C2;

	if(outfile){
	    fprintf(outfile, "%3s %25.16le %25.16le %25.16le %25.16le\n", "GR", jd_ref+t,
		    grx, gry, grz);
	    fflush(outfile);
	}

	dxdx += dterm7x_sumdx*over_C2 + dterm8x_sumdx*over_C2;
	dxdy += dterm7x_sumdy*over_C2 + dterm8x_sumdy*over_C2;
	dxdz += dterm7x_sumdz*over_C2 + dterm8x_sumdz*over_C2;	
	dxdvx += dterm7x_sumdvx*over_C2;
	dxdvy += dterm7x_sumdvy*over_C2;
	dxdvz += dterm7x_sumdvz*over_C2;

	dydx += dterm7y_sumdx*over_C2 + dterm8y_sumdx*over_C2;
	dydy += dterm7y_sumdy*over_C2 + dterm8y_sumdy*over_C2;
	dydz += dterm7y_sumdz*over_C2 + dterm8y_sumdz*over_C2;	
	dydvx += dterm7y_sumdvx*over_C2;
	dydvy += dterm7y_sumdvy*over_C2;
	dydvz += dterm7y_sumdvz*over_C2;

	dzdx += dterm7z_sumdx*over_C2 + dterm8z_sumdx*over_C2;
	dzdy += dterm7z_sumdy*over_C2 + dterm8z_sumdy*over_C2;
	dzdz += dterm7z_sumdz*over_C2 + dterm8z_sumdz*over_C2;	
	dzdvx += dterm7z_sumdvx*over_C2;
	dzdvy += dterm7z_sumdvy*over_C2;
	dzdvz += dterm7z_sumdvz*over_C2;
	
	//particles[i].ax += term7x_sum/C2 + term8x_sum/C2;
	//particles[i].ay += term7y_sum/C2 + term8y_sum/C2;
	//particles[i].az += term7z_sum/C2 + term8z_sum/C2;

	// Variational equation terms go here.
	for (int v=0; v < sim->var_config_N; v++){
	    struct reb_variational_configuration const vc = sim->var_config[v];
	    int tp = vc.testparticle;
	    struct reb_particle* const particles_var1 = particles + vc.index;		
	    if(tp == i){
	    
		// variational particle coords
		const double ddx = particles_var1[0].x;
		const double ddy = particles_var1[0].y;
		const double ddz = particles_var1[0].z;
		const double ddvx = particles_var1[0].vx;
		const double ddvy = particles_var1[0].vy;
		const double ddvz = particles_var1[0].vz;

		// Matrix multiplication
		const double dax =   ddx  * dxdx  + ddy  * dxdy  + ddz  * dxdz
		    +   ddvx * dxdvx + ddvy * dxdvy + ddvz * dxdvz;
		const double day =   ddx  * dydx  + ddy  * dydy  + ddz  * dydz
		    +   ddvx * dydvx + ddvy * dydvy + ddvz * dydvz;
		const double daz =   ddx  * dzdx  + ddy  * dzdy  + ddz  * dzdz
		    +   ddvx * dzdvx + ddvy * dzdvy + ddvz * dzdvz;

		// Accumulate acceleration terms
		particles_var1[0].ax += dax;
		particles_var1[0].ay += day;
		particles_var1[0].az += daz;
		
	    }
	}
    }
}

static void assist_additional_force_eih_GR_orig(struct reb_simulation* sim,
	    int eih_loop_limit,
	    double xo, double yo, double zo,
	    double vxo, double vyo, double vzo,
	    double axo, double ayo, double azo,	       	    
	    FILE *outfile,
	    FILE *eih_file){

    struct assist_extras* assist = (struct assist_extras*) sim->extras;
    struct assist_ephem* ephem = assist->ephem;
    const double jd_ref = ephem->jd_ref;

    // Einstein-Infeld-Hoffman PPN GR treatment
    // This is one of three options for GR.
    // This version is only rarely needed.

    const double au = JPL_EPHEM_CAU;    
    const double c = (JPL_EPHEM_CLIGHT/au)*86400;
    const double C2 = c*c;  // This could be stored as C2.
    const double over_C2 = 1./(c*c);    

    const unsigned int N = sim->N;  // N includes real+variational particles
    const unsigned int N_real = N - sim->N_var;

    const double t = sim->t;    

    struct reb_particle* const particles = sim->particles;

    double GM;
    
    // The Sun center is reference for these calculations.
    double xs, ys, zs, vxs, vys, vzs, axs, ays, azs;
    assist_all_ephem(ephem, assist->ephem_cache, ASSIST_BODY_SUN, t, &GM, &xs, &ys, &zs, &vxs, &vys, &vzs, &axs, &ays, &azs);

    double beta = 1.0;
    double gamma = 1.0;

    // Loop over test particles        
    for (int i=0; i<N_real; i++){

	double GMj, xj, yj, zj, vxj, vyj, vzj, axj, ayj, azj;    
	double GMk, xk, yk, zk, vxk, vyk, vzk, axk, ayk, azk;

	// Declare and initialize variational terms
	// Only do this if the variational terms are needed.
	double dxdx = 0.0;
	double dxdy = 0.0;
	double dxdz = 0.0;    
	double dxdvx = 0.0;
	double dxdvy = 0.0;
	double dxdvz = 0.0;    
	double dydx = 0.0;
	double dydy = 0.0;
	double dydz = 0.0;    
	double dydvx = 0.0;
	double dydvy = 0.0;
	double dydvz = 0.0;    
	double dzdx = 0.0;
	double dzdy = 0.0;
	double dzdz = 0.0;    
	double dzdvx = 0.0;
	double dzdvy = 0.0;
	double dzdvz = 0.0;    

	double term7x_sum = 0.0;
	double dterm7x_sumdx = 0.0;
	double dterm7x_sumdy = 0.0;
	double dterm7x_sumdz = 0.0;		
	double dterm7x_sumdvx = 0.0;
	double dterm7x_sumdvy = 0.0;
	double dterm7x_sumdvz = 0.0;		

	double term7y_sum = 0.0;
	double dterm7y_sumdx = 0.0;
	double dterm7y_sumdy = 0.0;
	double dterm7y_sumdz = 0.0;		
	double dterm7y_sumdvx = 0.0;
	double dterm7y_sumdvy = 0.0;
	double dterm7y_sumdvz = 0.0;		
	
	double term7z_sum = 0.0;
	double dterm7z_sumdx = 0.0;
	double dterm7z_sumdy = 0.0;
	double dterm7z_sumdz = 0.0;		
	double dterm7z_sumdvx = 0.0;
	double dterm7z_sumdvy = 0.0;
	double dterm7z_sumdvz = 0.0;		
	
	double term8x_sum = 0.0;
	double dterm8x_sumdx = 0.0;
	double dterm8x_sumdy = 0.0;
	double dterm8x_sumdz = 0.0;		

	double term8y_sum = 0.0;
	double dterm8y_sumdx = 0.0;
	double dterm8y_sumdy = 0.0;
	double dterm8y_sumdz = 0.0;		
	
	double term8z_sum = 0.0;
	double dterm8z_sumdx = 0.0;
	double dterm8z_sumdy = 0.0;
	double dterm8z_sumdz = 0.0;		

	double grx = 0.0;
	double gry = 0.0;
	double grz = 0.0;		

	for (int j=0; j<eih_loop_limit; j++){ // This is either 1 or ASSIST_BODY_NPLANETS

	    // Get position and mass of massive body j.
	    assist_all_ephem(ephem, assist->ephem_cache, j, t, &GMj,
		      &xj, &yj, &zj,
		      &vxj, &vyj, &vzj,
		      &axj, &ayj, &azj);

	    // Compute position vector of test particle i relative to massive body j.
	    const double dxij = particles[i].x + (xo - xj); 
	    const double dyij = particles[i].y + (yo - yj);
	    const double dzij = particles[i].z + (zo - zj);
	    const double rij2 = dxij*dxij + dyij*dyij + dzij*dzij;
	    const double _rij  = sqrt(rij2);
	    const double prefacij = GMj/(_rij*_rij*_rij);

	    const double dprefacijdx = -3.0*GMj/(_rij*_rij*_rij*_rij*_rij)*dxij;
	    const double dprefacijdy = -3.0*GMj/(_rij*_rij*_rij*_rij*_rij)*dyij;
	    const double dprefacijdz = -3.0*GMj/(_rij*_rij*_rij*_rij*_rij)*dzij;

	    // This is the place to do all the various i-j dot products
	    
	    const double vi2 = particles[i].vx*particles[i].vx +
		particles[i].vy*particles[i].vy +
		particles[i].vz*particles[i].vz;

	    const double term2 = gamma*over_C2*vi2;
	    const double dterm2dvx = 2.0*gamma*over_C2*particles[i].vx;
	    const double dterm2dvy = 2.0*gamma*over_C2*particles[i].vy;
	    const double dterm2dvz = 2.0*gamma*over_C2*particles[i].vz;	    

	    const double vj2 = (vxj-vxo)*(vxj-vxo) + (vyj-vyo)*(vyj-vyo) + (vzj-vzo)*(vzj-vzo);

	    const double term3 = (1+gamma)*over_C2*vj2;
	    // Variational equations do not depend on term3

	    const double vidotvj = particles[i].vx*(vxj-vxo) +
		particles[i].vy*(vyj-vyo) +
		particles[i].vz*(vzj-vzo);

	    const double term4 = -2*(1+gamma)*over_C2*vidotvj;
	    const double dterm4dvx = -2*(1+gamma)*over_C2*(vxj-vxo);
	    const double dterm4dvy = -2*(1+gamma)*over_C2*(vyj-vyo);
	    const double dterm4dvz = -2*(1+gamma)*over_C2*(vzj-vzo);	    	    
	    

	    const double rijdotvj = dxij*(vxj-vxo) + dyij*(vyj-vyo) + dzij*(vzj-vzo);

	    if(eih_file){
		fprintf(eih_file, " EIH_J%12d\n", j);	    
		fprintf(eih_file, "%25.16lE ", rijdotvj/_rij);
	    }

	    const double term5 = -1.5*over_C2*(rijdotvj*rijdotvj)/(_rij*_rij);
	    const double dterm5dx = -3.0*over_C2*rijdotvj/_rij*((vxj-vxo)/_rij - rijdotvj*dxij/(_rij*_rij*_rij));
	    const double dterm5dy = -3.0*over_C2*rijdotvj/_rij*((vyj-vyo)/_rij - rijdotvj*dyij/(_rij*_rij*_rij));
	    const double dterm5dz = -3.0*over_C2*rijdotvj/_rij*((vzj-vzo)/_rij - rijdotvj*dzij/(_rij*_rij*_rij));	    	    

	    double fx = (2+2*gamma)*particles[i].vx - (1+2*gamma)*(vxj-vxo);
	    double fy = (2+2*gamma)*particles[i].vy - (1+2*gamma)*(vyj-vyo);
	    double fz = (2+2*gamma)*particles[i].vz - (1+2*gamma)*(vzj-vzo);
	    double f = dxij*fx + dyij*fy + dzij*fz;

	    double dfdx = fx;
	    double dfdy = fy;
	    double dfdz = fz;	    	    
	    double dfdvx = dxij*(2+2*gamma);
	    double dfdvy = dyij*(2+2*gamma);
	    double dfdvz = dzij*(2+2*gamma);

	    double term7x = prefacij*f*(particles[i].vx-(vxj-vxo));
	    double term7y = prefacij*f*(particles[i].vy-(vyj-vyo));
	    double term7z = prefacij*f*(particles[i].vz-(vzj-vzo));

	    double dterm7xdx = dprefacijdx * f * (particles[i].vx-(vxj-vxo))
		+ prefacij * dfdx * (particles[i].vx-(vxj-vxo));
	    double dterm7xdy = dprefacijdy * f * (particles[i].vx-(vxj-vxo))
		+ prefacij * dfdy * (particles[i].vx-(vxj-vxo));
	    double dterm7xdz = dprefacijdz * f * (particles[i].vx-(vxj-vxo))
		+ prefacij * dfdz * (particles[i].vx-(vxj-vxo));
	    double dterm7xdvx = prefacij * dfdvx * (particles[i].vx-(vxj-vxo))
		+ prefacij * f;
	    double dterm7xdvy = prefacij * dfdvy * (particles[i].vx-(vxj-vxo));

	    double dterm7xdvz = prefacij * dfdvz * (particles[i].vx-(vxj-vxo));	    

	    double dterm7ydx = dprefacijdx * f * (particles[i].vy-(vyj-vyo))
		+ prefacij * dfdx * (particles[i].vy-(vyj-vyo));
	    double dterm7ydy = dprefacijdy * f * (particles[i].vy-(vyj-vyo))
		+ prefacij * dfdy * (particles[i].vy-(vyj-vyo));
	    double dterm7ydz = dprefacijdz * f * (particles[i].vy-(vyj-vyo))
		+ prefacij * dfdz * (particles[i].vy-(vyj-vyo));
	    double dterm7ydvx = prefacij * dfdvx * (particles[i].vy-(vyj-vyo));

	    double dterm7ydvy = prefacij * dfdvy * (particles[i].vy-(vyj-vyo))
		+ prefacij * f;		
	    double dterm7ydvz = prefacij * dfdvz * (particles[i].vy-(vyj-vyo));	    

	    double dterm7zdx = dprefacijdx * f * (particles[i].vz-(vzj-vzo))
		+ prefacij * dfdx * (particles[i].vz-(vzj-vzo));
	    double dterm7zdy = dprefacijdy * f * (particles[i].vz-(vzj-vzo))
		+ prefacij * dfdy * (particles[i].vz-(vzj-vzo));
	    double dterm7zdz = dprefacijdz * f * (particles[i].vz-(vzj-vzo))
		+ prefacij * dfdz * (particles[i].vz-(vzj-vzo));

	    double dterm7zdvx = prefacij * dfdvx * (particles[i].vz-(vzj-vzo));

	    double dterm7zdvy = prefacij * dfdvy * (particles[i].vz-(vzj-vzo));

	    double dterm7zdvz = prefacij * dfdvz * (particles[i].vz-(vzj-vzo))
		+ prefacij * f;
	    
	    term7x_sum += term7x;
	    term7y_sum += term7y;
	    term7z_sum += term7z;

	    dterm7x_sumdx += dterm7xdx;
	    dterm7x_sumdy += dterm7xdy;
	    dterm7x_sumdz += dterm7xdz;	    
	    dterm7x_sumdvx += dterm7xdvx;
	    dterm7x_sumdvy += dterm7xdvy;
	    dterm7x_sumdvz += dterm7xdvz;	    

	    dterm7y_sumdx += dterm7ydx;
	    dterm7y_sumdy += dterm7ydy;
	    dterm7y_sumdz += dterm7ydz;
	    dterm7y_sumdvx += dterm7ydvx;
	    dterm7y_sumdvy += dterm7ydvy;
	    dterm7y_sumdvz += dterm7ydvz;

	    dterm7z_sumdx += dterm7zdx;
	    dterm7z_sumdy += dterm7zdy;
	    dterm7z_sumdz += dterm7zdz;	    
	    dterm7z_sumdvx += dterm7zdvx;
	    dterm7z_sumdvy += dterm7zdvy;
	    dterm7z_sumdvz += dterm7zdvz;	    
	    

	    double term0 = 0.0;
	    double dterm0dx = 0.0;
	    double dterm0dy = 0.0;
	    double dterm0dz = 0.0;	    

	    double term1 = 0.0;
	    double dterm1dx = 0.0;
	    double dterm1dy = 0.0;
	    double dterm1dz = 0.0;	    
	    double dterm1dvx = 0.0;
	    double dterm1dvy = 0.0;
	    double dterm1dvz = 0.0;	    

	    axj = 0.0;
	    ayj = 0.0;
	    azj = 0.0;	    
	    
	    for (int k=0; k<ASSIST_BODY_NPLANETS; k++){

		// Get position and mass of massive body k.
		assist_all_ephem(ephem, assist->ephem_cache, k, t, &GMk,
			  &xk, &yk, &zk,
			  &vxk, &vyk, &vzk,
			  &axk, &ayk, &azk);

		// Compute position vector of test particle i relative to massive body k.
		const double dxik = particles[i].x + (xo - xk); 
		const double dyik = particles[i].y + (yo - yk);
		const double dzik = particles[i].z + (zo - zk);
		const double rik2 = dxik*dxik + dyik*dyik + dzik*dzik;
		const double _rik  = sqrt(rik2);

		// keep track of GM/rik sum
		term0 += GMk/_rik;

		dterm0dx -= GMk/(_rik*_rik*_rik) * dxik;
		dterm0dy -= GMk/(_rik*_rik*_rik) * dyik;
		dterm0dz -= GMk/(_rik*_rik*_rik) * dzik;				

		if(k != j){
		    // Compute position vector of massive body j relative to massive body k.
		    const double dxjk = xj - xk;
		    const double dyjk = yj - yk;
		    const double dzjk = zj - zk;
		    const double rjk2 = dxjk*dxjk + dyjk*dyjk + dzjk*dzjk;
		    const double _rjk  = sqrt(rjk2);

		    // keep track of GM/rjk sum
		    term1 += GMk/_rjk;

		    axj -= GMk*dxjk/(_rjk*_rjk*_rjk);
		    ayj -= GMk*dyjk/(_rjk*_rjk*_rjk);
		    azj -= GMk*dzjk/(_rjk*_rjk*_rjk);		    		    

		}

	    }

	    term0 *= -2*(beta+gamma)*over_C2;
	    dterm0dx *= -2*(beta+gamma)*over_C2;
	    dterm0dy *= -2*(beta+gamma)*over_C2;
	    dterm0dz *= -2*(beta+gamma)*over_C2;	    	    
	    
	    term1 *= -(2*beta-1)*over_C2;

	    const double rijdotaj = dxij*(axj-axo) + dyij*(ayj-ayo) + dzij*(azj-azo);
	    const double term6 = -0.5*over_C2*rijdotaj;
	    const double dterm6dx = -0.5*over_C2*(axj-axo);
	    const double dterm6dy = -0.5*over_C2*(ayj-ayo);	    
	    const double dterm6dz = -0.5*over_C2*(azj-azo);
	    
	    double term8x = GMj*axj/_rij*(3+4*gamma)/2;
	    double dterm8xdx = -GMj*axj/(_rij*_rij*_rij)*dxij*(3+4*gamma)/2;
	    double dterm8xdy = -GMj*axj/(_rij*_rij*_rij)*dyij*(3+4*gamma)/2;
	    double dterm8xdz = -GMj*axj/(_rij*_rij*_rij)*dzij*(3+4*gamma)/2;	    	    

	    double term8y = GMj*ayj/_rij*(3+4*gamma)/2;
	    double dterm8ydx = -GMj*ayj/(_rij*_rij*_rij)*dxij*(3+4*gamma)/2;
	    double dterm8ydy = -GMj*ayj/(_rij*_rij*_rij)*dyij*(3+4*gamma)/2;
	    double dterm8ydz = -GMj*ayj/(_rij*_rij*_rij)*dzij*(3+4*gamma)/2;	    	    

	    double term8z = GMj*azj/_rij*(3+4*gamma)/2;
	    double dterm8zdx = -GMj*azj/(_rij*_rij*_rij)*dxij*(3+4*gamma)/2;
	    double dterm8zdy = -GMj*azj/(_rij*_rij*_rij)*dyij*(3+4*gamma)/2;
	    double dterm8zdz = -GMj*azj/(_rij*_rij*_rij)*dzij*(3+4*gamma)/2;	    	    

	    term8x_sum += term8x;
	    term8y_sum += term8y;
	    term8z_sum += term8z;

	    dterm8x_sumdx += dterm8xdx;
	    dterm8x_sumdy += dterm8xdy;
	    dterm8x_sumdz += dterm8xdz;	    

	    dterm8y_sumdx += dterm8ydx;
	    dterm8y_sumdy += dterm8ydy;
	    dterm8y_sumdz += dterm8ydz;

	    dterm8z_sumdx += dterm8zdx;
	    dterm8z_sumdy += dterm8zdy;
	    dterm8z_sumdz += dterm8zdz;	    
	    
	    double factor = term0 + term1 + term2 + term3 + term4 + term5 + term6;

	    double dfactordx = dterm0dx + dterm1dx + dterm5dx + dterm6dx;
	    double dfactordy = dterm0dy + dterm1dy + dterm5dy + dterm6dy;
	    double dfactordz = dterm0dz + dterm1dz + dterm5dz + dterm6dz;	    
	    double dfactordvx = dterm1dvx + dterm2dvx + dterm4dvx;
	    double dfactordvy = dterm1dvy + dterm2dvy + dterm4dvy;
	    double dfactordvz = dterm1dvz + dterm2dvz + dterm4dvz;

	    if(eih_file){
		fprintf(eih_file, "%24.16lE ", -factor*C2);
		fprintf(eih_file, "%24.16lE %24.16lE %24.16lE %24.16lE ",
			-factor*C2*prefacij*dxij,
			-factor*C2*prefacij*dyij,
			-factor*C2*prefacij*dzij,
			f);	    
		fprintf(eih_file, "%24.16lE %24.16lE %24.16lE ",
			prefacij*f*(particles[i].vx-(vxj-vxo)),
			prefacij*f*(particles[i].vy-(vyj-vyo)),
			prefacij*f*(particles[i].vz-(vzj-vzo)));	    
		fprintf(eih_file, "%24.16lE %24.16lE %24.16lE ",
			term8x,
			term8y,
			term8z);
		fprintf(eih_file, "%24.16lE %24.16lE %24.16lE\n",
			axj, ayj, azj);

		fflush(eih_file);
	    }

	    grx += -prefacij*dxij*factor;
	    gry += -prefacij*dyij*factor;
	    grz += -prefacij*dzij*factor;
	    
	    particles[i].ax += -prefacij*dxij*factor;
	    particles[i].ay += -prefacij*dyij*factor;
	    particles[i].az += -prefacij*dzij*factor;

	    // Variational equation terms go here.

	    dxdx += -dprefacijdx*dxij*factor
		-prefacij*factor
		-prefacij*dxij*dfactordx;
	    
	    dxdy += -dprefacijdy*dxij*factor
		-prefacij*dxij*dfactordy;
	    
	    dxdz += -dprefacijdz*dxij*factor
		-prefacij*dxij*dfactordz;

	    dxdvx += 
		-prefacij*dxij*dfactordvx;

	    dxdvy += 
		-prefacij*dxij*dfactordvy;

	    dxdvz += 
		-prefacij*dxij*dfactordvz;

	    dydx += -dprefacijdx*dyij*factor
		-prefacij*dyij*dfactordx;
	    
	    dydy += -dprefacijdy*dyij*factor
		-prefacij*factor
		-prefacij*dyij*dfactordy;
	    
	    dydz += -dprefacijdz*dyij*factor
		-prefacij*dyij*dfactordz;

	    dydvx += 
		-prefacij*dyij*dfactordvx;

	    dydvy += 
		-prefacij*dyij*dfactordvy;

	    dydvz += 
		-prefacij*dyij*dfactordvz;
	    
	    dzdx += -dprefacijdx*dzij*factor
		-prefacij*dzij*dfactordx;

	    dzdy += -dprefacijdy*dzij*factor
		-prefacij*dzij*dfactordy;
	    
	    dzdz += -dprefacijdz*dzij*factor
		-prefacij*factor
		-prefacij*dzij*dfactordz;

	    dzdvx += 
		-prefacij*dzij*dfactordvx;

	    dzdvy += 
		-prefacij*dzij*dfactordvy;

	    dzdvz += 
		-prefacij*dzij*dfactordvz;

        }

	grx += term7x_sum*over_C2 + term8x_sum*over_C2;
	gry += term7y_sum*over_C2 + term8y_sum*over_C2;
	grz += term7z_sum*over_C2 + term8z_sum*over_C2;

	if(outfile){
	    fprintf(outfile, "%3s %25.16le %25.16le %25.16le %25.16le\n", "GR", jd_ref+t,
		    grx, gry, grz);
	    fflush(outfile);
	}

	dxdx += dterm7x_sumdx*over_C2 + dterm8x_sumdx*over_C2;
	dxdy += dterm7x_sumdy*over_C2 + dterm8x_sumdy*over_C2;
	dxdz += dterm7x_sumdz*over_C2 + dterm8x_sumdz*over_C2;	
	dxdvx += dterm7x_sumdvx*over_C2;
	dxdvy += dterm7x_sumdvy*over_C2;
	dxdvz += dterm7x_sumdvz*over_C2;

	dydx += dterm7y_sumdx*over_C2 + dterm8y_sumdx*over_C2;
	dydy += dterm7y_sumdy*over_C2 + dterm8y_sumdy*over_C2;
	dydz += dterm7y_sumdz*over_C2 + dterm8y_sumdz*over_C2;	
	dydvx += dterm7y_sumdvx*over_C2;
	dydvy += dterm7y_sumdvy*over_C2;
	dydvz += dterm7y_sumdvz*over_C2;

	dzdx += dterm7z_sumdx*over_C2 + dterm8z_sumdx*over_C2;
	dzdy += dterm7z_sumdy*over_C2 + dterm8z_sumdy*over_C2;
	dzdz += dterm7z_sumdz*over_C2 + dterm8z_sumdz*over_C2;	
	dzdvx += dterm7z_sumdvx*over_C2;
	dzdvy += dterm7z_sumdvy*over_C2;
	dzdvz += dterm7z_sumdvz*over_C2;
	
	particles[i].ax += term7x_sum*over_C2 + term8x_sum*over_C2;
	particles[i].ay += term7y_sum*over_C2 + term8y_sum*over_C2;
	particles[i].az += term7z_sum*over_C2 + term8z_sum*over_C2;

	// Variational equation terms go here.
	for (int v=0; v < sim->var_config_N; v++){
	    struct reb_variational_configuration const vc = sim->var_config[v];
	    int tp = vc.testparticle;
	    struct reb_particle* const particles_var1 = particles + vc.index;		
	    if(tp == i){
	    
		// variational particle coords
		const double ddx = particles_var1[0].x;
		const double ddy = particles_var1[0].y;
		const double ddz = particles_var1[0].z;
		const double ddvx = particles_var1[0].vx;
		const double ddvy = particles_var1[0].vy;
		const double ddvz = particles_var1[0].vz;

		// Matrix multiplication
		const double dax =   ddx  * dxdx  + ddy  * dxdy  + ddz  * dxdz
		    +   ddvx * dxdvx + ddvy * dxdvy + ddvz * dxdvz;
		const double day =   ddx  * dydx  + ddy  * dydy  + ddz  * dydz
		    +   ddvx * dydvx + ddvy * dydvy + ddvz * dydvz;
		const double daz =   ddx  * dzdx  + ddy  * dzdy  + ddz  * dzdz
		    +   ddvx * dzdvx + ddvy * dzdvy + ddvz * dzdvz;

		// Accumulate acceleration terms
		particles_var1[0].ax += dax;
		particles_var1[0].ay += day;
		particles_var1[0].az += daz;
		
	    }
	}
    }
}
