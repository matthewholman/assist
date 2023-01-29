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


enum {
    NO_ERR,        // no error
    ERR_JPL_EPHEM, // JPL ephemeris file not found
    ERR_JPL_AST,   // JPL asteroid file not found
    ERR_NAST,      // asteroid number out of range
    ERR_NEPH,      // planet number out of range
};

// Forward function declarations
static void assist_additional_force_direct(struct reb_simulation* sim, const struct assist_cache_item center, FILE *outfile);
static void assist_additional_force_solar_J2(struct reb_simulation* sim, const struct assist_cache_item center, FILE *outfile);
static void assist_additional_force_earth_J2J4(struct reb_simulation* sim, const struct assist_cache_item center, FILE *outfile);
static void assist_additional_force_non_gravitational(struct reb_simulation* sim, const struct assist_cache_item center, FILE *outfile);
static void assist_additional_force_potential_GR(struct reb_simulation* sim, const struct assist_cache_item center,  FILE *outfile);
static void assist_additional_force_simple_GR(struct reb_simulation* sim, const struct assist_cache_item center, FILE *outfile);
static void assist_additional_force_eih_GR(struct reb_simulation* sim, int eih_loop_limit, const struct assist_cache_item center, FILE *outfile, FILE *eih_file);

static const int N_ephem = 11;
static const int N_ast = 16;

void assist_additional_forces(struct reb_simulation* sim){

    // implement additional_forces here

    struct assist_extras* assist = (struct assist_extras*) sim->extras;
    struct assist_ephem* ephem = assist->ephem;
    int geo = assist->geocentric;

    const unsigned int N = sim->N;  // N includes real+variational particles
    const unsigned int N_real = N - sim->N_var;    

    // The limit of the EIH GR limit should be a free
    // parameter
    int eih_loop_limit = N_ephem; // 1;

    const double t = sim->t;

    struct assist_cache_item center = {0};

    // Check which center is used.
    // The current options are the barycenter (default) and geocenter.
    // We might consider adding the heliocenter.

    if(geo == 1){
	// geocentric
	// Get mass, position, velocity, and acceleration of the Earth for later use.
	// The offset position is used to adjust the particle positions.
	int flag = assist_all_ephem(ephem, assist->ephem_cache, 3, t, &center);
	if(flag != NO_ERR){
	    char outstring[50];
	    sprintf(outstring, "%s %d %d\n", "Ephemeris error a ", 3, flag);
	    reb_error(sim, outstring);
	}
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
        assist_additional_force_non_gravitational(sim, center, outfile);
    }
    if (assist->forces & ASSIST_FORCE_EARTH_HARMONICS){
        assist_additional_force_earth_J2J4(sim, center, outfile);
    }
    if (assist->forces & ASSIST_FORCE_SUN_HARMONICS){
        assist_additional_force_solar_J2(sim, center, outfile);
    }
    
    FILE *eih_file = NULL;
    // Uncomment this line and recompile for testing.
    //eih_file = fopen("eih_acc.out", "w");

    if (assist->forces & ASSIST_FORCE_GR_EIH){
        assist_additional_force_eih_GR(sim, eih_loop_limit, center, outfile, eih_file);
    }

    if (assist->forces & ASSIST_FORCE_GR_POTENTIAL){
        assist_additional_force_potential_GR(sim, center, outfile);
    }
    if (assist->forces & ASSIST_FORCE_GR_SIMPLE){
        assist_additional_force_simple_GR(sim, center, outfile);
    }
    if (assist->forces & (ASSIST_FORCE_SUN | ASSIST_FORCE_PLANETS | ASSIST_FORCE_ASTEROIDS)){
        assist_additional_force_direct(sim, center, outfile);
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
        assist_all_ephem(ephem, assist->ephem_cache, 3, t, &center);

        // This is the indirect term for geocentric equations
        // of motion.
        for (int j=0; j<N_real; j++){    
            sim->particles[j].ax -= center.ax;
            sim->particles[j].ay -= center.ay;
            sim->particles[j].az -= center.az;
        }
    }
}

static int planet_ephem(struct assist_ephem* ephem, const int i, const double jd_ref, const double t, struct assist_cache_item* const result){

    struct mpos_s now;

    // Calculate GM values for Earth and Moon
    // from Earth-moon ratio and sum.
#define    em_r JPL_EPHEM_EMRAT
#define    GMe (em_r/(1.+em_r) * JPL_EPHEM_GMB)
#define    GMm (1./(1.+em_r) * JPL_EPHEM_GMB)

    // The values below are G*mass.
    // Units are solar masses, au, days.
    // DE440/441 units: au^3 day^-2.
    const static double JPL_GM[11] =
	{
	    JPL_EPHEM_GMS, // 0 sun
	    JPL_EPHEM_GM1, // 1 mercury
	    JPL_EPHEM_GM2, // 2 venus
	    GMe,           // 3 earth
	    GMm,           // 4 moon
	    JPL_EPHEM_GM4, // 5 mars
	    JPL_EPHEM_GM5, // 6 jupiter
	    JPL_EPHEM_GM6, // 7 saturn
	    JPL_EPHEM_GM7, // 8 uranus
	    JPL_EPHEM_GM8, // 9 neptune
	    JPL_EPHEM_GM9, // 10 pluto
	};

    // Get position, velocity, and mass of body i in barycentric coords.

    result->m = JPL_GM[i];

    assist_jpl_calc(ephem->pl, &now, jd_ref, t, i); 

    // Convert to au/day and au/day^2
    vecpos_div(now.u, ephem->pl->cau);
    vecpos_div(now.v, ephem->pl->cau/86400.);
    vecpos_div(now.w, ephem->pl->cau/(86400.*86400.));

    result->x = now.u[0];
    result->y = now.u[1];
    result->z = now.u[2];
    result->vx = now.v[0];
    result->vy = now.v[1];
    result->vz = now.v[2];
    result->ax = now.w[0];
    result->ay = now.w[1];
    result->az = now.w[2];

    return(NO_ERR);
    
}

static int ast_ephem(struct assist_ephem* ephem, const int i, const double jd_ref, const double t, struct assist_cache_item* const result){

    //static int initialized = 0;

    //static struct spk_s *spl;
    struct mpos_s pos;

    // DE441
    // The values below are G*mass.
    // Units are solar masses, au, days.
    const static double JPL_GM[16] =    
    {
	JPL_EPHEM_MA0107, // 107 camilla
	JPL_EPHEM_MA0001, // 1 Ceres
	JPL_EPHEM_MA0065, // 65 cybele
	JPL_EPHEM_MA0511, // 511 davida
	JPL_EPHEM_MA0015, // 15 eunomia
	JPL_EPHEM_MA0031, // 31 euphrosyne	    
	JPL_EPHEM_MA0052, // 52 europa
	JPL_EPHEM_MA0010, // 10 hygiea
	JPL_EPHEM_MA0704, // 704 interamnia
	JPL_EPHEM_MA0007, // 7 iris
	JPL_EPHEM_MA0003, // 3 juno
	JPL_EPHEM_MA0002, // 2 pallas
	JPL_EPHEM_MA0016, // 16 psyche
	JPL_EPHEM_MA0087, // 87 sylvia
	JPL_EPHEM_MA0088, // 88 thisbe
	JPL_EPHEM_MA0004  // 4 vesta
    };

    if(i<0 || i>15){
	return(ERR_NAST);
    }

    if(ephem->spl==NULL){
	return(ERR_JPL_EPHEM);	
    }

    /*
    if (initialized == 0){
	char buf[] = "/Users/mholman/assist/data/sb441-n16.bsp";
	if ((spl = assist_spk_init(buf)) == NULL) {
	    printf("Couldn't find asteroid ephemeris file: %s\n", buf);
	    return(ERR_JPL_AST);
	}

	initialized = 1;

    }
    */
    
    // TODO: again, the units might be handled more
    // generally

    result->m = JPL_GM[i];

    assist_spk_calc(ephem->spl, i, jd_ref, t, &pos);

    result->x = pos.u[0];
    result->y = pos.u[1];
    result->z = pos.u[2];

    return(NO_ERR);

}


int assist_all_ephem(struct assist_ephem* ephem, struct assist_ephem_cache* ephem_cache, const int i, const double t, struct assist_cache_item* result){
    if (ephem_cache){
        const double* const ct = ephem_cache->t+7*i;
        for (int s=0; s<7; s+=1){
            if (t==ct[s]){
                *result = *(ephem_cache->items+i*7+s);
                return NO_ERR;
            }
        }
        // No match
    }

    const double jd_ref = ephem->jd_ref;

    // Get position and mass of massive body i.
    if(i < N_ephem){
        int flag = planet_ephem(ephem, i, jd_ref, t, result);
        if(flag != NO_ERR) return(flag);
    }else{
        // Get position and mass of asteroid i-N_ephem.
        int flag = ast_ephem(ephem, i-N_ephem, jd_ref, t, result);
        if(flag != NO_ERR) return(flag);

        struct assist_cache_item item;
        flag = assist_all_ephem(ephem, ephem_cache, 0, t, &item);
        if(flag != NO_ERR) return(flag);		    

        // Translate massive asteroids from heliocentric to barycentric.
        result->x += item.x; result->y += item.y; result->z += item.z;
        // velocities and accelerations are not needed for these
        // bodies
        //*vx = NAN; *vy = NAN; *vz = NAN;
        //*ax = NAN; *ay = NAN; *az = NAN;

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
        *(ephem_cache->items + i*7 +os) = *result;
    }
    return(NO_ERR);
}


static void assist_additional_force_direct(struct reb_simulation* sim, const struct assist_cache_item center, FILE *outfile){
    const unsigned int N = sim->N;  // N includes real+variational particles
    const unsigned int N_real = N - sim->N_var;

    struct assist_extras* assist = (struct assist_extras*) sim->extras;
    struct assist_ephem* ephem = assist->ephem;
    const double jd_ref = ephem->jd_ref;
    
    const double t = sim->t;    

    const int N_tot = N_ephem + N_ast;

    struct reb_particle* const particles = sim->particles;

    static const int order[] = {13, 16, 20, 25, 11, 23, 21, 15,
				24, 17, 19, 14, 18, 22, 26, 12,
				10,  4,  5,  1,  9,  8,  3,  2,  7, 6, 0};

    // Direct forces from massives bodies
    for (int k=0; k<N_tot; k++){
        int i = order[k];
        if (i==0 && !(assist->forces & ASSIST_FORCE_SUN)) continue;
        if (i>=1 && i<N_ephem && !(assist->forces & ASSIST_FORCE_PLANETS)) continue;
        if (i>=N_ephem && !(assist->forces & ASSIST_FORCE_ASTEROIDS)) continue;

        // Get position and mass of massive body i.
        // TOOD: make a version that returns the positions, velocities,
        // and accelerations for all the bodies at a given time.

        struct assist_cache_item item;
        int flag = assist_all_ephem(ephem, assist->ephem_cache, i, t, &item);

        if(flag != NO_ERR){
            char outstring[50];
            sprintf(outstring, "%s %d %d\n", "Ephemeris error b ", i, flag);	    
            reb_error(sim, outstring);
        }

        // Loop over test particles
        for (int j=0; j<N_real; j++){

            // Compute position vector of test particle j relative to massive body i.
            const double dx = particles[j].x + (center.x - item.x); 
            const double dy = particles[j].y + (center.y - item.y);
            const double dz = particles[j].z + (center.z - item.z);
            const double r2 = dx*dx + dy*dy + dz*dz;
            const double _r  = sqrt(r2);
            const double prefac = item.m/(_r*_r*_r);

            if(outfile){
                fprintf(outfile, "%3d %25.16le %25.16le %25.16le %25.16le %25.16le %25.16le %25.16le %25.16le\n", i, jd_ref+t, item.m, dx, dy, dz, -prefac*dx, -prefac*dy, -prefac*dz);
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

    for (int k=0; k<N_tot; k++){
	//for (int i=0; i<N_tot; i++){
	int i = order[k];
	//int i = k;

        // Get position and mass of massive body i.	
        struct assist_cache_item item;
	int flag = assist_all_ephem(ephem, assist->ephem_cache, i, t, &item);

	if(flag != NO_ERR){
	    char outstring[50];
	    sprintf(outstring, "%s %d %d\n", "Ephemeris error c ", i, flag);	    	    
	    reb_error(sim, outstring);
	}

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
	    const double dx = particles[j].x + (center.x - item.x);
	    const double dy = particles[j].y + (center.x - item.y);
	    const double dz = particles[j].z + (center.x - item.z);
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
		    particles_var1[0].ax += item.m * dax; 
		    particles_var1[0].ay += item.m * day; 
		    particles_var1[0].az += item.m * daz; 

		}
	    }
        }
    }
}

static void assist_additional_force_earth_J2J4(struct reb_simulation* sim, const struct assist_cache_item center, FILE *outfile){

    struct assist_extras* assist = (struct assist_extras*) sim->extras;
    struct assist_ephem* ephem = assist->ephem;
    const double jd_ref = ephem->jd_ref;
    
    //const double G = sim->G;
    const unsigned int N = sim->N;  // N includes real+variational particles
    const unsigned int N_real = N - sim->N_var;

    const double t = sim->t;    

    struct reb_particle* const particles = sim->particles;

    // 
    // Here is the treatment of the Earth's J2, J3, and J4.
    // J2 and J4 are based in part on gravitational_harmonics
    // example in reboundx.
    // Assumes the coordinates are geocentric.
    // Also assuming that Earth's pole is along the z
    // axis.  This is only precisely true at the J2000
    // epoch.

    // The geocenter is the reference for the Earth J2/J4 calculations.
    struct assist_cache_item iteme;
    assist_all_ephem(ephem, assist->ephem_cache, 3, t, &iteme);

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
        double dx = p.x + (center.x - iteme.x);
        double dy = p.y + (center.y - iteme.y);
        double dz = p.z + (center.z - iteme.z);

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

	double resx = iteme.m*J2e_prefac*J2e_fac*dx;
	double resy = iteme.m*J2e_prefac*J2e_fac*dy;
	double resz = iteme.m*J2e_prefac*(J2e_fac-2.)*dz;	

	// J3 terms
        const double J3e_prefac = 5.*J3e*Re_eq*Re_eq*Re_eq/r2/r2/r/2.;
        const double J3e_fac = 3.-7.*costheta2;
	
	resx += -iteme.m*J3e_prefac*(1./r2)*J3e_fac*dx*dz;
    resy += -iteme.m*J3e_prefac*(1./r2)*J3e_fac*dy*dz;
	resz += -iteme.m*J3e_prefac*(6.*costheta2 - 7.*costheta2*costheta2-0.6);
	
	// J4 terms
        const double J4e_prefac = 5.*J4e*Re_eq*Re_eq*Re_eq*Re_eq/r2/r2/r2/r/8.;
        const double J4e_fac = 63.*costheta2*costheta2-42.*costheta2 + 3.;

        resx += iteme.m*J4e_prefac*J4e_fac*dx;
        resy += iteme.m*J4e_prefac*J4e_fac*dy;
        resz += iteme.m*J4e_prefac*(J4e_fac+12.-28.*costheta2)*dz;

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

	const double dxdx = iteme.m*J2e_prefac*(J2e_fac-5.*J2e_fac2*dx*dx/r2);
	const double dydy = iteme.m*J2e_prefac*(J2e_fac-5.*J2e_fac2*dy*dy/r2);
	const double dzdz = iteme.m*J2e_prefac*(-1.)*J2e_fac3;
	const double dxdy = iteme.m*J2e_prefac*(-5.)*J2e_fac2*dx*dy/r2;
	const double dydz = iteme.m*J2e_prefac*(-5.)*(J2e_fac2-2.)*dy*dz/r2;
	const double dxdz = iteme.m*J2e_prefac*(-5.)*(J2e_fac2-2.)*dx*dz/r2;

	// J3 terms
        const double costheta = dz/r;
        const double J3e_fac2 = 21*(-3.*costheta2+1.)/r2;
        const double J3e_fac3 = 3*(-21.*costheta2*costheta2+14.*costheta2-1.)/r2;
        const double J3e_fac4 = (-63.*costheta2*costheta2+70.*costheta2-15.)*costheta/r;

	const double dxdxJ3 = iteme.m*J3e_prefac*costheta*(J3e_fac2*dx*dx-J3e_fac)/r;
	const double dydyJ3 = iteme.m*J3e_prefac*costheta*(J3e_fac2*dy*dy-J3e_fac)/r;
	const double dzdzJ3 = iteme.m*J3e_prefac*J3e_fac4;
	const double dxdyJ3 = iteme.m*J3e_prefac*J3e_fac2*costheta*dx*dy/r;
	const double dydzJ3 = iteme.m*J3e_prefac*J3e_fac3*dy;
	const double dxdzJ3 = iteme.m*J3e_prefac*J3e_fac3*dx;

	// J4 terms
        const double J4e_fac2= 33.*costheta2*costheta2-18.*costheta2 + 1.;
        const double J4e_fac3= 33.*costheta2*costheta2-30.*costheta2 + 5.;
        const double J4e_fac4= 231.*costheta2*costheta2*costheta2-315.*costheta2*costheta2+105.*costheta2 - 5.;
	
	const double dxdxJ4 = iteme.m*J4e_prefac*(J4e_fac-21.*J4e_fac2*dx*dx/r2);
	const double dydyJ4 = iteme.m*J4e_prefac*(J4e_fac-21.*J4e_fac2*dy*dy/r2);
	const double dzdzJ4 = iteme.m*J4e_prefac*(-3.)*J4e_fac4;
	const double dxdyJ4 = iteme.m*J4e_prefac*(-21.)*J4e_fac2*dx*dy/r2;
	const double dydzJ4 = iteme.m*J4e_prefac*(-21.)*J4e_fac3*dy*dz/r2;
	const double dxdzJ4 = iteme.m*J4e_prefac*(-21.)*J4e_fac3*dx*dz/r2;

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

static void assist_additional_force_solar_J2(struct reb_simulation* sim, const struct assist_cache_item center, FILE *outfile){

    struct assist_extras* assist = (struct assist_extras*) sim->extras;
    struct assist_ephem* ephem = assist->ephem;
    const double jd_ref = ephem->jd_ref;
    
    const unsigned int N = sim->N;  // N includes real+variational particles
    const unsigned int N_real = N - sim->N_var;

    const double t = sim->t;    

    struct reb_particle* const particles = sim->particles;

    // The Sun center is reference for these calculations.

    struct assist_cache_item itemr;

    assist_all_ephem(ephem, assist->ephem_cache, 0, t, &itemr); 
    const double GMsun = itemr.m;    

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
        double dx = p.x + (center.x - itemr.x);
        double dy = p.y + (center.y - itemr.y);
        double dz = p.z + (center.z - itemr.z);

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

static void assist_additional_force_non_gravitational(struct reb_simulation* sim, const struct assist_cache_item center, FILE *outfile){

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

    // Here is the treatment of non-gravitational forces.

    // The Sun center is reference for these calculations.
    struct assist_cache_item itemr;
    assist_all_ephem(ephem, assist->ephem_cache, 0, t, &itemr);

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
        double dx = p.x + (center.x - itemr.x);
        double dy = p.y + (center.y - itemr.y);
        double dz = p.z + (center.z - itemr.z);

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

	double dvx = p.vx + (center.vx - itemr.vx);
	double dvy = p.vy + (center.vy - itemr.vy);
	double dvz = p.vz + (center.vz - itemr.vz);

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

static void assist_additional_force_potential_GR(struct reb_simulation* sim, const struct assist_cache_item center,  FILE *outfile){
    
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
    
    // The Sun center is reference for these calculations.
    struct assist_cache_item items;
    assist_all_ephem(ephem, assist->ephem_cache, 0, t, &items);
    const double GMsun = items.m;

    //xs = ys = zs = 0.0;

    for (int j=0; j<N_real; j++){

        struct reb_particle p = particles[j];

	p.x += (center.x - items.x);
	p.y += (center.y - items.y);
	p.z += (center.z - items.z);

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

static void assist_additional_force_simple_GR(struct reb_simulation* sim, struct assist_cache_item center, FILE *outfile){
    
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


    
    // The Sun center is reference for these calculations.
    struct assist_cache_item items;
    assist_all_ephem(ephem, assist->ephem_cache, 0, t, &items);
    const double GMsun = items.m;

    for (int j=0; j<N_real; j++){

        struct reb_particle p = particles[j];

	p.x += (center.x - items.x);
	p.y += (center.y - items.y);
	p.z += (center.z - items.z);
	p.vx += (center.vx - items.vx);
	p.vy += (center.vy - items.vy);
	p.vz += (center.vz - items.vz);
	
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
	    const struct assist_cache_item center,
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

	double term7x_sum = 0.0;
	double term7y_sum = 0.0;
	double term7z_sum = 0.0;
	double term8x_sum = 0.0;
	double term8y_sum = 0.0;
	double term8z_sum = 0.0;

	double grx = 0.0;
	double gry = 0.0;
	double grz = 0.0;		

	for (int j=0; j<eih_loop_limit; j++){ // This is either 1 or N_ephem

	    // Get position and mass of massive body j.
        struct assist_cache_item itemj;
	    assist_all_ephem(ephem, assist->ephem_cache, j, t, &itemj);

	    // Compute position vector of test particle i relative to massive body j.
	    const double dxij = particles[i].x + (center.x - itemj.x); 
	    const double dyij = particles[i].y + (center.y - itemj.y);
	    const double dzij = particles[i].z + (center.z - itemj.z);
	    const double rij2 = dxij*dxij + dyij*dyij + dzij*dzij;
	    const double _rij  = sqrt(rij2);
	    const double prefacij = itemj.m/(rij2*_rij);

	    // This is the place to do all the various i-j dot products
	    
	    const double vi2 = particles[i].vx*particles[i].vx +
		particles[i].vy*particles[i].vy +
		particles[i].vz*particles[i].vz;

	    const double term2 = gamma*over_C2*vi2;

	    const double vj2 = (itemj.vx-center.vx)*(itemj.vx-center.vx) + (itemj.vy-center.vy)*(itemj.vy-center.vy) + (itemj.vz-center.vz)*(itemj.vz-center.vz);

	    const double term3 = (1+gamma)*over_C2*vj2;
	    // Variational equations do not depend on term3

	    const double vidotvj = particles[i].vx*(itemj.vx-center.vx) +
		particles[i].vy*(itemj.vy-center.vy) +
		particles[i].vz*(itemj.vz-center.vz);

	    const double term4 = -2*(1+gamma)*over_C2*vidotvj;

	    const double rijdotvj = dxij*(itemj.vx-center.vx) + dyij*(itemj.vy-center.vy) + dzij*(itemj.vz-center.vz);

	    if(eih_file){
		fprintf(eih_file, " EIH_J%12d\n", j);	    
		fprintf(eih_file, "%25.16lE ", rijdotvj/_rij);
	    }

	    const double term5 = -1.5*over_C2*(rijdotvj*rijdotvj)/(_rij*_rij);

	    const double fx = (2+2*gamma)*particles[i].vx - (1+2*gamma)*(itemj.vx-center.vx);
	    const double fy = (2+2*gamma)*particles[i].vy - (1+2*gamma)*(itemj.vy-center.vy);
	    const double fz = (2+2*gamma)*particles[i].vz - (1+2*gamma)*(itemj.vz-center.vz);
	    const double f = dxij*fx + dyij*fy + dzij*fz;

	    const double prefacij_f = prefacij*f;
	    const double term7x = prefacij_f*(particles[i].vx-(itemj.vx-center.vx));
	    const double term7y = prefacij_f*(particles[i].vy-(itemj.vy-center.vy));
	    const double term7z = prefacij_f*(particles[i].vz-(itemj.vz-center.vz));
	    
	    term7x_sum += term7x;
	    term7y_sum += term7y;
	    term7z_sum += term7z;

	    double term0 = 0.0;
	    double term1 = 0.0;

	    double axj = 0.0;
	    double ayj = 0.0;
	    double azj = 0.0;	    
	    
	    for (int k=0; k<N_ephem; k++){

		// Get position and mass of massive body k.
        struct assist_cache_item itemk;
	    assist_all_ephem(ephem, assist->ephem_cache, k, t, &itemk);

		// Compute position vector of test particle i relative to massive body k.
		const double dxik = particles[i].x + (center.x - itemk.x); 
		const double dyik = particles[i].y + (center.y - itemk.y);
		const double dzik = particles[i].z + (center.z - itemk.z);
		const double rik2 = dxik*dxik + dyik*dyik + dzik*dzik;
		const double _rik  = sqrt(rik2);

		// keep track of GM/rik sum
		term0 += itemk.m/_rik;

		if(k != j){
		    // Compute position vector of massive body j relative to massive body k.
		    const double dxjk = itemj.x - itemk.x;
		    const double dyjk = itemj.y - itemk.y;
		    const double dzjk = itemj.z - itemk.z;
		    const double rjk2 = dxjk*dxjk + dyjk*dyjk + dzjk*dzjk;
		    const double _rjk  = sqrt(rjk2);

		    // keep track of GM/rjk sum
		    term1 += itemk.m/_rjk;

		    const double fac = itemk.m/(rjk2*_rjk);
		    axj -= fac*dxjk;
		    ayj -= fac*dyjk;
		    azj -= fac*dzjk;

		}

	    }

	    term0 *= -2*(beta+gamma)*over_C2;
	    
	    term1 *= -(2*beta-1)*over_C2;

	    const double rijdotaj = dxij*(axj-center.ax) + dyij*(ayj-center.ay) + dzij*(azj-center.az);
	    const double term6 = -0.5*over_C2*rijdotaj;

	    const double term8_fac = itemj.m/_rij*(3+4*gamma)/2;
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
			prefacij*f*(particles[i].vx-(itemj.vx-center.vx)),
			prefacij*f*(particles[i].vy-(itemj.vy-center.vy)),
			prefacij*f*(particles[i].vz-(itemj.vz-center.vz)));	    
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
 
 	for (int j=0; j<eih_loop_limit; j++){ // This is either 1 or N_ephem
 
 	    // Get position and mass of massive body j.
        struct assist_cache_item itemj;
 	    assist_all_ephem(ephem, assist->ephem_cache, j, t, &itemj);
 
 	    // Compute position vector of test particle i relative to massive body j.
 	    const double dxij = particles[i].x + (center.x - itemj.x); 
 	    const double dyij = particles[i].y + (center.y - itemj.y);
 	    const double dzij = particles[i].z + (center.z - itemj.z);
 	    const double rij2 = dxij*dxij + dyij*dyij + dzij*dzij;
 	    const double _rij  = sqrt(rij2);
 	    const double prefacij = itemj.m/(_rij*_rij*_rij);
 
 	    const double dprefacijdx = -3.0*itemj.m/(_rij*_rij*_rij*_rij*_rij)*dxij;
 	    const double dprefacijdy = -3.0*itemj.m/(_rij*_rij*_rij*_rij*_rij)*dyij;
 	    const double dprefacijdz = -3.0*itemj.m/(_rij*_rij*_rij*_rij*_rij)*dzij;
 
 	    // This is the place to do all the various i-j dot products
 	    
 	    const double vi2 = particles[i].vx*particles[i].vx +
 		particles[i].vy*particles[i].vy +
 		particles[i].vz*particles[i].vz;
 
 	    const double term2 = gamma*over_C2*vi2;
 	    const double dterm2dvx = 2.0*gamma*over_C2*particles[i].vx;
 	    const double dterm2dvy = 2.0*gamma*over_C2*particles[i].vy;
 	    const double dterm2dvz = 2.0*gamma*over_C2*particles[i].vz;	    
 
 	    const double vj2 = (itemj.vx-center.vx)*(itemj.vx-center.vx) + (itemj.vy-center.vy)*(itemj.vy-center.vy) + (itemj.vz-center.vz)*(itemj.vz-center.vz);
 
 	    const double term3 = (1+gamma)*over_C2*vj2;
 	    // Variational equations do not depend on term3
 
 	    const double vidotvj = particles[i].vx*(itemj.vx-center.vx) +
 		particles[i].vy*(itemj.vy-center.vy) +
 		particles[i].vz*(itemj.vz-center.vz);
 
 	    const double term4 = -2*(1+gamma)*over_C2*vidotvj;
 	    const double dterm4dvx = -2*(1+gamma)*over_C2*(itemj.vx-center.vx);
 	    const double dterm4dvy = -2*(1+gamma)*over_C2*(itemj.vy-center.vy);
 	    const double dterm4dvz = -2*(1+gamma)*over_C2*(itemj.vz-center.vz);	    	    
 	    
 
 	    const double rijdotvj = dxij*(itemj.vx-center.vx) + dyij*(itemj.vy-center.vy) + dzij*(itemj.vz-center.vz);
 
 	    if(eih_file){
 		fprintf(eih_file, " EIH_J%12d\n", j);	    
 		fprintf(eih_file, "%25.16lE ", rijdotvj/_rij);
 	    }
 
 	    const double term5 = -1.5*over_C2*(rijdotvj*rijdotvj)/(_rij*_rij);
 	    const double term5_fac = 3.0*over_C2*rijdotvj/_rij;
 	    const double dterm5dx = -term5_fac*((itemj.vx-center.vx)/_rij - rijdotvj*dxij/(_rij*_rij*_rij));
 	    const double dterm5dy = -term5_fac*((itemj.vy-center.vy)/_rij - rijdotvj*dyij/(_rij*_rij*_rij));
 	    const double dterm5dz = -term5_fac*((itemj.vz-center.vz)/_rij - rijdotvj*dzij/(_rij*_rij*_rij));	    	    
 
 	    double fx = (2+2*gamma)*particles[i].vx - (1+2*gamma)*(itemj.vx-center.vx);
 	    double fy = (2+2*gamma)*particles[i].vy - (1+2*gamma)*(itemj.vy-center.vy);
 	    double fz = (2+2*gamma)*particles[i].vz - (1+2*gamma)*(itemj.vz-center.vz);
 	    double f = dxij*fx + dyij*fy + dzij*fz;
 
 	    double dfdx = fx;
 	    double dfdy = fy;
 	    double dfdz = fz;	    	    
 	    double dfdvx = dxij*(2+2*gamma);
 	    double dfdvy = dyij*(2+2*gamma);
 	    double dfdvz = dzij*(2+2*gamma);
 
 	    double term7x = prefacij*f*(particles[i].vx-(itemj.vx-center.vx));
 	    double term7y = prefacij*f*(particles[i].vy-(itemj.vy-center.vy));
 	    double term7z = prefacij*f*(particles[i].vz-(itemj.vz-center.vz));
 
 	    double dterm7xdx = dprefacijdx * f * (particles[i].vx-(itemj.vx-center.vx))
 		+ prefacij * dfdx * (particles[i].vx-(itemj.vx-center.vx));
 	    double dterm7xdy = dprefacijdy * f * (particles[i].vx-(itemj.vx-center.vx))
 		+ prefacij * dfdy * (particles[i].vx-(itemj.vx-center.vx));
 	    double dterm7xdz = dprefacijdz * f * (particles[i].vx-(itemj.vx-center.vx))
 		+ prefacij * dfdz * (particles[i].vx-(itemj.vx-center.vx));
 	    double dterm7xdvx = prefacij * dfdvx * (particles[i].vx-(itemj.vx-center.vx))
 		+ prefacij * f;
 	    double dterm7xdvy = prefacij * dfdvy * (particles[i].vx-(itemj.vx-center.vx));
 
 	    double dterm7xdvz = prefacij * dfdvz * (particles[i].vx-(itemj.vx-center.vx));	    
 
 	    double dterm7ydx = dprefacijdx * f * (particles[i].vy-(itemj.vy-center.vy))
 		+ prefacij * dfdx * (particles[i].vy-(itemj.vy-center.vy));
 	    double dterm7ydy = dprefacijdy * f * (particles[i].vy-(itemj.vy-center.vy))
 		+ prefacij * dfdy * (particles[i].vy-(itemj.vy-center.vy));
 	    double dterm7ydz = dprefacijdz * f * (particles[i].vy-(itemj.vy-center.vy))
 		+ prefacij * dfdz * (particles[i].vy-(itemj.vy-center.vy));
 	    double dterm7ydvx = prefacij * dfdvx * (particles[i].vy-(itemj.vy-center.vy));
 
 	    double dterm7ydvy = prefacij * dfdvy * (particles[i].vy-(itemj.vy-center.vy))
 		+ prefacij * f;		
 	    double dterm7ydvz = prefacij * dfdvz * (particles[i].vy-(itemj.vy-center.vy));	    
 
 	    double dterm7zdx = dprefacijdx * f * (particles[i].vz-(itemj.vz-center.vz))
 		+ prefacij * dfdx * (particles[i].vz-(itemj.vz-center.vz));
 	    double dterm7zdy = dprefacijdy * f * (particles[i].vz-(itemj.vz-center.vz))
 		+ prefacij * dfdy * (particles[i].vz-(itemj.vz-center.vz));
 	    double dterm7zdz = dprefacijdz * f * (particles[i].vz-(itemj.vz-center.vz))
 		+ prefacij * dfdz * (particles[i].vz-(itemj.vz-center.vz));
 
 	    double dterm7zdvx = prefacij * dfdvx * (particles[i].vz-(itemj.vz-center.vz));
 
 	    double dterm7zdvy = prefacij * dfdvy * (particles[i].vz-(itemj.vz-center.vz));
 
 	    double dterm7zdvz = prefacij * dfdvz * (particles[i].vz-(itemj.vz-center.vz))
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
 
 	    double axj = 0.0;
 	    double ayj = 0.0;
 	    double azj = 0.0;	    
 	    
 	    for (int k=0; k<N_ephem; k++){
 
 		// Get position and mass of massive body k.
        struct assist_cache_item itemk;
 	    assist_all_ephem(ephem, assist->ephem_cache, k, t, &itemk);
 
 		// Compute position vector of test particle i relative to massive body k.
 		const double dxik = particles[i].x + (center.x - itemk.x); 
 		const double dyik = particles[i].y + (center.y - itemk.y);
 		const double dzik = particles[i].z + (center.z - itemk.z);
 		const double rik2 = dxik*dxik + dyik*dyik + dzik*dzik;
 		const double _rik  = sqrt(rik2);
 
 		// keep track of GM/rik sum
 		term0 += itemk.m/_rik;
 
 		dterm0dx -= itemk.m/(_rik*_rik*_rik) * dxik;
 		dterm0dy -= itemk.m/(_rik*_rik*_rik) * dyik;
 		dterm0dz -= itemk.m/(_rik*_rik*_rik) * dzik;				
 
 		if(k != j){
 		    // Compute position vector of massive body j relative to massive body k.
 		    const double dxjk = itemj.x - itemk.x;
 		    const double dyjk = itemj.y - itemk.y;
 		    const double dzjk = itemj.z - itemk.z;
 		    const double rjk2 = dxjk*dxjk + dyjk*dyjk + dzjk*dzjk;
 		    const double _rjk  = sqrt(rjk2);
 
 		    // keep track of GM/rjk sum
 		    term1 += itemk.m/_rjk;
 
 		    axj -= itemk.m*dxjk/(_rjk*_rjk*_rjk);
 		    ayj -= itemk.m*dyjk/(_rjk*_rjk*_rjk);
 		    azj -= itemk.m*dzjk/(_rjk*_rjk*_rjk);		    		    
 
 		}
 
 	    }
 
 	    term0 *= -2*(beta+gamma)*over_C2;
 	    dterm0dx *= -2*(beta+gamma)*over_C2;
 	    dterm0dy *= -2*(beta+gamma)*over_C2;
 	    dterm0dz *= -2*(beta+gamma)*over_C2;	    	    
 	    
 	    term1 *= -(2*beta-1)*over_C2;
 
 	    const double rijdotaj = dxij*(axj-center.ax) + dyij*(ayj-center.ay) + dzij*(azj-center.az);
 	    const double term6 = -0.5*over_C2*rijdotaj;
 	    const double dterm6dx = -0.5*over_C2*(axj-center.ax);
 	    const double dterm6dy = -0.5*over_C2*(ayj-center.ay);	    
 	    const double dterm6dz = -0.5*over_C2*(azj-center.az);
 	    
 	    double term8x = itemj.m*axj/_rij*(3+4*gamma)/2;
 	    double dterm8xdx = -itemj.m*axj/(_rij*_rij*_rij)*dxij*(3+4*gamma)/2;
 	    double dterm8xdy = -itemj.m*axj/(_rij*_rij*_rij)*dyij*(3+4*gamma)/2;
 	    double dterm8xdz = -itemj.m*axj/(_rij*_rij*_rij)*dzij*(3+4*gamma)/2;	    	    
 
 	    double term8y = itemj.m*ayj/_rij*(3+4*gamma)/2;
 	    double dterm8ydx = -itemj.m*ayj/(_rij*_rij*_rij)*dxij*(3+4*gamma)/2;
 	    double dterm8ydy = -itemj.m*ayj/(_rij*_rij*_rij)*dyij*(3+4*gamma)/2;
 	    double dterm8ydz = -itemj.m*ayj/(_rij*_rij*_rij)*dzij*(3+4*gamma)/2;	    	    
 
 	    double term8z = itemj.m*azj/_rij*(3+4*gamma)/2;
 	    double dterm8zdx = -itemj.m*azj/(_rij*_rij*_rij)*dxij*(3+4*gamma)/2;
 	    double dterm8zdy = -itemj.m*azj/(_rij*_rij*_rij)*dyij*(3+4*gamma)/2;
 	    double dterm8zdz = -itemj.m*azj/(_rij*_rij*_rij)*dzij*(3+4*gamma)/2;	    	    
 
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
 			prefacij*f*(particles[i].vx-(itemj.vx-center.vx)),
 			prefacij*f*(particles[i].vy-(itemj.vy-center.vy)),
 			prefacij*f*(particles[i].vz-(itemj.vz-center.vz)));	    
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
 
// static void assist_additional_force_eih_GR_orig(struct reb_simulation* sim,
// 	    int eih_loop_limit,
// 	    double xo, double yo, double zo,
// 	    double vxo, double vyo, double vzo,
// 	    double axo, double ayo, double azo,	       	    
// 	    FILE *outfile,
// 	    FILE *eih_file){
// 
//     struct assist_extras* assist = (struct assist_extras*) sim->extras;
//     struct assist_ephem* ephem = assist->ephem;
//     const double jd_ref = ephem->jd_ref;
// 
//     // Einstein-Infeld-Hoffman PPN GR treatment
//     // This is one of three options for GR.
//     // This version is only rarely needed.
// 
//     const double au = JPL_EPHEM_CAU;    
//     const double c = (JPL_EPHEM_CLIGHT/au)*86400;
//     const double C2 = c*c;  // This could be stored as C2.
//     const double over_C2 = 1./(c*c);    
// 
//     const unsigned int N = sim->N;  // N includes real+variational particles
//     const unsigned int N_real = N - sim->N_var;
// 
//     const double t = sim->t;    
// 
//     struct reb_particle* const particles = sim->particles;
// 
//     double GM;
//     
//     // The Sun center is reference for these calculations.
//     double xs, ys, zs, vxs, vys, vzs, axs, ays, azs;
//     assist_all_ephem(ephem, assist->ephem_cache, 0, t, &GM, &xs, &ys, &zs, &vxs, &vys, &vzs, &axs, &ays, &azs);
// 
//     double beta = 1.0;
//     double gamma = 1.0;
// 
//     // Loop over test particles        
//     for (int i=0; i<N_real; i++){
// 
// 	double GMj, xj, yj, zj, vxj, vyj, vzj, axj, ayj, azj;    
// 	double GMk, xk, yk, zk, vxk, vyk, vzk, axk, ayk, azk;
// 
// 	// Declare and initialize variational terms
// 	// Only do this if the variational terms are needed.
// 	double dxdx = 0.0;
// 	double dxdy = 0.0;
// 	double dxdz = 0.0;    
// 	double dxdvx = 0.0;
// 	double dxdvy = 0.0;
// 	double dxdvz = 0.0;    
// 	double dydx = 0.0;
// 	double dydy = 0.0;
// 	double dydz = 0.0;    
// 	double dydvx = 0.0;
// 	double dydvy = 0.0;
// 	double dydvz = 0.0;    
// 	double dzdx = 0.0;
// 	double dzdy = 0.0;
// 	double dzdz = 0.0;    
// 	double dzdvx = 0.0;
// 	double dzdvy = 0.0;
// 	double dzdvz = 0.0;    
// 
// 	double term7x_sum = 0.0;
// 	double dterm7x_sumdx = 0.0;
// 	double dterm7x_sumdy = 0.0;
// 	double dterm7x_sumdz = 0.0;		
// 	double dterm7x_sumdvx = 0.0;
// 	double dterm7x_sumdvy = 0.0;
// 	double dterm7x_sumdvz = 0.0;		
// 
// 	double term7y_sum = 0.0;
// 	double dterm7y_sumdx = 0.0;
// 	double dterm7y_sumdy = 0.0;
// 	double dterm7y_sumdz = 0.0;		
// 	double dterm7y_sumdvx = 0.0;
// 	double dterm7y_sumdvy = 0.0;
// 	double dterm7y_sumdvz = 0.0;		
// 	
// 	double term7z_sum = 0.0;
// 	double dterm7z_sumdx = 0.0;
// 	double dterm7z_sumdy = 0.0;
// 	double dterm7z_sumdz = 0.0;		
// 	double dterm7z_sumdvx = 0.0;
// 	double dterm7z_sumdvy = 0.0;
// 	double dterm7z_sumdvz = 0.0;		
// 	
// 	double term8x_sum = 0.0;
// 	double dterm8x_sumdx = 0.0;
// 	double dterm8x_sumdy = 0.0;
// 	double dterm8x_sumdz = 0.0;		
// 
// 	double term8y_sum = 0.0;
// 	double dterm8y_sumdx = 0.0;
// 	double dterm8y_sumdy = 0.0;
// 	double dterm8y_sumdz = 0.0;		
// 	
// 	double term8z_sum = 0.0;
// 	double dterm8z_sumdx = 0.0;
// 	double dterm8z_sumdy = 0.0;
// 	double dterm8z_sumdz = 0.0;		
// 
// 	double grx = 0.0;
// 	double gry = 0.0;
// 	double grz = 0.0;		
// 
// 	for (int j=0; j<eih_loop_limit; j++){ // This is either 1 or N_ephem
// 
// 	    // Get position and mass of massive body j.
// 	    assist_all_ephem(ephem, assist->ephem_cache, j, t, &GMj,
// 		      &xj, &yj, &zj,
// 		      &vxj, &vyj, &vzj,
// 		      &axj, &ayj, &azj);
// 
// 	    // Compute position vector of test particle i relative to massive body j.
// 	    const double dxij = particles[i].x + (xo - xj); 
// 	    const double dyij = particles[i].y + (yo - yj);
// 	    const double dzij = particles[i].z + (zo - zj);
// 	    const double rij2 = dxij*dxij + dyij*dyij + dzij*dzij;
// 	    const double _rij  = sqrt(rij2);
// 	    const double prefacij = GMj/(_rij*_rij*_rij);
// 
// 	    const double dprefacijdx = -3.0*GMj/(_rij*_rij*_rij*_rij*_rij)*dxij;
// 	    const double dprefacijdy = -3.0*GMj/(_rij*_rij*_rij*_rij*_rij)*dyij;
// 	    const double dprefacijdz = -3.0*GMj/(_rij*_rij*_rij*_rij*_rij)*dzij;
// 
// 	    // This is the place to do all the various i-j dot products
// 	    
// 	    const double vi2 = particles[i].vx*particles[i].vx +
// 		particles[i].vy*particles[i].vy +
// 		particles[i].vz*particles[i].vz;
// 
// 	    const double term2 = gamma*over_C2*vi2;
// 	    const double dterm2dvx = 2.0*gamma*over_C2*particles[i].vx;
// 	    const double dterm2dvy = 2.0*gamma*over_C2*particles[i].vy;
// 	    const double dterm2dvz = 2.0*gamma*over_C2*particles[i].vz;	    
// 
// 	    const double vj2 = (vxj-vxo)*(vxj-vxo) + (vyj-vyo)*(vyj-vyo) + (vzj-vzo)*(vzj-vzo);
// 
// 	    const double term3 = (1+gamma)*over_C2*vj2;
// 	    // Variational equations do not depend on term3
// 
// 	    const double vidotvj = particles[i].vx*(vxj-vxo) +
// 		particles[i].vy*(vyj-vyo) +
// 		particles[i].vz*(vzj-vzo);
// 
// 	    const double term4 = -2*(1+gamma)*over_C2*vidotvj;
// 	    const double dterm4dvx = -2*(1+gamma)*over_C2*(vxj-vxo);
// 	    const double dterm4dvy = -2*(1+gamma)*over_C2*(vyj-vyo);
// 	    const double dterm4dvz = -2*(1+gamma)*over_C2*(vzj-vzo);	    	    
// 	    
// 
// 	    const double rijdotvj = dxij*(vxj-vxo) + dyij*(vyj-vyo) + dzij*(vzj-vzo);
// 
// 	    if(eih_file){
// 		fprintf(eih_file, " EIH_J%12d\n", j);	    
// 		fprintf(eih_file, "%25.16lE ", rijdotvj/_rij);
// 	    }
// 
// 	    const double term5 = -1.5*over_C2*(rijdotvj*rijdotvj)/(_rij*_rij);
// 	    const double dterm5dx = -3.0*over_C2*rijdotvj/_rij*((vxj-vxo)/_rij - rijdotvj*dxij/(_rij*_rij*_rij));
// 	    const double dterm5dy = -3.0*over_C2*rijdotvj/_rij*((vyj-vyo)/_rij - rijdotvj*dyij/(_rij*_rij*_rij));
// 	    const double dterm5dz = -3.0*over_C2*rijdotvj/_rij*((vzj-vzo)/_rij - rijdotvj*dzij/(_rij*_rij*_rij));	    	    
// 
// 	    double fx = (2+2*gamma)*particles[i].vx - (1+2*gamma)*(vxj-vxo);
// 	    double fy = (2+2*gamma)*particles[i].vy - (1+2*gamma)*(vyj-vyo);
// 	    double fz = (2+2*gamma)*particles[i].vz - (1+2*gamma)*(vzj-vzo);
// 	    double f = dxij*fx + dyij*fy + dzij*fz;
// 
// 	    double dfdx = fx;
// 	    double dfdy = fy;
// 	    double dfdz = fz;	    	    
// 	    double dfdvx = dxij*(2+2*gamma);
// 	    double dfdvy = dyij*(2+2*gamma);
// 	    double dfdvz = dzij*(2+2*gamma);
// 
// 	    double term7x = prefacij*f*(particles[i].vx-(vxj-vxo));
// 	    double term7y = prefacij*f*(particles[i].vy-(vyj-vyo));
// 	    double term7z = prefacij*f*(particles[i].vz-(vzj-vzo));
// 
// 	    double dterm7xdx = dprefacijdx * f * (particles[i].vx-(vxj-vxo))
// 		+ prefacij * dfdx * (particles[i].vx-(vxj-vxo));
// 	    double dterm7xdy = dprefacijdy * f * (particles[i].vx-(vxj-vxo))
// 		+ prefacij * dfdy * (particles[i].vx-(vxj-vxo));
// 	    double dterm7xdz = dprefacijdz * f * (particles[i].vx-(vxj-vxo))
// 		+ prefacij * dfdz * (particles[i].vx-(vxj-vxo));
// 	    double dterm7xdvx = prefacij * dfdvx * (particles[i].vx-(vxj-vxo))
// 		+ prefacij * f;
// 	    double dterm7xdvy = prefacij * dfdvy * (particles[i].vx-(vxj-vxo));
// 
// 	    double dterm7xdvz = prefacij * dfdvz * (particles[i].vx-(vxj-vxo));	    
// 
// 	    double dterm7ydx = dprefacijdx * f * (particles[i].vy-(vyj-vyo))
// 		+ prefacij * dfdx * (particles[i].vy-(vyj-vyo));
// 	    double dterm7ydy = dprefacijdy * f * (particles[i].vy-(vyj-vyo))
// 		+ prefacij * dfdy * (particles[i].vy-(vyj-vyo));
// 	    double dterm7ydz = dprefacijdz * f * (particles[i].vy-(vyj-vyo))
// 		+ prefacij * dfdz * (particles[i].vy-(vyj-vyo));
// 	    double dterm7ydvx = prefacij * dfdvx * (particles[i].vy-(vyj-vyo));
// 
// 	    double dterm7ydvy = prefacij * dfdvy * (particles[i].vy-(vyj-vyo))
// 		+ prefacij * f;		
// 	    double dterm7ydvz = prefacij * dfdvz * (particles[i].vy-(vyj-vyo));	    
// 
// 	    double dterm7zdx = dprefacijdx * f * (particles[i].vz-(vzj-vzo))
// 		+ prefacij * dfdx * (particles[i].vz-(vzj-vzo));
// 	    double dterm7zdy = dprefacijdy * f * (particles[i].vz-(vzj-vzo))
// 		+ prefacij * dfdy * (particles[i].vz-(vzj-vzo));
// 	    double dterm7zdz = dprefacijdz * f * (particles[i].vz-(vzj-vzo))
// 		+ prefacij * dfdz * (particles[i].vz-(vzj-vzo));
// 
// 	    double dterm7zdvx = prefacij * dfdvx * (particles[i].vz-(vzj-vzo));
// 
// 	    double dterm7zdvy = prefacij * dfdvy * (particles[i].vz-(vzj-vzo));
// 
// 	    double dterm7zdvz = prefacij * dfdvz * (particles[i].vz-(vzj-vzo))
// 		+ prefacij * f;
// 	    
// 	    term7x_sum += term7x;
// 	    term7y_sum += term7y;
// 	    term7z_sum += term7z;
// 
// 	    dterm7x_sumdx += dterm7xdx;
// 	    dterm7x_sumdy += dterm7xdy;
// 	    dterm7x_sumdz += dterm7xdz;	    
// 	    dterm7x_sumdvx += dterm7xdvx;
// 	    dterm7x_sumdvy += dterm7xdvy;
// 	    dterm7x_sumdvz += dterm7xdvz;	    
// 
// 	    dterm7y_sumdx += dterm7ydx;
// 	    dterm7y_sumdy += dterm7ydy;
// 	    dterm7y_sumdz += dterm7ydz;
// 	    dterm7y_sumdvx += dterm7ydvx;
// 	    dterm7y_sumdvy += dterm7ydvy;
// 	    dterm7y_sumdvz += dterm7ydvz;
// 
// 	    dterm7z_sumdx += dterm7zdx;
// 	    dterm7z_sumdy += dterm7zdy;
// 	    dterm7z_sumdz += dterm7zdz;	    
// 	    dterm7z_sumdvx += dterm7zdvx;
// 	    dterm7z_sumdvy += dterm7zdvy;
// 	    dterm7z_sumdvz += dterm7zdvz;	    
// 	    
// 
// 	    double term0 = 0.0;
// 	    double dterm0dx = 0.0;
// 	    double dterm0dy = 0.0;
// 	    double dterm0dz = 0.0;	    
// 
// 	    double term1 = 0.0;
// 	    double dterm1dx = 0.0;
// 	    double dterm1dy = 0.0;
// 	    double dterm1dz = 0.0;	    
// 	    double dterm1dvx = 0.0;
// 	    double dterm1dvy = 0.0;
// 	    double dterm1dvz = 0.0;	    
// 
// 	    axj = 0.0;
// 	    ayj = 0.0;
// 	    azj = 0.0;	    
// 	    
// 	    for (int k=0; k<N_ephem; k++){
// 
// 		// Get position and mass of massive body k.
// 		assist_all_ephem(ephem, assist->ephem_cache, k, t, &GMk,
// 			  &xk, &yk, &zk,
// 			  &vxk, &vyk, &vzk,
// 			  &axk, &ayk, &azk);
// 
// 		// Compute position vector of test particle i relative to massive body k.
// 		const double dxik = particles[i].x + (xo - xk); 
// 		const double dyik = particles[i].y + (yo - yk);
// 		const double dzik = particles[i].z + (zo - zk);
// 		const double rik2 = dxik*dxik + dyik*dyik + dzik*dzik;
// 		const double _rik  = sqrt(rik2);
// 
// 		// keep track of GM/rik sum
// 		term0 += GMk/_rik;
// 
// 		dterm0dx -= GMk/(_rik*_rik*_rik) * dxik;
// 		dterm0dy -= GMk/(_rik*_rik*_rik) * dyik;
// 		dterm0dz -= GMk/(_rik*_rik*_rik) * dzik;				
// 
// 		if(k != j){
// 		    // Compute position vector of massive body j relative to massive body k.
// 		    const double dxjk = xj - xk;
// 		    const double dyjk = yj - yk;
// 		    const double dzjk = zj - zk;
// 		    const double rjk2 = dxjk*dxjk + dyjk*dyjk + dzjk*dzjk;
// 		    const double _rjk  = sqrt(rjk2);
// 
// 		    // keep track of GM/rjk sum
// 		    term1 += GMk/_rjk;
// 
// 		    axj -= GMk*dxjk/(_rjk*_rjk*_rjk);
// 		    ayj -= GMk*dyjk/(_rjk*_rjk*_rjk);
// 		    azj -= GMk*dzjk/(_rjk*_rjk*_rjk);		    		    
// 
// 		}
// 
// 	    }
// 
// 	    term0 *= -2*(beta+gamma)*over_C2;
// 	    dterm0dx *= -2*(beta+gamma)*over_C2;
// 	    dterm0dy *= -2*(beta+gamma)*over_C2;
// 	    dterm0dz *= -2*(beta+gamma)*over_C2;	    	    
// 	    
// 	    term1 *= -(2*beta-1)*over_C2;
// 
// 	    const double rijdotaj = dxij*(axj-axo) + dyij*(ayj-ayo) + dzij*(azj-azo);
// 	    const double term6 = -0.5*over_C2*rijdotaj;
// 	    const double dterm6dx = -0.5*over_C2*(axj-axo);
// 	    const double dterm6dy = -0.5*over_C2*(ayj-ayo);	    
// 	    const double dterm6dz = -0.5*over_C2*(azj-azo);
// 	    
// 	    double term8x = GMj*axj/_rij*(3+4*gamma)/2;
// 	    double dterm8xdx = -GMj*axj/(_rij*_rij*_rij)*dxij*(3+4*gamma)/2;
// 	    double dterm8xdy = -GMj*axj/(_rij*_rij*_rij)*dyij*(3+4*gamma)/2;
// 	    double dterm8xdz = -GMj*axj/(_rij*_rij*_rij)*dzij*(3+4*gamma)/2;	    	    
// 
// 	    double term8y = GMj*ayj/_rij*(3+4*gamma)/2;
// 	    double dterm8ydx = -GMj*ayj/(_rij*_rij*_rij)*dxij*(3+4*gamma)/2;
// 	    double dterm8ydy = -GMj*ayj/(_rij*_rij*_rij)*dyij*(3+4*gamma)/2;
// 	    double dterm8ydz = -GMj*ayj/(_rij*_rij*_rij)*dzij*(3+4*gamma)/2;	    	    
// 
// 	    double term8z = GMj*azj/_rij*(3+4*gamma)/2;
// 	    double dterm8zdx = -GMj*azj/(_rij*_rij*_rij)*dxij*(3+4*gamma)/2;
// 	    double dterm8zdy = -GMj*azj/(_rij*_rij*_rij)*dyij*(3+4*gamma)/2;
// 	    double dterm8zdz = -GMj*azj/(_rij*_rij*_rij)*dzij*(3+4*gamma)/2;	    	    
// 
// 	    term8x_sum += term8x;
// 	    term8y_sum += term8y;
// 	    term8z_sum += term8z;
// 
// 	    dterm8x_sumdx += dterm8xdx;
// 	    dterm8x_sumdy += dterm8xdy;
// 	    dterm8x_sumdz += dterm8xdz;	    
// 
// 	    dterm8y_sumdx += dterm8ydx;
// 	    dterm8y_sumdy += dterm8ydy;
// 	    dterm8y_sumdz += dterm8ydz;
// 
// 	    dterm8z_sumdx += dterm8zdx;
// 	    dterm8z_sumdy += dterm8zdy;
// 	    dterm8z_sumdz += dterm8zdz;	    
// 	    
// 	    double factor = term0 + term1 + term2 + term3 + term4 + term5 + term6;
// 
// 	    double dfactordx = dterm0dx + dterm1dx + dterm5dx + dterm6dx;
// 	    double dfactordy = dterm0dy + dterm1dy + dterm5dy + dterm6dy;
// 	    double dfactordz = dterm0dz + dterm1dz + dterm5dz + dterm6dz;	    
// 	    double dfactordvx = dterm1dvx + dterm2dvx + dterm4dvx;
// 	    double dfactordvy = dterm1dvy + dterm2dvy + dterm4dvy;
// 	    double dfactordvz = dterm1dvz + dterm2dvz + dterm4dvz;
// 
// 	    if(eih_file){
// 		fprintf(eih_file, "%24.16lE ", -factor*C2);
// 		fprintf(eih_file, "%24.16lE %24.16lE %24.16lE %24.16lE ",
// 			-factor*C2*prefacij*dxij,
// 			-factor*C2*prefacij*dyij,
// 			-factor*C2*prefacij*dzij,
// 			f);	    
// 		fprintf(eih_file, "%24.16lE %24.16lE %24.16lE ",
// 			prefacij*f*(particles[i].vx-(vxj-vxo)),
// 			prefacij*f*(particles[i].vy-(vyj-vyo)),
// 			prefacij*f*(particles[i].vz-(vzj-vzo)));	    
// 		fprintf(eih_file, "%24.16lE %24.16lE %24.16lE ",
// 			term8x,
// 			term8y,
// 			term8z);
// 		fprintf(eih_file, "%24.16lE %24.16lE %24.16lE\n",
// 			axj, ayj, azj);
// 
// 		fflush(eih_file);
// 	    }
// 
// 	    grx += -prefacij*dxij*factor;
// 	    gry += -prefacij*dyij*factor;
// 	    grz += -prefacij*dzij*factor;
// 	    
// 	    particles[i].ax += -prefacij*dxij*factor;
// 	    particles[i].ay += -prefacij*dyij*factor;
// 	    particles[i].az += -prefacij*dzij*factor;
// 
// 	    // Variational equation terms go here.
// 
// 	    dxdx += -dprefacijdx*dxij*factor
// 		-prefacij*factor
// 		-prefacij*dxij*dfactordx;
// 	    
// 	    dxdy += -dprefacijdy*dxij*factor
// 		-prefacij*dxij*dfactordy;
// 	    
// 	    dxdz += -dprefacijdz*dxij*factor
// 		-prefacij*dxij*dfactordz;
// 
// 	    dxdvx += 
// 		-prefacij*dxij*dfactordvx;
// 
// 	    dxdvy += 
// 		-prefacij*dxij*dfactordvy;
// 
// 	    dxdvz += 
// 		-prefacij*dxij*dfactordvz;
// 
// 	    dydx += -dprefacijdx*dyij*factor
// 		-prefacij*dyij*dfactordx;
// 	    
// 	    dydy += -dprefacijdy*dyij*factor
// 		-prefacij*factor
// 		-prefacij*dyij*dfactordy;
// 	    
// 	    dydz += -dprefacijdz*dyij*factor
// 		-prefacij*dyij*dfactordz;
// 
// 	    dydvx += 
// 		-prefacij*dyij*dfactordvx;
// 
// 	    dydvy += 
// 		-prefacij*dyij*dfactordvy;
// 
// 	    dydvz += 
// 		-prefacij*dyij*dfactordvz;
// 	    
// 	    dzdx += -dprefacijdx*dzij*factor
// 		-prefacij*dzij*dfactordx;
// 
// 	    dzdy += -dprefacijdy*dzij*factor
// 		-prefacij*dzij*dfactordy;
// 	    
// 	    dzdz += -dprefacijdz*dzij*factor
// 		-prefacij*factor
// 		-prefacij*dzij*dfactordz;
// 
// 	    dzdvx += 
// 		-prefacij*dzij*dfactordvx;
// 
// 	    dzdvy += 
// 		-prefacij*dzij*dfactordvy;
// 
// 	    dzdvz += 
// 		-prefacij*dzij*dfactordvz;
// 
//         }
// 
// 	grx += term7x_sum*over_C2 + term8x_sum*over_C2;
// 	gry += term7y_sum*over_C2 + term8y_sum*over_C2;
// 	grz += term7z_sum*over_C2 + term8z_sum*over_C2;
// 
// 	if(outfile){
// 	    fprintf(outfile, "%3s %25.16le %25.16le %25.16le %25.16le\n", "GR", jd_ref+t,
// 		    grx, gry, grz);
// 	    fflush(outfile);
// 	}
// 
// 	dxdx += dterm7x_sumdx*over_C2 + dterm8x_sumdx*over_C2;
// 	dxdy += dterm7x_sumdy*over_C2 + dterm8x_sumdy*over_C2;
// 	dxdz += dterm7x_sumdz*over_C2 + dterm8x_sumdz*over_C2;	
// 	dxdvx += dterm7x_sumdvx*over_C2;
// 	dxdvy += dterm7x_sumdvy*over_C2;
// 	dxdvz += dterm7x_sumdvz*over_C2;
// 
// 	dydx += dterm7y_sumdx*over_C2 + dterm8y_sumdx*over_C2;
// 	dydy += dterm7y_sumdy*over_C2 + dterm8y_sumdy*over_C2;
// 	dydz += dterm7y_sumdz*over_C2 + dterm8y_sumdz*over_C2;	
// 	dydvx += dterm7y_sumdvx*over_C2;
// 	dydvy += dterm7y_sumdvy*over_C2;
// 	dydvz += dterm7y_sumdvz*over_C2;
// 
// 	dzdx += dterm7z_sumdx*over_C2 + dterm8z_sumdx*over_C2;
// 	dzdy += dterm7z_sumdy*over_C2 + dterm8z_sumdy*over_C2;
// 	dzdz += dterm7z_sumdz*over_C2 + dterm8z_sumdz*over_C2;	
// 	dzdvx += dterm7z_sumdvx*over_C2;
// 	dzdvy += dterm7z_sumdvy*over_C2;
// 	dzdvz += dterm7z_sumdvz*over_C2;
// 	
// 	particles[i].ax += term7x_sum*over_C2 + term8x_sum*over_C2;
// 	particles[i].ay += term7y_sum*over_C2 + term8y_sum*over_C2;
// 	particles[i].az += term7z_sum*over_C2 + term8z_sum*over_C2;
// 
// 	// Variational equation terms go here.
// 	for (int v=0; v < sim->var_config_N; v++){
// 	    struct reb_variational_configuration const vc = sim->var_config[v];
// 	    int tp = vc.testparticle;
// 	    struct reb_particle* const particles_var1 = particles + vc.index;		
// 	    if(tp == i){
// 	    
// 		// variational particle coords
// 		const double ddx = particles_var1[0].x;
// 		const double ddy = particles_var1[0].y;
// 		const double ddz = particles_var1[0].z;
// 		const double ddvx = particles_var1[0].vx;
// 		const double ddvy = particles_var1[0].vy;
// 		const double ddvz = particles_var1[0].vz;
// 
// 		// Matrix multiplication
// 		const double dax =   ddx  * dxdx  + ddy  * dxdy  + ddz  * dxdz
// 		    +   ddvx * dxdvx + ddvy * dxdvy + ddvz * dxdvz;
// 		const double day =   ddx  * dydx  + ddy  * dydy  + ddz  * dydz
// 		    +   ddvx * dydvx + ddvy * dydvy + ddvz * dydvz;
// 		const double daz =   ddx  * dzdx  + ddy  * dzdy  + ddz  * dzdz
// 		    +   ddvx * dzdvx + ddvy * dzdvy + ddvz * dzdvz;
// 
// 		// Accumulate acceleration terms
// 		particles_var1[0].ax += dax;
// 		particles_var1[0].ay += day;
// 		particles_var1[0].az += daz;
// 		
// 	    }
// 	}
//     }
// }
