/**
 * @file    assist.c
 * @brief   Central internal functions for ASSIST
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

#define FNAMESIZE 256
#define DEFAULT_JPL_SB_EPHEM "sb441-n16.bsp"

#include <stdint.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <limits.h>
#include <float.h>
#include "assist.h"
#include "rebound.h"

#include "spk.h"
#include "planets.h"

/**
 * @brief Struct containing pointers to intermediate values
 */
struct reb_dpconst7 {
    double* const restrict p0;  ///< Temporary values at intermediate step 0 
    double* const restrict p1;  ///< Temporary values at intermediate step 1 
    double* const restrict p2;  ///< Temporary values at intermediate step 2 
    double* const restrict p3;  ///< Temporary values at intermediate step 3 
    double* const restrict p4;  ///< Temporary values at intermediate step 4 
    double* const restrict p5;  ///< Temporary values at intermediate step 5 
    double* const restrict p6;  ///< Temporary values at intermediate step 6 
};

static struct reb_dpconst7 dpcast(struct reb_dp7 dp){
    struct reb_dpconst7 dpc = {
        .p0 = dp.p0, 
        .p1 = dp.p1, 
        .p2 = dp.p2, 
        .p3 = dp.p3, 
        .p4 = dp.p4, 
        .p5 = dp.p5, 
        .p6 = dp.p6, 
    };
    return dpc;
}

const int reb_max_messages_length = 1024;   // needs to be constant expression for array size
const int reb_max_messages_N = 10;

enum {
    NO_ERR,        // no error
    ERR_JPL_EPHEM, // JPL ephemeris file not found
    ERR_JPL_AST,   // JPL asteroid file not found
    ERR_NAST,      // asteroid number out of range
    ERR_NEPH,      // planet number out of range
};

int ebody[11] = {
        PLAN_SOL,                       // Sun (in barycentric)
        PLAN_MER,                       // Mercury center
        PLAN_VEN,                       // Venus center
        PLAN_EAR,                       // Earth center
        PLAN_LUN,                       // Moon center
        PLAN_MAR,                       // Mars center
        PLAN_JUP,                       // ...
        PLAN_SAT,
        PLAN_URA,
        PLAN_NEP,
        PLAN_PLU
};

#define STRINGIFY(s) str(s)
#define str(s) #s

const char* assist_build_str = __DATE__ " " __TIME__;   // Date and time build string. 
const char* assist_version_str = "1.0.0";         // **VERSIONLINE** This line gets updated automatically. Do not edit manually.
const char* assist_githash_str = STRINGIFY(ASSISTGITHASH);// This line gets updated automatically. Do not edit manually.

// TODO: this should be more general.  Perhaps
// this information could come directly from
// the JPL ephemerides.
static int number_bodies(int* N_ephem, int* N_ast){
    *N_ephem = 11;
    *N_ast = 16;

    return(*N_ephem + *N_ast);
}

// Added gravitational constant G (2020 Feb 26)
// Added vx, vy, vz for GR stuff (2020 Feb 27)
// Consolidated the routine, removing the if block.
//
static int ephem(const int i, const double jde, double* const GM,
	   double* const x, double* const y, double* const z,
	   double* const vx, double* const vy, double* const vz,
	   double* const ax, double* const ay, double* const az){

    static int initialized = 0;

    static struct _jpl_s *pl;
    struct mpos_s now;

    // The values below are G*mass.
    // Units are solar masses, au, days.
    // TODO: The units should probably be handled elsewhere.
    // DE440/441 units: au^3 day^-2.
    // TODO: These should be moved to an external source that
    // is easily modified.
    const static double JPL_GM[11] =
	{
	    0.2959122082841196e-03, // 0 sun
	    0.4912500194889318e-10, // 1 mercury
	    0.7243452332644119e-09, // 2 venus
	    0.8887692446707102e-09, // 3 earth
	    0.1093189462402435e-10, // 4 moon	    
	    0.9549548829725812e-10, // 5 mars
	    0.2825345825225792e-06, // 6 jupiter
	    0.8459705993376290e-07, // 7 saturn
	    0.1292026564968240e-07, // 8 uranus
	    0.1524357347885194e-07, // 9 neptune
	    0.2175096464893358e-11, // 10 pluto
	};

    if(i<0 || i>10){
	return(ERR_NEPH);
    }

    if (initialized == 0){
      if ((pl = jpl_init()) == NULL) {
	  return(ERR_JPL_EPHEM);	  
      }
      initialized = 1;
    }

    // Get position, velocity, and mass of body i in barycentric coords. 
    
    *GM = JPL_GM[i];

    jpl_calc(pl, &now, jde, ebody[i], PLAN_BAR); 

    // Convert to au/day and au/day^2
    // TODO: Consider making the units more flexible.
    vecpos_div(now.u, pl->cau);
    vecpos_div(now.v, pl->cau/86400.);
    vecpos_div(now.w, pl->cau/(86400.*86400.));

    *x = now.u[0];
    *y = now.u[1];
    *z = now.u[2];
    *vx = now.v[0];
    *vy = now.v[1];
    *vz = now.v[2];
    *ax = now.w[0];
    *ay = now.w[1];
    *az = now.w[2];

    return(NO_ERR);
    
}

static int ast_ephem(const int i, const double jde, double* const GM, double* const x, double* const y, double* const z){

    static int initialized = 0;

    static struct spk_s *spl;
    struct mpos_s pos;

    // The values below are G*mass.
    // Units are solar masses, au, days.
    // DE441
    // GMs were supplied by Davide.
    // TODO: these should be moved to an external
    // source that can be easily updated.
    const static double JPL_GM[16] =    
    {
	    3.2191392075878588e-15, // 107 camilla
	    1.3964518123081070e-13, // 1 ceres
	    2.0917175955133682e-15, // 65 cybele	
	    8.6836253492286545e-15, // 511 davida
	    4.5107799051436795e-15, // 15 eunomia
	    2.4067012218937576e-15, // 31 euphrosyne	    
	    5.9824315264869841e-15, // 52 europa
	    1.2542530761640810e-14, // 10 hygiea
	    6.3110343420878887e-15, // 704 interamnia
	    2.5416014973471498e-15, // 7 iris	    
	    4.2823439677995011e-15, // 3 juno
	    3.0471146330043200e-14, // 2 pallas
	    3.5445002842488978e-15, // 16 psyche
	    4.8345606546105521e-15, // 87 sylvia
	    2.6529436610356353e-15, // 88 thisbe
	    3.8548000225257904e-14, // 4 vesta          

    };

    if(i<0 || i>15){
	return(ERR_NAST);
    }

    if (initialized == 0){

	char buf[FNAMESIZE];

        /** use or environment-specified file, 
	 * or the default filename, in that order
         */
        if (getenv("JPL_SB_EPHEM")!=NULL)
	    strncpy(buf, getenv("JPL_SB_EPHEM"), FNAMESIZE-1);
        else
	    strncpy(buf, DEFAULT_JPL_SB_EPHEM, FNAMESIZE-1);

	//printf("%s\n", buf);
	FILE *file;
	file = fopen(buf, "r");
	if(file == NULL){
	    printf("couldn't open asteroid file.\n");
	}
	      

	if ((spl = spk_init(buf)) == NULL) {	
	    //if ((spl = spk_init("sb441-n16.bsp")) == NULL) {
	    return(ERR_JPL_AST);
      }
      
      initialized = 1;

    }

    // TODO: again, the units might be handled more
    // generally

    *GM = JPL_GM[i];
    spk_calc(spl, i, jde, &pos);          
    *x = pos.u[0];
    *y = pos.u[1];
    *z = pos.u[2];

    return(NO_ERR);

}

int all_ephem(const int i, const double t, double* const GM,
		      double* const x, double* const y, double* const z,
		      double* const vx, double* const vy, double* const vz,
		      double* const ax, double* const ay, double* const az
		      ){    

    int number_bodies(int* N_ephem, int* N_ast);
    static int N_ast = -1;
    static int N_ephem = -1;

    if(N_ast == -1 || N_ephem == -1){
	number_bodies(&N_ephem, &N_ast);
    }

    static double xs, ys, zs, vxs, vys, vzs, axs, ays, azs;
    static double GMs;
    static double t_last = -1e99;

    // For any given step, using the IAS15 integrator,
    // all_ephem will need to access positions and GM values
    // for 27 bodies at the times of 8 different substeps.
    // The integrator loops through through the same substeps
    // to convergence.
    // TODO: We can optimize by saving the positions
    // of the perturbers, rather than reevaluating them at
    // each iteration.
    
    // Get position and mass of massive body i.
    if(i < N_ephem){
	int flag = ephem(i, t, GM, x, y, z, vx, vy, vz, ax, ay, az);
	if(flag != NO_ERR) return(flag);
    }else{
	// Get position and mass of asteroid i-N_ephem.
	int flag = ast_ephem(i-N_ephem, t, GM, x, y, z);
	if(flag != NO_ERR) return(flag);	

	if(t != t_last){
	    flag = ephem(0, t, &GMs, &xs, &ys, &zs, &vxs, &vys, &vzs, &axs, &ays, &azs);
	    if(flag != NO_ERR) return(flag);		    
	    t_last = t;
	}

	// Translate massive asteroids from heliocentric to barycentric.
	*x += xs; *y += ys; *z += zs;
	*vx = NAN; *vy = NAN; *vz = NAN;
	*ax = NAN; *ay = NAN; *az = NAN;		
    }

    return(NO_ERR);
}

void assist_additional_forces(struct reb_simulation* sim){

    // Doesn't need to be hard-coded.
    const double c = 173.14463267424031;
    int* geo = 0;    
    sim->G = 0.295912208285591100E-03; // Gravitational constant (AU, solar masses, days)

    const double C2 = c*c;  // This could be stored as C2.

    struct assist_extras* assist = (struct assist_extras*) sim->extras;

    int nsubsteps = assist->nsubsteps;
    double* hg = assist->hg;

    // implement additional_forces here
    const double G = sim->G;
    const int N = sim->N;  // N includes real+variational particles
    const int N_real= N;

    struct reb_particle* const particles = sim->particles;

    int N_ephem, N_ast;
    int number_bodies(int* N_ephem, int* N_ast);

    const double t = sim->t;

    const int N_tot = number_bodies(&N_ephem, &N_ast);

    double GM;
    double x, y, z, vx, vy, vz, ax, ay, az;

    // TODO: eliminate these output files after testing.
    FILE *outfile;
    outfile = fopen("acc.out", "w");

    // Get mass, position, velocity, and acceleration of the Earth and Sun
    // for later use.
    // The hard-wired constants should be changed.
    
    // The offset position is used to adjust the particle positions.
    // The options are the barycenter (default) and geocenter.
    double xo, yo, zo, vxo, vyo, vzo, axo, ayo, azo;

    if(*geo == 1){
	// geocentric
	int flag = all_ephem(3, t, &GM, &xo, &yo, &zo, &vxo, &vyo, &vzo, &axo, &ayo, &azo);
	if(flag != NO_ERR){
	    char outstring[50];
	    sprintf(outstring, "%s %d %d\n", "Ephemeris error a ", 3, flag);
	    reb_error(sim, outstring);
	}
    }else{
	// barycentric
	xo = 0.0;  yo = 0.0;  zo = 0.0;
	vxo = 0.0; vyo = 0.0; vzo = 0.0;      
    }

    // Direct forces from massives bodies
    for (int i=0; i<N_tot; i++){
	
        // Get position and mass of massive body i.
	// TOOD: make a version that returns the positions, velocities,
	// and accelerations for all the bodies at a given time.
	int flag = all_ephem(i, t, &GM, &x, &y, &z, &vx, &vy, &vz, &ax, &ay, &az);
	if(flag != NO_ERR){
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

	    fprintf(outfile, "%3d %25.16le %25.16le %25.16le %25.16le %25.16le %25.16le %25.16le %25.16le\n", i, t, GM, dx, dy, dz, -prefac*dx, -prefac*dy, -prefac*dz);

/* This code is repeated below in the variational particle loop; Remove. AA 5/24/22
            // Values and cooefficients for variational equations
            // Looks like extra 
            const double r3inv = 1./(r2*_r);
            const double r5inv = 3.*r3inv/r2;

            const double dxdx = dx*dx*r5inv - r3inv;
            const double dydy = dy*dy*r5inv - r3inv;
            const double dzdz = dz*dz*r5inv - r3inv;
            const double dxdy = dx*dy*r5inv;
            const double dxdz = dx*dz*r5inv;
            const double dydz = dy*dz*r5inv;
*/

	    particles[j].ax -= prefac*dx;
	    particles[j].ay -= prefac*dy;
	    particles[j].az -= prefac*dz;

        }
    }

    fclose(outfile);

    
    // Acceleration of variational particles due to direct forces from massive bodies 
    // Loop over the perturbers
    for (int i=0; i<N_tot; i++){

        // Get position and mass of massive body i.	
	all_ephem(i, t, &GM, &x, &y, &z, &vx, &vy, &vz, &ax, &ay, &az);	
	
        for (int j=0; j<N_real; j++){ //loop over test particles

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
		    const double Gmi = GM;

		    // Matrix multiplication
		    const double dax =   ddx * dxdx + ddy * dxdy + ddz * dxdz;
		    const double day =   ddx * dxdy + ddy * dydy + ddz * dydz;
		    const double daz =   ddx * dxdz + ddy * dydz + ddz * dzdz;

		    // No variational mass contributions for test particles!

		    // Accumulate acceleration terms
		    particles_var1[0].ax += Gmi * dax; 
		    particles_var1[0].ay += Gmi * day; 
		    particles_var1[0].az += Gmi * daz; 

		}
	    }
        }
    }

    // We might move this into a somewhat separate part of the code,
    // similar to how different extra forces are typically handled in
    // reboundx
    // 
    // Here is the treatment of the Earth's J2 and J4.
    // Borrowed code from gravitational_harmonics example.
    // Assumes the coordinates are geocentric.
    // Also assuming that Earth's pole is along the z
    // axis.  This is only precisely true at the J2000
    // epoch.
    //

    // The geocenter is the reference for the Earth J2/J4 calculations.
    double xe, ye, ze, vxe, vye, vze, axe, aye, aze;    
    all_ephem(3, t, &GM, &xe, &ye, &ze, &vxe, &vye, &vze, &axe, &aye, &aze);

    double xr, yr, zr, vxr, vyr, vzr, axr, ayr, azr;
    xr = xe;  yr = ye;  zr = ze;

    // Hard-coded constants.  BEWARE!
    // Clean up on aisle 3!
    const double GMearth = 0.888769244512563400E-09;
    const double J2e = 0.00108262545;
    const double J4e = -0.000001616;
    const double au = 149597870.700;
    const double Re_eq = 6378.1263/au;
    // Unit vector to equatorial pole at the epoch
    // Clean this up!
    // Note also that the pole orientation is not changing during
    // the integration.

    double RAs =  359.87123273*M_PI/180.;
    double Decs =  89.88809752*M_PI/180.;

    //double xp = cos(Decs)*cos(RAs);
    //double yp = cos(Decs)*sin(RAs);
    //double zp = sin(Decs);

    double xp =  0.0019111736356920146;
    double yp = -1.2513100974355823e-05;
    double zp =   0.9999981736277104;

    //double xp =  0.0;
    //double yp =  0.0;
    //double zp =  1.0;
    
    double incl = acos(zp);
    double longnode;
    if(xp != 0.0 || yp !=0.0) {    
      longnode = atan2(xp, -yp);
    } else {
      longnode = 0.0;
    }

    // Rearrange this loop for efficiency
    for (int j=0; j<N_real; j++){

        const struct reb_particle p = particles[j];
        double dx = p.x + (xo - xr);
        double dy = p.y + (yo - yr);
        double dz = p.z + (zo - zr);

        const double r2 = dx*dx + dy*dy + dz*dz;
        const double r = sqrt(r2);
	
	// Rotate to Earth equatorial frame
	// This could be a single rotation

	// Rotate around z by RA
	double cosr = cos(-longnode);
	double sinr = sin(-longnode);

	double dxp =  dx * cosr - dy * sinr;
	double dyp =  dx * sinr + dy * cosr;
	double dzp =  dz;

	// Rotate around x by Dec
	double cosd = cos(-incl);
	double sind = sin(-incl);
	
	dx =  dxp;
	dy =  dyp * cosd - dzp * sind;
	dz =  dyp * sind + dzp * cosd;

	// Calculate acceleration in
	// Earth equatorial frame	

	// J2 terms
        const double costheta2 = dz*dz/r2;
        const double J2e_prefac = 3.*J2e*Re_eq*Re_eq/r2/r2/r/2.;
        const double J2e_fac = 5.*costheta2-1.;
        const double J2e_fac2 = 7.*costheta2-1.;
        const double J2e_fac3 = 35.*costheta2*costheta2-30.*costheta2+3.;

	double resx = GMearth*J2e_prefac*J2e_fac*dx;
	double resy = GMearth*J2e_prefac*J2e_fac*dy;
	double resz = GMearth*J2e_prefac*(J2e_fac-2.)*dz;	

	// J4 terms
        const double J4e_prefac = 5.*J4e*Re_eq*Re_eq*Re_eq*Re_eq/r2/r2/r2/r/8.;
        const double J4e_fac = 63.*costheta2*costheta2-42.*costheta2 + 3.;
        const double J4e_fac2= 33.*costheta2*costheta2-18.*costheta2 + 1.;
        const double J4e_fac3= 33.*costheta2*costheta2-30.*costheta2 + 5.;
        const double J4e_fac4= 231.*costheta2*costheta2*costheta2-315.*costheta2*costheta2+105.*costheta2 - 5.;

        resx += GMearth*J4e_prefac*J4e_fac*dx;
        resy += GMearth*J4e_prefac*J4e_fac*dy;
        resz += GMearth*J4e_prefac*(J4e_fac+12.-28.*costheta2)*dz;

	// Rotate back to original frame
	// Rotate around x by -Dec
	double resxp =  resx;
	double resyp =  resy * cosd + resz * sind;
	double reszp = -resy * sind + resz * cosd;
	
	// Rotate around z by -RA
	resx =  resxp * cosr + resyp * sinr;
	resy = -resxp * sinr + resyp * cosr;
	resz =  reszp;

	// Accumulate final acceleration terms
  	particles[j].ax += resx;
        particles[j].ay += resy; 
        particles[j].az += resz;

	// Constants for variational equations
	// J2 terms
	const double dxdx = GMearth*J2e_prefac*(J2e_fac-5.*J2e_fac2*dx*dx/r2);
	const double dydy = GMearth*J2e_prefac*(J2e_fac-5.*J2e_fac2*dy*dy/r2);
	const double dzdz = GMearth*J2e_prefac*(-1.)*J2e_fac3;
	const double dxdy = GMearth*J2e_prefac*(-5.)*J2e_fac2*dx*dy/r2;
	const double dydz = GMearth*J2e_prefac*(-5.)*(J2e_fac2-2.)*dy*dz/r2;
	const double dxdz = GMearth*J2e_prefac*(-5.)*(J2e_fac2-2.)*dx*dz/r2;
	// J4 terms
	const double dxdxJ4 = GMearth*J4e_prefac*(J4e_fac-21.*J4e_fac2*dx*dx/r2);
	const double dydyJ4 = GMearth*J4e_prefac*(J4e_fac-21.*J4e_fac2*dy*dy/r2);
	const double dzdzJ4 = GMearth*J4e_prefac*(-3.)*J4e_fac4;
	const double dxdyJ4 = GMearth*J4e_prefac*(-21.)*J4e_fac2*dx*dy/r2;
	const double dydzJ4 = GMearth*J4e_prefac*(-21.)*J4e_fac3*dy*dz/r2;
	const double dxdzJ4 = GMearth*J4e_prefac*(-21.)*J4e_fac3*dx*dz/r2;

	// Looping over variational particles
        for(int v = N_real + 6*j; v < N_real + 6*(j+1); v++){

	    double ddx = particles[v].x;
	    double ddy = particles[v].y;
	    double ddz = particles[v].z;

	    // Rotate to Earth equatorial frame
	    double ddxp =  ddx * cosr - ddy * sinr;
	    double ddyp =  ddx * sinr + ddy * cosr;
	    double ddzp =  ddz;
	    ddx =  ddxp;
	    ddy =  ddyp * cosd - ddzp * sind;
	    ddz =  ddyp * sind + ddzp * cosd;

	    // Matrix multiplication
	    double dax =   ddx * dxdx + ddy * dxdy + ddz * dxdz;
	    double day =   ddx * dxdy + ddy * dydy + ddz * dydz;
	    double daz =   ddx * dxdz + ddy * dydz + ddz * dzdz;

	    dax +=   ddx * dxdxJ4 + ddy * dxdyJ4 + ddz * dxdzJ4;
	    day +=   ddx * dxdyJ4 + ddy * dydyJ4 + ddz * dydzJ4;
	    daz +=   ddx * dxdzJ4 + ddy * dydzJ4 + ddz * dzdzJ4;

	    // Rotate back
	    double daxp =  dax;
	    double dayp =  day * cosd + daz * sind;
	    double dazp = -day * sind + daz * cosd;
	    dax =  daxp * cosr + dayp * sinr;
	    day = -daxp * sinr + dayp * cosr;
	    daz =  dazp;

	    // Accumulate acceleration terms
	    particles[v].ax += dax;
	    particles[v].ay += day;
	    particles[v].az += daz;

        }
    }

    // We might move this into a somewhat separate part of the code,
    // similar to how different extra forces are typically handled in
    // reboundx
    // Here is the treatment of the Sun's J2.
    // Borrowed code from gravitational_harmonics.

    // The Sun center is reference for these calculations.

    //all_ephem(i, t, &GM, &x, &y, &z, &vx, &vy, &vz, &ax, &ay, &az);	    
    all_ephem(0, t, &GM, &xr, &yr, &zr, &vxr, &vyr, &vzr, &axr, &ayr, &azr);	    

    // Hard-coded constants.  BEWARE!
    // Clean up on aisle 3!
    // Mass of sun in solar masses.    
    const double Msun = 1.0;  // hard-code parameter.
    const double Rs_eq = 696000.0/au;
    const double J2s = 2.1106088532726840e-07;
    //const double J2s = 100000000.*2.1106088532726840e-07; //boost J2 for testing

    RAs = 268.13*M_PI/180.;
    Decs = 63.87*M_PI/180.;

    xp = cos(Decs)*cos(RAs);
    yp = cos(Decs)*sin(RAs);
    zp = sin(Decs);

    incl = acos(zp);
    if(xp != 0.0 || yp !=0.0) {    
      longnode = atan2(xp, -yp);
    } else {
      longnode = 0.0;
    }
    
    for (int j=0; j<N_real; j++){

        const struct reb_particle p = particles[j];
        double dx = p.x + (xo - xr);
        double dy = p.y + (yo - yr);
        double dz = p.z + (zo - zr);

        const double r2 = dx*dx + dy*dy + dz*dz;
        const double r = sqrt(r2);

	// Rotate to solar equatorial frame

	// Rotate around z by RA
	double cosr = cos(-longnode);
	double sinr = sin(-longnode);
	
	// Rotate around z by RA
	double dxp =  dx * cosr - dy * sinr;
	double dyp =  dx * sinr + dy * cosr;
	double dzp =  dz;

	// Rotate around x by Dec
	double cosd = cos(-incl);
	double sind = sin(-incl);

	dx =  dxp;
	dy =  dyp * cosd - dzp * sind;
	dz =  dyp * sind + dzp * cosd;

	const double costheta2 = dz*dz/r2;
        const double J2s_prefac = 3.*J2s*Rs_eq*Rs_eq/r2/r2/r/2.;
        const double J2s_fac = 5.*costheta2-1.;
        const double J2s_fac2 = 7.*costheta2-1.;
        const double J2s_fac3 = 35.*costheta2*costheta2-30.*costheta2+3.;

	// Calculate acceleration
	double resx = G*Msun*J2s_prefac*J2s_fac*dx;
	double resy = G*Msun*J2s_prefac*J2s_fac*dy;
	double resz = G*Msun*J2s_prefac*(J2s_fac-2.)*dz;

        // Variational equations

	// Rotate back to original frame
	// Rotate around x by -Dec
	double resxp =  resx;
	double resyp =  resy * cosd + resz * sind;
	double reszp = -resy * sind + resz * cosd;
	
	// Rotate around z by -RA
	resx =  resxp * cosr + resyp * sinr;
	resy = -resxp * sinr + resyp * cosr;
	resz =  reszp;

        particles[j].ax += resx;
        particles[j].ay += resy;
        particles[j].az += resz;

	// Constants for variational equations
	const double dxdx = G*Msun*J2s_prefac*(J2s_fac-5.*J2s_fac2*dx*dx/r2);
	const double dydy = G*Msun*J2s_prefac*(J2s_fac-5.*J2s_fac2*dy*dy/r2);
	const double dzdz = G*Msun*J2s_prefac*(-1.)*J2s_fac3;
	const double dxdy = G*Msun*J2s_prefac*(-5.)*J2s_fac2*dx*dy/r2;
	const double dydz = G*Msun*J2s_prefac*(-5.)*(J2s_fac2-2.)*dy*dz/r2;
	const double dxdz = G*Msun*J2s_prefac*(-5.)*(J2s_fac2-2.)*dx*dz/r2;

	// Looping over variational particles
        for(int v = N_real + 6*j; v < N_real + 6*(j+1); v++){

	    double ddx = particles[v].x;
	    double ddy = particles[v].y;
	    double ddz = particles[v].z;

	    // Rotate to solar equatorial frame
	    double ddxp =  ddx * cosr - ddy * sinr;
	    double ddyp =  ddx * sinr + ddy * cosr;
	    double ddzp =  ddz;

	    ddx =  ddxp;
	    ddy =  ddyp * cosd - ddzp * sind;
	    ddz =  ddyp * sind + ddzp * cosd;

	    double dax =   ddx * dxdx + ddy * dxdy + ddz * dxdz;
	    double day =   ddx * dxdy + ddy * dydy + ddz * dydz;
	    double daz =   ddx * dxdz + ddy * dydz + ddz * dzdz;

	    // Rotate back to original frame
	    double daxp =  dax;
	    double dayp =  day * cosd + daz * sind;
	    double dazp = -day * sind + daz * cosd;
	    dax =  daxp * cosr + dayp * sinr;
	    day = -daxp * sinr + dayp * cosr;
	    daz =  dazp;

	    // Accumulate acceleration terms
	    particles[v].ax += dax;
	    particles[v].ay += day;
	    particles[v].az += daz;

        } 
       
    }

    // Here is the treatment of non-gravitational forces.

    // The Sun center is reference for these calculations.
    double xs, ys, zs, vxs, vys, vzs, axs, ays, azs;    
    all_ephem(0, t, &GM, &xs, &ys, &zs, &vxs, &vys, &vzs, &axs, &ays, &azs);
    xr = xs;  yr = ys;  zr = zs;
    vxr = vxs; vyr = vys; vzr = vzs;    

    // The non-grav parameters are specific to each object being
    // integrated.

    // Normal asteroids
    double A1 = 0.0;
    double A2 = 0.0;
    double A3 = 0.0;

    // 2020 CD3
    //double A1= 1.903810165823E-10;
    //double A2 = 0.0;
    //double A3 = 0.0;

    // Apophis
    //double A1 = 0.0;
    //double A2 = -5.592839897872E-14;
    //double A3 = 0.0;

    // 2020 SO
    //double A1 = 2.840852439404E-9; //0.0;
    //double A2 = -2.521527931094E-10;
    //double A3= 2.317289821804E-10;
    
    for (int j=0; j<N_real; j++){

        const struct reb_particle p = particles[j];
        double dx = p.x + (xo - xr);
        double dy = p.y + (yo - yr);
        double dz = p.z + (zo - zr);

        const double r2 = dx*dx + dy*dy + dz*dz;
        const double r = sqrt(r2);

	const double g = 1.0/r2;

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
        const double t = sqrt(t2);

	particles[j].ax += A1*g*dx/r + A2*g*tx/t + A3*g*hx/h;
        particles[j].ay += A1*g*dy/r + A2*g*ty/t + A3*g*hy/h;
        particles[j].az += A1*g*dz/r + A2*g*tz/t + A3*g*hz/h;

//      variational matrix elements

        const double r3    = r*r*r;
        const double v2    = dvx*dvx + dvy*dvy + dvz*dvz;
        const double rdotv = dx*dvx  + dy*dvy  + dz*dvz;
        const double vdott = dvx*tx  + dvy*ty  + dvz*tz;

        const double dgx  = -2.*dx*g*g;
        const double dgy  = -2.*dy*g*g;
        const double dgz  = -2.*dz*g*g;

        const double hxh3 = hx/h/h/h;
        const double hyh3 = hy/h/h/h;
        const double hzh3 = hz/h/h/h;

        const double txt3 = tx/t/t/t;
        const double tyt3 = ty/t/t/t;
        const double tzt3 = tz/t/t/t;

	const double dxdx = A1*(dgx*dx + g*(1/r - dx*dx/r3)) 
                          + A2*(dgx*hx + g*(-hxh3)*(v2*dx - rdotv*dvx)) 
                          + A3*(dgx*tx + g*((dx*dvx - rdotv)/t - txt3*(2.*dx*vdott - rdotv*tx)));
	const double dydy = A1*(dgy*dy + g*(1/r - dy*dy/r3)) 
                          + A2*(dgy*hy + g*(-hyh3)*(v2*dy - rdotv*dvy)) 
                          + A3*(dgy*ty + g*((dy*dvy - rdotv)/t - tyt3*(2.*dy*vdott - rdotv*ty)));
	const double dzdz = A1*(dgz*dz + g*(1/r - dz*dz/r3)) 
                          + A2*(dgz*hz + g*(-hzh3)*(v2*dz - rdotv*dvz)) 
                          + A3*(dgz*tz + g*((dz*dvz - rdotv)/t - tzt3*(2.*dz*vdott - rdotv*tz)));

	const double dxdy = A1*(dgy*dx + g*(-dx*dy/r3))
                          + A2*(dgy*hx + g*(dvz/h -hxh3*(v2*dy - rdotv*dvy)))
                          + A3*(dgy*tx + g*((2*dy*dvx - dx*dvy)/t - txt3*(2*dy*vdott - rdotv*ty)));
	const double dydx = A1*(dgx*dy + g*(-dx*dy/r3))
                          + A2*(dgx*hy + g*(-dvz/h -hyh3*(v2*dx - rdotv*dvx)))
                          + A3*(dgx*ty + g*((2*dx*dvy - dy*dvx)/t - tyt3*(2*dx*vdott - rdotv*tx)));
	const double dxdz = A1*(dgz*dx + g*(-dx*dz/r3))
                          + A2*(dgz*hx + g*(-dvy/h -hxh3*(v2*dz - rdotv*dvz)))
                          + A3*(dgz*tx + g*((2*dz*dvx - dx*dvz)/t - txt3*(2*dz*vdott - rdotv*tz)));

	const double dzdx = A1*(dgx*dz + g*(-dx*dz/r3))
                          + A2*(dgx*hz + g*(dvy/h -hzh3*(v2*dx - rdotv*dvx)))
                          + A3*(dgx*tz + g*((2*dx*dvz - dz*dvx)/t - tzt3*(2*dx*vdott - rdotv*tx)));
	const double dydz = A1*(dgz*dy + g*(-dy*dz/r3))
                          + A2*(dgz*hy + g*(dvx/h -hyh3*(v2*dz - rdotv*dvz)))
                          + A3*(dgz*ty + g*((2*dz*dvy - dy*dvz)/t - tyt3*(2*dz*vdott - rdotv*tz)));
	const double dzdy = A1*(dgy*dz + g*(-dy*dz/r3))
                          + A2*(dgy*hz + g*(-dvx/h -hzh3*(v2*dy - rdotv*dvy)))
                          + A3*(dgy*tz + g*((2*dy*dvz - dz*dvy)/t - tzt3*(2*dy*vdott - rdotv*ty)));

	const double dxdvx = A1*(0)
                           + A2*(-hxh3*(r2*dvx - dx*rdotv))
                           + A3*((dy*dy + dz*dz)/t - txt3*r2*tx);
 	const double dydvy = A1*(0)
                           + A2*(-hyh3*(r2*dvy - dy*rdotv))
                           + A3*((dx*dx + dz*dz)/t - tyt3*r2*ty);
	const double dzdvz = A1*(0)
                           + A2*(-hzh3*(r2*dvz - dz*rdotv))
                           + A3*((dx*dx + dy*dy)/t - tzt3*r2*tz);

	const double dxdvy = A1*(0)
                           + A2*(-dz/h - hxh3*(r2*dvy - dy*rdotv))
                           + A3*(-dy*dx/t - tyt3*r2*tx);
	const double dydvx = A1*(0)
                           + A2*(dz/h - hyh3*(r2*dvx - dx*rdotv))
                           + A3*(-dx*dy/t - txt3*r2*ty);
	const double dxdvz = A1*(0)
                           + A2*(dy/h - hxh3*(r2*dvz - dz*rdotv))
                           + A3*(-dz*dx/t - tzt3*r2*tx);

	const double dzdvx = A1*(0)
                           + A2*(-dy/h - hzh3*(r2*dvx - dx*rdotv))
                           + A3*(-dx*dz/t - txt3*r2*tz);
	const double dydvz = A1*(0)
                           + A2*(-x/h - hyh3*(r2*dvz - dz*rdotv))
                           + A3*(-dz*dy/t - tzt3*r2*ty);
	const double dzdvy = A1*(0)
                           + A2*(x/h - hzh3*(r2*dvy - dy*rdotv))
                           + A3*(-dy*dz/t - tyt3*r2*tz);

	// Looping over variational particles
        for(int v = N_real + 6*j; v < N_real + 6*(j+1); v++){

	    // variational particle coords -- transformed to appropriate coord system.
	    double ddx = particles[v].x + (xo - xr);
	    double ddy = particles[v].y + (yo - yr);
	    double ddz = particles[v].z + (zo - zr);
	    double ddvx = particles[v].vx + (vxo - vxr);
	    double ddvy = particles[v].vy + (vyo - vyr);
	    double ddvz = particles[v].vz + (vzo - vzr);

	    // Matrix multiplication
	    const double dax =   ddx  * dxdx  + ddy  * dxdy  + ddz  * dxdz
		+   ddvx * dxdvx + ddvy * dxdvy + ddvz * dxdvz;
	    const double day =   ddx  * dydx  + ddy  * dydy  + ddz  * dydz
		+   ddvx * dydvx + ddvy * dydvy + ddvz * dydvz;
	    const double daz =   ddx  * dzdx  + ddy  * dzdy  + ddz  * dzdz
		+   ddvx * dzdvx + ddvy * dzdvy + ddvz * dzdvz;

	    // Accumulate acceleration terms
	    particles[v].ax += dax;
	    particles[v].ay += day;
	    particles[v].az += daz;

        }
//  variational end
    }

    // Here is the Solar GR treatment
    // The Sun is the reference for these calculations.
    xr  = xs;  yr  = ys;  zr = zs;
    vxr = vxs; vyr = vys; vzr = vzs;

//  const double mu = G*Msun; 
    const double mu = 1.0*G*Msun; // for testing: boost GR effects
    const int max_iterations = 10; // hard-coded parameter.
    for (int j=0; j<N_real; j++){

        struct reb_particle p = particles[j];
        struct reb_vec3d vi;

	p.x += (xo - xr);
	p.y += (yo - yr);
	p.z += (zo - zr);
	p.vx += (vxo - vxr);
	p.vy += (vyo - vyr);
	p.vz += (vzo - vzr);
	
        vi.x = p.vx;
        vi.y = p.vy;
        vi.z = p.vz;
        double vi2=vi.x*vi.x + vi.y*vi.y + vi.z*vi.z;
        const double ri = sqrt(p.x*p.x + p.y*p.y + p.z*p.z);
	
        int q = 0;
        double A = (0.5*vi2 + 3.*mu/ri)/C2;
        struct reb_vec3d old_v;
        for(q=0; q<max_iterations; q++){
            old_v.x = vi.x;
            old_v.y = vi.y;
            old_v.z = vi.z;
            vi.x = p.vx/(1.-A);
            vi.y = p.vy/(1.-A);
            vi.z = p.vz/(1.-A);
            vi2 =vi.x*vi.x + vi.y*vi.y + vi.z*vi.z;
            A = (0.5*vi2 + 3.*mu/ri)/C2;
            const double dvx = vi.x - old_v.x;
            const double dvy = vi.y - old_v.y;
            const double dvz = vi.z - old_v.z;
            if ((dvx*dvx + dvy*dvy + dvz*dvz)/vi2 < DBL_EPSILON*DBL_EPSILON){
                break;
            }
        }
        const int default_max_iterations = 10;
        if(q==default_max_iterations){
            reb_warning(sim, "REBOUNDx Warning: 10 iterations in ephemeris forces failed to converge. This is typically because the perturbation is too strong for the current implementation.");
        }
  
        const double B = (mu/ri - 1.5*vi2)*mu/(ri*ri*ri)/C2;
        const double rdotrdot = p.x*p.vx + p.y*p.vy + p.z*p.vz;
        
        struct reb_vec3d vidot;
        vidot.x = p.ax + B*p.x;
        vidot.y = p.ay + B*p.y;
        vidot.z = p.az + B*p.z;
        
        const double vdotvdot = vi.x*vidot.x + vi.y*vidot.y + vi.z*vidot.z;
        const double D = (vdotvdot - 3.*mu/(ri*ri*ri)*rdotrdot)/C2;

        particles[j].ax += B*(1.-A)*p.x - A*p.ax - D*vi.x;
        particles[j].ay += B*(1.-A)*p.y - A*p.ay - D*vi.y;
        particles[j].az += B*(1.-A)*p.z - A*p.az - D*vi.z;

	const double prefac = mu/(ri*ri*ri)/C2;
	const double rdotv = p.x*p.vx+p.y*p.vy+p.z*p.vz;
	const double fac1 = mu/ri-vi2;
	const double fac2 = 3.*vi2/ri/ri-4.*mu/ri/ri/ri;;
	const double fac3 = 12.*rdotv/ri/ri;

	const double dxdx = prefac*(fac1+fac2*p.x*p.x+4.*p.vx*p.vx-fac3*p.vx*p.x);
	const double dydy = prefac*(fac1+fac2*p.y*p.y+4.*p.vy*p.vy-fac3*p.vy*p.y);
	const double dzdz = prefac*(fac1+fac2*p.z*p.z+4.*p.vz*p.vz-fac3*p.vz*p.z);

	const double dxdy = prefac*(fac2*p.x*p.y+4.*p.vx*p.vy-fac3*p.vx*p.y);
	const double dydx = prefac*(fac2*p.y*p.x+4.*p.vy*p.vx-fac3*p.vy*p.x);
	const double dxdz = prefac*(fac2*p.x*p.z+4.*p.vx*p.vz-fac3*p.vx*p.z);

	const double dzdx = prefac*(fac2*p.z*p.x+4.*p.vz*p.vx-fac3*p.vz*p.x);
	const double dydz = prefac*(fac2*p.y*p.z+4.*p.vy*p.vz-fac3*p.vy*p.z);
	const double dzdy = prefac*(fac2*p.z*p.y+4.*p.vz*p.vy-fac3*p.vz*p.y);

	const double dxdvx = prefac*(4.*rdotv-2.*p.x*p.vx+4.*p.x*p.vx);
	const double dydvy = prefac*(4.*rdotv-2.*p.y*p.vy+4.*p.y*p.vy);
	const double dzdvz = prefac*(4.*rdotv-2.*p.z*p.vz+4.*p.z*p.vz);

	const double dxdvy = prefac*(-2.*p.x*p.vy+4.*p.y*p.vx);
	const double dydvx = prefac*(-2.*p.y*p.vx+4.*p.x*p.vy);
	const double dxdvz = prefac*(-2.*p.x*p.vz+4.*p.z*p.vx);

	const double dzdvx = prefac*(-2.*p.z*p.vx+4.*p.x*p.vz);
	const double dydvz = prefac*(-2.*p.y*p.vz+4.*p.z*p.vy);
	const double dzdvy = prefac*(-2.*p.z*p.vy+4.*p.y*p.vz);

	// Looping over variational particles
        for(int v = N_real + 6*j; v < N_real + 6*(j+1); v++){

	    // variational particle coords
	    double ddx = particles[v].x;
	    double ddy = particles[v].y;
	    double ddz = particles[v].z;
	    double ddvx = particles[v].vx;
	    double ddvy = particles[v].vy;
	    double ddvz = particles[v].vz;

	    // Matrix multiplication
	    const double dax =   ddx  * dxdx  + ddy  * dxdy  + ddz  * dxdz
		+   ddvx * dxdvx + ddvy * dxdvy + ddvz * dxdvz;
	    const double day =   ddx  * dydx  + ddy  * dydy  + ddz  * dydz
		+   ddvx * dydvx + ddvy * dydvy + ddvz * dydvz;
	    const double daz =   ddx  * dzdx  + ddy  * dzdy  + ddz  * dzdz
		+   ddvx * dzdvx + ddvy * dzdvy + ddvz * dzdvz;

	    // Accumulate acceleration terms
 	    particles[v].ax += dax;
	    particles[v].ay += day;
	    particles[v].az += daz;

        }
    }

    if(*geo == 1){
	// geocentric
	all_ephem(3, t, &GM, &xo, &yo, &zo, &vxo, &vyo, &vzo, &axo, &ayo, &azo);

	//printf("%lf %le %le %le geo\n", t, axo, ayo, azo);

	// This is the indirect term for geocentric equations
	// of motion.
	for (int j=0; j<N_real; j++){    

	    //printf("%lf %le %le %le\n", t, particles[j].ax, particles[j].ay, particles[j].az);	    
	    particles[j].ax -= axo;
	    particles[j].ay -= ayo;
	    particles[j].az -= azo;

	}
    }

}


struct assist_extras* assist_attach(struct reb_simulation* sim){  
    if (sim == NULL){
        fprintf(stderr, "ASSIST Error: Simulation pointer passed to assist_attach was NULL.\n");
        return NULL;
    }
    struct assist_extras* assist = malloc(sizeof(*assist));
    // Initialization separate from memory allocation because python handles memory management
    assist_initialize(sim, assist); 
    return assist;
}


void assist_extras_cleanup(struct reb_simulation* r){
    struct assist_extras* assist =  r->extras;
    assist->sim = NULL;
}

void assist_initialize(struct reb_simulation* sim, struct assist_extras* assist){
    assist->sim = sim;
    sim->extras = assist;
    sim->extras_cleanup = assist_extras_cleanup;
    sim->additional_forces = assist_additional_forces;
}

void assist_free_pointers(struct assist_extras* assist){
    assist_detach(assist->sim, assist);
    //assist->c = NULL;
    //assist->ts = NULL;
    //assist->last_state = NULL;
    //assist->hg = NULL;                
}

void assist_free(struct assist_extras* assist){
    // Freeing pointers is separate because python handles memory management of structure itself.
    assist_free_pointers(assist);
    free(assist);
}

void assist_detach(struct reb_simulation* sim, struct assist_extras* assist){
    assist->sim = NULL;
}

void assist_error(struct assist_extras* assist, const char* const msg){
    if (assist->sim == NULL){
        fprintf(stderr, "ASSIST Error: A Simulation is no longer attached to the ASSIST extras instance. Most likely the Simulation has been freed.\n");
    } else{
        reb_error(assist->sim, msg);
    }
}

int nsubsteps = 10;    

// integration_function
// tstart: integration start time in tdb
// tstep: suggested initial time step (days)
// trange: amount of time to integrate (days)
// geocentric: 1==geocentric equations of motion, 0==heliocentric
// n_particles: number of input test particles
// instate: input states of test particles
// ts: output times and states.

int integration_function(double tstart, double tend, double tstep,
			 int geocentric,
			 double epsilon,
			 int n_particles,
			 double* instate,
			 int n_var,
			 int* invar_part,			 
			 double* invar,
			 int n_alloc,			 
			 int *n_out,
			 int nsubsteps,
			 double* hg,
			 double* outtime,
			 double* outstate,
			 double min_dt){
                         //double max_dt){			 

    struct reb_simulation* sim = reb_create_simulation();

    sim->t = tstart;
    sim->dt = tstep;    // time step in days, this is just an initial value.

    // Set up simulation constants
    // The gravitational constant should be set using the ephemeris routines,
    // so that it is ensured to consistent with the units used in those routines.
    sim->G = 0.295912208285591100E-03; // Gravitational constant (AU, solar masses, days)

    sim->integrator = REB_INTEGRATOR_IAS15;
    sim->save_messages = 1;
    sim->heartbeat = heartbeat;
    sim->display_data = NULL;
    sim->collision = REB_COLLISION_NONE;  // This is important and needs to be considered carefully.
    sim->collision_resolve = reb_collision_resolve_merge; // Not sure what this is for.
    sim->gravity = REB_GRAVITY_NONE;

    // These quantities are specific to IAS15.  Perhaps something more flexible could
    // be done so that other REBOUND integration routines could be explored.

    // Don't hard code this.
    sim->ri_ias15.min_dt = min_dt;  // to avoid very small time steps (default: 0.0, suggestion 1e-2)
    //sim->ri_ias15.max_dt = max_dt;  // to avoid very large time steps (default: inf, suggestion 32.0)
    //sim->ri_ias15.epsilon = 1e-8;   // to avoid convergence issue with geocentric orbits (default: 1e-9)   

    // This should be flexible.
    // However, it should be set to 1 in most cases.
    //sim->exact_finish_time = 0;
    sim->exact_finish_time = 1;

    // Add and initialize particles    
    for(int i=0; i<n_particles; i++){

	struct reb_particle tp = {0};

	tp.x  =  instate[6*i+0];
	tp.y  =  instate[6*i+1];
	tp.z  =  instate[6*i+2];
	tp.vx =  instate[6*i+3];
	tp.vy =  instate[6*i+4];
	tp.vz =  instate[6*i+5];

	reb_add(sim, tp);
    }

    // Add and initialize variational particles
    for(int i=0; i<n_var; i++){

	// invar_part[i] contains the index of the test particle that we vary.	
        int var_i = reb_add_var_1st_order(sim, invar_part[i]);
	
        sim->particles[var_i].x =  invar[6*i+0]; 
        sim->particles[var_i].y =  invar[6*i+1]; 
        sim->particles[var_i].z =  invar[6*i+2]; 
        sim->particles[var_i].vx = invar[6*i+3]; 
        sim->particles[var_i].vy = invar[6*i+4]; 
        sim->particles[var_i].vz = invar[6*i+5]; 
 
    }
    
    struct assist_extras* assist = assist_attach(sim);

    timestate *ts = (timestate*) malloc(sizeof(timestate));
    tstate* last_state = (tstate*) malloc(sim->N*sizeof(tstate));

    ts->t = outtime;
    ts->state = outstate;        

    assist->last_state = last_state;
    assist->ts = ts;

    assist->nsubsteps = nsubsteps;
    assist->hg = hg;

    // Need to read in parameters geocentric, c,

    int N = sim->N; // N includes real+variational particles

    ts->n_particles = n_particles;
    ts->n_alloc = n_alloc;

    reb_integrate(sim, tend);

    if (sim->messages){
	printf("error\n");
	fflush(stdout);
	for(int i=0; i<reb_max_messages_N; i++){
	    printf("mess: %d", i);
	    fflush(stdout);
	    if(sim->messages[i] != NULL){
		printf("%d %s", i, sim->messages[i]);
		printf("blah\n");
		fflush(stdout);
	    }
	}
    }

    //struct reb_particle* const particles = sim->particles;

    *n_out = sim->steps_done;

    int status = sim->status;

    printf("%p %p\n", assist, ts);
    fflush(stdout);

    // explicitly free all the memory allocated by ASSIST
    ts->t = NULL;
    ts->state = NULL;
    ts->last_state = NULL;
    free(ts);
    free(last_state);

    //assist_free(assist);    // this explicitly frees all the memory allocated by ASSIST
    
    reb_free_simulation(sim);

    return(status);
}

// This function is doing two related things:
// 1. Calculating the positions and velocities at the substeps
// 2. Storing the times and positions/velocities in the arrays
//    that are provided.
// For this to work, we need:
// * the last valid state for all particles,
// * the b coefficients for all the particles,
// * the last time step
//
// We need to adjust this so that it stores the positions
// and velocities at the substeps and the final computed
// state, rather than the previous computed state and
// the values at the substeps.

//static const double hg[11]   =   { 0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0};

void store_function(struct reb_simulation* sim){
    int N = sim->N;
    int N3 = 3*N;

    static int last_steps_done = 0;

    double s[9]; // Summation coefficients

    struct assist_extras* assist = (struct assist_extras*) sim->extras;

    int nsubsteps = assist->nsubsteps;
    double* hg = assist->hg;

    timestate* ts = ((struct assist_extras*) sim->extras)->ts;
    tstate* last_state = ((struct assist_extras*) sim->extras)->last_state;    

    static double* outtime;
    static double* outstate;

    int n_alloc;

    int step = sim->steps_done;

    outtime = ts->t;
    outstate = ts->state;
    n_alloc= ts->n_alloc;

    if(step==0){

	int state_offset = 0;
	int time_offset = 0;
	
	outtime[time_offset++] = sim->t;

	for(int j=0; j<N; j++){
	    last_state[j].t = sim->t;	
	    outstate[state_offset++] = sim->particles[j].x;
	    outstate[state_offset++] = sim->particles[j].y;
	    outstate[state_offset++] = sim->particles[j].z;
	    outstate[state_offset++] = sim->particles[j].vx;
	    outstate[state_offset++] = sim->particles[j].vy;
	    outstate[state_offset++] = sim->particles[j].vz;
	}

    }else if(sim->steps_done > last_steps_done){

	// Convenience variable.  The 'br' field contains the 
	// set of coefficients from the last completed step.
	const struct reb_dpconst7 b  = dpcast(sim->ri_ias15.br);

	double* x0 = malloc(sizeof(double)*N3);
	double* v0 = malloc(sizeof(double)*N3);
	double* a0 = malloc(sizeof(double)*N3);

	for(int j=0;j<N;j++) {

	    const int k0 = 3*j+0;
	    const int k1 = 3*j+1;
	    const int k2 = 3*j+2;

	    x0[k0] = last_state[j].x;
	    x0[k1] = last_state[j].y;
	    x0[k2] = last_state[j].z;

	    v0[k0] = last_state[j].vx;
	    v0[k1] = last_state[j].vy;
	    v0[k2] = last_state[j].vz;	

	    a0[k0] = last_state[j].ax;
	    a0[k1] = last_state[j].ay;
	    a0[k2] = last_state[j].az;
	    
	}
	int time_offset = (step-1)*nsubsteps+1;

	// Loop over intervals using specified spacings      
	for(int n=1;n<(nsubsteps+1);n++) {	    

	    // The hg[n] values here define the substeps used in the
	    // the integration, but they could be any values.
	    // A natural alternative would be the Chebyshev nodes.
	    // The degree would be altered, but the value would need
	    // to be included in the output.
	    // Another approach would be to output the ingredients for
	    // the evaluations below.
	    // x0, v0, a0, and 7 coefficients for each component
	    // This has the advantage of providing a complete state
	    // plus the accelerations at intervals.

	    s[0] = sim->dt_last_done * hg[n];

	    s[1] = s[0] * s[0] / 2.;
	    s[2] = s[1] * hg[n] / 3.;
	    s[3] = s[2] * hg[n] / 2.;
	    s[4] = 3. * s[3] * hg[n] / 5.;
	    s[5] = 2. * s[4] * hg[n] / 3.;
	    s[6] = 5. * s[5] * hg[n] / 7.;
	    s[7] = 3. * s[6] * hg[n] / 4.;
	    s[8] = 7. * s[7] * hg[n] / 9.;

	    double t = sim->t + sim->dt_last_done * (-1.0 + hg[n]);

	    outtime[time_offset++] = t;

	    // Predict positions at interval n using b values
	    // for all the particles
	    for(int j=0;j<N;j++) {  
		const int k0 = 3*j+0;
		const int k1 = 3*j+1;
		const int k2 = 3*j+2;

		double xx0 = x0[k0] + (s[8]*b.p6[k0] + s[7]*b.p5[k0] + s[6]*b.p4[k0] + s[5]*b.p3[k0] + s[4]*b.p2[k0] + s[3]*b.p1[k0] + s[2]*b.p0[k0] + s[1]*a0[k0] + s[0]*v0[k0] );
		double xy0 = x0[k1] + (s[8]*b.p6[k1] + s[7]*b.p5[k1] + s[6]*b.p4[k1] + s[5]*b.p3[k1] + s[4]*b.p2[k1] + s[3]*b.p1[k1] + s[2]*b.p0[k1] + s[1]*a0[k1] + s[0]*v0[k1] );
		double xz0 = x0[k2] + (s[8]*b.p6[k2] + s[7]*b.p5[k2] + s[6]*b.p4[k2] + s[5]*b.p3[k2] + s[4]*b.p2[k2] + s[3]*b.p1[k2] + s[2]*b.p0[k2] + s[1]*a0[k2] + s[0]*v0[k2] );

		// Store the results
		int offset = ((step-1)*nsubsteps + n)*6*N + 6*j;
		outstate[offset+0] = xx0;
		outstate[offset+1] = xy0;	  	  
		outstate[offset+2] = xz0;
	    }

	    s[0] = sim->dt_last_done * hg[n];
	    s[1] =      s[0] * hg[n] / 2.;
	    s[2] = 2. * s[1] * hg[n] / 3.;
	    s[3] = 3. * s[2] * hg[n] / 4.;
	    s[4] = 4. * s[3] * hg[n] / 5.;
	    s[5] = 5. * s[4] * hg[n] / 6.;
	    s[6] = 6. * s[5] * hg[n] / 7.;
	    s[7] = 7. * s[6] * hg[n] / 8.;

	    // Predict velocities at interval n using b values
	    // for all the particles
	    for(int j=0;j<N;j++) {

		const int k0 = 3*j+0;
		const int k1 = 3*j+1;
		const int k2 = 3*j+2;

		double vx0 = v0[k0] + s[7]*b.p6[k0] + s[6]*b.p5[k0] + s[5]*b.p4[k0] + s[4]*b.p3[k0] + s[3]*b.p2[k0] + s[2]*b.p1[k0] + s[1]*b.p0[k0] + s[0]*a0[k0];
		double vy0 = v0[k1] + s[7]*b.p6[k1] + s[6]*b.p5[k1] + s[5]*b.p4[k1] + s[4]*b.p3[k1] + s[3]*b.p2[k1] + s[2]*b.p1[k1] + s[1]*b.p0[k1] + s[0]*a0[k1];
		double vz0 = v0[k2] + s[7]*b.p6[k2] + s[6]*b.p5[k2] + s[5]*b.p4[k2] + s[4]*b.p3[k2] + s[3]*b.p2[k2] + s[2]*b.p1[k2] + s[1]*b.p0[k2] + s[0]*a0[k2];

		// Store the results
		int offset = ((step-1)*nsubsteps + n)*6*N + 6*j;				
		outstate[offset+3] = vx0;
		outstate[offset+4] = vy0;	  	  
		outstate[offset+5] = vz0;

	    }
	}

	free(x0);
	free(v0);
	free(a0);
    }
    last_steps_done = sim->steps_done;

    if((ts->n_alloc-step) < 1){
	sim->status = REB_EXIT_USER;
	return;
    }

}

void store_coefficients(struct reb_simulation* sim){
    int N = sim->N;
    int N3 = 3*N;

    static int last_steps_done = 0;

    double s[9]; // Summation coefficients

    struct assist_extras* assist = (struct assist_extras*) sim->extras;

    int nsubsteps = assist->nsubsteps;
    double* hg = assist->hg;

    timestate* ts = ((struct assist_extras*) sim->extras)->ts;
    tstate* last_state = ((struct assist_extras*) sim->extras)->last_state;    

    static double* outtime;
    static double* outstate;

    int n_alloc;

    int step = sim->steps_done;

    outtime = ts->t;
    outstate = ts->state;
    n_alloc= ts->n_alloc;

    if(step==0){

	//printf("initial step %d %lf\n", step, sim->t);
	
    }else if(step > last_steps_done){

	//printf("step %d\n", step);

	// Convenience variable.  The 'br' field contains the 
	// set of coefficients from the last completed step.
	const struct reb_dpconst7 b  = dpcast(sim->ri_ias15.br);

	double* x0 = malloc(sizeof(double)*N3);
	double* v0 = malloc(sizeof(double)*N3);
	double* a0 = malloc(sizeof(double)*N3);

	//double t = sim->t + sim->dt_last_done * (-1.0 + hg[n]);
	printf("x %lf %lf \n", sim->t-sim->dt_last_done, sim->dt_last_done);
	
	for(int j=0;j<N;j++) {

	    const int k0 = 3*j+0;
	    const int k1 = 3*j+1;
	    const int k2 = 3*j+2;

	    x0[k0] = last_state[j].x;
	    x0[k1] = last_state[j].y;
	    x0[k2] = last_state[j].z;

	    v0[k0] = last_state[j].vx;
	    v0[k1] = last_state[j].vy;
	    v0[k2] = last_state[j].vz;	

	    a0[k0] = last_state[j].ax;
	    a0[k1] = last_state[j].ay;
	    a0[k2] = last_state[j].az;

	    printf("%d %.16le %.16le %.16le\n",
		   j,
		   x0[k0], v0[k0], a0[k0]);
	    printf("%d %.16le %.16le %.16le %.16le %.16le %.16le %.16le\n",
		   j, b.p0[k0], b.p1[k0], b.p2[k0], b.p3[k0], b.p4[k0], b.p5[k0], b.p6[k0]);
	    printf("%d %.16le %.16le %.16le %.16le %.16le %.16le %.16le\n",
		   j, b.p0[k1], b.p1[k1], b.p2[k1], b.p3[k1], b.p4[k1], b.p5[k1], b.p6[k1]);
	    printf("%d %.16le %.16le %.16le %.16le %.16le %.16le %.16le\n",
		   j, b.p0[k2], b.p1[k2], b.p2[k2], b.p3[k2], b.p4[k2], b.p5[k2], b.p6[k2]);	    

	}

	free(x0);
	free(v0);
	free(a0);
    }

}

void store_last_state(struct reb_simulation* sim){

    //timestate* ts = ((struct assist_extras*) sim->extras)->ts;
    tstate* last_state = ((struct assist_extras*) sim->extras)->last_state;    
    
    int N = sim->N;    
    for(int j=0; j<N; j++){ 
	last_state[j].t = sim->t;	
	last_state[j].x = sim->particles[j].x;
	last_state[j].y = sim->particles[j].y;
	last_state[j].z = sim->particles[j].z;
	last_state[j].vx = sim->particles[j].vx;
	last_state[j].vy = sim->particles[j].vy;
	last_state[j].vz = sim->particles[j].vz;
	last_state[j].ax = sim->particles[j].ax;
	last_state[j].ay = sim->particles[j].ay;
	last_state[j].az = sim->particles[j].az;
    }
}

void heartbeat(struct reb_simulation* sim){

    void store_function(struct reb_simulation* sim);
    void store_last_state(struct reb_simulation* sim);
    void store_coefficients(struct reb_simulation* sim);        

    store_function(sim);
    store_coefficients(sim);    

    reb_update_acceleration(sim);

    store_last_state(sim);

}
