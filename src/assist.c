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

	    //fprintf(outfile, "%3d %25.16le %25.16le %25.16le %25.16le %25.16le %25.16le %25.16le %25.16le\n", i, t, GM, dx, dy, dz, -prefac*dx, -prefac*dy, -prefac*dz);

	    particles[j].ax -= prefac*dx;
	    particles[j].ay -= prefac*dy;
	    particles[j].az -= prefac*dz;

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
			 double min_dt,
			 double max_dt){			 

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
    sim->ri_ias15.max_dt = max_dt;  // to avoid very large time steps (default: inf, suggestion 32.0)
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

    // explicitly free all the memory allocated by REBOUNDx
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

	// Loop over intervals using Gauss-Radau spacings      
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

void heartbeat(struct reb_simulation* r){

    void store_function(struct reb_simulation* r);
    void store_last_state(struct reb_simulation* r);    

    store_function(r);

    reb_update_acceleration(r);

    store_last_state(r);

}
