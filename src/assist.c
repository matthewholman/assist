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

#include "const.h"

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <stdint.h>
#include <math.h>
#include <limits.h>
#include <float.h>
#include "assist.h"
#include "rebound.h"

#include "spk.h"
#include "planets.h"
#include "forces.h"

const int reb_max_messages_length = 1024;   // needs to be constant expression for array size
const int reb_max_messages_N = 10;

enum {
    NO_ERR,        // no error
    ERR_JPL_EPHEM, // JPL ephemeris file not found
    ERR_JPL_AST,   // JPL asteroid file not found
    ERR_NAST,      // asteroid number out of range
    ERR_NEPH,      // planet number out of range
};

#define STRINGIFY(s) str(s)
#define str(s) #s

const char* assist_build_str = __DATE__ " " __TIME__;   // Date and time build string. 
const char* assist_version_str = "1.0.1b7";         // **VERSIONLINE** This line gets updated automatically. Do not edit manually.
const char* assist_githash_str = STRINGIFY(ASSISTGITHASH);// This line gets updated automatically. Do not edit manually.
    
// Forward function declarations
static void store_function(struct reb_simulation* sim);
static void assist_heartbeat(struct reb_simulation* r);

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
    sim->force_is_velocity_dependent = 1;
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

// integration_function
// tstart: integration start time in tdb
// tend:   integration end time in tdb
// tstep:  suggested initial time step (days)
// geocentric:  1==geocentric equations of motion, 0==barycentric
// n_particles: number of input test particles
// instate:     input states of test particles
// part_params: input non-grav parameters for test particles
// n_var:       number of input variational particles
// invar_part:  index of host particle that each variational
//              particle refers to.
// invar:       input states of variational particles
// var_part_params: input non-grav parameters for variational particles
// n_alloc:     number of overall times steps for which there
//              is space allocated.
// n_out:       number of outputs
// nsubsteps:   number of substeps of output per overall
//              time step.
// hg:          array of output substep times as a fraction of
//              the step interval, i.e. values 0 to 1.
// outtime:     array of output times.
// outstate:    array of output states.
// min_dt:      minimum allowed time step.
int integration_function(double jd_ref,
			 double tstart, double tend, double tstep,
			 int geocentric,
			 double epsilon,
			 int n_particles,
			 double* instate,
			 double* part_params,			 
			 int n_var,
			 int* invar_part,			 
			 double* invar,
			 double* var_part_params, // particle constants for variational particles too			 
			 int n_alloc,			 
			 int *n_out,
			 int nsubsteps,
			 double* hg,
			 double* outtime,
			 double* outstate,
			 double min_dt){

    //const double au = JPL_EPHEM_CAU;    
    //const double c = (JPL_EPHEM_CLIGHT/au)*86400;
    
    struct reb_simulation* sim = reb_create_simulation();

    sim->t = tstart;
    sim->dt = tstep;    // time step in days, this is just an initial value.

    // TODO: decide how flexible these should be.
    sim->integrator = REB_INTEGRATOR_IAS15;
    sim->save_messages = 1;
    sim->heartbeat = assist_heartbeat;
    sim->display_data = NULL;
    sim->collision = REB_COLLISION_NONE;  // This is important and needs to be considered carefully.
    sim->collision_resolve = reb_collision_resolve_merge; // Not sure what this is for.
    sim->gravity = REB_GRAVITY_NONE;

    // This should be flexible.
    // However, it should be set to 1 in most cases.
    //sim->exact_finish_time = 0;
    sim->exact_finish_time = 1;
    
    // These quantities are specific to IAS15.  Perhaps something more flexible could
    // be done so that other REBOUND integration routines could be explored.

    sim->ri_ias15.min_dt = min_dt;    // to avoid very small time steps (default: 0.0, suggestion 1e-2)
    sim->ri_ias15.epsilon = epsilon;  // to avoid convergence issue with geocentric orbits (default: 1e-9)

    // Attach an assist struct to the simulation
    // and initialize some values.
    struct assist_extras* assist = assist_attach(sim);
    assist->particle_params = NULL;
    assist->N = 0;
    assist->geocentric = geocentric;
    assist->jd_ref = jd_ref; // 2451545.0; 

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

	// Could probably tie the next statements together
	// in one function with the previous call.  For each
	// real or variational particle added to rebound, there
	// needs to be a set of extra parameters.
	// The number of non-grav parameters should be flexible,
	// rather than fixed at 3.	
	assist->N++;
	// If part_params is NULL, skip this part
	if(part_params != NULL){
	    assist->particle_params = realloc(assist->particle_params, (assist->N)*3*sizeof(double));
	    assist->particle_params[3*i+0] = part_params[3*i+0];
	    assist->particle_params[3*i+1] = part_params[3*i+1];
	    assist->particle_params[3*i+2] = part_params[3*i+2];
	}
    }

    // Need to ensure that if n_var != 0 that
    // invar != NULL.
    // update the return value
    if(n_var != 0 && invar == NULL)
	return 0;

    // Add and initialize variational particles
    for(int i=0; i<n_var; i++){

	// invar_part[i] contains the index of the test particle that we vary.
        int var_i = reb_add_var_1st_order(sim, invar_part[i]);
	//assist_add(assist, params);	
	
        sim->particles[var_i].x =  invar[6*i+0]; 
        sim->particles[var_i].y =  invar[6*i+1]; 
        sim->particles[var_i].z =  invar[6*i+2]; 
        sim->particles[var_i].vx = invar[6*i+3]; 
        sim->particles[var_i].vy = invar[6*i+4]; 
        sim->particles[var_i].vz = invar[6*i+5];

	// Could probably tie the next statements together
	// in one function with the previous call.  For each
	// real or variational particle added to rebound, there
	// needs to be a set of extra parameters.
	// The number of non-grav parameters should be flexible,
	// rather than fixed at 3.	
	assist->N++;
	// If var_part_params is null, skip this part.
	if(var_part_params != NULL){	
	    assist->particle_params = realloc(assist->particle_params, (assist->N)*3*sizeof(double));
	    assist->particle_params[3*var_i+0] = var_part_params[3*i+0];
	    assist->particle_params[3*var_i+1] = var_part_params[3*i+1];
	    assist->particle_params[3*var_i+2] = var_part_params[3*i+2];
	}
    }

    // Allocate memory.
    timestate *ts = (timestate*) malloc(sizeof(timestate));
    tstate* last_state = (tstate*) malloc(sim->N*sizeof(tstate));

    ts->t = outtime;
    ts->state = outstate;        

    assist->last_state = last_state;
    assist->ts = ts;

    assist->nsubsteps = nsubsteps;
    assist->hg = hg;

    ts->n_particles = n_particles;
    ts->n_alloc = n_alloc;

    // Do the integration
    reb_integrate(sim, tend);

    // This should be handled different for python.
    if (sim->messages){
	for(int i=0; i<reb_max_messages_N; i++){
	    printf("mess: %d", i);
	    if(sim->messages[i] != NULL){
		printf("%d %s", i, sim->messages[i]);
	    }
	}
    }

    *n_out = sim->steps_done;

    int status = sim->status;

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

static void store_function(struct reb_simulation* sim){
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

    int step = sim->steps_done;

    outtime = ts->t;
    outstate = ts->state;

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

static void store_last_state(struct reb_simulation* sim){

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

static void assist_heartbeat(struct reb_simulation* sim){
    store_function(sim);
    reb_update_acceleration(sim);
    store_last_state(sim);
}



