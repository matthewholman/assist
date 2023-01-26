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
static void assist_heartbeat(struct reb_simulation* r);
static void assist_pre_timestep_modifications(struct reb_simulation* r);

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

//"path/to/planet/ephem", "path/to/smallbody/ephem"
struct assist_ephem* assist_ephem_init(char *user_planets_path, char *user_asteroids_path){

    char default_planets_path[] = "/data/linux_m13000p17000.441";
    char default_asteroids_path[] = "/data/sb441-n16.bsp";

    struct assist_ephem* ephem = calloc(1, sizeof(struct assist_ephem));
    ephem->jd_ref = 2451545.0; // Default jd_ref
    
    const int FNAMESIZE = 1024;
    char planets_path[FNAMESIZE];
    char asteroids_path[FNAMESIZE];        

    /** Use user-defined file or the default filename, 
     *  in that order.
     */

    if(user_planets_path == NULL && getenv("ASSIST_DIR")==NULL){
        fprintf(stderr, "No user or default planet ephemeris file\n");
        return(NULL);	  
    }

    if(user_planets_path == NULL){
        sprintf(planets_path, "%s%s", getenv("ASSIST_DIR"), default_planets_path);
    }else{
        strncpy(planets_path, user_planets_path, FNAMESIZE-1);	
    }

    if ((ephem->pl = assist_jpl_init(planets_path)) == NULL) {
        printf("Couldn't find planet ephemeris file: %s\n", planets_path);	  
        return(NULL);	  
    }

    if(user_asteroids_path == NULL && getenv("ASSIST_DIR")==NULL){
        fprintf(stderr, "No user or asteroid ephemeris file\n");
        return(NULL);	  
    }

    if(user_asteroids_path == NULL){
        sprintf(asteroids_path, "%s%s", getenv("ASSIST_DIR"), default_asteroids_path);
    }else{
        strncpy(asteroids_path, user_asteroids_path, FNAMESIZE-1);	
    }

    if ((ephem->spl = assist_spk_init(asteroids_path)) == NULL) {
        printf("Couldn't find asteroid ephemeris file: %s\n", asteroids_path);
        return(NULL);
    }

    return(ephem);
}

void assist_ephem_free(struct assist_ephem* ephem){
    assist_jpl_free(ephem->pl);
    assist_spk_free(ephem->spl);
    free(ephem);
}

struct assist_extras* assist_attach(struct reb_simulation* sim, struct assist_ephem* ephem){  
    if (sim == NULL){
        fprintf(stderr, "ASSIST Error: Simulation pointer passed to assist_attach was NULL.\n");
        return NULL;
    }
    int extras_should_free_ephem = 0;
    if (ephem == NULL){
        // Try default 
        ephem = assist_ephem_init(NULL, NULL);
        if (ephem == NULL){
            fprintf(stderr, "ASSIST Error: Ephemeris pointer passed to assist_attach was NULL. Initialization with default path failed.\n");
            return NULL;
        }
        extras_should_free_ephem = 1;
    }

    // Initialization separate from memory allocation because python handles memory management
    struct assist_extras* assist = calloc(1, sizeof(*assist));
    assist_initialize(sim, assist, ephem); 
    assist->extras_should_free_ephem = extras_should_free_ephem;
    
    return assist;
}


void assist_initialize(struct reb_simulation* sim, struct assist_extras* assist, struct assist_ephem* ephem){
    assist->sim = sim;
    assist->ephem_cache = calloc(1, sizeof(struct assist_ephem_cache));
    assist->ephem_cache->sun_t = NAN; 
    const int N_total = 16+11; // TODO
    assist->ephem_cache->items = calloc(N_total*7, sizeof(struct assist_cache_item));
    assist->ephem_cache->index = malloc(N_total*sizeof(int));
    assist->ephem_cache->t = malloc(N_total*7*sizeof(double));
    for (int i=0;i<N_total;i++){
        assist->ephem_cache->index[i] = 6;
    }
    for (int i=0;i<7*N_total;i++){
        assist->ephem_cache->t[i] = NAN;
    }

    assist->ephem = ephem;
    assist->particle_params = NULL;
    assist->forces = ASSIST_FORCE_SUN // default forces
                     | ASSIST_FORCE_PLANETS
                     | ASSIST_FORCE_ASTEROIDS
                     | ASSIST_FORCE_NON_GRAVITATIONAL
                     | ASSIST_FORCE_EARTH_HARMONICS
                     | ASSIST_FORCE_SUN_HARMONICS
                     | ASSIST_FORCE_GR_EIH;
    
    sim->integrator = REB_INTEGRATOR_IAS15;
    sim->gravity = REB_GRAVITY_NONE;
    sim->extras = assist;
    sim->additional_forces = assist_additional_forces;
    sim->force_is_velocity_dependent = 1;
}

void assist_free_pointers(struct assist_extras* assist){
    assist_detach(assist->sim, assist);
    free(assist->last_state_x);
    free(assist->last_state_v);
    free(assist->last_state_a);
    free(assist->ephem_cache);
    if (assist->extras_should_free_ephem){
        assist_ephem_free(assist->ephem);
    }
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


struct reb_particle assist_get_particle(struct assist_ephem* ephem, const int particle_id, const double t){
    struct reb_particle p = {0};
    double GM = 0;
    int flag = assist_all_ephem(ephem, NULL, particle_id, t, &GM, &p.x, &p.y, &p.z, &p.vx, &p.vy, &p.vz, &p.ax, &p.ay, &p.az);
    if (flag != NO_ERR){
        fprintf(stderr, "An error occured while trying to initialize particle from ephemeris data.\n");
    }
    p.m = GM; // Note this is GM, not M
    return p;
}


// assist_integrate
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
// output_n_alloc:     number of overall times steps for which there
//              is space allocated.
// n_out:       number of outputs
// nsubsteps:   number of substeps of output per overall
//              time step.
// hg:          array of output substep times as a fraction of
//              the step interval, i.e. values 0 to 1.
// output_t:      array of output times.
// output_state:  array of output states.
// min_dt:      minimum allowed time step.
int assist_integrate(struct assist_ephem* ephem,
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
		      int output_n_alloc,			 
		      int *n_out,
		      int nsubsteps,
		      double* hg,
		      double* output_t,
		      double* output_state,
		      double min_dt){

    //const double au = JPL_EPHEM_CAU;    
    //const double c = (JPL_EPHEM_CLIGHT/au)*86400;

    struct reb_simulation* sim = reb_create_simulation();

    sim->t = tstart;
    sim->dt = tstep;    // time step in days, this is just an initial value.

    // These quantities are specific to IAS15.  Perhaps something more flexible could
    // be done so that other REBOUND integration routines could be explored.

    sim->ri_ias15.min_dt = min_dt;    // to avoid very small time steps (default: 0.0, suggestion 1e-2)
    sim->ri_ias15.epsilon = epsilon;  // to avoid convergence issue with geocentric orbits (default: 1e-9)
    
    sim->heartbeat = assist_heartbeat;
    sim->pre_timestep_modifications = assist_pre_timestep_modifications;
    sim->save_messages = 1;

    // Attach an assist struct to the simulation
    // and initialize some values.
    struct assist_extras* assist = assist_attach(sim, ephem);
    assist->geocentric = geocentric;

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

    int N = sim->N;
	// If part_params is NULL, skip this part
	if(part_params != NULL){
	    assist->particle_params = realloc(assist->particle_params, N*3*sizeof(double));
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
    int N = sim->N;
	// If var_part_params is null, skip this part.
	if(var_part_params != NULL){	
	    assist->particle_params = realloc(assist->particle_params, N*3*sizeof(double));
	    assist->particle_params[3*var_i+0] = var_part_params[3*i+0];
	    assist->particle_params[3*var_i+1] = var_part_params[3*i+1];
	    assist->particle_params[3*var_i+2] = var_part_params[3*i+2];
	}
    }

    // Allocate memory.

    assist->output_t = output_t;
    assist->output_state = output_state;
    assist->output_n_alloc = output_n_alloc;
    assist->last_state_x = malloc(sim->N*3*sizeof(double));
    assist->last_state_v = malloc(sim->N*3*sizeof(double));
    assist->last_state_a = malloc(sim->N*3*sizeof(double));

    assist->nsubsteps = nsubsteps;
    assist->hg = hg;

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

    assist_free(assist);    // this explicitly frees all the memory allocated by ASSIST
    
    reb_free_simulation(sim);

    return(status);
}

// Fill array output with data using interpolation over the past timestep
// h=0 beginning of timestep, h=1 end of timestep
int assist_interpolate(struct reb_simulation* sim, double h, double* output){
    struct assist_extras* assist = (struct assist_extras*) sim->extras;
    int N = sim->N;

    // Convenience variable.  The 'br' field contains the
    // set of coefficients from the last completed step.
    const struct reb_dpconst7 b  = dpcast(sim->ri_ias15.br);

    double* x0 = assist->last_state_x;
    double* v0 = assist->last_state_v;
    double* a0 = assist->last_state_a;

    double s[9]; // Summation coefficients

    s[0] = sim->dt_last_done * h;

    s[1] = s[0] * s[0] / 2.;
    s[2] = s[1] * h / 3.;
    s[3] = s[2] * h / 2.;
    s[4] = 3. * s[3] * h / 5.;
    s[5] = 2. * s[4] * h / 3.;
    s[6] = 5. * s[5] * h / 7.;
    s[7] = 3. * s[6] * h / 4.;
    s[8] = 7. * s[7] * h / 9.;

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
        int offset = 6*j;
        output[offset+0] = xx0;
        output[offset+1] = xy0;
        output[offset+2] = xz0;
    }

    s[0] = sim->dt_last_done * h;
    s[1] =      s[0] * h / 2.;
    s[2] = 2. * s[1] * h / 3.;
    s[3] = 3. * s[2] * h / 4.;
    s[4] = 4. * s[3] * h / 5.;
    s[5] = 5. * s[4] * h / 6.;
    s[6] = 6. * s[5] * h / 7.;
    s[7] = 7. * s[6] * h / 8.;

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
        int offset = 6*j;
        output[offset+3] = vx0;
        output[offset+4] = vy0;
        output[offset+5] = vz0;

    }
    return 1;
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

static void assist_heartbeat(struct reb_simulation* sim){
    int N = sim->N;

    struct assist_extras* assist = (struct assist_extras*) sim->extras;

    int steps_done = assist->steps_done;
    int nsubsteps = assist->nsubsteps;

    double* output_t = assist->output_t;
    double* output_state = assist->output_state;

    int step = sim->steps_done;

    if(step==0){ // first output are just the initial conditions
        int state_offset = 0;

        output_t[0] = sim->t;

        for(int j=0; j<N; j++){
            output_state[state_offset++] = sim->particles[j].x;
            output_state[state_offset++] = sim->particles[j].y;
            output_state[state_offset++] = sim->particles[j].z;
            output_state[state_offset++] = sim->particles[j].vx;
            output_state[state_offset++] = sim->particles[j].vy;
            output_state[state_offset++] = sim->particles[j].vz;
        }

    }else if(sim->steps_done > steps_done){

        double* hg = assist->hg;

        // Loop over substeps
        for(int n=1;n<(nsubsteps+1);n++) {	    
            int offset = (step-1)*nsubsteps;
            output_t[offset+n] = sim->t + sim->dt_last_done * (-1.0 + hg[n]);
            assist_interpolate(sim, hg[n], output_state +  (offset + n)*6*N );
        }
    }
    steps_done = sim->steps_done;

    if((assist->output_n_alloc-step*nsubsteps) < 1){
        sim->status = REB_EXIT_USER;
        return;
    }

}

static void assist_pre_timestep_modifications(struct reb_simulation* sim){
    struct assist_extras* assist = sim->extras;
    assist->last_state_t = sim->t;

    reb_update_acceleration(sim); // This will later be recalculated. Could be optimized.

    for(int j=0; j<sim->N; j++){ 
        int offset = 3*j;
        assist->last_state_x[offset+0] = sim->particles[j].x;
        assist->last_state_x[offset+1] = sim->particles[j].y;
        assist->last_state_x[offset+2] = sim->particles[j].z;
        assist->last_state_v[offset+0] = sim->particles[j].vx;
        assist->last_state_v[offset+1] = sim->particles[j].vy;
        assist->last_state_v[offset+2] = sim->particles[j].vz;
        assist->last_state_a[offset+0] = sim->particles[j].ax;
        assist->last_state_a[offset+1] = sim->particles[j].ay;
        assist->last_state_a[offset+2] = sim->particles[j].az;
    }
}

