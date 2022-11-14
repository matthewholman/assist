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
#define DEFAULT_JPL_SB_EPHEM "../../data/sb441-n16.bsp"
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

static int ephem(const int i, const double jde, double* const GM,
	   double* const x, double* const y, double* const z,
	   double* const vx, double* const vy, double* const vz,
	   double* const ax, double* const ay, double* const az){

    static int initialized = 0;

    static struct _jpl_s *pl;
    struct mpos_s now;

    // Calculate GM values for Earth and Moon
    // from Earth-moon ratio and sum.
    const static double em_r = JPL_EPHEM_EMRAT;
    const static double GMe = em_r/(1.+em_r) * JPL_EPHEM_GMB;
    const static double GMm = 1./(1.+em_r) * JPL_EPHEM_GMB;

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

    /*
    const static double JPL_GM[16] =    
    {
	    3.2191392075878588e-15, // 107 camilla
	    JPL_EPHEM_MA0001, // 1 Ceres
	    //1.3964518123081070e-13, // 1 ceres
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
    */
    
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

	FILE *file;
	file = fopen(buf, "r");
	if(file == NULL){
	    fprintf(stderr, "Couldn't open asteroid file: %s\n", buf);
	}

	if ((spl = spk_init(buf)) == NULL) {
	    printf("Could find asteroid ephemeris file: %s\n", buf);
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

// TODO: this should be more general.  Perhaps
// this information could come directly from
// the JPL ephemerides.
static int number_bodies(int* N_ephem, int* N_ast){
    *N_ephem = 11;
    *N_ast = 16;

    return(*N_ephem + *N_ast);
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
    // all_ephem will need to access the positions and GM values
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

//void assist_add(struct assist_extras* assist, params);

void assist_additional_forces(struct reb_simulation* sim){

    // implement additional_forces here
    //printf("here\n");
    //fflush(stdout);

    struct assist_extras* assist = (struct assist_extras*) sim->extras;

    //printf("there\n");
    //fflush(stdout);

    // TODO: Constant should be hard-codes.
    // geo flag should be set from the outside.

    const unsigned int N = sim->N;  // N includes real+variational particles
    const unsigned int N_real = N - sim->N_var;    

    int geo = assist->geocentric;

    //printf("here geo: %d\n", geo);

    int N_ephem, N_ast;
    int number_bodies(int* N_ephem, int* N_ast);

    const int N_tot = number_bodies(&N_ephem, &N_ast);    

    // The limit of the EIH GR limit should be a free
    // parameter
    int eih_loop_limit = N_ephem; // 1; // N_ephem;

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
	axo = 0.0; ayo = 0.0; azo = 0.0;	
    }

    // TODO: eliminate the output files after testing
    // or make this more flexible
    FILE *outfile = NULL;
    outfile = fopen("acc.out", "w");

    // These should be executed in order from smallest
    // to largest

    FILE *eih_file = NULL;
    eih_file = fopen("eih_acc.out", "w");
    
    eih_GR(sim, eih_loop_limit,
	   xo, yo, zo, vxo, vyo, vzo, axo, ayo, azo,	   
	   outfile, eih_file);
    
    earth_J2J4(sim, xo, yo, zo, outfile);

    direct(sim, xo, yo, zo, outfile);
    
    solar_J2(sim, xo, yo, zo, outfile);

    non_gravs(sim, xo, yo, zo, vxo, vyo, vzo, outfile);

    // Pick one or the other of the next two routines

    //simple_GR(sim, xo, yo, zo, vxo, vyo, vzo, outfile);



    FILE *vfile = NULL;
    static int first=1;
    if(first==1){
	vfile = fopen("vary_acc.out", "w");
    }
    
    //test_vary(sim, vfile);
    test_vary_2nd(sim, vfile);    

    if(first==1){
	fclose(vfile);
	first=1;
    }

    fclose(eih_file);
    fflush(outfile);
    fclose(outfile);

    if(geo == 1){
	// geocentric
	// TODO: This part will need work for the variational equations
	// to work properly.
	all_ephem(3, t, &GM, &xo, &yo, &zo, &vxo, &vyo, &vzo, &axo, &ayo, &azo);

	// This is the indirect term for geocentric equations
	// of motion.
	for (int j=0; j<N_real; j++){    

	    sim->particles[j].ax -= axo;
	    sim->particles[j].ay -= ayo;
	    sim->particles[j].az -= azo;

	}
    }
}

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
int integration_function(double tstart, double tend, double tstep,
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

    struct reb_simulation* sim = reb_create_simulation();

    sim->t = tstart;
    sim->dt = tstep;    // time step in days, this is just an initial value.

    // TODO: decide how flexible these should be.
    sim->integrator = REB_INTEGRATOR_IAS15;
    sim->save_messages = 1;
    sim->heartbeat = heartbeat;
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

    int step = sim->steps_done;

    outtime = ts->t;
    outstate = ts->state;
    int n_alloc= ts->n_alloc;

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

    //int nsubsteps = assist->nsubsteps;
    //double* hg = assist->hg;

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

	// Convenience variable.  The 'br' field contains the 
	// set of coefficients from the last completed step.
	const struct reb_dpconst7 b  = dpcast(sim->ri_ias15.br);

	double* x0 = malloc(sizeof(double)*N3);
	double* v0 = malloc(sizeof(double)*N3);
	double* a0 = malloc(sizeof(double)*N3);

	//double t = sim->t + sim->dt_last_done * (-1.0 + hg[n]);
	//printf("%lf %lf \n", sim->t-sim->dt_last_done, sim->dt_last_done);
	
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

	    //printf("\n%d\n", j);
	    //printf("%.16le %.16le %.16le\n", x0[k0], v0[k0], a0[k0]);
	    //printf("%.16le %.16le %.16le %.16le %.16le %.16le %.16le\n",
	    //b.p0[k0], b.p1[k0], b.p2[k0], b.p3[k0], b.p4[k0], b.p5[k0], b.p6[k0]);
	    //printf("%.16le %.16le %.16le\n", x0[k1], v0[k1], a0[k1]);
	    //printf("%.16le %.16le %.16le %.16le %.16le %.16le %.16le\n",
	    //b.p0[k1], b.p1[k1], b.p2[k1], b.p3[k1], b.p4[k1], b.p5[k1], b.p6[k1]);
	    //printf("%.16le %.16le %.16le\n", x0[k2], v0[k2], a0[k2]);	    
	    //printf("%.16le %.16le %.16le %.16le %.16le %.16le %.16le\n",
	    //b.p0[k2], b.p1[k2], b.p2[k2], b.p3[k2], b.p4[k2], b.p5[k2], b.p6[k2]);	    

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

void direct(struct reb_simulation* sim, double xo, double yo, double zo, FILE *outfile){

    //const double G = sim->G;
    const unsigned int N = sim->N;  // N includes real+variational particles
    const unsigned int N_real = N - sim->N_var;

    const double t = sim->t;    

    int N_ephem, N_ast;
    int number_bodies(int* N_ephem, int* N_ast);

    const int N_tot = number_bodies(&N_ephem, &N_ast);

    struct reb_particle* const particles = sim->particles;

    double GM;
    double x, y, z, vx, vy, vz, ax, ay, az;

    static const order[] = {26, 25, 24, 23, 22, 21, 20, 19, 18, 17,
			    16, 15, 14, 13, 12, 11, 10, 9, 8, 5, 4, 
			    1, 2, 7, 6, 3, 0};

    // Direct forces from massives bodies
    //for (int i=0; i<N_tot; i++){
    for (int k=0; k<N_tot; k++){    
	int i = order[k];
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

	    if(outfile){
		fprintf(outfile, "%3d %25.16le %25.16le %25.16le %25.16le %25.16le %25.16le %25.16le %25.16le\n", i, t, GM, dx, dy, dz, -prefac*dx, -prefac*dy, -prefac*dz);
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
	//for (int i=N_tot-1; i>=0; i--){
	//for (int i=0; i<N_tot; i++){

	int i = order[k];

        // Get position and mass of massive body i.	
	int flag = all_ephem(i, t, &GM, &x, &y, &z, &vx, &vy, &vz, &ax, &ay, &az);
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

void earth_J2J4(struct reb_simulation* sim, double xo, double yo, double zo, FILE *outfile){

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
    all_ephem(3, t, &GM, &xe, &ye, &ze, &vxe, &vye, &vze, &axe, &aye, &aze);
    const double GMearth = GM;

    double xr, yr, zr; //, vxr, vyr, vzr, axr, ayr, azr;
    xr = xe;  yr = ye;  zr = ze;

    // Hard-coded constants.  BEWARE!
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
	    fprintf(outfile, "%3s %25.16le %25.16le %25.16le %25.16le\n", "J24", t, resx, resy, resz);
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

void solar_J2(struct reb_simulation* sim, double xo, double yo, double zo, FILE *outfile){

    const unsigned int N = sim->N;  // N includes real+variational particles
    const unsigned int N_real = N - sim->N_var;

    const double t = sim->t;    

    struct reb_particle* const particles = sim->particles;

    double GM;

    // The Sun center is reference for these calculations.

    double xr, yr, zr, vxr, vyr, vzr, axr, ayr, azr;

    all_ephem(0, t, &GM, &xr, &yr, &zr, &vxr, &vyr, &vzr, &axr, &ayr, &azr);
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
	    fprintf(outfile, "%3s %25.16le %25.16le %25.16le %25.16le\n", "J2", t, resx, resy, resz);
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

void non_gravs(struct reb_simulation* sim,
	       double xo, double yo, double zo,
	       double vxo, double vyo, double vzo,	       
	       FILE *outfile){

    const unsigned int N = sim->N;  // N includes real+variational particles
    const unsigned int N_real = N - sim->N_var;
    const unsigned int N_var = sim->N_var;  // N includes real+variational particles    

    const double t = sim->t;    

    struct reb_particle* const particles = sim->particles;

    struct assist_extras* assist = (struct assist_extras*) sim->extras;

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
    all_ephem(0, t, &GMsun, &xs, &ys, &zs, &vxs, &vys, &vzs, &axs, &ays, &azs);

    xr = xs;  yr = ys;  zr = zs;
    vxr = vxs; vyr = vys; vzr = vzs;    

    // The non-grav parameters are specific to each object being
    // integrated.

    // Normal asteroids
    //double A1 = 0.0;
    //double A2 = 0.0;
    //double A3 = 0.0;

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

	//printf(" A123: %lf %lf %lf\n", A1, A2, A3);
	

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

	// 'Oumuamua
	double ALN = 0.04083733261;
	double NK = 2.6;
	double NM = 2.0;
	double NN = 3.0;
	double r0 = 5.0;
	
	const double g = ALN*pow(r/r0, -NM)*pow(1.0+pow(r/r0, NN), -NK);

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

	particles[j].ax += A1*g*dx/r + A2*g*tx/_t + A3*g*hx/h;
        particles[j].ay += A1*g*dy/r + A2*g*ty/_t + A3*g*hy/h;
        particles[j].az += A1*g*dz/r + A2*g*tz/_t + A3*g*hz/h;

	// variational matrix elements
	// Only evaluate the constants if there are variational particles

        const double r3    = r*r*r;
        const double v2    = dvx*dvx + dvy*dvy + dvz*dvz;
        const double rdotv = dx*dvx  + dy*dvy  + dz*dvz;
        const double vdott = dvx*tx  + dvy*ty  + dvz*tz;

	// Need to update this for the new g(r) function.
	const double dgdr = -2.*g/r;
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
		//double dA1 = part_params[N_real+v].A1;
		//double dA2 = part_params[N_real+v].A2;
		//double dA3 = part_params[N_real+v].A3;
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

void simple_GR(struct reb_simulation* sim,
	       double xo, double yo, double zo,
	       double vxo, double vyo, double vzo,	       
	       FILE *outfile){

    // Damour and Deruelle solar GR treatment

    const double au = JPL_EPHEM_CAU;    
    const double c = (JPL_EPHEM_CLIGHT/au)*86400;
    const double C2 = c*c;  // This could be stored as C2.
    
    const unsigned int N = sim->N;  // N includes real+variational particles
    const unsigned int N_real = N - sim->N_var;

    const double t = sim->t;    

    struct reb_particle* const particles = sim->particles;

    double GM;

    double xr, yr, zr, vxr, vyr, vzr, axr, ayr, azr;
    
    // The Sun center is reference for these calculations.
    double xs, ys, zs, vxs, vys, vzs, axs, ays, azs;

    all_ephem(0, t, &GM, &xs, &ys, &zs, &vxs, &vys, &vzs, &axs, &ays, &azs);
    const double GMsun = GM;    

    xr  = xs;  yr  = ys;  zr = zs;
    vxr = vxs; vyr = vys; vzr = vzs;
    axr = axs; ayr = ays; azr = azs;    

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
	    fprintf(outfile, "%3s %25.16le %25.16le %25.16le %25.16le\n", "GR", t,
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

void eih_GR(struct reb_simulation* sim,
	    int eih_loop_limit,
	    double xo, double yo, double zo,
	    double vxo, double vyo, double vzo,
	    double axo, double ayo, double azo,	       	    
	    FILE *outfile,
	    FILE *eih_file){

    // Einstein-Infeld-Hoffman PPN GR treatment
    // This is one of two options for GR.
    // This one version is only rarely needed.

    const double au = JPL_EPHEM_CAU;    
    const double c = (JPL_EPHEM_CLIGHT/au)*86400;
    const double C2 = c*c;  // This could be stored as C2.
    
    // Doesn't need to be hard-coded.
    //const double c = 173.14463267424031;
    //const double C2 = c*c;  // This could be stored as C2.
    
    const unsigned int N = sim->N;  // N includes real+variational particles
    const unsigned int N_real = N - sim->N_var;

    const double t = sim->t;    

    struct reb_particle* const particles = sim->particles;

    double GM;
    //double x, y, z, vx, vy, vz, ax, ay, az;

    double xr, yr, zr, vxr, vyr, vzr, axr, ayr, azr;
    
    // The Sun center is reference for these calculations.
    double xs, ys, zs, vxs, vys, vzs, axs, ays, azs;
    all_ephem(0, t, &GM, &xs, &ys, &zs, &vxs, &vys, &vzs, &axs, &ays, &azs);
    const double GMsun = GM;    

    xr  = xs;  yr  = ys;  zr = zs;
    vxr = vxs; vyr = vys; vzr = vzs;
    axr = axs; ayr = ays; azr = azs;

    int N_ephem, N_ast;
    int number_bodies(int* N_ephem, int* N_ast);

    const int N_tot = number_bodies(&N_ephem, &N_ast);

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

	for (int j=0; j<eih_loop_limit; j++){ // This is either 1 or N_ephem
	//for (int j=0; j<1; j++){	
	//for (int j=0; j<N_ephem; j++){

	    // Get position and mass of massive body j.
	    all_ephem(j, t, &GMj,
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

	    const double term2 = gamma/C2*vi2;
	    const double dterm2dvx = 2.0*gamma/C2*particles[i].vx;
	    const double dterm2dvy = 2.0*gamma/C2*particles[i].vy;
	    const double dterm2dvz = 2.0*gamma/C2*particles[i].vz;	    

	    const double vj2 = (vxj-vxo)*(vxj-vxo) + (vyj-vyo)*(vyj-vyo) + (vzj-vzo)*(vzj-vzo);

	    const double term3 = (1+gamma)/C2*vj2;
	    // Variational equations do not depend on term3

	    const double vidotvj = particles[i].vx*(vxj-vxo) +
		particles[i].vy*(vyj-vyo) +
		particles[i].vz*(vzj-vzo);

	    const double term4 = -2*(1+gamma)/C2*vidotvj;
	    const double dterm4dvx = -2*(1+gamma)/C2*(vxj-vxo);
	    const double dterm4dvy = -2*(1+gamma)/C2*(vyj-vyo);
	    const double dterm4dvz = -2*(1+gamma)/C2*(vzj-vzo);	    	    
	    

	    const double rijdotvj = dxij*(vxj-vxo) + dyij*(vyj-vyo) + dzij*(vzj-vzo);

	    if(eih_file){
		fprintf(eih_file, " EIH_J%12d\n", j);	    
		fprintf(eih_file, "%25.16lE ", rijdotvj/_rij);
	    }

	    const double term5 = -1.5/C2*(rijdotvj*rijdotvj)/(_rij*_rij);
	    const double dterm5dx = -3.0/C2*rijdotvj/_rij*((vxj-vxo)/_rij - rijdotvj*dxij/(_rij*_rij*_rij));
	    const double dterm5dy = -3.0/C2*rijdotvj/_rij*((vyj-vyo)/_rij - rijdotvj*dyij/(_rij*_rij*_rij));
	    const double dterm5dz = -3.0/C2*rijdotvj/_rij*((vzj-vzo)/_rij - rijdotvj*dzij/(_rij*_rij*_rij));	    	    

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
	    
	    for (int k=0; k<N_ephem; k++){

		// Get position and mass of massive body k.
		all_ephem(k, t, &GMk,
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

	    term0 *= -2*(beta+gamma)/C2;
	    dterm0dx *= -2*(beta+gamma)/C2;
	    dterm0dy *= -2*(beta+gamma)/C2;
	    dterm0dz *= -2*(beta+gamma)/C2;	    	    
	    
	    term1 *= -(2*beta-1)/C2;

	    const double rijdotaj = dxij*(axj-axo) + dyij*(ayj-ayo) + dzij*(azj-azo);
	    const double term6 = -0.5/C2*rijdotaj;
	    const double dterm6dx = -0.5/C2*(axj-axo);
	    const double dterm6dy = -0.5/C2*(ayj-ayo);	    
	    const double dterm6dz = -0.5/C2*(azj-azo);
	    
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

	grx += term7x_sum/C2 + term8x_sum/C2;
	gry += term7y_sum/C2 + term8y_sum/C2;
	grz += term7z_sum/C2 + term8z_sum/C2;

	if(outfile){
	    fprintf(outfile, "%3s %25.16le %25.16le %25.16le %25.16le\n", "GR", t,
		    grx, gry, grz);
	    fflush(outfile);
	}

	dxdx += dterm7x_sumdx/C2 + dterm8x_sumdx/C2;
	dxdy += dterm7x_sumdy/C2 + dterm8x_sumdy/C2;
	dxdz += dterm7x_sumdz/C2 + dterm8x_sumdz/C2;	
	dxdvx += dterm7x_sumdvx/C2;
	dxdvy += dterm7x_sumdvy/C2;
	dxdvz += dterm7x_sumdvz/C2;

	dydx += dterm7y_sumdx/C2 + dterm8y_sumdx/C2;
	dydy += dterm7y_sumdy/C2 + dterm8y_sumdy/C2;
	dydz += dterm7y_sumdz/C2 + dterm8y_sumdz/C2;	
	dydvx += dterm7y_sumdvx/C2;
	dydvy += dterm7y_sumdvy/C2;
	dydvz += dterm7y_sumdvz/C2;

	dzdx += dterm7z_sumdx/C2 + dterm8z_sumdx/C2;
	dzdy += dterm7z_sumdy/C2 + dterm8z_sumdy/C2;
	dzdz += dterm7z_sumdz/C2 + dterm8z_sumdz/C2;	
	dzdvx += dterm7z_sumdvx/C2;
	dzdvy += dterm7z_sumdvy/C2;
	dzdvz += dterm7z_sumdvz/C2;
	
	particles[i].ax += term7x_sum/C2 + term8x_sum/C2;
	particles[i].ay += term7y_sum/C2 + term8y_sum/C2;
	particles[i].az += term7z_sum/C2 + term8z_sum/C2;

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

void test_vary(struct reb_simulation* sim, FILE *vfile){

    static int first=1;
    
    const unsigned int N = sim->N;  // N includes real+variational particles
    const unsigned int N_real = N - sim->N_var;

    const double t = sim->t;    

    struct reb_particle* const particles = sim->particles;

    // What is happening is that a set of test particles with
    // small initial displacements from the 0th test particles
    // are integrated with the 0th test particle.  The differences
    // between the trajectories can be compared to results from
    // the variational particles.
    
    double delt = 1e-8;
    if(first==1){
	for (int j=1; j<N_real; j++){
	    double dx = particles[j].ax - particles[0].ax;
	    double dy = particles[j].ay - particles[0].ay;
	    double dz = particles[j].az - particles[0].az;
	    if(vfile){
		fprintf(vfile, "%3d %25.16le %25.16le %25.16le %25.16le\n", j, t, dx/delt, dy/delt, dz/delt);
	    }
	}
    
	for (int j=0; j<N_real; j++){ //loop over test particles
	    for (int v=0; v < sim->var_config_N; v++){
		struct reb_variational_configuration const vc = sim->var_config[v];
		int tp = vc.testparticle;
		struct reb_particle* const particles_var1 = particles + vc.index;
		if(vfile){		
		    if(tp == j){
			fprintf(vfile, "%3d %25.16le %25.16le %25.16le %25.16le\n",
				j, t, particles_var1[0].ax, particles_var1[0].ay, particles_var1[0].az);
		    }
		}
	    }
	}
    }
}

void test_vary_2nd(struct reb_simulation* sim, FILE *vfile){

    static int first=1;
    
    const unsigned int N = sim->N;  // N includes real+variational particles
    const unsigned int N_real = N - sim->N_var;

    const double t = sim->t;    

    struct reb_particle* const particles = sim->particles;

    // What is happening is that a set of test particles with
    // small initial displacements from the 0th test particles
    // are integrated with the 0th test particle.
    // This assumes that there is a pair of test particles for
    // each of the six dimensions of the initial conditions.
    // The differences between the trajectories can be used to
    // get 2nd order numerical derivatives to compared to result
    // from the variational particles.
    
    double delt = 1e-8;
    if(first==1){
	for (int j=1; j<7; j++){
	    double dx = particles[j+6].ax - particles[j].ax;
	    double dy = particles[j+6].ay - particles[j].ay;
	    double dz = particles[j+6].az - particles[j].az;
	    if(vfile){
		fprintf(vfile, "%3d %25.16le %25.16le %25.16le %25.16le\n", j, t, 0.5*dx/delt, 0.5*dy/delt, 0.5*dz/delt);
	    }
	}
    
	for (int j=0; j<N_real; j++){ //loop over test particles
	    for (int v=0; v < sim->var_config_N; v++){
		struct reb_variational_configuration const vc = sim->var_config[v];
		int tp = vc.testparticle;
		struct reb_particle* const particles_var1 = particles + vc.index;
		if(vfile){		
		    if(tp == j){
			fprintf(vfile, "%3d %25.16le %25.16le %25.16le %25.16le\n",
				j, t, particles_var1[0].ax, particles_var1[0].ay, particles_var1[0].az);
		    }
		}
	    }
	}
    }
}

