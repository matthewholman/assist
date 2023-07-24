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

#define STRINGIFY(s) str(s)
#define str(s) #s

const char* assist_build_str = __DATE__ " " __TIME__;   // Date and time build string. 
const char* assist_version_str = "1.1.3";         // **VERSIONLINE** This line gets updated automatically. Do not edit manually.
const char* assist_githash_str = STRINGIFY(ASSISTGITHASH);// This line gets updated automatically. Do not edit manually.


// These correspond to ASSIST_STATUS enum.

const char* assist_error_messages[] = {
    "No error has occured.", // ASSIST_SUCCESS
    "The JPL planet ephemeris file has not been found.", // ASSIST_ERROR_EPHEM_FILE
    "The JPL asteroid ephemeris file has not been found. Asteroid forces have been disabled.", // ASSIST_ERROR_AST_FILE
    "The requested asteroid ID has not been found.", // ASSIST_ERROR_NAST
    "The requested planet ID has not been found.", // ASSIST_ERROR_NEPHEM
    "The requested time is outside the coverage provided by the ephemeris file.", // ASSIST_ERROR_COVERAGE
};
const int assist_error_messages_N = ASSIST_ERROR_N;

    
// Forward function declarations
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

int assist_ephem_init(struct assist_ephem* ephem, char *user_planets_path, char *user_asteroids_path){

    char default_planets_path[] = "/data/linux_m13000p17000.441";
    char default_asteroids_path[] = "/data/sb441-n16.bsp";

    ephem->jd_ref = 2451545.0; // Default jd_ref
    
    const int FNAMESIZE = 1024;
    char planets_path[FNAMESIZE];
    char asteroids_path[FNAMESIZE];        

    /** Use user-defined file or the default filename, 
     *  in that order.
     */

    if(user_planets_path == NULL && getenv("ASSIST_DIR")==NULL){
        return ASSIST_ERROR_EPHEM_FILE;	  
    }

    if(user_planets_path == NULL){
        sprintf(planets_path, "%s%s", getenv("ASSIST_DIR"), default_planets_path);
    }else{
        strncpy(planets_path, user_planets_path, FNAMESIZE-1);	
    }

    if ((ephem->jpl = assist_jpl_init(planets_path)) == NULL) {
        return ASSIST_ERROR_EPHEM_FILE;	  
    }

    int asteroids_path_not_found = 0;


    if(user_asteroids_path == NULL){
        if(getenv("ASSIST_DIR")==NULL){
            asteroids_path_not_found = 1;
        }else{
            sprintf(asteroids_path, "%s%s", getenv("ASSIST_DIR"), default_asteroids_path);
        }
    }else{
        strncpy(asteroids_path, user_asteroids_path, FNAMESIZE-1);	
    }

    if (asteroids_path_not_found != 1){
        if ((ephem->spl = assist_spk_init(asteroids_path)) == NULL) {
            asteroids_path_not_found = 1;
        }
    }
            
    if (asteroids_path_not_found != 1){
        // Try to find masses of bodies in spk file in ephemeris constants
        for(int n=0; n<ephem->spl->num; n++){ // loop over all asteroids
            int found = 0;
            for(int c=0; c<ephem->jpl->num; c++){ // loop over all constants
                if (strncmp(ephem->jpl->str[c], "MA", 2) == 0) {
                    int cid = atoi(ephem->jpl->str[c]+2);
                    int offset = 2000000;
                    if (cid==ephem->spl->targets[n].code-offset){
                        ephem->spl->targets[n].mass = ephem->jpl->con[c];
                        found = 1;
                        break;
                    }
                }
            }
            // Use lookup table for new KBO objects in DE440/441
            // Source: https://ssd.jpl.nasa.gov/ftp/eph/planets/bsp/README.txt
            int massmap[] = {
                // ID, SPK_ID
                8001,  2136199,
                8002,  2136108,
                8003,  2090377,
                8004,  2136472,
                8005,  2050000,
                8006,  2084522,
                8007,  2090482,
                8008,  2020000,
                8009,  2055637,
                8010,  2028978,
                8011,  2307261,
                8012,  2174567,
                8013,  3361580,
                8014,  3308265,
                8015,  2055565,
                8016,  2145452,
                8017,  2090568,
                8018,  2208996,
                8019,  2225088,
                8020,  2019521,
                8021,  2120347,
                8022,  2278361,
                8023,  3525142,
                8024,  2230965,
                8025,  2042301,
                8026,  2455502,
                8027,  3545742,
                8028,  2523639,
                8029,  2528381,
                8030,  3515022,
            };
            if (found==0){
                int mapped = -1;
                for (int m=0; m<sizeof(massmap); m+=2){
                    if (massmap[m+1]==ephem->spl->targets[n].code){
                        mapped = massmap[m];
                        break;
                    }
                }
                if (mapped != -1){
                    for(int c=0; c<ephem->jpl->num; c++){ // loop over all constants (again)
                        if (strncmp(ephem->jpl->str[c], "MA", 2) == 0) {
                            int cid = atoi(ephem->jpl->str[c]+2);
                            if (cid==mapped){
                                ephem->spl->targets[n].mass = ephem->jpl->con[c];
                                found = 1;
                                break;
                            }
                        }
                    }
                }
            }
            if (found==0){
                fprintf(stderr,"WARNING: Cannot find mass for asteroid %d (NAIF ID Number %d).\n", n, ephem->spl->targets[n].code );
            }

        }
    }else{
        fprintf(stderr, "(ASSIST) %s\n", assist_error_messages[ASSIST_ERROR_AST_FILE]);
    }

    return ASSIST_SUCCESS;
}

struct assist_ephem* assist_ephem_create(char *user_planets_path, char *user_asteroids_path){
    struct assist_ephem* ephem = calloc(1, sizeof(struct assist_ephem));
    int error = assist_ephem_init(ephem, user_planets_path, user_asteroids_path);
    if (error != ASSIST_SUCCESS){
        fprintf(stderr, "(ASSIST) An error occured while trying to initialize the ephemeris structure.\n");
        fprintf(stderr, "(ASSIST) %s\n", assist_error_messages[error]);
        assist_ephem_free(ephem);
        return NULL;
    }
    return ephem;
}

void assist_ephem_free_pointers(struct assist_ephem* ephem){
    if (ephem->jpl){
        assist_jpl_free(ephem->jpl);
    }
    if (ephem->spl){
        assist_spk_free(ephem->spl);
    }
}
void assist_ephem_free(struct assist_ephem* ephem){
    assist_ephem_free_pointers(ephem);
    free(ephem);
}

struct assist_extras* assist_attach(struct reb_simulation* sim, struct assist_ephem* ephem){  
    if (sim == NULL){
        fprintf(stderr, "(ASSIST) Error: Simulation pointer passed to assist_attach was NULL.\n");
        return NULL;
    }
    int extras_should_free_ephem = 0;
    if (ephem == NULL){
        // Try default 
        ephem = assist_ephem_create(NULL, NULL);
        if (ephem == NULL){
            fprintf(stderr, "(ASSIST) Error: Ephemeris pointer passed to assist_attach was NULL. Initialization with default path failed.\n");
            return NULL;
        }
        extras_should_free_ephem = 1;
    }

    // Initialization separate from memory allocation because python handles memory management
    struct assist_extras* assist = calloc(1, sizeof(*assist));
    assist_init(assist, sim, ephem); 
    assist->extras_should_free_ephem = extras_should_free_ephem;
    
    return assist;
}

void assist_extras_cleanup(struct reb_simulation* sim){
    struct assist_extras* assist = sim->extras;
    assist->sim = NULL;
}

void assist_init(struct assist_extras* assist, struct reb_simulation* sim, struct assist_ephem* ephem){
    assist->sim = sim;
    assist->ephem_cache = calloc(1, sizeof(struct assist_ephem_cache));
    int N_total = ASSIST_BODY_NPLANETS;
    if (ephem->spl){
        N_total += ephem->spl->num;
    }
    assist->gr_eih_sources = 1; // Only include Sun by default
    assist->ephem_cache->items = calloc(N_total*7, sizeof(struct assist_cache_item));
    assist->ephem_cache->t = malloc(N_total*7*sizeof(double));
    for (int i=0;i<7*N_total;i++){
        assist->ephem_cache->t[i] = -1e306;
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
    assist->last_state = NULL; 
    assist->current_state = NULL; 
    sim->integrator = REB_INTEGRATOR_IAS15;
    sim->gravity = REB_GRAVITY_NONE;
    sim->extras = assist;
    sim->extras_cleanup = assist_extras_cleanup;
    sim->additional_forces = assist_additional_forces;
    sim->force_is_velocity_dependent = 1;
}

void assist_free_pointers(struct assist_extras* assist){
    if (assist->sim){
        assist_detach(assist->sim, assist);
        assist->sim = NULL;
    }
    if (assist->last_state){
        free(assist->last_state);
        assist->last_state = NULL;
    }
    if (assist->current_state){
        free(assist->current_state);
        assist->current_state = NULL;
    }
    if (assist->ephem_cache){
        if (assist->ephem_cache->items){
            free(assist->ephem_cache->items);
        }
        if (assist->ephem_cache->t){
            free(assist->ephem_cache->t);
        }
        free(assist->ephem_cache);
        assist->ephem_cache = NULL;
    }
    if (assist->extras_should_free_ephem && assist->ephem){
        assist_ephem_free(assist->ephem);
        assist->ephem = NULL;
    }
}


void assist_free(struct assist_extras* assist){
    // Freeing pointers is separate because python handles memory management of structure itself.
    assist_free_pointers(assist);
    free(assist);
}

void assist_detach(struct reb_simulation* sim, struct assist_extras* assist){
    if (assist->sim){
        sim->extras = NULL;
        sim->extras_cleanup = NULL;
        sim->additional_forces = NULL;
        sim->pre_timestep_modifications = NULL;
    }
    assist->sim = NULL;
}

void assist_error(struct assist_extras* assist, const char* const msg){
    if (assist->sim == NULL){
        fprintf(stderr, "(ASSIST) Error: A Simulation is no longer attached to the ASSIST extras instance. Most likely the Simulation has been freed.\n");
    } else{
        reb_error(assist->sim, msg);
    }
}


struct reb_particle assist_get_particle_with_error(struct assist_ephem* ephem, const int particle_id, const double t, int* error){
    struct reb_particle p = {0};
    double GM = 0;
    int flag = assist_all_ephem(ephem, NULL, particle_id, t, &GM, &p.x, &p.y, &p.z, &p.vx, &p.vy, &p.vz, &p.ax, &p.ay, &p.az);
    *error = flag;
    p.m = GM; // Note this is GM, not M
    return p;
}


struct reb_particle assist_get_particle(struct assist_ephem* ephem, const int particle_id, const double t){
    int error = 0;
    struct reb_particle p = assist_get_particle_with_error(ephem, particle_id, t, &error);
    if (error != ASSIST_SUCCESS){
        fprintf(stderr, "(ASSIST) An error occured while trying to initialize particle from ephemeris data.\n");
        fprintf(stderr, "(ASSIST) %s\n", assist_error_messages[error]);
    }
    return p;
}

void assist_interpolate(const struct reb_particle* const last_state, const struct reb_dp7 b_coeff, double dt_last_done, double h, int N, struct reb_particle* output){
    const struct reb_dpconst7 b  = dpcast(b_coeff);

    double s[9]; // Summation coefficients
    double sv[9]; // Summation coefficients

    s[0] = dt_last_done * h;

    s[1] = s[0] * s[0] / 2.;
    s[2] = s[1] * h / 3.;
    s[3] = s[2] * h / 2.;
    s[4] = 3. * s[3] * h / 5.;
    s[5] = 2. * s[4] * h / 3.;
    s[6] = 5. * s[5] * h / 7.;
    s[7] = 3. * s[6] * h / 4.;
    s[8] = 7. * s[7] * h / 9.;
    
    sv[0] = dt_last_done * h;
    sv[1] =      sv[0] * h / 2.;
    sv[2] = 2. * sv[1] * h / 3.;
    sv[3] = 3. * sv[2] * h / 4.;
    sv[4] = 4. * sv[3] * h / 5.;
    sv[5] = 5. * sv[4] * h / 6.;
    sv[6] = 6. * sv[5] * h / 7.;
    sv[7] = 7. * sv[6] * h / 8.;

    // Predict positions and velocities at interval n using b values
    // for all the particles
    for(int j=0;j<N;j++) {
        const int k0 = 3*j+0;
        const int k1 = 3*j+1;
        const int k2 = 3*j+2;

        output[j].x = last_state[j].x + (s[8]*b.p6[k0] + s[7]*b.p5[k0] + s[6]*b.p4[k0] + s[5]*b.p3[k0] + s[4]*b.p2[k0] + s[3]*b.p1[k0] + s[2]*b.p0[k0] + s[1]*last_state[j].ax + s[0]*last_state[j].vx );
        output[j].y = last_state[j].y + (s[8]*b.p6[k1] + s[7]*b.p5[k1] + s[6]*b.p4[k1] + s[5]*b.p3[k1] + s[4]*b.p2[k1] + s[3]*b.p1[k1] + s[2]*b.p0[k1] + s[1]*last_state[j].ay + s[0]*last_state[j].vy );
        output[j].z = last_state[j].z + (s[8]*b.p6[k2] + s[7]*b.p5[k2] + s[6]*b.p4[k2] + s[5]*b.p3[k2] + s[4]*b.p2[k2] + s[3]*b.p1[k2] + s[2]*b.p0[k2] + s[1]*last_state[j].az + s[0]*last_state[j].vz );

        output[j].vx = last_state[j].vx + sv[7]*b.p6[k0] + sv[6]*b.p5[k0] + sv[5]*b.p4[k0] + sv[4]*b.p3[k0] + sv[3]*b.p2[k0] + sv[2]*b.p1[k0] + sv[1]*b.p0[k0] + sv[0]*last_state[j].ax;
        output[j].vy = last_state[j].vy + sv[7]*b.p6[k1] + sv[6]*b.p5[k1] + sv[5]*b.p4[k1] + sv[4]*b.p3[k1] + sv[3]*b.p2[k1] + sv[2]*b.p1[k1] + sv[1]*b.p0[k1] + sv[0]*last_state[j].ay;
        output[j].vz = last_state[j].vz + sv[7]*b.p6[k2] + sv[6]*b.p5[k2] + sv[5]*b.p4[k2] + sv[4]*b.p3[k2] + sv[3]*b.p2[k2] + sv[2]*b.p1[k2] + sv[1]*b.p0[k2] + sv[0]*last_state[j].az;
    }
}

struct reb_simulation* assist_create_interpolated_simulation(struct reb_simulationarchive* sa, double t){
    if (sa==NULL) return NULL;

    // Find blob just after time t
    if (t <= sa->t[1]){ // Note cannot use first blob because accelerations are missing!
        printf("Requested time outside range of SimulationArchive.\n");
        return NULL;
    }
    if (t >= sa->t[sa->nblobs-1]){
        printf("Requested time outside range of SimulationArchive.\n");
        return NULL;
    }
    long blob = 0;
    for (long i=1; i<sa->nblobs; i++){
        if (sa->t[i] >= t){
            blob = i;
            break;
        }
    }
    
    //Very hackish solutions. Should be improved!
    enum reb_input_binary_messages warnings = REB_INPUT_BINARY_WARNING_NONE;
    struct reb_simulation* r2 = reb_create_simulation();
    reb_create_simulation_from_simulationarchive_with_messages(r2, sa, blob-1, &warnings);
    // r2 = reb_input_process_warnings(r2, warnings); Ignoring warnings for now
    
    struct reb_simulation* r3 = reb_create_simulation();
    reb_create_simulation_from_simulationarchive_with_messages(r3, sa, blob, &warnings);
    // r3 = reb_input_process_warnings(r3, warnings);
    
    double h = (t - r2->t)/(r3->dt_last_done);
    assist_interpolate_simulation(r2, r3, h);
    reb_free_simulation(r3);
    return r2;
}

void assist_swap_particles(struct reb_simulation* sim){
    struct assist_extras* ax = sim->extras;
    struct reb_particle* p = sim->particles;
    sim->particles = ax->current_state;
    ax->current_state = p; 
}

void assist_integrate_or_interpolate(struct assist_extras* ax, double t){
    struct reb_simulation* sim = ax->sim;
    
    sim->pre_timestep_modifications = assist_pre_timestep_modifications;
    sim->exact_finish_time = 0;

    if (ax->current_state==NULL){
        ax->current_state = malloc(sizeof(struct reb_particle)*6*sim->N);
        ax->last_state = malloc(sizeof(struct reb_particle)*6*sim->N);
    }else{
        assist_swap_particles(sim);
    }

    double dts = copysign(1., sim->dt_last_done);
    if ( !(dts*(sim->t-sim->dt_last_done)  <  dts*t &&  dts*t < dts*sim->t) ){
        // Integrate if requested time not in interval of last timestep
        reb_integrate(sim, t);
    }
    double h = 1.0-(sim->t -t) / sim->dt_last_done; 
    if (sim->t - t==0.){
        memcpy(ax->current_state, sim->particles, sizeof(struct reb_particle)*sim->N);
    }else if (h<0.0 || h>=1.0 || !isnormal(h)){
        printf("Error: cannot interpolate beyond timestep bounds (h=%e).\n",h);
    }else if (sim->ri_ias15.br.p0 == NULL) {
        printf("Error: cannot interpolate before first timestep is complete (h=%e).\n",h);
    }else{
        assist_interpolate(ax->last_state, sim->ri_ias15.br, sim->dt_last_done, h, sim->N, ax->current_state);
    }
    assist_swap_particles(sim);
}

int assist_interpolate_simulation(struct reb_simulation* sim1, struct reb_simulation* sim2, double h){
    int N = sim1->N;

    // Convenience variable.  The 'br' field contains the
    // set of coefficients from the last completed step.
    const struct reb_dpconst7 b  = dpcast(sim2->ri_ias15.br);

    double* x0 = sim1->ri_ias15.x0;
    double* v0 = sim1->ri_ias15.v0;
    double* a0 = sim2->ri_ias15.a0; // Note: sim2 !

    double s[9]; // Summation coefficients

    s[0] = sim2->dt_last_done * h;

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
        sim1->particles[j].x = xx0;
        sim1->particles[j].y = xy0;
        sim1->particles[j].z = xz0;
    }

    s[0] = sim2->dt_last_done * h;
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
        sim1->particles[j].vx = vx0;
        sim1->particles[j].vy = vy0;
        sim1->particles[j].vz = vz0;

    }
    sim1->t += s[0];
    return 1;
}

static void assist_pre_timestep_modifications(struct reb_simulation* sim){
    struct assist_extras* assist = sim->extras;
    reb_update_acceleration(sim); // This will later be recalculated. Could be optimized.
    memcpy(assist->last_state, sim->particles, sizeof(struct reb_particle)*sim->N);
}

