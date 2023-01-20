/**
 * @file    assist.h
 * @brief   ASSIST API definition.
 * @author  Hanno Rein 
 * 
 * @section     LICENSE
 * Copyright (c) 2022 Matthew Holman, Hanno Rein
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
#ifndef _ASSIST_ASSIST_H_H
#define _ASSIST_ASSIST_H_H

#ifndef M_PI 
#define M_PI 3.1415926535879323846
#endif

#include <stdint.h>
#include <limits.h>
#include "rebound.h"
#ifndef ASSISTGITHASH
#define ASSISTGITHASH notavailable0000000000000000000000000001 
#endif // ASSISTGITHASH

extern const char* assist_build_str;      ///< Date and time build string.
extern const char* assist_version_str;    ///< Version string.
extern const char* assist_githash_str;    ///< Current git hash.

typedef struct {
    double A1;
    double A2;
    double A3;        
} particle_params;


// ENUM to enable/disable different forces
enum ASSIST_FORCES { 
    ASSIST_FORCE_NONE               = 0,
    ASSIST_FORCE_SUN                = 0x01,
    ASSIST_FORCE_PLANETS            = 0x02,
    ASSIST_FORCE_ASTEROIDS          = 0x04,
    ASSIST_FORCE_NON_GRAVITATIONAL  = 0x08, 
    ASSIST_FORCE_EARTH_HARMONICS    = 0x10,
    ASSIST_FORCE_SUN_HARMONICS      = 0x20,
    ASSIST_FORCE_GR_EIH             = 0x40,
    ASSIST_FORCE_GR_SIMPLE          = 0x80,
    ASSIST_FORCE_GR_POTENTIAL       = 0x100,
};

struct assist_extras {
    struct reb_simulation* sim;
    int geocentric;
    double last_state_t;
    double* last_state_x;
    double* last_state_v;
    double* last_state_a;
    int nsubsteps;
    double* hg;
    //particle_params* particle_params;
    double* particle_params;
    double* output_t; 
    double* output_state;
    int output_n_alloc;
    int steps_done;
    int forces;
    double jd_ref;
};

/**
 * @brief Adds ASSIST functionality to a passed REBOUND simulation.
 * @param sim Pointer to the reb_simulation on which to add ASSIST functionality.
 * @return Pointer to an assist_extras structure.
 */
struct assist_extras* assist_attach(struct reb_simulation* sim);

/**
 * @brief Frees all memory allocated by ASSIST instance.
 * @details Should be called after simulation is done if memory is a concern.
 * @param assist The assist_extras pointer returned from the initial call to assist_attach.
 */
void assist_free(struct assist_extras* assist);

/**
 * @brief Detaches ASSIST from simulation, resetting all the simulation's function pointers that ASSIST has set.
 * @details This does not free the memory allocated by ASSIST (call assist_free).
 * @param sim Pointer to the simulation from which to remove ASSIST
 */
void assist_detach(struct reb_simulation* sim, struct assist_extras* assist);

/**
 * @brief Output an error message.
 * @details This function should be used if an error occurs rather than simply using print. 
 *          The message will be passed to python.
 * @param assist The assist_extras pointer.
 * @param msg The error message.
 */
void assist_error(struct assist_extras* assist, const char* const msg);


// Find particle position and velocity based on ephemeris data
struct reb_particle assist_get_particle(struct assist_extras* assist, const int particle_id, const double t);

// Functions called from python:
void assist_initialize(struct reb_simulation* sim, struct assist_extras* assist); // Initializes all pointers and values.
void assist_free_pointers(struct assist_extras* assist);

void assist_interpolate(struct reb_simulation* sim, double h, double* output);

int assist_integrate(double jd_ref,
		     double tstart, double tend, double tstep,
		     int geocentric,
		     double epsilon,
		     int n_particles,
		     double* instate,
		     //particle_params* part_params,
		     double* part_params,
		     int n_var,
		     int* invar_part,			 
		     double* invar,
		     //particle_params* var_part_params,
		     double* var_part_params,			 
		     int n_alloc,			 
		     int *n_out,
		     int nsubsteps,
		     double* hg,
		     double* outtime,
		     double* outstate,
		     double min_dt);


void test_vary(struct reb_simulation* sim, FILE *vfile);

void test_vary_2nd(struct reb_simulation* sim, FILE *vfile);



#endif
