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

enum ASSIST_STATUS{
    ASSIST_SUCCESS,         // no error
    ASSIST_ERROR_EPHEM_FILE,   // JPL ephemeris file not found
    ASSIST_ERROR_AST_FILE,     // JPL asteroid file not found
    ASSIST_ERROR_NAST,         // asteroid number out of range
    ASSIST_ERROR_NEPHEM,      // planet number out of range
};

enum ASSIST_BODY {
    ASSIST_BODY_SUN         = 0,
    ASSIST_BODY_MERCURY     = 1,
    ASSIST_BODY_VENUS       = 2,
    ASSIST_BODY_EARTH       = 3,
    ASSIST_BODY_MOON        = 4,
    ASSIST_BODY_MARS        = 5,
    ASSIST_BODY_JUPITER     = 6,
    ASSIST_BODY_SATURN      = 7,
    ASSIST_BODY_URANUS      = 8,
    ASSIST_BODY_NEPTUNE     = 9,
    ASSIST_BODY_PLUTO       = 10,

    ASSIST_BODY_NPLANETS    = 11,

    ASSIST_BODY_CAMILLA     = 11,
    ASSIST_BODY_CERES       = 12,
    ASSIST_BODY_CYBELE      = 13,
    ASSIST_BODY_DAVIDA      = 14,
    ASSIST_BODY_EUNOMIA     = 15,
    ASSIST_BODY_EUPHROSYNE  = 16,
    ASSIST_BODY_EUROPA      = 17,
    ASSIST_BODY_HYGIEA      = 18,
    ASSIST_BODY_INTERAMNIA  = 19,
    ASSIST_BODY_IRIS        = 20,
    ASSIST_BODY_JUNO        = 21,
    ASSIST_BODY_PALLAS      = 22,
    ASSIST_BODY_PSYCHE      = 23,
    ASSIST_BODY_SYLVIA      = 24,
    ASSIST_BODY_THISBE      = 25,
    ASSIST_BODY_VESTA       = 26,

    ASSIST_BODY_NASTEROIDS  = 16,
};

struct assist_ephem {
    double jd_ref;
    struct jpl_s* jpl;
    struct spk_s* spl;
};

struct assist_cache_item {
    double GM;
    double x;
    double y;
    double z;
    double vx;
    double vy;
    double vz;
    double ax;
    double ay;
    double az;
};


struct assist_ephem_cache {
    double* t;
    int* index;
    struct assist_cache_item* items;
};

struct assist_extras {
    struct reb_simulation* sim;
    struct assist_ephem* ephem;
    struct assist_ephem_cache* ephem_cache;
    int extras_should_free_ephem;   // Internal use only. Set to 1 if extras allocated memory for ephem.
    int geocentric;
    struct reb_particle* last_state;
    struct reb_particle* current_state;
    //particle_params* particle_params;
    double* particle_params;
    int steps_done;
    int forces;
};

/**
 * @brief Adds ASSIST functionality to a passed REBOUND simulation.
 * @param sim Pointer to the reb_simulation on which to add ASSIST functionality.
 * @return Pointer to an assist_extras structure.
 */
struct assist_extras* assist_attach(struct reb_simulation* sim, struct assist_ephem* ephem);

/**
 * @brief Frees all memory allocated by ASSIST instance.
 * @details Should be called after simulation is done if memory is a concern.
 * @param assist The assist_extras pointer returned from the initial call to assist_attach.
 */
void assist_free(struct assist_extras* assist);

void assist_ephem_free(struct assist_ephem* ephem);

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


int assist_interpolate_simulation(struct reb_simulation* sim1, struct reb_simulation* sim2, double h);
struct reb_simulation* assist_create_interpolated_simulation(struct reb_simulationarchive* sa, double t);
void assist_integrate_or_interpolate(struct assist_extras* ax, double t);

// Find particle position and velocity based on ephemeris data
struct reb_particle assist_get_particle(struct assist_ephem* ephem, const int particle_id, const double t);

// Functions called from python:
void assist_init(struct assist_extras* assist, struct reb_simulation* sim, struct assist_ephem* ephem);
void assist_free_pointers(struct assist_extras* assist);


void test_vary(struct reb_simulation* sim, FILE *vfile);

void test_vary_2nd(struct reb_simulation* sim, FILE *vfile);

struct assist_ephem* assist_ephem_create(char *planets_file_name, char *asteroids_file_name);

#endif
