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

struct assist_extras {
    struct reb_simulation* sim;
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
 * @brief Output an error message.
 * @details This function should be used if an error occurs rather than simply using print. 
 *          The message will be passed to python.
 * @param assist The assist_extras pointer.
 * @param msg The error message.
 */
void assist_error(struct assist_extras* assist, const char* const msg);

#endif
