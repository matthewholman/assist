/**
 * @file    assist.c
 * @brief   Central internal functions for ASSIST
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

/* Main routines called each timestep. */
#include <stdint.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <limits.h>
#include <float.h>
#include "assist.h"
#include "rebound.h"

#define STRINGIFY(s) str(s)
#define str(s) #s

const char* assist_build_str = __DATE__ " " __TIME__; // Date and time build string. 
const char* assist_version_str = "3.4.1";         // **VERSIONLINE** This line gets updated automatically. Do not edit manually.
const char* assist_githash_str = STRINGIFY(REBXGITHASH);             // This line gets updated automatically. Do not edit manually.



struct assist_extras* assist_attach(struct reb_simulation* sim){  
    if (sim == NULL){
        fprintf(stderr, "ASSIST Error: Simulation pointer passed to assist_attach was NULL.\n");
        return NULL;
    }
    struct assist_extras* assist = malloc(sizeof(*assist));
    // Initialization seperate from memory allocation because python handles memory management
    assist_initialize(sim, assist); 
    return assist;
}

void assist_initialize(struct reb_simulation* sim, struct assist_extras* assist){
    // Initialize assist, connect to rebound, ...
    assist->sim = sim;
}

void assist_free_pointers(struct assist_extras* assist){
    assist_detach(assist->sim, assist);
}

void assist_free(struct assist_extras* assist){
    // Freeing pointers is seperate because python handles memory management of structure itself.
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
