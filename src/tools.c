/**
 * @file    tools.c
 * @brief   Various tools for ASSIST
 * @author  Hanno Rein
 * 
 * @section     LICENSE
 * Copyright (c) 2025 Hanno Rein
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
#include <stdint.h>
#include <math.h>
#include "assist.h"
#include "rebound.h"

// This function returns a new rebound simulation with all particles added. 
// Return simulation must be freed by caller.
// if merge_moon = 1, then the Earth and Moon are added as a single particle.
struct reb_simulation* assist_simulation_convert_to_rebound(struct reb_simulation* r, struct assist_ephem* ephem, int merge_moon){  
    struct reb_simulation* r2 = reb_simulation_create();
    r2->t = r->t;
    r2->dt = r->dt;
    r2->ri_ias15.epsilon = r->ri_ias15.epsilon;
    r2->ri_ias15.adaptive_mode = r->ri_ias15.adaptive_mode;

    for (int i=0; i<11; i++){
        if (!merge_moon || i!=ASSIST_BODY_EARTH || i!=ASSIST_BODY_MOON){
            int error=0;
            struct reb_particle p = assist_get_particle_with_error(ephem, i, r->t,&error);
            if (error != ASSIST_SUCCESS){
                fprintf(stderr, "(ASSIST) An error occured while trying to initialize particle from ephemeris data.\n");
            }else{
                reb_simulation_add(r2, p);
            }
        }
    }
    if (merge_moon){
        int error =0;
        struct reb_particle p1 = assist_get_particle_with_error(ephem, ASSIST_BODY_EARTH, r->t,&error);
        struct reb_particle p2 = assist_get_particle_with_error(ephem, ASSIST_BODY_MOON, r->t,&error);
        if (error != ASSIST_SUCCESS){
            fprintf(stderr, "(ASSIST) An error occured while trying to initialize particle from ephemeris data.\n");
        }else{
            struct reb_particle p = reb_particle_com_of_pair(p1, p2);
            reb_simulation_add(r2, p);
        }
    }
    r2->N_active = r2->N;
    for (int i=0; i<r->N; i++){
        reb_simulation_add(r2, r->particles[i]);
    }
    return r2;
}
