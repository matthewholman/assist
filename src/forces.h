/**
 * @file    forces.h
 * @brief   Implementation of physical fources.
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
#ifndef _ASSIST_FORCES_H
#define _ASSIST_FORCES_H

// This is the force routine called by REBOUND. 
// It includes all forces, including gravity.
void assist_additional_forces(struct reb_simulation* sim);

int assist_all_ephem(struct assist_ephem* ephem, struct assist_ephem_cache* cache, const int i, const double t, double* const GM,
		      double* const x, double* const y, double* const z,
		      double* const vx, double* const vy, double* const vz,
		      double* const ax, double* const ay, double* const az
		      );

#endif
