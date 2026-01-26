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
#include <fcntl.h>
#include <unistd.h>
#include "assist.h"
#include "rebound.h"

#include "spk.h"
#include "forces.h"
#include "ascii_ephem.h"

#define STRINGIFY(s) str(s)
#define str(s) #s

const char* assist_build_str = __DATE__ " " __TIME__;   // Date and time build string. 
const char* assist_version_str = "1.1.9";         // **VERSIONLINE** This line gets updated automatically. Do not edit manually.
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

// -----------------------------
// Ephemeris format + discovery
// -----------------------------

int assist_detect_ascii_bin_signature(int fd) {
    // ASCII-derived binary ephemeris files have constant names at offset 0x00FC (252 bytes).
    // We check whether the first few 6-byte constant names look plausible.

    char const_names[6 * 3];  // Read first 3 constant names (6 chars each)
    if (lseek(fd, 0x00FC, SEEK_SET) != 0x00FC) {
        return 0;
    }

    if (read(fd, const_names, sizeof(const_names)) != (ssize_t)sizeof(const_names)) {
        return 0;
    }

    int valid_names = 0;
    for (int i = 0; i < 3; i++) {
        const char* name = &const_names[i * 6];

        // Names should start with a letter and contain at least 2 alphanumerics,
        // with the remainder being spaces or (rarely) null padding.
        if ((name[0] >= 'A' && name[0] <= 'Z') || (name[0] >= 'a' && name[0] <= 'z')) {
            int has_letters_or_digits = 0;
            int has_only_valid_chars = 1;
            for (int j = 0; j < 6; j++) {
                unsigned char c = (unsigned char)name[j];
                if ((c >= 'A' && c <= 'Z') || (c >= 'a' && c <= 'z') ||
                    (c >= '0' && c <= '9')) {
                    has_letters_or_digits++;
                } else if (c == ' ' || c == '\0') {
                    // padding is fine
                } else {
                    has_only_valid_chars = 0;
                    break;
                }
            }

            if (has_only_valid_chars && has_letters_or_digits >= 2) {
                valid_names++;
            }
        }
    }

    // If at least 2 out of 3 names look like ephemeris constants, it's likely a .440/.441.
    return (valid_names >= 2);
}

ephemeris_file_format_t assist_detect_ephemeris_file_format(int fd) {
    char buf[1024];
    ssize_t bytes_read;

    // Read first chunk of file
    lseek(fd, 0, SEEK_SET);
    bytes_read = read(fd, buf, sizeof(buf));
    if (bytes_read <= 0) {
        return FILE_FORMAT_UNKNOWN;
    }

    // Check for valid .bsp SPK file (DAF/SPK header) first
    if (bytes_read >= 8 && strncmp(buf, "DAF/SPK ", 8) == 0) {
        return FILE_FORMAT_VALID_BSP;
    }

    // Check for ASCII-derived binary ephemeris format signature
    if (assist_detect_ascii_bin_signature(fd)) {
        // Reset file position after signature check
        lseek(fd, 0, SEEK_SET);
        return FILE_FORMAT_ASCII_BIN;
    }

    return FILE_FORMAT_UNKNOWN;
}

int assist_discover_planets_path(char* out_path, size_t out_path_size, const char* assist_dir) {
    if (out_path == NULL || out_path_size == 0 || assist_dir == NULL) {
        return 0;
    }
    const char* candidates[] = {
        "/data/de441.bsp",
        "/data/de440.bsp",
        "/data/linux_m13000p17000.441",
        "/data/linux_p1550p2650.440",
    };
    for (int c = 0; c < 4; c++) {
        snprintf(out_path, out_path_size, "%s%s", assist_dir, candidates[c]);
        if (access(out_path, R_OK) == 0) {
            return 1;
        }
    }
    return 0;
}

    
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

    char default_planets_path[] = "/data/de441.bsp";
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
        // Implement discovery order under ASSIST_DIR
        const char* base = getenv("ASSIST_DIR");
        if (base){
            if (!assist_discover_planets_path(planets_path, FNAMESIZE, base)){
                // fall back to default path (de441.bsp) for error path consistency
                snprintf(planets_path, FNAMESIZE, "%s%s", base, default_planets_path);
            }
        }else{
            return ASSIST_ERROR_EPHEM_FILE;
        }
    }else{
        strncpy(planets_path, user_planets_path, FNAMESIZE-1);	
    }

    // Detect planets file format and load appropriate kernel
    ephem->spk_planets = NULL;
    ephem->ascii_planets = NULL;
    ephemeris_file_format_t planets_format = FILE_FORMAT_UNKNOWN;
    {
        int fd = open(planets_path, O_RDONLY);
        if (fd >= 0){
            planets_format = assist_detect_ephemeris_file_format(fd);
            close(fd);
        }
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
        if ((ephem->spk_asteroids = assist_spk_init(asteroids_path)) == NULL) {
            asteroids_path_not_found = 1;
        }
    }
    
    if (planets_format == FILE_FORMAT_VALID_BSP){
        // SPK planets path
        ephem->spk_planets = assist_spk_init(planets_path);
        if (ephem->spk_planets == NULL){
            return ASSIST_ERROR_EPHEM_FILE;
        }
        ephem->planets_source = FILE_FORMAT_VALID_BSP;
        // Precompute SPK planet indices
        for (int k=0; k<ASSIST_BODY_NPLANETS; k++) ephem->spk_target_index[k] = -1;
        ephem->spk_emb_index = -1;
        // Match the NAIF target codes used by `assist_spk_calc_planets_by_assist` (and the
        // ASCII-derived (.440/.441) ephemeris column meanings): planet *barycenters* for all planets except
        // that Earth/Moon are handled explicitly via 399/301 + EMB=3.
        static const int naif_by_assist[] = { 10, 1, 2, 399, 301, 4, 5, 6, 7, 8, 9 };
        for (int k=0; k<ASSIST_BODY_NPLANETS; k++){
            struct spk_target* t = assist_spk_find_target(ephem->spk_planets, naif_by_assist[k]);
            if (t){
                ephem->spk_target_index[k] = (int)(t - ephem->spk_planets->targets);
            }
        }
        // EMB index
        {
            struct spk_target* emb = assist_spk_find_target(ephem->spk_planets, 3);
            if (emb){ ephem->spk_emb_index = (int)(emb - ephem->spk_planets->targets); }
        }
        // SPK planets path: load constants/masses from SPK comments and join
        struct spk_constants_and_masses data = assist_load_spk_constants_and_masses(planets_path);
        assist_apply_spk_constants(ephem, &data);
        if (ephem->spk_planets) {
            assist_spk_join_masses(ephem->spk_planets, &data.masses, ephem->EMRAT);
        }
        if (ephem->spk_asteroids) {
            assist_spk_join_masses(ephem->spk_asteroids, &data.masses, ephem->EMRAT);
        }
        assist_free_spk_constants_and_masses(&data);
        ephem->planets_calc = assist_spk_calc_planets_by_assist;
    }else if (planets_format == FILE_FORMAT_ASCII_BIN){
        // Try ASCII-derived binary ephemeris (.440/.441)
        ephem->ascii_planets = assist_ascii_init(planets_path);
        if (ephem->ascii_planets == NULL){
            return ASSIST_ERROR_EPHEM_FILE;
        }
        ephem->planets_source = FILE_FORMAT_ASCII_BIN;
        // Copy constants from the ASCII-derived binary ephemeris
        ephem->J2E = ephem->ascii_planets->J2E;
        ephem->J3E = ephem->ascii_planets->J3E;
        ephem->J4E = ephem->ascii_planets->J4E;
        ephem->J2SUN = ephem->ascii_planets->J2SUN;
        ephem->AU = ephem->ascii_planets->AU;
        ephem->RE = ephem->ascii_planets->RE;
        ephem->CLIGHT = ephem->ascii_planets->CLIGHT;
        ephem->ASUN = ephem->ascii_planets->ASUN;
        ephem->EMRAT = ephem->ascii_planets->cem;
        ephem->Re_eq = ephem->RE / ephem->AU;
        ephem->Rs_eq = ephem->ASUN / ephem->AU;
        ephem->c_AU_per_day = (ephem->CLIGHT / ephem->AU) * 86400.0;
        ephem->c_squared = ephem->c_AU_per_day * ephem->c_AU_per_day;
        ephem->over_c_squared = 1.0 / ephem->c_squared;
        ephem->planets_calc = assist_ascii_calc_from_ephem;
        // If asteroids SPK loaded, set asteroid masses from JPL constants
        if (ephem->spk_asteroids){
            for (int n=0; n<ephem->spk_asteroids->num; n++){
                int code = ephem->spk_asteroids->targets[n].code;
                // Asteroid SPK targets use NAIF codes like 2000000 + (asteroid number).
                // Map to MAxxxx constants in the .440/.441 constants table.
                if (code >= 2000000){
                    char key[7];
                    snprintf(key, sizeof(key), "MA%04d", code - 2000000);
                    double gm = 0.0;
                    if (assist_ascii_find_constant(ephem->ascii_planets, key, &gm)){
                        ephem->spk_asteroids->targets[n].mass = gm;
                    }
                } else if (code == 399) {
                    // (Unlikely in sb*.bsp, but keep consistent if present.)
                    ephem->spk_asteroids->targets[n].mass = ephem->ascii_planets->mass[ASSIST_BODY_EARTH];
                } else if (code == 301) {
                    ephem->spk_asteroids->targets[n].mass = ephem->ascii_planets->mass[ASSIST_BODY_MOON];
                } else if (code == 10) {
                    ephem->spk_asteroids->targets[n].mass = ephem->ascii_planets->mass[ASSIST_BODY_SUN];
                }
            }
        }
    }else{
        // Provide a helpful error message describing what went wrong and how to fix it.
        // We deliberately do not try to interpret ASCII source or unknown formats for planets.
        fprintf(stderr, "(ASSIST) Error: Failed to initialize planets ephemeris from '%s'.\n", planets_path);
        fprintf(stderr, "(ASSIST) ASSIST supports NAIF SPK kernels (.bsp, e.g. de440.bsp) and JPL binary ephemerides (.440/.441, e.g. linux_p1550p2650.440).\n");
        fprintf(stderr, "(ASSIST) See the README section 'Planet ephemeris formats (DE440)' for download instructions.\n");
        return ASSIST_ERROR_EPHEM_FILE;
    }
            
    if (asteroids_path_not_found == 1){
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
    if (ephem->spk_planets != NULL){
        assist_spk_free(ephem->spk_planets);
    }
    if (ephem->spk_asteroids != NULL){
        assist_spk_free(ephem->spk_asteroids);
    }
    if (ephem->ascii_planets != NULL){
        assist_ascii_free(ephem->ascii_planets);
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
    if (ephem->spk_asteroids){
        N_total += ephem->spk_asteroids->num;
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
    // Default values for non-gravitational parameters
	assist->alpha = 1.0;
	assist->nk = 0.0;
	assist->nm = 2.0;
	assist->nn = 5.093;
	assist->r0 = 1.0;

    sim->integrator = REB_INTEGRATOR_IAS15;
    sim->gravity = REB_GRAVITY_NONE;
    sim->extras = assist;
    sim->extras_cleanup = assist_extras_cleanup;
    sim->additional_forces = assist_additional_forces;
    sim->force_is_velocity_dependent = 1;
    sim->ri_ias15.adaptive_mode = 1; // Use legacy IAS15 timestepping mode
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
        reb_simulation_error(assist->sim, msg);
    }
}


struct reb_particle assist_get_particle_with_error(const struct assist_ephem* ephem, const int particle_id, const double t, int* error){
    struct reb_particle p = {0};
    double GM = 0;
    int flag = assist_all_ephem(ephem, NULL, particle_id, t, &GM, &p.x, &p.y, &p.z, &p.vx, &p.vy, &p.vz, &p.ax, &p.ay, &p.az);
    *error = flag;
    p.m = GM; // Note this is GM, not M
    return p;
}


struct reb_particle assist_get_particle(const struct assist_ephem* ephem, const int particle_id, const double t){
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
    enum reb_simulation_binary_error_codes warnings = REB_SIMULATION_BINARY_WARNING_NONE;
    struct reb_simulation* r2 = reb_simulation_create();
    reb_simulation_create_from_simulationarchive_with_messages(r2, sa, blob-1, &warnings);
    // r2 = reb_input_process_warnings(r2, warnings); Ignoring warnings for now
    
    struct reb_simulation* r3 = reb_simulation_create();
    reb_simulation_create_from_simulationarchive_with_messages(r3, sa, blob, &warnings);
    // r3 = reb_input_process_warnings(r3, warnings);
    
    double h = (t - r2->t)/(r3->dt_last_done);
    assist_interpolate_simulation(r2, r3, h);
    reb_simulation_free(r3);
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
        ax->current_state = malloc(sizeof(struct reb_particle)*sim->N);
        ax->last_state = malloc(sizeof(struct reb_particle)*sim->N);
        // Initialize new arrays with sim->particles (sets mass, radius, hash, if users use those) 
        memcpy(ax->current_state, sim->particles, sizeof(struct reb_particle)*sim->N);
        memcpy(ax->last_state, sim->particles, sizeof(struct reb_particle)*sim->N);
    }else{
        assist_swap_particles(sim);
    }

    double dts = copysign(1., sim->dt_last_done);
    if ( !(dts*(sim->t-sim->dt_last_done)  <  dts*t &&  dts*t < dts*sim->t) ){
        // Integrate if requested time not in interval of last timestep
        
        reb_simulation_integrate(sim, t);

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
    reb_simulation_update_acceleration(sim); // This will later be recalculated. Could be optimized.
    memcpy(assist->last_state, sim->particles, sizeof(struct reb_particle)*sim->N);
}


