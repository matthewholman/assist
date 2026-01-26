#ifndef _ASSIST_ASCII_EPHEM_H
#define _ASSIST_ASCII_EPHEM_H
// ASCII-derived ephemeris backend (.440/.441)
//
// These files are produced by converting the upstream JPL "ASCII blocks"
// ephemerides into a machine-dependent binary format.
//
// ASSIST uses these files as an alternative planets provider to SPK kernels
// (.bsp). The public-facing selector for this backend is `FILE_FORMAT_ASCII_BIN`.

#include "assist.h"

struct ascii_s * assist_ascii_init(char *path);
void assist_ascii_free(struct ascii_s *ascii);
void assist_ascii_work(double *P, int ncm, int ncf, int niv, double t0, double t1, double *u, double *v, double *w);
enum ASSIST_STATUS assist_ascii_calc(struct ascii_s *pl, double jd_ref, double jd_rel, int body, 
		 double* const GM,
		 double* const x, double* const y, double* const z,
		 double* const vx, double* const vy, double* const vz,
         double* const ax, double* const ay, double* const az);

enum ASSIST_STATUS assist_ascii_calc_from_ephem(const struct assist_ephem* ephem, double jd_ref, double jd_rel, int body, 
		 double* const GM,
		 double* const x, double* const y, double* const z,
		 double* const vx, double* const vy, double* const vz,
         double* const ax, double* const ay, double* const az);

// Lookup a 6-character constant (e.g. "GMS", "GM1", "MA0001") in a .440/.441 file.
// Returns 1 if found, 0 otherwise.
int assist_ascii_find_constant(const struct ascii_s* ascii, const char* name, double* out_value);

// Order of columns in the ASCII-derived ephemeris file
enum {
         ASCII_MER,                        // Mercury
         ASCII_VEN,                        // Venus
         ASCII_EMB,                        // Earth-Moon barycenter
         ASCII_MAR,                        // Mars
         ASCII_JUP,                        // Jupiter
         ASCII_SAT,                        // Saturn
         ASCII_URA,                        // Uranus
         ASCII_NEP,                        // Neptune
         ASCII_PLU,                        // Pluto
         ASCII_LUN,                        // Moon (geocentric)
         ASCII_SUN,                        // the Sun
         ASCII_NUT,                        // nutations
         ASCII_LIB,                        // lunar librations
         ASCII_MAN,                        // lunar mantle
         ASCII_TDB,                        // TT-TDB (< 2 ms)

         ASCII_N,                          // Number of columns
 };

struct ascii_s {
        double beg, end;                // begin and end times
        double inc;                     // time step size
        double cau;                     // definition of AU
        double cem;                     // Earth/Moon mass ratio
        int32_t num;                    // number of constants
        int32_t ver;                    // ephemeris version
        int32_t off[ASCII_N];           // indexing offset
        int32_t ncf[ASCII_N];           // number of chebyshev coefficients
        int32_t niv[ASCII_N];           // number of interpolation intervals
        int32_t ncm[ASCII_N];           // number of components / dimension
        double mass[ASCII_N];           // G*mass for all bodies (ASSIST_BODY indexing for 0..ASSIST_BODY_NPLANETS-1)
        double J2E;                     // Other constant names follow JPL
        double J3E;
        double J4E;
        double J2SUN;
        double AU;
        double RE;
        double CLIGHT;
        double ASUN;
        size_t len, rec;                // file and record sizes
        void *map;                      // memory mapped location
	    double *con;			        // constant values
	    char **str;			            // constant names
};

// From Weryk's code
/////// private interface :

static inline void vecpos_off(double *u, const double *v, const double w)
        { u[0] += v[0] * w; u[1] += v[1] * w; u[2] += v[2] * w; }
static inline void vecpos_set(double *u, const double *v)
        { u[0] = v[0]; u[1] = v[1]; u[2] = v[2]; }
static inline void vecpos_nul(double *u)
        { u[0] = u[1] = u[2] = 0.0; }
static inline void vecpos_div(double *u, double v)
        { u[0] /= v; u[1] /= v; u[2] /= v; }

#endif // _ASSIST_ASCII_EPHEM_H


