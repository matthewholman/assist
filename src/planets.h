#include "assist.h"
#ifndef _JPL_EPHEM_H
#define _JPL_EPHEM_H

struct _jpl_s * assist_jpl_init(char *str);
int assist_jpl_free(struct _jpl_s *jpl);
void assist_jpl_work(double *P, int ncm, int ncf, int niv, double t0, double t1, double *u, double *v, double *w);
int assist_jpl_calc(struct assist_ephem* ephem, const int body, const double jd_ref, const double jd_rel, struct assist_cache_item* const result);
double assist_jpl_mass(struct _jpl_s *pl, int tar);

// Order of columns in JPL Ephemeris file
enum {
         JPL_MER,                        // Mercury
         JPL_VEN,                        // Venus
         JPL_EMB,                        // Earth
         JPL_MAR,                        // Mars
         JPL_JUP,                        // Jupiter
         JPL_SAT,                        // Saturn
         JPL_URA,                        // Uranus
         JPL_NEP,                        // Neptune
         JPL_PLU,                        // Pluto
         JPL_LUN,                        // Moon (geocentric)
         JPL_SUN,                        // the Sun
         JPL_NUT,                        // nutations
         JPL_LIB,                        // lunar librations
         JPL_MAN,                        // lunar mantle
         JPL_TDB,                        // TT-TDB (< 2 ms)

         JPL_N,                          // Number of columns
 };

struct _jpl_s {
        double beg, end;                // begin and end times
        double inc;                     // time step size
        double cau;                     // definition of AU
        double cau_over;                // inverse of cau
        double cem;                     // Earth/Moon mass ratio
        int32_t num;                    // number of constants
        int32_t ver;                    // ephemeris version
        int32_t off[JPL_N];                // indexing offset
        int32_t ncf[JPL_N];                // number of chebyshev coefficients
        int32_t niv[JPL_N];                // number of interpolation intervals
        int32_t ncm[JPL_N];                // number of components / dimension
///
        size_t len, rec;                // file and record sizes
        void *map;                      // memory mapped location
	double *con;			// constant values
	char **str;			// constant names
};

#endif // _JPL_EPHEM_H

