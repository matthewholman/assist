
#ifndef _JPL_EPHEM_H
#define _JPL_EPHEM_H

struct _jpl_s * assist_jpl_init(void);
int assist_jpl_free(struct _jpl_s *jpl);
void assist_jpl_work(double *P, int ncm, int ncf, int niv, double t0, double t1, double *u, double *v, double *w);
int assist_jpl_calc(struct _jpl_s *pl, struct mpos_s *pos, double jd_ref, double jd_rel, int body);
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

#endif // _JPL_EPHEM_H

