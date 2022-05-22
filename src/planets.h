struct _jpl_s * jpl_init(void);
int jpl_free(struct _jpl_s *jpl);
void jpl_work(double *P, int ncm, int ncf, int niv, double t0, double t1, double *u, double *v, double *w);
int jpl_calc(struct _jpl_s *jpl, struct mpos_s *now, double jde, int n, int m);

// these are the body codes for the user to specify
enum {
        PLAN_BAR,                       // <0,0,0>
        PLAN_SOL,                       // Sun (in barycentric)
        PLAN_EAR,                       // Earth centre
        PLAN_EMB,                       // Earth-Moon barycentre
        PLAN_LUN,                       // Moon centre
        PLAN_MER,                       // ... plus the rest
        PLAN_VEN,
        PLAN_MAR,
        PLAN_JUP,
        PLAN_SAT,
        PLAN_URA,
        PLAN_NEP,
        PLAN_PLU,	

        _NUM_TEST,
};

int body[11];

// these are array indices for the internal interface
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

        _NUM_JPL,
};

struct _jpl_s {
        double beg, end;                // begin and end times
        double inc;                     // time step size
        double cau;                     // definition of AU
        double cem;                     // Earth/Moon mass ratio
        int32_t num;                    // number of constants
        int32_t ver;                    // ephemeris version
        int32_t off[_NUM_JPL];          // indexing offset
        int32_t ncf[_NUM_JPL];          // number of chebyshev coefficients
        int32_t niv[_NUM_JPL];          // number of interpolation intervals
        int32_t ncm[_NUM_JPL];          // number of components / dimension
///
        size_t len, rec;                // file and record sizes
        void *map;                      // memory mapped location
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
