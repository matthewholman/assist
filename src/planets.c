#include <sys/mman.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <fcntl.h>
#include <string.h>
#include <unistd.h>

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "spk.h"
#include "planets.h"
#include "const.h"
#include "assist.h"

#ifndef JPL_EPHEM_FILE
#define JPL_EPHEM_FILE "../../data/linux_m13000p17000.441"
#endif

/*
 *  assist_jpl_work
 *
 *  Interpolate the appropriate Chebyshev polynomial coefficients.
 *
 *      ncf - number of coefficients per component
 *      ncm - number of components (ie: 3 for most)
 *      niv - number of intervals / sets of coefficients
 *
 */

void assist_jpl_work(double *P, int ncm, int ncf, int niv, double t0, double t1, double *u, double *v, double *w)
{
        double T[24], S[24];
        double U[24];
        double t, c;
        int p, m, n, b;

        // adjust to correct interval
        t = t0 * (double)niv;
        t0 = 2.0 * fmod(t, 1.0) - 1.0;
        c = (double)(niv * 2) / t1 / 86400.0;

        b = (int)t;

        // set up Chebyshev polynomials and derivatives
        T[0] = 1.0; T[1] = t0;
        S[0] = 0.0; S[1] = 1.0;
        U[0] = 0.0; U[1] = 0.0;	U[2] = 4.0;

        for (p = 2; p < ncf; p++) {
                T[p] = 2.0 * t0 * T[p-1] - T[p-2];
                S[p] = 2.0 * t0 * S[p-1] + 2.0 * T[p-1] - S[p-2];
        }
        for (p = 3; p < ncf; p++) {
                U[p] = 2.0 * t0 * U[p-1] + 4.0 * S[p-1] - U[p-2];
        }

        // compute the position/velocity
        for (m = 0; m < ncm; m++) {
                u[m] = v[m] = w[m] = 0.0;
                n = ncf * (m + b * ncm);

                for (p = 0; p < ncf; p++) {
                        u[m] += T[p] * P[n+p];
                        v[m] += S[p] * P[n+p] * c;
                        w[m] += U[p] * P[n+p] * c * c;
                }
        }
}
 
/*
 *  assist_jpl_init
 *
 *  Initialise everything needed ... probaly not be compatible with a non-430 file.
 *
 */

struct jpl_s * assist_jpl_init(char *str)
{
        struct jpl_s *jpl;
	struct stat sb;
	//char *str;
	ssize_t ret;
	off_t off;
        int fd, p;

        /** use or environment-specified file, 
	 * or the default filename, in that order
         */
	//if ((str = getenv("JPL_PLANET_EPHEM")) == NULL)
	//str = JPL_EPHEM_FILE;

        if ((fd = open(str, O_RDONLY)) < 0)
                return NULL;

        jpl = malloc(sizeof(struct jpl_s));
        memset(jpl, 0, sizeof(struct jpl_s));

        if (fstat(fd, &sb) < 0)
                goto err;

	// FIXME : probably should ensure the file is sized corrrectly
	// FIXME : also could read 3*84 bytes to see if this is a JPL file

	// skip the header and constant names for now
        if (lseek(fd, 0x0A5C, SEEK_SET) < 0)
                goto err;

        // read header
        ret  = read(fd, &jpl->beg, sizeof(double));
        ret += read(fd, &jpl->end, sizeof(double));
        ret += read(fd, &jpl->inc, sizeof(double));
        ret += read(fd, &jpl->num, sizeof(int32_t));
        ret += read(fd, &jpl->cau, sizeof(double));
        ret += read(fd, &jpl->cem, sizeof(double));

	//printf("%lf %lf %lf %d\n", jpl->beg, jpl->end, jpl->inc, jpl->num);

        // number of coefficients for all components
        for (p = 0; p < JPL_N; p++)
                jpl->ncm[p] = 3;

        // exceptions:
        jpl->ncm[JPL_NUT] = 2; // nutations
        jpl->ncm[JPL_TDB] = 1; // TT-TDB

        for (p = 0; p < 12; p++) {
                ret += read(fd, &jpl->off[p], sizeof(int32_t));
                ret += read(fd, &jpl->ncf[p], sizeof(int32_t));
                ret += read(fd, &jpl->niv[p], sizeof(int32_t));
        }

        ret += read(fd, &jpl->ver,     sizeof(int32_t));
        ret += read(fd, &jpl->off[12], sizeof(int32_t));
        ret += read(fd, &jpl->ncf[12], sizeof(int32_t));
        ret += read(fd, &jpl->niv[12], sizeof(int32_t));

	// get all the constant names, from two lcoations
	jpl->str = calloc(jpl->num, sizeof(char *));
	off = lseek(fd, 0, SEEK_CUR);
	lseek(fd, 0x00FC, SEEK_SET);

	// retrieve the names of the first 400 constants
	for (p = 0; p < 400; p++) {
		jpl->str[p] = calloc(1, 8);
		read(fd, jpl->str[p], 6);
	}

	lseek(fd, off, SEEK_SET);

	// read the remaining constant names
	for (p = 400; p < jpl->num; p++) {
		jpl->str[p] = calloc(1, 8);
		read(fd, jpl->str[p], 6);
	}

        // finishing reading
        for (p = 13; p < 15; p++) {
                ret += read(fd, &jpl->off[p], sizeof(int32_t));
                ret += read(fd, &jpl->ncf[p], sizeof(int32_t));
                ret += read(fd, &jpl->niv[p], sizeof(int32_t));
        }

        // adjust for correct indexing (ie: zero based)
        for (p = 0; p < JPL_N; p++)
                jpl->off[p] -= 1;

        // save file size, and determine 'kernel size'
        jpl->len = sb.st_size;
        jpl->rec = sizeof(double) * 2;

        for (p = 0; p < JPL_N; p++)
                jpl->rec += sizeof(double) * jpl->ncf[p] * jpl->niv[p] * jpl->ncm[p];

        // memory map the file, which makes us thread-safe with kernel caching
        jpl->map = mmap(NULL, jpl->len, PROT_READ, MAP_SHARED, fd, 0);

        if (jpl->map == NULL)
                goto err;

	// now read the constant values after seeking to where they are
	if (lseek(fd, jpl->rec, SEEK_SET) < 0)
		goto err;

	jpl->con = calloc(jpl->num, sizeof(double));

	for (p = 0; p < jpl->num; p++)
		read(fd, &jpl->con[p], sizeof(double));

        // this file descriptor is no longer needed since we are memory mapped
        if (close(fd) < 0)
                { ; } // perror ...
#if defined(MADV_RANDOM)
        if (madvise(jpl->map, jpl->len, MADV_RANDOM) < 0)
                { ; } // perror ...
#endif

        return jpl;

err:    close(fd);
        free(jpl);

        return NULL;
}

/*
 *  assist_jpl_free
 *
 */
int assist_jpl_free(struct jpl_s *jpl)
{
	int p;

        if (jpl == NULL)
                return -1;

        if (munmap(jpl->map, jpl->len) < 0)
                { ; } // perror...

	for (p = 0; p < jpl->num; p++)
		free(jpl->str[p]);

	free(jpl->str);
	free(jpl->con);
        memset(jpl, 0, sizeof(struct jpl_s));
        free(jpl);
        return 0;
}

/*
 *  jpl_calc
 *
 *  Caculate the position+velocity in _equatorial_ coordinates.
 *  Assumes pos is initially zero.
 */
enum ASSIST_STATUS assist_jpl_calc(struct jpl_s *jpl, double jd_ref, double jd_rel, int body, 
		 double* const GM,
		 double* const out_x, double* const out_y, double* const out_z,
		 double* const out_vx, double* const out_vy, double* const out_vz,
         double* const out_ax, double* const out_ay, double* const out_az){
    double t, *z;
    u_int32_t blk;

    if (jpl == NULL || jpl->map == NULL)
        return ASSIST_ERROR_EPHEM_FILE;
    if(body<0 || body >= ASSIST_BODY_NPLANETS)
	    return(ASSIST_ERROR_NEPHEM);
    
    struct mpos_s pos;

    // The values below are G*mass.
    // Units are solar masses, au, days.
    // DE440/441 units: au^3 day^-2.
    const static double JPL_GM[ASSIST_BODY_NPLANETS] =
	{
	    JPL_EPHEM_GMS, // 0 sun
	    JPL_EPHEM_GM1, // 1 mercury
	    JPL_EPHEM_GM2, // 2 venus
	    (JPL_EPHEM_EMRAT/(1.+JPL_EPHEM_EMRAT) * JPL_EPHEM_GMB),  // 3 earth    Calculate GM values for Earth and Moon
	    (1./(1.+JPL_EPHEM_EMRAT) * JPL_EPHEM_GMB),               // 4 moon     from Earth-moon ratio and sum.
	    JPL_EPHEM_GM4, // 5 mars
	    JPL_EPHEM_GM5, // 6 jupiter
	    JPL_EPHEM_GM6, // 7 saturn
	    JPL_EPHEM_GM7, // 8 uranus
	    JPL_EPHEM_GM8, // 9 neptune
	    JPL_EPHEM_GM9, // 10 pluto
	};


    // Get position, velocity, and mass of body i in barycentric coords.

    *GM = JPL_GM[body];

        // check if covered by this file
        if (jd_ref + jd_rel < jpl->beg || jd_ref + jd_rel > jpl->end)
            return ASSIST_ERROR_EPHEM_FILE;

        // compute record number and 'offset' into record
        blk = (u_int32_t)((jd_ref + jd_rel - jpl->beg) / jpl->inc);
        z = (double*)jpl->map + (blk + 2) * jpl->rec/sizeof(double);
        t = ((jd_ref - jpl->beg - (double)blk * jpl->inc) + jd_rel) / jpl->inc;

        switch (body) { // The indices in the jpl-> arrays match the JPL component index for the body
            case ASSIST_BODY_SUN:
                assist_jpl_work(&z[jpl->off[10]], jpl->ncm[10], jpl->ncf[10], jpl->niv[10], t, jpl->inc, pos.u, pos.v, pos.w);
                break;
            case ASSIST_BODY_MERCURY:
                assist_jpl_work(&z[jpl->off[JPL_MER]], jpl->ncm[JPL_MER], jpl->ncf[JPL_MER], jpl->niv[JPL_MER], t, jpl->inc, pos.u, pos.v, pos.w);
                break;
            case ASSIST_BODY_VENUS:
                assist_jpl_work(&z[jpl->off[JPL_VEN]], jpl->ncm[JPL_VEN], jpl->ncf[JPL_VEN], jpl->niv[JPL_VEN], t, jpl->inc, pos.u, pos.v, pos.w);
                break;
            case ASSIST_BODY_EARTH:
                {
                    struct mpos_s emb, lun;
                    assist_jpl_work(&z[jpl->off[JPL_EMB]], jpl->ncm[JPL_EMB], jpl->ncf[JPL_EMB], jpl->niv[JPL_EMB], t, jpl->inc, emb.u, emb.v, emb.w); // earth moon barycenter
                    assist_jpl_work(&z[jpl->off[JPL_LUN]], jpl->ncm[JPL_LUN], jpl->ncf[JPL_LUN], jpl->niv[JPL_LUN], t, jpl->inc, lun.u, lun.v, lun.w);

                    vecpos_set(pos.u, emb.u);
                    vecpos_off(pos.u, lun.u, -1.0 / (1.0 + jpl->cem));

                    vecpos_set(pos.v, emb.v);
                    vecpos_off(pos.v, lun.v, -1.0 / (1.0 + jpl->cem));

                    vecpos_set(pos.w, emb.w);
                    vecpos_off(pos.w, lun.w, -1.0 / (1.0 + jpl->cem));
                }
                break;
            case ASSIST_BODY_MOON: 
                {
                    struct mpos_s emb, lun;
                    assist_jpl_work(&z[jpl->off[JPL_EMB]], jpl->ncm[JPL_EMB], jpl->ncf[JPL_EMB], jpl->niv[JPL_EMB], t, jpl->inc, emb.u, emb.v, emb.w);
                    assist_jpl_work(&z[jpl->off[JPL_LUN]], jpl->ncm[JPL_LUN], jpl->ncf[JPL_LUN], jpl->niv[JPL_LUN], t, jpl->inc, lun.u, lun.v, lun.w);

                    vecpos_set(pos.u, emb.u);
                    vecpos_off(pos.u, lun.u, jpl->cem / (1.0 + jpl->cem));

                    vecpos_set(pos.v, emb.v);
                    vecpos_off(pos.v, lun.v, jpl->cem / (1.0 + jpl->cem));

                    vecpos_set(pos.w, emb.w);
                    vecpos_off(pos.w, lun.w, jpl->cem / (1.0 + jpl->cem));
                }
                break;
            case ASSIST_BODY_MARS:
                assist_jpl_work(&z[jpl->off[JPL_MAR]], jpl->ncm[JPL_MAR], jpl->ncf[JPL_MAR], jpl->niv[JPL_MAR], t, jpl->inc, pos.u, pos.v, pos.w);
                break;
            case ASSIST_BODY_JUPITER:
                assist_jpl_work(&z[jpl->off[JPL_JUP]], jpl->ncm[JPL_JUP], jpl->ncf[JPL_JUP], jpl->niv[JPL_JUP], t, jpl->inc, pos.u, pos.v, pos.w);
                break;
            case ASSIST_BODY_SATURN:
                assist_jpl_work(&z[jpl->off[JPL_SAT]], jpl->ncm[JPL_SAT], jpl->ncf[JPL_SAT], jpl->niv[JPL_SAT], t, jpl->inc, pos.u, pos.v, pos.w);
                break;
            case ASSIST_BODY_URANUS:
                assist_jpl_work(&z[jpl->off[JPL_URA]], jpl->ncm[JPL_URA], jpl->ncf[JPL_URA], jpl->niv[JPL_URA], t, jpl->inc, pos.u, pos.v, pos.w);
                break;
            case ASSIST_BODY_NEPTUNE:
                assist_jpl_work(&z[jpl->off[JPL_NEP]], jpl->ncm[JPL_NEP], jpl->ncf[JPL_NEP], jpl->niv[JPL_NEP], t, jpl->inc, pos.u, pos.v, pos.w);
                break;
            case ASSIST_BODY_PLUTO:
                assist_jpl_work(&z[jpl->off[JPL_PLU]], jpl->ncm[JPL_PLU], jpl->ncf[JPL_PLU], jpl->niv[JPL_PLU], t, jpl->inc, pos.u, pos.v, pos.w);
                break;
            default:
                return ASSIST_ERROR_NEPHEM; // body not found
                break;
        }

    // Convert to au/day and au/day^2
    vecpos_div(pos.u, jpl->cau);
    vecpos_div(pos.v, jpl->cau/86400.);
    vecpos_div(pos.w, jpl->cau/(86400.*86400.));

    *out_x = pos.u[0];
    *out_y = pos.u[1];
    *out_z = pos.u[2];
    *out_vx = pos.v[0];
    *out_vy = pos.v[1];
    *out_vz = pos.v[2];
    *out_ax = pos.w[0];
    *out_ay = pos.w[1];
    *out_az = pos.w[2];

    return(ASSIST_SUCCESS);

}

/*
 *  assist_jpl_mass
 *
 */
double assist_jpl_mass(struct jpl_s *jpl, int tar)
{
	char buf[14];
	int n;

	if (jpl == NULL)
		return 0.0;

	if (tar >= 10000)
		return 0.0;

	snprintf(buf, sizeof(buf), "MA%04d", tar);

	for (n = 0; n < jpl->num; n++) {
		if (strncmp(jpl->str[n], buf, 6) == 0)
			return jpl->con[n];
	}

	// not found
	return 0.0;
}

