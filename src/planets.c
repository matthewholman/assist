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

struct _jpl_s * assist_jpl_init(char *str)
{
        struct _jpl_s *jpl;
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

        jpl = malloc(sizeof(struct _jpl_s));
        memset(jpl, 0, sizeof(struct _jpl_s));

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
int assist_jpl_free(struct _jpl_s *jpl)
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
        memset(jpl, 0, sizeof(struct _jpl_s));
        free(jpl);
        return 0;
}

/*
 *  jpl_calc
 *
 *  Caculate the position+velocity in _equatorial_ coordinates.
 *  Assumes pos is initially zero.
 */
int assist_jpl_calc(struct _jpl_s *pl, struct mpos_s *pos, double jd_ref, double jd_rel, int body) {
        double t, *z;
        u_int32_t blk;

        if (pl == NULL || pl->map == NULL || pos == NULL)
                return -1;

        // check if covered by this file
        if (jd_ref + jd_rel < pl->beg || jd_ref + jd_rel > pl->end)
                return -1;

        // compute record number and 'offset' into record
        blk = (u_int32_t)((jd_ref + jd_rel - pl->beg) / pl->inc);
        z = (double*)pl->map + (blk + 2) * pl->rec/sizeof(double);
        t = ((jd_ref - pl->beg - (double)blk * pl->inc) + jd_rel) / pl->inc;

        switch (body) { // The indices in the pl-> arrays match the JPL component index for the body
            case 0: // SUN
                assist_jpl_work(&z[pl->off[10]], pl->ncm[10], pl->ncf[10], pl->niv[10], t, pl->inc, pos->u, pos->v, pos->w);
                break;
            case 1: // MER
                assist_jpl_work(&z[pl->off[JPL_MER]], pl->ncm[JPL_MER], pl->ncf[JPL_MER], pl->niv[JPL_MER], t, pl->inc, pos->u, pos->v, pos->w);
                break;
            case 2: // VEN
                assist_jpl_work(&z[pl->off[JPL_VEN]], pl->ncm[JPL_VEN], pl->ncf[JPL_VEN], pl->niv[JPL_VEN], t, pl->inc, pos->u, pos->v, pos->w);
                break;
            case 3: // EAR
                {
                    struct mpos_s emb, lun;
                    assist_jpl_work(&z[pl->off[JPL_EMB]], pl->ncm[JPL_EMB], pl->ncf[JPL_EMB], pl->niv[JPL_EMB], t, pl->inc, emb.u, emb.v, emb.w); // earth moon barycenter
                    assist_jpl_work(&z[pl->off[JPL_LUN]], pl->ncm[JPL_LUN], pl->ncf[JPL_LUN], pl->niv[JPL_LUN], t, pl->inc, lun.u, lun.v, lun.w);

                    vecpos_set(pos->u, emb.u);
                    vecpos_off(pos->u, lun.u, -1.0 / (1.0 + pl->cem));

                    vecpos_set(pos->v, emb.v);
                    vecpos_off(pos->v, lun.v, -1.0 / (1.0 + pl->cem));

                    vecpos_set(pos->w, emb.w);
                    vecpos_off(pos->w, lun.w, -1.0 / (1.0 + pl->cem));
                }
                break;
            case 4: // LUN 
                {
                    struct mpos_s emb, lun;
                    assist_jpl_work(&z[pl->off[JPL_EMB]], pl->ncm[JPL_EMB], pl->ncf[JPL_EMB], pl->niv[JPL_EMB], t, pl->inc, emb.u, emb.v, emb.w);
                    assist_jpl_work(&z[pl->off[JPL_LUN]], pl->ncm[JPL_LUN], pl->ncf[JPL_LUN], pl->niv[JPL_LUN], t, pl->inc, lun.u, lun.v, lun.w);

                    vecpos_set(pos->u, emb.u);
                    vecpos_off(pos->u, lun.u, pl->cem / (1.0 + pl->cem));

                    vecpos_set(pos->v, emb.v);
                    vecpos_off(pos->v, lun.v, pl->cem / (1.0 + pl->cem));

                    vecpos_set(pos->w, emb.w);
                    vecpos_off(pos->w, lun.w, pl->cem / (1.0 + pl->cem));
                }
                break;
            case 5: // MAR
                assist_jpl_work(&z[pl->off[JPL_MAR]], pl->ncm[JPL_MAR], pl->ncf[JPL_MAR], pl->niv[JPL_MAR], t, pl->inc, pos->u, pos->v, pos->w);
                break;
            case 6: // JUP
                assist_jpl_work(&z[pl->off[JPL_JUP]], pl->ncm[JPL_JUP], pl->ncf[JPL_JUP], pl->niv[JPL_JUP], t, pl->inc, pos->u, pos->v, pos->w);
                break;
            case 7: // SAT
                assist_jpl_work(&z[pl->off[JPL_SAT]], pl->ncm[JPL_SAT], pl->ncf[JPL_SAT], pl->niv[JPL_SAT], t, pl->inc, pos->u, pos->v, pos->w);
                break;
            case 8: // URA
                assist_jpl_work(&z[pl->off[JPL_URA]], pl->ncm[JPL_URA], pl->ncf[JPL_URA], pl->niv[JPL_URA], t, pl->inc, pos->u, pos->v, pos->w);
                break;
            case 9: // NEP
                assist_jpl_work(&z[pl->off[JPL_NEP]], pl->ncm[JPL_NEP], pl->ncf[JPL_NEP], pl->niv[JPL_NEP], t, pl->inc, pos->u, pos->v, pos->w);
                break;
            case 10: // PLU
                assist_jpl_work(&z[pl->off[JPL_PLU]], pl->ncm[JPL_PLU], pl->ncf[JPL_PLU], pl->niv[JPL_PLU], t, pl->inc, pos->u, pos->v, pos->w);
                break;
            case 11: // BAR
                     // Nothing needs to be done
                break;
            default:
                return -1; // body not found
                break;
        }

        pos->jde = jd_ref + jd_rel;
        return 0;
}

/*
 *  assist_jpl_mass
 *
 */
double assist_jpl_mass(struct _jpl_s *pl, int tar)
{
	char buf[14];
	int n;

	if (pl == NULL)
		return 0.0;

	if (tar >= 10000)
		return 0.0;

	snprintf(buf, sizeof(buf), "MA%04d", tar);

	for (n = 0; n < pl->num; n++) {
		if (strncmp(pl->str[n], buf, 6) == 0)
			return pl->con[n];
	}

	// not found
	return 0.0;
}

