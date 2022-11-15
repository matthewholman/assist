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

int body[11] = {
        PLAN_SOL,                       // Sun (in barycentric)
        PLAN_MER,                       // Mercury center
        PLAN_VEN,                       // Venus center
        PLAN_EAR,                       // Earth center
        PLAN_LUN,                       // Moon center
        PLAN_MAR,                       // Mars center
        PLAN_JUP,                       // ...
        PLAN_SAT,
        PLAN_URA,
        PLAN_NEP,
        PLAN_PLU
};

/*
 *  jpl_work
 *
 *  Interpolate the appropriate Chebyshev polynomial coefficients.
 *
 *      ncf - number of coefficients per component
 *      ncm - number of components (ie: 3 for most)
 *      niv - number of intervals / sets of coefficients
 *
 */

void jpl_work(double *P, int ncm, int ncf, int niv, double t0, double t1, double *u, double *v, double *w)
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
 *  jpl_init
 *
 *  Initialise everything needed ... probaly not be compatible with a non-430 file.
 *
 */

struct _jpl_s * jpl_init(void)
{
        struct _jpl_s *jpl;
	struct stat sb;
	char *str;
	ssize_t ret;
	off_t off;
        int fd, p;

        /** use or environment-specified file, 
	 * or the default filename, in that order
         */
	if ((str = getenv("JPL_PLANET_EPHEM")) == NULL)
		str = JPL_EPHEM_FILE;

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

	printf("%lf %lf %lf %d\n", jpl->beg, jpl->end, jpl->inc, jpl->num);

        // number of coefficients is assumed
        for (p = 0; p < _NUM_JPL; p++)
                jpl->ncm[p] = 3;

        jpl->ncm[JPL_NUT] = 2;
        jpl->ncm[JPL_TDB] = 1;

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
        for (p = 0; p < _NUM_JPL; p++)
                jpl->off[p] -= 1;

        // save file size, and determine 'kernel size'
        jpl->len = sb.st_size;
        jpl->rec = sizeof(double) * 2;

        for (p = 0; p < _NUM_JPL; p++)
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
        if (madvise(jpl->map, jpl->len, MADV_RANDOM) < 0)
                { ; } // perror ...

        return jpl;

err:    close(fd);
        free(jpl);

        return NULL;
}

/*
 *  jpl_free
 *
 */
int jpl_free(struct _jpl_s *jpl)
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
 *
 */

static void _bar(struct _jpl_s *jpl, double *z, double t, struct mpos_s *pos)
        { vecpos_nul(pos->u); vecpos_nul(pos->v); vecpos_nul(pos->w); }
static void _sun(struct _jpl_s *jpl, double *z, double t, struct mpos_s *pos)
        { jpl_work(&z[jpl->off[JPL_SUN]], jpl->ncm[JPL_SUN], jpl->ncf[JPL_SUN], jpl->niv[JPL_SUN], t, jpl->inc, pos->u, pos->v, pos->w); }
static void _emb(struct _jpl_s *jpl, double *z, double t, struct mpos_s *pos)
        { jpl_work(&z[jpl->off[JPL_EMB]], jpl->ncm[JPL_EMB], jpl->ncf[JPL_EMB], jpl->niv[JPL_EMB], t, jpl->inc, pos->u, pos->v, pos->w); }
static void _mer(struct _jpl_s *jpl, double *z, double t, struct mpos_s *pos)
        { jpl_work(&z[jpl->off[JPL_MER]], jpl->ncm[JPL_MER], jpl->ncf[JPL_MER], jpl->niv[JPL_MER], t, jpl->inc, pos->u, pos->v, pos->w); }
static void _ven(struct _jpl_s *jpl, double *z, double t, struct mpos_s *pos)
        { jpl_work(&z[jpl->off[JPL_VEN]], jpl->ncm[JPL_VEN], jpl->ncf[JPL_VEN], jpl->niv[JPL_VEN], t, jpl->inc, pos->u, pos->v, pos->w); }
static void _mar(struct _jpl_s *jpl, double *z, double t, struct mpos_s *pos)
        { jpl_work(&z[jpl->off[JPL_MAR]], jpl->ncm[JPL_MAR], jpl->ncf[JPL_MAR], jpl->niv[JPL_MAR], t, jpl->inc, pos->u, pos->v, pos->w); }
static void _jup(struct _jpl_s *jpl, double *z, double t, struct mpos_s *pos)
        { jpl_work(&z[jpl->off[JPL_JUP]], jpl->ncm[JPL_JUP], jpl->ncf[JPL_JUP], jpl->niv[JPL_JUP], t, jpl->inc, pos->u, pos->v, pos->w); }
static void _sat(struct _jpl_s *jpl, double *z, double t, struct mpos_s *pos)
        { jpl_work(&z[jpl->off[JPL_SAT]], jpl->ncm[JPL_SAT], jpl->ncf[JPL_SAT], jpl->niv[JPL_SAT], t, jpl->inc, pos->u, pos->v, pos->w); }
static void _ura(struct _jpl_s *jpl, double *z, double t, struct mpos_s *pos)
        { jpl_work(&z[jpl->off[JPL_URA]], jpl->ncm[JPL_URA], jpl->ncf[JPL_URA], jpl->niv[JPL_URA], t, jpl->inc, pos->u, pos->v, pos->w); }
static void _nep(struct _jpl_s *jpl, double *z, double t, struct mpos_s *pos)
        { jpl_work(&z[jpl->off[JPL_NEP]], jpl->ncm[JPL_NEP], jpl->ncf[JPL_NEP], jpl->niv[JPL_NEP], t, jpl->inc, pos->u, pos->v, pos->w); }
static void _plu(struct _jpl_s *jpl, double *z, double t, struct mpos_s *pos)
        { jpl_work(&z[jpl->off[JPL_PLU]], jpl->ncm[JPL_PLU], jpl->ncf[JPL_PLU], jpl->niv[JPL_PLU], t, jpl->inc, pos->u, pos->v, pos->w); }

static void _ear(struct _jpl_s *jpl, double *z, double t, struct mpos_s *pos)
{
        struct mpos_s emb, lun;
        jpl_work(&z[jpl->off[JPL_EMB]], jpl->ncm[JPL_EMB], jpl->ncf[JPL_EMB], jpl->niv[JPL_EMB], t, jpl->inc, emb.u, emb.v, emb.w);
        jpl_work(&z[jpl->off[JPL_LUN]], jpl->ncm[JPL_LUN], jpl->ncf[JPL_LUN], jpl->niv[JPL_LUN], t, jpl->inc, lun.u, lun.v, lun.w);

        vecpos_set(pos->u, emb.u);
        vecpos_off(pos->u, lun.u, -1.0 / (1.0 + jpl->cem));

        vecpos_set(pos->v, emb.v);
        vecpos_off(pos->v, lun.v, -1.0 / (1.0 + jpl->cem));

        vecpos_set(pos->w, emb.w);
        vecpos_off(pos->w, lun.w, -1.0 / (1.0 + jpl->cem));
}

/* This was not fully tested */
static void _lun(struct _jpl_s *jpl, double *z, double t, struct mpos_s *pos)
{
        struct mpos_s emb, lun;

        jpl_work(&z[jpl->off[JPL_EMB]], jpl->ncm[JPL_EMB], jpl->ncf[JPL_EMB], jpl->niv[JPL_EMB], t, jpl->inc, emb.u, emb.v, emb.w);
        jpl_work(&z[jpl->off[JPL_LUN]], jpl->ncm[JPL_LUN], jpl->ncf[JPL_LUN], jpl->niv[JPL_LUN], t, jpl->inc, lun.u, lun.v, lun.w);

        vecpos_set(pos->u, emb.u);
        vecpos_off(pos->u, lun.u, jpl->cem / (1.0 + jpl->cem));

        vecpos_set(pos->v, emb.v);
        vecpos_off(pos->v, lun.v, jpl->cem / (1.0 + jpl->cem));

        vecpos_set(pos->w, emb.w);
        vecpos_off(pos->w, lun.w, jpl->cem / (1.0 + jpl->cem));

}


// function pointers are used to avoid a pointless switch statement
// Added _lun here (2020 Feb 26)
static void (* _help[_NUM_TEST])(struct _jpl_s *, double *, double, struct mpos_s *)
    = { _bar, _sun, _ear, _emb, _lun, _mer, _ven, _mar, _jup, _sat, _ura, _nep, _plu};

int jpl_calc(struct _jpl_s *pl, struct mpos_s *now, double jde, double rel, int n, int m)
{
        struct mpos_s pos;
        struct mpos_s ref;
        double t, *z;
        u_int32_t blk;
        int p;

        if (pl == NULL || pl->map == NULL || now == NULL)
                return -1;

        // check if covered by this file
        if (jde + rel < pl->beg || jde + rel > pl->end)
                return -1;

        // compute record number and 'offset' into record
        blk = (u_int32_t)((jde + rel - pl->beg) / pl->inc);
        z = pl->map + (blk + 2) * pl->rec;
        t = ((jde - pl->beg - (double)blk * pl->inc) + rel) / pl->inc;

        // the magick of function pointers
        _help[n](pl, z, t, &pos);
        _help[m](pl, z, t, &ref);

        for (p = 0; p < 3; p++) {
                now->u[p] = pos.u[p] - ref.u[p];
                now->v[p] = pos.v[p] - ref.v[p];
                now->w[p] = pos.w[p] - ref.w[p];
        }

        now->jde = jde + rel;
        return 0;
}

/*
 *  jpl_mass
 *
 */
double jpl_mass(struct _jpl_s *pl, int tar)
{
	char buf[8];
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

