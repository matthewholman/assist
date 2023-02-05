
// spk.c - code to handle spice kernel position files

// https://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/req/daf.html
// https://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/req/spk.html
// https://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/req/naif_ids.html

#include <sys/types.h>
#include <sys/stat.h>
#include <sys/mman.h>
#include <math.h>
#include <errno.h>
#include <fcntl.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <time.h>
#include "assist.h"
#include "spk.h"
#include "const.h"


struct sum_s {
	double beg;		// begin epoch, seconds since J2000.0
	double end;		// ending epoch
	int tar;		// target code
	int cen;		// centre code (10 = sun)
	int ref;		// reference frame (1 = J2000.0)
	int ver;		// type of ephemeris (2 = chebyshev)
	int one;		// initial array address
	int two;		// final array address
};


/*
 *  spk_free
 *
 */
int assist_spk_free(struct spk_s *pl)
{
	int m;

	if (pl == NULL)
		return -1;

	for (m = 0; m < pl->num; m++) {
		free(pl->one[m]);
		free(pl->two[m]);
	}

	munmap(pl->map, pl->len);
	memset(pl, 0, sizeof(struct spk_s));
	free(pl);
	return 0;
}


/*
 *  assist_spk_init
 *
 */

// convert SPK epoch (since J2000.0) to julian day number
static double inline _jul(double eph)
	{ return 2451545.0 + eph / 86400.0; }

// check for any non-7bit ascii characters
static int _com(const char *buf)
{
	for (int n = 0; n < 1024; n++)
		{ if (buf[n] < 0) return 0; }

	return 1;
}

// display output strings to console
//static void _sho(const char *buf)
//{
//	int n;
//
//	// this doesn't always handle the newlines correctly
//	for (n = 0; buf[n+1] != '\0';)
//		n += fprintf(stdout, "%s\n", &buf[n]) - 1;
//}

struct spk_s * assist_spk_init(const char *path)
{
	struct spk_s *pl;
	struct stat sb;
	char buf[1024];
	struct sum_s *sum;
	double *val;
	int fd, nd, ni, nc;
	int m, n, c, b, B;
	off_t off;

	if ((fd = open(path, O_RDONLY)) < 0)
		return NULL;

	pl = malloc(sizeof(struct spk_s));
	memset(pl, 0, sizeof(struct spk_s));
	val = (double *)buf;
	sum = (struct sum_s *)buf;

	if (fstat(fd, &sb) < 0)
		goto err;

	// LOCIDW
	if (read(fd, buf, 8) != 8)
		goto err;

	if (strncmp(buf, "DAF/SPK", 7) != 0) {
		errno = EILSEQ;
		goto err;
	}

	// ND, NI
	read(fd, &nd, sizeof(int));
	read(fd, &ni, sizeof(int));

	// length of each segment, must match our sum_s struct
	nc = 8 * ( nd + (ni + 1) / 2 );

	if (nc != sizeof(struct sum_s)) {
		errno = EILSEQ;
		goto err;
	}

	// could check the other headers, but we really don't care
	// so find the first summary record (after potential comments)
	off = lseek(fd, 1024, SEEK_SET);
	read(fd, buf, 1024);

	while (_com(buf) > 0) {
	//	_sho(buf);
		off = lseek(fd, 0, SEEK_CUR);
		read(fd, buf, 1024);
	}

	// we are at the first summary block, validate
	if (val[1] != 0.0) {
		errno = EILSEQ;
		goto err;
	}

	// okay, let's go
	m = 0;
next:	n = (int)val[0] - 1;
	B = (int)val[2];

	for (b = 0; b < B; b++) {
		sum = (struct sum_s *)&buf[24 + b * sizeof(struct sum_s)];

//		fprintf(stdout, "beg %.1f end %.1f tar %d cen %d ref %d ver %d one %d two %d\n",
//				_jul(sum->beg), _jul(sum->end), sum->tar, sum->cen,
//				sum->ref, sum->ver, sum->one, sum->two);

		// pick out new target!
		if (sum->tar != pl->tar[m]) {
			m = pl->num++;
			pl->tar[m] = sum->tar;
			pl->cen[m] = sum->cen;
			pl->beg[m] = _jul(sum->beg);
			pl->res[m] = _jul(sum->end) - pl->beg[m];
			pl->one[m] = calloc(32768, sizeof(int));
			pl->two[m] = calloc(32768, sizeof(int));
		}

		// add index
		c = pl->ind[m]++;
		pl->one[m][c] = sum->one;
		pl->two[m][c] = sum->two;
	}

	if (n >= 0) {
		// this could probably be more elegant if the mmap happened sooner
		off = lseek(fd, n * 1024, SEEK_SET);
		read(fd, buf, 1024);
		goto next;
	}

	// memory map : kernel caching and thread safe
	pl->len = sb.st_size;
	pl->map = mmap(NULL, pl->len, PROT_READ, MAP_SHARED, fd, 0);

	if (pl->map == NULL)
		goto err;

	if (close(fd) < 0)
		{ ; }
#if defined(MADV_RANDOM)
	if (madvise(pl->map, pl->len, MADV_RANDOM) < 0)
		{ ; }
#endif

	return pl;

err:	perror(path);
	free(pl);

	return NULL;
}


/*
 *  assist_spk_find
 *
 *  See if we have the given body.
 *
 */
int assist_spk_find(struct spk_s *pl, int tar)
{
	int n;

	if (pl == NULL)
		return -1;

	for (n = 0; n < pl->num; n++)
		if (pl->tar[n] == tar)
			{ return n; }

	return -1;
}


/*
 *  assist_spk_calc
 *
 *  Compute the position and velocity after fetching the chebyshev polynomials.
 *
 */

enum ASSIST_STATUS assist_spk_calc(struct spk_s *pl, double jde, double rel, int m, double* GM, double* out_x, double* out_y, double* out_z)
{
    const static double JPL_GM[16] =    
    {
	JPL_EPHEM_MA0107, // 107 camilla
	JPL_EPHEM_MA0001, // 1 Ceres
	JPL_EPHEM_MA0065, // 65 cybele
	JPL_EPHEM_MA0511, // 511 davida
	JPL_EPHEM_MA0015, // 15 eunomia
	JPL_EPHEM_MA0031, // 31 euphrosyne	    
	JPL_EPHEM_MA0052, // 52 europa
	JPL_EPHEM_MA0010, // 10 hygiea
	JPL_EPHEM_MA0704, // 704 interamnia
	JPL_EPHEM_MA0007, // 7 iris
	JPL_EPHEM_MA0003, // 3 juno
	JPL_EPHEM_MA0002, // 2 pallas
	JPL_EPHEM_MA0016, // 16 psyche
	JPL_EPHEM_MA0087, // 87 sylvia
	JPL_EPHEM_MA0088, // 88 thisbe
	JPL_EPHEM_MA0004  // 4 vesta
    };

    if(pl==NULL){
        return(ASSIST_ERROR_AST_FILE);	
    }

    if(m<0 || m > pl->num){
        return(ASSIST_ERROR_NAST);
    }

    // TODO: again, the units might be handled more
    // generally

    *GM = JPL_GM[m];

	int n, b, p, P, R;
	double T[32];
    // double S[32]; // Not used at the moment
	double *val, z;
    struct mpos_s pos = {0};

	for (n = 0; n < 3; n++)
		pos.u[n] = pos.v[n] = 0.0;

	// find location of 'directory' describing the data records
	n = (int)((jde + rel - pl->beg[m]) / pl->res[m]);
	//val = (double*)pl->map + sizeof(double) * (pl->two[m][n] - 1);
	val = (double *)pl->map + pl->two[m][n] - 1;	

	// record size and number of coefficients per coordinate
	R = (int)val[-1];
	P = (R - 2) / 3; // must be < 32 !!

	// pick out the precise record
	b = (int)(((jde - _jul(val[-3])) + rel) / (val[-2] / 86400.0));
	//val = (double*)pl->map + sizeof(double) * (pl->one[m][n] - 1)
	//+ sizeof(double) * b * R;
	val = (double *)pl->map + (pl->one[m][n] - 1) + b * R;

	// scale to interpolation units
	z = ((jde - _jul(val[0])) + rel) / (val[1] / 86400.0);

	// set up Chebyshev polynomials
	T[0] = 1.0; T[1] = z;   
    // Not used at the moment:
    // S[1] = 1.0; S[0] = 0.0;

	for (p = 2; p < P; p++) {
		T[p] = 2.0 * z * T[p-1] - T[p-2];
        // Not used at the moment:
		// S[p] = 2.0 * z * S[p-1] + 2.0 * T[p-1] - S[p-2];
	}

	for (n = 0; n < 3; n++) {
		b = 2 + n * P;

		// sum interpolation stuff
		for (p = 0; p < P; p++) {
		    pos.u[n] += val[b + p] * T[p];
            // Not used at the moment:
		    //pos.v[n] += val[b + p] * S[p];
			
		}

		// restore units to [AU] and [AU/day]
		pos.u[n] /= 149597870.7;
        // Not used at the moment:
		//pos.v[n] /= 149597870.7 / 86400.0;  
		//pos.v[n] /= val[1];
	}
    
    *out_x = pos.u[0];
    *out_y = pos.u[1];
    *out_z = pos.u[2];

	return ASSIST_SUCCESS;
}

