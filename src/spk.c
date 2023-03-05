
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


#define record_length 1024
// check for any non-7bit ascii characters
static int _com(const char *record) {
	for (int n = 0; n < record_length; n++)
		{ if (record[n] < 0) return 0; }

	return 1;
}

struct spk_s * assist_spk_init(const char *path) {
    // For file format information, see https://naif.jpl.nasa.gov/pub/naif/toolkit_docs/FORTRAN/req/daf.html

    // Format for one summary
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

    // File is split into records. We read one record at a time.
    union {
	    char buf[record_length];
        struct {
            double next;    // The record number of the next summary record in the file. Zero if this is the final summary record.
            double prev;    // The record number of the previous summary record in the file. Zero if this is the initial summary record.
            double nsum;    // Number of summaries in this record
            struct sum_s s[25]; // Summaries (25 is the maximum)
        } summary;          // Summary record
        struct {
            char locidw[8]; // An identification word
            int nd;         // The number of double precision components in each array summary.
            int ni;         // The number of integer components in each array summary.
        } file;             // File record
    } record;

    // Try opening file.
	int fd = open(path, O_RDONLY);
	if (fd < 0){
		return NULL;
    }

	// Read the file record. 
    read(fd, &record, 1024);
    // Check if the file is a valid Double Precision Array File
	if (strncmp(record.file.locidw, "DAF/SPK", 7) != 0) {
        fprintf(stderr,"Error parsing DAF/SPK file. Incorrect header.\n");
		close(fd);
		return NULL;
	}

    // Check that the size of a summary record is equal to the size of our struct.
	int nc = 8 * ( record.file.nd + (record.file.ni + 1) / 2 );
	if (nc != sizeof(struct sum_s)) {
        fprintf(stderr,"Error parsing DAF/SPK file. Wrong size of summary record.\n");
		close(fd);
		return NULL;
	}
    
    // Continue reading file until we find a non-ascii record.
    do {
		read(fd, record.buf, 1024);
    } while (_com(record.buf) > 0);

	// We are at the first summary block, validate
	if ((int64_t)record.buf[8] != 0) {
        fprintf(stderr, "Error parsing DAF/SPL file. Cannot find summary block.\n");
		close(fd);
        return NULL; 
	}

	// okay, let's go
	struct spk_s* pl = calloc(1, sizeof(struct spk_s));
	int m = 0; 
    while (1){ // Loop over records 
        for (int b = 0; b < (int)record.summary.nsum; b++) { // Loop over summaries
            struct sum_s* sum = &record.summary.s[b]; // get current summary

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
            pl->one[m][pl->ind[m]] = sum->one;
            pl->two[m][pl->ind[m]] = sum->two;
            pl->end[m] = _jul(sum->end);
            pl->ind[m]++;
        }

        // Location of next record
        int n = (int)record.summary.next - 1;
        if (n<0){ // this is already the last record.
            break;
        }
        // Find and read next record
        lseek(fd, n * 1024, SEEK_SET);
        read(fd, record.buf, 1024);
    }

    // Get file size
	struct stat sb;
	if (fstat(fd, &sb) < 0){
        fprintf(stderr, "Error calculating size for DAF/SPL file.\n");
        return NULL; 
    }
	pl->len = sb.st_size;

    // Memory map
	pl->map = mmap(NULL, pl->len, PROT_READ, MAP_SHARED, fd, 0);
	if (pl->map == NULL){
        fprintf(stderr, "Error creating memory map.\n");
        return NULL; // Will leak memory
    }

#if defined(MADV_RANDOM)
	if (madvise(pl->map, pl->len, MADV_RANDOM) < 0){
        fprintf(stderr, "Error while calling madvise().\n");
        return NULL; // Will leak memory
    }
#endif

	close(fd);
	return pl;
}


enum ASSIST_STATUS assist_spk_calc(struct spk_s *pl, double jde, double rel, int m, double* GM, double* out_x, double* out_y, double* out_z)
{
    if(pl==NULL){
        return(ASSIST_ERROR_AST_FILE);	
    }

    if(m<0 || m > pl->num){
        return(ASSIST_ERROR_NAST);
    }
        
    if (jde + rel < pl->beg[m] || jde + rel > pl->end[m]){
        return ASSIST_ERROR_COVERAGE;
    }

    *GM = pl->mass[m]; // Note mass constants defined in DE440/441 ephemeris files. If not found mass of 0 is used.

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

