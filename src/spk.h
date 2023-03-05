
// spk.h

// https://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/req/daf.html
// https://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/req/spk.html
// https://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/req/naif_ids.html

#ifndef _SPK_H
#define _SPK_H

#include "assist.h"

struct mpos_s {
	double u[3];
	double v[3];
	double w[3];  
	double jde;
};


struct spk_target {
    int code;       // Target code
    int cen;     // Centre target
    double mass;    // Mass. Set to 0 if not found in ephemeris file.
    double beg;     // Begin epoch
    double end;     // End epoch
    double res;     // Epoch step
	int *one;		// Record index
	int *two;		// ... ditto
	int ind;		// Length of index

};

struct spk_s {
    struct spk_target* targets;
	int num;			// number of targets
	int allocated_num;	// space allocated for this many targets
	void *map;			// memory map
	size_t len;			// map length
};


int assist_spk_free(struct spk_s *pl);
struct spk_s * assist_spk_init(const char *path);
enum ASSIST_STATUS assist_spk_calc(struct spk_s *pl, double jde, double rel, int m, double* GM, double* x, double* y, double* z);

#endif // _SPK_H

