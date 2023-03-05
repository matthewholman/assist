
// spk.h

// https://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/req/daf.html
// https://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/req/spk.html
// https://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/req/naif_ids.html

#ifndef _SPK_H
#define _SPK_H

#include "assist.h"
#define _SPK_MAX	32	// maximum body count

struct mpos_s {
	double u[3];
	double v[3];
	double w[3];  
	double jde;
};

struct spk_s {

	int tar[_SPK_MAX];		// target code
	double mass[_SPK_MAX];		// mass (0 unless id is found in jpl file)
	int cen[_SPK_MAX];		// centre target
	double beg[_SPK_MAX];		// begin epoch
	double end[_SPK_MAX];		// begin epoch
	double res[_SPK_MAX];		// epoch step
	int *one[_SPK_MAX];		// record index
	int *two[_SPK_MAX];		// ... ditto
	int ind[_SPK_MAX];		// length of index

	int num;			// number of targets
	void *map;			// memory map
	size_t len;			// map length
};


int assist_spk_free(struct spk_s *pl);
struct spk_s * assist_spk_init(const char *path);
int assist_spk_find(struct spk_s *pl, int m);
enum ASSIST_STATUS assist_spk_calc(struct spk_s *pl, double jde, double rel, int m, double* GM, double* x, double* y, double* z);

#endif // _SPK_H

