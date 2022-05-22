
// spk.h

// https://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/req/daf.html
// https://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/req/spk.html
// https://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/req/naif_ids.html

#ifndef _SPK_H
#define _SPK_H

#define _SPK_MAX	32	// maximum body count

enum {
	SPK_NAIF_SSB		= 0,
	SPK_NAIF_MER		= 1,
	SPK_NAIF_VEN		= 2,
	SPK_NAIF_EMB		= 3,
	SPK_NAIF_MAR		= 4,
	SPK_NAIF_JUP		= 5,
	SPK_NAIF_SAT		= 6,
	SPK_NAIF_URA		= 7,
	SPK_NAIF_NEP		= 8,
	SPK_NAIF_PLU		= 9,
	SPK_NAIF_SUN		= 10,
	SPK_NAIF_MER___		= 199,		// ?
	SPK_NAIF_VEN___		= 299,		// ?
	SPK_NAIF_EAR___		= 399,		// ?
	SPK_NAIF_MOON		= 301,
	SPK_NAIF_CAMLLA         = 2000107,      //
	SPK_NAIF_CERES		= 2000001,	//
	SPK_NAIF_CYBELE		= 2000065,	//
	SPK_NAIF_DAVIDA		= 2000511,	//
	SPK_NAIF_EUNOMIA	= 2000015,	//
	SPK_NAIF_EUPHROSYNE	= 2000031,	//
	SPK_NAIF_EUROPA		= 2000052,	//
	SPK_NAIF_HYGIEA		= 2000010,	//
	SPK_NAIF_INTERAMNIA	= 2000704,	//
	SPK_NAIF_IRIS           = 2000007,      //
	SPK_NAIF_JUNO		= 2000003,	//
	SPK_NAIF_PALLAS		= 2000002,	//
	SPK_NAIF_PSYCHE		= 2000016,	//
	SPK_NAIF_SYLVIA		= 2000087,	//
	SPK_NAIF_THISBE		= 2000088,	//
	SPK_NAIF_VESTA		= 2000004,	//

};

struct mpos_s {
	double u[3];
	double v[3];
	double w[3];  
	double jde;
};

struct spk_s {

	int tar[_SPK_MAX];		// target code
	int cen[_SPK_MAX];		// centre target
	double beg[_SPK_MAX];		// begin epoch
	double res[_SPK_MAX];		// epoch step
	int *one[_SPK_MAX];		// record index
	int *two[_SPK_MAX];		// ... ditto
	int ind[_SPK_MAX];		// length of index

	int num;			// number of targets
	void *map;			// memory map
	size_t len;			// map length
};


int spk_free(struct spk_s *pl);
struct spk_s * spk_init(const char *path);
int spk_find(struct spk_s *pl, int m);
int spk_calc(struct spk_s *pl, int tar, double jde, struct mpos_s *pos);

#endif // _SPK_H

