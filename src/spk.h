
// spk.h

// https://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/req/daf.html
// https://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/req/spk.html
// https://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/req/naif_ids.html

#ifndef _SPK_H
#define _SPK_H

#include "assist.h"

#define RECORD_LENGTH 1024
#define MAX_COMMENT_LENGTH 1000000

struct mpos_s {
	double u[3];
	double v[3];
	double w[3];  
	double jde;
};

struct spk_constants {
    double AU;                     // definition of AU
    double EMRAT;                     // Earth/Moon mass ratio
    double J2E;                     // Other constant names follow JPL
    double J3E;
    double J4E;
    double J2SUN;
    double RE;
    double CLIGHT;
    double ASUN;
};

struct mass_data {
    char **names;
    double *values;
    int count;
};

// Represents global values of constants and masses
// fetched from the comments of an SPK file (used with DE440)
struct spk_global {
    struct spk_constants con;
    struct mass_data masses;

};

// Represents available data for a target in a SPK file
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

// Represents a collection of targets in an SPK file
struct spk_s {
    struct spk_target* targets;
	int num;			// number of targets
	int allocated_num;	// space allocated for this many targets
	void *map;			// memory map
	size_t len;			// map length
};

    // Format for one summary
struct sum_s {
    double beg;        // begin epoch, seconds since J2000.0
    double end;        // ending epoch
    int tar;        // target code
    int cen;        // centre code (10 = sun)
    int ref;        // reference frame (1 = J2000.0)
    int ver;        // type of ephemeris (2 = chebyshev)
    int one;        // initial array address
    int two;        // final array address
};

// File is split into records. We read one record at a time.
union record_t {
    char buf[RECORD_LENGTH];
    struct {
        double next;    // The record number of the next summary record in the file. Zero if this is the final summary record.
        double prev;    // The record number of the previous summary record in the file. Zero if this is the initial summary record.
        double nsum;    // Number of summaries in this record
        struct sum_s s[25]; // Summaries (25 is the maximum)
    } summary;          // Summary record
    // See: https://naif.jpl.nasa.gov/pub/naif/toolkit_docs/FORTRAN/req/daf.html#The%20File%20Record
    struct {
        char locidw[8]; // An identification word
        int nd;         // The number of double precision components in each array summary.
        int ni;         // The number of integer components in each array summary.
        char locifn[60];// The internal name or description of the array file.
        int fward;      // The record number of the initial summary record in the file.
        int bward;      // The record number of the final summary record in the file.
        int free;        // Next available free record
    } file;             // File record
};





int assist_spk_free(struct spk_s *pl);
int assist_free_spk_constants(struct spk_global *sg);
struct spk_s * assist_spk_init(const char *path);
union record_t * assist_load_spk_file_record(int fd);
struct spk_global * assist_load_spk_constants(const char *path);
void parse_comments(int fd, int first_summary_record, char **comments);
void assist_spk_join_masses(struct spk_s *sp, struct spk_global *sg);
enum ASSIST_STATUS assist_spk_calc(struct spk_s *pl, double jde, double rel, int m, double* GM, double* x, double* y, double* z);
enum ASSIST_STATUS assist_spk_calc_planets(struct assist_ephem* ephem, double jde, double rel, int m, double* GM, double* x, double* y, double* z, double* vx, double* vy, double* vz, double* ax, double* ay, double* az);

#endif // _SPK_H

