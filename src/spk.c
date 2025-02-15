// spk.c - code to handle spice kernel position files

// https://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/req/daf.html
// https://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/req/spk.html
// https://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/req/naif_ids.html

#include <sys/types.h>
#include <sys/stat.h>
#include <sys/mman.h>
#include <math.h>
#include <ctype.h>
#include <errno.h>
#include <fcntl.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <time.h>
#include "assist.h"
#include "spk.h"


/*
 *  spk_free
 *
 */
int assist_spk_free(struct spk_s *pl)
{

	if (pl == NULL)
		return -1;

    if (pl->targets){
        for (int m = 0; m < pl->num; m++) {
            free(pl->targets[m].one);
            free(pl->targets[m].two);
        }
        free(pl->targets);
    }

	munmap(pl->map, pl->len);
	memset(pl, 0, sizeof(struct spk_s));
	free(pl);
	return 0;
}

/*
 *  assist_free_spk_constants
 *
 */
int assist_free_spk_constants(struct spk_global *sg) {

    if (sg->masses.names) {
        for (int i = 0; i < sg->masses.count; i++) {
            free(sg->masses.names[i]);
        }
        free(sg->masses.names);
    }

    free(sg->masses.values);
    free(sg);

    return ASSIST_SUCCESS;
}

void parse_comments(int fd, int first_summary_record, char **comments) {
    // Calculate the number of records in the comment section
    int num_records = first_summary_record - 2;  // Number of comment records
    if (num_records <= 0) {
        *comments = NULL;
        return;
    }

    // Allocate a buffer for the comments
    int comment_length = num_records * RECORD_LENGTH;
    char *buffer = malloc(comment_length + 1);  // +1 for null-terminator
    if (!buffer) {
        perror("malloc");
        exit(EXIT_FAILURE);
    }

    // Read the comment records
    ssize_t bytes_read;
    size_t total_bytes_read = 0;
    // Seek to beginning of file
    if (lseek(fd, 0, SEEK_SET) == -1) {
        perror("lseek");
        free(buffer);
        exit(EXIT_FAILURE);
    }
    // Read in each comment record
    for (int i = 0; i < num_records; i++) {
        if (lseek(fd, (1 + i) * RECORD_LENGTH, SEEK_SET) == -1) {
            perror("lseek");
            free(buffer);
            exit(EXIT_FAILURE);
        }
        bytes_read = read(fd, buffer + total_bytes_read, RECORD_LENGTH);
        if (bytes_read == -1) {
            perror("read");
            free(buffer);
            exit(EXIT_FAILURE);
        }

        // Remove all null character and end of comment markers from the end of the buffer
        // by deleting the bytes. These pad the ends and create extra newlines.
        while (buffer[total_bytes_read + bytes_read - 1] == '\0' || buffer[total_bytes_read + bytes_read - 1] == '\4') {
            bytes_read--;
        }

        total_bytes_read += bytes_read;
    }


    // DAF comments use the null character to indicate end of a line
    // replace with newline character
    for (char *p = buffer; p < buffer + total_bytes_read; p++) {
        if (*p == '\0') {
            *p = '\n';
        }
    }

    // Assign the buffer to the comments pointer
    *comments = buffer;


}

// sscanf doesn't know the 'D' scientific notation, so replace
// it with 'E' when surrounded by digit and sign characters
void replace_d_with_e(char *str) {
    char *d_char = strchr(str, 'D');
    if (d_char != NULL) {
        // Ensure it's a scientific notation before replacing
        if ((d_char > str && isdigit(*(d_char - 1))) && (isdigit(*(d_char + 1)) || (*(d_char + 1) == '+' || *(d_char + 1) == '-'))) {
            *d_char = 'E';
        }
    }
}


// Reads in the constants and masses from the comments of the DE440.bsp file
struct spk_global * assist_load_spk_constants(const char *path) {
    // Try opening file.
    int fd = open(path, O_RDONLY);
    if (fd < 0){
        return NULL;
    }

    // Load the file record
    union record_t * record = assist_load_spk_file_record(fd);

    char *comments = NULL;
    parse_comments(fd, record->file.fward, &comments);

    if (comments == NULL) {
        printf("No comments found.\n");
        return NULL;
    }

    struct spk_global* sg = malloc(sizeof(struct spk_global));
    if (!sg) {
        // Handle memory allocation failure
        close(fd);
        return NULL;
    }

    // Initialize the spk_global struct using designated initializers
    *sg = (struct spk_global){
        .con = {
            .AU = 0,
            .EMRAT = 0,
            .J2E = 0,
            .J3E = 0,
            .J4E = 0,
            .J2SUN = 0,
            .RE = 0,
            .CLIGHT = 0,
            .ASUN = 0
        },
        .masses = {
            .names = NULL,
            .values = NULL,
            .count = 0
        }
    };

    char *line = strtok(comments, "\n");
    char key[64];
    char value_str[64];
    double value;
    int in_constants_section = 0;

    while (line) {
        if (strstr(line, "Initial conditions and constants used for integration:")) {
            in_constants_section = 1;
        }

        if (in_constants_section) {
            replace_d_with_e(line);
            if (sscanf(line, "%63s %63s", key, value_str) == 2) {
                // Convert the value string to double
                value = strtod(value_str, NULL);
                if (strcmp(key, "cau") == 0) sg->con.AU = value;
                else if (strcmp(key, "EMRAT") == 0) sg->con.EMRAT = value;
                else if (strcmp(key, "J2E") == 0) sg->con.J2E = value;
                else if (strcmp(key, "J3E") == 0) sg->con.J3E = value;
                else if (strcmp(key, "J4E") == 0) sg->con.J4E = value;
                else if (strcmp(key, "J2SUN") == 0) sg->con.J2SUN = value;
                else if (strcmp(key, "AU") == 0) sg->con.AU = value;
                else if (strcmp(key, "RE") == 0) sg->con.RE = value;
                else if (strcmp(key, "CLIGHT") == 0) sg->con.CLIGHT = value;
                else if (strcmp(key, "ASUN") == 0) sg->con.ASUN = value;
                else if (strncmp(key, "GM", 2) == 0 || strncmp(key, "MA", 2) == 0) {
                    sg->masses.names = realloc(sg->masses.names, (sg->masses.count + 1) * sizeof(char *));
                    sg->masses.values = realloc(sg->masses.values, (sg->masses.count + 1) * sizeof(double));
                    sg->masses.names[sg->masses.count] = strdup(key);
                    sg->masses.values[sg->masses.count] = value;
                    sg->masses.count++;
                }
            }
        }

        line = strtok(NULL, "\n");
    }

    free(comments);
    close(fd);
    return sg;
}

/*
 *  assist_spk_init
 *
 */

// convert SPK epoch (since J2000.0) to julian day number
static double inline _jul(double eph)
	{ return 2451545.0 + eph / 86400.0; }

// Populate mass data for spk targets from the global masses array
void assist_spk_join_masses(struct spk_s *sp, struct spk_global *sg) {
    if (sp == NULL || sg == NULL) {
        return;
    }

    // Create an array based mapping of the GMX and target code formats
    struct {
        const char *name;
        int code;
    } planet_codes[] = {
        {"GMS", 10}, // Sun
        {"GM1", 1}, // Mercury
        {"GM2", 2}, // Venus
        {"GMB", 399}, // Earth
        {"GMB", 3}, // Earth-Moon Barycenter
        {"GMB", 301}, // Moon
        {"GM4", 4}, // Mars
        {"GM5", 5}, // Jupiter
        {"GM6", 6}, // Saturn
        {"GM7", 7}, // Uranus
        {"GM8", 8}, // Neptune
        {"GM9", 9} // Pluto
    };

    // Join the mass data by iterating through the targets
    for (int m = 0; m < sp->num; m++) {
        // Determine the label we are matching for in the masses array
        // If it is a planet, use the lookup value from planet_codes
        // If it is an asteroid, we format it as "MA" + code
        // Allocate the label and make sure it starts as a null string
        char *mass_label = calloc(64, sizeof(char));

        for (int i = 0; i < sizeof(planet_codes) / sizeof(planet_codes[0]); i++) {
            if (sp->targets[m].code == planet_codes[i].code) {
                strncpy(mass_label, planet_codes[i].name, 63);
                mass_label[63] = '\0'; // Ensure null termination
                break;
            }
        }
        // If mass label is still empty, it is an asteroid
        if (strlen(mass_label) == 0) {
            // Format the asteroid code to be MAdddd where dddd is 4 digit
            // 0 masked target code - 2000000
            sprintf(mass_label, "MA%04d", sp->targets[m].code - 2000000);
        }

        // Find the mass in the masses array
        for (int i = 0; i < sg->masses.count; i++) {
            if (strcmp(sg->masses.names[i], mass_label) == 0) {
                // Earth and moon mass is stored as one value
                // Use the ratio to split it using pl->con.EMRAT
                if (sp->targets[m].code == 399) {
                    sp->targets[m].mass = sg->masses.values[i] * (sg->con.EMRAT / (1. + sg->con.EMRAT));
                } else if (sp->targets[m].code == 301) {
                    sp->targets[m].mass = sg->masses.values[i] * (1./(1.+sg->con.EMRAT));
                } else {
                    sp->targets[m].mass = sg->masses.values[i];
                }
                break;
            }
        }

        if (sp->targets[m].mass == 0 && sp->targets[m].code != 199 && sp->targets[m].code != 299) {
            printf("Mass not found for target code: %d\n", sp->targets[m].code);
        }
    }
}

// Load the file record of an spk file
    // For file format information, see https://naif.jpl.nasa.gov/pub/naif/toolkit_docs/FORTRAN/req/daf.html
union record_t * assist_load_spk_file_record(int fd) {
    // Allocate memory for the record
    union record_t *record = calloc(1, sizeof(union record_t));
    if (!record) {
        fprintf(stderr, "Memory allocation failed.\n");
        close(fd);
        return NULL;
    }

    // Seek to the beginning of the file
    if (lseek(fd, 0, SEEK_SET) == -1) {
        perror("lseek");
        free(record);
        close(fd);
        return NULL;
    }

    // Read the file record
    ssize_t bytesRead = read(fd, record, RECORD_LENGTH);
    if (bytesRead != RECORD_LENGTH) {
        if (bytesRead == -1) {
            perror("read");
        } else {
            fprintf(stderr, "Incomplete read. Expected %d bytes, got %zd bytes.\n", RECORD_LENGTH, bytesRead);
        }
        free(record);
        close(fd);
		return NULL;
    }

    // Check if the file is a valid Double Precision Array File
    if (strncmp(record->file.locidw, "DAF/SPK", 7) != 0) {
        fprintf(stderr,"Error parsing DAF/SPK file. Incorrect header.\n");
        free(record);
		close(fd);
		return NULL;
	}

    // Check that the size of a summary record is equal to the size of our struct
    int nc = 8 * (record->file.nd + (record->file.ni + 1) / 2);
	if (nc != sizeof(struct sum_s)) {
        fprintf(stderr,"Error parsing DAF/SPK file. Wrong size of summary record.\n");
        free(record);
		close(fd);
		return NULL;
	}
    
    return record;
}


// Initialize the targets of a single spk file
// Note that target masses will not be populated until assist_spk_join_masses
// is called.
struct spk_s * assist_spk_init(const char *path) {
    // Try opening file.
    int fd = open(path, O_RDONLY);
    if (fd < 0){
        return NULL;
    }

    // Load the file record
    union record_t * record = assist_load_spk_file_record(fd);

    // Seek until the first summary record using the file record's fward pointer.
    // Record numbers start from 1 not 0 so we subtract 1 to get to the correct record.
    lseek(fd, (record->file.fward - 1) * RECORD_LENGTH, SEEK_SET);
    read(fd, record->buf, RECORD_LENGTH);

	// We are at the first summary block, validate
    if ((int64_t)record->buf[8] != 0) {
        fprintf(stderr, "Error parsing DAF/SPL file. Cannot find summary block.\n");
		close(fd);
        return NULL; 
	}

	// okay, let's go
	struct spk_s* pl = calloc(1, sizeof(struct spk_s));
    while (1){ // Loop over records 
        for (int b = 0; b < (int)record->summary.nsum; b++) { // Loop over summaries
            struct sum_s* sum = &record->summary.s[b]; // get current summary
            
            // Index in our arrays for current target
            int m = -1;

            // Check to see if we are adding to an existing target
            for (int i = 0; i < pl->num; i++) {
                if (pl->targets[i].code == sum->tar) {
                    m = i;
                    break;
                }
            }

            // If this is a new target, add it to the list
            if (m == -1) {
                m = pl->num;
                if (pl->num <= pl->allocated_num){
                    pl->allocated_num += 32; // increase space in batches of 32
                    pl->targets = realloc(pl->targets, pl->allocated_num*sizeof(struct spk_target));
                }
                pl->num++;
                pl->targets[m].code = sum->tar;
                pl->targets[m].cen = sum->cen;
                pl->targets[m].beg = _jul(sum->beg);
                pl->targets[m].res = _jul(sum->end) - pl->targets[m].beg;
                pl->targets[m].one = calloc(32768, sizeof(int));
                pl->targets[m].two = calloc(32768, sizeof(int));
                pl->targets[m].ind = -1;
                // Set default of mass to 0
                pl->targets[m].mass = 0;
            }

            // add index for target
            pl->targets[m].ind++;
            pl->targets[m].one[pl->targets[m].ind] = sum->one;
            pl->targets[m].two[pl->targets[m].ind] = sum->two;
            pl->targets[m].end = _jul(sum->end);
        }

        // Location of next record
        long n = (long)record->summary.next - 1;
        if (n<0){
            // this is already the last record.
            break;
        }else{
            // Find and read next record
            lseek(fd, n * RECORD_LENGTH, SEEK_SET);
            read(fd, record->buf, RECORD_LENGTH);
        }
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
    struct spk_target* target = &(pl->targets[m]);
    if (jde + rel < target->beg || jde + rel > target->end){
        return ASSIST_ERROR_COVERAGE;
    }

    *GM = target->mass; // Note mass constants defined in DE440/441 ephemeris files. If not found mass of 0 is used.

	int n, b, p, P, R;
	double T[32];
    // double S[32]; // Not used at the moment
	double *val, z;
    struct mpos_s pos = {0};

	for (n = 0; n < 3; n++)
		pos.u[n] = pos.v[n] = 0.0;

	// find location of 'directory' describing the data records
	n = (int)((jde + rel - target->beg) / target->res);
	val = (double *)pl->map + target->two[n] - 1;	

	// record size and number of coefficients per coordinate
	R = (int)val[-1];
	P = (R - 2) / 3; // must be < 32 !!

	// pick out the precise record
	b = (int)(((jde - _jul(val[-3])) + rel) / (val[-2] / 86400.0));
	//+ sizeof(double) * b * R;
	val = (double *)pl->map + (target->one[n] - 1) + b * R;

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

struct spk_target* assist_spk_find_target(struct spk_s *pl, int code) {
    for (int i = 0; i < pl->num; i++) {
        if (pl->targets[i].code == code) {
            return &(pl->targets[i]);
        }
    }
    return NULL;
}

struct mpos_s assist_spk_target_pos(struct spk_s *pl, struct spk_target* target, double jde, double rel)
{
    int n, b, p, P, R;
    double T[32];
    double S[32];
    double U[32];
    double *val, z;
    struct mpos_s pos = {0};

    // find location of 'directory' describing the data records
    n = (int)((jde + rel - target->beg) / target->res);
    val = (double *)pl->map + target->two[n] - 1;

    // record size and number of coefficients per coordinate
    R = (int)val[-1];
    P = (R - 2) / 3; // must be < 32 !!

    // pick out the precise record
    b = (int)(((jde - _jul(val[-3])) + rel) / (val[-2] / 86400.0));
    val = (double *)pl->map + (target->one[n] - 1) + b * R;

    // scale to interpolation units
    z = ((jde - _jul(val[0])) + rel) / (val[1] / 86400.0);

    // Calculate the scaling factor 'c'
    double c = 1.0 / val[1];

    // set up Chebyshev polynomials
    T[0] = 1.0; T[1] = z;
    S[0] = 0.0; S[1] = 1.0;
    U[0] = 0.0; U[1] = 0.0; U[2] = 4.0;

    for (p = 2; p < P; p++) {
		T[p] = 2.0 * z * T[p-1] - T[p-2];
        S[p] = 2.0 * z * S[p-1] + 2.0 * T[p-1] - S[p-2];
	}
    for (p = 3; p < P; p++) {
        U[p] = 2.0 * z * U[p-1] + 4.0 * S[p-1] - U[p-2];
    }


    for (n = 0; n < 3; n++) {
        b = 2 + n * P;
        pos.u[n] = pos.v[n] = pos.w[n] = 0.0;

        // sum interpolation stuff
        for (p = 0; p < P; p++) {
            double coeff = val[b + p];
            pos.u[n] += coeff * T[p];
            pos.v[n] += coeff * S[p] * c;
            pos.w[n] += coeff * U[p] * c * c;
        }
    }

    return pos;
}


// Calculate the position and velocity of planets from DE440 SPK file
enum ASSIST_STATUS assist_spk_calc_planets(struct assist_ephem* ephem, double jde, double rel, int code, double* GM, double* out_x, double* out_y, double* out_z, double* out_vx, double* out_vy, double* out_vz, double* out_ax, double* out_ay, double* out_az)
{

    struct spk_s* pl = ephem->spk_planets;


    struct spk_target* target = assist_spk_find_target(pl, code);

    if (target == NULL) {
        return(ASSIST_ERROR_NEPHEM);
    }

    if (jde + rel < target->beg || jde + rel > target->end){
        return ASSIST_ERROR_COVERAGE;
    }

    *GM = target->mass; // Note mass constants defined in DE440/441 ephemeris files. If not found mass of 0 is used.

    struct mpos_s pos = assist_spk_target_pos(pl, target, jde, rel);

    // Earth and Moon must be translated from EMB to SSB frame
    if (code == 301 || code == 399) {
        struct spk_target* emb = assist_spk_find_target(pl, 3);
        struct mpos_s emb_pos = assist_spk_target_pos(pl, emb, jde, rel);

        for (int i = 0; i < 3; i++) {
            pos.u[i] += emb_pos.u[i];
            pos.v[i] += emb_pos.v[i];
            pos.w[i] += emb_pos.w[i];
        }

    }

    // Convert to AU and AU/day
    const double au = assist_get_constant(ephem, "AU");
    const double seconds_per_day = 86400.;

    for (int i = 0; i < 3; i++) {
        pos.u[i] /= au;
        pos.v[i] /= au / seconds_per_day;
        pos.w[i] /= au / (seconds_per_day * seconds_per_day);
    }

    *out_x = pos.u[0];
    *out_y = pos.u[1];
    *out_z = pos.u[2];
    *out_vx = pos.v[0];
    *out_vy = pos.v[1];
    *out_vz = pos.v[2];
    *out_ax = pos.w[0];
    *out_ay = pos.w[1];
    *out_az = pos.w[2];

    return ASSIST_SUCCESS;
}
