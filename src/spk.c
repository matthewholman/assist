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
#include <stdint.h>
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
    for (int i = 0; str[i] != '\0'; i++) {
        if (str[i] == 'D' || str[i] == 'd') {
            str[i] = 'e';
        }
    }
}






/*
 *  assist_spk_init
 *
 */

// convert SPK epoch (since J2000.0) to julian day number
static double inline _jul(double eph)
	{ return 2451545.0 + eph / 86400.0; }

// Populate mass data for spk targets from mass data structure
void assist_spk_join_masses(struct spk_s *sp, const struct mass_data* masses, double emrat) {
    if (sp == NULL || masses == NULL || masses->names == NULL) {
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
        // Skip if mass is already set
        if (sp->targets[m].mass != 0) {
            continue;
        }

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

        // Find the mass in the mass arrays
        for (int i = 0; i < masses->count; i++) {
            if (strcmp(masses->names[i], mass_label) == 0) {
                // Earth and moon mass is stored as one value
                // Use the ratio to split it using emrat
                if (sp->targets[m].code == 399) {
                    sp->targets[m].mass = masses->values[i] * (emrat / (1. + emrat));
                } else if (sp->targets[m].code == 301) {
                    sp->targets[m].mass = masses->values[i] * (1./(1.+emrat));
                } else {
                    sp->targets[m].mass = masses->values[i];
                }
                break;
            }
        }

        if (sp->targets[m].mass == 0 && sp->targets[m].code != 199 && sp->targets[m].code != 299) {
            printf("Mass not found for target code: %d\n", sp->targets[m].code);
        }
        
        free(mass_label);
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

// Detect legacy ephemeris file formats and provide helpful error messages
// (enum ephemeris_file_format_t is now defined in spk.h)

// Simple signature check for legacy binary format
int assist_detect_legacy_binary_signature(int fd) {
    // Legacy binary files have constant names at offset 0x00FC (252 bytes)
    // Check for typical ephemeris constant names like "GMS   ", "GM1   ", etc.
    
    char const_names[6 * 3];  // Read first 3 constant names (6 chars each)
    if (lseek(fd, 0x00FC, SEEK_SET) != 0x00FC) {
        return 0;  // Can't seek to expected position
    }
    
    if (read(fd, const_names, sizeof(const_names)) != sizeof(const_names)) {
        return 0;  // Can't read constant names
    }
    
    // Check if names look like typical ephemeris constants
    // They should be printable ASCII with trailing spaces, like "GMS   ", "GM1   "
    int valid_names = 0;
    for (int i = 0; i < 3; i++) {
        char* name = &const_names[i * 6];
        
        // Check if it looks like a constant name: starts with letters, may have trailing spaces
        if ((name[0] >= 'A' && name[0] <= 'Z') || (name[0] >= 'a' && name[0] <= 'z')) {
            int has_letters = 0;
            int has_spaces = 0;
            for (int j = 0; j < 6; j++) {
                if ((name[j] >= 'A' && name[j] <= 'Z') || (name[j] >= 'a' && name[j] <= 'z') || 
                    (name[j] >= '0' && name[j] <= '9')) {
                    has_letters++;
                } else if (name[j] == ' ') {
                    has_spaces++;
                } else if (name[j] == '\0') {
                    // Null termination is ok
                    break;
                }
            }
            
            // Valid constant name should have at least 2 letters and some spaces for padding
            if (has_letters >= 2 && (has_spaces > 0 || strlen(name) < 6)) {
                valid_names++;
            }
        }
    }
    
    // If at least 2 out of 3 names look like ephemeris constants, it's likely a legacy file
    return (valid_names >= 2);
}

ephemeris_file_format_t assist_detect_ephemeris_file_format(int fd) {
    char buf[1024];
    ssize_t bytes_read;
    
    // Read first chunk of file
    lseek(fd, 0, SEEK_SET);
    bytes_read = read(fd, buf, sizeof(buf));
    if (bytes_read <= 0) {
        return FILE_FORMAT_UNKNOWN;
    }
    
    // Check for valid .bsp SPK file (DAF/SPK header) first
    if (bytes_read >= 8 && strncmp(buf, "DAF/SPK ", 8) == 0) {
        return FILE_FORMAT_VALID_BSP;
    }
    
    // Check for legacy binary format signature
    if (assist_detect_legacy_binary_signature(fd)) {
        // Reset file position after signature check
        lseek(fd, 0, SEEK_SET);
        return FILE_FORMAT_BINARY_LEGACY;
    }
    
    // Check for ASCII source file markers
    // These files typically contain header information about the ephemeris
    // and start with comment blocks
    if (strstr(buf, "JPL") && strstr(buf, "EPHEMERIS") && 
        (strstr(buf, "ASCII") || strstr(buf, "ascii"))) {
        return FILE_FORMAT_ASCII_SOURCE;
    }
    
    return FILE_FORMAT_UNKNOWN;
}

static void print_legacy_format_error(ephemeris_file_format_t format, const char* filepath) {
    if (format == FILE_FORMAT_ASCII_SOURCE) {
        fprintf(stderr, "(ASSIST) Error: ASCII ephemeris source file detected. ASSIST requires binary SPK (.bsp) format.\n");
    } else if (format == FILE_FORMAT_BINARY_LEGACY) {
        fprintf(stderr, "(ASSIST) Error: Legacy binary ephemeris file detected. ASSIST requires binary SPK (.bsp) format.\n");
    } else {
        fprintf(stderr, "(ASSIST) Error: Unsupported ephemeris file format. ASSIST requires binary SPK (.bsp) format.\n");
    }
    
    fprintf(stderr, "(ASSIST) Download from: https://naif.jpl.nasa.gov/pub/naif/generic_kernels/spk/planets/de440.bsp\n");
    
    // Check if we're running in Python environment
    if (getenv("PYTHONPATH") || getenv("VIRTUAL_ENV") || getenv("CONDA_DEFAULT_ENV")) {
        fprintf(stderr, "(ASSIST) Or install via: pip install naif-de440 jpl-small-bodies-de441-n16\n");
    }
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

    // Check for legacy ephemeris file formats
    ephemeris_file_format_t format = assist_detect_ephemeris_file_format(fd);
    if (format != FILE_FORMAT_VALID_BSP) {
        print_legacy_format_error(format, path);
        close(fd);
        return NULL;
    }

    // Load the file record
    union record_t * record = assist_load_spk_file_record(fd);
    if (!record) {
        close(fd);
        return NULL;
    }

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


enum ASSIST_STATUS assist_spk_calc(const struct spk_s *pl, double jde, double rel, int m, double* GM, double* out_x, double* out_y, double* out_z)
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

struct spk_target* assist_spk_find_target(const struct spk_s *pl, int code) {
    for (int i = 0; i < pl->num; i++) {
        if (pl->targets[i].code == code) {
            return &(pl->targets[i]);
        }
    }
    return NULL;
}

struct mpos_s assist_spk_target_pos(const struct spk_s *pl, const struct spk_target* target, double jde, double rel)
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
enum ASSIST_STATUS assist_spk_calc_planets(const struct assist_ephem* ephem, double jde, double rel, int code, double* GM, double* out_x, double* out_y, double* out_z, double* out_vx, double* out_vy, double* out_vz, double* out_ax, double* out_ay, double* out_az)
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
    const double au = ephem->AU;
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

// Load both constants and masses from SPK file comments in one pass
struct spk_constants_and_masses assist_load_spk_constants_and_masses(const char *path) {
    struct spk_constants_and_masses data = {0};
    
    // Try opening file.
    int fd = open(path, O_RDONLY);
    if (fd < 0){
        return data;
    }

    // Load the file record
    union record_t * record = assist_load_spk_file_record(fd);
    if (!record) {
        close(fd);
        return data;
    }

    char *comments = NULL;
    parse_comments(fd, record->file.fward, &comments);
    free(record);

    if (comments == NULL) {
        close(fd);
        return data;
    }

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
                
                // Parse constants
                if (strcmp(key, "cau") == 0) data.AU = value;
                else if (strcmp(key, "EMRAT") == 0) data.EMRAT = value;
                else if (strcmp(key, "J2E") == 0) data.J2E = value;
                else if (strcmp(key, "J3E") == 0) data.J3E = value;
                else if (strcmp(key, "J4E") == 0) data.J4E = value;
                else if (strcmp(key, "J2SUN") == 0) data.J2SUN = value;
                else if (strcmp(key, "AU") == 0) data.AU = value;
                else if (strcmp(key, "RE") == 0) data.RE = value;
                else if (strcmp(key, "CLIGHT") == 0) data.CLIGHT = value;
                else if (strcmp(key, "ASUN") == 0) data.ASUN = value;
                // Parse masses
                else if (strncmp(key, "GM", 2) == 0 || strncmp(key, "MA", 2) == 0) {
                    data.masses.names = realloc(data.masses.names, (data.masses.count + 1) * sizeof(char *));
                    data.masses.values = realloc(data.masses.values, (data.masses.count + 1) * sizeof(double));
                    data.masses.names[data.masses.count] = strdup(key);
                    data.masses.values[data.masses.count] = value;
                    data.masses.count++;
                }
            }
        }

        line = strtok(NULL, "\n");
    }

    free(comments);
    close(fd);
    return data;
}

// Apply constants from parsed data to ephem structure
void assist_apply_spk_constants(struct assist_ephem* ephem, const struct spk_constants_and_masses* data) {
    if (!ephem || !data) return;
    
    // Apply base constants
    ephem->AU = data->AU;
    ephem->EMRAT = data->EMRAT;
    ephem->J2E = data->J2E;
    ephem->J3E = data->J3E;
    ephem->J4E = data->J4E;
    ephem->J2SUN = data->J2SUN;
    ephem->RE = data->RE;
    ephem->CLIGHT = data->CLIGHT;
    ephem->ASUN = data->ASUN;
    
    // Calculate derived constants
    ephem->Re_eq = ephem->RE / ephem->AU;                    // Earth radius in AU
    ephem->Rs_eq = ephem->ASUN / ephem->AU;                  // Sun radius in AU
    ephem->c_AU_per_day = (ephem->CLIGHT / ephem->AU) * 86400; // Speed of light in AU/day
    ephem->c_squared = ephem->c_AU_per_day * ephem->c_AU_per_day; // c^2 in (AU/day)^2
    ephem->over_c_squared = 1.0 / ephem->c_squared;        // 1/c^2 in (day/AU)^2
}

// Free constants and masses data structure
void assist_free_spk_constants_and_masses(struct spk_constants_and_masses* data) {
    if (data && data->masses.names) {
        for (int i = 0; i < data->masses.count; i++) {
            free(data->masses.names[i]);
        }
        free(data->masses.names);
        free(data->masses.values);
        data->masses.names = NULL;
        data->masses.values = NULL;
        data->masses.count = 0;
    }
}
