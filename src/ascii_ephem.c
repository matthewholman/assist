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
#include "ascii_ephem.h"
#include "assist.h"

/*
 *  assist_ascii_work
 *
 *  Interpolate the appropriate Chebyshev polynomial coefficients.
 *
 *      ncf - number of coefficients per component
 *      ncm - number of components (ie: 3 for most)
 *      niv - number of intervals / sets of coefficients
 *
 */

void assist_ascii_work(double *P, int ncm, int ncf, int niv, double t0, double t1, double *u, double *v, double *w)
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
 *  assist_ascii_init
 *
 *  Initialise everything needed to read an ASCII-derived binary ephemeris
 *  file (.440/.441).
 *
 */

static double getConstant(struct ascii_s* ascii, char* name){
    for (int p = 0; p < ascii->num; p++) {
        if (strncmp(name, ascii->str[p], 6) == 0){
            return ascii->con[p];
        }
    }
    fprintf(stderr,"WARNING: Constant [%s] not found in ephemeris file.\n",name); 
    return 0;
}

int assist_ascii_find_constant(const struct ascii_s* ascii, const char* name, double* out_value){
    if (out_value) *out_value = 0.0;
    if (ascii == NULL || ascii->str == NULL || ascii->con == NULL || name == NULL || out_value == NULL){
        return 0;
    }

    // Constant names are exactly 6 bytes, space-padded.
    char key[6];
    for (int i=0; i<6; i++) key[i] = ' ';
    for (int i=0; i<6 && name[i] != '\0'; i++) key[i] = name[i];

    for (int p=0; p<ascii->num; p++){
        if (memcmp(ascii->str[p], key, 6) == 0){
            *out_value = ascii->con[p];
            return 1;
        }
    }
    return 0;
}

struct ascii_s * assist_ascii_init(char *str)
{
    struct stat sb;
    ssize_t ret;
    int fd;

    if ((fd = open(str, O_RDONLY)) < 0){
        return NULL;
    }

    if (fstat(fd, &sb) < 0){
        close(fd);
        fprintf(stderr, "Error while trying to determine filesize.\n");
        return NULL;
    }
    

    // skip the header and constant names for now
    if (lseek(fd, 0x0A5C, SEEK_SET) < 0){
        close(fd);
        fprintf(stderr, "Error while seeking to header.\n");
        return NULL;
    }

    struct ascii_s* ascii = calloc(1, sizeof(struct ascii_s));

    // read header
    ret  = read(fd, &ascii->beg, sizeof(double));     // Start JD
    ret += read(fd, &ascii->end, sizeof(double));     // End JD
    ret += read(fd, &ascii->inc, sizeof(double));     // Days per block
    ret += read(fd, &ascii->num, sizeof(int32_t));    // Number of constants
    ret += read(fd, &ascii->cau, sizeof(double));     // AU to km 
    ret += read(fd, &ascii->cem, sizeof(double));     // Ratio between Earth/Moon

    // number of coefficients for all components
    for (int p = 0; p < ASCII_N; p++){
        ascii->ncm[p] = 3;
    }
    // exceptions:
    ascii->ncm[ASCII_NUT] = 2; // nutations
    ascii->ncm[ASCII_TDB] = 1; // TT-TDB

    for (int p = 0; p < 12; p++) {                      // Columns 1-12 of Group 1050
        ret += read(fd, &ascii->off[p], sizeof(int32_t));
        ret += read(fd, &ascii->ncf[p], sizeof(int32_t));
        ret += read(fd, &ascii->niv[p], sizeof(int32_t));
    }

    ret += read(fd, &ascii->ver,     sizeof(int32_t));    // Version. e.g. 440
    ret += read(fd, &ascii->off[12], sizeof(int32_t));    // Columns 13 of Group 1050
    ret += read(fd, &ascii->ncf[12], sizeof(int32_t));
    ret += read(fd, &ascii->niv[12], sizeof(int32_t));

    // Get all the constant names
    ascii->str = calloc(ascii->num, sizeof(char *));

    // retrieve the names of the first 400 constants
    lseek(fd, 0x00FC, SEEK_SET);    
    for (int p = 0; p < 400; p++) {     // Group 1040
        ascii->str[p] = calloc(8, sizeof(char));
        read(fd, ascii->str[p], 6);
    }

    // read the remaining constant names
    lseek(fd, 0x0B28, SEEK_SET);
    for (int p = 400; p < ascii->num; p++) {
        ascii->str[p] = calloc(8, sizeof(char));
        read(fd, ascii->str[p], 6);
    }

    for (int p = 13; p < 15; p++) {                     // Columns 14 and 15 of Group 1050
        ret += read(fd, &ascii->off[p], sizeof(int32_t));
        ret += read(fd, &ascii->ncf[p], sizeof(int32_t));
        ret += read(fd, &ascii->niv[p], sizeof(int32_t));
    }
    (void)ret; // silence unused-but-set warnings (read errors are handled via fstat/mmap usage)

    // adjust for correct indexing (ie: zero based)
    for (int p = 0; p < ASCII_N; p++){
        ascii->off[p] -= 1;
    }

    // save file size, and determine 'kernel size' or 'block size' (=8144 bytes for DE440/441)
    ascii->len = sb.st_size;
    ascii->rec = sizeof(double) * 2;

    for (int p = 0; p < ASCII_N; p++){
        ascii->rec += sizeof(double) * ascii->ncf[p] * ascii->niv[p] * ascii->ncm[p];
    }

    // memory map the file, which makes us thread-safe with kernel caching
    ascii->map = mmap(NULL, ascii->len, PROT_READ, MAP_SHARED, fd, 0);

    if (ascii->map == NULL){ 
        close(fd);
        free(ascii); // note constants leak
        fprintf(stderr, "Error while calling mmap().\n");
        return NULL;
    }

    // Read constants
    ascii->con = calloc(ascii->num, sizeof(double));
    lseek(fd, ascii->rec, SEEK_SET); // Starts at offset of 1 block size
    for (int p = 0; p < ascii->num; p++){
        read(fd, &ascii->con[p], sizeof(double));
        //printf("%6d  %s   %.5e\n",p,ascii->str[p],ascii->con[p]);
    }


    // Find GM values (ASSIST_BODY indexing)
    ascii->mass[ASSIST_BODY_SUN]     = getConstant(ascii, "GMS   ");  // Sun
    ascii->mass[ASSIST_BODY_MERCURY] = getConstant(ascii, "GM1   ");  // Mercury
    ascii->mass[ASSIST_BODY_VENUS]   = getConstant(ascii, "GM2   ");
	double emrat = getConstant(ascii, "EMRAT  "); // Earth Moon Ratio
    double gmb = getConstant(ascii, "GMB   ");    // Earth Moon combined
    ascii->mass[ASSIST_BODY_EARTH]   = (emrat/(1.+emrat)) * gmb;    // Earth
    ascii->mass[ASSIST_BODY_MOON]    = 1./(1+emrat) * gmb;          // Moon
    ascii->mass[ASSIST_BODY_MARS]    = getConstant(ascii, "GM4   ");  // Mars
    ascii->mass[ASSIST_BODY_JUPITER] = getConstant(ascii, "GM5   ");  // Jupiter
    ascii->mass[ASSIST_BODY_SATURN]  = getConstant(ascii, "GM6   ");
    ascii->mass[ASSIST_BODY_URANUS]  = getConstant(ascii, "GM7   ");
    ascii->mass[ASSIST_BODY_NEPTUNE] = getConstant(ascii, "GM8   ");
    ascii->mass[ASSIST_BODY_PLUTO]   = getConstant(ascii, "GM9   "); // Pluto


    // Other constants
    ascii->J2E = getConstant(ascii, "J2E   ");
    ascii->J3E = getConstant(ascii, "J3E   ");
    ascii->J4E = getConstant(ascii, "J4E   ");
    ascii->J2SUN = getConstant(ascii, "J2SUN ");
    ascii->AU = getConstant(ascii, "AU    ");
    ascii->RE = getConstant(ascii, "RE    ");
    ascii->CLIGHT = getConstant(ascii, "CLIGHT");
    ascii->ASUN = getConstant(ascii, "ASUN  ");

    // this file descriptor is no longer needed since we are memory mapped
    if (close(fd) < 0) { 
        fprintf(stderr, "Error while closing file.\n");
    }
#if defined(MADV_RANDOM)
    if (madvise(ascii->map, ascii->len, MADV_RANDOM) < 0){
        fprintf(stderr, "Error during madvise.\n");
    }
#endif

    return ascii;

}

void assist_ascii_free(struct ascii_s *ascii) {
    if (ascii == NULL){
        return;
    }
    if (munmap(ascii->map, ascii->len) < 0){ 
        fprintf(stderr, "Error during munmap().\n");
    }
    for (int p = 0; p < ascii->num; p++){
        free(ascii->str[p]);
    }
    free(ascii->str);
    free(ascii->con);
    free(ascii);
}

/*
 *  assist_ascii_calc
 *
 *  Caculate the position+velocity in _equatorial_ coordinates.
 *  Assumes pos is initially zero.
 */
enum ASSIST_STATUS assist_ascii_calc(struct ascii_s *ascii, double jd_ref, double jd_rel, int body, 
		 double* const GM,
		 double* const out_x, double* const out_y, double* const out_z,
		 double* const out_vx, double* const out_vy, double* const out_vz,
         double* const out_ax, double* const out_ay, double* const out_az){
    double t, *z;
    uint32_t blk;

    if (ascii == NULL || ascii->map == NULL)
        return ASSIST_ERROR_EPHEM_FILE;
    if(body<0 || body >= ASSIST_BODY_NPLANETS)
	    return(ASSIST_ERROR_NEPHEM);
    
    struct mpos_s pos;

    // Get mass, position, velocity, and mass of body i in barycentric coords.
    *GM = ascii->mass[body];

        // check if covered by this file
        if (jd_ref + jd_rel < ascii->beg || jd_ref + jd_rel > ascii->end)
            return ASSIST_ERROR_COVERAGE;

        // compute record number and 'offset' into record
        blk = (uint32_t)((jd_ref + jd_rel - ascii->beg) / ascii->inc);
        z = (double*)ascii->map + (blk + 2) * ascii->rec/sizeof(double);
        t = ((jd_ref - ascii->beg - (double)blk * ascii->inc) + jd_rel) / ascii->inc;

        switch (body) { // The indices in the ascii-> arrays match the column index for the body
            case ASSIST_BODY_SUN:
                assist_ascii_work(&z[ascii->off[ASCII_SUN]], ascii->ncm[ASCII_SUN], ascii->ncf[ASCII_SUN], ascii->niv[ASCII_SUN], t, ascii->inc, pos.u, pos.v, pos.w);
                break;
            case ASSIST_BODY_MERCURY:
                assist_ascii_work(&z[ascii->off[ASCII_MER]], ascii->ncm[ASCII_MER], ascii->ncf[ASCII_MER], ascii->niv[ASCII_MER], t, ascii->inc, pos.u, pos.v, pos.w);
                break;
            case ASSIST_BODY_VENUS:
                assist_ascii_work(&z[ascii->off[ASCII_VEN]], ascii->ncm[ASCII_VEN], ascii->ncf[ASCII_VEN], ascii->niv[ASCII_VEN], t, ascii->inc, pos.u, pos.v, pos.w);
                break;
            case ASSIST_BODY_EARTH:
                {
                    struct mpos_s emb, lun;
                    assist_ascii_work(&z[ascii->off[ASCII_EMB]], ascii->ncm[ASCII_EMB], ascii->ncf[ASCII_EMB], ascii->niv[ASCII_EMB], t, ascii->inc, emb.u, emb.v, emb.w); // earth moon barycenter
                    assist_ascii_work(&z[ascii->off[ASCII_LUN]], ascii->ncm[ASCII_LUN], ascii->ncf[ASCII_LUN], ascii->niv[ASCII_LUN], t, ascii->inc, lun.u, lun.v, lun.w);

                    vecpos_set(pos.u, emb.u);
                    vecpos_off(pos.u, lun.u, -1.0 / (1.0 + ascii->cem));

                    vecpos_set(pos.v, emb.v);
                    vecpos_off(pos.v, lun.v, -1.0 / (1.0 + ascii->cem));

                    vecpos_set(pos.w, emb.w);
                    vecpos_off(pos.w, lun.w, -1.0 / (1.0 + ascii->cem));
                }
                break;
            case ASSIST_BODY_MOON: 
                {
                    struct mpos_s emb, lun;
                    assist_ascii_work(&z[ascii->off[ASCII_EMB]], ascii->ncm[ASCII_EMB], ascii->ncf[ASCII_EMB], ascii->niv[ASCII_EMB], t, ascii->inc, emb.u, emb.v, emb.w);
                    assist_ascii_work(&z[ascii->off[ASCII_LUN]], ascii->ncm[ASCII_LUN], ascii->ncf[ASCII_LUN], ascii->niv[ASCII_LUN], t, ascii->inc, lun.u, lun.v, lun.w);

                    vecpos_set(pos.u, emb.u);
                    vecpos_off(pos.u, lun.u, ascii->cem / (1.0 + ascii->cem));

                    vecpos_set(pos.v, emb.v);
                    vecpos_off(pos.v, lun.v, ascii->cem / (1.0 + ascii->cem));

                    vecpos_set(pos.w, emb.w);
                    vecpos_off(pos.w, lun.w, ascii->cem / (1.0 + ascii->cem));
                }
                break;
            case ASSIST_BODY_MARS:
                assist_ascii_work(&z[ascii->off[ASCII_MAR]], ascii->ncm[ASCII_MAR], ascii->ncf[ASCII_MAR], ascii->niv[ASCII_MAR], t, ascii->inc, pos.u, pos.v, pos.w);
                break;
            case ASSIST_BODY_JUPITER:
                assist_ascii_work(&z[ascii->off[ASCII_JUP]], ascii->ncm[ASCII_JUP], ascii->ncf[ASCII_JUP], ascii->niv[ASCII_JUP], t, ascii->inc, pos.u, pos.v, pos.w);
                break;
            case ASSIST_BODY_SATURN:
                assist_ascii_work(&z[ascii->off[ASCII_SAT]], ascii->ncm[ASCII_SAT], ascii->ncf[ASCII_SAT], ascii->niv[ASCII_SAT], t, ascii->inc, pos.u, pos.v, pos.w);
                break;
            case ASSIST_BODY_URANUS:
                assist_ascii_work(&z[ascii->off[ASCII_URA]], ascii->ncm[ASCII_URA], ascii->ncf[ASCII_URA], ascii->niv[ASCII_URA], t, ascii->inc, pos.u, pos.v, pos.w);
                break;
            case ASSIST_BODY_NEPTUNE:
                assist_ascii_work(&z[ascii->off[ASCII_NEP]], ascii->ncm[ASCII_NEP], ascii->ncf[ASCII_NEP], ascii->niv[ASCII_NEP], t, ascii->inc, pos.u, pos.v, pos.w);
                break;
            case ASSIST_BODY_PLUTO:
                assist_ascii_work(&z[ascii->off[ASCII_PLU]], ascii->ncm[ASCII_PLU], ascii->ncf[ASCII_PLU], ascii->niv[ASCII_PLU], t, ascii->inc, pos.u, pos.v, pos.w);
                break;
            default:
                return ASSIST_ERROR_NEPHEM; // body not found
                break;
        }

    // Convert to au/day and au/day^2
    vecpos_div(pos.u, ascii->cau);
    vecpos_div(pos.v, ascii->cau/86400.);
    vecpos_div(pos.w, ascii->cau/(86400.*86400.));

    *out_x = pos.u[0];
    *out_y = pos.u[1];
    *out_z = pos.u[2];
    *out_vx = pos.v[0];
    *out_vy = pos.v[1];
    *out_vz = pos.v[2];
    *out_ax = pos.w[0];
    *out_ay = pos.w[1];
    *out_az = pos.w[2];

    return(ASSIST_SUCCESS);

}

enum ASSIST_STATUS assist_ascii_calc_from_ephem(const struct assist_ephem* ephem, double jd_ref, double jd_rel, int body, 
		 double* const GM,
		 double* const x, double* const y, double* const z,
		 double* const vx, double* const vy, double* const vz,
         double* const ax, double* const ay, double* const az){
    return assist_ascii_calc(ephem->ascii_planets, jd_ref, jd_rel, body, GM, x, y, z, vx, vy, vz, ax, ay, az);
}


