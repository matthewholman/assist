// Ephem slicer
// (c) Hanno Rein 2024
// A tool to trim down the DE440 ephemeris file, covering a smaller baseline
//
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <fcntl.h>

#include "planets.h" // for jpl_s struct

// https://gist.github.com/dgoguerra/7194777
static const char *humanSize(uint64_t bytes) {
    char *suffix[] = {"B", "KB", "MB", "GB", "TB"};
    char length = sizeof(suffix) / sizeof(suffix[0]);

    int i = 0;
    double dblBytes = bytes;

    if (bytes > 1024) {
        for (i = 0; (bytes / 1024) > 0 && i<length-1; i++, bytes /= 1024)
            dblBytes = bytes / 1024.0;
    }

    static char output[200];
    sprintf(output, "%.02lf %s", dblBytes, suffix[i]);
    return output;
}


// NOVAS. Fliegel, H. & Van Flandern, T.  Comm. of the ACM, Vol. 11, No. 10, October 1968, p. 657.
void cal_date (double tjd, short int *year, short int *month, short int *day, double *hour) {
   long int jd, k, m, n;

   double djd;

   djd = tjd + 0.5;
   jd = (long int) djd;

   *hour = fmod (djd,1.0) * 24.0;

   k     = jd + 68569L;
   n     = 4L * k / 146097L;

   k     = k - (146097L * n + 3L) / 4L;
   m     = 4000L * (k + 1L) / 1461001L;
   k     = k - 1461L * m / 4L + 31L;

   *month = (short int) (80L * k / 2447L);
   *day   = (short int) (k - 2447L * (long int) *month / 80L);
   k      = (long int) *month / 11L;

   *month = (short int) ((long int) *month + 2L - 12L * k);
   *year  = (short int) (100L * (n - 49L) + m + k);

   return;
}

void printf_jd(double tjd){
    short int year, month, day;
    double hour;
    cal_date(tjd, &year, &month, &day, &hour);
    printf("%d-%d-%d", year, month, day);
}



int main(int argc, char **argv) {
    ssize_t ret;
    int fd;

    if (argc<2){
        printf("Usage: ./ephem_slicer input_file [first_block] [number_of_blocks] \n\n");
        printf("If first_block and number_of_blocks not given, then the input file is analyzed but no new file is created.\n\n");
        printf("Examples: ./ephem_slicer ../../data/linux_p1550p2650.440 \n");
        printf("          ./ephem_slicer ../../data/linux_p1550p2650.440 10 3 \n");
        exit(1);
    }


    if ((fd = open(argv[1], O_RDONLY)) < 0){
        fprintf(stderr, "Cannot open input file.\n");
        exit(1);
    }

    // skip the header and constant names for now
    if (lseek(fd, 0x0A5C, SEEK_SET) < 0){
        close(fd);
        fprintf(stderr, "Error while seeking to header.\n");
        exit(1);
    }

    struct jpl_s* jpl = calloc(1, sizeof(struct jpl_s));

    // read header
    ret  = read(fd, &jpl->beg, sizeof(double));     // Start JD
    ret += read(fd, &jpl->end, sizeof(double));     // End JD
    ret += read(fd, &jpl->inc, sizeof(double));     // Days per block
    ret += read(fd, &jpl->num, sizeof(int32_t));    // Number of constants
    ret += read(fd, &jpl->cau, sizeof(double));     // AU to km 
    ret += read(fd, &jpl->cem, sizeof(double));     // Ratio between Earth/Moon

    // number of coefficients for all components
    for (int p = 0; p < JPL_N; p++){
        jpl->ncm[p] = 3;
    }
    // exceptions:
    jpl->ncm[JPL_NUT] = 2; // nutations
    jpl->ncm[JPL_TDB] = 1; // TT-TDB

    for (int p = 0; p < 12; p++) {                      // Columns 1-12 of Group 1050
        ret += read(fd, &jpl->off[p], sizeof(int32_t));
        ret += read(fd, &jpl->ncf[p], sizeof(int32_t));
        ret += read(fd, &jpl->niv[p], sizeof(int32_t));
    }

    ret += read(fd, &jpl->ver,     sizeof(int32_t));    // Version. e.g. 440
    ret += read(fd, &jpl->off[12], sizeof(int32_t));    // Columns 13 of Group 1050
    ret += read(fd, &jpl->ncf[12], sizeof(int32_t));
    ret += read(fd, &jpl->niv[12], sizeof(int32_t));

    // Get all the constant names
    jpl->str = calloc(jpl->num, sizeof(char *));

    // retrieve the names of the first 400 constants
    lseek(fd, 0x00FC, SEEK_SET);    
    for (int p = 0; p < 400; p++) {     // Group 1040
        jpl->str[p] = calloc(8, sizeof(char));
        read(fd, jpl->str[p], 6);
    }

    // read the remaining constant names
    lseek(fd, 0x0B28, SEEK_SET);
    for (int p = 400; p < jpl->num; p++) {
        jpl->str[p] = calloc(8, sizeof(char));
        read(fd, jpl->str[p], 6);
    }

    for (int p = 13; p < 15; p++) {                     // Columns 14 and 15 of Group 1050
        ret += read(fd, &jpl->off[p], sizeof(int32_t));
        ret += read(fd, &jpl->ncf[p], sizeof(int32_t));
        ret += read(fd, &jpl->niv[p], sizeof(int32_t));
    }

    // adjust for correct indexing (ie: zero based)
    for (int p = 0; p < JPL_N; p++){
        jpl->off[p] -= 1;
    }

    // save file size, and determine 'kernel size' or 'block size' (=8144 bytes for DE440/441)
    jpl->rec = sizeof(double) * 2;

    for (int p = 0; p < JPL_N; p++){
        jpl->rec += sizeof(double) * jpl->ncf[p] * jpl->niv[p] * jpl->ncm[p];
    }

    unsigned int nblocks = (jpl->end-jpl->beg)/jpl->inc;

    printf("Original file (%s):\n", argv[1]);
    printf("Days per block:     %.4f\n", jpl->inc);
    printf("Number of blocks:   %u\n", nblocks);
    printf("Block size (bytes): %zu\n", jpl->rec);
    printf("Start (JD):  %.4f   ", jpl->beg);
    printf_jd(jpl->beg);
    printf("\nEnd   (JD):  %.4f   ", jpl->end);
    printf_jd(jpl->end);
    printf("\nHeader size: %s\n", humanSize((uint64_t)(jpl->rec*2)));
    printf("Data size:   %s\n", humanSize((uint64_t)(jpl->rec*nblocks)));

    if (argc<4){
        printf("No new file created.\n");
        exit(0);
    }
    int first_block = atoi(argv[2]);    
    int num_blocks = atoi(argv[3]);    
    char* outfilename = malloc(sizeof(char)*1024);
    strcpy(outfilename, "linux_");
    double beg = jpl->beg+jpl->inc*first_block;
    double end = jpl->beg+jpl->inc*(first_block+num_blocks);
        

    
    short int beg_year, beg_month, beg_day;
    double beg_hour;
    short int end_year, end_month, end_day;
    double end_hour;
    cal_date(beg, &beg_year, &beg_month, &beg_day, &beg_hour);
    cal_date(end, &end_year, &end_month, &end_day, &end_hour);
    
   
    if (beg_year<0) {
        sprintf(outfilename + strlen(outfilename), "m%04d", -beg_year);
    }else{
        sprintf(outfilename + strlen(outfilename), "p%04d", beg_year);
    }
    if (end_year<0) {
        sprintf(outfilename + strlen(outfilename), "m%04d.440", -end_year);
    }else{
        sprintf(outfilename + strlen(outfilename), "p%04d.440", end_year);
    }


    FILE* fdo;
    if ((fdo = fopen(outfilename, "w")) < 0){
        fprintf(stderr, "Cannot open output file.\n");
        exit(1);
    }

    // copy header
    lseek(fd, 0, SEEK_SET);
    char* buf = malloc(2*jpl->rec);
    read(fd, buf, 2*jpl->rec);
    fwrite(buf, 2*jpl->rec, 1, fdo); 

    printf("\n\nNew file (%s):\n", outfilename);
    printf("Start (JD):  %.4f   ", beg);
    printf_jd(beg);
    printf("\nEnd   (JD):  %.4f   ", end);
    printf_jd(end);
    printf("\nHeader size: %s\n", humanSize((uint64_t)(jpl->rec*2)));
    printf("Data size:   %s\n", humanSize((uint64_t)(jpl->rec*num_blocks)));
    
    // Copy data blocks
    for (int blk = first_block; blk<first_block + num_blocks; blk++){
        lseek(fd, (blk+2)*jpl->rec, SEEK_SET);
        read(fd, buf, jpl->rec);
        fwrite(buf, jpl->rec, 1, fdo); 
    }

    // Overwrite header dates
    fseek(fdo, 0x0A5C, SEEK_SET);
    fwrite(&beg, sizeof(double), 1, fdo);
    fwrite(&end, sizeof(double), 1, fdo);

    if (close(fd) < 0) { 
        fprintf(stderr, "Error while closing input file.\n");
        exit(1);
    }
    if (fclose(fdo) < 0) { 
        fprintf(stderr, "Error while closing output file.\n");
        exit(1);
    }

    return 0;

}

