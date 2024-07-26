#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>
#include "assist.h"
#include "forces.h"

/*
* Compare the planet ephemeris for the ascii/binary and bsp versions of DE440.
*
* Note that it appears the actual coefficient values in the files vary
* in the range of 1-2 units in the last decimal place.
* 
*/


int main(int argc, char *argv[]) {

    char *linux_binary_path = "../../data/linux_p1550p2650.440";
    char *bsp_path = "../../data/de440.bsp";

    struct assist_ephem *binary_ephem = assist_ephem_create(linux_binary_path, NULL);
    if (!binary_ephem) {
        fprintf(stderr, "Failed to initialize ephemeris data.\n");
        return EXIT_FAILURE;
    }

    struct assist_ephem *bsp_ephem = assist_ephem_create(bsp_path, NULL);
    if (!bsp_ephem) {
        fprintf(stderr, "Failed to initialize ephemeris data.\n");
        return EXIT_FAILURE;
    }

    const char *constant_names[] = {
        "AU",
        "EMRAT",
        "J2E",
        "J3E",
        "J4E",
        "J2SUN",
        "RE",
        "CLIGHT",
        "ASUN"
    };

    const double threshold = 1e-20;

    for (int i = 0; i < 9; i++) {
        double binary_constant = assist_get_constant(binary_ephem, constant_names[i]);
        double bsp_constant = assist_get_constant(bsp_ephem, constant_names[i]);
        printf("%s: %f %f\n", constant_names[i], binary_constant, bsp_constant);
        assert(fabs(binary_constant - bsp_constant) < threshold);
    };

    assist_ephem_free(binary_ephem);
    assist_ephem_free(bsp_ephem);
  
    return EXIT_SUCCESS;
}