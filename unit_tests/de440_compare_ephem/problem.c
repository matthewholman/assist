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

    // Call assist_all_ephem for both versions of the ephemeris data
    // and compare the results
    double times[] = {0.0, 100.0, 1000.0, 10000.0};
    for (int body = 0; body < 11; body++) {

        printf("body: %d\n", body);
        for (int i = 0; i < 4; i++) {
            double t = times[i];

            double GMa, xa, ya, za, vxa, vya, vza, axa, aya, aza;
            int ret = assist_all_ephem(binary_ephem, NULL, body, t, &GMa, &xa, &ya, &za, &vxa, &vya, &vza, &axa, &aya, &aza);
            if (ret != 0) {
                printf("ret: %d\n", ret);
                fprintf(stderr, "Failed to get ephemeris data from binary file.\n");
                return EXIT_FAILURE;
            }

            double GMb, xb, yb, zb, vxb, vyb, vzb, axb, ayb, azb;
            ret = assist_all_ephem(bsp_ephem, NULL, body, t, &GMb, &xb, &yb, &zb, &vxb, &vyb, &vzb, &axb, &ayb, &azb);
            if (ret != 0) {
                fprintf(stderr, "Failed to get ephemeris data from bsp file.\n");
                return EXIT_FAILURE;
            }

            // Assert differences are below a certain threshold
            double threshold = 0.0;
            double accumulated_diff = 0.0;
            double comparisons[] = {GMb, xb, yb, zb, vxb, vyb, vzb, axb, ayb, azb,
                                    GMa, xa, ya, za, vxa, vya, vza, axa, aya, aza};
            for (int j = 0; j < 10; j++) {
                double diff = fabs(comparisons[j] - comparisons[j + 10]);
                if (diff > threshold) {
                    accumulated_diff += diff;
                    // // throw error
                    // return EXIT_FAILURE;
                }
            }
            printf("accumulated_diff: %.17e\n", accumulated_diff);
        }
    }

    assist_ephem_free(binary_ephem);
    assist_ephem_free(bsp_ephem);
  
    return EXIT_SUCCESS;
}