#include <stdio.h>
#include <stdlib.h>
#include "spk.h"

/*
* Test loading combinations of different ephem and spk files.
*/


int main(int argc, char *argv[]) {
    if (argc != 3) {
        fprintf(stderr, "Usage: %s <path to planets .bsp file> <path to asteroids .bsp file>\n", argv[0]);
        return EXIT_FAILURE;
    }

    char *planets_path = argv[1];
    char *asteroids_path = argv[2];
    struct assist_ephem *ephem = assist_ephem_create(planets_path, asteroids_path);
    if (!ephem) {
        fprintf(stderr, "Failed to initialize ephemeris data.\n");
        return EXIT_FAILURE;
    }

    return EXIT_SUCCESS;
}