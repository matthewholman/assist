#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include "assist.h"
#include "spk.h"


void assert_populated_masses(struct spk_s *spk_data) {
    if (!spk_data) {
        printf("SPK data is NULL.\n");
        return;
    }

    // Ensure we have targets loaded
    assert(spk_data->num > 0);
    printf("Found %d targets in SPK data.\n", spk_data->num);

    for (int i = 0; i < spk_data->num; i++) {
        struct spk_target *target = &spk_data->targets[i];
        // We do not calculate mass for mercury and venus barycenters
        if (target->code != 199 && target->code != 299) {
            // Assert that the masses are not zero
            assert(target->mass != 0.0);
        }
    }
}


int main(int argc, char *argv[]) {

    const char *planets_path = "../../data/de440.bsp";
    const char *asteroids_path = "../../data/sb441-n16.bsp";

    struct assist_ephem ephem = {0};
    
    // Load planets
    ephem.spk_planets = assist_spk_init(planets_path);
    if (!ephem.spk_planets) {
        fprintf(stderr, "Failed to initialize planet SPK data.\n");
        return EXIT_FAILURE;
    }

    // Load asteroids
    ephem.spk_asteroids = assist_spk_init(asteroids_path);
    if (!ephem.spk_asteroids) {
        fprintf(stderr, "Failed to initialize asteroid SPK data.\n");
        return EXIT_FAILURE;
    }
    
    // Load constants and masses, then apply them
    struct spk_constants_and_masses data = assist_load_spk_constants_and_masses(planets_path);
    assist_apply_spk_constants(&ephem, &data);
    assist_spk_join_masses(ephem.spk_planets, &data.masses, ephem.EMRAT);
    assist_spk_join_masses(ephem.spk_asteroids, &data.masses, ephem.EMRAT);
    assist_free_spk_constants_and_masses(&data);
    
    assert_populated_masses(ephem.spk_planets);
    assert_populated_masses(ephem.spk_asteroids);

    assist_spk_free(ephem.spk_planets);
    assist_spk_free(ephem.spk_asteroids);

    return EXIT_SUCCESS;
}