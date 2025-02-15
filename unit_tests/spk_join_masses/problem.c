#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include "spk.h"


void assert_populated_masses(struct spk_s *spk_data) {
    if (!spk_data) {
        printf("SPK data is NULL.\n");
        return;
    }

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

    struct spk_global *sg = assist_load_spk_constants(planets_path);
    if (!sg) {
        fprintf(stderr, "Failed to initialize SPK data.\n");
        return EXIT_FAILURE;
    }
    struct spk_s *planet_data = assist_spk_init(planets_path);
    if (!planet_data) {
        fprintf(stderr, "Failed to initialize SPK data.\n");
        return EXIT_FAILURE;
    }

    assist_spk_join_masses(planet_data, sg);
    assert_populated_masses(planet_data);

    struct spk_s *asteroid_data = assist_spk_init(asteroids_path);
    if (!asteroid_data) {
        fprintf(stderr, "Failed to initialize SPK data.\n");
        return EXIT_FAILURE;
    }

    assist_spk_join_masses(asteroid_data, sg);
    assert_populated_masses(asteroid_data);

    assist_spk_free(planet_data);
    assist_spk_free(asteroid_data);

    assist_free_spk_constants(sg);

    return EXIT_SUCCESS;
}