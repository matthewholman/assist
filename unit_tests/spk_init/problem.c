#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include "spk.h"

int main(int argc, char *argv[]) {
    char *planet_path = "../../data/de440.bsp";
    char *asteroid_path = "../../data/sb441-n16.bsp";

    struct spk_s *planet_spk = assist_spk_init(planet_path);
    if (!planet_spk) {
        fprintf(stderr, "Failed to initialize SPK data.\n");
        return EXIT_FAILURE;
    }

    // Assert that there are 14 targets in the planet SPK data
    // 1 : Mercury
    // 2 : Venus
    // 3 : Earth-Moon barycenter
    // 4 : Mars
    // 5 : Jupiter
    // 6 : Saturn
    // 7 : Uranus
    // 8 : Neptune
    // 9 : Pluto
    // 10 : Sun
    // 301 : Moon
    // 399 : Earth
    // 199 : Mercury (Mercury barycenter)
    // 299 : Venus (Venus barycenter)

    assert(planet_spk->num == 14);
    int counted_targets = 0;

    for (int i = 0; i < 16; i++) {
        // Assuming a valid target has a non-zero code, adjust the condition as necessary
        if (planet_spk->targets[i].code != 0) {
            counted_targets++;
        }
    }
    assert(counted_targets == 14);

    assist_spk_free(planet_spk);

    struct spk_s *asteroid_spk = assist_spk_init(asteroid_path);
    if (!asteroid_spk) {
        fprintf(stderr, "Failed to initialize SPK data.\n");
        return EXIT_FAILURE;
    }

    // Assert that there are 16 targets in the asteroid SPK data
    assert(asteroid_spk->num == 16);
    counted_targets = 0;
    
    for (int i = 0; i < 20; i++) {
        // Assuming a valid target has a non-zero code, adjust the condition as necessary
        if (asteroid_spk->targets[i].code != 0) {
            counted_targets++;
        }
    }

    assist_spk_free(asteroid_spk);  

    return EXIT_SUCCESS;
}