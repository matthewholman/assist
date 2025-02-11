#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include "spk.h"


void assert_constants_exist(struct spk_global *sg) {
    struct spk_constants *sc = &sg->con;
    assert(sc->AU != 0.0);
    assert(sc->EMRAT != 0.0);
    assert(sc->J2E != 0.0);
    assert(sc->J3E != 0.0);
    assert(sc->J4E != 0.0);
    assert(sc->J2SUN != 0.0);
    assert(sc->RE != 0.0);
    assert(sc->CLIGHT != 0.0);
    assert(sc->ASUN != 0.0);
    assert(sg->masses.count = 420);
    for (int i = 0; i < sg->masses.count; i++) {
        assert(sg->masses.names[i] != NULL);
        assert(sg->masses.values[i] != 0.0);
    }
}

int main(int argc, char *argv[]) {

    const char *planets_path = "../../data/de440.bsp";

    struct spk_global *sg = assist_load_spk_constants(planets_path);
    if (!sg) {
        fprintf(stderr, "Failed to initialize SPK data.\n");
        return EXIT_FAILURE;
    }

    assert_constants_exist(sg);

    return EXIT_SUCCESS;
}