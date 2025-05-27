#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include "assist.h"
#include "spk.h"


void assert_constants_exist(struct assist_ephem *ephem) {
    assert(ephem->AU != 0.0);
    assert(ephem->EMRAT != 0.0);
    assert(ephem->J2E != 0.0);
    assert(ephem->J3E != 0.0);
    assert(ephem->J4E != 0.0);
    assert(ephem->J2SUN != 0.0);
    assert(ephem->RE != 0.0);
    assert(ephem->CLIGHT != 0.0);
    assert(ephem->ASUN != 0.0);
    printf("All constants exist.\n");
}

int main(int argc, char *argv[]) {

    const char *planets_path = "../../data/de440.bsp";

    struct assist_ephem ephem = {0};
    struct spk_constants_and_masses data = assist_load_spk_constants_and_masses(planets_path);
    assist_apply_spk_constants(&ephem, &data);
    assist_free_spk_constants_and_masses(&data);

    assert_constants_exist(&ephem);

    return EXIT_SUCCESS;
}