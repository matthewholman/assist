// compile with
// gcc -o get_mass get_mass.c planets.c spk.c -lm -O3 -ggdb

#ifndef _GNU_SOURCE
#define _GNU_SOURCE
#endif

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "planets.h"
#include "spk.h"

int main(int argc, char **argv)
{
    struct jpl_s *pl;
    double m;
    int n;

    if ((pl = jpl_init()) == NULL)
        abort();

    for (n = 0; n < 1000; n++) {
        m = jpl_mass(pl, n);

        if (m > 0.0)
            fprintf(stdout, "%d\t%e\n", n, m);
    }

    jpl_free(pl);
    return 0;
}
