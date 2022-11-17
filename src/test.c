// compile with
// gcc -o test test.c spk.c -lm -O3 -ggdb

#ifndef _GNU_SOURCE
#define _GNU_SOURCE
#endif

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "spk.h"

#define C_AU	149597870.7	// [km]

static double _rad(double *u)
{
	return sqrt(u[0] * u[0] + u[1] * u[1] + u[2] * u[2]);
}

int main(int argc, char **argv)
{
	struct spk_s *pl;
	struct mpos_s pos;
	double jde, rel, r;
	int p, n;

	memset(&pos, 0, sizeof(struct mpos_s));
	jde = 2450000.0;

	for (p = 1; p < argc; p++) {
		if ((pl = spk_init(argv[p])) == NULL)
			{ perror(argv[p]); continue; }

		rel = 10000.0 * drand48();
		fprintf(stdout, "# %s %d %.6f\n", argv[0], pl->num, jde + rel);

		// use spk_find to look for a specific body

		for (n = 0; n < pl->num; n++) {
			spk_calc(pl, n, jde, rel, &pos);
			r = _rad(pos.u);

			fprintf(stdout, "%02d %3d %3d %+12.12e %+12.12e %+12.12e %+12.12e %+12.12e %+12.12e -> %.3f\n",
				n, pl->tar[n], pl->cen[n],
				pos.u[0], pos.u[1], pos.u[2],
				pos.v[0], pos.v[1], pos.v[2],
				r);
		}

		spk_free(pl);
	}

	return 0;
}


