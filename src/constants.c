// compile with
// gcc -o constants constants.c planets.c spk.c -lm -O3 -ggdb

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
	struct _jpl_s *pl;
	FILE *fp;
	int p;

	if ((fp = fopen("const.h", "w")) == NULL) {
		perror("const.h");
		exit(EXIT_FAILURE);
	}

	if ((pl = jpl_init()) == NULL)
		abort();

	fprintf(fp, "// const.h\n\n");
	fprintf(fp, "#ifndef _JPL_CONST_H\n");
	fprintf(fp, "#define _JPL_CONST_H\n\n");

	// write all the constants out ... we could save a default filename
	fprintf(fp, "#define JPL_EPHEM_%-20s \"%s\"\n", "FILE", "FIXME");
	fprintf(fp, "#define JPL_EPHEM_%-20s %02d\n", "VER", pl->ver);
	fprintf(fp, "#define JPL_EPHEM_%-20s %02d\n", "NUM", pl->num);
	fprintf(fp, "#define JPL_EPHEM_%-20s %.6f\n", "BEG", pl->beg);
	fprintf(fp, "#define JPL_EPHEM_%-20s %.6f\n", "END", pl->end);
	fprintf(fp, "#define JPL_EPHEM_%-20s %.16e\n", "CEM", pl->cem);
	fprintf(fp, "#define JPL_EPHEM_%-20s %.16e\n", "CAU", pl->cau);

	for (p = 0; p < pl->num; p++)
		fprintf(fp, "#define JPL_EPHEM_%-20s %.16e\t// %d\n", pl->str[p], pl->con[p], p);

	// we could write the masses into a constant array
//	fprintf(fp, "const double jpl_mass[] = {\n");
//	...

	fprintf(fp, "\n#endif // _JPL_CONST_H\n\n");
	jpl_free(pl);
	fclose(fp);
	return 0;
}

