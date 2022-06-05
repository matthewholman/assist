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
	int p, n, max, num;
	char *str;

	if (argc < 2) {
		fprintf(stderr, "usage: ./constants jpl_ephem.bin\n");
		exit(EXIT_FAILURE);
	}

	if ((fp = fopen("const.h", "w")) == NULL) {
		perror("const.h");
		exit(EXIT_FAILURE);
	}

	if (setenv("JPL_PLANET_EPHEM", argv[1], 1) < 0) {
		perror("setenv");
		exit(EXIT_FAILURE);
	}

	if ((pl = jpl_init()) == NULL)
		abort();

	fprintf(fp, "// const.h\n\n");
	fprintf(fp, "#ifndef _JPL_CONST_H\n");
	fprintf(fp, "#define _JPL_CONST_H\n\n");

	// write all the constants out ... we could save a default filename
	fprintf(fp, "//#define JPL_EPHEM_%-20s \"%s\"\n", "FILE", argv[1]);
	fprintf(fp, "#define JPL_EPHEM_%-20s %02d\n", "VER", pl->ver);
	fprintf(fp, "#define JPL_EPHEM_%-20s %02d\n", "NUM", pl->num);
	fprintf(fp, "#define JPL_EPHEM_%-20s %.6f\n", "BEG", pl->beg);
	fprintf(fp, "#define JPL_EPHEM_%-20s %.6f\n", "END", pl->end);
	fprintf(fp, "#define JPL_EPHEM_%-20s %.16e\n", "CEM", pl->cem);
	fprintf(fp, "#define JPL_EPHEM_%-20s %.16e\n", "CAU", pl->cau);

	for (p = 0; p < pl->num; p++)
		fprintf(fp, "#define JPL_EPHEM_%-20s %.16e\t// %d\n", pl->str[p], pl->con[p], p);
	/*
	for (num = max = p = 0; p < pl->num; p++) {
		// print all non-mass terms
	    
		if (pl->str[p][0] != 'M' || pl->str[p][1] != 'A')
			fprintf(fp, "#define JPL_EPHEM_%-20s %.16e\t// %d\n", pl->str[p], pl->con[p], p);
		else {
		    // FIXME : these could be used to show how many mass values are available from jpl_mass()
		    sscanf(pl->str[p], "MA%d", &max);
		    fprintf(fp, "\n#define JPL_EPHEM_AST %d\n", num);
		    //fprintf(fp, "#define JPL_EPHEM_MAX %d\n", max);
		    //printf("max: %d\n", max);
		    num += 1;
		}
	}
	*/


	fprintf(fp, "\n#endif // _JPL_CONST_H\n\n");
	jpl_free(pl);
	fclose(fp);
	return 0;
}

