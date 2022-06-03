
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

int main(int argc, char **argv)
{
	FILE *fp;
	char buf[256];
	char obj[64];
	char *str;
	off_t off;
	int p, tar;

	for (p = 1; p < argc; p++) {
		if ((fp = fopen(argv[p], "r")) == NULL)
			{ perror(argv[p]); continue; }

		off = 2048;

		while (!feof(fp)) {
			fseek(fp, off, SEEK_SET);
			fgets(buf, sizeof(buf), fp);
			off += 1 + strlen(buf);

			if (buf[0] == '\0')
				continue;
			//if (buf[0] == 0x04)
			//break;

			fprintf(stdout, "'%s' %d\n", buf, off);

			// name comes before number
			if ((str = strstr(buf, "Short Name:")) != NULL)
				sscanf(str, "%*[^:]: %s", obj);

			if ((str = strstr(buf, "ID  Number:")) != NULL) {
				sscanf(str, "%*[^:]: %d", &tar);
				fprintf(stdout, "#define SPK_NAIF_%-20s %d\n", obj, tar);
			}
		}

		if (fclose(fp) < 0)
			perror(argv[p]);
	}

	return 0;
}

