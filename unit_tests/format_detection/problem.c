#include <stdio.h>
#include <stdlib.h>
#include <fcntl.h>
#include <unistd.h>
#include <string.h>
#include "assist.h"

int main(){
    // Try detect real files if present
    int failures = 0;

    // Detect SPK (.bsp) if available
    int fd;
    fd = open("../../data/de440.bsp", O_RDONLY);
    if (fd>=0){
        ephemeris_file_format_t f = assist_detect_ephemeris_file_format(fd);
        if (f != FILE_FORMAT_VALID_BSP){
            fprintf(stderr, "Expected SPK (.bsp) to be valid. Got %d\n", (int)f);
            failures++;
        }
        close(fd);
    } else {
        // Not fatal; skip if file missing
    }

    // Detect ASCII-derived binary ephemeris (.440/.441) if available
    fd = open("../../data/linux_p1550p2650.440", O_RDONLY);
    if (fd>=0){
        ephemeris_file_format_t f = assist_detect_ephemeris_file_format(fd);
        if (f != FILE_FORMAT_ASCII_BIN){
            fprintf(stderr, "Expected .440/.441 ASCII-derived binary. Got %d\n", (int)f);
            failures++;
        }
        close(fd);
    }

    // Create a small ASCII source-like temp file (unsupported; should be UNKNOWN)
    char tmp[] = "/tmp/assist_ascii_XXXXXX";
    int tfd = mkstemp(tmp);
    if (tfd>=0){
        const char* content = "JPL EPHEMERIS ASCII SOURCE\nConstants...\n";
        write(tfd, content, strlen(content));
        lseek(tfd, 0, SEEK_SET);
        ephemeris_file_format_t f = assist_detect_ephemeris_file_format(tfd);
        if (f != FILE_FORMAT_UNKNOWN){
            fprintf(stderr, "Expected UNKNOWN for ASCII source. Got %d\n", (int)f);
            failures++;
        }
        close(tfd);
        unlink(tmp);
    }

    if (failures){
        fprintf(stderr, "format_detection: %d failures\n", failures);
        return 1;
    }
    printf("format_detection: OK\n");
    return 0;
}



