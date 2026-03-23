#include <stdio.h>
#include <math.h>
#include "assist.h"
#include "ascii_ephem.h"

static void assert_range(const char* name, double v, double lo, double hi, int* failures){
    if (!(v>lo && v<hi)){
        fprintf(stderr, "%s out of range: %g\n", name, v);
        (*failures)++;
    }
}

int main(){
    int failures = 0;
    struct ascii_s* a = assist_ascii_init("../../data/linux_p1550p2650.440");
    if (!a){
        fprintf(stderr, "Failed to load ASCII-derived binary ephemeris (.440/.441).\n");
        return 1;
    }
    assert_range("AU", a->AU, 1.49e8, 1.50e8, &failures);
    assert_range("CLIGHT", a->CLIGHT, 2.9e5, 3.1e5, &failures);
    assert_range("EMRAT", a->cem, 70.0, 100.0, &failures);

    double GM, x,y,z, vx,vy,vz, ax,ay,az;
    int flag = assist_ascii_calc(a, 2451545.0, 0.0, ASSIST_BODY_SUN, &GM, &x,&y,&z, &vx,&vy,&vz, &ax,&ay,&az);
    if (flag!=ASSIST_SUCCESS || !isfinite(x) || !isfinite(y) || !isfinite(z)){
        fprintf(stderr, "assist_ascii_calc failed or non-finite outputs.\n");
        failures++;
    }
    assist_ascii_free(a);

    if (failures){
        fprintf(stderr, "ascii_init_constants: %d failures\n", failures);
        return 1;
    }
    printf("ascii_init_constants: OK\n");
    return 0;
}



