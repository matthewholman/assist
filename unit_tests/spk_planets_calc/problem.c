#include <stdio.h>
#include <math.h>
#include "assist.h"
#include "spk.h"

static int finite3(double a, double b, double c){
    return isfinite(a) && isfinite(b) && isfinite(c);
}

int main(){
    int failures = 0;
    struct spk_s* sp = assist_spk_init("../../data/de440.bsp");
    if (!sp){
        fprintf(stderr, "Failed to load SPK de440.bsp\n");
        return 1;
    }
    struct spk_constants_and_masses data = assist_load_spk_constants_and_masses("../../data/de440.bsp");
    struct assist_ephem ephem = (struct assist_ephem){0};
    assist_apply_spk_constants(&ephem, &data);
    ephem.spk_planets = sp;

    double GM,x,y,z,vx,vy,vz,ax,ay,az;
    int flag = assist_spk_calc_planets(&ephem, 2451545.0, 0.0, 10, &GM, &x,&y,&z, &vx,&vy,&vz, &ax,&ay,&az);
    if (flag!=ASSIST_SUCCESS || !finite3(x,y,z) || !finite3(vx,vy,vz)){
        fprintf(stderr, "assist_spk_calc_planets failed or non-finite for Sun.\n");
        failures++;
    }
    flag = assist_spk_calc_planets(&ephem, 2451545.0, 0.0, 399, &GM, &x,&y,&z, &vx,&vy,&vz, &ax,&ay,&az);
    if (flag!=ASSIST_SUCCESS){
        fprintf(stderr, "assist_spk_calc_planets failed for Earth.\n");
        failures++;
    }
    flag = assist_spk_calc_planets(&ephem, 2451545.0, 0.0, 301, &GM, &x,&y,&z, &vx,&vy,&vz, &ax,&ay,&az);
    if (flag!=ASSIST_SUCCESS){
        fprintf(stderr, "assist_spk_calc_planets failed for Moon.\n");
        failures++;
    }

    assist_free_spk_constants_and_masses(&data);
    assist_spk_free(sp);

    if (failures){
        fprintf(stderr, "spk_planets_calc: %d failures\n", failures);
        return 1;
    }
    printf("spk_planets_calc: OK\n");
    return 0;
}
