// Clean up list

/**
 * This example demonstrates how to use ASSIST to integrate a test particle in 
 * the field of the Sun, planets, moon, and a set of massive asteroids, whose
 * positions come from JPL's DE440/441 ephemeris.
 */
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <math.h>
#include "rebound.h"
#include "assist.h"

void read_inputs(char *filename, double* tepoch, double* tstart, double* tend, double* tstep, 
		 int *geocentric,
		 double *epsilon,
		 int *n_particles,		 
		 double **instate,
		 particle_const **part_cons,
		 //double **part_cons,		 
		 int *n_var,
		 int **invar_part,		 
		 double **invar,
		 particle_const **invar_part_cons,		 		 
		 //double **invar_part_cons,		 
		 double **cov_mat);

/*
void read_inputs(char *filename, double* tepoch, double* tstart, double* tend, double* tstep, 
		 int *geocentric,
		 double *epsilon,
		 int *n_particles,		 
		 double **instate,
		 int *n_var,
		 int **invar_part,		 
		 double **invar,
		 double **cov_mat);
*/

int main(int argc, char* argv[]){

    // These variables are used to capture the results of read_inputs()
    double tepoch, tstart, tend, tstep;
    int geocentric;
    double epsilon;
    int n_particles;    
    double *instate;
    int n_var;    
    int *invar_part; // This stores the particle index that goes with the variational particle
    double *invar;
    particle_const* invar_part_cons;
    double* cov_mat;
    particle_const* part_const = NULL;

    //void heartbeat(struct reb_simulation* r);

    if(argc ==2){
	//read_inputs("holman_ic", &tepoch, &tstart, &tend, &tstep,
	read_inputs(argv[1], &tepoch, &tstart, &tend, &tstep,    
		    &geocentric, &epsilon,
		    &n_particles,		    
		    &instate,
		    &part_const,		    
		    &n_var,
		    &invar_part,
		    &invar,
		    &invar_part_cons,
		    &cov_mat);
    }else{
        printf("No Input File\n");
        exit(EXIT_FAILURE);
    }
	
    
    // This is for allocating memory to store the steps and
    // substeps of the integration
    int n_alloc = ((tend-tstart)/tstep + 1)*10;
    int nsubsteps = 10;
    double* outstate = NULL;
    double* outtime = NULL;
    outstate = (double *) malloc((8*n_alloc+1)*nsubsteps*6*sizeof(double));
    outtime  = (double *) malloc((8*n_alloc+1)*sizeof(double));

    //particle_const* part_cons = NULL;
    //part_cons = (particle_const*) malloc(n_particles*sizeof(particle_const));

    //double hg[11]   =   { 0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0};
    // Set substeps for output.
    double hg[nsubsteps+1];
    for(int i=0; i<=nsubsteps; i++){
	hg[i]=(1.0/nsubsteps)*i;
    }

    int n_steps;
    int status = integration_function(tstart, tend, tstep, 0, 1e-9,
				      n_particles, instate,
				      n_var, invar_part, invar,
				      part_const,
				      invar_part_cons,				      
				      n_alloc,
				      &n_steps,
				      nsubsteps,
				      hg,
				      outtime, outstate,
				      0.0); //, INFINITY);//  0.01, 32.0);

    int nouts = n_steps*nsubsteps + 1;

    for (int j=0; j<nouts; j++){
	for (int k=0; k<n_particles; k++){
	    int offset = j*(k+1)*6;
	    printf("%d %lf %.12lf %.12lf %.12lf %.12lf %.12lf %.12lf\n",
		   j, outtime[j],
		   outstate[offset+0], outstate[offset+1], outstate[offset+2],
		   outstate[offset+3], outstate[offset+4], outstate[offset+5]);
	}
    }
}

void read_inputs(char *filename, double* tepoch, double* tstart, double* tend, double* tstep, 
		 int *geocentric,
		 double *epsilon,
		 int *n_particles,		 
		 double **instate,
		 particle_const **inpart_cons,
		 //double **inpart_cons,		 
		 int *n_var,
		 int **invar_part,		 
		 double **invar,
		 particle_const **invar_part_cons,		 		 
		 //double **invar_part_cons,
		 double **cov_mat){

     char label[100]; /* hardwired for length */
     FILE* fp;

     int np = 0;
     int nvar = 0;

     int n_alloc = 1;
     double* state = malloc(n_alloc*6*sizeof(double));

     int n_var_alloc = 1;
     double* var = malloc(n_var_alloc*6*sizeof(double));
     int* var_part = malloc(n_var_alloc*sizeof(int));

     int var_flag = 0;

     particle_const* part_cons = malloc(n_alloc*sizeof(particle_const));
     //double* part_cons = malloc(n_alloc*3*sizeof(double));

     particle_const* var_part_cons = malloc(n_var_alloc*sizeof(particle_const));     

     double* cov = malloc(36*sizeof(double));

     int non_gravs = 0;
     double A1, A2, A3;
     
     if((fp = fopen(filename, "r")) != NULL){

      while(fscanf(fp, "%s", label) != EOF){
	  if(!strcmp(label, "tepoch")){
	      fscanf(fp, "%lf", tepoch);     
	  } else if(!strcmp(label, "tstart")){
	      fscanf(fp, "%lf", tstart);
	  } else if(!strcmp(label, "tstep")){
	      fscanf(fp, "%lf", tstep);
	  } else if(!strcmp(label, "tend")){
	      fscanf(fp, "%lf", tend);
	  } else if(!strcmp(label, "epsilon")){
	      fscanf(fp, "%lf", epsilon);
	  } else if(!strcmp(label, "geocentric")){
	      fscanf(fp, "%d", geocentric);
	  } else if(!strcmp(label, "state")){
	      fscanf(fp, "%lf%lf%lf", &state[6*np+0], &state[6*np+1], &state[6*np+2]);
	      fscanf(fp, "%lf%lf%lf", &state[6*np+3], &state[6*np+4], &state[6*np+5]);
	      // Set default values
	      part_cons[np].A1 = 0.0;
	      part_cons[np].A2 = 0.0;
	      part_cons[np].A3 = 0.0;
	      //part_cons[3*np+0] = 0.0;
	      //part_cons[3*np+1] = 0.0;
	      //part_cons[3*np+2] = 0.0;
	 
	      np++;

	      // Resize the array, if needed.
	      if(np==n_alloc){
		  n_alloc *= 2;
		  state = realloc(state, n_alloc*6*sizeof(double));
		  part_cons = realloc(part_cons, n_alloc*sizeof(particle_const));
		  //part_cons = realloc(part_cons, n_alloc*3*sizeof(double));
	      }

	  } else if(!strcmp(label, "A1A2A3")){

	      printf("A1A2A3\n");
	      fflush(stdout);
	      
	      fscanf(fp, "%lf %lf %lf", &A1, &A2, &A3);
	      if(var_flag){

		  printf("ngvar\n");
		  fflush(stdout);
		  var_part_cons[nvar-1].A1 = A1;
		  var_part_cons[nvar-1].A2 = A2;
		  var_part_cons[nvar-1].A3 = A3;
	      }else{
		  part_cons[np-1].A1 = A1;
		  part_cons[np-1].A2 = A2;
		  part_cons[np-1].A3 = A3;
	      }
		  
	      //part_cons[3*(np-1)+0] = A1;
	      //part_cons[3*(np-1)+1] = A2;
	      //part_cons[3*(np-1)+2] = A3;	      

	      if(A1!=0. || A2!=0. || A3!=0.){
		  non_gravs = 1;
	      }

	      // Check if any of the values is non-zero.
	      // If so, then prepare to pass the array.
	 
	  } else if(!strcmp(label, "var")){
	      var_flag=1;

	      printf("var here\n");
	      fflush(stdout);
		
	      fscanf(fp, "%d", &var_part[nvar]);
	      fscanf(fp, "%lf%lf%lf", &var[6*nvar+0], &var[6*nvar+1], &var[6*nvar+2]);
	      fscanf(fp, "%lf%lf%lf", &var[6*nvar+3], &var[6*nvar+4], &var[6*nvar+5]);	 
	      nvar++;

	      printf("var here 2\n");
	      fflush(stdout);
	      

	      // Resize the array, if needed.
	      if(nvar==n_var_alloc){
		  
		  n_var_alloc *= 2;
		  var = realloc(var, n_var_alloc*6*sizeof(double));
		  var_part = realloc(var_part, n_var_alloc*sizeof(int));
		  var_part_cons = realloc(var_part_cons, n_var_alloc*sizeof(particle_const));

		  printf("var there\n");
		  fflush(stdout);
		  
	      }
	 
	  } else if(!strcmp(label, "covariance")){
	      fscanf(fp, "%lf%lf%lf%lf%lf%lf", &cov[0], &cov[1], &cov[2],&cov[3], &cov[4], &cov[5]);
	      fscanf(fp, "%lf%lf%lf%lf%lf%lf", &cov[6], &cov[7], &cov[8],&cov[9], &cov[10], &cov[11]);
	      fscanf(fp, "%lf%lf%lf%lf%lf%lf", &cov[12], &cov[13], &cov[14],&cov[15], &cov[16], &cov[17]);
	      fscanf(fp, "%lf%lf%lf%lf%lf%lf", &cov[18], &cov[19], &cov[20],&cov[21], &cov[22], &cov[23]);
	      fscanf(fp, "%lf%lf%lf%lf%lf%lf", &cov[24], &cov[25], &cov[26],&cov[27], &cov[28], &cov[29]);
	      fscanf(fp, "%lf%lf%lf%lf%lf%lf", &cov[30], &cov[31], &cov[32],&cov[33], &cov[34], &cov[35]);
 
	  } else {
	      printf("No label: %s\n", label);
	      exit(EXIT_FAILURE);
	  }
      }

      printf("%d\n", non_gravs);
      fflush(stdout);
      
      if(non_gravs){
	  for(int i=0; i<np; i++){
	      printf("A1A2A3: %lf %lf %lf\n", part_cons[i].A1, part_cons[i].A2, part_cons[i].A3);
	      //printf("A1A2A3: %lf %lf %lf\n", part_cons[3*i+0], part_cons[3*i+1], part_cons[3*i+2]);	      
	  }

	  part_cons = realloc(part_cons, np*sizeof(particle_const));
	  //part_cons = realloc(part_cons, np*3*sizeof(double));	  
	  *inpart_cons = part_cons;
	  //var_part_cons = realloc(var_part_cons, np*3*sizeof(double));
	  var_part_cons = realloc(var_part_cons, nvar*sizeof(particle_const));	  	  
	  *inpart_cons = part_cons;
	  *invar_part_cons = var_part_cons;
	  
      }else{
	  free(part_cons);
	  *inpart_cons = NULL;
      }	  
      
      //deallocate unused space
      state = realloc(state, np*6*sizeof(double));
      var   = realloc(var, nvar*6*sizeof(double));
      var_part   = realloc(var_part, nvar*6*sizeof(int));            

      *n_particles = np;
      *n_var = nvar;      
      *instate = state;
      *invar = var;
      *invar_part = var_part;            
      *cov_mat = cov;

      fclose(fp);

     }else{
	 printf("no such file %s\n", filename); 
	 exit(EXIT_FAILURE);
     }

     return;
}




