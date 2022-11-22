/******************************************************************************
 * E1code4
 ******************************************************************************
 * Routine that runs the velocity verlet algorithm
 * Use as template to construct your program!
 */

/*
 * Calculate the acceleration
 * @a - vector that is filled with acceleration
 * @u - vector with the current positions
 * @m - vector with masses
 * @kappa - Spring constant
 * @size_of_u - the size of the position, acceleration and mass array
 */
#include <stdio.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "lattice.h"
#include "potential.h"
#include "lattice.h"
#include "potential.h"
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <time.h>

void calc_acc(double *a, double *u, double *m, double kappa, int size_of_u)
{
    /* Declaration of variables */
    int i;
    
    /* Calculating the acceleration on the boundaries */
    a[0] = kappa*(- 2*u[0] + u[1])/m[0];
    a[size_of_u - 1] = kappa*(u[size_of_u - 2] - 2*u[size_of_u - 1])/m[size_of_u - 1];
    
    /* Calculating the acceleration of the inner points */
    for (i = 1; i < size_of_u - 1; i++){
        a[i] = kappa*(u[i - 1] - 2*u[i] + u[i + 1])/m[i];
    }
}

/*
 * Perform the velocity verlet alogrithm 
 * @n_timesteps - The number of time steps to be performed
 * @n_particles - number of particles in the system
 * @v - array of velocity (Empty allocated array) : sizeof(v) = n_particles
 * @q_n - position of the n'th atom : sizeof(q_n) = n_timesteps+1
 * @dt - timestep
 * @m - vector with masses of atoms sizeof(n_particles)
 * @kappa - Spring constant
 */

void velocity_verlet(int n_timesteps, int n_particles, double *v,double* position_array,double dt, double *m,
		     double kappa)
{
    double q[n_particles][3]; //måste generalisera så att den kan hålla alla atomer i lattice,  256x3
    double a[n_particles][3]; //samma här,  256x3
    
    for(int i =0;i<n_particles){
    	for(int j = 0;j<3;j++){ //den här nested loopen är till för att spara alla startvärden i q[i][j]
    		q[i][j] = position_array[i][j];  //försöker spara mina atom positioner i q matrisen så att dem kan uppdateras enligt velocity verlet nedan.
    		}
    	}
    calc_acc(a, q, m, kappa, n_particles);  //vi måste generalisera accelerationsfunktionen till 3 dimensioner också.
    for (int i = 1; i < n_timesteps + 1; i++) {
        /* v(t+dt/2) */
        for (int j = 0; j < n_particles; j++) {
        	for(k=0; k<COLS;k++){ //vi måste loopa igenom varenda atom i våran lattice, [256][3]
            	    v[j][k] += dt * 0.5 * a[j]; //uppdaterar hastghet för alla 256x3 atomer (a[j] måste generaliseras till 3 dimensioner)
        }
        }
        
        /* q(t+dt) */
        for (int j = 0; j < n_particles; j++) {
        	for(k=0; k<COLS;k++){
            	    q[j][k] += dt * v[j][k]; //uppdaterar position för alla 256x3 atomer
            
        }
        }
        /* a(t+dt) */
        calc_acc(a, q, m, kappa, n_particles);
        
        /* v(t+dt) */
        for (int j = 0; j < n_particles; j++) {
            for(k=0; k<COLS;k++){
                    v[j][k] += dt * 0.5 * a[j]; 

        }
        /* Save the displacement of the 256 atoms */
        position_array[i][0] = q[i][0];
        position_array[i][1] = q[i][1];  //lattice uppdaterade atompositioner sparas för varje iteration/timestep
        position_array[i][2] = q[i][2];     

    }    
}


int main() {

int n_particles = 256;

#define ROWS 256
#define COLS 3
	double position_array[ROWS][COLS]; //defining a 2D array to fill with the fcc lattice. 
	
	for(int i=0; i<ROWS; i++){
		for(int j=0; j<COLS;j++){
		position_array[i][j] = 0;
	}
	}
int N = 4; //number of unit cells
double a0 = 4; //lattice parameter
init_fcc(position_array, N, a0); //skapar FCC lattice så att vi kan lägga den i velocity_verlet().

	double v[ROWS][COLS]; //defining a 2D empty velocity array that we can use in the 3D velocity verlet algoritm.
	for(int i=0; i<ROWS; i++){
		for(int j=0; j<COLS;j++){
		v[i][j] = 0;
	}
	}
	double m[ROWS]; //mass array 256x1 array.
	for(int i=0; i<ROWS; i++){
		m[i] = 1;
	}
		double kappa[ROWS]; //kappa array 256x1 array
	for(int i=0; i<ROWS; i++){
		kappa[i] = 1; //fyll kappa array med kappa värden..
	}

velocity_verlet(1000,256,v,position_array,25,m,kappa);
//v är en 2d array av 256x3 nu.
//m är en array av 256x1 nu.
//kappa är en array av 256x1 nu.

   FILE *fp = fopen("velocityverlet_FCC_lattice.csv", "w");
    for(int i = 0; i < 256; ++i){
    	for(int j = 0; j<3;j++){
	    fprintf(fp, "%f\n", position_array[i][j]);
    }
    }
    fclose(fp);
return 0;



}








