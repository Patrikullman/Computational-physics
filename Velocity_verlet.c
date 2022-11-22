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
//
void velocity_verlet(int n_timesteps, int n_particles, double *v,double* position_array,double* forces,double dt, double *m,
		     double kappa)
{
    double q[n_particles][3]; //måste generalisera så att den kan hålla alla atomer i lattice,  256x3
    double a[n_particles][3]; //samma här,  256x3
   
    for(int i =0;i<n_particles){
    	for(int j = 0;j<3;j++){ //den här nested loopen är till för att spara alla startvärden i q[i][j]
    		q[i][j] = position_array[i][j];  //försöker spara mina atom positioner i q matrisen så att dem kan uppdateras enligt velocity verlet nedan.
    		}
   
    	}
    	get_forces_AL(forces,position_array,cell_lenght, n_particles) //basically just replacing acceleration function call with get_forces_Al().
    for (int i = 1; i < n_timesteps + 1; i++) {
        /* v(t+dt/2) */
        for (int j = 0; j < n_particles; j++) {
        	for(k=0; k<COLS;k++){ //vi måste loopa igenom varenda atom i våran lattice, [256][3]
            	    v[j][k] += dt * 0.5 * forces[j][k]/m; //lägg in a = F/m här.
        }
        }
        /* q(t+dt) */
        for (int j = 0; j < n_particles; j++) {
        	for(k=0; k<COLS;k++){
            	    q[j][k] += dt * v[j][k]; //uppdaterar position för alla 256x3 atomer
        
        }
        }
        /* get the forces instead of the accelerations */
        get_forces_AL(forces,position_array,cell_lenght, n_particles)
        
        /* v(t+dt) */
        for (int j = 0; j < n_particles; j++) {
            for(k=0; k<COLS;k++){
                    v[j][k] += dt * 0.5 * forces[j][k]/m;; 

        }
        /* Save the displacement of the 256 atoms */
        position_array[i][0] = q[i][0];
        position_array[i][1] = q[i][1];  //lattice uppdaterade atompositioner sparas för varje iteration/timestep
        position_array[i][2] = q[i][2];     

    }
}
int main() {
double m = 1; //what is the mass of Al in correct units?
double kappa = 1; //what is kappa in correct units?
int n_particles = 256;
//Below im defining the vectors and matrices that will go into the verlet algoritm
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
init_fcc(position_array, N, a0); //constructing the FCC lattice. 

	double v[ROWS][COLS]; //defining an empty 2D velocity array that we can use in the velocity verlet algoritm.
	for(int i=0; i<ROWS; i++){
		for(int j=0; j<COLS;j++){
		v[i][j] = 0;
	}
	}
	double forces[ROWS][COLS]; //defining a 2D array of size [256][3] to fill with forces, by using "get_forces_AL()" function from H1potential.c
	for(int i=0; i<ROWS; i++){
		for(int j=0; j<COLS;j++){
		position_array[i][j] = 0;
	}
	}

velocity_verlet(1000,n_particles,v,position_array,forces,dt,m,kappa);

//once we get the algoritm working,we can extract velocities and positions to study the time-evolution of E_pot, E_kin and E_tot. 
   FILE *fp = fopen("velocityverlet_FCC_lattice.csv", "w");
    for(int i = 0; i < 256; ++i){
    	for(int j = 0; j<3;j++){
	    fprintf(fp, "%f\n",position_array[i][j]);
    }
    }
    fclose(fp);
return 0;



}








