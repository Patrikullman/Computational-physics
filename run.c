#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "lattice.h"
#include "potential.h"
#include "lattice.h"
#include "potential.h"

// these allows us to use gsl functions. 
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <time.h>

int
run(
    int argc,
    char *argv[]
   )
{
/////using init_fcc to create the fcc lattice.

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
init_fcc(position_array, N, a0); //Function takes a matrix of size [4*N*N*N][3] as input and stores a fcc lattice in it. 

for(int i=0; i<ROWS; i++){
	for(int j=0; j<COLS;j++){
		printf("position_array %f\n", position_array[i][j]);  //just printing out the fcc lattice to see.
		}
	}
/////Below we are calculating E_pot for various lattice parameters and plotting them in python.
double E_pot;
double a0_array[] = {3.95,3.96,3.97,3.98,3.99,4.0,4.1,4.2}; //unit Ångström.
double natoms = 256;
    for (int i = 0; i<8;i++){     //Calculate the lattice and energy for many structures, i.e varying a0 so we can easily plot it and make a fit in python. 
            init_fcc(position_array, N, a0_array[i]); 
            E_pot = get_energy_AL(position_array, N * a0_array[i], natoms)/64; //divide by 64 because we need to divide by N*N*N = 4*4*4 = 64
    	    FILE *FP;
            FP = fopen("lattice_energies.csv", "a"); //Once we get the a0 and corresponding E_pot's  we write it to a file and plot E_pot as a function of a0_array in python
            fprintf(FP, "%f,%f\n", a0_array[i],E_pot); 
          }
          

//////now, we want to introduce small deviations in the lattice. 


const gsl_rng_type * T;
    gsl_rng * r;

    int i, n = 256;
    gsl_rng_env_setup();

    T = gsl_rng_default;
    r = gsl_rng_alloc(T); 

    // typically the seed is set to the current time.
    //time_t seed = time(NULL);
    int seed = 42;
    gsl_rng_set(r, seed);
	double* ptr; //initializing array ptr to store our random uniform numbers. 
    for (i = 0; i < ROWS; i++){   //for loop for our uniform numbers
	double u = gsl_rng_uniform(r); //uniform numbers are stored in variable u
	ptr = (double*)malloc(1000 * sizeof(double)); //dynamically allocating an array, to store all u values.
	ptr[i] = u*(0.26); //saving all u values in the ptr array. (this is a uniform array between 0 and 0.26 that i can use to deviate my atom positions.)
	for (int i = 0;i<ROWS;i++){
		for(int j=0;j<COLS;j++){ //here im looping through each element, giving it 50% chance to the position getting either + och - 6.5% deviation. 
			if (ptr[i] <= 0.26/2) {  //this if-statement creates a 50% chance.  (0.26 comes from 6.5% of a0=4 as lattice parameter)
				position_array[i][j] = position_array[i][j] + ptr[i]; 
				printf("p %f\n",position_array[i][j]);
				}
			else {
			position_array[i][j] = position_array[i][j] - ptr[i];
}
}
}
	//position_array[i][0] = position_array[i][0] + ptr[i];  //just introducing a deviation +- 0.26 in the original positions. 
	//position_array[i][1] = position_array[i][1] - ptr[i];  //this works to give diviations in the lattice but its not completely random. 
	//position_array[i][2] = position_array[i][2] + ptr[i];
	//printf(" %f\n %f\n %f\n",position_array[i][0],position_array[i][1],position_array[i][2]);

    }
    gsl_rng_free (r);

double pos_x[ROWS][COLS]; 
double pos_y[ROWS][COLS];
double pos_z[ROWS][COLS];
for(int i = 0; i<ROWS;i++){
	pos_x[i][0] = position_array[i][0];
	pos_y[i][1] = position_array[i][1];
	pos_z[i][2] = position_array[i][2];
}

for(int i=0; i<ROWS; i++){
	for(int j=0; j<COLS;j++){
	//printf("%f\n", array[i][j]);
	FILE *fp;
	fp = fopen("lattice_positions.csv", "a");  // writing positions of the fcc atoms to a file.
	fprintf(fp, "%f,%f,%f\n", pos_x[i][0],pos_y[i][1],pos_z[i][2] ); 
	}
	}






    return 0;
}
