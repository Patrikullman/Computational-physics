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

#include <stdlib.h>
#include <stdio.h>

void calc_acc(double *a, double *u, double *m, double kappa, int size_of_u)
{
    /* Declaration of variables */
    int i;
    
    /* Calculating the acceleration on the boundaries */
    a[0] = kappa*(- 1*u[0] + u[1])/m[0];
    a[size_of_u - 1] = kappa*(u[size_of_u - 2] - 1*u[size_of_u - 1])/m[size_of_u - 1];
    
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
void velocity_verlet(int n_timesteps, int n_particles, double *v, double *q_1,
		     double *q_2, double *q_3, double dt, double *m,
		     double kappa, double *v1, double *v2, double *v3)
{
    double q[n_particles];
    double a[n_particles];
    

    q[0] = q_1[0];
    q[1] = q_2[0];
    q[2] = q_3[0];

    v[0] = v1[0];
    v[1] = v2[0];
    v[2] = v3[0];

    calc_acc(a, q, m, kappa, n_particles);
    for (int i = 1; i < n_timesteps + 1; i++) {
    
        /* v(t+dt/2) */
        for (int j = 0; j < n_particles; j++) {
            v[j] += dt * 0.5 * a[j];
        }
        
        /* q(t+dt) */
        for (int j = 0; j < n_particles; j++) {
            q[j] += dt * v[j];
        }
        
        /* a(t+dt) */
        calc_acc(a, q, m, kappa, n_particles);
        
        /* v(t+dt) */
        for (int j = 0; j < n_particles; j++) {
            v[j] += dt * 0.5 * a[j];
        }
        
        /* Save the displacement of the three atoms */
        q_1[i] = q[0];
        q_2[i] = q[1];
        q_3[i] = q[2];

        /* Save the velocity of the three atoms */
        v1[i] = v[0];
        v2[i] = v[1];
        v3[i] = v[2];

    }
}

void write_to_file(char *fname, double *frequencies,
		   double *spectrum, int n_points)
{
    FILE *fp = fopen(fname, "w");
    fprintf(fp, "time, signal\n");
    for(int i = 0; i < n_points; ++i){
	    fprintf(fp, "%f,%f\n", frequencies[i], spectrum[i]);
    }
    fclose(fp);
}

int main()
{
    /* Defining constants and parameters */
    int n_timesteps = 1000; 
    int n_particles = 3;
    double dt = 0.25/n_timesteps;
    double kappa = 1600/16.0218; // kappa expressed in terms of asu 
    int row_size = n_timesteps; 
    int column_size = n_particles+2; // add 2 because of the walls
    double m[3] = {16.0/9649, 12.0/9649, 16.0/9649}; // mass in asu

    /* Allocating memory for position arrays and storing initial position */
    double* q_1 = (double*)malloc(n_timesteps * sizeof(double));  
    double* q_2 = (double*)malloc(n_timesteps * sizeof(double)); 
    double* q_3 = (double*)malloc(n_timesteps * sizeof(double));
    q_1[0] = 0.01;
    q_2[0] = 0.005;
    q_3[0] = -0.005;

    /* Allocating memory for velocity arrays and storing initial velocities for each particle */
    double* v = (double*)malloc(n_particles * sizeof(double));
    for (int i = 0; i < n_particles; ++i){
        v[i] = 0.;
    }
    
    double* v1 = (double*)malloc(n_timesteps * sizeof(double));  
    double* v2 = (double*)malloc(n_timesteps * sizeof(double)); 
    double* v3 = (double*)malloc(n_timesteps * sizeof(double));
    v1[0] = 0.; //initial velocities
    v2[0] = 0.;
    v3[0] = 0.;

   
    

    /* Running velocity verlet algorithm */
    velocity_verlet(n_timesteps, n_particles, v, q_1, q_2, q_3, dt, m, kappa, v1, v2, v3);
        
    /* Writing the positions at each time step to file */
    FILE *fp = fopen("velocityverlet_CO2.csv", "w");
    for(int i = 0; i < n_timesteps; ++i){
        fprintf(fp, "%f,%f,%f\n", q_1[i], q_2[i], q_3[i]); 
    }
    fclose(fp);

    /* Allocating memory for 2D arrays to store position and velocities */
    double** position_array = (double**)malloc(row_size * sizeof(double*));
    for (int i = 0; i < row_size; i++)
        position_array[i] = (double*)malloc(n_particles * sizeof(double));

    double** velocity_array = (double**)malloc(row_size * sizeof(double*));
    for (int i = 0; i < row_size; i++)
        velocity_array[i] = (double*)malloc(n_particles * sizeof(double)); 
    

    /* Storing the positions and velocities in 2D arrays */
    for(int i = 0; i < n_timesteps; ++i){
            //position_array[i][0] = 0;
            position_array[i][1] = q_1[i];
            position_array[i][2] = q_2[i];
            position_array[i][3] = q_3[i];
            //position_array[i][4] = 0;

            velocity_array[i][0] = v1[i];
            velocity_array[i][1] = v2[i];
            velocity_array[i][2] = v3[i];
    }
    

    /* Allocating memory to store energy data */
    double* Ek = (double*)malloc(n_timesteps * sizeof(double));
    double* Ep = (double*)malloc(n_timesteps * sizeof(double));
    double* Etot = (double*)malloc(n_timesteps * sizeof(double));

    /* Calculating energy for each time step and writing data to file */
    FILE *fp2 = fopen("energy_data_CO2.csv", "w");
    for (int i = 0; i < n_timesteps; ++i){

            Ep[i] += kappa/2 * (
                (position_array[i][2]-position_array[i][1])*(position_array[i][2]-position_array[i][1])+
                (position_array[i][3]-position_array[i][2])*(position_array[i][3]-position_array[i][2]));

        
        for(int j = 0; j < 3; ++j){
            Ek[i] += m[j] * velocity_array[i][j] * velocity_array[i][j] / 2;
        }
        
        Etot[i] = Ek[i] + Ep[i];
        fprintf(fp2, "%f,%f,%f\n", Ep[i], Ek[i], Etot[i]);
    }
    fclose(fp2);


    return 0;
}
