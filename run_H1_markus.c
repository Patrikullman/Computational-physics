#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include "lattice.h"
#include "potential.h"
#include "tools.h"



int
run(
    int argc,
    char *argv[]
   ) 
{
    int N = 4; // Number of unit cells in each direction

    /* Creating array with slightly varying lattice constants */
    double a0_arr[10] = {3.95, 3.97, 3.98, 3.99, 4, 4.01, 4.03, 4.05, 4.08, 4.09}; 

    
    double position_array[4*N*N*N][3]; // initializing an array with coordinates/positions of the atoms
    int natoms = sizeof(position_array) / (3 * sizeof(double)); // number of atoms 
    
    FILE *fp = fopen("lattice_para_vs_Ep.csv", "w");
    for(int i = 0; i < 10; ++i){
        
        init_fcc(position_array, N, a0_arr[i]); 
        double E_pot = get_energy_AL(position_array, N * a0_arr[i], natoms) / (N*N*N); 
        fprintf(fp, "%f, %f\n", a0_arr[i], E_pot); 
    }
    fclose(fp);

    FILE *fp2 = fopen("positions.csv", "w");
    for(int i = 0; i < 256; ++i){
        fprintf(fp2, "%f, %f, %f\n", position_array[i][0], position_array[i][1], position_array[i][2]);
    }
    fclose(fp2);



    /* Now follows some code used to study time evolution of energies */

    double m = 26*0.103625e-3; // Mass of Al in asu
    double a0 = 4.0303522; // The lattice constant
    int n_timesteps = 1000;
    double dt = 0.01;
    double velocity_array[4*N*N*N][3]; 
    
    /* Initialize the velocity array with all velocities equal to 0 */
    for(int i = 0; i < 2*N; ++i){
        for(int j = 0; j < 2*N; ++j){
            for(int k = 0; k < N; ++k){
                velocity_array[i * N * 2 * N + j *N + k][0] = 0;
                velocity_array[i * N * 2 * N + j *N + k][1] = 0;
                velocity_array[i * N * 2 * N + j *N + k][2] = 0;
            }
        }
    }

    double pos_array[4*N*N*N][3]; // Array with positions of each atom in ONE lattice
    double pos_evolution[n_timesteps][4*N*N*N][3]; // Array for time evolution of lattice
    double force[4*N*N*N][3]; // Array with the forces between the atoms in the lattice
    init_fcc(pos_array, N, a0); // Initialize the lattice
    
    /* Code to give the atoms a deviation from the equilibrium */
    const gsl_rng_type * T;
    gsl_rng * r;

    gsl_rng_env_setup();

    T = gsl_rng_default;
    r = gsl_rng_alloc(T);

    int seed = 42;
    gsl_rng_set(r, seed);

    for(int i = 0; i < 2*N; ++i){
        for(int j = 0; j < 2*N; ++j){
            for(int k = 0; k < N; ++k){
                pos_array[i * N * 2 * N + j *N + k][0] += 2*(gsl_rng_uniform(r)-0.5)*a0*0.065;
                pos_array[i * N * 2 * N + j *N + k][1] += 2*(gsl_rng_uniform(r)-0.5)*a0*0.065;
                pos_array[i * N * 2 * N + j *N + k][2] += 2*(gsl_rng_uniform(r)-0.5)*a0*0.065;
            }
        }
    }

    FILE *fp3 = fopen("positions_deviated.csv", "w");
    for(int i = 0; i < 256; ++i){
        fprintf(fp3, "%f, %f, %f\n", position_array[i][0], position_array[i][1], position_array[i][2]);
    }
    fclose(fp3);

    for(int i = 0; i < n_timesteps; ++i){
        for(int j = 0; j < 2*N; ++j){
            for(int k = 0; k < 2*N; ++k){
                for(int l = 0; l < N; ++l){

                    pos_evolution[i][j * N * 2 * N + k * N + l][0] = pos_array[j * N * 2 * N + k * N + l][0];
                    pos_evolution[i][j * N * 2 * N + k * N + l][1] = pos_array[j * N * 2 * N + k * N + l][1];
                    pos_evolution[i][j * N * 2 * N + k * N + l][2] = pos_array[j * N * 2 * N + k * N + l][2];
                    //printf("pos ev:%f\n", pos_evolution[i][j * N * 2 * N + k * N + l][0]);
                }
            }
        }
    }

    double *Ek = (double*)malloc(sizeof(double)*(n_timesteps+ 1));
    velocity_verlet(n_timesteps, 4*N*N*N, velocity_array, pos_array, pos_evolution, dt, m, force, N*a0, Ek);
    //printf("Ek:%f\n", Ek[250]);

    //printf("pos:%f\n", pos_array[250][1]);

    double *Ep = (double*)malloc(sizeof(double)*(n_timesteps+1));
    //double *Ek = (double*)malloc(sizeof(double)*n_timesteps);
    //double Ek[n_timesteps];
    double *Etot = (double*)malloc(sizeof(double)*(n_timesteps+1));
    double *times = (double*)malloc(sizeof(double)*(n_timesteps+1));

    for(int i = 0; i < n_timesteps; ++i){
       // printf("pos_ev:%f\n", pos_evolution[i][7][0]);
    }

    for(int i = 0; i < n_timesteps; ++i){
        //printf("pos_ev:%f\n", pos_evolution[i][7][0]);
        Ep[i] = get_energy_AL(pos_evolution[i], N*a0, 4*N*N*N); // /(N*N*N);
        //printf("Ep:%f\n", Ep[i]);
        //Ek[i] = get_virial_AL(pos_evolution[i], N*a0, 4*N*N*N)/(N*N*N);
        Etot[i] = Ep[i] + Ek[i];
        times[i] = i*dt;
        //printf("Ek, Ep:%f,%f\n", Ek[i], Ep[i]);
    }
    
    write_to_file("energy_arrays.csv", times, Ep, Ek, Etot, n_timesteps);


   
    double kB = 1.380649e-23; // Boltzmann's constant
    double temp_inst = 0; // The instantanous temperature
    double temp = 0; // The final temperature

    /* Loop through and add up the inst temperatures, finally calculate average*/
    for(int i = 0; i < n_timesteps; ++i){
        temp_inst += (2*Ek[i]*1.60218e-19) / (3 * kB * (4*N*N*N) );
    }
    temp = temp_inst / n_timesteps;
    printf("T: %f\n", temp);

    

    /* Now follows some code used for equilibrating the system at specified temp and pressure */
    
    double *Ek_scaled = (double*)malloc(sizeof(double)*(n_timesteps+1));
   // double *C_V = (double*)malloc(n_timesteps * sizeof(double));
    double *temp_inst_verlet_scaled = (double*)malloc(n_timesteps * sizeof(double));

    
    velocity_verlet_scaled(n_timesteps, 4*N*N*N, velocity_array, pos_array, pos_evolution, dt, m, force, N*a0, Ek_scaled, temp_inst_verlet_scaled);
    
    double *Ep_scaled = (double*)malloc(sizeof(double)*(n_timesteps+1));
    double *Etot_scaled = (double*)malloc(sizeof(double)*(n_timesteps+1));
    double *times_scaled = (double*)malloc(sizeof(double)*(n_timesteps+1));
    double Ek_mean = 0;
    double Ek_vari = 0;

    for(int i = 600; i < n_timesteps; ++i){
        Ek_mean += Ek_scaled[i];
        // Ek_vari += (Ek_scaled[i]*Ek_scaled[i]);
    }
    Ek_mean = Ek_mean / 400;

    double* EK_centered = (double*)malloc(400*sizeof(double));
    for(int i = 600; i < n_timesteps; ++i){
        EK_centered[i-600] = Ek_scaled[i-600]-Ek_mean;
    }
    for(int i = 600; i < n_timesteps; ++i){
        Ek_vari += EK_centered[i-600]*EK_centered[i-600];
    }
    Ek_vari =  Ek_vari/400;

    FILE *fp4 = fopen("particle_evol_x.csv", "w");
    FILE *fp5 = fopen("heat_capa.csv", "w");
    for(int i = 0; i < n_timesteps; ++i){
       
        Ep_scaled[i] = get_energy_AL(pos_evolution[i], N*a0, 4*N*N*N);
        Etot_scaled[i] = Ep[i] + Ek[i];
        times_scaled[i] = i*dt;
        fprintf(fp4, "%f,%f,%f,%f,%f\n", times[i], pos_evolution[i][0][0], pos_evolution[i][101][0], 
        pos_evolution[i][150][0], pos_evolution[i][199][0]);

       
        
        
    }
    
    fclose(fp4);
    fclose(fp5);
    double kBeV = 8.617e-5;
    double C_V = (3*(4*N*N*N)*kBeV / 2) /
        (1 - (2 / (3*(4*N*N*N)*kBeV*kBeV*temp_inst_verlet_scaled[n_timesteps-1]*temp_inst_verlet_scaled[n_timesteps-1])) * (Ek_vari));

    printf("Ek vari: %0.30f, CV: %0.30f\n", Ek_vari, C_V);

    double temp_inst_scaled = 0; // The instantanous temperature
    double temp_scaled = 0; // The final temperature

    /* Loop through and add up the inst temperatures, finally calculate average*/
    for(int i = 0; i < n_timesteps; ++i){
        temp_inst_scaled += (2*Ek_scaled[i]*1.60218e-19) / (3 * kB * (4*N*N*N) );
    }
    temp_scaled = temp_inst_scaled / n_timesteps;
    printf("T scaled: %f\n", temp_scaled);
    //printf("x0: %f, x1: %f, x2: %f, x3: %f", pos_evolution[0][0][0], pos_evolution[0][101][0], 
       // pos_evolution[0][150][0], pos_evolution[0][199][0]);

    double** distance_array = (double**)malloc(256 * sizeof(double*));
    for (int i = 0; i < 256; i++){
        distance_array[i] = (double*)malloc(256 * sizeof(double));
    }

    FILE *fp6 = fopen("dist_between_points.csv", "w");
    for(int i = 0; i < 256; ++i){
        for(int j = 0; j < 256; ++j){
            //if(j = i){
              //distance_array[i][j] = 0;
            //}
            //else{
                distance_array[i][j] = distance(pos_evolution[n_timesteps-1][i], pos_evolution[n_timesteps-1][j]);
            //}
            fprintf(fp6, "%d,%f\n", i, distance_array[i][j]);
            //printf("dist: %f\n", distance_array[i][j]);
        }
    }
    fclose(fp6);

    return 0;
}