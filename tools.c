#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "potential.h"


void velocity_verlet(int n_timesteps, int n_particles, double v[][3], double position[][3],
double position_evolution[][n_particles][3], double dt, double m, double force[][3], 
double cell_length, double E_kin[n_timesteps])
{
    for(int i = 1;  i < n_timesteps + 1; ++i){
        int N = 4;
        for(int j = 0; j < 2*N; ++j){
            for(int k = 0; k < 2*N; ++k){
                for(int l = 0; l < N; ++l){

                    
                    v[j * N * 2 * N + k * N + l][0] += dt * 0.5 * force[j * N * 2 * N + k * N + l][0]/m;
                    v[j * N * 2 * N + k * N + l][1] += dt * 0.5 * force[j * N * 2 * N + k * N + l][1]/m;
                    v[j * N * 2 * N + k * N + l][2] += dt * 0.5 * force[j * N * 2 * N + k * N + l][2]/m;

                    position[j * N * 2 * N + k * N + l][0] += dt * v[j * N * 2 * N + k * N + l][0];
                    position[j * N * 2 * N + k * N + l][1] += dt * v[j * N * 2 * N + k * N + l][1];
                    position[j * N * 2 * N + k * N + l][2] += dt * v[j * N * 2 * N + k * N + l][2];

                }
            }
        }
        //printf("pos in verlet:%f\n", position[250][1]);

        get_forces_AL(force, position, cell_length, n_particles);
        //printf("force:%f\n", force[250][1]);

        for(int j = 0; j < 2*N; ++j){
            for(int k = 0; k < 2*N; ++k){
                for(int l = 0; l < N; ++l){

                    v[j * N * 2 * N + k * N + l][0] += dt * 0.5 * force[j * N * 2 * N + k * N + l][0]/m;
                    v[j * N * 2 * N + k * N + l][1] += dt * 0.5 * force[j * N * 2 * N + k * N + l][1]/m;
                    v[j * N * 2 * N + k * N + l][2] += dt * 0.5 * force[j * N * 2 * N + k * N + l][2]/m;

                }
            }
        }

        for(int j = 0; j < 2*N; ++j){
            for(int k = 0; k < 2*N; ++k){
                for(int l = 0; l < N; ++l){

                    position_evolution[i][j * N * 2 * N + k * N + l][0] = position[j * N * 2 * N + k * N + l][0];
                    position_evolution[i][j * N * 2 * N + k * N + l][1] = position[j * N * 2 * N + k * N + l][1];
                    position_evolution[i][j * N * 2 * N + k * N + l][2] = position[j * N * 2 * N + k * N + l][2];
                    
                }
            }
        }

        E_kin[i] = 0;
        for(int j = 0; j < 2*N; ++j){
            for(int k = 0; k < 2*N; ++k){
                for(int l = 0; l < N; ++l){

                    E_kin[i] += (m/2) * (v[j * N * 2 * N + k * N + l][0] * v[j * N * 2 * N + k * N + l][0] +
                                         v[j * N * 2 * N + k * N + l][1] * v[j * N * 2 * N + k * N + l][1] +
                                         v[j * N * 2 * N + k * N + l][2] * v[j * N * 2 * N + k * N + l][2]);
                    
                }
            }
        }
    }
}

void write_to_file(char *fname, double *times,
            double *pot_energy, double *kin_energy, 
            double *tot_energy, int n_points)
{
    FILE *fp = fopen(fname, "w");
    fprintf(fp, "time, pot energy, kin energy, tot energy\n");
    for(int i = 0; i < n_points; ++i){
	    fprintf(fp, "%f,%f,%f,%f\n", times[i],
        pot_energy[i], kin_energy[i], tot_energy[i]);
    }
    fclose(fp);
}
