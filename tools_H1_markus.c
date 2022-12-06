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

void velocity_verlet_scaled(int n_timesteps, int n_particles, double v[][3], double position[][3],
double position_evolution[][n_particles][3], double dt, double m, double force[][3], 
double cell_length, double E_kin[n_timesteps], double temp_inst_verlet_scaled[n_timesteps])
{
    FILE *fp = fopen("inst_temp.csv", "w");
    FILE *fp2 = fopen("inst_pres.csv", "w");
    FILE *fp3 = fopen("a0_evolution.csv", "w");
    int N = 4;
    
    for(int i = 1;  i < n_timesteps + 1; ++i){
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

        if(i < 300){ // Limiting time step for equilibration at higher temp  before cooling down to actual desired temp
            int T_eq = 2000 + 273; // The higher temp to ensure melting
            double kB = 1.380649e-23;
            double alpha_T[n_timesteps];
            double tau_T = 100 * dt;
            double temp_inst2 = 0;

            
            temp_inst2 = (2*E_kin[i]*1.60218e-19) / (3 * kB * (4*N*N*N) );
            alpha_T[i] = 1 + (((2*dt) / tau_T) * ((T_eq - temp_inst2) / temp_inst2));
            temp_inst_verlet_scaled[i] = temp_inst2;

            /* Scale the velocity and turn off scaling at time step 600 at which equilibrium is considered to be reached*/
            if(i < 600){
                for(int j = 0; j < 2*N; ++j){
                    for(int k = 0; k < 2*N; ++k){
                        for(int l = 0; l < N; ++l){

                            v[j * N * 2 * N + k * N + l][0] = sqrt(alpha_T[i]) * v[j * N * 2 * N + k * N + l][0];
                            v[j * N * 2 * N + k * N + l][1] = sqrt(alpha_T[i]) * v[j * N * 2 * N + k * N + l][1];
                            v[j * N * 2 * N + k * N + l][2] = sqrt(alpha_T[i]) * v[j * N * 2 * N + k * N + l][2];

                        }
                    }
                }
            }
            
            fprintf(fp, "%f, %f\n", temp_inst2, i*dt);
            
            //FILE *fp2 = fopen("inst_pres.csv", "w");
            //FILE *fp3 = fopen("a0_evolution.csv", "w");
            int P_eq = 1/1.60218e+6; //100000 / 1.60218e+11; // Pressure in asu at which we wish to equilibrate 
            double alpha_P[n_timesteps];
            double tau_P = 100*dt;
            double kappa_T = 2.32; //0.04e+9 / 1.60218e+11; //0.04GPa^-1
            double P_inst = 0;

            P_inst = ( E_kin[i] + 3 * get_virial_AL(position, cell_length, n_particles) ) / 
            (3 * cell_length * cell_length * cell_length);
            //printf("p inst: %f\n", P_inst);
            alpha_P[i] = 1 - ( ( kappa_T * dt / tau_P ) * (P_eq - P_inst) ); 

            /* Scale the position while time step less than 600 */
            if(i < 600){
                for(int j = 0; j < 2*N; ++j){
                    for(int k = 0; k < 2*N; ++k){
                        for(int l = 0; l < N; ++l){

                            position[j * N * 2 * N + k * N + l][0] = cbrt(alpha_P[i]) * position[j * N * 2 * N + k * N + l][0];
                            position[j * N * 2 * N + k * N + l][1] = cbrt(alpha_P[i]) * position[j * N * 2 * N + k * N + l][1];
                            position[j * N * 2 * N + k * N + l][2] = cbrt(alpha_P[i]) * position[j * N * 2 * N + k * N + l][2];

                        }
                    }
                }
            }
            fprintf(fp2, "%f, %f\n", P_inst, i*dt);
            fprintf(fp3, "%f, %f\n", cell_length, i*dt);
            cell_length = cbrt(alpha_P[i]) * cell_length;
            //printf("alpha_P: %f\n", alpha_P[i]);
            //printf("a0: %f\n", cell_length / N);
        }
        else{
            int T_eq = 700 + 273; // The actual temp at which we wish to equilibrate
            double kB = 1.380649e-23;
            double alpha_T[n_timesteps];
            double tau_T = 100 * dt;
            double temp_inst2 = 0;

            
            temp_inst2 = (2*E_kin[i]*1.60218e-19) / (3 * kB * (4*N*N*N) );
            alpha_T[i] = 1 + (((2*dt) / tau_T) * ((T_eq - temp_inst2) / temp_inst2));
            temp_inst_verlet_scaled[i] = temp_inst2;

            /* Scale the velocity and turn off scaling at time step 600 at which equilibrium is considered to be reached*/
            if(i < 600){
                for(int j = 0; j < 2*N; ++j){
                    for(int k = 0; k < 2*N; ++k){
                        for(int l = 0; l < N; ++l){

                            v[j * N * 2 * N + k * N + l][0] = sqrt(alpha_T[i]) * v[j * N * 2 * N + k * N + l][0];
                            v[j * N * 2 * N + k * N + l][1] = sqrt(alpha_T[i]) * v[j * N * 2 * N + k * N + l][1];
                            v[j * N * 2 * N + k * N + l][2] = sqrt(alpha_T[i]) * v[j * N * 2 * N + k * N + l][2];

                        }
                    }
                }
            }
            
            fprintf(fp, "%f, %f\n", temp_inst2, i*dt);
            
            //FILE *fp2 = fopen("inst_pres.csv", "w");
            //FILE *fp3 = fopen("a0_evolution.csv", "w");
            int P_eq = 1/1.60218e+6; //100000 / 1.60218e+11; // Pressure in asu at which we wish to equilibrate 
            double alpha_P[n_timesteps];
            double tau_P = 100*dt;
            double kappa_T = 2.32; //0.04e+9 / 1.60218e+11; //0.04GPa^-1
            double P_inst = 0;

            P_inst = ( E_kin[i] + 3 * get_virial_AL(position, cell_length, n_particles) ) / 
            (3 * cell_length * cell_length * cell_length);
            //printf("p inst: %f\n", P_inst);
            alpha_P[i] = 1 - ( ( kappa_T * dt / tau_P ) * (P_eq - P_inst) ); 

            /* Scale the position while time step less than 600 */
            if(i < 600){
                for(int j = 0; j < 2*N; ++j){
                    for(int k = 0; k < 2*N; ++k){
                        for(int l = 0; l < N; ++l){

                            position[j * N * 2 * N + k * N + l][0] = cbrt(alpha_P[i]) * position[j * N * 2 * N + k * N + l][0];
                            position[j * N * 2 * N + k * N + l][1] = cbrt(alpha_P[i]) * position[j * N * 2 * N + k * N + l][1];
                            position[j * N * 2 * N + k * N + l][2] = cbrt(alpha_P[i]) * position[j * N * 2 * N + k * N + l][2];

                        }
                    }
                }
            }
            fprintf(fp2, "%f, %f\n", P_inst, i*dt);
            fprintf(fp3, "%f, %f\n", cell_length, i*dt);
            cell_length = cbrt(alpha_P[i]) * cell_length;
              
            //printf("alpha_P: %f\n", alpha_P[i]);
            //printf("a0: %f\n", cell_length / N);
        }

         

    }
    fclose(fp);
    fclose(fp3);
    fclose(fp2);

    /*
    fclose(fp);

    FILE *fp2 = fopen("inst_pres.csv", "w");
    FILE *fp3 = fopen("a0_evolution.csv", "w");
    for(int i = 0; i < n_timesteps; ++i){ 
        int P_eq = 100000 / 1.60218e+11; // Pressure in asu at which we wish to equilibrate 
        double alpha_P[n_timesteps];
        double tau_P = 100*dt;
        double kappa_T = 0.04e+9 / 1.60218e+11; //0.04GPa^-1
        double P_inst = 0;

        P_inst = ( E_kin[i] + 3 * get_virial_AL(position, cell_length, n_particles) ) / 
        (3 * cell_length * cell_length * cell_length);
        alpha_P[i] = 1 - ( ( kappa_T * dt / tau_P ) * (P_eq - P_inst) ); 

        for(int j = 0; j < 2*N; ++j){
            for(int k = 0; k < 2*N; ++k){
                for(int l = 0; l < N; ++l){

                    position[j * N * 2 * N + k * N + l][0] = cbrt(alpha_P[i]) * position[j * N * 2 * N + k * N + l][0];
                    position[j * N * 2 * N + k * N + l][1] = cbrt(alpha_P[i]) * position[j * N * 2 * N + k * N + l][1];
                    position[j * N * 2 * N + k * N + l][2] = cbrt(alpha_P[i]) * position[j * N * 2 * N + k * N + l][2];

                }
            }
        }
        fprintf(fp2, "%f, %f\n", P_inst, i*dt);
        fprintf(fp3, "%f, %f\n", cell_length, i*dt);
        cell_length = cbrt(alpha_P[i]) * cell_length;
        printf("a0: %f\n", cell_length / N);
        
        
    }
    
    fclose(fp3);
    fclose(fp2);
    */
}


double distance(double point1[3], double point2[3]){
    double diff_x = point1[0]-point2[0];
    double diff_y = point1[1]-point2[1];
    double diff_z = point1[2]-point2[2];

    double dist = sqrt(diff_x*diff_x + diff_y*diff_y + diff_z*diff_z);
    return dist;
}