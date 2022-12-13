#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include "tools.h"

int
run(
    int argc,
    char *argv[]
   )
{
    int n_T = 100;
    int n_P = 1000;
    
    double* T = (double*)malloc(n_T * sizeof(double));
    double* P = (double*)malloc(n_P * sizeof(double));

    int T_start = 300;
    for(int i = 0; i < n_T; ++i){
        T[i] = T_start + 10*i;
    }

    int P_start = 0.01;
    for(int i = 0; i < n_P; ++i){
        P[i] = P_start + 0.001*i;
    }

    int N = 256;
    double Ecucu = -436e-3; // Binding energy Cu-Cu
    double Ezizi = -113e-3; // Binding energy Zi-Zi
    double Ecuzi = -294e-3; // Binding energy Cu-Zi
    double E0 = 2*N*(Ecucu + Ezizi + 2*Ecuzi);
    double dE = Ecucu + Ezizi - 2*Ecuzi;
    double kB = 8.617333e-5; // Boltzmann const in eV
    double* F = (double*)malloc(n_P * sizeof(double)); // Allocate array for free energy values
    double* min_P = (double*)malloc(n_T * sizeof(double)); // Allocate array for P values corresponding to minimum free energy at specific temp
    double* E_mf = (double*)malloc(n_T * sizeof(double)); // Allocate array for the mean field energy
    
    /* Calculate the order param, free energy, and mean field energy at different temps */
    FILE *fp = fopen("f_values.csv", "w");
    FILE *fp2 = fopen("PT_values.csv", "w");
    for(int j = 0; j < n_T; ++j){
        for(int i = 1; i < n_P; ++i){
                F[i] = E0 - 2*N*P[i]*P[i]*dE - 2*N*kB*T[j]*log(2) + N*kB*T[j]*((1+P[i])*log(1+P[i])+(1-P[i])*log(1-P[i]));
                fprintf(fp, "%f,%f,%f\n", P[i], T[j], F[i]);
        }
        
        int k = 0;
        while (F[k] > F[k+1])
        {
            min_P[j] = P[k+1];
            k += 1;
        }   
        
        E_mf[j] = E0 - 2*N*min_P[j]*min_P[j]*dE;
        fprintf(fp2, "%f,%f,%f\n", T[j], min_P[j], E_mf[j]);
    }


    
    fclose(fp2);
    fclose(fp);
    free(F);
    free(min_P);

    /* Calculate heat capacity at different temps */
    
    double* dx = (double*)malloc((n_T-1) * sizeof(double)); // Allocate array for delta_x step
    double* dy = (double*)malloc((n_T-1) * sizeof(double)); // Allocate array for delta_y step
    double* C = (double*)malloc((n_T-1) * sizeof(double)); // Allocate array for heat capacity
    double* T_deri = (double*)malloc((n_T-1) * sizeof(double)); // Allocate array for the temps between which the slope is calculated

    /* Set values to T_deri */
    for(int i = 0; i < n_T-1; ++i){
        T_deri[i] = T[i];
    }

    /* Calculate the slopes/derivatives and hence also the heat capacity */
    FILE *fp3 = fopen("heat_capa.csv", "w");
    for(int i = 0; i < n_T-1; ++i){
        dx[i] = T[i] - T[i+1];
        dy[i] = E_mf[i] - E_mf[i+1];
        C[i] = dy[i] / dx[i];
        fprintf(fp3, "%f,%f\n", T_deri[i], C[i]);
    }
    fclose(fp3);

    int* temperatures = (int*)malloc(100 * sizeof(int));
    int no_of_iterations = 6000; // Number of Monte Carlo steps


    FILE *fp4 = fopen("temp_energy.csv", "w");
    for(int t = 1; t < 101; ++t){

        temperatures[t] = 10*t; // put temperature t in the temperature array
        double temp = 10*t; // Set temperature in kelvin

        double lattice_a[10][10][10]; // The initial Cu atoms
        double lattice_b[10][10][10]; // The initial Zi atoms
        double lattice_a_prop[10][10][10];
        double lattice_b_prop[10][10][10];


        // Set initial values 
        for(int i = 0; i < 10; ++i){
            for(int j = 0; j < 10; ++j){
                for(int k = 0; k < 10; ++k){
                    lattice_a[i][j][k] = 0; // Lattice a, originally containing Cu atoms
                    lattice_b[i][j][k] = 1; // Lattice b, originally containing Zi atoms
                    lattice_a_prop[i][j][k] = 0; // A proposed a lattice when atoms has been swaped
                    lattice_b_prop[i][j][k] = 1; // A proposed b lattice when atoms has been swaped
                }
            }
        }

        int Naa = 0; // Number of Cu-Cu bindings 
        int Nbb = 0; // Number of Zi-Zi bindings
        int Nab = 0; // Number of Cu-Zi bindings

        /* Loop through the two lattice and calculate an initial energy, with boundary conditions taken into a count */
        for(int i = 0; i < 10; ++i){
            for(int j = 0; j < 10; ++j){
                for(int k = 0; k < 10; ++k){
                    if(lattice_a[i][j][k] == lattice_b[i][j][k]){
                        if(lattice_a[i][j][k] == 0){
                            Naa += 1;
                        }
                        else{
                            Nbb += 1;
                        }
                    }
                    else{
                        Nab += 1;
                    }

                    if(i != 9){
                        if(lattice_a[i][j][k] == lattice_b[i+1][j][k]){
                            if(lattice_a[i][j][k] == 0){
                                Naa += 1;
                            }
                            else{
                                Nbb += 1;
                            }
                        }
                        else{
                            Nab += 1;
                        }
                    }
                    else{
                        if(lattice_a[i][j][k] == lattice_b[0][j][k]){
                            if(lattice_a[i][j][k] == 0){
                                Naa += 1;
                            }
                            else{
                                Nbb += 1;
                            }
                        }
                        else{
                            Nab += 1;
                        }
                    }

                    if(j != 9){
                        if(lattice_a[i][j][k] == lattice_b[i][j+1][k]){
                            if(lattice_a[i][j][k] == 0){
                                Naa += 1;
                            }
                            else{
                                Nbb += 1;
                            }
                        }
                        else{
                            Nab += 1;
                        }
                    }
                    else{
                        if(lattice_a[i][j][k] == lattice_b[i][0][k]){
                            if(lattice_a[i][j][k] == 0){
                                Naa += 1;
                            }
                            else{
                                Nbb += 1;
                            }
                        }
                        else{
                            Nab += 1;
                        }
                    }

                    if(k != 9){
                        if(lattice_a[i][j][k] == lattice_b[i][j][k+1]){
                            if(lattice_a[i][j][k] == 0){
                                Naa += 1;
                            }
                            else{
                                Nbb += 1;
                            }
                        }
                        else{
                            Nab += 1;
                        }
                    }
                    else{
                        if(lattice_a[i][j][k] == lattice_b[i][j][0]){
                            if(lattice_a[i][j][k] == 0){
                                Naa += 1;
                            }
                            else{
                                Nbb += 1;
                            }
                        }
                        else{
                            Nab += 1;
                        }
                    }

                    if(i != 9 && j != 9){
                        if(lattice_a[i][j][k] == lattice_b[i+1][j+1][k]){
                            if(lattice_a[i][j][k] == 0){
                                Naa += 1;
                            }
                            else{
                                Nbb += 1;
                            }
                        }
                        else{
                            Nab += 1;
                        }
                    }
                    else{
                        if(lattice_a[i][j][k] == lattice_b[0][0][k]){
                            if(lattice_a[i][j][k] == 0){
                                Naa += 1;
                            }
                            else{
                                Nbb += 1;
                            }
                        }
                        else{
                            Nab += 1;
                        }
                    }

                    if(j != 9 && k != 9){
                        if(lattice_a[i][j][k] == lattice_b[i][j+1][k+1]){
                            if(lattice_a[i][j][k] == 0){
                                Naa += 1;
                            }
                            else{
                                Nbb += 1;
                            }
                        }
                        else{
                            Nab += 1;
                        }
                    }
                    else{
                        if(lattice_a[i][j][k] == lattice_b[i][0][0]){
                            if(lattice_a[i][j][k] == 0){
                                Naa += 1;
                            }
                            else{
                                Nbb += 1;
                            }
                        }
                        else{
                            Nab += 1;
                        }
                    }

                    if(i != 9 && k != 9){
                        if(lattice_a[i][j][k] == lattice_b[i+1][j][k+1]){
                            if(lattice_a[i][j][k] == 0){
                                Naa += 1;
                            }
                            else{
                                Nbb += 1;
                            }
                        }
                        else{
                            Nab += 1;
                        }
                    }
                    else{
                        if(lattice_a[i][j][k] == lattice_b[0][j][0]){
                            if(lattice_a[i][j][k] == 0){
                                Naa += 1;
                            }
                            else{
                                Nbb += 1;
                            }
                        }
                        else{
                            Nab += 1;
                        }
                    }

                    if(i != 9 && j != 9 && k != 9){
                        if(lattice_a[i][j][k] == lattice_b[i+1][j+1][k+1]){
                            if(lattice_a[i][j][k] == 0){
                                Naa += 1;
                            }
                            else{
                                Nbb += 1;
                            }
                        }
                        else{
                            Nab += 1;
                        }
                    }
                    else{
                        if(lattice_a[i][j][k] == lattice_b[0][0][0]){
                            if(lattice_a[i][j][k] == 0){
                                Naa += 1;
                            }
                            else{
                                Nbb += 1;
                            }
                        }
                        else{
                            Nab += 1;
                        }
                    }
                }
            }
        }

        //printf("Naa: %d, Nbb: %d, Nab: %d\n", Naa, Nbb, Nab);

        double E_start = Naa * Ecucu + Nbb * Ezizi + Nab * Ecuzi; // Initial energy
        double E_current = E_start; // Set the current energy to the initial
        int accept_count = 0;

        const gsl_rng_type * TT;
        gsl_rng * r;

        gsl_rng_env_setup();

        TT = gsl_rng_default;
        r = gsl_rng_alloc(TT);

        int seed = 42;
        gsl_rng_set(r, seed);


        /* Monte Carlo with no_of_iterations number of steps */
        for(int l = 0; l < no_of_iterations; ++l){
        

            // Get random index of atom in lattice_a to swap with neighbour from lattice_b
            int x_index = 9 * gsl_rng_uniform(r);
            int y_index = 9 * gsl_rng_uniform(r);
            int z_index = 9 * gsl_rng_uniform(r);

            //printf("x index: %d, index y: %d, index z: %d\n", x_index, y_index, z_index);

            // Get which neighbour in lattice_b to swap with the lattice_a atom (based on the position of the lattice_a atom)
            int which_NN_in_b_x = x_index;
            int which_NN_in_b_y = y_index;
            int which_NN_in_b_z = z_index;
            if(gsl_rng_uniform(r) > 0.5){
                which_NN_in_b_x = x_index + 1;
            }
            if(gsl_rng_uniform(r) > 0.5){
                which_NN_in_b_y = z_index + 1;
            }
            if(gsl_rng_uniform(r) > 0.5){
                which_NN_in_b_z = z_index + 1;
            }

            

            // Value from lattice_b to put in lattice_a
            int to_a = lattice_b[which_NN_in_b_x][which_NN_in_b_y][which_NN_in_b_z];
            // Replace position in the proposed lattice_b with atom from lattice_a
            lattice_b_prop[which_NN_in_b_x][which_NN_in_b_y][which_NN_in_b_z] = lattice_a[x_index][y_index][z_index];
            // Insert value from lattice_b in the proposed lattice_a
            lattice_a_prop[x_index][y_index][z_index] = to_a;

            int Naa_MC = 0; // Number of Cu-Cu bindings 
            int Nbb_MC = 0; // Number of Zi-Zi bindings
            int Nab_MC = 0; // Number of Cu-Zi bindings

            // Check number of nearest neighbour of each atom type
            for(int i = 0; i < 10; ++i){
                for(int j = 0; j < 10; ++j){
                    for(int k = 0; k < 10; ++k){
                        if(lattice_a_prop[i][j][k] == lattice_b_prop[i][j][k]){
                            if(lattice_a_prop[i][j][k] == 0){
                                Naa_MC += 1;
                            }
                            else{
                                Nbb_MC += 1;
                            }
                        }
                        else{
                            Nab_MC += 1;
                        }

                        if(i != 9){
                            if(lattice_a_prop[i][j][k] == lattice_b_prop[i+1][j][k]){
                                if(lattice_a_prop[i][j][k] == 0){
                                    Naa_MC += 1;
                                }
                                else{
                                    Nbb_MC += 1;
                                }
                            }
                            else{
                                Nab_MC += 1;
                            }
                        }
                        else{
                            if(lattice_a_prop[i][j][k] == lattice_b_prop[0][j][k]){
                                if(lattice_a_prop[i][j][k] == 0){
                                    Naa_MC += 1;
                                }
                                else{
                                    Nbb_MC += 1;
                                }
                            }
                            else{
                                Nab_MC += 1;
                            }
                        }

                        if(j != 9){
                            if(lattice_a_prop[i][j][k] == lattice_b_prop[i][j+1][k]){
                                if(lattice_a_prop[i][j][k] == 0){
                                    Naa_MC += 1;
                                }
                                else{
                                    Nbb_MC += 1;
                                }
                            }
                            else{
                                Nab_MC += 1;
                            }
                        }
                        else{
                            if(lattice_a_prop[i][j][k] == lattice_b_prop[i][0][k]){
                                if(lattice_a_prop[i][j][k] == 0){
                                    Naa_MC += 1;
                                }
                                else{
                                    Nbb_MC += 1;
                                }
                            }
                            else{
                                Nab_MC += 1;
                            }
                        }

                        if(k != 9){
                            if(lattice_a_prop[i][j][k] == lattice_b_prop[i][j][k+1]){
                                if(lattice_a_prop[i][j][k] == 0){
                                    Naa_MC += 1;
                                }
                                else{
                                    Nbb_MC += 1;
                                }
                            }
                            else{
                                Nab_MC += 1;
                            }
                        }
                        else{
                            if(lattice_a_prop[i][j][k] == lattice_b_prop[i][j][0]){
                                if(lattice_a_prop[i][j][k] == 0){
                                    Naa_MC += 1;
                                }
                                else{
                                    Nbb_MC += 1;
                                }
                            }
                            else{
                                Nab_MC += 1;
                            }
                        }

                        if(i != 9 && j != 9){
                            if(lattice_a_prop[i][j][k] == lattice_b_prop[i+1][j+1][k]){
                                if(lattice_a_prop[i][j][k] == 0){
                                    Naa_MC += 1;
                                }
                                else{
                                    Nbb_MC += 1;
                                }
                            }
                            else{
                                Nab_MC += 1;
                            }
                        }
                        else{
                            if(lattice_a_prop[i][j][k] == lattice_b_prop[0][0][k]){
                                if(lattice_a_prop[i][j][k] == 0){
                                    Naa_MC += 1;
                                }
                                else{
                                    Nbb_MC += 1;
                                }
                            }
                            else{
                                Nab_MC += 1;
                            }
                        }

                        if(j != 9 && k != 9){
                            if(lattice_a_prop[i][j][k] == lattice_b_prop[i][j+1][k+1]){
                                if(lattice_a_prop[i][j][k] == 0){
                                    Naa_MC += 1;
                                }
                                else{
                                    Nbb_MC += 1;
                                }
                            }
                            else{
                                Nab_MC += 1;
                            }
                        }
                        else{
                            if(lattice_a_prop[i][j][k] == lattice_b_prop[i][0][0]){
                                if(lattice_a_prop[i][j][k] == 0){
                                    Naa_MC += 1;
                                }
                                else{
                                    Nbb_MC += 1;
                                }
                            }
                            else{
                                Nab_MC += 1;
                            }
                        }

                        if(i != 9 && k != 9){
                            if(lattice_a_prop[i][j][k] == lattice_b_prop[i+1][j][k+1]){
                                if(lattice_a_prop[i][j][k] == 0){
                                    Naa_MC += 1;
                                }
                                else{
                                    Nbb_MC += 1;
                                }
                            }
                            else{
                                Nab_MC += 1;
                            }
                        }
                        else{
                            if(lattice_a_prop[i][j][k] == lattice_b_prop[0][j][0]){
                                if(lattice_a_prop[i][j][k] == 0){
                                    Naa_MC += 1;
                                }
                                else{
                                    Nbb_MC += 1;
                                }
                            }
                            else{
                                Nab_MC += 1;
                            }
                        }

                        if(i != 9 && j != 9 && k != 9){
                            if(lattice_a_prop[i][j][k] == lattice_b_prop[i+1][j+1][k+1]){
                                if(lattice_a_prop[i][j][k] == 0){
                                    Naa_MC += 1;
                                }
                                else{
                                    Nbb_MC += 1;
                                }
                            }
                            else{
                                Nab_MC += 1;
                            }
                        }
                        else{
                            if(lattice_a_prop[i][j][k] == lattice_b_prop[0][0][0]){
                                if(lattice_a_prop[i][j][k] == 0){
                                    Naa_MC += 1;
                                }
                                else{
                                    Nbb_MC += 1;
                                }
                            }
                            else{
                                Nab_MC += 1;
                            }
                        }
                    }
                }
            }

            double E_prop = Naa_MC * Ecucu + Nbb_MC * Ezizi + Nab_MC * Ecuzi; // Energy with the proposed lattices
            //printf("Naa MC: %d, Nbb MC: %d, Nab MC: %d\n", Naa_MC, Nbb_MC, Nab_MC);
            printf("E prop: %f, E current: %f\n", E_prop, E_current);

            double delta_E = E_prop - E_current; // Energy difference
            
            double condition = exp(-delta_E / (kB*temp)); // Condition to except new energy and lattice configuration

        

            /* If condition > rand uniform number, update current energy and set the a and b lattices equal to the proposed ones*/
            
            if(E_prop < E_current){
                accept_count += 1;
                E_current = E_prop;
                
                for(int i = 0; i < 10; ++i){
                    for(int j = 0; j < 10; ++j){
                        for(int k = 0; k < 10; ++k){
                            lattice_a[i][j][k] = lattice_a_prop[i][j][k];
                            lattice_b[i][j][k] = lattice_b_prop[i][j][k];
                        }
                    }
                }
            }
            else{
                if(condition >= gsl_rng_uniform(r)){
                    accept_count += 1;

                    for(int i = 0; i < 10; ++i){
                        for(int j = 0; j < 10; ++j){
                            for(int k = 0; k < 10; ++k){
                                lattice_a[i][j][k] = lattice_a_prop[i][j][k];
                                lattice_b[i][j][k] = lattice_b_prop[i][j][k];
                            }
                        }
                    }
                }
            }
            
            printf("Accept count: %d\n", accept_count);
        }
        
        fprintf(fp4, "%f, %f\n", temp, E_current);
    

    }

    return 0;

}	