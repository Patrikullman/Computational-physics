#include <stdio.h>
#include <math.h>
#include <stdlib.h>

int
run(
    int argc,
    char *argv[]
   )
{
    int n_T = 81;
    int n_P = 100;
    
    int* T = (int*)malloc(n_T * sizeof(int));
    double* P = (double*)malloc(n_P * sizeof(double));

    int T_start = 300;
    for(int i = 0; i < n_T; ++i){
        T[i] = T_start + 10*i;
    }

    int P_start = 0.01;
    for(int i = 0; i < n_P; ++i){
        P[i] = P_start + 0.01*i;
    }

    //double T[26] = {390, 420, 450, 500, 530, 560, 580, 600, 630, 660, 680, 700, 730, 760, 780, 800, 830, 860, 880, 900, 930, 960, 980, 1000, 1030, 1100}; //kör scriptet och ändra här, mellan 390-->900K
    //double P[60] = {0.99,0.98,0.96,0.94,0.92,0.90,0.88,0.86,0.84,0.82,0.8,0.78,0.76,0.74,0.72,0.7,0.68,0.66,0.64,0.62,0.60,0.58,0.56,0.54,0.52,0.50,0.48,0.46,0.44,0.42,0.40,0.38,0.36,0.34,0.32,0.30,0.28,0.26,0.024,0.22,0.20,0.18,0.16,0.14,0.12,0.10,0.08,0.06};
    int N = 256;
    double Ecucu = -436e-3;
    double Ezizi = -113e-3;
    double Ecuzi = -294e-3;
    double E0 = 2*N*(Ecucu + Ezizi + 2*Ecuzi);
    double dE = Ecucu + Ezizi - 2*Ecuzi;
    double kB = 8.617333e-5;
    double* F = (double*)malloc(n_P * sizeof(double));
    double* min_P = (double*)malloc(n_T * sizeof(double));
    double* E_mf = (double*)malloc(n_T * sizeof(double)); // The mean field energy
    
    
    FILE *fp = fopen("f_values.csv", "w");
    FILE *fp2 = fopen("PT_values.csv", "w");
    for(int j = 0; j < n_T; ++j){
        for(int i = 1; i < n_P; ++i){
                F[i] = E0 - 2*N*P[i]*P[i]*dE - 2*N*kB*T[j]*log(2) + N*kB*T[j]*((1+P[i])*log(1+P[i])+(1-P[i])*log(1-P[i]));
                fprintf(fp, "%f,%d,%f\n", P[i], T[j], F[i]);
        }
        
        int k = 0;
        while (F[k] > F[k+1])
        {
            min_P[j] = P[k+1];
            k += 1;
        }   
        
        E_mf[j] = E0 - 2*N*min_P[j]*min_P[j]*dE;
        fprintf(fp2, "%d,%f,%f\n", T[j], min_P[j], E_mf[j]);
    }


    
    fclose(fp2);
    fclose(fp);
    free(F);
    free(min_P);

    int* dx = (int*)malloc((n_T-1) * sizeof(int));
    double* dy = (double*)malloc((n_T-1) * sizeof(double));
    double* C = (double*)malloc((n_T-1) * sizeof(double));
    int* T_deri = (int*)malloc((n_T-1) * sizeof(int));

    for(int i = 0; i < n_T-1; ++i){
        T_deri[i] = T[i];
    }

    FILE *fp3 = fopen("heat_capa.csv", "w");
    for(int i = 0; i < n_T-1; ++i){
        dx[i] = T[i] - T[i+1];
        dy[i] = E_mf[i] - E_mf[i+1];
        C[i] = dy[i] / dx[i];
        fprintf(fp3, "%d,%f\n", T_deri[i], C[i]);
    }
    fclose(fp3);

    return 0;
}	