#include <stdio.h>
#include <math.h>
#include <stdlib.h>


int
run(
    int argc,
    char *argv[]
   )
{
    
    double P[11] = {0.9999, 0.9, 0.8, 0.7, 0.6, 0.5, 0.4, 0.3, 0.2, 0.1, 0.0};
    int N = 256;
    double Ecucu = -436e-3;
    double Ezizi = -113e-3;
    double Ecuzi = -294e-3;
    double E0 = 2*N*(Ecucu + Ezizi + 2*Ecuzi);
    double dE = Ecucu + Ezizi - 2*Ecuzi;
    double kB = 8.617333e-5;
    //double F_energy[11][1500];
    double** F = (double**)malloc(11 * sizeof(double*));
    for (int i = 0; i < 11; i++)
        F[i] = (double*)malloc(1500 * sizeof(double));
    
    FILE *fp = fopen("f_values.csv", "w");
    for(int i = 1; i < 2; ++i){
        for(int j = 0; j < 1500; ++j){
            int T = j;
            F[i][j] = E0 - 2*N*P[i]*P[i]*dE - 2*N*kB*T*log(2) + N*kB*T*((1+P[i])*log(1+P[i])+(1-P[i])*log(1-P[i]));
            fprintf(fp, "%f,%d,%f\n", P[i], T, F[i][j]);
        }
    }
    fclose(fp);

    return 0;
}
