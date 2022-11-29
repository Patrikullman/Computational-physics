#include <stdio.h>
#include <math.h>
#include <stdlib.h>

int
run(
    int argc,
    char *argv[]
   )
{
    
    double T = 400; //kör scriptet och ändra här, mellan 390-->900K
    double P[60] = {0.99,0.98,0.96,0.94,0.92,0.90,0.88,0.86,0.84,0.82,0.8,0.78,0.76,0.74,0.72,0.7,0.68,0.66,0.64,0.62,0.60,0.58,0.56,0.54,0.52,0.50,0.48,0.46,0.44,0.42,0.40,0.38,0.36,0.34,0.32,0.30,0.28,0.26,0.024,0.22,0.20,0.18,0.16,0.14,0.12,0.10,0.08,0.06};
    int N = 256;
    double Ecucu = -436e-3;
    double Ezizi = -113e-3;
    double Ecuzi = -294e-3;
    double E0 = 2*N*(Ecucu + Ezizi + 2*Ecuzi);
    double dE = Ecucu + Ezizi - 2*Ecuzi;
    double kB = 8.617333e-5;
    double* F = (double*)malloc(60 * sizeof(double*));
    //or (int i = 0; i < 11; i++)
     //   F[i] = (double*)malloc(1500 * sizeof(double));
    
    
    FILE *fp = fopen("f_values.csv", "w");
    for(int i = 1; i < 60; ++i){
            F[i] = E0 - 2*N*P[i]*P[i]*dE - 2*N*kB*T*log(2) + N*kB*T*((1+P[i])*log(1+P[i])+(1-P[i])*log(1-P[i]));
            fprintf(fp, "%f,%f,%f\n", P[i], T, F[i]);
        }
    
    fclose(fp);
free(F);

    return 0;
}	
