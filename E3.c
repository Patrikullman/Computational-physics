#include <stdio.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <time.h>
#include <math.h>

int
main()
{
/*
    const gsl_rng_type * T;
    gsl_rng * r;

    int i, n = 10000; //perform the integral for N= 10^1, 10^2, 10^3,10^4.
    gsl_rng_env_setup();

    T = gsl_rng_default;
    r = gsl_rng_alloc(T);

    int seed = 42;
    gsl_rng_set(r, seed); 
	double f_x;
	double f_x_squared;
	double* x;
	double variance_squared;
	double variance;
	x = (double*)malloc(10000 * sizeof(double));
    for (i = 0; i < n; i++){  
	double u = gsl_rng_uniform(r); //uniform numbers are stored in variable u
	x[i] = u; //storing random numbers in x array.
	
	f_x += x[i]*(1-x[i]); //using monte carlo to calculate the 1D integral, i.e summing up the integral with random numbers as input.
	f_x_squared += (x[i]*(1-x[i]))*(x[i]*(1-x[i])); //each term squared, for variance calc.
	printf("%f\n",x[i]);
	
	}

FILE *fp1 = fopen("uniform.csv", "w");
    for(int i = 0; i < n; ++i){
        fprintf(fp1,"%f\n", x[i]);
    }
    fclose(fp1);


double integral = f_x/n; //taking average, to get the value of the integral.
printf("monte carlo integral: %f\n",integral);
 
variance_squared = f_x_squared/n - ((f_x/n)*(f_x/n)); //calculating standard dev for each iteration
variance = sqrt(variance_squared); //variance
printf("Variance: %f\n",variance);

double I_plus = integral + variance/sqrt(n);
double I_minus = integral - variance/sqrt(n);
printf("I = I_n + var / sqrt(n):%f\n ", I_plus);
printf("I = I_n - var / sqrt(n):%f\n ", I_minus);


//for 10^1 : integral = 0.182055
//for 10^2 : integral = 0.160331
//for 10^3 : integral = 0.164759
//for 10^4 : integral = 0.166982
//true value according to symbolab is 0.16666....

*/

//                    Här börjar 2.
//------------------------------------------------------------------------------------------------------------------------------------------------------------
//We now want to do the monte carlo integration but with importance sampling and transformation method so we reduce the variance. 

//f = sin(PI*x); //p(x)	
//primitive_f = -(1/PI)*cos(PI*x) //primitive of sin(pi*x) this is our F(x).
//inverse_F = acos(-pi*x)/pi // F^-1
//#define PI 3.141592
double PI = 3.14159262;
    const gsl_rng_type * T;
    gsl_rng * r;

    int i, n = 1000; //perform the integral for N= 10^1, 10^2, 10^3,10^4.
    gsl_rng_env_setup();

    T = gsl_rng_default;
    r = gsl_rng_alloc(T);

    // typically the seed is set to the current time.
    //time_t seed = time(NULL);
    int seed = 4112;
    gsl_rng_set(r, seed); 
	double* x;
	double* eta;
	double g;
	x = (double*)malloc(10000 * sizeof(double)); //dynamically allocating an array, to store all u values.
	eta = (double*)malloc(10000 * sizeof(double));
    for (i = 0; i < n; i++){
        double u = gsl_rng_uniform(r); //uniform numbers are stored in variable u  
	x[i] = u; //our uniform random numbers
	eta [i] = acos(1-2*x[i])/ PI; //acos(PI * x[i]/PI) ; //putting our random numbers into the inverse function.
	printf("%f\n",eta[i]);
	g += (eta[i]*(1-eta[i])) / (PI*sin(eta[i])); //sin(PI*eta[i]);  // g(x) = f(x) / p(x)
	printf("%f\n",g);

	}

FILE *fp = fopen("eta.csv", "w");
    for(int i = 0; i < n; ++i){
        fprintf(fp,"%f\n", eta[i]);
    }
    fclose(fp);


double importance_sampling_integral = g/n;
printf("monte carlo integral importance sampling: %f\n",importance_sampling_integral);



    return 0;
}



