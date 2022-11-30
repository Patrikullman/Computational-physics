#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

/* Defining the function to be integrated */
double func(double x, double y, double z){
    double PI = 3.14159265;
    
    double f = pow(PI, (double) -3/2) * ( x*x + x*x*y*y + x*x*y*y*z*z ) * exp( -(x*x + y*y + z*z) );
    return f;
    }

/* Defining the weight function */
double weight(double x, double y, double z){
    double PI = 3.14159265;
    
    double f = pow(PI, (double) -3/2) * exp( -(x*x + y*y + z*z) );
    return f;
}

int main()
{

    const gsl_rng_type * T;
    gsl_rng * r;

    gsl_rng_env_setup();

    T = gsl_rng_default;
    r = gsl_rng_alloc(T);

    int seed = 42;
    gsl_rng_set(r, seed);
    
   double delta = 2.3; // choose a delta so slightly less than 50% of the proposed steps get accepted
   double x_current = 0.; // start position x
   double y_current = 0.; // start position y
   double z_current = 0.; // start position z

   double inputs[10000][3]; // array to store inputs 
   int accept = 0; // acceptance count
   int tot_count = 0; // total loop count
   double sum = 0;
   double mean = 0;

   for(int i = 0; i < 10000; ++i){
        /* Write the current position to the input array */
        inputs[i][0] = x_current;
        inputs[i][1] = y_current;
        inputs[i][2] = z_current;

        /* Propose a new step for x, y and z */
        double x_prop = x_current + delta * (gsl_rng_uniform(r) - 0.5);
        double y_prop = y_current + delta * (gsl_rng_uniform(r) - 0.5);
        double z_prop = z_current + delta * (gsl_rng_uniform(r) - 0.5);

        /* Check the current and proposed likelihood for the steps and compare the */
        double current_lik = weight(x_current, y_current, z_current);
        //printf("current lik: %f\n", current_lik);
        double prop_lik = weight(x_prop, y_prop, z_prop);
        //printf("prop lik: %f\n", prop_lik);
        double ratio  = prop_lik / current_lik;

        /* Generate uniform number [0,1] */
        double event = gsl_rng_uniform(r);

        /* Accept new step if ratio >= event*/
        if (event <= ratio)
        {
            x_current = x_prop;
            y_current = y_prop;
            z_current = z_prop;
            accept += 1;
        }

        tot_count += 1;
        //printf("x: %f, y: %f, z: %f\n", x_current,y_current,z_current);
   }
   printf("acceptance count: %d\n", accept);
   printf("total count: %d\n", tot_count);
   printf("acceptance ratio: %f\n", (double)accept/tot_count);

    
    for(int i = 0; i < 10000; ++i){
        sum += func(inputs[i][0], inputs[i][1], inputs[i][2]);
        //double x = inputs[i][0];
        //double y = inputs[i][1];
        //double z = inputs[i][2];
        //double f = func(inputs[i][0], inputs[i][1], inputs[i][2]);
        //printf("x: %f, y: %f, z: %f sum: %f\n", x,y,z,f);
    }

    /* Calculate integral */
    double I = (double) sum / 10000;
    printf("Integral value: %f\n", I);



    return 0;
}
