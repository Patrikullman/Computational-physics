#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>



int main()
{
    /* Open the MC.txt file and read it to file pointer */
    FILE *fp = fopen("MC.txt", "r");

    /* Create array to store values from the MC.txt file */
    double* numArray = (double*)malloc(1e6 * sizeof(double));
    
    /* Loop through the values written to the file pointer and put in array */
    for (int i = 0; i < 1e6; i++)
    {
        fscanf(fp, "%lf", &numArray[i]);
    }

    double sum = 0;
    double sum_squared = 0;
    double mean = 0;
    double vari = 0;

    for(int i = 0; i < 1e6; ++i){
        sum += numArray[i]; // sum up all values in the array
        sum_squared += ( numArray[i] * numArray[i] ); // sum up the squared of all values in array
    }
    mean = sum / 1e6; // mean value of values in array
    vari = ( sum_squared / 1e6 ) - (mean * mean); // variance of values in array

    double* phi = (double*)malloc(1e6 * sizeof(double)); // allocate memeory to store phi values for each k

    /* Defining lower and upper limit for reasonable candidate to be the index for the statistical inefficiency s*/
    double lower = exp(-2)-0.001; 
    double upper = exp(-2)+0.001;
    /* Loop through all k and check if a given k in phi could be a candidate for the statistical inefficiency */
    for(int k = 0; k < 1000; ++k){ // note that i think k should go to 1e6 but that takes forever...
        int t = k + 1; // set initial index to loop from in numerator
        double num = 0; // store value of numerator in phi
        double den = 0; // store value of denominator in phi
        for(int i = t; i < 1e6; ++i){ // sum up numerator in phi
            num += (numArray[i]-mean)*(numArray[i-k]-mean);  
        }
        for (int i = 0; i < 1e6; ++i){ // sum up denominator in phi
            den += ( (numArray[i]-mean)*(numArray[i]-mean) );
        }
        
        phi[k] = num / den; // calculate value of phi for index k
        /* Check if phi[k] close to exp(-2), in that case k could be the statistical inefficiency s
           If thats the case, print the value of phi and the corresponding k */
        if ( (lower < phi[k]) && (phi[k] < upper) ) 
        {
            printf("Based on autocorrelation, phi: %f, k: %d, potential s: %d\n", phi[k], k, k);
        }
        
    }
    

    


    int NB = 100; // Number of blocks
    int B = 1e6 / NB; // The block size

    /* Allocate memory for each block average */
    double* Fj = (double*)malloc(NB * sizeof(double));
    
    /* Calculate block average */
    for(int j = 0; j < NB; ++j){
        double block_sum = 0;
        for(int i = 0; i < B; ++i){
            block_sum += numArray[j*B + i];    
        }
        Fj[j] = block_sum / B;
    }

    double Fj_sum = 0;
    double Fj_sum_squared = 0;

    /* Sum up the block averages and the square of all the block averages */
    for(int i = 0; i < NB; ++i){
        Fj_sum += Fj[i];
        Fj_sum_squared += ( Fj[i] * Fj[i] );
    }

    /* Calculate mean and variance of block averages */
    double Fj_mean = Fj_sum / NB;
    double Fj_vari = (Fj_sum_squared / NB) - (Fj_mean * Fj_mean);

    /* Calculate the statistical inefficiency based on block averaging */
    double s_block = ( B * Fj_vari ) / vari;
    printf("Based on block averaging, s = %f\n", s_block);
    

    return 0;
}


