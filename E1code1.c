/******************************************************************************
 * E1code1
 ******************************************************************************
 * generate h(t) = a*cos(2*pi*f*t + phi)   ans then   + a2*cos(2*pi*f2*t + phi2)
 * and writes the result to a h_t.csv
 *
 * Compile me as:
 * clang E1code1.c -o <executable name> -lm
 */
/******************************************************************************
 * Includes
 *****************************************************************************/
#include <stdio.h>   //fopen, fprintf
#include <stdlib.h>  //malloc
#include <math.h>    //cos
#include <stdint.h>  //uint64_t

/******************************************************************************
 * Constants
 *****************************************************************************/
#define PI 3.14159


/******************************************************************************
 * Helper functions
 *****************************************************************************/
/*
 * constructs the signal
 * @signal - array to be filled with signal values
 * @t - time array filled with discrete time stamps
 * @len_t - the length of the time array
 * @a - amplitude of signal
 * @f - frequency of signal
 * @phi - phase of signal
*/
double *generate_signal(double *signal, double *t, uint64_t len_t, double a,
			double f, double phi,double a2, double f2, double phi2) 
{
    for(int i = 0; i < len_t; i++){
	signal[i] = a * cos(2 * PI * f * t[i] + phi) + a2 * cos(2 * PI * f2 * t[i] + phi2);
    }
    return signal;
}

/*
 * constructs time array
 * @array - array to be filled with time values
 * @start - start value
 * @len_t - number of times stamps in array
 * @dt - time step between two consecutive times
*/
void arange(double *array, double start, int len_t, double dt){
    for(int i = 0; i < len_t; i++){
	array[i] = start + i*dt;
    }
} 

/*
 * constructs time array
 * @fname - File name 
 * @time_array - array of time values
 * @signal - array with signal values
 * @n_points - number of points
*/
void write_to_file(char *fname, double *time_array,
		   double *signal, int n_points)
{
    FILE *fp = fopen(fname, "w");
    fprintf(fp, "time, signal\n");
    for(int i = 0; i < n_points; ++i){
	    fprintf(fp, "%f,%f\n", time_array[i], signal[i]);
    }
    fclose(fp);
}

int main()
{
    int N = 250; double dt = 0.083;    // we choose the sampling interval according to the nyquist critical frequency. 1/(2*6) = 0.083 see page 66 MD pdf canvas. 
    double a = 1; double f = 2; double phi = 0;
    double a2 =1; double f2 = 6; double phi2 = 0;
    double time_array[N];
    arange(time_array, 0, N, dt);
    
    double signal[N];
    generate_signal(signal, time_array, N, a, f, phi,a2, f2, phi2);
    write_to_file("signal.csv", time_array, signal, N);
    return 0;
}
