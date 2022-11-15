/******************************************************************************
 * E1code4
 ******************************************************************************
 * Routine that runs the velocity verlet algorithm
 * Use as template to construct your program!
 */

/*
 * Calculate the acceleration
 * @a - vector that is filled with acceleration
 * @u - vector with the current positions
 * @m - vector with masses
 * @kappa - Spring constant
 * @size_of_u - the size of the position, acceleration and mass array
 */
 #include <stdio.h>
 #include <stdlib.h>
void calc_acc(double *a, double *u, double *m, double kappa, int size_of_u)
{
    /* Declaration of variables */
    int i;
    
    /* Calculating the acceleration on the boundaries */
    a[0] = kappa*(- 2*u[0] + u[1])/m[0];
    a[size_of_u - 1] = kappa*(u[size_of_u - 2] - 2*u[size_of_u - 1])/m[size_of_u - 1];
    
    /* Calculating the acceleration of the inner points */
    for (i = 1; i < size_of_u - 1; i++){
        a[i] = kappa*(u[i - 1] - 2*u[i] + u[i + 1])/m[i];
    }
}

/*
 * Perform the velocity verlet alogrithm 
 * @n_timesteps - The number of time steps to be performed
 * @n_particles - number of particles in the system
 * @v - array of velocity (Empty allocated array) : sizeof(v) = n_particles
 * @q_n - position of the n'th atom : sizeof(q_n) = n_timesteps+1
 * @dt - timestep
 * @m - vector with masses of atoms sizeof(n_particles)
 * @kappa - Spring constant
 */
void velocity_verlet(int n_timesteps, int n_particles, double *v, double *q_1,
		     double *q_2, double *q_3, double dt, double *m,
		     double kappa)
{
    double q[n_particles]; //skapar en array med 3 platser
    double a[n_particles];
    q[0] = q_1[0];
    q[1] = q_2[0];
    q[2] = q_3[0];
    calc_acc(a, q, m, kappa, n_particles);
    for (int i = 1; i < n_timesteps + 1; i++) {
        /* v(t+dt/2) */
        for (int j = 0; j < n_particles; j++) {
            v[j] += dt * 0.5 * a[j];
        }
        
        /* q(t+dt) */
        for (int j = 0; j < n_particles; j++) {
            q[j] += dt * v[j];  //se gitlab sida 2 av readme, där visar dem ekvationerna fr velocity verlet och hur vi får dem. 
            //printf("%f",q[j]);
        }
        
        /* a(t+dt) */
        calc_acc(a, q, m, kappa, n_particles);
        
        /* v(t+dt) */
        for (int j = 0; j < n_particles; j++) {
            v[j] += dt * 0.5 * a[j];
        }
        
        /* Save the displacement of the three atoms */
        q_1[i] = q[0];
        q_2[i] = q[1];
        q_3[i] = q[2];
        printf("q_1[i]: %f,%f,%f\n",q_1[i],q_2[i], q_3[i]);
//de vi får ut här är q_1, q_2 och q_3 och vi har ju en nested loop som går 100 varv (100 timesteps) och varje varv så sparas ju positionerna i q[] arrayn,  där q[q_1, q_2, q_3]. 
    }    
}


int main() {

double* q;
double* q_1;   //here im trying to allocate three arrays for q1, q2 q3,   the positions of the three atoms. the size of q_n should be timesteps+1
double* q_2;   
double* q_3;
double* v;
int n_particles = 3;
double kappa= 1;
q = (double*)malloc(1000 * sizeof(int));
q_1 = (double*)malloc(1000 * sizeof(int));  //imorgon ska jag testa att använda C2 2D-array för att skapa matris för q1q2q3 
q_2 = (double*)malloc(1000 * sizeof(int)); 
q_3 = (double*)malloc(1000 * sizeof(int));
q_1[0] = 0.01;//initial positions 
q_2[0] = 0.005;
q_3[0] = -0.005;
double m[3] = {6.0, 6.0, 6.0}; //masses of the three particles size of n_particles, så size  = 3
v = (double*)calloc(n_particles,n_particles * sizeof(int)); //empty allocated array, size n_particles,  så size är 3


velocity_verlet(500,3,v,q_1,q_2,q_3,25, m, 0.0001);  //im wondering what the optimal parameters are... do we need to convert ? 


   FILE *fp = fopen("velocityverlet.csv", "w");
    for(int i = 0; i < 500; ++i){
	    fprintf(fp, "%f,%f,%f\n", q_1[i], q_2[i], q_3[i]);  //trying to add positions and time to a file so we can plot.  we get q[i] from velocity verlet. 
    }
    fclose(fp);
    









// how can we plot the energies? potential, kinetic and total energy.
//Total energy is the hamiltonian,  it should be a straight line. 

// H = sum(m*v[i])²+sum(kappa/2*(q[i+1]-q[i])
//double E_pot;
//for(int i=0;i<3;i++){
//	E_pot = 0.0001/2*(q_2[i]-q_1[i]);
//	printf("%f\n",E_pot);
//}





return 0;











}





/*


//double q_1[1] = {0.2};
//double q_2[1] = {0.1};
//double q_3[1] = {0.3};
//double q_1 = 0.1;
//double q_2 = 0.3;
//double q_3 = 0.2;


void arange(double *array, double start, int len_t, double dt){ //the function arange constructs a time array,  i took it from E1code1.c 
    for(int i = 0; i < len_t; i++){
	array[i] = start + i*dt;
    }
} 



*/












