#include<stdlib.h>
#include<stdio.h>
#include<math.h>
//contains functions for using complex number
//I =  the imaginary function
#include<complex.h>
#include<malloc.h>
double pi = 3.14159265359;

// Function prototype for second derivative
// dfdt is going to take the wavefunction, take the second
// derivative at each point multiplied by the imaginary unit 
// and return the time-derivative of the wavefunction.
// Arguements: arg 1 - dim, number of points in the wf 
//             arg 2 - psivec - array of wave function values
//             arg 3 - dpsi - array of d/dt wave function values
//             arg 4 - dx  - increment along x axis
void dfdt(int dim, double complex *psivec, double complex *dpsi, double dx   );  





int main() {


  
  int dim = 10;
  double *x;
  double complex *wfn, *dpsi;
  double L = 1.;

  double dx;
  int i;

 dx = L/dim;

  x = (double *)malloc(dim*sizeof(double));
 wfn = (double complex *)malloc(dim*sizeof(double complex));
 dpsi= (double complex *)malloc(dim*sizeof(double complex));
  for (i=0; i<=dim; i++)  {


    x[i] = i*dx;
    wfn[i] = sqrt(2./L)*sin(pi*x[i]/L) + 0.*I;
    printf(" %f %f %f\n",x[i], creal(wfn[i]), cimag(wfn[i]));

  }

dfdt( dim, wfn, dpsi, dx);

  for (i=0; i<=dim; i++) {
  printf(" %f (%f, %f) (%f, %f)\n",x[i],creal(dpsi[i]),cimag(wfn[i]) *pi*pi/2., cimag(dpsi[i]), -1*creal(wfn[i])*pi*pi/2.);


  }

  
}




void dfdt(int dim, double complex *psivec, double complex *dpsi, double dx   ) {



//write code here
  int j;
  dpsi[0] = 0. + 0.*I;
  dpsi[dim] = 0. + 0.*I;

  for (j=1; j<dim; j++) {

    dpsi[j] = (I/2.)*(psivec[j+1] - 2*psivec[j] + psivec[j-1])/(dx*dx);




}

}

