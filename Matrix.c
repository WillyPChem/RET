#include<math.h>
#include<stdio.h>
#include<stdlib.h>
#include<malloc.h>
#include<complex.h>

//void Commutator (int dim, double H[dim][dim], double D[dim][dim], double P[dim][dim]);

void Commutator(int dim, double *H, double complex *D, double complex *P);
double pi = 3.14159625;
int main() {

int dim = 2.;
double *H;
double complex *D, *P;

H = (double *)malloc(dim*dim*sizeof(double));
D = (double complex *)malloc(dim*dim*sizeof(double complex));
P = (double complex *)malloc(dim*dim*sizeof(double complex));

// Element DM(0,0) = D(0*dim+0)
// Element DM(1,0) = D(1*dim+0)
// Element DM(0,1) = D(0*dim+1)
// In general, DM(i,j) = D(i*dim+j) where DM is the matrix
// form of D and D is a vector of length dim*dim
//

/*D[1][1] = 3;
D[0][1] = 1;
H[1][1] = 1.0;
H[0][0] = 1.5;
*/

D[0*dim + 0] = 1/2.;
D[1*dim + 1] = 1/2.;

H[0*dim + 0] = (pi*pi/2);
H[1*dim + 1] = 2*pi*pi;


Commutator( dim, H, D,  P);

for (int i=0; i<dim; i++) {
  for (int j=0; j<dim; j++) {

    printf("  %f  ",P[i*dim+j]);

  }
  printf(" \n ");

}

/*
for (int i=0; i<dim; i++) {

  for (int j=0; j<dim; j++) {

    double sum2 = 0;
    double sum1 = 0; 
    for (int k=0; k<dim; k++) {
   
      sum2 += H[i][k]*D[k][j];
      sum1 += D[i][k]*H[k][j];
    }
    Pa[i][j] = sum1;
      printf("  Pa[%i][%i] is %f\n",i,j,sum1);
   
     Pb[i][j] = sum2;
      printf(" Pb[%i][%i] is %f\n",i,j,sum2);  
 }
   
}
*/

}


/*void Commutator (int dim, double H[dim][dim], double D[dim][dim], double P[dim][dim]) {

// write code here!
 for (int i=0; i<dim; i++) {

   for (int j=0; j<dim; j++) {

   for(int k=0; k<dim; ++k)
			{
      sum2 += H[i][k]*D[k][j];
      sum1 += D[i][k]*H[k][j];
    }
    Pa[i][j] = sum1;
      printf("  Pa[%i][%i] is %f\n",i,j,sum1);

     Pb[i][j] = sum2;
      printf(" Pb[%i][%i] is %f\n",i,j,sum2);
*/      


void Commutator(int dim, double *H, double complex *D, double complex *P) {

// write code here!
for (int i=0; i<dim; i++) {

  for (int j=0; j<dim; j++) {


    double complex sum2 = 0+0*I;
    double complex sum1 = 0+0*I;
    for (int k=0; k<dim; k++) {

      sum1 += H[i*dim+k]*D[k*dim+j]*I;
      sum2 -= D[i*dim+k]*H[k*dim+j]*I;
    }
    P[i*dim+j] = sum1 - sum2;
    printf(" Pb[%i][%i] is %f %f\n",i,j,creal(sum1-sum2),cimag(sum1-sum2));
}
}
}



