#include<math.h>
#include<stdio.h>
#include<stdlib.h>
#include<malloc.h>

<<<<<<< HEAD
void Commutator (int dim, double H[dim][dim], double D[dim][dim], double P[dim][dim]);
=======
void Commutator(int dim, double *H, double *D, double *P);
>>>>>>> 6429e597c5a2793ad4511947767feccf9a440da7
int main() {

int dim = 2;
double *H, *D, *P;

H = (double *)malloc(dim*dim*sizeof(double));
D = (double *)malloc(dim*dim*sizeof(double));
P = (double *)malloc(dim*dim*sizeof(double));

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

D[1*dim + 1] = 3;
D[0*dim + 1] = 1;

H[1*dim + 1] = 1;
H[0*dim + 0] = 1.5;


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

<<<<<<< HEAD
void Commutator (int dim, double H[dim][dim], double D[dim][dim], double P[dim][dim]) {

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
      

=======
void Commutator(int dim, double *H, double *D, double *P) {

// write code here!
for (int i=0; i<dim; i++) {

  for (int j=0; j<dim; j++) {
>>>>>>> 6429e597c5a2793ad4511947767feccf9a440da7

    double sum2 = 0;
    double sum1 = 0;
    for (int k=0; k<dim; k++) {

      sum1 += H[i*dim+k]*D[k*dim+j];
      sum2 += D[i*dim+k]*H[k*dim+j];
    }
    P[i*dim+j] = sum1 - sum2;
    printf(" Pb[%i][%i] is %f\n",i,j,sum1-sum2);
}
}
}
}
}

