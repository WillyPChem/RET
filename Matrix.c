#include<math.h>
#include<stdio.h>
#include<stdlib.h>

void Commutator(int dim, double H[dim][dim], double D[dim][dim], double P[dim][dim]);
int main() {

int dim = 5;
double H[dim][dim], D[dim][dim], Pa[dim][dim], Pb[dim][dim];

D[1][1] = 3;
D[0][1] = 1;

H[1][1] = 1.0;
H[0][0] = 1.5;

for (int i=0; i<dim; i++) {

  for (int j=0; j<dim; j++) {

    double sum2 = 0;
    double sum1 = 0; 
    for (int k=0; k<2; k++) {
   
      sum2 += H[i][k]*D[k][j];
      sum1 += D[i][k]*H[k][j];
    }
    Pa[i][j] = sum1;
      printf("  Pa[%i][%i] is %f\n",i,j,sum1);
   
     Pb[i][j] = sum2;
      printf(" Pb[%i][%i] is %f\n",i,j,sum2);  
 }
   
}

}

void Commutator(int dim, double H[dim][dim], double D[dim][dim], double P[dim][dim]) {

// write code here!

}
