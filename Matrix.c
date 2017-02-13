#include<math.h>
#include<stdio.h>
#include<stdlib.h>

int main() {

double H[2][2], D[2][2], Pa[2][2], Pb[2][2];

D[1][1] = 3;
D[0][1] = 1;

H[1][1] = 1.0;
H[0][0] = 1.5;

for (int i=0; i<2; i++) {

  for (int j=0; j<2; j++) {

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

