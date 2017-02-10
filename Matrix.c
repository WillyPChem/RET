# include<math.h>
#include<stdio.h>
#include<stdlib.h>

int main() {

double H[2][2], D[2][2];

D[1][1] = 3;
D[0][1] = 1;

H[1][1] = 1.0;
H[0][0] = 1.5;

for (int i=0; i<2; i++) {

  for (int j=0; j<2; j++) {

    printf("  %f  ",H[i][j]);

  }
  printf(" \n ");
}


}

