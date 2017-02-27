#include<stdio.h>
#include<stdlib.h>

int dim = 5;
int array[dim][dim];

void print2(int n, int m);
int main() {

return 0;

}

void print2(int n, int m) {
  int i, j;
  for (i=0; i<n; i++) {

    for (j=0; j<m; j++) {

      printf("  %i \n",array[i][j]);

    }
  }


}
