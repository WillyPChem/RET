#include<math.h>
#include<stdio.h>
#include<stdlib.h>
#include<malloc.h>
#include<complex.h>

//void Commutator (int dim, double H[dim][dim], double D[dim][dim], double P[dim][dim]);
void RK3(int dim, double *xvec, double complex *wfn, double dx, double dt);
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


Commutator( dim, H, D, P);

for (int i=0; i<dim; i++) {
  for (int j=0; j<dim; j++) {

    printf("  %f  ",P[i*dim+j]);

  }
  printf(" \n ");

}

}

void RK3(int dim, double *xvec, double complex *wfn, double dx, double dt) {

  int i;
  double complex *wfn_dot, *wfn2, *wfn3, *wfn_np1, *k1, *k2, *k3;

  wfn_dot = (double complex *)malloc((dim+1)*sizeof(double complex));
  wfn2 = (double complex *)malloc((dim+1)*sizeof(double complex));
  wfn3 = (double complex *)malloc((dim+1)*sizeof(double complex));
  wfn_np1 = (double complex *)malloc((dim+1)*sizeof(double complex));
  k1 = (double complex *)malloc((dim+1)*sizeof(double complex));
  k2 = (double complex *)malloc((dim+1)*sizeof(double complex));
  k3 = (double complex *)malloc((dim+1)*sizeof(double complex));
  
  // Must zero out all elements of these arrays
  for (i=0; i<=dim; i++) {
    wfn_dot[i] = 0. + 0.*I;
    wfn2[i] = 0. + 0.*I;
    wfn3[i] = 0. + 0.*I;
    wfn_np1[i] = 0. + 0.*I;
    k1[i] = 0. + 0.*I;
    k2[i] = 0. + 0.*I;
    k3[i] = 0. + 0.*I;

  }

  // Get dPsi(n)/dt at initial time!
  void Commutator (int dim, double complex wfn, double complex wfn_dot, double dx);
  // Compute approximate wfn update with Euler step
  for (i=0; i<=dim; i++) {
    k1[i] = dt*wfn_dot[i];
    wfn2[i] = wfn[i] + k1[i]/2.;
  }
  // Get dPsi(n+k1/2)/dt
  void Commutator (int dim, double complex wfn2, double complex wfn_dot, double dx);
  // Compute approximate wfn update with Euler step
  for (i=0; i<=dim; i++) {
    k2[i] = dt*wfn_dot[i];
    wfn3[i] = wfn[i] + k2[i]/2.;
  }
  // Get dPsi(n+k2/2)/dt
  void Commutator (int dim, double complex  wfn3, double complex wfn_dot, double dx);
  // Compute approximate update with Euler step
  for (i=0; i<=dim; i++) {
    k3[i] = dt*wfn_dot[i];
    wfn_np1[i] = wfn[i] + k1[i]/6. + 2.*k2[i]/3. + k3[i]/6.;
    wfn[i] = wfn_np1[i];
}
  free(wfn_dot);
  free(wfn2);
  free(wfn3);
  free(wfn_np1);
  free(k1);
  free(k2);
  free(k3);

}





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



