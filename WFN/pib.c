#include<stdlib.h>
#include<malloc.h>
#include<stdio.h>
#include<math.h>
#include<complex.h>

double pi = 4.*atan(1.0);

void FiniteDifference_2D(int dim, double complex *f, double complex *dfdt, double dx);
void PrintVec(int dim, double *xvec, double complex *fvec);

void RK3(int dim, double *xvec, double complex *wfn, double dx, double dt);

int main() {

  int dim = 100;
  double L = 1.;
  double *x;
  double complex *wfn, *wfn_dot, *an;

  wfn = (double complex *)malloc((dim+1)*sizeof(double complex));
  wfn_dot = (double complex *)malloc((dim+1)*sizeof(double complex));
  an  = (double complex *)malloc((dim+1)*sizeof(double complex));
  x   = (double *)malloc((dim+1)*sizeof(double));
  double dx = L/dim;
  double dt = 0.000005;

  for (int i=0; i<=dim; i++) {

    x[i] = dx*i;
    wfn[i] = sqrt(2./L)*sin(5*pi*dx*i/L) + 0.*I;
  }

  int piter=1;
  for (int i=0; i<100000; i++) {
    RK3(dim, x, wfn, dx, dt); 
    if (i%20==0) {
      printf("\n\n#%i\n",piter);
      PrintVec(dim, x, wfn);
      piter++;
    }
  }
  
}

void PrintVec(int dim, double *xvec, double complex *fvec) {

  for (int i=0; i<=dim; i++) {

    printf("  %f  %12.10f  %12.10f\n",xvec[i],creal(fvec[i]),cimag(fvec[i]));

  }

}

void FiniteDifference_2D(int dim, double complex *f, double complex *dfdt, double dx) {

  // Take care of first and last points of fpp
  for (int i=1; i<dim; i++) {
    dfdt[i] = 0. + 0.*I;
  }
  //  Loop through points 1-dim-2
  for (int i=1; i<dim; i++) {

    dfdt[i] = (I/2.)*(f[i+1] - 2.*f[i] + f[i-1])/(dx*dx);

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
  FiniteDifference_2D(dim, wfn, wfn_dot, dx);
  // Compute approximate wfn update with Euler step
  for (i=0; i<=dim; i++) {
    k1[i] = dt*wfn_dot[i];
    wfn2[i] = wfn[i] + k1[i]/2.;
  }
  // Get dPsi(n+k1/2)/dt
  FiniteDifference_2D(dim, wfn2, wfn_dot, dx);
  // Compute approximate wfn update with Euler step
  for (i=0; i<=dim; i++) {
    k2[i] = dt*wfn_dot[i];
    wfn3[i] = wfn[i] + k2[i]/2.;
  }
  // Get dPsi(n+k2/2)/dt
  FiniteDifference_2D(dim, wfn3, wfn_dot, dx);
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
