#include<math.h>
#include<stdio.h>
#include<stdlib.h>
#include<malloc.h>
#include<complex.h>

//void Commutator (int dim, double H[dim][dim], double D[dim][dim], double P[dim][dim]);
void RK3(int dim, double *xvec, double complex *wfn, double dx, double dt);
void Commutator(int dim, double *H, double complex *D, double complex *P);
void PrintComplexMatrix(int dim, double complex *M);
void PrintRealMatrix(int dim, double *M);
void FormDM(int dim, double complex *C, double complex *Ct, double complex *DM);

double pi = 3.14159625;
int main() {

int dim = 2.;
double *H;
double complex *D, *P, *C, *Ct, *Dp;
double dt = .001;
double dx = 1.;

H = (double *)malloc(dim*dim*sizeof(double));
D = (double complex *)malloc(dim*dim*sizeof(double complex));
P = (double complex *)malloc(dim*dim*sizeof(double complex));
Dp= (double complex *)malloc(dim*dim*sizeof(double complex));
C = (double complex *)malloc(dim*sizeof(double complex));
Ct= (double complex *)malloc(dim*sizeof(double complex));


// Element DM(0,0) = D(0*dim+0)
// Element DM(1,0) = D(1*dim+0)
// Element DM(0,1) = D(0*dim+1)
// In general, DM(i,j) = D(i*dim+j) where DM is the matrix
// form of D and D is a vector of length dim*dim
//

C[0] = sqrt(1/2.) + 0.*I;
C[1] = sqrt(1/2.) + 0.*I;
Ct[0] = C[0];
Ct[1] = C[1];

D[0*dim + 0] = 1/2. + 0.*I;
D[0*dim + 1] = 1/2. + 0.*I;
D[1*dim + 0] = 1/2. + 0.*I;
D[1*dim + 1] = 1/2. + 0.*I;

H[0*dim + 0] = (pi*pi/(2*100));
H[1*dim + 1] = 4*pi*pi/(2*100);


//Commutator( dim, H, D, P);
printf("  hard coded DM \n");
PrintComplexMatrix(dim, D);
FormDM(dim, C, Ct, Dp);
printf("  wfn DM \n");
PrintComplexMatrix(dim, Dp);


for (int i=1; i<500; i++) {

  RK3(dim, H, D, dx, dt);

  for (int j=0; j<dim; j++) {

    C[j]  = (sqrt(1/2.) + 0.*I)*cexp(-I*H[j*dim+j]*dt*i);
    Ct[j] = (sqrt(1/2.) + 0.*I)*cexp(I*H[j*dim+j]*dt*i); 

  }
  FormDM(dim, C, Ct, Dp);
  printf("  Updated from Commutator \n");
  PrintComplexMatrix(dim, D);
  printf("  Updated Analytically \n");
  PrintComplexMatrix(dim, Dp);

}
}

void PrintRealMatrix(int dim, double *M) {

  printf("\n");
  for (int i=0; i<dim; i++) {

    for (int j=0; j<dim; j++) {

      printf(" %f ",M[i*dim+j]);

    }
    printf("\n");
  }

  printf("\n");
}

void PrintComplexMatrix(int dim, double complex *M) {
 
  printf("\n");
  for (int i=0; i<dim; i++) {

    for (int j=0; j<dim; j++) {

      printf(" (%f,%f) ",creal(M[i*dim+j]),cimag(M[i*dim+j]));

    }
    printf("\n");
  }
  printf("\n");
}
void RK3(int dim, double *H, double complex *wfn, double dx, double dt) {

  int i, j;
  double complex *wfn_dot, *wfn2, *wfn3, *wfn_np1, *k1, *k2, *k3;
  
  wfn_dot = (double complex *)malloc((dim*dim)*sizeof(double complex));
  wfn2 = (double complex *)malloc((dim*dim)*sizeof(double complex));
  wfn3 = (double complex *)malloc((dim*dim)*sizeof(double complex));
  wfn_np1 = (double complex *)malloc((dim*dim)*sizeof(double complex));
  k1 = (double complex *)malloc((dim*dim)*sizeof(double complex));
  k2 = (double complex *)malloc((dim*dim)*sizeof(double complex));
  k3 = (double complex *)malloc((dim*dim)*sizeof(double complex));
  
  // Must zero out all elements of these arrays
  for (i=0; i<dim; i++) {
   for (j=0; j<dim; j++) {
    wfn_dot[i*dim+j] = 0. + 0.*I;
    wfn2[i*dim+j] = 0. + 0.*I;
    wfn3[i*dim+j] = 0. + 0.*I;
    wfn_np1[i*dim+j] = 0. + 0.*I;
    k1[i*dim+j] = 0. + 0.*I;
    k2[i*dim+j] = 0. + 0.*I;
    k3[i*dim+j] = 0. + 0.*I;
    }
  }

  // Get dPsi(n)/dt at initial time!
  Commutator (dim, H, wfn, wfn_dot);
  // Compute approximate wfn update with Euler step
  for (i=0; i<dim; i++) {
    for (j=0; j<dim; j++) {
    k1[i*dim+j] = dt*wfn_dot[i*dim+j];
    wfn2[i*dim+j] = wfn[i*dim+j] + k1[i*dim+j]/2.;
   }

  }
  // Get dPsi(n+k1/2)/dt
  Commutator (dim, H, wfn2, wfn_dot);
  // Compute approximate wfn update with Euler step
  for (i=0; i<dim; i++) {
    for (j=0; j<dim; j++) {
    k2[i*dim+j] = dt*wfn_dot[i*dim+j];
    wfn3[i*dim+j] = wfn[i*dim+j] + k2[i*dim+j]/2.;
  }
  }
  // Get dPsi(n+k2/2)/dt
  Commutator (dim,H, wfn3, wfn_dot);
  // Compute approximate update with Euler step
  for (i=0; i<dim; i++) {
    for (j=0; j<dim; j++) {
    k3[i*dim+j] = dt*wfn_dot[i*dim+j];
    wfn_np1[i*dim+j] = wfn[i*dim+j] + k1[i*dim+j]/6. + 2.*k2[i*dim+j]/3. + k3[i*dim+j]/6.;
    wfn[i*dim+j] = wfn_np1[i*dim+j];
  }
}
  free(wfn_dot);
  free(wfn2);
  free(wfn3);
  free(wfn_np1);
  free(k1);
  free(k2);
  free(k3);

}

void FormDM(int dim, double complex *C, double complex *Ct, double complex *DM) {

  for (int i=0; i<dim; i++) {

    for (int j=0; j<dim; j++) {

      DM[i*dim+j] = C[i]*Ct[j];

    }

  }

}



void Commutator(int dim, double *H, double complex *D, double complex *P) {

// write code here!
for (int i=0; i<dim; i++) {

  for (int j=0; j<dim; j++) {


    double complex sum2 = 0+0*I;
    double complex sum1 = 0+0*I;
    for (int k=0; k<dim; k++) {

      sum1 -= H[i*dim+k]*D[k*dim+j]*I;
      sum2 += D[i*dim+k]*H[k*dim+j]*I;
    }
    P[i*dim+j] = sum1 + sum2;
    //printf(" Pb[%i][%i] is %f %f\n",i,j,creal(sum1-sum2),cimag(sum1-sum2));
}
}
}



