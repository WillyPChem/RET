#include<math.h>
#include<stdio.h>
#include<stdlib.h>
#include<malloc.h>
#include<complex.h>

//void Commutator (int dim, double H[dim][dim], double D[dim][dim], double P[dim][dim]);
void RK3(int Nlevel, double time, double *E, double *Mu, double *Dis, double complex *D, double dt);
void Liouville(int dim, double complex *H, double complex *D, double complex *P);
void AntiCommutator(int dim, double *H, double complex *D, double complex *P);
void PrintComplexMatrix(int dim, double complex *M);
void PrintRealMatrix(int dim, double *M);
void FormDM(int dim, double complex *C, double complex *Ct, double complex *DM);
double E_Field(double time);

// NOTE!!!  You need three global variables for the rates associated with 
// The Lindblad operators - gamma, beta, and alpha should be defined here according
// to their values in the journal of physical chemistry letters paper
int Nlevel = 3;
int dim = Nlevel*Nlevel;

double pi = 4*atan(1.0);
double wn_to_au = 4.55633528e-6; 

int main() {

  
  int dim = Nlevel*Nlevel;
  double *E, *Mu, *Dis;
  double complex *H, *D, *P;
  double dt = 0.01;

  H = (double complex*)malloc(dim*sizeof(double complex));
  D = (double complex *)malloc(dim*sizeof(double complex));
  P = (double complex *)malloc(dim*sizeof(double complex));
  E  = (double *)malloc(dim*sizeof(double));
  Mu = (double *)malloc(dim*sizeof(double));
  Dis = (double *)malloc(dim*sizeof(double));
  
  // Initialize Denistry matrix as a superposition state of energy eigenstate 1 and energy eigenstate 2
  // Density matrix element D(i,j) is accessed as D[i*Nlevel+j];
  D[0] = 1. + 0.*I;
  for (int i=1; i<dim; i++) D[i] = 0. + 0.*I;

  // Variable that is a file pointer for each data file:
  FILE *Efp, *mufp, *disfp;

  // Open each file for reading
  Efp = fopen("Energy.txt","r");
  mufp = fopen("Dipole.txt","r");
  disfp = fopen("Dissipation.txt","r");  
  
  double val;
  for (int i=0; i<dim; i++) {


       // Read from energy file and store to the energy matrix
       fscanf(Efp,"%lf",&val);
       E[i] = val;

       fscanf(mufp,"%lf",&val);
       Mu[i] = val;

       fscanf(disfp,"%lf",&val);
       Dis[i] = val;

  }

  printf("\nE\n");
  PrintRealMatrix(Nlevel, E);
  printf("\nMu\n");
  PrintRealMatrix(Nlevel,Mu);
  printf("\nDiss\n");
  PrintRealMatrix(Nlevel,Dis);

  double tr=0.;
  for (int i=1; i<400000; i++) {

    //void RK3(int Nlevel, double time, double *E, double *Mu, double *Dis, double complex *D, double dt)
    RK3(Nlevel, dt*i, E, Mu, Dis, D, dt);

    printf("\n %f ",dt*i);
    tr=0.;
    for (int j=0; j<Nlevel; j++) {

      printf(" %12.10e",creal(D[j*Nlevel+j]));
      tr+=creal(D[j*Nlevel+j]);
    }
    printf("\n #Trace is %12.10e",tr);
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

      printf(" (%12.10e,%12.10e) ",creal(M[i*dim+j]),cimag(M[i*dim+j]));

    }
    printf("\n");
  }
  printf("\n");
}
void RK3(int Nlevel, double time, double *E, double *Mu, double *Dis, double complex *D, double dt) {

  int i, j;
  double complex *D_dot, *D2, *D3, *D_np1, *k1, *k2, *k3;
  double complex *H;  // Contribution to Ddot from Lindblad dissipation
  int dim = Nlevel*Nlevel; 
  double Efield;

  D_dot = (double complex *)malloc(dim*sizeof(double complex));
  D2    = (double complex *)malloc(dim*sizeof(double complex));
  D3    = (double complex *)malloc(dim*sizeof(double complex));
  D_np1 = (double complex *)malloc(dim*sizeof(double complex));
  k1    = (double complex *)malloc(dim*sizeof(double complex));
  k2    = (double complex *)malloc(dim*sizeof(double complex));
  k3    = (double complex *)malloc(dim*sizeof(double complex));
  H     = (double complex *)malloc(dim*sizeof(double complex));

  
  // Must zero out all elements of these arrays
  for (i=0; i<dim; i++) {
    D_dot[i] = 0. + 0.*I;
    D2[i] = 0. + 0.*I;
    D3[i] = 0. + 0.*I;
    D_np1[i] = 0. + 0.*I;
    k1[i] = 0. + 0.*I;
    k2[i] = 0. + 0.*I;
    k3[i] = 0. + 0.*I;
    H[i] = 0. + 0.*I;
   
  }

  Efield = E_Field(time);

  // Compute full Hamiltonian at current time t 
  for (i=0; i<dim; i++) {


      H[i] = E[i] + Efield*Mu[i] + I*Dis[i];

  } 

  //PrintComplexMatrix(Nlevel, H);

  // Get dPsi(n)/dt at initial time!
  Liouville(Nlevel, H, D, D_dot);

  //PrintComplexMatrix(Nlevel, D);
  //PrintComplexMatrix(Nlevel, D_dot);
  // Compute approximate wfn update with Euler step
  for (i=0; i<dim; i++) {
    k1[i] = dt*D_dot[i];
    D2[i] = D[i] + k1[i]/2.;
  }

  // Update Field!
  Efield = E_Field(time+dt/2.);

  // Compute full Hamiltonian at partially updated time t 
  for (i=0; i<dim; i++) {

      H[i] = E[i] + Efield*Mu[i] + I*Dis[i];

  }

  //PrintComplexMatrix(Nlevel, H);
  // Get dPsi(n+k1/2)/dt
  Liouville (Nlevel, H, D2, D_dot);
  // Compute approximate wfn update with Euler step
  for (i=0; i<dim; i++) {
    k2[i] = dt*D_dot[i];
    D3[i] = D[i] + k2[i]/2.;
  }

  // Get dPsi(n+k2/2)/dt
  Liouville (Nlevel, H, D3, D_dot);
  // Compute approximate update with Euler step
  for (i=0; i<dim; i++) {
    k3[i] = dt*D_dot[i];
    D_np1[i] = D[i] + k1[i]/6. + 2.*k2[i]/3. + k3[i]/6.;
    D[i] = D_np1[i];
}


  free(D_dot);
  free(D2);
  free(D3);
  free(D_np1);
  free(k1);
  free(k2);
  free(k3);
  free(H);
}

void FormDM(int dim, double complex *C, double complex *Ct, double complex *DM) {

  for (int i=0; i<dim; i++) {

    for (int j=0; j<dim; j++) {

      DM[i*dim+j] = C[i]*Ct[j];

    }

  }

}

void Liouville(int dim, double complex *H, double complex *D, double complex *P) {


  // write code here!
  for (int i=0; i<dim; i++) {

    for (int j=0; j<dim; j++) {


      double complex sum2 = 0.+0.*I;
      double complex sum1 = 0.+0.*I;
      for (int k=0; k<dim; k++) {

        sum1 -= H[i*dim+k]*D[k*dim+j]*I;
        sum2 += D[i*dim+k]*H[k*dim+j]*I;
      }
      P[i*dim+j] = sum1 + sum2;
    }
  }
}



void AntiCommutator(int dim, double *H, double complex *D, double complex *P) {

// write code here!
for (int i=0; i<dim; i++) {

  for (int j=0; j<dim; j++) {


    double complex sum2 = 0.+0.*I;
    double complex sum1 = 0.+0.*I;
    for (int k=0; k<dim; k++) {

      sum1 += H[i*dim+k]*D[k*dim+j];
      sum2 += D[i*dim+k]*H[k*dim+j];
    }
    P[i*dim+j] = sum1 + sum2;
    //printf(" Pb[%i][%i] is %f %f\n",i,j,creal(sum1-sum2),cimag(sum1-sum2));
}
}
}


double E_Field(double time) {

  double Ef;
  if (time<415.) {

    Ef = 0.001*sin(time*pi/415)*sin(time*pi/415)*sin(0.07423*time);

  }
  else Ef = 0.;
    
  return Ef;


}

