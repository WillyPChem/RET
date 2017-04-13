#include<math.h>
#include<stdio.h>
#include<stdlib.h>
#include<malloc.h>
#include<complex.h>

//void Commutator (int dim, double H[dim][dim], double D[dim][dim], double P[dim][dim]);
void RK3(int dim, double *xvec, double complex *wfn, double dx, double dt);
void Commutator(int dim, double *H, double complex *D, double complex *P);
void AntiCommutator(int dim, double *H, double complex *D, double complex *P);
void L_Deph(int dim, double alpha, double complex *D, double *bas, double complex *P);
void L_Diss(int dim, int g, double beta, double complex *D, double *bas, double complex *P);
void L_Sink(int dim, int chan, int s, double gamma, double complex *D, double *bas, double complex *P);
double pi = 3.14159625;
int main() {

int dim = 9.;
double *H;
double complex *D, *P;
double dt = .1;
double dx = 1.;

H = (double *)malloc(dim*dim*sizeof(double));
D = (double complex *)malloc(dim*dim*sizeof(double complex));
P = (double complex *)malloc(dim*dim*sizeof(double complex));



// Element DM(0,0) = D(0*dim+0)
// Element DM(1,0) = D(1*dim+0)
// Element DM(0,1) = D(0*dim+1)
// In general, DM(i,j) = D(i*dim+j) where DM is the matrix
// form of D and D is a vector of length dim*dim
//

// Initialize Denistry matrix as a superposition state of energy eigenstate 1 and energy eigenstate 2
// D comes from Psi = sqrt(1/2) psi_1 + sqrt(1/2) psi_2
D[0*dim + 0] = 0 + 0.*I;
D[0*dim + 1] = 0 + 0.*I;
D[1*dim + 0] = 0 + 0.*I;
D[1*dim + 1] = 1. + 0.*I;

/*
D[0*dim + 0] = 0.;
D[1*dim + 1] = 0.;
D[2*dim + 2] = 0.; 
D[3*dim + 3] = 0.;
D[4*dim + 4] = 0.;
D[5*dim + 5] = 0.;
D[6*dim + 6] = 0.;
D[7*dim + 7] = 0.;
D[8*dim + 8] = 0.;
D[9*dim + 9] = 0.;
*/
H[0*dim + 0] = 0.;
H[0*dim + 1] = 0.;
H[0*dim + 2] = 0.;
H[0*dim + 3] = 0.;
H[0*dim + 4] = 0.; 
H[0*dim + 5] = 0.;
H[0*dim + 6] = 0.;
H[0*dim + 7] = 0.;
H[0*dim + 8] = 0.;
H[0*dim + 9] = 0.;

H[1*dim + 0] =  0.;
H[1*dim + 1] =  0.0012757;
H[1*dim + 2] = -0.0004829;
H[1*dim + 3] =  0.0000364;
H[1*dim + 4] = -0.0000227;
H[1*dim + 5] =  0.0000273;
H[1*dim + 6] = -0.0000364;
H[1*dim + 7] = -0.0000182;
H[1*dim + 8] =  0.;


H[2*dim + 0] = 0.;
H[2*dim + 1] = -0.0004829;
H[2*dim + 2] =  0.0019136;
H[2*dim + 3] =  0.0001275;
H[2*dim + 4] =  0.0000273;
H[2*dim + 5] =  0.0000091;
H[2*dim + 6] =  0.0000592;
H[2*dim + 7] =  0.0000045;
H[2*dim + 8] =  0.;

H[3*dim + 0] =  0.;
H[3*dim + 1] =  0.0000364;
H[3*dim + 2] =  0.0001275;
H[3*dim + 3] =  0.;
H[3*dim + 4] = -0.0002824; 
H[3*dim + 5] = -0.0000045;
H[3*dim + 6] = -0.0000410;
H[3*dim + 7] =  0.0000774;
H[3*dim + 8] =  0.; 

H[4*dim + 0] =  0.;
H[4*dim + 1] = -0.0000227;
H[4*dim + 2] =  0.0000273;
H[4*dim + 3] = -0.0002824;
H[4*dim + 4] =  0.0007973;
H[4*dim + 5] = -0.0003189;
H[4*dim + 6] = -0.0000865;
H[4*dim + 7] = -0.0002597;
H[4*dim + 8] =  0.;

H[5*dim + 0] =  0.;
H[5*dim + 1] =  0.0000273;
H[5*dim + 2] =  0.0000091;
H[5*dim + 3] = -0.0000045;
H[5*dim + 4] = -0.0003189;
H[5*dim + 5] =  0.0014580;
H[5*dim + 6] =  0.0001822;
H[5*dim + 7] = -0.0000091;
H[5*dim + 8] =  0.;

H[6*dim + 0] =  0.;
H[6*dim + 1] = -0.0000364;
H[6*dim + 2] =  0.0000592;
H[6*dim + 3] = -0.0000410;
H[6*dim + 4] = -0.0000865;
H[6*dim + 5] =  0.0001822;
H[6*dim + 6] =  0.0016402;
H[6*dim + 7] =  0.0001458;
H[6*dim + 8] = 0.;

H[7*dim + 0] =  0.;
H[7*dim + 1] = -0.0000182;
H[7*dim + 2] =  0.0000045;
H[7*dim + 3] =  0.0000774;
H[7*dim + 4] = -0.0002597;
H[7*dim + 5] = -0.0000091;
H[7*dim + 6] =  0.0001458;
H[7*dim + 7] =  0.0011846;
H[7*dim + 8] =  0.;

H[8*dim + 0] = 0.;
H[8*dim + 1] = 0.;
H[8*dim + 2] = 0.;
H[8*dim + 3] = 0.;
H[8*dim + 4] = 0.;
H[8*dim + 5] = 0.;
H[8*dim + 6] = 0.;
H[8*dim + 7] = 0.;
H[8*dim + 8] = 0.;

//RK3( dim, H, D, P);
// Hamiltonian for particle in a box of length 10 Bohr radii
//H[0*dim + 0] = pi*pi/(2*100);
//H[1*dim + 1] = 4*pi*pi/(2*100);


//Commutator( dim, H, D, P);
/*
for (int i=0; i<dim; i++) {
  for (int j=0; j<dim; j++) {

  printf(" (%f , %f) ",creal(D[i*dim+j]), cimag(D[i*dim+j]));
    
  }
  printf("\n");
} */
//void RK3(int dim, double *xvec, double complex *wfn, double dx, double dt)


for (int j=0; j<500000; j++) {

  

  RK3(dim, H, D, dx, dt);
  printf(" %f ",j*dt);
  for (int i=0; i<dim; i++) {
  //for (int j=0; j<dim; j++) {

    printf(" %e ",creal(D[i*dim+i]));
//}
  }
  printf("\n");
 }

/*
RK3(dim, H, D, dx, dt);
for (int i=0; i<dim; i++) {
  for (int j=0; j<dim; j++) {

  printf(" (%f , %f) ",creal(D[i*dim+j]), cimag(D[i*dim+j]));

  }
  printf("\n");
}
*/
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

// Function that takes the anticommutator of H and D -> i.e. {H,D} = HD + DH
// In this case, H is taken to be a real matrix and D is taken to be a complex matrix
// The result of the commutator will be a complex matrix
void AntiCommutator(int dim, double *H, double complex *D, double complex *P) {

  for (int i=0; i<dim; i++) {

    for (int j=0; j<dim; j++) {

      double complex sum2 = 0+0*I;
      double complex sum1 = 0+0*I;
      for (int k=0; k<dim; k++) {

        sum1 += H[i*dim+k]*D[k*dim+j];
        sum2 += D[i*dim+k]*H[k*dim+j];
    
      }
    
      P[i*dim+j] = sum1 + sum2;
  
    }
  }
}

// Lioville operator that controls the transfer of exciton popolution to the reaction center.
// The argument of this functions are as follows:
// int dim -> number of states including ground state, exciton states, and sink state... full
//            FMO model + reaction center will have dim = 9
// int chan -> the index of the exciton state the couples directly to the reaction center, aka, the channel.  In the FMO model
//             described in J Chem Phys Lett, chromophore 3 couples to the reaction center, so 
//             so chan=3 in that case 
// int s    -> the index of the state the corresponds to the reaction center (aka sink)... in the FMO model
//             described in the J Chem Phys Lett, the sink has index = 8
// double gamma -> a rate constant which describes rate of transfer from the chan chromophore to the sink state
// double complex *D -> the current density matrix
// double *bas -> the basis states for the density matrix
// double complex *P -> The contribution to the time-derivative of the Density matrix coming from coupling to the sink
void L_Sink(int dim, int chan, int s, double gamma, double complex *D, double *bas, double complex *P) {

  int i, j, k;
  double *chan_bas, *s_bas;
  double complex *temp_t1, *temp_t2, *LD;
  chan_bas = (double *)malloc(dim*dim*sizeof(double));
  s_bas    = (double *)malloc(dim*dim*sizeof(double));
  temp_t1  = (double complex *)malloc(dim*dim*sizeof(double complex));
  temp_t2  = (double complex *)malloc(dim*dim*sizeof(double complex));
  LD       = (double complex *)malloc(dim*dim*sizeof(double complex));



  // Form |chan><chang| and |s><s| matrices
  for (i=0; i<dim; i++) {
    for (j=0; j<dim; j++) {
      chan_bas[i*dim+j] = bas[chan*dim+i]*bas[j*dim+chan];
      s_bas[i*dim+j] = bas[s*dim+i]*bas[j*dim+s];
    }
  }



  for (i=0; i<dim; i++) {

    for (j=0; j<dim; j++) {

      temp_t1[i*dim+j] = 2*D[chan*dim+chan]*s_bas[i*dim+j];
    }
  }
  AntiCommutator(dim, chan_bas, D, temp_t2);

  for (i=0; i<dim; i++) {
    for (j=0; j<dim; j++) {
      P[i*dim+j] = gamma*temp_t1[i*dim+j] - gamma*temp_t2[i*dim+j];
    }
  }

 free(chan_bas);
 free(s_bas);
 free(temp_t1);
 free(temp_t2);
 free(LD);

}


