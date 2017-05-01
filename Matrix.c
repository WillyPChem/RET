#include<math.h>
#include<stdio.h>
#include<stdlib.h>
#include<malloc.h>
#include<complex.h>

//void Commutator (int dim, double H[dim][dim], double D[dim][dim], double P[dim][dim]);
void RK3(int dim, double *bas, double *xvec, double complex *D, double dx, double dt);
void Commutator(int dim, double *H, double complex *D, double complex *P);
void AntiCommutator(int dim, double *H, double complex *D, double complex *P);
void L_Deph(int dim, double alpha, double complex *D, double *bas, double complex *P);
void L_Diss(int dim, int g, double beta, double complex *D, double *bas, double complex *P);
void L_Sink(int dim, int chan, int s, double gamma, double complex *D, double *bas, double complex *P);
double gam = 1.21e-8/2;
double bet = 7.26e-5/2;
double alph = 1.52e-4/2;
double pi = 3.14159625;

// NOTE!!!  You need three global variables for the rates associated with 
// The Lindblad operators - gamma, beta, and alpha should be defined here according
// to their values in the journal of physical chemistry letters paper
int dim = 9;
int g = 0;
int s = 8;
int chan = 3;

int main() {

double *H, *bas;
double complex *D, *P;
double dt = 1.;
double dx = 1.;

H = (double *)malloc(dim*dim*sizeof(double));
D = (double complex *)malloc(dim*dim*sizeof(double complex));
P = (double complex *)malloc(dim*dim*sizeof(double complex));
bas = (double *)malloc(dim*dim*sizeof(double));

for (int i=0; i<dim; i++) {
  for (int j=0; j<dim; j++) {
    if (i==j) bas[i*dim+i] = 1.;
    else bas[i*dim+j] = 0.;
  }
}

// Element DM(1,0) = D(1*dim+0)
// D comes from Psi = sqrt(1/2) psi_1 + sqrt(1/2) psi_2
// Element DM(0,1) = D(0*dim+1)
// In general, DM(i,j) = D(i*dim+j) where DM is the matrix
// form of D and D is a vector of length dim*dim
//

// Initialize Denistry matrix as a superposition state of energy eigenstate 1 and energy eigenstate 2
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


for (int j=0; j<50000; j++) {

  

  RK3(dim, bas, H, D, dx, dt);
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

// Function that Solves the Liouville equation to update the densisty matrix D from D(t) to D(t+dt)
// Using 3rd Order Runge-Kutta (RK3)... Must compute the time-derivative of D(t) and two different
// partially updated D(ts) [D(t + k1/2) and D(t+k2/2)]
void RK3(int dim, double *bas, double *H, double complex *D, double dx, double dt) {

  int i, j;
  double complex *D_dot, *D2, *D3, *D_np1, *k1, *k2, *k3, *LDiss, *LDeph, *LSink;
  // NOTE!!!!  Need to declare three additional complex double arrays that will store the three contributions 
  // to Ddot (time derivative of the Density matrix) that arise from the Lindblad operator (L_Diss, L_Deph, and L_Sink)
  // Declare those arrays here!
  
  // Time-derivative of Density Matrix from Commutator
  D_dot = (double complex *)malloc((dim*dim)*sizeof(double complex));
  // Partially updated Density Matrix from first RK step
  D2 = (double complex *)malloc((dim*dim)*sizeof(double complex));
  // Partially updated Density Matrix from second RK step
  D3 = (double complex *)malloc((dim*dim)*sizeof(double complex));
  // Temporary Array for fully updated Density Matrix
  D_np1 = (double complex *)malloc((dim*dim)*sizeof(double complex));
  // Three "stride" arrays for the RK3 updates
  k1 = (double complex *)malloc((dim*dim)*sizeof(double complex));
  k2 = (double complex *)malloc((dim*dim)*sizeof(double complex));
  k3 = (double complex *)malloc((dim*dim)*sizeof(double complex));
  
  // NOTE!!! Need to allocate Meomry for the three arrays associated with the Lindblad operators
  // Allocate Memory here!
  LDiss = (double complex *)malloc((dim*dim)*sizeof(double complex));
  LDeph = (double complex *)malloc((dim*dim)*sizeof(double complex));
  LSink = (double complex *)malloc((dim*dim)*sizeof(double complex));
  // Must zero out all elements of these arrays
  for (i=0; i<dim; i++) {
   for (j=0; j<dim; j++) {
    D_dot[i*dim+j] = 0. + 0.*I;
    D2[i*dim+j] = 0. + 0.*I;
    D3[i*dim+j] = 0. + 0.*I;
    D_np1[i*dim+j] = 0. + 0.*I;
    k1[i*dim+j] = 0. + 0.*I;
    k2[i*dim+j] = 0. + 0.*I;
    k3[i*dim+j] = 0. + 0.*I;
    LDiss[i*dim+j] = 0. + 0.*I;
    LDeph[i*dim+j] = 0. + 0.*I;
    LSink[i*dim+j] = 0. + 0.*I;
  }
  }

  // Get d/dt of D(t)
   L_Deph(dim, alph, D, bas, LDeph);
   L_Diss(dim, g, bet, D, bas, LDiss);
   L_Sink(dim, chan, s, gam, D, bas, LSink);
   Commutator (dim, H, D, D_dot);
  // NOTE!!!  Need to get the Lindblad contributions to D_dot as well
  // Call the three Lindblad Functions with the appropriate arguments here!!!
  // Compute approximate wfn update with Euler step
  for (i=0; i<dim; i++) {
    for (j=0; j<dim; j++) {
    // NOTE!!! Partial update needs to include contributions from Lindblad terms!
    // rewrite the right hand side of the "k1[i*dim+j] =" line so that these terms are included 
    k1[i*dim+j] = dt*D_dot[i*dim+j] + dt*LDeph[i*dim+j] + dt*LDiss[i*dim+j] + dt*LSink[i*dim+j];
    D2[i*dim+j] = D[i*dim+j] + k1[i*dim+j]/2.;
   }

  }
  // Get d/dt of partially-updated D(t)+k1/2
   L_Deph(dim, alph, D2, bas, LDeph);
   L_Diss(dim, g, bet, D2, bas, LDiss);
   L_Sink(dim, chan, s, gam, D2, bas, LSink);
   Commutator (dim, H, D2, D_dot);
  // NOTE!!!  Need to get the Lindblad contributions to D_dot as well
  // Call the three Lindblad Functions with the appropriate arguments here!!!
  // Compute approximate wfn update with Euler step
  for (i=0; i<dim; i++) {
    for (j=0; j<dim; j++) {
    // NOTE!!! Partial update needs to include contributions from Lindblad terms!
    // rewrite the right hand side of the "k2[i*dim+j] =" line so that these terms are included 
    k2[i*dim+j] = dt*D_dot[i*dim+j]+ dt*LDeph[i*dim+j] + dt*LDiss[i*dim+j] + dt*LSink[i*dim+j];
    D3[i*dim+j] = D[i*dim+j] + k2[i*dim+j]/2.;
  }
  }
  // Get d/dt of partially-updated D(t)+k2/2
 L_Deph(dim, alph, D3, bas, LDeph);
 L_Diss(dim, g, bet, D3, bas, LDiss);
 L_Sink(dim, chan, s, gam, D3, bas, LSink);
 Commutator (dim,H, D3, D_dot);
  // NOTE!!!  Need to get the Lindblad contributions to D_dot as well
  // Call the three Lindblad Functions with the appropriate arguments here!!!
  // Compute approximate update with Euler step
  for (i=0; i<dim; i++) {
    for (j=0; j<dim; j++) {
    // NOTE!!! Partial update needs to include contributions from Lindblad terms!
    // rewrite the right hand side of the "k3[i*dim+j] =" line so that these terms are included 
    k3[i*dim+j] = dt*D_dot[i*dim+j]+ dt*LDeph[i*dim+j] + dt*LDiss[i*dim+j] + dt*LSink[i*dim+j];
    // Complete RK3 update of D(t) -> D(t+dt)
    D_np1[i*dim+j] = D[i*dim+j] + k1[i*dim+j]/6. + 2.*k2[i*dim+j]/3. + k3[i*dim+j]/6.;
    D[i*dim+j] = D_np1[i*dim+j];
  }
}
  // Free Memory for all temporary arrays!
  // NOTE!!! You need to free memory for the new arrays you created for the Lindblad operator
  free(D_dot);
  free(D2);
  free(D3);
  free(D_np1);
  free(k1);
  free(k2);
  free(k3);
  free(LDiss);
  free(LDeph);
  free(LSink);
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

// Lioville operator that controls the transfer of exciton popolution back to the ground state
// The argument of this functions are as follows:
// int dim -> number of states including ground state, exciton states, and sink state... full
//            FMO model + reaction center will have dim = 9
// int g -> the index of the ground state - in this model, it is g = 0
// double beta -> a rate constant which describes rate of transfer from exciton states on chromophores to the ground state
// double complex *D -> the current density matrix
// double *bas -> the basis states for the density matrix
// double complex *P -> The contribution to the time-derivative of the Density matrix coming from relaxation to the ground state
void L_Diss(int dim, int g, double beta, double complex *D, double *bas, double complex *P) {

  int i, j, k;
  double *temp_bas, *g_bas;
  double complex *temp_t1, *temp_t2, *LD;
  temp_bas = (double *)malloc(dim*dim*sizeof(double));
  g_bas    = (double *)malloc(dim*dim*sizeof(double));
  temp_t1  = (double complex *)malloc(dim*dim*sizeof(double complex));
  temp_t2  = (double complex *)malloc(dim*dim*sizeof(double complex));
  LD       = (double complex *)malloc(dim*dim*sizeof(double complex));

  // Form |g><g| matrix
  for (i=0; i<dim; i++) {
    for (j=0; j<dim; j++) {
      g_bas[i*dim+j] = bas[g*dim+i]*bas[j*dim+g];
      LD[i*dim+j] = 0. + 0.*I;
    }
  }


  for (k=1; k<dim-1; k++) {

    for (i=0; i<dim; i++) {

      for (j=0; j<dim; j++) {

        temp_bas[i*dim+j] = bas[k*dim+i]*bas[j*dim+k];
        temp_t1[i*dim+j] = 2*D[k*dim+k]*g_bas[i*dim+j];
      }
    }
    AntiCommutator(dim, temp_bas, D, temp_t2);
    for (i=0; i<dim; i++) {
      for (j=0; j<dim; j++) {
        LD[i*dim+j] += beta*temp_t1[i*dim+j] - beta*temp_t2[i*dim+j];
      }
    }
 }
 for (i=0; i<dim; i++) {
   for (j=0; j<dim; j++) {
     P[i*dim+j] = LD[i*dim+j];
   }
 }


 free(temp_bas);
 free(g_bas);
 free(temp_t1);
 free(temp_t2);
 free(LD);
}


// double *bas -> the basis states for the density matrix
// double complex *P -> The contribution to the time-derivative of the Density matrix coming from dephasing
void L_Deph(int dim, double alpha, double complex *D, double *bas, double complex *P) {

  int i, j, k;
  double *temp_bas;
  double complex *temp_t1, *temp_t2, *LD;
  temp_bas = (double *)malloc(dim*dim*sizeof(double));
  temp_t1  = (double complex *)malloc(dim*dim*sizeof(double complex));
  temp_t2  = (double complex *)malloc(dim*dim*sizeof(double complex));
  LD       = (double complex *)malloc(dim*dim*sizeof(double complex));

  for (i=0; i<dim; i++) {
    for (j=0; j<dim; j++) {
      LD[i*dim+j] = 0. + 0.*I;
    }
  }

  for (k=1; k<dim-1; k++) {

    for (i=0; i<dim; i++) {

      for (j=0; j<dim; j++) {

        temp_bas[i*dim+j] = bas[k*dim+i]*bas[j*dim+k];
        temp_t1[i*dim+j] = 2*D[k*dim+k]*temp_bas[i*dim+j];
      }
    }
    AntiCommutator(dim, temp_bas, D, temp_t2);
    for (i=0; i<dim; i++) {
      for (j=0; j<dim; j++) {
        LD[i*dim+j] += alpha*temp_t1[i*dim+j] - alpha*temp_t2[i*dim+j];
      }
    }
 }

 for (i=0; i<dim; i++) {
   for (j=0; j<dim; j++) {
     P[i*dim+j] = LD[i*dim+j];
   }
 }
 
 free(temp_bas);
 free(temp_t1);
 free(temp_t2);
 free(LD);

}

