#include<math.h>
#include<stdio.h>
#include<stdlib.h>
#include<malloc.h>
#include<complex.h>

//void Commutator (int dim, double H[dim][dim], double D[dim][dim], double P[dim][dim]);
void RK3(int dim, double *bas, double *xvec, double complex *wfn, double dx, double dt);
void Liouville(int dim, double *H, double complex *D, double complex *P);
void AntiCommutator(int dim, double *H, double complex *D, double complex *P);
void PrintComplexMatrix(int dim, double complex *M);
void PrintRealMatrix(int dim, double *M);
void FormDM(int dim, double complex *C, double complex *Ct, double complex *DM);
void L_Deph(int dim, double alpha, double complex *D, double *bas, double complex *P);
void L_Diss(int dim, int g, double beta, double complex *D, double *bas, double complex *P);
void L_Sink(int dim, int chan, int s, double gamma, double complex *D, double *bas, double complex *P);

double gam = 1.21e-8*0;
double bet = 7.26e-5;
double alph = 1.52e-4;
double pi = 3.14159625;

// NOTE!!!  You need three global variables for the rates associated with 
// The Lindblad operators - gamma, beta, and alpha should be defined here according
// to their values in the journal of physical chemistry letters paper
int dim = 9;
int g = 0;
int s = 8;
int chan = 3;

double pi = 4*atan(1.0);
double wn_to_au = 4.55633528e-6; 
// Beta for L_Diss from Mazziotti paper
// alpha, beta, gamma - 1.52e-4, 7.26e-5, 1.21e-8 
//double be = 7.26e-5;

double be = 1.21e-8;
double gam = 1.52e-4;
//double gam = 1.21e-8;
double alph = 7.26e-5;
int main() {

  int dim = 9;
  double *H, *bas;
  double complex *D, *P, *C, *Ct, *Dp;
  double dt = 0.1;
  double dx = 1.;

  H = (double *)malloc(dim*dim*sizeof(double));
  D = (double complex *)malloc(dim*dim*sizeof(double complex));
  P = (double complex *)malloc(dim*dim*sizeof(double complex));
  Dp= (double complex *)malloc(dim*dim*sizeof(double complex));
  C = (double complex *)malloc(dim*sizeof(double complex));
  Ct= (double complex *)malloc(dim*sizeof(double complex));
  bas = (double *)malloc(dim*dim*sizeof(double));

  // Element DM(0,0) = D(0*dim+0)
  // Element DM(1,0) = D(1*dim+0)
  // Element DM(0,1) = D(0*dim+1)
  // In general, DM(i,j) = D(i*dim+j) where DM is the matrix
  // form of D and D is a vector of length dim*dim


  //  Build basis vectors |k> for excitonic states and ground state...
  //  e.g. |g> = (1, 0, 0, 0, 0, 0, 0, 0,0)
  //  and  |1> = (0, 1, 0, 0, 0, 0, 0, 0,0)

  for (int i=0; i<dim; i++) {
    for (int j=0; j<dim; j++) {
       bas[i*dim+j] = 0.;
       D[i*dim+j] = 0. * I*0.;
    }
  }
  for (int i=0; i<dim; i++) {
    bas[i*dim+i] = 1.;
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
=======


  //  The density matrix... initialize in some excitonic state, i.e. D[1][1] = 1.0 
  //  initializes in the first excitonic state
  D[1*dim+1] = 1. + 0.*I;
  // Below is a superposition btw ground and first excitonic state
  /*D[0*dim + 0] = 1/2. + 0.*I;
  D[0*dim + 1] = 1/2. + 0.*I;
  D[1*dim + 0] = 1/2. + 0.*I;
  D[1*dim + 1] = 1/2. + 0.*I;
  */

  // Build Hamiltonian for Liouville equation in atomic units - note that ground state will have energy 0... coupling 
  // between excited-states and ground state will occur through dissipation term in Liouville equation, gs -> excited
  // states will have no coupling ... 
  FILE *Hfp;
  Hfp = fopen("Ham.txt","r");
  char *crap;
  crap = (char *)malloc(1000*sizeof(char));
  double val;

  for (int i=1; i<dim-1; i++) {
    for (int j=1; j<dim-1; j++) {
      fscanf(Hfp,"%s",crap);
      val = atof(crap);

      H[i*dim+j] = val*wn_to_au;
      //H[i*dim+j] = val;
    }
  }
  fclose(Hfp);
  PrintRealMatrix(dim, H);



  /* Bact. Chlorophor Hamiltonian from 
    "Efficient energy transfer in light-harvesting systems, I: 
    optimal temperature, reorganization
    energy and spatial–temporal correlations"
    2010 New J. Phys. 12 105012
    Jianlan Wu, Fan Liu, Young Shen, Jianshu Cao, and Robert J Silbey
 
  H = [280, −106, 8, −5, 6, −8, −4,
 −106, 420, 28, 6, 2, 13, 1,
  8, 28, 0, −62, −1, −9, 17,
 −5, 6, −62, 175, −70, −19, −57,
  6, 2, −1, −70, 320, 40, −2,
 −8, 13, −9, −19, 40, 360, 32,
 −4, 1, 17, −57, −2, 32, 260]
 */

  double tr=0.;
  for (int i=1; i<500000; i++) {

    RK3(dim, bas, H, D, dx, dt);

    printf("\n %f ",dt*i*0.02418);
    tr=0.;
    for (int j=0; j<dim; j++) {

      printf(" %12.10e",creal(D[j*dim+j]));
      tr+=creal(D[j*dim+j]);
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
void RK3(int dim, double *bas, double *H, double complex *wfn, double dx, double dt) {

  int i, j;
  double complex *wfn_dot, *wfn2, *wfn3, *wfn_np1, *k1, *k2, *k3;
  double complex *LDiss, *LSink, *LDeph;  // Contribution to Ddot from Lindblad dissipation
 
  wfn_dot = (double complex *)malloc((dim*dim)*sizeof(double complex));
  LDiss = (double complex *)malloc((dim*dim)*sizeof(double complex));
  LSink = (double complex *)malloc((dim*dim)*sizeof(double complex));
  LDeph = (double complex *)malloc((dim*dim)*sizeof(double complex));
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
    LDiss[i*dim+j] = 0. + 0.*I;
    LSink[i*dim+j] = 0. + 0.*I;
    LDeph[i*dim+j] = 0. + 0.*I;
    }
  }

  // Get dPsi(n)/dt at initial time!
  Liouville (dim, H, wfn, wfn_dot);
  L_Diss(dim, 0, be, wfn, bas, LDiss);
  L_Deph(dim, alph, wfn, bas, LDeph);
  L_Sink(dim, 3, 8, gam, wfn, bas, LSink); 

  // Compute approximate wfn update with Euler step
  for (i=0; i<dim; i++) {
    for (j=0; j<dim; j++) {
    k1[i*dim+j] = dt*(wfn_dot[i*dim+j] + LDiss[i*dim+j] + LDeph[i*dim+j] + LSink[i*dim+j] );
    wfn2[i*dim+j] = wfn[i*dim+j] + k1[i*dim+j]/2.;
   }

  }
  // Get dPsi(n+k1/2)/dt
  Liouville (dim, H, wfn2, wfn_dot);
  L_Diss(dim, 0, be, wfn2, bas, LDiss);
  L_Deph(dim, alph, wfn2, bas, LDeph);
  L_Sink(dim, 3, 8, gam, wfn2, bas, LSink);
  // Compute approximate wfn update with Euler step
  for (i=0; i<dim; i++) {
    for (j=0; j<dim; j++) {
    k2[i*dim+j] = dt*(wfn_dot[i*dim+j] + LDiss[i*dim+j] + LDeph[i*dim+j] + LSink[i*dim+j] );
    wfn3[i*dim+j] = wfn[i*dim+j] + k2[i*dim+j]/2.;
  }
  }
  // Get dPsi(n+k2/2)/dt
  Liouville (dim,H, wfn3, wfn_dot);
  L_Diss(dim, 0, be, wfn3, bas, LDiss);
  L_Deph(dim, alph, wfn3, bas, LDeph);
  L_Sink(dim, 3, 8, gam, wfn3, bas, LSink);
  // Compute approximate update with Euler step
  for (i=0; i<dim; i++) {
    for (j=0; j<dim; j++) {
    k3[i*dim+j] = dt*(wfn_dot[i*dim+j] + LDiss[i*dim+j] + LDeph[i*dim+j] + LSink[i*dim+j] );
    wfn_np1[i*dim+j] = wfn[i*dim+j] + k1[i*dim+j]/6. + 2.*k2[i*dim+j]/3. + k3[i*dim+j]/6.;
    //wfn_np1[i*dim+j] = wfn[i*dim+j] + k1[i*dim+j];
    wfn[i*dim+j] = wfn_np1[i*dim+j];
  }
}
  free(wfn_dot);
  free(LDiss);
  free(LDeph);
  free(LSink);
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



void Liouville(int dim, double *H, double complex *D, double complex *P) {

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
 //PrintComplexMatrix(dim, LD);       
 //void PrintRealMatrix(int dim, double *M);
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
      LD[i*dim+j] = 0. + 0.*I;
      temp_bas[i*dim+j] = 0.;
      temp_t1[i*dim+j] = 0. + 0.*I;
      temp_t2[i*dim+j] = 0. + 0.*I;
      chan_bas[i*dim+j] = bas[chan*dim+i]*bas[j*dim+chan];
      s_bas[i*dim+j] = bas[s*dim+i]*bas[j*dim+s];
    }
  }

  


  for (i=0; i<dim; i++) {

    for (j=0; j<dim; j++) {

      temp_t1[i*dim+j] = D[chan*dim+chan]*s_bas[i*dim+j];
    }
  }
  AntiCommutator(dim, chan_bas, D, temp_t2);

  for (i=0; i<dim; i++) {
    for (j=0; j<dim; j++) {
      P[i*dim+j] = 2*gamma*temp_t1[i*dim+j] - gamma*temp_t2[i*dim+j];
    }
  }

 free(chan_bas);
 free(s_bas);
 free(temp_t1);
 free(temp_t2);
 free(LD);

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


