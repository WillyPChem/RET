#include<math.h>
#include<stdio.h>
#include<stdlib.h>
#include<malloc.h>
#include<complex.h>

//void Commutator (int dim, double H[dim][dim], double D[dim][dim], double P[dim][dim]);
void RK3(int Nlevel, double time, double *bas, double *E, double *Mu, double *Dis, double complex *D, double dt);
void Liouville(int dim, double complex *H, double complex *D, double complex *P);
void AntiCommutator(int dim, double *H, double complex *D, double complex *P);
void PrintComplexMatrix(int dim, double complex *M);
void PrintRealMatrix(int dim, double *M);
void FormDM(int dim, double complex *C, double complex *Ct, double complex *DM);
void L_Diss(int Nlevel, double *gamma, double complex *D, double *bas, double complex *P);
void Fourier (double complex *dm, int n, double dt);
void TrMuMol (int dim, double complex *H, double complex *D, double *Mu);
double E_Field(double time);
double complex TrMuD(int Nlevel, double *Mu, double complex *D);

// NOTE!!!  You need three global variables for the rates associated with 
// The Lindblad operators - gamma, beta, and alpha should be defined here according
// to their values in the journal of physical chemistry letters paper
//int Nlevel = 3;
//int dim = Nlevel*Nlevel;

double pi = 4*atan(1.0);
double wn_to_au = 4.55633528e-6; 

int main() {
//Nanoparticles Variables here
  int Nlevel = 3;
  int numTime = 400000;
  int zeropad = 100000;
  int dim = Nlevel*Nlevel;  
  double *E, *Mu, *Dis, *bas;
  double complex *H, *D, *P, *dipole;
  double dt = 0.01;

  dipole = (double complex *)malloc((numTime+zeropad)*sizeof(double complex));
  H = (double complex*)malloc(dim*sizeof(double complex));
  D = (double complex *)malloc(dim*sizeof(double complex));
  P = (double complex *)malloc(dim*sizeof(double complex));
  E  = (double *)malloc(dim*sizeof(double));
  Mu = (double *)malloc(dim*sizeof(double));
  Dis = (double *)malloc(dim*sizeof(double));
  bas = (double *)malloc(dim*sizeof(double));

//Molecule Variables here
 // int Nlevel = 3;
 // int dim = Nlevel*Nlevel;
  double *EMG, *MuMG, *DisMG, *basMG;
  double complex *HMG, *DMG, *PMG, *dipoleMG;
 // double dt = 0.01;

  dipoleMG = (double complex *)malloc((numTime+zeropad)*sizeof(double complex));
  HMG = (double complex*)malloc(dim*sizeof(double complex));
  DMG = (double complex *)malloc(dim*sizeof(double complex));
  PMG = (double complex *)malloc(dim*sizeof(double complex));
  EMG  = (double *)malloc(dim*sizeof(double));
  MuMG = (double *)malloc(dim*sizeof(double));
  DisMG = (double *)malloc(dim*sizeof(double));
  basMG = (double *)malloc(dim*sizeof(double));

  
  // Initialize Denistry matrix as a superposition state of energy eigenstate 1 and energy eigenstate 2
  // Density matrix element D(i,j) is accessed as D[i*Nlevel+j];
  D[0] = 1. + 0.*I;
  for (int i=1; i<dim; i++) D[i] = 0. + 0.*I;


  // BUILD DM BASIS
  for (int i=0; i<Nlevel; i++) {
    for (int j=0; j<Nlevel; j++) {

      if (i==j) bas[i*Nlevel+j] = 1.0;

      else bas[i*Nlevel+j] = 0.;

    }
  }

  // Variable that is a file pointer for each data file:
  FILE *Efp, *mufp, *disfp;

  // Open each file for reading
  Efp = fopen("EnergyAu.txt","r");
  mufp = fopen("DipoleAu.txt","r");
  disfp = fopen("DissipationAu.txt","r");  
  
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
  double complex dipole_moment;
  FILE *dfp;
  dfp = fopen("DipoleMoment.txt","w");

  dipole[0] = TrMuD(Nlevel, Mu, D);
  for (int i=1; i<numTime; i++) {

    RK3(Nlevel, dt*i, bas, E, Mu, Dis, D, dt);

    printf("\n %f ",dt*i);
    tr=0.;
    for (int j=0; j<Nlevel; j++) {

      printf(" %12.10e",creal(D[j*Nlevel+j]));
      tr+=creal(D[j*Nlevel+j]);
    }
    printf("\n #Trace is %12.10e",tr);
    dipole_moment = TrMuD(Nlevel, Mu, D);
    dipole[i] = dipole_moment;
    fprintf(dfp," %f  %12.10e  %12.10e\n",dt*i,creal(dipole_moment),cimag(dipole_moment));
  }

  for (int i=numTime; i<zeropad; i++) {

    dipole[i] = 0. + 0.*I;

  }

  Fourier(dipole, numTime+zeropad, dt);

  fclose(dfp);
  return 0;
}

 FILE *EfpMG, *mufpMG, *disfpMG;

  // Open each file for reading
  EfpMG = fopen("Energy.txt","r");
  mufpMG = fopen("Dipole.txt","r");
  disfpMG = fopen("Dissipation.txt","r");

 double val;
  for (int i=0; i<dim; i++) {


       // Read from energy file and store to the energy matrix
       fscanf(Efp,"%lf",&val);
       EMG[i] = val;

       fscanf(mufp,"%lf",&val);
       MuMG[i] = val;

       fscanf(disfp,"%lf",&val);
       DisMG[i] = val;

  }

  printf("\nE\n");
  PrintRealMatrix(Nlevel, EMG);
  printf("\nMu\n");
  PrintRealMatrix(Nlevel,MuMG);
  printf("\nDiss\n");
  PrintRealMatrix(Nlevel,DisMG);

  double tr=0.;
  double complex dipole_momentMG;
  FILE *dfp;
  dfp = fopen("DipoleMomentMG.txt","w");

  dipoleMG[0] = TrMuD(Nlevel, MuMG, D);
  for (int i=1; i<numTime; i++) {

    RK3(Nlevel, dt*i, basMG, EMG, MuMG, DisMG, DMG, dt);

    printf("\n %f ",dt*i);
    tr=0.;
    for (int j=0; j<Nlevel; j++) {

      printf(" %12.10e",creal(D[j*Nlevel+j]));
      tr+=creal(D[j*Nlevel+j]);
    }
    printf("\n #Trace is %12.10e",tr);
    dipole_momentMG = TrMuDMG(Nlevel, MuMG, DMG);
    dipoleMG[i] = dipole_momentMG;
    fprintf(dfp," %f  %12.10e  %12.10e\n",dt*i,creal(dipole_momentMG),cimag(dipole_momentMG));
  }
  for (int i=numTime; i<zeropad; i++) {

    dipoleMG[i] = 0. + 0.*I;

  }

  Fourier(dipole, numTime+zeropad, dt);

  fclose(dfp);
  return 0;
 

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
void RK3(int Nlevel, double time, double *bas, double *E, double *Mu, double *Dis, double complex *D, double dt) {

  int i, j;
  double complex *D_dot, *D2, *D3, *D_np1, *k1, *k2, *k3;
  double complex *H, *LD;  // Contribution to Ddot from Lindblad dissipation
  double *gamma;
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
  LD    = (double complex *)malloc(dim*sizeof(double complex));
  gamma = (double *)malloc(Nlevel*sizeof(double));
  
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

  gamma[0] = 0.;
  gamma[1] = Dis[1*Nlevel+1];
  gamma[2] = Dis[2*Nlevel+2];

  Efield = E_Field(time);

  // Compute full Hamiltonian at current time t 
  for (i=0; i<dim; i++) {


      H[i] = E[i] + Efield*Mu[i]; // - I*Dis[i];

  } 

  //PrintComplexMatrix(Nlevel, H);

  // Get dPsi(n)/dt at initial time!
  Liouville(Nlevel, H, D, D_dot);
  L_Diss(Nlevel, gamma, D, bas, LD);
  //PrintComplexMatrix(Nlevel, D);
  //PrintComplexMatrix(Nlevel, D_dot);


  // Compute approximate wfn update with Euler step
  for (i=0; i<dim; i++) {
    k1[i] = dt*(D_dot[i]+LD[i]);
    D2[i] = D[i] + k1[i]/2.;
  }

  // Update Field!
  Efield = E_Field(time+dt/2.);

  // Compute full Hamiltonian at partially updated time t 
  for (i=0; i<dim; i++) {

      H[i] = E[i] + Efield*Mu[i]; // - I*Dis[i];

  }

  //PrintComplexMatrix(Nlevel, H);
  // Get dPsi(n+k1/2)/dt
  Liouville (Nlevel, H, D2, D_dot);
  L_Diss(Nlevel, gamma, D2, bas, LD);
  
  // Compute approximate wfn update with Euler step
  for (i=0; i<dim; i++) {
    k2[i] = dt*(D_dot[i] + LD[i]);
    D3[i] = D[i] + k2[i]/2.;
  }

  // Get dPsi(n+k2/2)/dt
  Liouville (Nlevel, H, D3, D_dot);
  L_Diss(Nlevel, gamma, D3, bas, LD);

  // Compute approximate update with Euler step
  for (i=0; i<dim; i++) {
    k3[i] = dt*(D_dot[i] + LD[i]);
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
  free(gamma);
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
  double tau = 75.;
  if (time<tau) {

    Ef = 0.01*sin(time*pi/tau)*sin(time*pi/tau)*sin(0.07423*time);

  }
  else Ef = 0.;
    
  return Ef;


}

void L_Diss(int Nlevel, double *gamma, double complex *D, double *bas, double complex *P) {

  int i, j, k;
  double *temp_bas, *g_bas;
  double complex *temp_t1, *temp_t2, *LD;
  temp_bas = (double *)malloc(Nlevel*Nlevel*sizeof(double));
  g_bas    = (double *)malloc(Nlevel*Nlevel*sizeof(double));
  temp_t1  = (double complex *)malloc(Nlevel*Nlevel*sizeof(double complex));
  temp_t2  = (double complex *)malloc(Nlevel*Nlevel*sizeof(double complex));
  LD       = (double complex *)malloc(Nlevel*Nlevel*sizeof(double complex));

  // Form |g><g| matrix
  for (i=0; i<Nlevel; i++) {
    for (j=0; j<Nlevel; j++) {
      g_bas[i*Nlevel+j] = bas[0*Nlevel+i]*bas[j*Nlevel+0];
      LD[i*Nlevel+j] = 0. + 0.*I;
    }
  }
  //printf("  |g><g|  \n");
  //PrintRealMatrix(Nlevel, g_bas);

  for (k=1; k<Nlevel; k++) {

    for (i=0; i<Nlevel; i++) {

      for (j=0; j<Nlevel; j++) {

        temp_bas[i*Nlevel+j] = bas[k*Nlevel+i]*bas[j*Nlevel+k];
        temp_t1[i*Nlevel+j] = 2*D[k*Nlevel+k]*g_bas[i*Nlevel+j];
      }
    }
   // printf("   |%i><%i| \n",k,k);
   // PrintRealMatrix(Nlevel, temp_bas);

    AntiCommutator(Nlevel, temp_bas, D, temp_t2);
    for (i=0; i<Nlevel; i++) {
      for (j=0; j<Nlevel; j++) {
        LD[i*Nlevel+j] += gamma[2]*temp_t1[i*Nlevel+j] - gamma[2]*temp_t2[i*Nlevel+j];
      }
    }
 }
 for (i=0; i<Nlevel; i++) {
   for (j=0; j<Nlevel; j++) {
     P[i*Nlevel+j] = LD[i*Nlevel+j];
   }
 }

 //PrintComplexMatrix(Nlevel, P);

 free(temp_bas);
 free(g_bas);
 free(temp_t1);
 free(temp_t2);
 free(LD);
}

void Fourier (double complex *dm, int n, double dt){
  FILE *fp;
  fp = fopen("Absorption_SpectrumAu.txt","w");
  double wmin=0.5*0.07;
  double wmax=2*0.07;
  int maxk = 500;
  double dw = (wmax-wmin)/maxk;
  double time;
 
  for (int k = 0; k <=maxk; k++) {
    double sumreal = 0;
    double sumimag = 0;
    double  w = wmin+k*dw;

    for (int t = 0; t < n; t++){
      time = dt*t;
      double angle = time*w;
      sumreal += creal(dm[t]) * cos(angle) + 0. * sin(angle);
      sumimag  += -creal(dm[t]) * sin(angle) + 0. * cos(angle);
    }
    fprintf(fp," %12.10e  %12.10e\n",w*(27.2114),sumreal*sumreal+sumimag*sumimag);
  }
  fclose(fp);
}

double complex TrMuD(int Nlevel, double *Mu, double complex *D) {

  double complex tr = 0. + 0.*I;
  for (int i=0; i<Nlevel; i++) {

    double complex sum = 0. + 0.*I;
    for (int k=0; k<Nlevel; k++) {

      sum += Mu[i*Nlevel+k]*D[k*Nlevel+i];

    }
    tr += sum;

  } 

  return tr;

}

double complex Hint(int Nlevel, double complex *D, double *MU, double complex *DMG, double *MuMG, double R);

double R = 1.;
double complex TrMumol = 0. + 0.*I;
double complex TrMuNp = 0. + 0.*I;


int i,j,k;
 for (int i=0; i<Nlevel; i++){

     for (int j=0, j<Nlevel; j++){

      for (int k=0, k<Nlevel; k++){

  sum += (1/R^3)*(MuNp*MuMol- ((MuNp*R)(R*MuMol)/R^2))

}
  Hint += sum;
}

return Hint;

}
