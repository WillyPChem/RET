#include <fftw3.h>
#include<math.h>
#include<stdio.h>
#include<malloc.h>
#include<complex.h>

#define REAL 0
#define IMAG 1
//void Commutator (int dim, double H[dim][dim], double D[dim][dim], double P[dim][dim]);
void RK3(int Nlevel, double time, double *bas, double *E, double *Hint, double *Mu, double *Dis, double complex *D, double dt);
void Liouville(int dim, double complex *H, double complex *D, double complex *P);
void AntiCommutator(int dim, double *H, double complex *D, double complex *P);
void PrintComplexMatrix(int dim, double complex *M);
void PrintRealMatrix(int dim, double *M);
void FormDM(int dim, double complex *C, double complex *Ct, double complex *DM);
void L_Diss(int Nlevel, double *gamma, double complex *D, double *bas, double complex *P);
void Fourier (double complex *dm, int n, double dt, double complex *ftout, double *freqvec);
double complex TrMuD(int Nlevel, double *Mu, double complex *D);
double E_Field(double time);
void FillDFTArray(int fftw_iter, double real_val, double imag_val, fftw_complex* ar);

// Function prototype for H_interaction
void H_interaction(int dim, double *Hint, double *mu, double dpm, double *R);


// NOTE!!!  You need three global variables for the rates associated with 
// The Lindblad operators - gamma, beta, and alpha should be defined here according
// to their values in the journal of physical chemistry letters paper
//int Nlevel = 3;
//int dim = Nlevel*Nlevel;

double pi = 4*atan(1.0);
double wn_to_au = 4.55633528e-6; 
double mu_au_to_si = 8.47835326e-30; // 1 a.u. = 8.47835326e-30 C m
double E_au_to_si = 5.14220652e11;  // 1 a.u. = 5.14220652e11 V/m
double omega_au = 4.134147e+16;;
int main() {

  //Nanoparticles Variables here
  int numTime = 1000000;
  int zeropad = 1000000;
  double *E, *Mu, *Dis, *bas, *Hint;
  double complex *H, *D, *P;
  double dt = 0.01;
  int Nlevel, dim;
  
  // NP levels can be variable in principle
  printf("  How many states are in your NP system? \n");
  scanf("%i",&Nlevel);
  dim = Nlevel*Nlevel;

  int dft_dim = numTime+zeropad;
  // FFTW variables here -> inputs to fft
  fftw_complex *dipole;
  fftw_complex *efield;
  fftw_complex *nps;
  fftw_complex *efs;

  // Allocate memory for FFTW arrays
  dipole = (fftw_complex*)malloc(dft_dim*sizeof(fftw_complex));
  efield = (fftw_complex*)malloc(dft_dim*sizeof(fftw_complex));
  nps = (fftw_complex*)malloc(dft_dim*sizeof(fftw_complex));
  efs = (fftw_complex*)malloc(dft_dim*sizeof(fftw_complex));

  fftw_plan npp = fftw_plan_dft_1d(dft_dim,
                                      dipole,
                                      nps,
                                      FFTW_BACKWARD,
                                      FFTW_ESTIMATE);

  fftw_plan efp = fftw_plan_dft_1d(dft_dim,
                                      efield,
                                      efs,
                                      FFTW_BACKWARD,
                                      FFTW_ESTIMATE);

  // Allocate memory for all other arrays
  // NP
  H = (double complex*)malloc(dim*sizeof(double complex));
  D = (double complex *)malloc(dim*sizeof(double complex));
  P = (double complex *)malloc(dim*sizeof(double complex));
  E  = (double *)malloc(dim*sizeof(double));
  Mu = (double *)malloc(dim*sizeof(double));
  Dis = (double *)malloc(dim*sizeof(double));
  bas = (double *)malloc(dim*sizeof(double));
  Hint = (double *)malloc(dim*sizeof(double));
  // Variables for instantaneous quantities  
  double tr;
  double complex dipole_moment;
  FILE *dfp;
  FILE *popfp;

  FILE *Efp, *Mufp, *Disfp;

  // Open each file for reading
  Efp = fopen("Matrices/SMA_PEAK1/Energy.txt","r");
  Mufp = fopen("Matrices/SMA_PEAK1/Dipole.txt","r");
  Disfp = fopen("Matrices/SMA_PEAK1/Dissipation.txt","r");

  // Density matrix element D(i,j) is accessed as D[i*Nlevel+j];
  D[0] = 1. + 0.*I;
  // NP
  for (int i=1; i<dim; i++){
    D[i] = 0. + 0.*I;
  }

  // BUILD DM BASIS - this comes into play in Lindblad operator
  // NP
  for (int i=0; i<Nlevel; i++) {
    for (int j=0; j<Nlevel; j++) {
      if (i==j){
        bas[i*Nlevel+j] = 1.0;
      }     
      else{
        bas[i*Nlevel+j] = 0.;
      }
    }
  }
  // Get parameters for NP and MG from files
  double val;
  // NP
  for (int i=0; i<dim; i++) {


       // Read from energy file and store to the energy matrix
       fscanf(Efp,"%lf",&val);
       E[i] = val;

       fscanf(Mufp,"%lf",&val);
       Mu[i] = val/1.;

       fscanf(Disfp,"%lf",&val);
       Dis[i] = val;

       Hint[i] = 0.;
  }
  // Print parameters to screen
  printf("\nE\n");
  PrintRealMatrix(Nlevel, E);
  printf("\nMu\n");
  PrintRealMatrix(Nlevel,Mu);
  printf("\nDis\n");
  PrintRealMatrix(Nlevel,Dis);
  printf("\nBas\n");
  PrintRealMatrix(Nlevel,bas);
  printf("\nDM\n");
  PrintComplexMatrix(Nlevel,D);
  // Data files for printing instantaneous data
  dfp = fopen("DATA/SMA_PEAK1/DipoleMoment.dat","w");
  popfp = fopen("DATA/SMA_PEAK1/Population.dat","w");

  // Get initial dipole moments
  dipole_moment = TrMuD(Nlevel, Mu, D)*mu_au_to_si;

  FillDFTArray(0, creal(dipole_moment), cimag(dipole_moment), dipole);
  FillDFTArray(0, 0., 0., efield);

  

  for (int i=1; i<numTime; i++) {

    // Calculate Hint now! 
    RK3(Nlevel, dt*i, bas, E, Hint, Mu, Dis, D, dt);
    
    dipole_moment = TrMuD(Nlevel, Mu, D)*mu_au_to_si; 
    FillDFTArray(i, creal(dipole_moment), cimag(dipole_moment), dipole);
    

    fprintf(popfp,"\n %f ",dt*i);
    tr=0.;
    for (int j=0; j<Nlevel; j++) {

      fprintf(popfp," %12.10e",creal(D[j*Nlevel+j]));
      tr+=creal(D[j*Nlevel+j]);

    }
    fprintf(popfp," %12.10e",tr);

    // Uncomment if you want a file with dipole moment data in it for the nanoparticle
    fprintf(dfp," %f  %12.10e  %12.10e\n",dt*i,creal(dipole_moment),cimag(dipole_moment));
 
    FillDFTArray(i,  E_au_to_si*E_Field(dt*i), 0, efield);
  }

  // Done with iterations, now zero-pad the dipole and efield arrays!
  for (int i=numTime; i<zeropad; i++) {

    FillDFTArray(i, 0., 0., dipole);
    FillDFTArray(i, 0., 0., efield);

  }

  // This actually performns the Fourier transform 
  fftw_execute(npp);
  fftw_execute(efp);
  

/*My slow DFT function! 
  Fourier(dipole, numTime+zeropad, dt, NPSpectrum, EV);
  Fourier(dipoleMG, numTime+zeropad, dt, MGSpectrum, EV);
  Fourier(efield, numTime+zeropad, dt, LaserSpectrum, EV);
*/

  
  FILE *absfp; 
  absfp = fopen("DATA/SMA_PEAK1/AbsorptionSpectrum.dat","w");
  fprintf(absfp, "#  Energy (ev)    Scattering NP      Absorption NP\n");
  
  int nfreq = 5001;
  // ~4.05 eV is max energy/ max freq
  double maxfreq = 50*0.08188379587298;
  double df = maxfreq / (nfreq - 1);
  double eps_0 = 1.0 / (4.0 * M_PI);
  double dw = 2*pi/((numTime+zeropad)*dt);
  double area = 0.;
  double sig_max = 0.;
  for (int i=1; i<(numTime+zeropad); i++) {

    // This is omega in atomic units - same as energy in atomic units
    double omega = 2*pi*i/((numTime+zeropad)*dt);
    if (omega * 27.211>maxfreq) break;

    double omega_si = omega*omega_au;
    double eev = omega*27.211;
    double k = omega_si/2.99792458e+8;
    double pre_scat = k*k*k*k/(6*pi*8.854187e-12*8.854187e-12); 
    double pre_abs = k/(8.854187e-12);
    double complex npr = nps[i][0]/numTime;
    double complex npi = nps[i][1]/numTime;

    double complex efr = efs[i][0]/numTime;
    double complex efi = efs[i][1]/numTime;

    //double pre = k/(8.854187e-12);
    double complex alphaNP = (npr+I*npi)/(efr+I*efi);

    double sig_scat_NP = pre_scat * creal(alphaNP*conj(alphaNP));
    double sig_abs_NP = pre_abs * cimag(alphaNP);

    if (sig_abs_NP>sig_max) { sig_max = sig_abs_NP; }
    area += sig_abs_NP*dw;
    // Going to print absorption cross section in abs/micron^2
    fprintf(absfp, "  %12.10e  %12.10e  %12.10e\n",eev,sig_scat_NP, sig_abs_NP);
  }
  
  printf("  AREA IS     %12.10e\n", area);
  printf("  SIG_MAX is  %12.10e\n", sig_max);
  fclose(absfp);
  
  //fclose(dfp);
  return 0;
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
void RK3(int Nlevel, double time, double *bas, double *E, double *Hint, double *Mu, double *Dis, double complex *D, double dt) {

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

  for (i=0; i<Nlevel; i++) {
    gamma[i] = Dis[i*Nlevel+i];
  }
  Efield = E_Field(time);

  // Compute full Hamiltonian at current time t 
  for (i=0; i<dim; i++) {


      H[i] = E[i] + Hint[i] - Efield*Mu[i]; // - I*Dis[i];

  } 

  //PrintComplexMatrix(Nlevel, H);

  // Get dPsi(n)/dt at initial time!
  // Two main changes needed to couple the molecule and nanoparticle:
  // (1) Liouville function needs to include H_interaction
  // (2) We need to use Liouville/L_Diss to update both the molecule and the nanoparticle density matrix
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

      H[i] = E[i] + Hint[i] - Efield*Mu[i]; // - I*Dis[i];

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
  double tau = 50.;

  //Ef = 0.01*sin(pi*time/tau)*sin(pi*time/tau)*exp(-0.005*time)*(sin(0.07423*time)+sin(0.1*time)+sin(0.5*time));
  if (time<tau) {

    Ef = 0.0001*sin(time*pi/tau)*sin(time*pi/tau)*sin(0.07423*time);

  }
  else Ef = 0.;
  //Ef = 0.000001*sin(0.05*time);  
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
  double gk;
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

    gk = gamma[k];
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
        LD[i*Nlevel+j] += gk*temp_t1[i*Nlevel+j] - gk*temp_t2[i*Nlevel+j];
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

void Fourier(double complex *dm, int n, double dt, double complex *ftout, double *freqvec){
  //FILE *fp;
  //fp = fopen("Absorption_SpectrumAu.txt","w");
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
      sumreal += creal(cexp(-I*angle)*dm[t])*dt;
      sumimag += cimag(cexp(-I*angle)*dm[t])*dt;
      //sumreal += creal(dm[t]) * cos(angle) + cimag(dm[t]) * sin(angle);
      //sumimag  += creal(dm[t]) * sin(angle) + cimag(dm[t]) * cos(angle);
    }
    ftout[k] = sumreal + sumimag*I;
    // Energy in eV
    freqvec[k] = w*27.2114;
    //fprintf(fp," %12.10e  %12.10e\n",w*(27.2114),sumreal*sumreal+sumimag*sumimag);
  }
  //fclose(fp);
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


void H_interaction(int dim, double *Hint, double *mu, double dpm, double *R) {
  
  int i; 
 // double *tmp1, *tmp2;
  double oer2, oer3;
  double scal_R = sqrt(R[0]*R[0] + R[1]*R[1] + R[2]*R[2]);
 
  oer2 = pow(scal_R,-2.);
  oer3 = pow(scal_R,-3.);
 
  // Write code between here!
 
 for (i=0; i<dim*dim; i++){

   // Very important!  Currently assuming z-polarized light, so only the <mu>_z and mu_z terms
   //                  are non-zero, hence we consider only R[2]*mu and <mu>*R[2] 
   Hint[i] = oer3*(dpm*mu[i] -3*R[2]*mu[i]*R[2]*dpm*oer2);
 } 

}

void FillDFTArray(int fftw_iter, double real_val, double imag_val, fftw_complex* ar) {

  ar[fftw_iter][REAL] = real_val;
  ar[fftw_iter][IMAG] = imag_val;

}
