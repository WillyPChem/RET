#include<math.h>
#include<stdio.h>
#include<stdlib.h>
#include<malloc.h>
#include<complex.h>
#include<fftw3.h>


//void Commutator (int dim, double H[dim][dim], double D[dim][dim], double P[dim][dim]);
void RK3(int Nlevel, double time, double *bas, double *E, double *Hint, double *Mu, double *Dis, double complex *D, double dt);
void Liouville(int dim, double complex *H, double complex *D, double complex *P);
void AntiCommutator(int dim, double *H, double complex *D, double complex *P);
void PrintComplexMatrix(int dim, double complex *M);
void PrintRealMatrix(int dim, double *M);
void FormDM(int dim, double complex *C, double complex *Ct, double complex *DM);
void L_Diss(int Nlevel, double *gamma, double complex *D, double *bas, double complex *P);
void Fourier (double complex *dm, int n, double dt, double complex *ftout, double *freqvec);
void TrMuMol (int dim, double complex *H, double complex *D, double *Mu);
double E_Field(double time);
double complex TrMuD(int Nlevel, double *Mu, double complex *D);

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

  // MG variables here!
  int NlevelMG = 3;
  int dimMG = NlevelMG*NlevelMG;
  double *EMG, *MuMG, *MuZERO, *DisMG, *basMG, *HintMG;
  double complex *HMG, *DMG, *PMG;


  // FFTW variables here!
  fftw_complex *dipole;
  fftw_complex *efield;
  fftw_complex *dipoleMG;

  // Allocate memory for FFTW arrays
  dipole = (fftw_complex*)malloc((numTime+zeropad)*sizeof(fftw_complex));
  efield = (fftw_complex*)malloc((numTime+zeropad)*sizeof(fftw_complex));
  dipoleMG = (fftw_complex*)malloc((numTime+zeropad)*sizeof(fftw_complex));

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
  // MG
  HMG = (double complex*)malloc(dimMG*sizeof(double complex));
  DMG = (double complex *)malloc(dimMG*sizeof(double complex));
  PMG = (double complex *)malloc(dimMG*sizeof(double complex));
  EMG  = (double *)malloc(dimMG*sizeof(double));
  MuMG = (double *)malloc(dimMG*sizeof(double));
  MuZERO = (double *)malloc(dimMG*sizeof(double));
  DisMG = (double *)malloc(dimMG*sizeof(double));
  basMG = (double *)malloc(dimMG*sizeof(double));
  HintMG = (double *)malloc(dimMG*sizeof(double));

  // Variables for instantaneous quantities  
  double tr, trMG;
  double complex dipole_moment, dipole_momentMG;
  FILE *dfp, *dfpMG;
  FILE *popfp, *popMGfp;

  // Separation vector
  double *r;
  r = (double *)malloc(3*sizeof(double));
  r[0] = 0.;
  r[1] = 0.;
  r[2] = 2000000.;


  // Files for stuff:
  // Variable that is a file pointer for each data file:
  /*char *Efn, *Mufn, *Disfn, *EMGfn, *MuMGfn, *DisMGfn;
  Efn = (char *)malloc(1000*sizeof(char));
  Mufn = (char *)malloc(1000*sizeof(char));
  Disfn = (char *)malloc(1000*sizeof(char));
  EMGfn = (char *)malloc(1000*sizeof(char));
  MuMGfn = (char *)malloc(1000*sizeof(char));
  DisMGfn = (char *)malloc(1000*sizeof(char));
  */
  FILE *Efp, *Mufp, *Disfp, *EfpMG, *MufpMG, *DisfpMG;

  // Open each file for reading
  Efp = fopen("Matrices/PRB_JPCC/EnergyAu.txt","r");
  Mufp = fopen("Matrices/PRB_JPCC/DipoleAu.txt","r");
  Disfp = fopen("Matrices/PRB_JPCC/DissipationAu.txt","r");

  EfpMG = fopen("Matrices/PRB_JPCC/Energy.txt","r");
  MufpMG = fopen("Matrices/PRB_JPCC/Dipole.txt","r");
  DisfpMG = fopen("Matrices/PRB_JPCC/Dissipation.txt","r");


  // Density matrix element D(i,j) is accessed as D[i*Nlevel+j];
  D[0] = 1. + 0.*I;
  DMG[0] = 1. + 0.*I;
  // NP
  for (int i=1; i<dim; i++){
    D[i] = 0. + 0.*I;
  }
  // MG
  for (int i=1; i<dimMG; i++){
    DMG[i] = 0. + 0.*I;
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
  // MG
  for (int i=0; i<NlevelMG; i++) {
    for (int j=0; j<NlevelMG; j++) {
      if (i==j){
        basMG[i*Nlevel+j] = 1.0;
      }
      else{
        basMG[i*Nlevel+j] = 0.;
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
       Mu[i] = val;

       fscanf(Disfp,"%lf",&val);
       Dis[i] = val;
  }
  // MG
  for (int i=0; i<dimMG; i++) {

       fscanf(EfpMG,"%lf",&val);
       EMG[i] = val;

       fscanf(MufpMG,"%lf",&val);
       MuMG[i] = val;
       MuZERO[i] = 0.;
       fscanf(DisfpMG,"%lf",&val);
       DisMG[i] = val;


  }

  // Print parameters to screen
  printf("\nE\n");
  PrintRealMatrix(Nlevel, E);
  printf("\nMu\n");
  PrintRealMatrix(Nlevel,Mu);
  printf("\nDis\n");
  PrintRealMatrix(Nlevel,Dis);

  printf("\nEMG\n");
  PrintRealMatrix(NlevelMG, EMG);
  printf("\nMuMG\n");
  PrintRealMatrix(NlevelMG, MuMG);
  printf("\nDisMG\n");
  PrintRealMatrix(NlevelMG, DisMG);

  // Data files for printing instantaneous data
  dfp = fopen("DATA/PRB_JPCC/DipoleMoment.txt","w");
  dfpMG = fopen("DATA/PRB_JPCC/DipoleMomentMG.txt", "w");
  popfp = fopen("DATA/PRB_JPCC/Population.txt","w");
  popMGfp = fopen("DATA/PRB_JPCC/PopulationMG.txt","w");

  // Get initial dipole moments
  dipole_moment = TrMuD(Nlevel, Mu, D);
  dipole_momentMG = TrMuD(NlevelMG, MuMG, DMG);

  // Store dm to vectors along with initial field
  dipole[0]= dipole_moment; 

  dipoleMG[0] = dipole_momentMG;

  efield[0] = 0. + 0.*I;


  //void H_interaction(int dim, double *Hint, double *mu, double dpm, double R) 
  H_interaction(Nlevel, Hint, Mu, creal(dipole_momentMG), r); 
  H_interaction(NlevelMG, HintMG, MuMG, creal(dipole_moment), r);

  

  for (int i=1; i<numTime; i++) {

    // Calculate Hint now!
    
    RK3(Nlevel, dt*i, bas, E, Hint, Mu, Dis, D, dt);
    
    dipole_moment = TrMuD(Nlevel, Mu, D)*mu_au_to_si; 
    
    H_interaction(NlevelMG, HintMG, MuMG, creal(dipole_moment), r);
    
    RK3(NlevelMG, dt*i, basMG, EMG, HintMG, MuMG, DisMG, DMG, dt);
    //RK3(Nlevel, dt*i, basMG, EMG, HintMG, MuZERO, DisMG, DMG, dt);
    
    dipole_momentMG = TrMuD(NlevelMG, MuMG, DMG)*mu_au_to_si;
    
    H_interaction(Nlevel, Hint, Mu, creal(dipole_momentMG), r); 
   

    fprintf(popfp,"\n %f ",dt*i);
    fprintf(popMGfp,"\n %f ",dt*i);
    tr=0.;
    trMG = 0.;
    for (int j=0; j<Nlevel; j++) {

      fprintf(popfp," %12.10e",creal(D[j*Nlevel+j]));
      tr+=creal(D[j*Nlevel+j]);

    }
    for (int j=0; j<NlevelMG; j++) {


      fprintf(popMGfp,"  %12.10e",creal(DMG[j*Nlevel+j]));
      trMG+=creal(DMG[j*Nlevel+j]);

    }
    fprintf(popfp," %12.10e",tr);
    fprintf(popMGfp," %12.10e",trMG);

    dipole[i] = dipole_moment;
    //dipole[i][1] = 0.;

    // Uncomment if you want a file with dipole moment data in it for the nanoparticle
    fprintf(dfp," %f  %12.10e  %12.10e\n",dt*i,creal(dipole_moment),cimag(dipole_moment));
  
    dipoleMG[i] = dipole_momentMG;
    //dipoleMG[i][1] = 0.;
    // Uncomment if you want a file with dipole moment data in it for molecule
    fprintf(dfpMG,"%f %12.10e %12.10e\n",dt*i,creal(dipole_momentMG),cimag(dipole_momentMG));

    efield[i] = E_au_to_si*E_Field(dt*i) + 0.*I;
    //efield[i][1] = 0.;
  }

  for (int i=numTime; i<zeropad; i++) {

    dipole[i] = 0.+0.*I;
    //dipole[i][1] = 0.;

    dipoleMG[i] = 0.+0.*I;
    //dipoleMG[i][1] = 0.;

    efield[i] = 0.+0.*I;
    //efield[i][1] = 0.;
  }
 
   
  // double complex *ftout, double *freqvec
  double complex *NPSpectrum, *MGSpectrum, *LaserSpectrum;
  double *EV;
  // Currently FT over 500 different energies, hence the arrays for the various spectra
  // need to have length 500... will make them longer for good measure
  NPSpectrum = (double complex *)malloc(1000*sizeof(double complex));
  MGSpectrum = (double complex *)malloc(1000*sizeof(double complex));
  LaserSpectrum = (double complex *)malloc(1000*sizeof(double complex));
  EV = (double *)malloc(1000*sizeof(double));

  // Now calculate the spectra
  fftw_plan np, mg, ef;
  np = fftw_plan_dft_1d((numTime+zeropad), dipole, dipole, FFTW_BACKWARD, FFTW_ESTIMATE);
  mg = fftw_plan_dft_1d((numTime+zeropad), dipoleMG, dipoleMG, FFTW_BACKWARD, FFTW_ESTIMATE);
  ef = fftw_plan_dft_1d((numTime+zeropad), efield, efield, FFTW_BACKWARD, FFTW_ESTIMATE);
 
  fftw_execute(np);
  fftw_execute(mg);
  fftw_execute(ef);
 
/*My slow DFT function! 
  Fourier(dipole, numTime+zeropad, dt, NPSpectrum, EV);
  Fourier(dipoleMG, numTime+zeropad, dt, MGSpectrum, EV);
  Fourier(efield, numTime+zeropad, dt, LaserSpectrum, EV);
*/

  
  FILE *absfp; 
  absfp = fopen("DATA/PRB_JPCC/MAbsorptionSpectrum.txt","w");
  fprintf(absfp, "#  Energy (ev)    Absorption NP      Absorption MG\n");
  
  int nfreq = 5001;
  // ~4.05 eV is max energy/ max freq
  double maxfreq = 50*0.08188379587298;
  double df = maxfreq / (nfreq - 1);
  double eps_0 = 1.0 / (4.0 * M_PI);
  printf("  numTime+zeropad %i\n",numTime+zeropad);
  printf("  min freq in au:  %12.10e\n",2*pi/(dt*(numTime+zeropad)));
  printf("  min freq in si:  %12.10e\n",2*pi/(dt*(numTime+zeropad))*omega_au);
  printf("  max freq in au:  %12.10e\n",2*pi/dt);
  printf("  max freq in si:  %12.10e\n",(2*pi/dt)*omega_au);
  for (int i=1; i<(numTime+zeropad); i++) {

    // This is omega in atomic units - same as energy in atomic units
    double omega = 2*pi*i/((numTime+zeropad)*dt);
    if (omega * 27.211>maxfreq) break;

    double omega_si = omega*omega_au;
    double eev = omega*27.211;
    double k = omega_si/2.99792458e+8;
    double pre_ext = k*k*k*k/(6*pi*8.854187e-12*8.854187e-12); 
    double pre_abs = k/(8.854187e-12);
    double complex npr = creal(dipole[i])/numTime;
    double complex npi = cimag(dipole[i])/numTime;

    double complex mgr = creal(dipoleMG[i])/numTime;
    double complex mgi = cimag(dipoleMG[i])/numTime;

    double complex efr = creal(efield[i])/numTime;
    double complex efi = cimag(efield[i])/numTime;

    //double pre = k/(8.854187e-12);
    double complex alphaNP = (npr+npi)/(efr+efi);
    double complex alphaMG = (mgr+mgi)/(efr+efi);
    //double sig_NP = pre*cimag(alphaNP);
    //double sig_MG = pre*cimag(alphaMG);
    double sig_ext_NP = pre_ext * creal(alphaNP*conj(alphaNP));
    double sig_ext_MG = pre_ext * creal(alphaMG*conj(alphaMG));
    double sig_abs_NP = pre_abs * cimag(alphaNP);
    double sig_abs_MG = pre_abs * cimag(alphaMG);

    // Going to print absorption cross section in abs/micron^2
    fprintf(absfp, "  %12.10e  %12.10e  %12.10e  %12.10e  %12.10e %12.10e  %12.10e\n",eev,npr, npi, mgr, mgi, efr, efi);
  }
  
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
