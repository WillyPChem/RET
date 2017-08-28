#include <fftw3.h>
#include<math.h>
#include<stdio.h>
#include<malloc.h>
#include<complex.h>

#define REAL 0
#define IMAG 1
//void Commutator (int dim, double H[dim][dim], double D[dim][dim], double P[dim][dim]);
void RK3(int Nlevel, double time, double *bas, double *E, double *Hint, double *Mu_x, double *Mu_y, double *Mu_z, double *Dis, double complex *D, double dt);
void Liouville(int dim, double complex *H, double complex *D, double complex *P);
void AntiCommutator(int dim, double *H, double complex *D, double complex *P);
void PrintComplexMatrix(int dim, double complex *M);
void PrintRealMatrix(int dim, double *M);
void FormDM(int dim, double complex *C, double complex *Ct, double complex *DM);
void L_Diss(int Nlevel, double *gamma, double complex *D, double *bas, double complex *P);
void Fourier (double complex *dm, int n, double dt, double complex *ftout, double *freqvec);
double complex TrMuD(int Nlevel, double *Mu, double complex *D);
double E_Field(double time, int comp);
void FillDFTArray(int fftw_iter, double real_val, double imag_val, fftw_complex* ar);
double D_Error(int dim, double complex *D);
// Function prototype for H_interaction
void H_interaction(int dim, double *Hint, double *mu_x, double *mu_y, double *mu_z, double dpm_x, double dpm_y, double dpm_z, double *R);
void DipoleAcceleration(int dim, double dt, fftw_complex* dp, fftw_complex* dpa);

// NOTE!!!  You need three global variables for the rates associated with 
// The Lindblad operators - gamma, beta, and alpha should be defined here according
// to their values in the journal of physical chemistry letters paper
//int Nlevel = 3;
//int dim = Nlevel*Nlevel;
// Some important constants
double pi = 4*atan(1.0);
double wn_to_au = 4.55633528e-6; 
double mu_au_to_si = 8.47835326e-30; // 1 a.u. = 8.47835326e-30 C m
double E_au_to_si = 5.14220652e11;  // 1 a.u. = 5.14220652e11 V/m
double omega_au = 4.134147e+16;;
int main() {

  // Timestep for propagation (in atomic units)
  double dt = 0.01;
  // Number of timesteps for propagation
  int numTime = 2000000;
  // Number of additional timesteps that will be "zero-padded" in FT vectors
  int zeropad = 2000000;
  // Total length of FT vectors
  int dft_dim = numTime+zeropad;

  // Hamiltonian terms for NP (donor)
  double *E, *Mu_x, *Mu_y, *Mu_z, *Dis, *bas, *Hint;
  double complex *H, *D, *P;

  // Nlevel -> number of states for NP, dim -> Dimension of Nanoparticle arrays
  int Nlevel, dim;
  
  // NP levels can be variable in principle
  printf("  How many states are in your NP system? \n");
  scanf("%i",&Nlevel);
  dim = Nlevel*Nlevel;
  // MG variables here!

  // NlevelMG -> number of states for MG, dimMG -> Dimension of MG arrays
  int NlevelMG = 3;
  int dimMG = NlevelMG*NlevelMG;

  // Hamiltonian terms for MG (acceptor)
  double *EMG, *MuMG_x, *MuMG_y, *MuMG_z, *MuZERO, *DisMG, *basMG, *HintMG;
  double complex *HMG, *DMG, *PMG;

  // Expectation value variables
  double tr, trMG;
  double complex dipole_moment_x, dipole_moment_y, dipole_moment_z, dipole_moment;
  double complex dipole_momentMG_x, dipole_momentMG_y, dipole_momentMG_z, dipole_momentMG;

  // Separation vector
  double *r;

  // Electric field variables
  double Efx, Efy, Efz, Enorm;

  // File pointers!
  // For printing dipole expectation value data (scalar quantity)
  FILE *dfp, *dfpMG;
  // For printing Population and Energy Transfer data 
  FILE *popfp, *popMGfp, *ecumfp;
  // For reading information about energy eigenstates
  FILE *Efp, *EfpMG;
  // For reading information about excited-state lifetimes
  FILE *Disfp, *DisfpMG;
  // For reading information about transition dipole moment components
  FILE *Muxfp, *Muyfp, *Muzfp, *MuMGxfp, *MuMGyfp, *MuMGzfp;
  // For writing absorption, scattering, and emission spectra
  FILE *absfp, *emsfp;
  // FFTW variables here -> inputs to fft
  fftw_complex *dipole;
  fftw_complex *dipoleA;
  fftw_complex *efield;
  fftw_complex *dipoleMG;
  fftw_complex *dipoleMGA;
  fftw_complex *nps;
  fftw_complex *mgs;
  fftw_complex *efs;
  fftw_complex *emiss;
  fftw_complex *emissMG;

  // Allocate memory for FFTW arrays
  dipole = (fftw_complex*)malloc(dft_dim*sizeof(fftw_complex));
  dipoleA= (fftw_complex*)malloc(dft_dim*sizeof(fftw_complex));
  efield = (fftw_complex*)malloc(dft_dim*sizeof(fftw_complex));
  dipoleMG = (fftw_complex*)malloc(dft_dim*sizeof(fftw_complex));
  dipoleMGA= (fftw_complex*)malloc(dft_dim*sizeof(fftw_complex));
  nps = (fftw_complex*)malloc(dft_dim*sizeof(fftw_complex));
  mgs = (fftw_complex*)malloc(dft_dim*sizeof(fftw_complex));
  efs = (fftw_complex*)malloc(dft_dim*sizeof(fftw_complex));
  emiss = (fftw_complex*)malloc(dft_dim*sizeof(fftw_complex));
  emissMG = (fftw_complex*)malloc(dft_dim*sizeof(fftw_complex));

  // This determines the type of FT that will be executed later
  fftw_plan npp = fftw_plan_dft_1d(dft_dim,
                                      dipole,
                                      nps,
                                      FFTW_BACKWARD,
                                      FFTW_ESTIMATE);

  fftw_plan mgp = fftw_plan_dft_1d(dft_dim,
                                      dipoleMG,
                                      mgs,
                                      FFTW_BACKWARD,
                                      FFTW_ESTIMATE);


  fftw_plan efp = fftw_plan_dft_1d(dft_dim,
                                      efield,
                                      efs,
                                      FFTW_BACKWARD,
                                      FFTW_ESTIMATE);

  fftw_plan emp = fftw_plan_dft_1d(dft_dim, 
 					dipoleA,
					emiss,
					FFTW_BACKWARD,
					FFTW_ESTIMATE);

  fftw_plan emgp = fftw_plan_dft_1d(dft_dim,
					dipoleMGA,
					emissMG,
					FFTW_BACKWARD,
					FFTW_ESTIMATE); 


  // Allocate memory for all other arrays
  // NP
  H = (double complex*)malloc(dim*sizeof(double complex));
  D = (double complex *)malloc(dim*sizeof(double complex));
  P = (double complex *)malloc(dim*sizeof(double complex));
  E  = (double *)malloc(dim*sizeof(double));
  Mu_x = (double *)malloc(dim*sizeof(double));
  Mu_y = (double *)malloc(dim*sizeof(double));
  Mu_z = (double *)malloc(dim*sizeof(double));
  Dis = (double *)malloc(dim*sizeof(double));
  bas = (double *)malloc(dim*sizeof(double));
  Hint = (double *)malloc(dim*sizeof(double));
  // MG
  HMG = (double complex*)malloc(dimMG*sizeof(double complex));
  DMG = (double complex *)malloc(dimMG*sizeof(double complex));
  PMG = (double complex *)malloc(dimMG*sizeof(double complex));
  EMG  = (double *)malloc(dimMG*sizeof(double));
  MuMG_x = (double *)malloc(dimMG*sizeof(double));
  MuMG_y = (double *)malloc*(dimMG*sizeof(double));
  MuMG_z = (double *)malloc*(dimMG*sizeof(double));
  MuZERO = (double *)malloc(dimMG*sizeof(double));
  DisMG = (double *)malloc(dimMG*sizeof(double));
  basMG = (double *)malloc(dimMG*sizeof(double));
  HintMG = (double *)malloc(dimMG*sizeof(double));

  // Allocate memory for separation vector
  r = (double *)malloc(3*sizeof(double));
  // Define the separation between NP and MG
  r[0] = 500000000.;
  r[1] = 0.;
  r[2] = 0.;

  // INPUT FILES!
  // Nanoparticle inputfiles
  //Efp = fopen("Matrices/PLASMON/Energy_Au.txt","r");
  //Mufp = fopen("Matrices/PLASMON/Dipole_Au.txt","r");
  //Disfp = fopen("Matrices/PLASMON/Dissipation_Au.txt","r");
  Efp = fopen("Matrices/SMA_PEAK1/Energy8s.txt","r");
  Muxfp = fopen("Matrices/SMA_PEAK1/Dipole8s_x.txt","r");
  Muyfp = fopen("Matrices/SMA_PEAK1/Dipole8s_y.txt","r");
  Muzfp = fopen("Matrices/SMA_PEAK1/Dipole8s_z.txt","r");
  Disfp = fopen("Matrices/SMA_PEAK1/Dissipation8s.txt","r");

  // MG inputfiles
  EfpMG = fopen("Matrices/SMA_PEAK1/Energy.txt","r");
  MuMGxfp = fopen("Matrices/SMA_PEAK1/Dipole_x.txt","r");
  MuMGyfp = fopen("Matrices/SMA_PEAK1/Dipole_y.txt","r");
  MuMGzfp = fopen("Matrices/SMA_PEAK1/Dipole_z.txt","r");
  DisfpMG = fopen("Matrices/SMA_PEAK1/Dissipation.txt","r");

  // OUTPUT FILES!  MAKE SURE NAMES ARE CONSISTENT
  dfp = fopen("DATA/SMA_PEAK1/DipoleMoment_SMAAir_inf.dat","w");
  dfpMG = fopen("DATA/SMA_PEAK1/DipoleMomentMG_SMAAir_inf.dat", "w");
  popfp = fopen("DATA/SMA_PEAK1/Population_SMAAir_inf.dat","w");
  popMGfp = fopen("DATA/SMA_PEAK1/PopulationMG_SMAAir_inf.dat","w");
  ecumfp = fopen("DATA/SMA_PEAK1/CumulativeEnergyTransfer_SMAAir_inf.dat","w");
  emsfp = fopen("DATA/SMA_PEAK1/EmissionSpectru_SMAAir_inf.dat","w");
  absfp = fopen("DATA/SMA_PEAK1/AbsorptionSpectrum_SMAAir_inf.dat","w");


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
        basMG[i*NlevelMG+j] = 1.0;
      }
      else{
        basMG[i*NlevelMG+j] = 0.;
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

       fscanf(Muxfp,"%lf",&val);
       Mu_x[i] = val;

       fscanf(Muyfp,"%lf",&val);
       Mu_y[i] = val;
  
       fscanf(Muzfp,"%lf",&val);
       Mu_z[i] = val;

       fscanf(Disfp,"%lf",&val);
       Dis[i] = val;
  }
  // MG
  for (int i=0; i<dimMG; i++) {

       fscanf(EfpMG,"%lf",&val);
       EMG[i] = val;

       fscanf(MuMGxfp,"%lf",&val);
       MuMG_x[i] = val;

       fscanf(MuMGyfp,"%lf",&val);
       MuMG_y[i] = val;
  
       fscanf(MuMGzfp,"%lf",&val);
       MuMG_z[i] = val;

       MuZERO[i] = 0.;
       fscanf(DisfpMG,"%lf",&val);
       DisMG[i] = val;


  }

  // Print parameters to screen
  printf("\nE\n");
  PrintRealMatrix(Nlevel, E);
  printf("\nMu_x\n");
  PrintRealMatrix(Nlevel,Mu_x);
  printf("\nMu_y\n");
  PrintRealMatrix(Nlevel,Mu_y);
  printf("\nMu_z\n");
  PrintRealMatrix(Nlevel,Mu_z);
  printf("\nDis\n");
  PrintRealMatrix(Nlevel,Dis);
  printf("\nBas\n");
  PrintRealMatrix(Nlevel,bas);
  printf("\nDM\n");
  PrintComplexMatrix(Nlevel,D);

  printf("\nEMG\n");
  PrintRealMatrix(NlevelMG, EMG);
  printf("\nMuMG_x\n");
  PrintRealMatrix(NlevelMG, MuMG_x);
  printf("\nMuMG_y\n");
  PrintRealMatrix(NlevelMG, MuMG_y);
  printf("\nMuMG_z\n");
  PrintRealMatrix(NlevelMG, MuMG_z);
  printf("\nDisMG\n");
  PrintRealMatrix(NlevelMG, DisMG);
  printf("\nBasMG\n");
  PrintRealMatrix(NlevelMG, basMG);
  printf("\nDMG\n");
  PrintComplexMatrix(NlevelMG, DMG);


  // Get initial dipole moments
  dipole_moment_x = TrMuD(Nlevel, Mu_x, D)*mu_au_to_si;
  dipole_moment_y = TrMuD(Nlevel, Mu_y, D)*mu_au_to_si;
  dipole_moment_z = TrMuD(Nlevel, Mu_z, D)*mu_au_to_sil
  dipole_moment = csqrt(dipole_moment_x*dipole_moment_x + dipole_moment_y*dipole_moment_y + dipole_moment_z*dipole_moment_z);

  dipole_momentMG_x = TrMuD(NlevelMG, MuMG_x, DMG)*mu_au_to_si;
  dipole_momentMG_y = TrMuD(NlevelMG, MuMG_y, DMG)*mu_au_to_si;
  dipole_momentMG_z = TrMuD(NlevelMG, MuMG_z, DMG)*mu_au_to_si;
  dipole_momentMG = csqrt(dipole_momentMG_x*dipole_momentMG_x + dipole_momentMG_y*dipole_momentMG_y + dipole_momentMG_z*dipole_momentMG_z);

  FillDFTArray(0, creal(dipole_moment), cimag(dipole_moment), dipole);
  FillDFTArray(0, creal(dipole_momentMG), cimag(dipole_momentMG), dipoleMG);
  FillDFTArray(0, 0., 0., efield);


  //void H_interaction(int dim, double *Hint, double *mu, double dpm, double R) 
  H_interaction(Nlevel, Hint, Mu_x, Mu_y, Mu_z, creal(dipole_momentMG_x), creal(dipole_momentMG_y), creal(dipole_momentMG_z), r); 
  H_interaction(NlevelMG, HintMG, MuMG_x, MuMG_y, MuMG_z, creal(dipole_moment_x), creal(dipole_moment_y), creal(dipole_moment_z), r);

  double max_MG_Error, MG_Error;
  double EnMG, E_Transfer, TransferTime;
  double e1_curr, e1_prev, e1_cum, e2_curr, e2_prev, e2_cum;
  e1_curr = 0.;
  e1_prev = 0.;
  e1_cum = 0.;
  e2_curr = 0.;
  e2_prev = 0.;
  e2_cum = 0.;
  max_MG_Error = D_Error(NlevelMG*NlevelMG, DMG);
  for (int i=1; i<numTime; i++) {

    // Get amount of energy transferred into the excited states of MG as of last step
    e1_prev = creal(DMG[1*NlevelMG+1])*EMG[1*NlevelMG+1];
    e2_prev = creal(DMG[2*NlevelMG+2])*EMG[2*NlevelMG+2];

    // How perturbed is the Density Matrix on MG?
    MG_Error = D_Error(NlevelMG*NlevelMG, DMG);

    // What is the energy of MG system now?
    EnMG = creal(TrMuD(NlevelMG, EMG, DMG));

    // If more perturbed than before, update max_MG_Error
    // And E_Transfer
    if (MG_Error > max_MG_Error) {

      max_MG_Error = MG_Error;
      E_Transfer = EnMG;
      TransferTime = i*dt;

    }
    
    // Update the NP 
    RK3(Nlevel, dt*i, bas, E, Hint, Mu, Dis, D, dt);

    // Calculate dipole moment of NP
    dipole_moment = TrMuD(Nlevel, Mu, D); 
    FillDFTArray(i, creal(dipole_moment*mu_au_to_si), cimag(dipole_moment*mu_au_to_si), dipole);

    // Calculate interaction matrix between NP and MG: H_int^{np->mg}    
    H_interaction(NlevelMG, HintMG, MuMG_x, MuMG_y, MuMG_z, creal(dipole_moment_x), creal(dipole_moment_y), creal(dipole_moment_z), r);
    
    // Update the MG
    RK3(NlevelMG, dt*i, basMG, EMG, HintMG, MuMG, DisMG, DMG, dt);
    //RK3(NlevelMG, dt*i, basMG, EMG, HintMG, MuZERO, DisMG, DMG, dt);
    
    // Calculate the dipole moment of MG
    dipole_momentMG = TrMuD(NlevelMG, MuMG, DMG);
    FillDFTArray(i, creal(dipole_momentMG*mu_au_to_si), cimag(dipole_momentMG*mu_au_to_si), dipoleMG);
    
    // Calculate interaction matrix between MG and NP: H_int^{mg->np} 
    H_interaction(Nlevel, Hint, Mu_x, Mu_y, Mu_z, creal(dipole_momentMG_x), creal(dipole_momentMG_y), creal(dipole_momentMG_z), r); 
   
    // Get amount of energy tranfserred into the excited states of MG as of now
    e1_curr = creal(DMG[1*NlevelMG+1])*EMG[1*NlevelMG+1];
    e2_curr = creal(DMG[2*NlevelMG+2])*EMG[2*NlevelMG+2];

    // Is there more than before?  If so, incrementally increase cumulative energy transferred to that state
    if (e1_curr>e1_prev) e1_cum+=(e1_curr-e1_prev);
    if (e2_curr>e2_prev) e2_cum+=(e2_curr-e2_prev);

    // Print population data!
    fprintf(popfp,"\n %f ",dt*i);
    fprintf(popMGfp,"\n %f ",dt*i);
    tr=0.;
    trMG = 0.;
    for (int j=0; j<Nlevel; j++) {

      fprintf(popfp," %12.10e",creal(D[j*Nlevel+j]));
      tr+=creal(D[j*Nlevel+j]);

    }
    for (int j=0; j<NlevelMG; j++) {


      fprintf(popMGfp,"  %12.10e",creal(DMG[j*NlevelMG+j]));
      trMG+=creal(DMG[j*NlevelMG+j]);

    }
    fprintf(popfp," %12.10e",tr);
    fprintf(popMGfp," %12.10e",trMG);
  
    // Print data on energy transferred to MG (from field and from NP)
    fprintf(ecumfp, "%f   %12.10e  %12.10e  %12.10e  %12.10e\n",dt*i,e1_curr, e1_cum,e2_curr, e2_cum);

    // Print dipole moment of NP
    fprintf(dfp," %f  %12.10e  %12.10e\n",dt*i,creal(dipole_moment),cimag(dipole_moment));
  
    // Print dipole moment of MG
    fprintf(dfpMG,"%f %12.10e %12.10e\n",dt*i,creal(dipole_momentMG),cimag(dipole_momentMG));

    Efx = E_Field(dt*i, 0);
    Efy = E_Field(dt*i, 1);
    Efz = E_Field(dt*i, 2);
    Enorm = sqrt(Efx*Efx+Efy*Efy+Efz*Efz); 
    FillDFTArray(i,  E_au_to_si*Enorm, 0, efield);
  }
  
  printf("  E_mg_1 is %12.10f: \n", EMG[3*1+1]*27.211);
  printf("  E_mg_2 is %12.10f: \n", EMG[3*2+2]*27.211); 
  printf("  Max(||D(0)-D(t)||) is %12.10e at t=%12.10e\n",max_MG_Error,TransferTime);

  printf("  r_x          r_y         r_z          E_PEAK(eV)       E_e1(eV)          E_e2(eV)          E_tot(eV)\n");
  printf("  %6.3e    %6.3e   %6.3e    %12.10e %12.10e  %12.10e  %12.10e\n",r[0],r[1],r[2],E_Transfer*27.211, e1_cum*27.211, e2_cum*27.211, (e1_cum+e2_cum)*27.211);


  // Compute dipole accelration for NP and MG system 
  DipoleAcceleration(numTime, dt, dipole, dipoleA);
  DipoleAcceleration(numTime, dt, dipoleMG, dipoleMGA);

  // Now to the spectra!
  for (int i=numTime; i<zeropad; i++) {

    FillDFTArray(i, 0., 0., dipole);
    FillDFTArray(i, 0., 0., dipoleMG);
    FillDFTArray(i, 0., 0., efield);
    FillDFTArray(i, 0., 0., dipoleA);
    FillDFTArray(i, 0., 0., dipoleMGA);

  }
 
  fftw_execute(npp);
  fftw_execute(mgp);
  fftw_execute(efp);
  fftw_execute(emp);
  fftw_execute(emgp);
 
  
  fprintf(absfp, "#  Energy (ev)    SCAT NP      SCAT MG       ABS NP       ABS MG\n");
  
  int nfreq = 5001;
  // ~6.48 eV is max energy/ max freq
  double maxfreq = 70*0.08188379587298;
  double df = maxfreq / (nfreq - 1);
  double eps_0 = 1.0 / (4.0 * M_PI);
  for (int i=1; i<(numTime+zeropad); i++) {

    // This is omega in atomic units - same as energy in atomic units
    double omega = 2*pi*i/((numTime+zeropad)*dt);
    if (omega * 27.211>maxfreq) break;

    double omega_si = omega*omega_au;
    double eev = omega*27.211;
    double k = omega_si/2.99792458e+8;
    double pre_scat = k*k*k*k/(6*pi*8.854187e-12*8.854187e-12); 
    double pre_abs = k/(8.854187e-12);
    double npr = nps[i][0]/numTime;
    double npi = nps[i][1]/numTime;

    double mgr = mgs[i][0]/numTime;
    double mgi = mgs[i][1]/numTime;

    double efr = efs[i][0]/numTime;
    double efi = efs[i][1]/numTime;

    double emnpr = emiss[i][0]/numTime;
    double emnpi = emiss[i][1]/numTime;
    double emmgr = emissMG[i][0]/numTime;
    double emmgi = emissMG[i][1]/numTime;

    double sig_emiss_np = creal( emnpr*emnpr+emnpi*emnpi );
    double sig_emiss_mg = creal( emmgr*emmgr+emmgi*emmgi );

    double complex alphaNP = (npr+I*npi)/(efr+I*efi);
    double complex alphaMG = (mgr+I*mgi)/(efr+I*efi);

    double sig_scat_NP = pre_scat * creal(alphaNP*conj(alphaNP));
    double sig_scat_MG = pre_scat * creal(alphaMG*conj(alphaMG));
    double sig_abs_NP = pre_abs * cimag(alphaNP);
    double sig_abs_MG = pre_abs * cimag(alphaMG);

    // Going to print absorption and scattering cross section in m^2
    fprintf(absfp, "  %12.10e  %12.10e  %12.10e  %12.10e  %12.10e %12.10e  %12.10e\n",eev,sig_scat_NP, sig_scat_MG, sig_abs_NP, sig_abs_MG, efr, efi);
    fprintf(emsfp, "  %12.10e  %12.10e  %12.10e  \n",eev, sig_emiss_np, sig_emiss_mg);
  }
  
  fclose(absfp);
  
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
void RK3(int Nlevel, double time, double *bas, double *E, double *Hint, double *Mu_x, double *Mu_y, double *Mu_z, double *Dis, double complex *D, double dt) {

  int i, j;
  double complex *D_dot, *D2, *D3, *D_np1, *k1, *k2, *k3;
  double complex *H, *LD;  // Contribution to Ddot from Lindblad dissipation
  double *gamma;
  int dim = Nlevel*Nlevel; 
  double Efield_x, Efield_y, Efield_z;

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
  Efield_x = E_Field(time, 0);
  Efield_y = E_Field(time, 1);
  Efield_z = E_Field(time, 2);

  // Compute full Hamiltonian at current time t 
  for (i=0; i<dim; i++) {


      H[i] = E[i] + Hint[i] - Efield_x*Mu_x[i] - Efield_y*Mu_y[i] - Efield_z*Mu_z[i]; // - I*Dis[i];

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
  Efield_x = E_Field(time+dt/2., 0);
  Efield_y = E_Field(time+dt/2., 1);
  Efield_z = E_Field(time+dt/2., 2);

  // Compute full Hamiltonian at partially updated time t 
  for (i=0; i<dim; i++) {

      H[i] = E[i] + Hint[i] - Efield_x*Mu_x[i] - Efield_y*Mu_y[i] - Efield_z*Mu_z[i]; // - I*Dis[i];

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
  free(LD);
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


double E_Field(double time, int comp) {

  double Ef;
  double tau = 75.;

  // Define polarization of source here!
  // Conventionally, light travels along x-direction so has no polarization along x
  double p_x, p_y, p_z;
  p_x = 0;
  p_y = sqrt(1./2);
  p_z = sqrt(1./2);

  if (time<tau) {

    if (comp==0) {

      Ef = p_x*0.0003*sin(time*pi/tau)*sin(time*pi/tau)*sin(0.07423*time);
    
    }
    else if (comp==1) {

      Ef = p_y*0.0003*sin(time*pi/tau)*sin(time*pi/tau)*sin(0.07423*time);

    }
    else if (comp==2) {

      Ef = p_z*0.0003*sin(time*pi/tau)*sin(time*pi/tau)*sin(0.07423*time);

    }
    // Defaults to z-polarized
    else {

      Ef = p_z*0.0003*sin(time*pi/tau)*sin(time*pi/tau)*sin(0.07423*time);

    }
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

//  Needs x, y, and z-components of dipole array and dipole moment 
void H_interaction(int dim, double *Hint, double *mu_x, double *mu_y, double *mu_z, double dpm_x, double dpm_y, double dpm_z, double *R) {
  
  int i; 
 // double *tmp1, *tmp2;
  double oer2, oer3;
  double scal_R = sqrt(R[0]*R[0] + R[1]*R[1] + R[2]*R[2]);
  double t1, t2; 
  oer2 = pow(scal_R,-2.);
  oer3 = pow(scal_R,-3.);
 
  // Write code between here!
 
 for (i=0; i<dim*dim; i++){

   // Using components of dipole operator and dipole moment now!
   t1 = dpm_x*mu_x[i] + dpm_y*mu_y[i] + dpm_z*mu_z[i];
   t2 = (R[0]*mu_x[i] + R[1]*mu_y[i] + R[2]*mu_z[i])*(R[0]*dpm_x + R[1]*dpm_y + R[2]*dpm_z);
   Hint[i] = oer3*(t1 - 3*t2/oer2);
 } 

}

void FillDFTArray(int fftw_iter, double real_val, double imag_val, fftw_complex* ar) {

  ar[fftw_iter][REAL] = real_val;
  ar[fftw_iter][IMAG] = imag_val;

}

// Computes the Frobenius norm of the difference between the current density matrix and the initial density matrix
double D_Error(int dim, double complex *D) {

  double complex e;
  double complex one = 1. + 0.*I;
  double complex zero = 0. + 0.*I;
  double FN;

  e = 0.+0*I;
 
  e = (one - D[0])*conj(one-D[0]);

  for (int i=1; i<dim; i++) {

    e += (zero - D[i])*conj(zero-D[i]);  

  }

  FN = creal(csqrt(e));

  return FN;
}

void DipoleAcceleration(int dim, double dt, fftw_complex* dp, fftw_complex* dpa) {

  double complex st1, st2, st3, st4, st5, accel;

  // initialize stencil variables to zero
  st1 = 0.+0.*I;
  st2 = 0.+0.*I;
  st3 = 0.+0.*I;
  st4 = 0.+0.*I;
  st5 = 0.+0.*I;  
  for (int i=0; i<dim; i++) {

    st1 = st2;
    st2 = st3;
    st3 = st4;
    st5 = dp[i][REAL] + dp[i][IMAG]*I;
    accel = (-1./12.)*(st1+st5) + (4./3.)*(st2+st4)-(5./2.)*st3;
    accel /= (dt*dt);
    dpa[i][REAL] = creal(accel);
    dpa[i][IMAG] = cimag(accel);

    //printf("  %12.10e  %12.10e  %12.10e  %12.10e  %12.10e\n",dt*i,creal(st5),cimag(st5),creal(accel),cimag(accel));    

  }

}

