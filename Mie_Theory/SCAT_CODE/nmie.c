//**********************************************************************************//
// This library implements the algorithm for a multilayered sphere described by:    //
//   [1] Wen Yang, "Improved recursive algorithm for light scattering by a          //
//       multilayered sphere", Appl. Opt., 42 (2003) 1710.                          //
//                                                                                  //
// You can find the description of all the used equations in:                       //
//   [2] Ovidio Peña and Umapada Pal, "Scattering of electromagnetic radiation by   //
//       a multilayered sphere", Computer Physics Communications, Under review.     //
//                                                                                  //
// Hereinafter all equations numbers refer to [2]                                   //
//**********************************************************************************//
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include "ucomplex.h"
#include "nmie.h"

#define round(x) ((x) >= 0 ? (int)((x) + 0.5):(int)((x) - 0.5))

// Calculate Nstop - equation (17)
int Nstop(double xL) {
  int result;

  if (xL <= 8) {
    result = round(xL + 4*pow(xL, 1/3) + 1);
  } else if (xL <= 4200) {
    result = round(xL + 4.05*pow(xL, 1/3) + 2);
  } else {
    result = round(xL + 4*pow(xL, 1/3) + 2);
  }

  return result;
}

int Nmax(int L, double x[], complex m[]) {
  int i, result;

  result = Nstop(x[L + 1]);
  for (i = 1; i <= L; i++) {
    if (result < Cabs(RCmul(x[i], m[i]))) {
      result = round(Cabs(RCmul(x[i], m[i])));
    }
    if (result < Cabs(RCmul(x[i - 1], m[i]))) {
      result = round(Cabs(RCmul(x[i - 1], m[i])));
    }
  }

  return result + 15;
}

// Calculate an - equation (5)
complex calc_an(int n, double XL, complex Ha, complex mL, complex PsiXL, complex ZetaXL, complex PsiXLM1, complex ZetaXLM1) {
  complex Num = Csub(Cmul(Cadd(Cdiv(Ha, mL), Complex(n/XL, 0)), PsiXL), PsiXLM1);
  complex Denom = Csub(Cmul(Cadd(Cdiv(Ha, mL), Complex(n/XL, 0)), ZetaXL), ZetaXLM1);

  return Cdiv(Num, Denom);
}

// Calculate bn - equation (6)
complex calc_bn(int n, double XL, complex Hb, complex mL, complex PsiXL, complex ZetaXL, complex PsiXLM1, complex ZetaXLM1) {
  complex Num = Csub(Cmul(Cadd(Cmul(Hb, mL), Complex(n/XL, 0)), PsiXL), PsiXLM1);
  complex Denom = Csub(Cmul(Cadd(Cmul(Hb, mL), Complex(n/XL, 0)), ZetaXL), ZetaXLM1);

  return Cdiv(Num, Denom);
}

// Calculates S1_n - equation (25a)
complex calc_S1_n(int n, complex an, complex bn, double Pin, double Taun) {
  return RCmul((double)(n + n + 1)/(double)(n*n + n), Cadd(RCmul(Pin, an), RCmul(Taun, bn)));
}

// Calculates S2_n - equation (25b) (it's the same as (25a), just switches Pin and Taun)
complex calc_S2_n(int n, complex an, complex bn, double Pin, double Taun) {
  return calc_S1_n(n, an, bn, Taun, Pin);
}

//**********************************************************************************//
// This function calculates the actual scattering parameters and amplitudes         //
//                                                                                  //
// Input parameters:                                                                //
//   L: Number of layers                                                            //
//   x: Array containing the size parameters of the layers [1..L]                   //
//   m: Array containing the relative refractive indexes of the layers [1..L]       //
//   nTheta: Number of scattering angles                                            //
//   Theta: Array containing all the scattering angles where the scattering         //
//          amplitudes will be calculated                                           //
//                                                                                  //
// Output parameters:                                                               //
//   Qext: Efficiency factor for extinction                                         //
//   Qsca: Efficiency factor for scattering                                         //
//   Qabs: Efficiency factor for absorption (Qabs = Qext - Qsca)                    //
//   Qbk: Efficiency factor for backscattering                                      //
//   Qpr: Efficiency factor for the radiation pressure                              //
//   g: Asymmetry factor (g = (Qext-Qpr)/Qsca)                                      //
//   Albedo: Single scattering albedo (Albedo = Qsca/Qext)                          //
//   S1, S2: Complex scattering amplitudes                                          //
//                                                                                  //
// Return value:                                                                    //
//   Number of multipolar expansion terms used for the calculations                 //
//**********************************************************************************//

int nMie(int L, double x[], complex m[], int nTheta, double Theta[], double *Qext, double *Qsca, double *Qabs, double *Qbk, double *Qpr, double *g, double *Albedo, complex S1[], complex S2[]) {
  int n_max = Nmax(L, x, m);

  complex an, bn, anP1, bnP1, Qbktmp;

  complex D1_lmlx[n_max + 2][L + 1], D1_lmlxM1[n_max + 2][L + 1];
  complex D3_lmlx[n_max + 1][L + 1], D3_lmlxM1[n_max + 1][L + 1];
  complex D1XL[n_max + 2], D3XL[n_max + 1];
  complex PsiZeta_lmlx[n_max + 1][L + 1], PsiZeta_lmlxM1[n_max + 1][L + 1];
  complex PsiXL[n_max + 1], ZetaXL[n_max + 1], PsiZetaXL[n_max + 1];
  complex Q[n_max + 1][L + 1];
  complex Ha[n_max + 1][L + 1], Hb[n_max + 1][L + 1];
  double Pi[n_max + 1][nTheta], Tau[n_max + 1][nTheta];
  complex z1, z2;
  complex Num, Denom;
  complex G1, G2;
  complex Temp;
  double Tmp, x2;

  int n, l, t;

  // Initialize the scattering parameters
  *Qext = 0;
  *Qsca = 0;
  *Qabs = 0;
  *Qbk = 0;
  Qbktmp = Complex(0, 0);
  *Qpr = 0;
  *g = 0;
  *Albedo = 0;

  // Initialize Pi, Tau and the scattering amplitudes
  for (t = 0; t < nTheta; t++) {
    Pi[0][t] = 0.0;
    Tau[0][t] = 0.0;
    S1[t] = Complex(0, 0);
    S2[t] = Complex(0, 0);
  }

  //********************************************************//
  // Calculate D1, D3 and PsiZeta for z1 in the first layer //
  //********************************************************//
  z1 = RCmul(x[1], m[1]);

  // Downward recurrence for D1 - equations (16a) and (16b)
  D1_lmlx[n_max + 1][1] = Complex(0, 0);
  for (n = n_max + 1; n > 0; n--) {
    D1_lmlx[n - 1][1] = Csub(Cdiv(Complex(n, 0), z1), Cdiv(Complex(1, 0), Cadd(D1_lmlx[n][1], Cdiv(Complex(n, 0), z1))));
  }

  // Upward recurrence for PsiZeta and D3 - equations (18a) - (18d)
  PsiZeta_lmlx[0][1] = RCmul(0.5, Csub(Complex(1, 0), Cmul(Complex(cos(2*z1.r), sin(2*z1.r)), Complex(exp(-2*z1.i), 0))));
  D3_lmlx[0][1] = Complex(0, 1);
  for (n = 1; n <= n_max; n++) {
    PsiZeta_lmlx[n][1] = Cmul(PsiZeta_lmlx[n - 1][1], Cmul(Csub(Cdiv(Complex(n, 0), z1), D1_lmlx[n - 1][1]), Csub(Cdiv(Complex(n, 0), z1), D3_lmlx[n - 1][1])));

    D3_lmlx[n][1] = Cadd(D1_lmlx[n][1], Cdiv(Complex(0, 1), PsiZeta_lmlx[n][1]));
  }

  //******************************************************************//
  // Calculate Ha and Hb in the first layer - equations (7a) and (8a) //
  //******************************************************************//
  for (n = 1; n <= n_max; n++) {
    Ha[n][1] = D1_lmlx[n][1];
    Hb[n][1] = D1_lmlx[n][1];
  }

  //*******************************************//
  // Iteration from the layer 2 to the layer L //
  //*******************************************//
  for (l = 2; l <= L; l++) {
    //**************************************************************//
    //Calculate D1, D3 and PsiZeta for z1 and z2 in the layers 2..L //
    //**************************************************************//
    z1 = RCmul(x[l], m[l]);
    z2 = RCmul(x[l - 1], m[l]);

    // Downward recurrence for D1 - equations (16a) and (16b)
    D1_lmlx[n_max + 1][l] = Complex(0, 0);
    D1_lmlxM1[n_max + 1][l] = Complex(0, 0);
    for (n = n_max + 1; n > 0; n--) {
      D1_lmlx[n - 1][l] = Csub(Cdiv(Complex(n, 0), z1), Cdiv(Complex(1, 0), Cadd(D1_lmlx[n][l], Cdiv(Complex(n, 0), z1))));
      D1_lmlxM1[n - 1][l] = Csub(Cdiv(Complex(n, 0), z2), Cdiv(Complex(1, 0), Cadd(D1_lmlxM1[n][l], Cdiv(Complex(n, 0), z2))));
    }

    // Upward recurrence for PsiZeta and D3 - equations (18a) - (18d)
    PsiZeta_lmlx[0][l] = RCmul(0.5, Csub(Complex(1, 0), Cmul(Complex(cos(2*z1.r), sin(2*z1.r)), Complex(exp(-2*z1.i), 0))));
    PsiZeta_lmlxM1[0][l] = RCmul(0.5, Csub(Complex(1, 0), Cmul(Complex(cos(2*z2.r), sin(2*z2.r)), Complex(exp(-2*z2.i), 0))));

    D3_lmlx[0][l] = Complex(0, 1);
    D3_lmlxM1[0][l] = Complex(0, 1);

    for (n = 1; n <= n_max; n++) {
      PsiZeta_lmlx[n][l] = Cmul(PsiZeta_lmlx[n - 1][l], Cmul(Csub(Cdiv(Complex(n, 0), z1), D1_lmlx[n - 1][l]), Csub(Cdiv(Complex(n, 0), z1), D3_lmlx[n - 1][l])));
      PsiZeta_lmlxM1[n][l] = Cmul(PsiZeta_lmlxM1[n - 1][l], Cmul(Csub(Cdiv(Complex(n, 0), z2), D1_lmlxM1[n - 1][l]), Csub(Cdiv(Complex(n, 0), z2), D3_lmlxM1[n - 1][l])));

      D3_lmlx[n][l] = Cadd(D1_lmlx[n][l], Cdiv(Complex(0, 1), PsiZeta_lmlx[n][l]));
      D3_lmlxM1[n][l] = Cadd(D1_lmlxM1[n][l], Cdiv(Complex(0, 1), PsiZeta_lmlxM1[n][l]));
    }

    //******************************************//
    //Calculate Q, Ha and Hb in the layers 2..L //
    //******************************************//

    // Upward recurrence for Q - equations (19a) and (19b)
    Num = RCmul(exp(-2*(z1.i - z2.i)), Complex(cos(-2*z2.r) - exp(-2*z2.i), sin(-2*z2.r)));
    Denom = Complex(cos(-2*z1.r) - exp(-2*z1.i), sin(-2*z1.r));
    Q[0][l] = Cdiv(Num, Denom);

    for (n = 1; n <= n_max; n++) {
      Num = Cmul(Cadd(Cmul(z1, D1_lmlx[n][l]), Complex(n, 0)), Csub(Complex(n, 0), Cmul(z1, D3_lmlx[n - 1][l])));
      Denom = Cmul(Cadd(Cmul(z2, D1_lmlxM1[n][l]), Complex(n, 0)), Csub(Complex(n, 0), Cmul(z2, D3_lmlxM1[n - 1][l])));

      Tmp = (x[l - 1]*x[l - 1])/(x[l]*x[l]);

      Q[n][l] = Cdiv(Cmul(RCmul(Tmp, Q[n - 1][l]), Num), Denom);
    }

    // Upward recurrence for Ha and Hb - equations (7b), (8b) and (12) - (15)
    for (n = 1; n <= n_max; n++) {
      //Ha
      G1 = Csub(Cmul(m[l], Ha[n][l - 1]), Cmul(m[l - 1], D1_lmlxM1[n][l]));
      G2 = Csub(Cmul(m[l], Ha[n][l - 1]), Cmul(m[l - 1], D3_lmlxM1[n][l]));

      Temp = Cmul(Q[n][l], G1);

      Num = Csub(Cmul(G2, D1_lmlx[n][l]), Cmul(Temp, D3_lmlx[n][l]));
      Denom = Csub(G2, Temp);

      Ha[n][l] = Cdiv(Num, Denom);

      //Hb
      G1 = Csub(Cmul(m[l - 1], Hb[n][l - 1]), Cmul(m[l], D1_lmlxM1[n][l]));
      G2 = Csub(Cmul(m[l - 1], Hb[n][l - 1]), Cmul(m[l], D3_lmlxM1[n][l]));

      Temp = Cmul(Q[n][l], G1);

      Num = Csub(Cmul(G2, D1_lmlx[n][l]), Cmul(Temp, D3_lmlx[n][l]));
      Denom = Csub(G2, Temp);

      Hb[n][l] = Cdiv(Num, Denom);
    }
  }

  //************************************//
  //Calculate D1, D3 and PsiZeta for XL //
  //************************************//
  z1 = Complex(x[L], 0);

  // Downward recurrence for D1XL - equations (16a) and (16b)
  D1XL[n_max + 1] = Complex(0, 0);
  for (n = n_max + 1; n > 0; n--) {
    D1XL[n - 1] = Csub(Cdiv(Complex(n, 0), z1), Cdiv(Complex(1, 0), Cadd(D1XL[n], Cdiv(Complex(n, 0), z1))));
  }

  //Upward recurrence for PsiXL, ZetaXL and D3XL - equations (18b), (18d) and (20a) - (21b)
  PsiXL[0] = Complex(sin(z1.r), 0);
  ZetaXL[0] = Complex(sin(z1.r), -cos(z1.r));

  PsiZetaXL[0] = RCmul(0.5, Csub(Complex(1, 0), Cmul(Complex(cos(2*z1.r), sin(2*z1.r)), Complex(exp(-2*z1.i), 0))));
  D3XL[0] = Complex(0, 1);
  for (n = 1; n <= n_max; n++) {
    PsiXL[n] = Cmul(PsiXL[n - 1], Csub(Cdiv(Complex(n, 0), z1), D1XL[n - 1]));
    ZetaXL[n] = Cmul(ZetaXL[n - 1], Csub(Cdiv(Complex(n, 0), z1), D3XL[n - 1]));

    PsiZetaXL[n] = Cmul(PsiZetaXL[n - 1], Cmul(Csub(Cdiv(Complex(n, 0), z1), D1XL[n - 1]), Csub(Cdiv(Complex(n, 0), z1), D3XL[n - 1])));
    D3XL[n] = Cadd(D1XL[n], Cdiv(Complex(0, 1), PsiZetaXL[n]));
  }

  //*********************************************************************//
  // Finally, we calculate the an and bn coefficients and the resulting  //
  // scattering parameters                                               //
  //*********************************************************************//
  x2 = x[L]*x[L];

  anP1 = calc_an(1, x[L], Ha[1][L], m[L], PsiXL[1], ZetaXL[1], PsiXL[0], ZetaXL[0]);
  bnP1 = calc_bn(1, x[L], Hb[1][L], m[L], PsiXL[1], ZetaXL[1], PsiXL[0], ZetaXL[0]);
  for (n = 1; n <= n_max; n++) {
    an = anP1;
    bn = bnP1;

    anP1 = calc_an(n + 1, x[L], Ha[n + 1][L], m[L], PsiXL[n + 1], ZetaXL[n + 1], PsiXL[n], ZetaXL[n]);
    bnP1 = calc_bn(n + 1, x[L], Hb[n + 1][L], m[L], PsiXL[n + 1], ZetaXL[n + 1], PsiXL[n], ZetaXL[n]);

    // Equation (27)
    *Qext = *Qext + (double)(n + n + 1)*(an.r + bn.r);
    // Equation (28)
    *Qsca = *Qsca + (double)(n + n + 1)*(an.r*an.r + an.i*an.i + bn.r*bn.r + bn.i*bn.i);
    // Equation (29)
    *Qpr = *Qpr + ((n*(n + 2)/(n + 1))*((Cadd(Cmul(an, Conjg(anP1)), Cmul(bn, Conjg(bnP1)))).r) + ((double)(n + n + 1)/(n*(n + 1)))*(Cmul(an, Conjg(bn)).r));

    // Equation (33)
    Qbktmp = Cadd(Qbktmp, RCmul((double)((n + n + 1)*(1 - 2*(n % 2))), Csub(an, bn)));

    //****************************************************//
    // Calculate Pi_n and Tau_n for all values of Theta   //
    // Equations (26a) - (26c)                            //
    //****************************************************//
    for (t = 0; t < nTheta; t++) {
      Pi[n][t] = ((n == 1) ? 1.0 : (((double)(n + n - 1)*cos(Theta[t])*Pi[n - 1][t] - (double)n*Pi[n - 2][t])/((double)(n - 1))));
      Tau[n][t] = (double)n*cos(Theta[t])*Pi[n][t] - (double)(n + 1)*Pi[n - 1][t];

      S1[t] = Cadd(S1[t], calc_S1_n(n, an, bn, Pi[n][t], Tau[n][t]));
      S2[t] = Cadd(S2[t], calc_S2_n(n, an, bn, Pi[n][t], Tau[n][t]));
    }
  }

  *Qext = 2*(*Qext)/x2;                                 // Equation (27)
  *Qsca = 2*(*Qsca)/x2;                                 // Equation (28)
  *Qpr = *Qext - 4*(*Qpr)/x2;                           // Equation (29)

  *Qabs = *Qext - *Qsca;                                // Equation (30)
  *Albedo = *Qsca / *Qext;                              // Equation (31)
  *g = (*Qext - *Qpr) / *Qsca;                          // Equation (32)

  *Qbk = (Qbktmp.r*Qbktmp.r + Qbktmp.i*Qbktmp.i)/x2;    // Equation (33)

  return n_max;
}
