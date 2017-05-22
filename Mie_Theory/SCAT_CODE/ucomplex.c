#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include "ucomplex.h"

complex Cadd(complex z1, complex z2) {
  complex zr;
  zr.r = z1.r + z2.r;
  zr.i = z1.i + z2.i;
  return zr;
}

complex Csub(complex z1, complex z2) {
  complex zr;
  zr.r = z1.r - z2.r;
  zr.i = z1.i - z2.i;
  return zr;
}

complex Cmul(complex z1, complex z2) {
  complex zr;
  zr.r = (z1.r*z2.r) - (z1.i*z2.i);
  zr.i = (z1.r*z2.i) + (z1.i*z2.r);
  return zr;
}

complex RCmul(double r, complex z) {
  complex zr;
  zr.r = z.r*r;
  zr.i = z.i*r;
  return zr;
}

complex Cdiv(complex z1, complex z2) {
/* The following algorithm is used to properly handle
   denominator overflow:

              |  a + b(d/c)   c - a(d/c)
              |  ---------- + ----------- I     if |d| < |c|
   a + b I    |  c + d(d/c)   a + d(d/c)
   -------  = |
   c + d I    |  b + a(c/d)   -a + b(c/d)
              |  ---------- + ----------- I     if |d| >= |c|
              |  d + c(c/d)   d + c(c/d)
*/
  complex zr;
  double tmp, denom;

  if (fabs(z2.r) > fabs(z2.i)) {
    tmp = z2.i/z2.r;
    denom = z2.r + z2.i*tmp;
    zr.r = (z1.r + z1.i*tmp)/denom;
    zr.i = (z1.i - z1.r*tmp)/denom;
  } else {
    tmp = z2.r/z2.i;
    denom = z2.i + z2.r*tmp;
    zr.r = (z1.i + z1.r*tmp)/denom;
    zr.i = (-z1.r + z1.i*tmp)/denom;
  }
  return zr;
}    

complex Complex(double re, double im) {
  complex zr;
  zr.r = re;
  zr.i = im;
  return zr;
}

double Cabs (complex z) {
  double r;
  r = sqrt((z.r*z.r) + (z.i*z.i));
  return r;
}

double Carc (complex z) {
  double r;
  r = atan2(z.i, z.r);
  return r;
}

complex Conjg (complex z) {
  complex zr;
  zr.r = z.r;
  zr.i = -z.i;
  return zr;
}

complex Cinv (complex z) {
  complex zr;
  double denom;

  denom = (z.r*z.r) + (z.i*z.i);
  zr.r = z.r/denom;
  zr.i = -z.i/denom;
  return zr;
}

complex Cexp (complex z) {
/* exp(x + iy) = exp(x).exp(iy) = exp(x).[cos(y) + i sin(y)] */
  complex zr;
  double expz;

  expz = exp(z.r);
  zr.r = expz*cos(z.i);
  zr.i = expz*sin(z.i);
  return zr;
}

complex Clog (complex z) {
/* ln( p exp(i0)) = ln(p) + i0 + 2kpi */
  complex zr;
  zr.r = log(Cabs(z));
  zr.i = atan2(z.i, z.r);
  return zr;
}

complex Csqrt (complex z) {
  complex zr;
  double root, q;

  if ((z.r != 0.0) || (z.i != 0.0)) {
    root = sqrt(0.5*(fabs(z.r) + Cabs(z)));
    q = z.i/(2.0*root);
    if (z.r >= 0.0) {
      zr.r = root;
      zr.i = q;
    } else if (z.i < 0.0) {
             zr.r = -q;
             zr.i = -root;
           } else {
             zr.r =  q;
             zr.i =  root;
           }
    } else zr = z;
  return zr;
}

/* complex trigonometric functions */

complex Ccos (complex z) {
/* complex cosinus */
/* cos(x+iy) = cos(x).cos(iy) - sin(x).sin(iy) */
/* cos(ix) = cosh(x) et sin(ix) = i.sinh(x) */
  complex zr;
  zr.r = cos(z.r)*cosh(z.i);
  zr.i = -sin(z.r)*sinh(z.i);
  return zr;
}

complex Csin (complex z) {
/* sinus complex */
/* sin(x+iy) = sin(x).cos(iy) + cos(x).sin(iy) */
/* cos(ix) = cosh(x) et sin(ix) = i.sinh(x) */
  complex zr;
  zr.r = sin(z.r)*cosh(z.i);
  zr.i = cos(z.r)*sinh(z.i);
  return zr;
}

complex Ctan (complex z) {
/* tangente */
  return Cdiv(Csin(z), Ccos(z));
}

/* inverse trigonometric functions */

complex Carc_cos (complex z) {
/* arc cosinus complex */
/* arccos(z) = -i.argch(z) */
  return Cmul(Complex(0.0, -1.0), Carc_ch(z));
}

complex Carc_sin (complex z) {
/* arc sinus complex */
/* arcsin(z) = -i.argsh(i.z) */
  return Cmul(Complex(0.0, -1.0), Carc_sh(Cmul(Complex(0.0, 1.0), z)));
}

complex Carc_tan (complex z) {
/* arc tangente complex */
/* arctg(z) = -i.argth(i.z) */
  return Cmul(Complex(0.0, -1.0), Carc_th(Cmul(Complex(0.0, 1.0), z)));
}

 /* hyberbolic complex functions */

complex Cch (complex z) {
/* hyberbolic cosinus */
/* cosh(x+iy) = cosh(x).cosh(iy) + sinh(x).sinh(iy) */
/* cosh(iy) = cos(y) et sinh(iy) = i.sin(y) */
  complex zr;
  zr.r = cosh(z.r)*cos(z.i);
  zr.i = sinh(z.r)*sin(z.i);
  return zr;
}

complex Csh (complex z) {
/* hyberbolic sinus */
/* sinh(x+iy) = sinh(x).cosh(iy) + cosh(x).sinh(iy) */
/* cosh(iy) = cos(y) et sinh(iy) = i.sin(y) */
  complex zr;
  zr.r = sinh(z.r)*cos(z.i);
  zr.i = cosh(z.r)*sin(z.i);
  return zr;
}

complex Cth (complex z) {
/* hyberbolic complex tangent */
/* th(x) = sinh(x)/cosh(x) */
/* cosh(x) > 1 qq x */
  complex temp, zr;
  temp = Cch(z);
  zr = Csh(z);
  return Cdiv(zr, temp);
}

/* inverse complex hyperbolic functions */

complex Carc_ch (complex z) {
/*   hyberbolic arg cosinus */
/*                          _________  */
/* argch(z) = -/+ ln(z + i.V 1 - z.z)  */
  complex zr;
  zr = Clog(Cadd(z, Cmul(Complex(0.0, 1.0), Csqrt(Csub(Cmul(z, z), Complex(1.0, 0.0))))));
  zr.r = -zr.r;
  zr.i = -zr.i;
  return zr;
}

complex Carc_sh (complex z) {
/*   hyperbolic arc sinus       */
/*                    ________  */
/* argsh(z) = ln(z + V 1 + z.z) */
  return Clog(Cadd(z, Csqrt(Cadd(Cmul(z, z), Complex(1.0, 0.0)))));
}

complex Carc_th (complex z) {
/* hyperbolic arc tangent */
/* argth(z) = 1/2 ln((z + 1)/(1 - z)) */
  return Cdiv(Clog(Cdiv(Cadd(z, Complex(1.0, 0.0)), Csub(z, Complex(1.0, 0.0)))), Complex(2.0, 0.0));
}

