#include<stdio.h>
#include<stdlib.h>
#include<math.h>
<<<<<<< HEAD
int main() { 
int i, MAX_I;
 double x, dx, fx,fpx,fxf,xf;

//rise = fxf-fx
//run = xf-x
//slope = (fxf-fx)/(xf-x)


 double pi = 4*atan(1.0);
//double pi = 3.141592653589793;

MAX_I = 8; 
dx = 10./MAX_I;
 // evaluate and print f(x) = -(x-4)^2 + 4
 for (i=0; i<MAX_I; i++ ){

  x = i*dx;
  xf = (i+1)*dx;

  
 // fx = -1*(x-4)*(x-4) + 4;
  fx = sqrt(2/10.)*sin(pi*x/10.);
  fxf = sqrt(2/10.)*sin(pi*xf/10.);

   ror = (fxf-fx)/(xf-x);

//calculate analytical derivatives

fpx = sqrt(2/10.)*(pi/10)*cos(pi*x/10);

double rms=

  printf(" %i %f %f \n ",i,x, ror, fpx);


}
=======

//calculating derivative of psi_1(x)
int main() {

  int i, MAX_I;
  double x, xf, dx, fx, fpx, fxf, ror;
  //rise = fxf-fx
  //run = xf - x
  //slope = (fxf-fx)/(xf-x)
 
  //double pi = 4*atan(1.0)
  double pi = 3.141592653589793;
 

  MAX_I =1000;
  dx = 10./MAX_I; 
  // evaluate and print f(x) = -(x-4)^2 + 4

  for (i=0; i<=MAX_I; i++) { 
  	// calculate numerical derivative 
  	x = i*dx;
	xf = (i+1)*dx;
	
	fx= sqrt(2/10.)*sin(pi*x/10.);
	fxf= sqrt (2/10.)*sin(pi*xf/10.);
	ror= (fxf-fx)/(xf-x);
	// calculate analytical derivative 
	
  	fpx = sqrt(2/10.)*(pi/10)*cos(pi*x/10);

        double rms = sqrt( (fpx-ror)*(fpx-ror));
 	printf("%i %f %f %f %f\n ",i, x,ror, fpx, rms);
    }
>>>>>>> 28caba614b148dd8375462107b8ec4cda4db370c
}
