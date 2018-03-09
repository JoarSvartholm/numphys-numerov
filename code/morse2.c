#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <gsl/gsl_const_mksa.h>
#include "shoot.h"





int main(int argc, char const *argv[]) {
  simpars s;
  parameters p;

  double mH = 1.67374e-27;
  double mCl = 5.88715e-26;
  p.EB = 7.392e-19;
  p.a = 1.812e10;
  p.Re = 1.275e-10;
  p.mu = mH*mCl/(mH+mCl);
  double hbar = GSL_CONST_MKSA_PLANCKS_CONSTANT_HBAR;

  s.x0 = 0.5e-10;
  s.x1 = 3e-10;
  s.Y0 = 0;
  s.Y01 = 1e-12;
  s.N = 800;
  s.h = (s.x1-s.x0)/(2*s.N);


  double psi[2*s.N+1],psivec[4];



  //Input control
  p.E = hbar*p.a*sqrt(0.5*p.EB/p.mu)*0.5;
  //if((int)p.E%2)  a = -1; else a = 1;

  shoot(&p,&s,psi,psivec);

  double F = pow(psivec[0]-psivec[2],2) + pow(psivec[1]-psivec[3],2);
  //printf("%e %e %e %e %e\n",F,psivec[0],psivec[1],psivec[2],psivec[3] );

  //Solution match
  for(int i=0;i<3;i++){
  //  if(abs((psiLeft[N-i]-psiRight[N-2+i])/psiLeft[N-i])>1e-8){
    //  fprintf(stderr, "Solution mismatch\n");
    //  exit(-1);
    //}
  }


  //Normalization
  double norm = 0.5*psi[0]*psi[0];
  for(int i=1;i<2*s.N;i++){
    norm +=psi[i]*psi[i];
  }
  norm+=psi[2*s.N]*psi[2*s.N];
  for(int i=0;i<2*s.N+1;i++){
    psi[i]/=sqrt(norm);
  }

  //Print result
  for(int i=0;i<2*s.N+1;i++){
    printf("%e\n", psi[i] );
  }



  return 0;
}
