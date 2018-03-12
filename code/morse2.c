#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <gsl/gsl_const_mksa.h>
#include "shoot.h"





int main(int argc, char const *argv[]) {
  simpars s;
  parameters p;
  double F,E;
  int iter=0,maxiter = 100;

  double mH = 1.67374e-27;
  double mCl = 5.88715e-26;
  p.EB = 7.392e-19;
  p.a = 1.812e10;
  p.Re = 1.275e-10;
  p.mu = mH*mCl/(mH+mCl);
  double hbar = GSL_CONST_MKSA_PLANCKS_CONSTANT_HBAR;
  double e = GSL_CONST_MKSA_ELECTRON_CHARGE;

  s.x0 = 0.5e-10;
  s.x1 = 3e-10;
  s.Y0 = 0;
  s.Y01 = 1e-12;
  s.N = 500;
  s.h = (s.x1-s.x0)/(2*s.N);
  s.NL = 350;
  s.NR = 650;

  double psi[2*s.N+1],psivec[4];

  //Input control
  double n = atof(argv[1]);
  E = hbar*p.a*sqrt(2.*p.EB/p.mu)*(n+0.5) - pow(hbar*p.a*sqrt(2.*p.EB/p.mu)*(n+0.5),2)/(4.*p.EB);
  double dE = 0;

do{

  E = E+dE;
  shoot(&p,&s,psi,psivec,E);
  F = pow(psivec[0]-psivec[2],2) + pow(psivec[1]-psivec[3],2);

  //printf("%e %e  %d   %e %e %e %e\n",F,E,iter,  psivec[0],psivec[1],psivec[2],psivec[3]);
  shoot(&p,&s,psi,psivec,E+0.01*E);
  double Fp = pow(psivec[0]-psivec[2],2) + pow(psivec[1]-psivec[3],2);
  shoot(&p,&s,psi,psivec,E-0.01*E);
  double Fm = pow(psivec[0]-psivec[2],2) + pow(psivec[1]-psivec[3],2);
  double FD = (Fp-Fm)/(0.02*E);
  dE = -F/FD;

  iter++;
}while(F>1e-10 && iter<maxiter);

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
  if(argc<3){
  for(int i=0;i<2*s.N+1;i++){
    if(psi[1]<0) printf("%e\n", -psi[i] );
    else printf("%e\n", psi[i] );
  }
}else printf("n = %d , Eigenenergy = %e eV\n", (int)n,E/e );

  return 0;
}
