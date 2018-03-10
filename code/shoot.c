#include <stdio.h>
#include <math.h>
#include <gsl/gsl_const_mksa.h>
#include "shoot.h"

double V(double R,void *params){
  parameters pars = *(parameters *) params;
  //return pars.a*pars.a*pars.EB*(R-pars.Re)*(R-pars.Re);
  return pars.EB*(1-exp(-pars.a*(R-pars.Re)))*(1-exp(-pars.a*(R-pars.Re)));
}


void shoot(parameters *p,simpars *s, double *psi,double *psivec,double E){
  double psiRight[s->NR +2], psiLeft[s->NL+2];
  double Vj,Y,Y1,Ytemp,xj,kj;
  double h = s->h;
  double hbar = GSL_CONST_MKSA_PLANCKS_CONSTANT_HBAR;


  //Left shot
  psiLeft[0] = s->Y0/s->x0;
  psiLeft[1] = s->Y01/(s->x0+h);
  Y = s->Y0; Y1 = s->Y01;
  for(int j = 1;j<=s->NL;j++){
    xj = s->x0 + (double)j*h;

    Vj = V(xj,p);
    kj = 2.0*p->mu*(E-Vj)/(hbar*hbar);
    Ytemp = Y1*(2.0 - h*h*kj/(1.+h*h*kj/12.)) - Y;

    Y=Y1;
    Y1=Ytemp;
    //kj = 2.0*p->mu/(hbar*hbar)*(E-V(xj+h,p));
    //psiLeft[j+1] = Y1/(1 + h*h*kj/12.);
    psiLeft[j+1] = Y1/((1.0 + p->mu/(hbar*hbar)*h*h*(E-V(xj+h,p))/6.0));
    psiLeft[j+1] /= (xj+h);
    //if(j==1) printf("%d %e %e\n", j, xj, psiLeft[j+1] );
  }
  //Right shot
  psiRight[0] = s->Y0/s->x1;
  psiRight[1] = s->Y01/(s->x1-h);
  Y = s->Y0; Y1 = s->Y01;
  for(int j = 1;j<=s->NR;j++){
    xj = s->x1 - (double)j*h;

    Vj = V(xj,p);
    kj = 2.0*p->mu*(E-Vj)/(hbar*hbar);
    Ytemp = Y1*(2.0 - h*h*kj/(1.+h*h*kj/12.)) - Y;

    Y=Y1;
    Y1=Ytemp;
    psiRight[j+1] = Y1/((1.0 + p->mu/(hbar*hbar)*h*h*(E-V(xj-h,p))/6.0));
    psiRight[j+1] /= (xj-h);

    //printf("%e \n",psiRight[j+1]);
  }

  double psim = psiRight[s->NR];
  for(int i=0;i<s->NR+2;i++){
    psiRight[i]=psiRight[i]/psim;
  }
  psim = psiLeft[s->NL];
  for(int i=0;i<s->NL+2;i++){
    psiLeft[i]=psiLeft[i]/psim;
  }

  //Construct solution vector
  for(int i=0;i<s->NL+1;i++){
    psi[i] = psiLeft[i];
  }
  for(int i=0;i<s->NR+1;i++){
    psi[2*s->N-i] = psiRight[i];
  }


  psivec[0] = psiLeft[s->NL-1];
  psivec[1] = psiLeft[s->NL+1];
  psivec[2] = psiRight[s->NR+1];
  psivec[3] = psiRight[s->NR-1];


}
