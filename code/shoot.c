#include <stdio.h>
#include <math.h>
#include <gsl/gsl_const_mksa.h>
#include "shoot.h"

double V(double R,void *params){
  parameters pars = *(parameters *) params;
  //return pars.a*pars.a*pars.EB*(R-pars.Re)*(R-pars.Re);
  return pars.EB*(1-exp(-pars.a*(R-pars.Re)))*(1-exp(-pars.a*(R-pars.Re)));
}


void shoot(parameters *p,simpars *s, double *psi,double *psivec){
  double psiRight[2*s->N +2], psiLeft[s->N+2];
  double Vj,Y,Y1,Ytemp,xj,kj;
  double h = s->h;
  double E = p->E;
  double hbar = GSL_CONST_MKSA_PLANCKS_CONSTANT_HBAR;


  //Left shot
  psiLeft[0] = s->Y0/s->x0;
  psiLeft[1] = s->Y01/(s->x0+h);
  Y = s->Y0; Y1 = s->Y01;
  for(int j = 1;j<=0.7*s->N;j++){
    xj = s->x0 + (double)j*h;

    Vj = V(xj,p);
    kj = 2.0*p->mu*(p->E-Vj)/(hbar*hbar);
    Ytemp = Y1*(2.0 - h*h*kj/(1.+h*h*kj/12.)) - Y;

    Y=Y1;
    Y1=Ytemp;
    //kj = 2.0*p->mu/(hbar*hbar)*(p->E-V(xj+h,p));
    //psiLeft[j+1] = Y1/(1 + h*h*kj/12.);
    psiLeft[j+1] = Y1/((1.0 + p->mu/(hbar*hbar)*h*h*(E-V(xj+h,p))/6.0));
    psiLeft[j+1] /= (xj+h);
    //if(j==1) printf("%d %e %e\n", j, xj, psiLeft[j+1] );
    printf("%e \n",psiLeft[j+1]);
    //if(j>=s->N-2) printf("%e\n", xj+h );
  }
  //Right shot
  psiRight[0] = s->Y0/s->x1;
  psiRight[1] = s->Y01/(s->x1-h);
  Y = s->Y0; Y1 = s->Y01;
  for(int j = 1;j<=1.5*s->N;j++){
    xj = s->x1 - (double)j*h;

    Vj = V(xj,p);
    kj = 2.0*p->mu*(p->E-Vj)/(hbar*hbar);
    Ytemp = Y1*(2.0 - h*h*kj/(1.+h*h*kj/12.)) - Y;

    Y=Y1;
    Y1=Ytemp;
    psiRight[j+1] = Y1/((1.0 + p->mu/(hbar*hbar)*h*h*(E-V(xj-h,p))/6.0));
    psiRight[j+1] /= (xj-h);

    //printf("%e \n",psiRight[j+1]);
    //if(j>=s->N-2) printf("%e\n", xj-h );
  }

  for(int i=0;i<s->N+2;i++){
    psiRight[i]/=psiLeft[s->N];
    psiLeft[i]/=psiLeft[s->N];
  }

  //Construct solution vector
  for(int i=0;i<s->N+1;i++){
    psi[i] = psiLeft[i];
    psi[2*s->N-i] = psiRight[i];
  }

  psivec[0] = psiLeft[s->N-1];
  psivec[1] = psiLeft[s->N+1];
  psivec[2] = psiRight[s->N+1];
  psivec[3] = psiRight[s->N-1];

  for(double i=0;i<2*s->N;i++){
    xj = s->x0 + h*i;
    Vj = V(xj,p);
  //  printf("%e\n",Vj );
  }


}
