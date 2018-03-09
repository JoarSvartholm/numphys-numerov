#include <stdio.h>
#include <stdlib.h>
#include <math.h>

double V(double x){
  return 0.5*x*x;
}

int main(int argc, char const *argv[]) {
  int a,i,j,N=500;
  double psiLeft[N+1],psiRight[N+1],psi[2*N+1];
  double h=16./(2*N);
  double x,x0 = -8,x1 = 8;
  double Y0 = 0, Y01 = 1e-12;
  double Vj,Y,Y1,Ytemp;

  //Input control
  double E = atof(argv[1]);
  if((int)E%2)  a = -1; else a = 1;

  //Left shot
  psiLeft[0] = Y0;
  psiLeft[1] = Y01;
  Y = Y0; Y1 = Y01;
  for(j = 1;j<N;j++){
    x = x0 + (double)j*h;

    Vj = V(x);
    Ytemp = Y1*(2.0 - h*h*2.0*(E-Vj)/(1+h*h*(E-Vj)/6.0)) - Y;

    Y=Y1;
    Y1=Ytemp;
    psiLeft[j+1] = Y1/(1.0 + h*h*(E-V(x+h))/6.0);
  }
  //Right shot
  psiRight[0] = Y0;
  psiRight[1] = (double)a*Y01;
  Y = Y0; Y1 = (double)a*Y01;
  for(j = 1;j<N;j++){
    x = x1 - (double)j*h;

    Vj = V(x);
    Ytemp = Y1*(2.0 - h*h*2.0*(E-Vj)/(1+h*h*(E-Vj)/6.0)) - Y;

    Y=Y1;
    Y1=Ytemp;
    psiRight[j+1] = Y1/(1.0 + h*h*(E-V(x-h))/6.0);
  }

  //Solution match
  for(i=0;i<3;i++){
    if(abs((psiLeft[N-i]-psiRight[N-2+i])/psiLeft[N-i])>1e-8){
      fprintf(stderr, "Solution mismatch\n");
      exit(-1);
    }
  }

  //Construct solution vector
  for(i=0;i<N+1;i++){
    psi[i] = psiLeft[i];
    psi[2*N-i] = psiRight[i];
  }

  //Normalization
  double norm = 0.5*psi[0]*psi[0];
  for(i=1;i<2*N;i++){
    norm +=psi[i]*psi[i];
  }
  norm+=psi[2*N]*psi[2*N];
  for(i=0;i<2*N+1;i++){
    psi[i]/=sqrt(norm);
  }

  //Print result
  for(i=0;i<2*N+1;i++){
    printf("%e\n", psi[i] );
  }


  return 0;
}
