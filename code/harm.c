#include <stdio.h>
#include <stdlib.h>

double V(double x){
  return x*x*0.5;
}

int main(int argc, char const *argv[]) {

  if(argc<2) {fprintf(stderr, "Wrong number of parameters\n" ); exit(-1);}
  double E = atof(argv[1]);

  double Y0 = 0,Y1 = 1e-12, x0 = -8;
  double psi,Vj,Yt,x;
  double h = 16./1000;

  if(argc<3){
    printf("%e\n %e\n", Y0,Y1 );
    for(int j = 2;j<=500;j++){
      x = x0 + (double)j*h;

      Vj = V(x);
      Yt = Y1*(2.0 - h*h*2.0*(E-Vj)/(1+h*h*(E-Vj)/6.0)) - Y0;

      Y0=Y1;
      Y1=Yt;

      psi = Y1/(1.0 + h*h*(E-V(x+h))/6.0);

      printf("%e\n",psi );

    }

    Y0=0,Y1=1e-12;
    printf("%e\n %e\n", Y0,Y1 );
    for(int j = 2;j<=500;j++){
      x = x0 + 16 - (double)j*h;

      Vj = V(x);
      Yt = Y1*(2.0 - h*h*2.0*(E-Vj)/(1+h*h*(E-Vj)/6.0)) - Y0;

      Y0=Y1;
      Y1=Yt;

      psi = Y1/(1.0 + h*h*(E-V(x+h))/6.0);

      printf("%e\n",psi );

    }
  }else{
    printf("%e\n %e\n", Y0,Y1 );
    for(int j = 2;j<=1000;j++){
      x = x0 + (double)j*h;

      Vj = V(x);
      Yt = Y1*(2.0 - h*h*2.0*(E-Vj)/(1+h*h*(E-Vj)/6.0)) - Y0;

      Y0=Y1;
      Y1=Yt;

      psi = Y1/(1.0 + h*h*(E-V(x+h))/6.0);

      printf("%e\n",psi );
    }
  }

  return 0;
}
