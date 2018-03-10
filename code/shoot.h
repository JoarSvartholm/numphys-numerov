
typedef struct {
  double EB,a,Re,mu,E;
}parameters;

typedef struct{
  double Y0,Y01,x0,x1,h;
  int N,NL,NR;
}simpars;


void shoot(parameters *p,simpars *s, double *psi,double *psivec);
