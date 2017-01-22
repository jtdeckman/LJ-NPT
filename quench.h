#include <cmath>
#define SIGN(a,b) ((b) > 0.0 ? fabs(a) : -fabs(a))
#define SHFT(a,b,c,d) (a)=(b);(b)=(c);(c)=(d);
#define MAX(a,b) ((a) > (b) ? (a) : (b))

typedef unsigned int uint;

class Optim {

    int ITMAX;
    static const int ITMAXB=100;
    static const int TOO_MANY_ITERATIONS=4;
    static const double EPS=1.0e-10;
    static const double ZEPS=1.0e-10; 
    static const double TINY=1.0e-20;
    static const double TOL=2.0e-4;
    static const double CGOLD=0.3819660;
    static const double GOLD=1.618034;
    static const double GLIMIT=100.0;

    int fitstatus, ncom, nmolcom, N, Nm, Nm3;
    double coe, eps, sigma, sigma2, sigma6, sigma12, roe, rc, k;
    double coe6, coe12;

    double *pcom, *xicom, *g, *h, *xi, *xt;
 
    double (Optim::*nrfunc)(double [], uint N);

    double f1dim(double x);
    void linmin(double p[], double xi[], int n, double *fret, double (Optim::*func)(double [], uint));
    void mnbrak(double *ax,double *bx,double *cx,double *fa,double *fb,double *fc, double (Optim::*func)(double));
    double brent(double ax, double bx, double cx, double (Optim::*f)(double), double tol, double *xmin);
    int getfitstatus(void);
    void  setfitstatus(int thestatus);

  public:

    Optim(double, double, int, int);
    ~Optim() {};
    void init(uint size, double sig);
 
    void frprmn1(double p[], int n,double ftol, double &grm, int *iter,double *fret, double (Optim::*func)(double [], uint),
    void (Optim::*dfunc)(double [], double [], uint)); 
    
    void gradLJ(double*, double*, uint); 
    double energyLJ(double*, uint);

};


Optim::Optim(double sig = 1.0, double ep = 1.0, int n = 10, int maxi=1000) : sigma(sig), eps(ep), Nm(n), Nm3(3*n), N(n), ITMAX(maxi) {  
      
       pcom = new double[Nm3];
       xicom = new double[Nm3];
       g = new double[Nm3];
       h = new double[Nm3];
       xi = new double[Nm3];
       xt = new double[Nm3];

       coe = 4*eps;
       sigma2 = sigma*sigma;
       sigma6 = sigma2*sigma2*sigma2;
       sigma12 = sigma6*sigma6;

       coe6 = coe*24.0*sigma6;
       coe12 = coe*48.0*sigma12;

       fitstatus=0;
       ncom=0; 
       nmolcom=0;
   
       for(int i=0; i<Nm3; i++) {
          pcom[i] = 0;
          xicom[i] = 0;
          g[i] = 0;
          h[i] = 0;
          xi[i] = 0;
          xt[i] = 0;
       }

}      

double Optim::energyLJ(double* Q, uint N ) {

  double rsq=0, E=0;
  double r6, r12;

  for(int i=0; i<N; i++)
    for(int j=i+1; j<N; j++) {
      for(int k=0; k<3; k++) 
        rsq += (Q[i*3+k]-Q[j*3+k])*(Q[i*3+k]-Q[j*3+k]);
      
      r6 = rsq*rsq*rsq;
      r12 = r6*r6;
      E +=  coe*(sigma12/r12  - sigma6/r6);
      rsq = 0;
    }

   
  return(E);
}

void Optim::gradLJ(double* r, double* gr, uint n) {

  double rsq, rsq4, rsq7, gfact, rvec[3];
 
  for(uint i=0; i<3*n; i++)
    gr[i]=0;

  for(uint i=0; i < n-1; i++) {  

    for(uint j=i+1; j < n; j++) {  

      rsq = 0.0;

      for(uint k=0; k<3; k++) {
        rvec[k] = r[3*i+k] - r[3*j+k];
        rsq += rvec[k]*rvec[k];
      }

      rsq4 = rsq*rsq*rsq*rsq;
      rsq7 = (rsq4*rsq4)/rsq;

      gfact = coe6/rsq4 - coe12/rsq7;

      for(uint k=0; k<3; k++) {
        rsq = gfact*rvec[k];
        gr[3*i+k] -= rsq;
        gr[3*j+k] += rsq;
      }
    }
  }
}	

int Optim::getfitstatus(void)
 {
    return(fitstatus);
 }
 
 void  Optim::setfitstatus(int thestatus)
 {
   fitstatus = thestatus;
 }

 double Optim::brent(double ax, double bx, double cx, double (Optim::*f)(double), 
 double tol, double *xmin)
{
  int iter;
  double a,b,d=0,etemp,fu,fv,fw,fx,p,q,r,tol1,tol2,u,v,w,x,xm;
  double e=0.0;

   a=((ax < cx) ? ax : cx);
   b=((ax > cx) ? ax : cx);
   x=w=v=bx;
   fw=fv=fx=(this->*f)(x);
   for (iter=1;iter<=ITMAXB;iter++) {
     xm=0.5*(a+b);
     tol2=2.0*(tol1=tol*fabs(x)+ZEPS);
     if (fabs(x-xm) <= (tol2-0.5*(b-a))) {
       *xmin=x;
       return fx;
     }
     if (fabs(e) > tol1) {
       r=(x-w)*(fx-fv);
       q=(x-v)*(fx-fw);
       p=(x-v)*q-(x-w)*r;
       q=2.0*(q-r);
       if (q > 0.0) p = -p;
       q=fabs(q);
       etemp=e;
      e=d;
       if (fabs(p) >= fabs(0.5*q*etemp) || p <= q*(a-x) || p >= q*(b-x))
      d=CGOLD*(e=(x >= xm ? a-x : b-x));
     else {
        d=p/q;
         u=x+d;
        if (u-a < tol2 || b-u < tol2)
           d=SIGN(tol1,xm-x);
       }
     } else {
       d=CGOLD*(e=(x >= xm ? a-x : b-x));
     }
     u=(fabs(d) >= tol1 ? x+d : x+SIGN(tol1,d));
    fu=(this->*f)(u);
     if (fu <= fx) {
       if (u >= x) a=x; else b=x;
       SHFT(v,w,x,u);
       SHFT(fv,fw,fx,fu);
     } else {
       if (u < x) a=u; else b=u;
       if (fu <= fw || w == x) {
         v=w;
         w=u;
         fv=fw;
         fw=fu;
       } else if (fu <= fv || v == x || v == w) {
         v=u;
        fv=fu;
       }
     }
   }
   setfitstatus(getfitstatus() + TOO_MANY_ITERATIONS);
   *xmin=x;
   return fx;
 }


void Optim::mnbrak(double *ax,double *bx,double *cx,double *fa,double *fb,double *fc, double (Optim::*func)(double))
 {
   double ulim,u,r,q,fu,dum;
   *fa=(this->*func)(*ax);
   *fb=(this->*func)(*bx);
   if (*fb > *fa) {
   SHFT(dum,*ax,*bx,dum);
     SHFT(dum,*fb,*fa,dum);
   }
  *cx=(*bx)+GOLD*(*bx-*ax);
   *fc=(this->*func)(*cx);
   while (*fb > *fc) {
     r=(*bx-*ax)*(*fb-*fc);
     q=(*bx-*cx)*(*fb-*fa);
    u=(*bx)-((*bx-*cx)*q-(*bx-*ax)*r)/
       (2.0*SIGN(MAX(fabs(q-r),TINY),q-r));
     ulim=(*bx)+GLIMIT*(*cx-*bx);
     if ((*bx-u)*(u-*cx) > 0.0) {
      fu=(this->*func)(u);
       if (fu < *fc) {
         *ax=(*bx);
        *bx=u;
        *fa=(*fb);
         *fb=fu;
         return;
      } else if (fu > *fb) {
        *cx=u;
        *fc=fu;
       return;
       }
       u=(*cx)+GOLD*(*cx-*bx);
       fu=(this->*func)(u);
     } else if ((*cx-u)*(u-ulim) > 0.0) {
       fu=(this->*func)(u);
       if (fu < *fc) {
        SHFT(*bx,*cx,u,*cx+GOLD*(*cx-*bx))
        SHFT(*fb,*fc,fu,(this->*func)(u))
       }
     } else if ((u-ulim)*(ulim-*cx) >= 0.0) {
       u=ulim;
      fu=(this->*func)(u);
     } else {
       u=(*cx)+GOLD*(*cx-*bx);
       fu=(this->*func)(u);
     }
     SHFT(*ax,*bx,*cx,u);
    SHFT(*fa,*fb,*fc,fu);
   }
 }
 

void Optim::linmin(double p[], double xi[], int n, double *fret, double (Optim::*func)(double [], uint N))
 {
  int j;
  double xx,xmin,fx,fb,fa,bx,ax;
  ncom=n;
  nrfunc=func;
   for (j=0;j<3*n;j++) {
    pcom[j]=p[j];
     xicom[j]=xi[j];
   }
   ax=0.0;
   xx=1.0;
   bx=2.0;
  mnbrak(&ax,&xx,&bx,&fa,&fx,&fb,&Optim::f1dim);
  *fret= brent(ax,xx,bx,&Optim::f1dim,TOL,&xmin);
   for (j=0;j<3*n;j++) {
     xi[j] *= xmin;
     p[j] += xi[j];
   }
 }

void Optim::frprmn1(double p[], int n,double ftol, double& grm, int *iter,double *fret,       
 double (Optim::*func)(double [], uint n), void (Optim::*dfunc)(double [], double [], uint n)) // Energy func. (energyLJ(), Gradient func. (gradLJ())
 {
   int j,its;
   double gg,gam,fp,dgg;
  
   fp=(this->*func)(p,n);
   (this->*dfunc)(p,xi,n);
   for (j=0;j<3*n;j++) {
     g[j] = -xi[j];
    xi[j]=h[j]=g[j];
   }
  for (its=1;its<=ITMAX;its++) {
     *iter=its;
     linmin(p,xi,n,fret,func);
     if (2.0*fabs(*fret-fp) <= ftol*(fabs(*fret)+fabs(fp)+EPS)) {
       (this->*dfunc)(p,xi,n);
       grm=fabs(xi[0]);
       for(int i=1;i<3*n;i++) { if(fabs(xi[i]) > grm){grm=xi[i];} }
       return;
     }
     fp=(this->*func)(p,n);
    (this->*dfunc)(p,xi,n);
    dgg=gg=0.0;
    for (j=0; j<3*n; j++) {
      gg += g[j]*g[j];
      dgg += (xi[j]+g[j])*xi[j];
     }
    if (gg == 0.0) {
       return;
     }
    gam=dgg/gg;
    for (j=0;j<3*n;j++) {
       g[j] = -xi[j];
       xi[j]=h[j]=g[j]+gam*h[j];
     }
   }
   setfitstatus(getfitstatus() + TOO_MANY_ITERATIONS);
 }

double Optim::f1dim(double x)
 {
   int j;
   double f; 
   for (j=0;j<3*ncom;j++) xt[j]=pcom[j]+x*xicom[j];
   f=(this->*nrfunc)(xt,ncom);
   return f;
 }
