#include "quench.h"

extern "C" {
double c_quench(double *q, int n, int maxi, int *nit) {
  
  double fret, grm=10000;
  int numit, ITMAX;
  const double ftol = 1e-14;

  Optim qc(1.0, 1.0, n, maxi);

  qc.frprmn1(q, n, ftol, grm, &numit, &fret, &Optim::energyLJ, &Optim::gradLJ);

  *nit = numit;

  return grm;
}
}
