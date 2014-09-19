#include <R.h>
#include <Rinternals.h>

SEXP rowWhichMaxC(SEXP x) {
  int m = nrows(x);
  int n = ncols(x);
  double *pout, *px;
  SEXP out = PROTECT(allocMatrix(REALSXP, m, 2));
  
  pout = REAL(out);
  px = REAL(x);
  
  for (int i = 0; i < m; i++) {
    pout[0 + i] = px[0 + i];           //out(i, 0)
    pout[m + i] = 1;                  //out(i, 1)
    for (int j = 1; j < n; j++) {
        if (px[j*m + i] > pout[i]) {
            pout[i] = px[j*m + i];
    		pout[i + m] = j + 1;
        }
    }
  }
  
  UNPROTECT(1);
  return out;
}
