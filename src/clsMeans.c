#include <R.h>
#include <Rinternals.h>

SEXP clsMeans(SEXP mM, SEXP vCls, SEXP vKcls) {
    int m = nrows(mM);
    int p = ncols(mM);
    int k = length(vKcls);
    SEXP clsMean = PROTECT(allocMatrix(REALSXP, k, p));
    SEXP clsSize = PROTECT(allocVector(INTSXP, k));
	double *pclsMean = REAL(clsMean);
	int *pclsSize = INTEGER(clsSize);
    memset(pclsMean, 0, k * p * sizeof(double));
	memset(pclsSize, 0, k * sizeof(int));
	int *pCls = INTEGER(vCls);
	int *pKcls = INTEGER(vKcls);
	double *pM = REAL(mM);
	
    for(int i = 0; i < m; i++) {
        if(!R_IsNA(pCls[i])) {
            for(int j = 0; j < k; j++) {
                if(pCls[i] == pKcls[j]) {
                    for(int a = 0; a < p; a++) {
                        pclsMean[j + k*a] += pM[i + m*a];
                    }
                    pclsSize[j]++;
                }
            }
        }
    }
    
    for(int j = 0; j < k; j++) {
        for(int a = 0; a < p; a++) {
            pclsMean[j + k*a] = pclsMean[j + k*a]/pclsSize[j];
        }
    }
    
	UNPROTECT(2);
    return clsMean;
}