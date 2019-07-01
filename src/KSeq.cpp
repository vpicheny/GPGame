#include <Rcpp.h>
using namespace Rcpp;

//' @param Z matrix ([npts x nobj] of objective values)
//' @param Nadir,Shadow vectors defining the line intersecting the Pareto front
//' @noRd
// [[Rcpp::export]]
int getKS_cpp(NumericMatrix Z, NumericVector Nadir, NumericVector Shadow) {
  NumericVector SmN(Nadir.length());
  SmN = Shadow - Nadir;
  int nobj = Z.ncol();
  int nr = Z.nrow();
  int iKS = 0;
  double maxmin = -INFINITY;
  double tmp;
  double* ptrZ = &Z(0,0);
  double* ptrN = &Nadir(0);
  double* ptrSmN = &SmN(0);
  
  for(int i = 0; i < nr; i++, ptrZ++){
    tmp = INFINITY;
    for(int j = 0; j < nobj; j++, ptrN++, ptrSmN++){
      if((*ptrZ - *ptrN) / *ptrSmN < tmp){
        tmp = (*ptrZ - *ptrN) / *ptrSmN;
      }
      ptrZ += nr;
    }
    ptrZ -= nr * nobj;
    ptrN -= nobj;
    ptrSmN -= nobj;
    
    if(maxmin < tmp){
      maxmin = tmp;
      iKS = i;
    }
  }
  
  return(iKS + 1);
}
