#include <Rcpp.h>
using namespace Rcpp;
// localLoss
// 
// @keywords internal
// 
// @description Calculates the local loss
// 
// @return a list containing two matrices containing the local loss
//' @import Rcpp
// [[Rcpp::export]]
List localLoss(const NumericMatrix &L,const NumericMatrix &R, 
               const NumericVector &is, const NumericVector &js, 
               const NumericMatrix &error_matrix) {
    
    NumericMatrix dL(L.nrow(), L.ncol());
    NumericMatrix dR(R.nrow(), R.ncol());
    for(int n = 0; n < is.size(); n++){
        int row = is[n] - 1;
        int column  = js[n] - 1;
        double x = error_matrix(row, column);
        for(int lIndex = 0; lIndex < dL.row(row).size(); lIndex++){
            dL.row(row)[lIndex] = dL.row(row)[lIndex] + x * R.column(column)[lIndex];
            dR.column(column)[lIndex] = dR.column(column)[lIndex] + x * L.row(row)[lIndex];
        }

    }
    return List::create(  _["dL"]  = dL, _["dR"]  = dR );
}


