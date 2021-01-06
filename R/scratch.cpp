#include <Rcpp.h>
#include <Rmath.h>
using namespace Rcpp;
 
// [[Rcpp::export]]
void test() { 
  //NumericVector normSam = rnorm(1e4);
  IntegerVector data = seq(1,1e4);

  IntegerVector fold = seq(1,1e2);

//  for(int i = 0; i < fold.size(); i++) {
//    fold[i] = i*4;
//  }

  IntegerVector subset = data[fold];

  IntegerVector data2 = seq(1,1e4);

  IntegerVector eq = ifelse(data >= data2, 1, 0);



  LogicalVector negMask(fold.size(), true);

  Rcout << "Length of negMask: " << negMask.size() << "and first two elements" << negMask[0] << " " << negMask[1] << "\n";

  eq = eq - 1; 

  //subset = normSam[fold];
}
