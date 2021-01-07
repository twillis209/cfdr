#include <Rcpp.h>
#include <Rmath.h>

using namespace Rcpp;

// [[Rcpp::export]]
NumericVector wrap_cummin() {
  NumericVector x = rnorm(1e3);
  NumericVector res = cummin(x);
  return res;
}
