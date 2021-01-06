#include <Rcpp.h>
#include <Rmath.h>

using namespace Rcpp;

// From https://github.com/dmbates/ecdfExample

// [[Rcpp::export]]
NumericVector ecdf_cpp(NumericVector reference, NumericVector sample) {
  NumericVector sortedRef = clone(reference);

  std::sort(sortedRef.begin(), sortedRef.end());

  NumericVector estimatedQuantiles(sample.size());

  for(int i = 0; i < sample.size(); ++i) {
    int j = 0;
    while(sortedRef[j] <= sample[i] && j < sortedRef.size()) ++j;
    estimatedQuantiles[i] =  j/((double) sortedRef.size());
  }

  return estimatedQuantiles;
}


/*



https://github.com/stepcie/sslcov/blob/master/src/ecdf.cpp

*/

// [[Rcpp::export]]
NumericVector ecdf_cpp2(NumericVector reference, NumericVector sample) {
  NumericVector sortedRef = clone(reference);

  std::sort(sortedRef.begin(), sortedRef.end());

  NumericVector estimatedQuantiles(sample.size());

  for(int i = 0; i < sample.size(); ++i) {
    estimatedQuantiles[i] = (std::upper_bound(sortedRef.begin(), sortedRef.end(), sample[i])-sortedRef.begin());
  }

  return estimatedQuantiles/((double) sortedRef.size());
}


// Method which works if there are no duplicates

// [[Rcpp::export]]
NumericVector ecdf_cpp_noDup(NumericVector reference, NumericVector sample) {
  NumericVector sortedRef = clone(reference);
  NumericVector sortedSam = clone(sample);

  std::sort(sortedRef.begin(), sortedRef.end());
  std::sort(sortedSam.begin(), sortedSam.end());

  IntegerVector ord = match(sortedSam, sample);

  NumericVector estimatedQuantiles(sample.size());

  for(int i = 0, j = 0; i < sortedSam.size(); i++) {
    while(sortedRef[j] <= sortedSam[i] && j < sortedRef.size()) ++j;
    estimatedQuantiles[ord[i]-1] =  j/((double) sortedRef.size());
  }

  return estimatedQuantiles;
}


/*

//[[Rcpp::export]]
IntegerVector cppcp(NumericVector samp, NumericVector ref, IntegerVector ord) {
  int nobs = samp.size();
  IntegerVector ans(nobs);
  for (int i = 0, j = 0; i < nobs; ++i) {
    int ind(ord[i] - 1); // C++ uses 0-based indices
    double ssampi(samp[ind]);
    while (ref[j] < ssampi && j < ref.size()) ++j;
    ans[ind] = j;     // j is the 1-based index of the lower bound
  }
  return ans;
}

*/
