#include <Rcpp.h>
#include <Rmath.h>

using namespace Rcpp;

//' Empirical cdf
//'
//' This function estimates an empirical cdf using the reference parameter and then 
//' evaluates the estimated empirical cdf at the points specified in the sample parameter. 
//' 
//' @param reference NumericVector
//' @param sample NumericVector
//' @return NumericVector of estimated quantiles of sample values
//' 
//' @author Tom Willis
// [[Rcpp::export]]
NumericVector ecdf_cpp(NumericVector reference, NumericVector sample) {
  NumericVector sortedRef = clone(reference);

  std::sort(sortedRef.begin(), sortedRef.end());

  NumericVector estimatedQuantiles(sample.size());

  for(int i = 0; i < sample.size(); ++i) {
    estimatedQuantiles[i] = (std::upper_bound(sortedRef.begin(), sortedRef.end(), sample[i])-sortedRef.begin());
  }

  return estimatedQuantiles/((double) sortedRef.size());
}

//' Linear interpolation function
//'
//' This function carries out linear interpolation, reproducing a subset of the behaviour of stats::approx in R. Linear interpolation is performed with f(x)=f(x0)+[(f(x1)-f(x0))/(x1-x0)]*(x-x0)
//'
//' x is assumed to be sorted; if sorted in descending order, the order of elements in x and y is reversed. If a value of xout lies outside the range [min(x), max(x)], the value is interpolated using the nearest extremum (this corresponds to the behaviour of stats::approx with method=2).
//'
//' @param x NumericVector of x coordinates of points to be interpolated
//' @param y NumericVector of y coordinates of points to be interpolated
//' @param xout NumericVector of points at which to interpolate
//'
//' @return NumericVector of interpolated values 
//' 
//' @author Tom Willis
// [[Rcpp::export]]
NumericVector approx_cpp(NumericVector x, NumericVector y, NumericVector xout) {
  NumericVector yout(xout.size());

  NumericVector xc = clone(x);
  NumericVector yc = clone(y);

  if(!std::is_sorted(x.begin(), x.end())) {
    if(!std::is_sorted(x.begin(), x.end(), std::greater<double>())) {
        stop("Input vector x is not sorted in ascending or descending order");
      } else {
      std::reverse(xc.begin(), xc.end());
      std::reverse(yc.begin(), yc.end());
      }
  } 

  double x1_index, x1, x0, y1, y0;

  // For our purposes, x is strictly increasing
  for(int i = 0; i < xout.size(); ++i) {
    // std::lower_bound finds the smallest element *not less than* xout[i]
    auto upper = std::lower_bound(xc.begin(), xc.end(), xout[i]);

    // TODO This can be optimised as well
    if(upper == xc.begin()) {
      yout[i] = yc[0];
      continue;
    } else if(upper == xc.end()) {
      yout[i] = yc[yc.size()-1];
      continue;
    }

    x1_index = std::distance(xc.begin(), upper);
    x1 = *upper;
    x0 = xc[x1_index-1];
    y1 = yc[x1_index];
    y0 = yc[x1_index-1];

    yout[i] = y0+((y1-y0)/(x1-x0))*(xout[i]-x0);
  }
 
  return yout;

}

//' Returns coordinates of L-regions for cFDR method. Drop-in replacement for the cfdr::vl function in the use case wherein indices and fold are supplied and mode=2.
//' 
//' @param p principal p-values
//' @param q conditional p-values
//' @param adj adjust cFDR values and hence curves L using estimate of Pr(H0|Pj<pj)
//' @param indices indices of points at which to compute v(L)
//' @param at cfdr cutoff/cutoffs. Defaults to null
//' @param fold indices of points to exclude when from calculation of L-curves 
//' @param p_threshold if H0 is to be rejected automatically whenever p<p_threshold, include this in all regions L
//' @param nt number of test points in x-direction, default 5000
//' @param nv resolution for constructing L-curves, default 1000
//' @param scale return curves on the p- or z- plane. Y values are equally spaced on the z-plane.
//' @param closed determines whether curves are closed polygons encircling regions L (closed=T), or lines indicating the rightmost border of regions L
//' @param verbose print progress 
//' @return list containing elements x, y. Assuming n curves are calculated (where n=length(indices) or length(at)) and closed=T, x is a matrix of dimension n x (4+nv), y ix a vector of length (4+nv).
//' 
//' @author Tom Willis
// [[Rcpp::export]]
List vl_mode2(NumericVector p, NumericVector q, IntegerVector indices, IntegerVector fold,  bool adj=true, Nullable<NumericVector> at=R_NilValue, int nt=5000, int nv=1000, double p_threshold=0, CharacterVector scale=CharacterVector::create("p", "z"), bool closed=true, bool verbose=false, double gx=0.00001) {

  NumericVector zp = -qnorm(p/2);

  NumericVector zq = -qnorm(q/2);

  if (any(!is_finite(zp+zq))) {
    stop("P-values p,q must be in [1e-300,1]");
  }

  // maximum limits for integration
  double mx = max(NumericVector::create(10.0, max(abs(-qnorm(p/2)))));
  double my = max(NumericVector::create(10.0, max(abs(-qnorm(q/2)))));

  NumericVector ccut = NumericVector(indices.size(),0.0);
  
  NumericVector yval2(nv+1);

  for(int i = 0; i < yval2.size(); i++) {
    yval2[i] = my * ((double) i/(yval2.size()-1)); 
  }
  
  yval2 = yval2[Range(0,nv-1)];

  NumericMatrix xval2(ccut.size(), yval2.size());

  for(int i = 0; i < ccut.size(); i++) {
    for(int j = 0; j < yval2.size(); j++) {
      xval2(i, j) = 1.0 * yval2[j]; 
    }
  }

  NumericVector pval2 = 2.0*pnorm(-yval2);

  NumericVector xtest(nt,0.0);

  for(int i = 0; i < nt; i++) {
    xtest[i] = mx * ((double) i/(nt-1)); 
  }

  NumericVector ptest=2.0*pnorm(-xtest);

  // warning: fold is a 1-based indexing mask and we are computing a 0-based mask
  LogicalVector negFoldMask(zq.size(), true);

  for(int j = 0; j < fold.size(); j++) {
    negFoldMask[(fold[j]-1)]=false;
   }

  NumericVector zp_negFold = zp[negFoldMask];
  NumericVector zq_negFold = zq[negFoldMask];
  double zp_index, zq_index;

  for(int i = 0; i < indices.size(); i++){

    zp_index = zp[(indices[i]-1)];
    zq_index = zq[(indices[i]-1)];

    LogicalVector zq_negFold_w_bool = zq_negFold >= zq_index;
    NumericVector zp_negFold_w = zp_negFold[zq_negFold_w_bool];
    int w_size =  zp_negFold_w.size();

    if(w_size >=1 ) {

      NumericVector cfsub = (1.0+(1.0/w_size))*ptest/(1.0+(1.0/w_size)-ecdf_cpp(zp_negFold_w, xtest));

      cfsub=NumericVector(cummin(cfsub));

      NumericVector y_in = cfsub-gx*xtest + gx*mx;

      ccut[i]=approx_cpp(xtest,y_in,NumericVector::create(zp_index))[0];
    } else {
      ccut[i] = p[(indices[i]-1)];
    }
  }

  ccut=ccut*(1.0+ 1e-6); // ccut + 1e-8 // prevent floating-point comparison errors

  NumericVector correct, correct_ccut;
  
    if (adj) {

      NumericVector p_negFold = p[negFoldMask];
      NumericVector q_negFold = q[negFoldMask];

      LogicalVector p_negFold_w_bool = p_negFold > 0.5;
      NumericVector p_negFold_w = p_negFold[p_negFold_w_bool];
      int w_size = p_negFold_w.size();

      correct = (1.0+ecdf_cpp(q_negFold[p_negFold_w_bool], pval2)*p_negFold.size())/(1.0+ecdf_cpp(q_negFold, pval2)*p_negFold.size());

      correct=NumericVector(cummin(correct));

      correct_ccut=approx_cpp(pval2,correct,q[(indices-1)]);

      } else {
        correct = NumericVector(pval2.size(), 1.0);
        correct_ccut = NumericVector(ccut.size(), 1.0);
      }

    ccut=ccut*correct_ccut;

    NumericVector zp_ind = zp[(indices-1)];
    NumericVector zq_ind = zq[(indices-1)];

    zp_ind = zp_ind*(nt/mx);
    zq_ind = zq_ind*(nv/my);

    zp_ind = ceiling(zp_ind);
    zq_ind = ceiling(zq_ind);

    zp_ind = pmax(1.0, pmin(zp_ind,nt));
    zq_ind = pmax(1.0, pmin(zq_ind,nv));

    for(int i = 0; i < yval2.size(); i++) {
      LogicalVector zq_negFold_yval2_bool = zq_negFold > yval2[i];
      NumericVector zp_negFold_yval2 = zp_negFold[zq_negFold_yval2_bool];
      int w_size = zp_negFold_yval2.size();

      if(w_size >=1) {
        
        NumericVector cfsub = (1.0+(1.0/w_size))*ptest/(1.0+(1.0/w_size)-ecdf_cpp(zp_negFold_yval2, xtest));

        cfsub=NumericVector(cummin(cfsub));

        xval2.column(i) = approx_cpp((cfsub-gx*xtest + gx*mx)*correct[i],xtest,ccut);
      } else {
        xval2.column(i)=-qnorm((ccut/correct_ccut)/2.0);
      }
    }


    int xval2_nrow = xval2.nrow();
    int xval2_ncol = xval2.ncol();

    double comparator = -R::qnorm(p_threshold/2.0, 0.0, 1.0, true, false);

    LogicalMatrix xval2_w(xval2_nrow, xval2_ncol);
    for(int i = 0; i < xval2_nrow; ++i) {
      for(int j = 0; j < xval2_ncol; ++j) {
        xval2(i,j) = xval2(i,j) > comparator ? comparator : xval2(i,j); 
      }
    }

    if(closed) {
      int yval2_length = yval2.size();

      NumericVector yval2_mod = NumericVector(4+yval2_length);

      yval2_mod[0] = R_PosInf;
      yval2_mod[1] = 0;
      yval2_mod[yval2_length+2] = R_PosInf;
      yval2_mod[yval2_length+3] = R_PosInf;

      for(int i = 0; i < yval2_length; i++) {
        yval2_mod[2+i] = yval2[i];
      }
      yval2 = yval2_mod;


      NumericMatrix xval2_mod(xval2_nrow, xval2_ncol+4);
      NumericVector r_posinf_vec(xval2_nrow, R_PosInf);
      xval2_mod.column(0) = r_posinf_vec;
      xval2_mod.column(1) = r_posinf_vec;

      xval2_mod.column(xval2_ncol+2) = xval2.column(nv-1);
      xval2_mod.column(xval2_ncol+3) = r_posinf_vec;

      for(int i = 0; i < xval2_nrow; ++i) {
        for(int j = 0; j < xval2_ncol; ++j) {
          xval2_mod(i,j+2) = xval2(i,j);
        }
      }

      xval2 = xval2_mod;

    }

    NumericMatrix X;
    NumericVector Y;

    if(scale[0]=="p") {
      X = NumericMatrix(xval2.nrow(), xval2.ncol());

      for(int i = 0; i < xval2_nrow; ++i) {
        X.row(i) = 2.0*Rcpp::pnorm(-abs(xval2.row(i)));
      }

      Y = 2.0*Rcpp::pnorm(-abs(yval2));
    } else {
      X = xval2;
      Y = yval2;
    }

    return List::create(Named("x")=X, Named("y")=Y);
}
