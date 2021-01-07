#include <Rcpp.h>
#include <Rmath.h>

using namespace Rcpp;

/*

  I first got the idea of this implementation from a comment on a Reddit post. That code used `std::lower_bound`, which seemed not to work. I found the implementation in the following link when trying to fix this problem and switched to using `std::upper_bound` after viewing it (https://github.com/stepcie/sslcov/blob/master/src/ecdf.cpp). 

*/

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

/*
        ccut[i]=approx(xtest,cfsub-gx*xtest + gx*mx,zp[indices[i]],rule=2)$y

 TODO other rules; for now, we just have hard-coded a default of rule=2

 TODO just returning interpolated function values for now

 TODO if we knew there were no duplicates, we could speed things up by sorting xout, retaining the order, interpolating, then sorting the interpolated values back into the original order of xout

 TODO edge cases where xout lies outside extrema or xout is equal to an element of x

 Linear interpolation is performed with f(x)=f(x0)+[(f(x1)-f(x0))/(x1-x0)]*(x-x0)

 xtest is monotonic

 zp[indices[i]] is quantiles of a vector of p-values and not sorted

 rule=2 means that if an xout value falls outside [min(xtest), max(xtest)], then we use the closest of these two extrema
*/

// [[Rcpp::export]]
NumericVector approx_cpp(NumericVector x, NumericVector y, NumericVector xout) {
  NumericVector yout(xout.size());


  double x1_index, x1, x0, y1, y0;

  // For our purposes, x is strictly increasing
  for(int i = 0; i < xout.size(); i++) {
    // std::lower_bound finds the smallest element *not less than* xout[i]
    auto upper = std::lower_bound(x.begin(), x.end(), xout[i]);

    // TODO This can be optimised as well
    if(upper == x.begin()) {
      yout[i] = y[0];
      continue;
    } else if(upper == x.end()) {
      yout[i] = y[y.size()-1];
      continue;
    }

    // TODO not using this at the moment
    // Gets lower bound x0 (don't confuse this lower bound with the function lower_bound, which returns a value *not less than* xout[i])
    //auto lower = std::max_element(x.begin(), upper);

    x1_index = std::distance(x.begin(), upper);
    x1 = *upper;
    x0 = x[x1_index-1];
    y1 = y[x1_index];
    y0 = y[x1_index-1];

    yout[i] = y0+((y1-y0)/(x1-x0))*(xout[i]-x0);
  }
 
  return yout;

}

// TODO maybe we should correct the 1-indexed vectors at the outset

// [[Rcpp::export]]
List vl_mode2(NumericVector p, NumericVector q, IntegerVector indices, IntegerVector fold, bool adj=true, int nt=5000, int nv=1000, double p_threshold=0, CharacterVector scale=CharacterVector::create("p", "z"), bool closed=true, bool verbose=false, double gx=10e-5, Nullable<NumericVector> at=R_NilValue) {
//void vl_mode2(const NumericVector& p, const NumericVector& q, const IntegerVector& indices, const IntegerVector& fold, bool adj = TRUE, const NumericVector& at = R_NilValue, int nt = 5000, int nv = 1000, double p_threshold = 0, CharacterVector scale = CharacterVector({"p", "z"}), bool closed = TRUE, bool verbose = FALSE, double gx=10^-5) {


  NumericVector zp = qnorm(p/2);
  NumericVector zq = -1.0*qnorm(q/2);


  if (any(!is_finite(zp+zq))) {
    stop("P-values p,q must be in [1e-300,1]");
  }


  // maximum limits for integration
  double mx = max(NumericVector::create(10.0, max(abs(-qnorm(p/2)))));
  double my = max(NumericVector::create(10.0, max(abs(-qnorm(q/2)))));


  //gx=0 //1/1000  // use for 'tiebreaks'- if a point is on a curve with nonzero area, force the L-curve through that point

  NumericVector ccut;
  
  if (indices==R_NilValue) {
    if (at==R_NilValue) {
      stop("One of the parameters 'indices', 'at' must be set");
    }
    ccut=at;
    int mode=0;
  } else {
    // Note that the vector constructor is mimicking rep(0, length(indices))
    ccut = NumericVector(indices.size());
  }

  
  // yval2=seq(from=0,to=my,length.out=nv+1)[1:nv];
  NumericVector yval2(nv+1);

  for(int i = 0; i < yval2.size(); i++) {
    yval2[i] = my * (i/(yval2.size()-1)); 
  }
  
  yval2 = yval2[Range(0,nv-1)];

  // xval2=outer(rep(1,length(ccut)),yval2);
  NumericMatrix xval2(ccut.size(), yval2.size());

  for(int i = 0; i < ccut.size(); i++) {
    for(int j = 0; j < yval2.size(); j++) {
      xval2(i, j) = 1.0 * yval2[j]; 
    }
  }

  // pval2=2*pnorm(-yval2);
  NumericVector pval2 = 2*pnorm(-yval2);

  // xtest=seq(0,mx,length.out=nt);
  NumericVector xtest(nt);

  for(int i = 0; i < xtest.size(); i++) {
    xtest[i] = mx * (i/(xtest.size()-1)); 
  }

  // ptest=2*pnorm(-xtest);
  NumericVector ptest=2*pnorm(-xtest);

  // Omitting the mode 0 and 1 blocks from this section of code

  // warning: fold is a 1-based indexing mask and we are computing a 0-based mask
  LogicalVector negFoldMask(zq.size(), true);

  for(int j = 0; j < fold.size(); j++) {
    negFoldMask[(fold[j]-1)]=false;
  }

  NumericVector zp_negFold, zq_negFold, zq_indices;

  if(indices != R_NilValue) {

    ccut = NumericVector(indices.size());

    for(int i = 0; i < indices.size(); i++){
      // TODO slow not to preallocate the vector w here, but we don't know ahead of time; surely it would *not* be quicker to calculate the indices first?

      // The following code simulates the R code below
      //      w=which(zq[-fold] >= zq[indices[i]])


      zp_negFold = zp[negFoldMask];
      zq_negFold = zq[negFoldMask];
      zq_indices = zq[(indices[i]-1)];

      // 0-based mask
      IntegerVector w;

      for(int j = 0; j < fold.size(); j++) {
        if(zq_negFold[j] >= zq_indices[j]) {
          w.push_back(j);
        }
      }

      if(w.size() >=1 ) {
        // TODO incomplete line, needs to be finished to match the R line below
        NumericVector cfsub = (1+(1.0/w.size()))*ptest/(1+(1/w.size()))-ecdf_cpp(zp_negFold[w], xtest);

        // Compiler does not like the following line and error message is completely opaque to me, so I resort to a more verbose construction
        //cfsub=cummin(cfsub);

        NumericVector cummin_cfsub=cummin(cfsub);

        // approx_cpp returns a singleton vector here as zp[(indices[i]-1)] is a scalar
        ccut[i]=approx_cpp(xtest,cummin_cfsub-gx*xtest + gx*mx,zp[(indices[i]-1)])[0];
      } else {
        ccut[i] = p[(indices[i]-1)];
      }
    }
  }

  
  ccut=ccut*(1+ 1e-6); // ccut + 1e-8 // prevent floating-point comparison errors
  
    NumericVector out(ccut.size());

    NumericVector correct, correct_ccut;
  
    if (adj) {

      NumericVector p_negFold = p[negFoldMask];
      NumericVector q_negFold = q[negFoldMask];

      // 0-based mask
      IntegerVector w;

      for(int j = 0; j < p_negFold.size(); j++) {
        if(p_negFold[j] > 0.5) {
          w.push_back(j);
        }
      }

      //correct=cummin((1+ecdf(q[-fold][which(p[-fold]>0.5)])(pval2)*length(p[-fold]))/
      NumericVector p_negFold_w = p_negFold[w];

      NumericVector correct = cummin((1+ecdf_cpp(q_negFold[p_negFold_w], pval2)*p_negFold.size()));

        if(indices != R_NilValue) correct_ccut=approx_cpp(pval2,correct,q[indices]);
      //cummin((1+ecdf(pc[which(p>0.5)])(pc)*length(p))/(1+ rank(pc)))
      } else {
        correct(pval2.size(), 1);
        correct_ccut(ccut.size(), 1);
      }

    if(indices != R_NilValue) ccut=ccut*correct_ccut;

    NumericVector zp_ind = zp[indices];
    NumericVector zq_ind = zq[indices];

    // TODO I'm not sure if nt/mx will behave as it should
    zp_ind = zp_ind*(nt/mx);
    zq_ind = zq_ind*(nv/my);

    zp_ind = ceiling(zp_ind);
    zq_ind = ceiling(zq_ind);

    zp_ind = pmax(1, pmin(zp_ind,nt));
    zq_ind = pmax(1, pmin(zq_ind,nv));

    for(int i = 0; i < yval2.size(); i++) {
      IntegerVector w;

      for(int j = 0; j < zq_negFold.size(); j++) {
        if(zq_negFold[j] > yval2[i]) {
          w.push_back(j);
        }
      }

      if(w.size() >=1) {
        
        NumericVector cfsub = (1+(1.0/w.size()))*ptest/(1+(1/w.size()))-ecdf_cpp(zp_negFold[w], xtest);

        // Compiler does not like the following line and error message is completely opaque to me, so I resort to a more verbose construction
        //cfsub=cummin(cfsub);

        NumericVector cummin_cfsub=cummin(cfsub);

        // TODO not sure if f=1 matters if method = 'linear'
        //        xval2.column(i) = approx_cpp((cfsub-gx*xtest + gx*mx)*correct[i],xtest,ccut,rule=2,f=1);
        xval2.column(i) = approx_cpp((cfsub-gx*xtest + gx*mx)*correct[i],xtest,ccut);
      } else {
        // TODO truncation due to integer division?
        xval2.column(i)=-qnorm((ccut/correct_ccut)/2.0);
      }
    }


    IntegerVector w;

    for(int j = 0; j < fold.size(); j++) {
      if(zq_negFold[j] >= zq_indices[j]) {
        w.push_back(j);
      }
    }

        
  
 // https://stackoverflow.com/questions/49026407/how-can-i-do-logical-operations-on-rcppnumericmatrix-using-a-sugar-manner

    int xval2_nrow = xval2.nrow();
    int xval2_ncol = xval2.ncol();

    // The probability distribution functions in R:: deal with scalars and those in Rcpp:: with vectors, apparently the R::qnorm function has no defaults so one needs to call it verbosely like this
    double comparator = -R::qnorm(p_threshold/2.0, 0.0, 1.0, true, false);

    LogicalMatrix xval2_w(xval2_nrow, xval2_ncol);
    for(int i = 0; i < xval2_nrow; ++i) {
      for(int j = 0; j < xval2_ncol; ++j) {
        //xval2[which(xval2> -qnorm(p_threshold/2))]=-qnorm(p_threshold/2)
        xval2(i,j) = xval2(i,j) > comparator ? comparator : xval2(i,j); 
      }
    }

    // TODO there must be some sugar for this copying operation
    if(closed) {
      // yval2=c(Inf,0,yval2,Inf,Inf)
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

      //xval2=cbind(Inf,Inf,xval2,xval2[,nv],Inf)
      NumericMatrix xval2_mod(xval2_nrow, xval2_ncol+4);
      NumericVector r_posinf_vec(xval2_nrow, R_PosInf);
      xval2_mod.column(0) = r_posinf_vec;
      xval2_mod.column(1) = r_posinf_vec;
      xval2_mod.column(xval2_ncol+2) = xval2.column(nv);
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

    // TODO might be able to use Rcpp::pnorm for these calls

    if(scale[1]=="p") {
      X = NumericMatrix(xval2.nrow(), xval2.ncol());

      for(int i = 0; i < xval2_nrow; ++i) {
        for(int j = 0; j < xval2_ncol; ++j) {
          X(i,j) = 2.0*R::pnorm(-abs(xval2(i,j)), 0.0,1.0,true,false);
        }
      }

      Y = NumericVector(yval2.size());

      for(int i = 0; i < yval2.size(); ++i) {
        Y[i] = 2.0*R::pnorm(-abs(yval2[i]), 0.0,1.0,true,false);
      }

    } else {
      X = xval2;
      Y = yval2;
    }

    return List::create(Named("x")=X, Named("y")=Y);
}

// TODO internal interpolation function, not used in mode 2
//double interpol()

// vl(p,q, indices=ind[[i]], fold=folds[[i]])
// TODO mode argument hasn't been dealt with
// TODO const declarations for arguments
// TODO noNA optimisation
// TODO all subsetting with 1-based vectors!
