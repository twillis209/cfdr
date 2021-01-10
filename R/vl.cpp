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

 TODO TODO the call to approx_cpp with x=pval2, where pval2 is sorted largest to smallest

 Linear interpolation is performed with f(x)=f(x0)+[(f(x1)-f(x0))/(x1-x0)]*(x-x0)

 rule=2 means that if an xout value falls outside [min(xtest), max(xtest)], then we use the closest of these two extrema
*/

// [[Rcpp::export]]
NumericVector approx_cpp(NumericVector x, NumericVector y, NumericVector xout) {
  NumericVector yout(xout.size());

  NumericVector xc = clone(x);
  NumericVector yc = clone(y);

  if(!std::is_sorted(x.begin(), x.end())) {
    if(!std::is_sorted(x.begin(), x.end(), std::greater<double>())) {
        stop("Input vector x is not sorted in ascending or descending order");
      } else {
      //warning("Reversing order as x is sorted in descending order");
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


// [[Rcpp::export]]
List vl_mode2_part1(NumericVector p, NumericVector q, IntegerVector indices, IntegerVector fold,  bool adj=true, Nullable<NumericVector> at=R_NilValue, int nt=5000, int nv=1000, double p_threshold=0, CharacterVector scale=CharacterVector::create("p", "z"), bool closed=true, bool verbose=false, double gx=0.00001) {

  NumericVector zp = -qnorm(p/2);

  NumericVector zq = -qnorm(q/2);

  if (any(!is_finite(zp+zq))) {
    stop("P-values p,q must be in [1e-300,1]");
  }

  // maximum limits for integration
  double mx = max(NumericVector::create(10.0, max(abs(-qnorm(p/2)))));
  double my = max(NumericVector::create(10.0, max(abs(-qnorm(q/2)))));

  //Rprintf("mx: %.5f my: %.5f\n", mx, my);
  
    // Note that the vector constructor is mimicking rep(0, length(indices))
  NumericVector ccut = NumericVector(indices.size(),0.0);
  
  // yval2=seq(from=0,to=my,length.out=nv+1)[1:nv];
  NumericVector yval2(nv+1);

  for(int i = 0; i < yval2.size(); i++) {
    yval2[i] = my * ((double) i/(yval2.size()-1)); 
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
  NumericVector xtest(nt,0.0);

  // nt == xtest.size()
  for(int i = 0; i < nt; i++) {
    xtest[i] = mx * ((double) i/(nt-1)); 
  }

  // // Rcout << "test" << std::endl;
  // ptest=2*pnorm(-xtest);
  NumericVector ptest=2*pnorm(-xtest);

  /*
   zp
   zq
   mx
   my
   ccut
   xval2
   pval2
   ptest

   */

  return List::create(Named("ptest")=ptest, Named("pval2")=pval2, Named("xval2")=xval2, Named("my")=my, Named("mx")=mx, Named("zp")=zp, Named("zq")=zq, Named("yval2")=yval2, Named("xtest")=xtest);
}


// [[Rcpp::export]]
List vl_mode2_part2(NumericVector p, NumericVector q, IntegerVector indices, IntegerVector fold,  bool adj=true, Nullable<NumericVector> at=R_NilValue, int nt=5000, int nv=1000, double p_threshold=0, CharacterVector scale=CharacterVector::create("p", "z"), bool closed=true, bool verbose=false, double gx=0.00001) {

  NumericVector zp = -qnorm(p/2);

  NumericVector zq = -qnorm(q/2);

  if (any(!is_finite(zp+zq))) {
    stop("P-values p,q must be in [1e-300,1]");
  }

  // maximum limits for integration
  double mx = max(NumericVector::create(10.0, max(abs(-qnorm(p/2)))));
  double my = max(NumericVector::create(10.0, max(abs(-qnorm(q/2)))));

  //Rprintf("mx: %.5f my: %.5f\n", mx, my);
  
    // Note that the vector constructor is mimicking rep(0, length(indices))
  NumericVector ccut = NumericVector(indices.size(),0.0);
  
  // yval2=seq(from=0,to=my,length.out=nv+1)[1:nv];
  NumericVector yval2(nv+1);

  for(int i = 0; i < yval2.size(); i++) {
    yval2[i] = my * ((double) i/(yval2.size()-1)); 
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
  NumericVector xtest(nt,0.0);

  // nt == xtest.size()
  for(int i = 0; i < nt; i++) {
    xtest[i] = mx * ((double) i/(nt-1)); 
  }

  // // Rcout << "test" << std::endl;
  // ptest=2*pnorm(-xtest);
  NumericVector ptest=2.0*pnorm(-xtest);

  // warning: fold is a 1-based indexing mask and we are computing a 0-based mask
  LogicalVector negFoldMask(zq.size(), true);

  for(int j = 0; j < fold.size(); j++) {
    negFoldMask[(fold[j]-1)]=false;
   }

  NumericVector zp_negFold = zp[negFoldMask];
  NumericVector zq_negFold = zq[negFoldMask];
  double zp_index, zq_index;
  NumericMatrix cfsub(indices.size(), nt);
  NumericMatrix ecdf_results(indices.size(), nt);

  // // Rcout << "test2" << std::endl;

  for(int i = 0; i < indices.size(); i++){
    // The following code simulates the R code below
    //      w=which(zq[-fold] >= zq[indices[i]])

    zp_index = zp[(indices[i]-1)];
    zq_index = zq[(indices[i]-1)];

    LogicalVector zq_negFold_w_bool = zq_negFold >= zq_index;
    NumericVector zp_negFold_w = zp_negFold[zq_negFold_w_bool];
    int w_size =  zp_negFold_w.size();

    //Rprintf("w_size: %d \n", w_size);

    // // Rcout << "test3b" << std::endl;
    if(w_size >=1 ) {
      // TODO incomplete line, needs to be finished to match the R line below
      //      cfsub.row(i) = (1.0+(1.0/w_size))*ptest/(1.0+(1.0/w_size))-ecdf_cpp(zp_negFold_w, xtest);

      cfsub.row(i) = (1.0+(1.0/w_size))*ptest/(1.0+(1.0/w_size)-ecdf_cpp(zp_negFold_w, xtest));

      ecdf_results.row(i) = ecdf_cpp(zp_negFold_w, xtest);

      // // Rcout << "test3c" << std::endl;
      cfsub.row(i)=NumericVector(cummin(cfsub.row(i)));

      NumericVector y_in = cfsub.row(i)-gx*xtest + gx*mx;

      // // Rcout << "y_in.size(): " << y_in.size() << std::endl;

      // // Rcout << "zp[(indices[i]-1)]: " << zp[(indices[i]-1)] << std::endl;
      // // Rcout << "test3d" << std::endl;

      // TODO It must be tremendously inefficient to keep creating a vector here
      // approx_cpp returns a singleton vector here as zp[(indices[i]-1)] is a scalar
      ccut[i]=approx_cpp(xtest,y_in,NumericVector::create(zp_index))[0];

      // // Rcout << "test3e" << std::endl;
    } else {

      // // Rcout << "test3f" << std::endl;
      ccut[i] = p[(indices[i]-1)];

      // // Rcout << "test3g" << std::endl;
    }
  }

  return List::create(Named("negFold")=negFoldMask, Named("zp_negFold")=zp_negFold, Named("ccut")=ccut, Named("cfsub")=cfsub, Named("ecdf_results")=ecdf_results);
}

// [[Rcpp::export]]
List vl_mode2_part3(NumericVector p, NumericVector q, IntegerVector indices, IntegerVector fold,  bool adj=true, Nullable<NumericVector> at=R_NilValue, int nt=5000, int nv=1000, double p_threshold=0, CharacterVector scale=CharacterVector::create("p", "z"), bool closed=true, bool verbose=false, double gx=0.00001) {

  NumericVector zp = -qnorm(p/2);

  NumericVector zq = -qnorm(q/2);

  if (any(!is_finite(zp+zq))) {
    stop("P-values p,q must be in [1e-300,1]");
  }

  // maximum limits for integration
  double mx = max(NumericVector::create(10.0, max(abs(-qnorm(p/2)))));
  double my = max(NumericVector::create(10.0, max(abs(-qnorm(q/2)))));

  //Rprintf("mx: %.5f my: %.5f\n", mx, my);
  
    // Note that the vector constructor is mimicking rep(0, length(indices))
  NumericVector ccut = NumericVector(indices.size(),0.0);
  
  // yval2=seq(from=0,to=my,length.out=nv+1)[1:nv];
  NumericVector yval2(nv+1);

  for(int i = 0; i < yval2.size(); i++) {
    yval2[i] = my * ((double) i/(yval2.size()-1)); 
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
  NumericVector pval2 = 2.0*pnorm(-yval2);

  // xtest=seq(0,mx,length.out=nt);
  NumericVector xtest(nt,0.0);

  // nt == xtest.size()
  for(int i = 0; i < nt; i++) {
    xtest[i] = mx * ((double) i/(nt-1)); 
  }

  // // Rcout << "test" << std::endl;
  // ptest=2*pnorm(-xtest);
  NumericVector ptest=2.0*pnorm(-xtest);

  // warning: fold is a 1-based indexing mask and we are computing a 0-based mask
  LogicalVector negFoldMask(zq.size(), true);

  for(int j = 0; j < fold.size(); j++) {
    negFoldMask[(fold[j]-1)]=false;
   }

  NumericVector zp_negFold = zp[negFoldMask];
  NumericVector zq_negFold = zq[negFoldMask];
  double zp_index, zq_index;

  // // Rcout << "test2" << std::endl;

  for(int i = 0; i < indices.size(); i++){
    // The following code simulates the R code below
    //      w=which(zq[-fold] >= zq[indices[i]])

    zp_index = zp[(indices[i]-1)];
    zq_index = zq[(indices[i]-1)];

    LogicalVector zq_negFold_w_bool = zq_negFold >= zq_index;
    NumericVector zp_negFold_w = zp_negFold[zq_negFold_w_bool];
    int w_size =  zp_negFold_w.size();

    //Rprintf("w_size: %d \n", w_size);

    // // Rcout << "test3b" << std::endl;
    if(w_size >=1 ) {
      // TODO incomplete line, needs to be finished to match the R line below
      //      cfsub.row(i) = (1.0+(1.0/w_size))*ptest/(1.0+(1.0/w_size))-ecdf_cpp(zp_negFold_w, xtest);

      NumericVector cfsub = (1.0+(1.0/w_size))*ptest/(1.0+(1.0/w_size)-ecdf_cpp(zp_negFold_w, xtest));

      // // Rcout << "test3c" << std::endl;
      cfsub=NumericVector(cummin(cfsub));

      NumericVector y_in = cfsub-gx*xtest + gx*mx;

      // // Rcout << "y_in.size(): " << y_in.size() << std::endl;

      // // Rcout << "zp[(indices[i]-1)]: " << zp[(indices[i]-1)] << std::endl;
      // // Rcout << "test3d" << std::endl;

      // TODO It must be tremendously inefficient to keep creating a vector here
      // approx_cpp returns a singleton vector here as zp[(indices[i]-1)] is a scalar
      ccut[i]=approx_cpp(xtest,y_in,NumericVector::create(zp_index))[0];

      // // Rcout << "test3e" << std::endl;
    } else {

      // // Rcout << "test3f" << std::endl;
      ccut[i] = p[(indices[i]-1)];

      // // Rcout << "test3g" << std::endl;
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

      // TODO explain why we wrap here, for now, just note that not doing this leads to a compiler error
      correct=NumericVector(cummin(correct));

      correct_ccut=approx_cpp(pval2,correct,q[(indices-1)]);

      } else {
        correct = NumericVector(pval2.size(), 1.0);
        correct_ccut = NumericVector(ccut.size(), 1.0);
      }

    // Rcout << "test5" << std::endl;

    ccut=ccut*correct_ccut;

    /*
    NumericVector zp_ind = zp[(indices-1)];
    NumericVector zq_ind = zq[(indices-1)];

    // TODO I'm not sure if nt/mx will behave as it should
    zp_ind = zp_ind*(nt/mx);
    zq_ind = zq_ind*(nv/my);

    zp_ind = ceiling(zp_ind);
    zq_ind = ceiling(zq_ind);

    zp_ind = pmax(1, pmin(zp_ind,nt));
    zq_ind = pmax(1, pmin(zq_ind,nv));

    for(int i = 0; i < yval2.size(); i++) {
      LogicalVector zq_negFold_yval2_bool = zq_negFold > yval2[i];
      NumericVector zp_negFold_yval2 = zp_negFold[zq_negFold_yval2_bool];
      int w_size = zp_negFold_yval2.size();

      if(w_size >=1) {
        
        NumericVector cfsub = (1+(1.0/w_size))*ptest/(1+(1/w_size))-ecdf_cpp(zp_negFold_yval2, xtest);

        cfsub=NumericVector(cummin(cfsub));

        // TODO not sure if f=1 matters if method = 'linear'
        //        xval2.column(i) = approx_cpp((cfsub-gx*xtest + gx*mx)*correct[i],xtest,ccut,rule=2,f=1);
        xval2.column(i) = approx_cpp((cfsub-gx*xtest + gx*mx)*correct[i],xtest,ccut);
      } else {
        // TODO truncation due to integer division?
        xval2.column(i)=-qnorm((ccut/correct_ccut)/2.0);
      }
    }

    // Rcout << "test6" << std::endl;
  
 // https://stackoverflow.com/questions/49026407/how-can-i-do-logical-operations-on-rcppnumericmatrix-using-a-sugar-manner

    int xval2_nrow = xval2.nrow();
    int xval2_ncol = xval2.ncol();

    // The probability distribution functions in R:: deal with scalars and those in Rcpp:: with vectors, apparently the R::qnorm function has no defaults so one needs to call it verbosely like this
    double comparator = -R::qnorm(p_threshold/2.0, 0.0, 1.0, true, false);

    LogicalMatrix xval2_w(xval2_nrow, xval2_ncol);
    for(int i = 0; i < xval2_nrow; ++i) {
      for(int j = 0; j < xval2_ncol; ++j) {
        xval2(i,j) = xval2(i,j) > comparator ? comparator : xval2(i,j); 
      }
    }

    // Rcout << "test7" << std::endl;

    // TODO there must be some sugar for this copying operation
    if(closed) {
      // yval2=c(Inf,0,yval2,Inf,Inf)
      int yval2_length = yval2.size();
      NumericVector yval2_mod = NumericVector(4+yval2_length);
      yval2_mod[0] = R_PosInf;
      yval2_mod[1] = 0;
      yval2_mod[yval2_length+2] = R_PosInf;
      yval2_mod[yval2_length+3] = R_PosInf;

      // Rcout << "test7a" << std::endl;

      for(int i = 0; i < yval2_length; i++) {
        yval2_mod[2+i] = yval2[i];
      }
      yval2 = yval2_mod;


      // Rcout << "test7b" << std::endl;

      //xval2=cbind(Inf,Inf,xval2,xval2[,nv],Inf)
      NumericMatrix xval2_mod(xval2_nrow, xval2_ncol+4);
      NumericVector r_posinf_vec(xval2_nrow, R_PosInf);
      xval2_mod.column(0) = r_posinf_vec;
      xval2_mod.column(1) = r_posinf_vec;

      // Rcout << "test7bi" << std::endl;
      xval2_mod.column(xval2_ncol+2) = xval2.column(nv-1);
      xval2_mod.column(xval2_ncol+3) = r_posinf_vec;

      // Rcout << "test7c" << std::endl;

      for(int i = 0; i < xval2_nrow; ++i) {
        for(int j = 0; j < xval2_ncol; ++j) {
          xval2_mod(i,j+2) = xval2(i,j);
        }
      }

      xval2 = xval2_mod;

    }

    // Rcout << "test8" << std::endl;

    NumericMatrix X;
    NumericVector Y;

    // TODO might be able to use Rcpp::pnorm for these calls

    if(scale[0]=="p") {
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

    // Rcout << "test9" << std::endl;

    return List::create(Named("x")=X, Named("y")=Y);
    */
    return List::create(Named("correct")=correct, Named("ccut")=ccut, Named("correct_ccut")=correct_ccut, Named("pval2")=pval2, Named("q[indices-1]")=q[(indices-1)]);
}

// [[Rcpp::export]]
List vl_mode2_part4(NumericVector p, NumericVector q, IntegerVector indices, IntegerVector fold,  bool adj=true, Nullable<NumericVector> at=R_NilValue, int nt=5000, int nv=1000, double p_threshold=0, CharacterVector scale=CharacterVector::create("p", "z"), bool closed=true, bool verbose=false, double gx=0.00001) {

  NumericVector zp = -qnorm(p/2);

  NumericVector zq = -qnorm(q/2);

  if (any(!is_finite(zp+zq))) {
    stop("P-values p,q must be in [1e-300,1]");
  }

  // maximum limits for integration
  double mx = max(NumericVector::create(10.0, max(abs(-qnorm(p/2)))));
  double my = max(NumericVector::create(10.0, max(abs(-qnorm(q/2)))));

  //Rprintf("mx: %.5f my: %.5f\n", mx, my);
  
    // Note that the vector constructor is mimicking rep(0, length(indices))
  NumericVector ccut = NumericVector(indices.size(),0.0);
  
  // yval2=seq(from=0,to=my,length.out=nv+1)[1:nv];
  NumericVector yval2(nv+1);

  for(int i = 0; i < yval2.size(); i++) {
    yval2[i] = my * ((double) i/(yval2.size()-1)); 
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
  NumericVector pval2 = 2.0*pnorm(-yval2);

  // xtest=seq(0,mx,length.out=nt);
  NumericVector xtest(nt,0.0);

  // nt == xtest.size()
  for(int i = 0; i < nt; i++) {
    xtest[i] = mx * ((double) i/(nt-1)); 
  }

  // // Rcout << "test" << std::endl;
  // ptest=2*pnorm(-xtest);
  NumericVector ptest=2.0*pnorm(-xtest);

  // warning: fold is a 1-based indexing mask and we are computing a 0-based mask
  LogicalVector negFoldMask(zq.size(), true);

  for(int j = 0; j < fold.size(); j++) {
    negFoldMask[(fold[j]-1)]=false;
   }

  NumericVector zp_negFold = zp[negFoldMask];
  NumericVector zq_negFold = zq[negFoldMask];
  double zp_index, zq_index;

  // // Rcout << "test2" << std::endl;

  for(int i = 0; i < indices.size(); i++){
    // The following code simulates the R code below
    //      w=which(zq[-fold] >= zq[indices[i]])

    zp_index = zp[(indices[i]-1)];
    zq_index = zq[(indices[i]-1)];

    LogicalVector zq_negFold_w_bool = zq_negFold >= zq_index;
    NumericVector zp_negFold_w = zp_negFold[zq_negFold_w_bool];
    int w_size =  zp_negFold_w.size();

    //Rprintf("w_size: %d \n", w_size);

    // // Rcout << "test3b" << std::endl;
    if(w_size >=1 ) {
      // TODO incomplete line, needs to be finished to match the R line below
      //      cfsub.row(i) = (1.0+(1.0/w_size))*ptest/(1.0+(1.0/w_size))-ecdf_cpp(zp_negFold_w, xtest);

      NumericVector cfsub = (1.0+(1.0/w_size))*ptest/(1.0+(1.0/w_size)-ecdf_cpp(zp_negFold_w, xtest));

      // // Rcout << "test3c" << std::endl;
      cfsub=NumericVector(cummin(cfsub));

      NumericVector y_in = cfsub-gx*xtest + gx*mx;

      // // Rcout << "y_in.size(): " << y_in.size() << std::endl;

      // // Rcout << "zp[(indices[i]-1)]: " << zp[(indices[i]-1)] << std::endl;
      // // Rcout << "test3d" << std::endl;

      // TODO It must be tremendously inefficient to keep creating a vector here
      // approx_cpp returns a singleton vector here as zp[(indices[i]-1)] is a scalar
      ccut[i]=approx_cpp(xtest,y_in,NumericVector::create(zp_index))[0];

      // // Rcout << "test3e" << std::endl;
    } else {

      // // Rcout << "test3f" << std::endl;
      ccut[i] = p[(indices[i]-1)];

      // // Rcout << "test3g" << std::endl;
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

      // TODO explain why we wrap here, for now, just note that not doing this leads to a compiler error
      correct=NumericVector(cummin(correct));

      correct_ccut=approx_cpp(pval2,correct,q[(indices-1)]);

      } else {
        correct = NumericVector(pval2.size(), 1.0);
        correct_ccut = NumericVector(ccut.size(), 1.0);
      }

    // Rcout << "test5" << std::endl;

    ccut=ccut*correct_ccut;

    // Testing the code from this point onwards in part 4

    NumericVector zp_ind = zp[(indices-1)];
    NumericVector zq_ind = zq[(indices-1)];

    // TODO I'm not sure if nt/mx will behave as it should
    zp_ind = zp_ind*(nt/mx);
    zq_ind = zq_ind*(nv/my);

    zp_ind = ceiling(zp_ind);
    zq_ind = ceiling(zq_ind);

    zp_ind = pmax(1.0, pmin(zp_ind,nt));
    zq_ind = pmax(1.0, pmin(zq_ind,nv));

    NumericMatrix cfsub_mat(yval2.size(), xtest.size());
    NumericMatrix ecdf_mat(yval2.size(), xtest.size());

    for(int i = 0; i < yval2.size(); i++) {
      LogicalVector zq_negFold_yval2_bool = zq_negFold > yval2[i];
      NumericVector zp_negFold_yval2 = zp_negFold[zq_negFold_yval2_bool];
      int w_size = zp_negFold_yval2.size();

      if(w_size >=1) {
        
        NumericVector cfsub = (1.0+(1.0/w_size))*ptest/(1.0+(1.0/w_size)-ecdf_cpp(zp_negFold_yval2, xtest));

        ecdf_mat.row(i) = ecdf_cpp(zp_negFold_yval2, xtest);

        cfsub=NumericVector(cummin(cfsub));

        cfsub_mat.row(i) = cfsub;

        xval2.column(i) = approx_cpp((cfsub-gx*xtest + gx*mx)*correct[i],xtest,ccut);
      } else {
        // TODO truncation due to integer division?
        xval2.column(i)=-qnorm((ccut/correct_ccut)/2.0);
      }
    }

    return(List::create(Named("xval2")=xval2, Named("zp_ind")=zp_ind, Named("zq_ind")=zq_ind, Named("cfsub_mat")=cfsub_mat, Named("ecdf_mat")=ecdf_mat));

      /*

    // Rcout << "test6" << std::endl;
  
 // https://stackoverflow.com/questions/49026407/how-can-i-do-logical-operations-on-rcppnumericmatrix-using-a-sugar-manner

    int xval2_nrow = xval2.nrow();
    int xval2_ncol = xval2.ncol();

    // The probability distribution functions in R:: deal with scalars and those in Rcpp:: with vectors, apparently the R::qnorm function has no defaults so one needs to call it verbosely like this
    double comparator = -R::qnorm(p_threshold/2.0, 0.0, 1.0, true, false);

    LogicalMatrix xval2_w(xval2_nrow, xval2_ncol);
    for(int i = 0; i < xval2_nrow; ++i) {
      for(int j = 0; j < xval2_ncol; ++j) {
        xval2(i,j) = xval2(i,j) > comparator ? comparator : xval2(i,j); 
      }
    }

    // Rcout << "test7" << std::endl;

    // TODO there must be some sugar for this copying operation
    if(closed) {
      // yval2=c(Inf,0,yval2,Inf,Inf)
      int yval2_length = yval2.size();
      NumericVector yval2_mod = NumericVector(4+yval2_length);
      yval2_mod[0] = R_PosInf;
      yval2_mod[1] = 0;
      yval2_mod[yval2_length+2] = R_PosInf;
      yval2_mod[yval2_length+3] = R_PosInf;

      // Rcout << "test7a" << std::endl;

      for(int i = 0; i < yval2_length; i++) {
        yval2_mod[2+i] = yval2[i];
      }
      yval2 = yval2_mod;


      // Rcout << "test7b" << std::endl;

      //xval2=cbind(Inf,Inf,xval2,xval2[,nv],Inf)
      NumericMatrix xval2_mod(xval2_nrow, xval2_ncol+4);
      NumericVector r_posinf_vec(xval2_nrow, R_PosInf);
      xval2_mod.column(0) = r_posinf_vec;
      xval2_mod.column(1) = r_posinf_vec;

      // Rcout << "test7bi" << std::endl;
      xval2_mod.column(xval2_ncol+2) = xval2.column(nv-1);
      xval2_mod.column(xval2_ncol+3) = r_posinf_vec;

      // Rcout << "test7c" << std::endl;

      for(int i = 0; i < xval2_nrow; ++i) {
        for(int j = 0; j < xval2_ncol; ++j) {
          xval2_mod(i,j+2) = xval2(i,j);
        }
      }

      xval2 = xval2_mod;

    }

    // Rcout << "test8" << std::endl;

    NumericMatrix X;
    NumericVector Y;

    // TODO might be able to use Rcpp::pnorm for these calls

    if(scale[0]=="p") {
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
      */
}

// [[Rcpp::export]]
List vl_mode2_part5(NumericVector p, NumericVector q, IntegerVector indices, IntegerVector fold,  bool adj=true, Nullable<NumericVector> at=R_NilValue, int nt=5000, int nv=1000, double p_threshold=0, CharacterVector scale=CharacterVector::create("p", "z"), bool closed=true, bool verbose=false, double gx=0.00001) {

  NumericVector zp = -qnorm(p/2);

  NumericVector zq = -qnorm(q/2);

  if (any(!is_finite(zp+zq))) {
    stop("P-values p,q must be in [1e-300,1]");
  }

  // maximum limits for integration
  double mx = max(NumericVector::create(10.0, max(abs(-qnorm(p/2)))));
  double my = max(NumericVector::create(10.0, max(abs(-qnorm(q/2)))));

  //Rprintf("mx: %.5f my: %.5f\n", mx, my);
  
    // Note that the vector constructor is mimicking rep(0, length(indices))
  NumericVector ccut = NumericVector(indices.size(),0.0);
  
  // yval2=seq(from=0,to=my,length.out=nv+1)[1:nv];
  NumericVector yval2(nv+1);

  for(int i = 0; i < yval2.size(); i++) {
    yval2[i] = my * ((double) i/(yval2.size()-1)); 
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
  NumericVector pval2 = 2.0*pnorm(-yval2);

  // xtest=seq(0,mx,length.out=nt);
  NumericVector xtest(nt,0.0);

  // nt == xtest.size()
  for(int i = 0; i < nt; i++) {
    xtest[i] = mx * ((double) i/(nt-1)); 
  }

  // // Rcout << "test" << std::endl;
  // ptest=2*pnorm(-xtest);
  NumericVector ptest=2.0*pnorm(-xtest);

  // warning: fold is a 1-based indexing mask and we are computing a 0-based mask
  LogicalVector negFoldMask(zq.size(), true);

  for(int j = 0; j < fold.size(); j++) {
    negFoldMask[(fold[j]-1)]=false;
   }

  NumericVector zp_negFold = zp[negFoldMask];
  NumericVector zq_negFold = zq[negFoldMask];
  double zp_index, zq_index;

  // // Rcout << "test2" << std::endl;

  for(int i = 0; i < indices.size(); i++){
    // The following code simulates the R code below
    //      w=which(zq[-fold] >= zq[indices[i]])

    zp_index = zp[(indices[i]-1)];
    zq_index = zq[(indices[i]-1)];

    LogicalVector zq_negFold_w_bool = zq_negFold >= zq_index;
    NumericVector zp_negFold_w = zp_negFold[zq_negFold_w_bool];
    int w_size =  zp_negFold_w.size();

    //Rprintf("w_size: %d \n", w_size);

    // // Rcout << "test3b" << std::endl;
    if(w_size >=1 ) {
      // TODO incomplete line, needs to be finished to match the R line below
      //      cfsub.row(i) = (1.0+(1.0/w_size))*ptest/(1.0+(1.0/w_size))-ecdf_cpp(zp_negFold_w, xtest);

      NumericVector cfsub = (1.0+(1.0/w_size))*ptest/(1.0+(1.0/w_size)-ecdf_cpp(zp_negFold_w, xtest));

      // // Rcout << "test3c" << std::endl;
      cfsub=NumericVector(cummin(cfsub));

      NumericVector y_in = cfsub-gx*xtest + gx*mx;

      // // Rcout << "y_in.size(): " << y_in.size() << std::endl;

      // // Rcout << "zp[(indices[i]-1)]: " << zp[(indices[i]-1)] << std::endl;
      // // Rcout << "test3d" << std::endl;

      // TODO It must be tremendously inefficient to keep creating a vector here
      // approx_cpp returns a singleton vector here as zp[(indices[i]-1)] is a scalar
      ccut[i]=approx_cpp(xtest,y_in,NumericVector::create(zp_index))[0];

      // // Rcout << "test3e" << std::endl;
    } else {

      // // Rcout << "test3f" << std::endl;
      ccut[i] = p[(indices[i]-1)];

      // // Rcout << "test3g" << std::endl;
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

      // TODO explain why we wrap here, for now, just note that not doing this leads to a compiler error
      correct=NumericVector(cummin(correct));

      correct_ccut=approx_cpp(pval2,correct,q[(indices-1)]);

      } else {
        correct = NumericVector(pval2.size(), 1.0);
        correct_ccut = NumericVector(ccut.size(), 1.0);
      }

    // Rcout << "test5" << std::endl;

    ccut=ccut*correct_ccut;

    // Testing the code from this point onwards in part 4

    NumericVector zp_ind = zp[(indices-1)];
    NumericVector zq_ind = zq[(indices-1)];

    // TODO I'm not sure if nt/mx will behave as it should
    zp_ind = zp_ind*(nt/mx);
    zq_ind = zq_ind*(nv/my);

    zp_ind = ceiling(zp_ind);
    zq_ind = ceiling(zq_ind);

    zp_ind = pmax(1.0, pmin(zp_ind,nt));
    zq_ind = pmax(1.0, pmin(zq_ind,nv));

    NumericMatrix cfsub_mat(yval2.size(), xtest.size());
    NumericMatrix ecdf_mat(yval2.size(), xtest.size());

    for(int i = 0; i < yval2.size(); i++) {
      LogicalVector zq_negFold_yval2_bool = zq_negFold > yval2[i];
      NumericVector zp_negFold_yval2 = zp_negFold[zq_negFold_yval2_bool];
      int w_size = zp_negFold_yval2.size();

      if(w_size >=1) {
        
        NumericVector cfsub = (1.0+(1.0/w_size))*ptest/(1.0+(1.0/w_size)-ecdf_cpp(zp_negFold_yval2, xtest));

        ecdf_mat.row(i) = ecdf_cpp(zp_negFold_yval2, xtest);

        cfsub=NumericVector(cummin(cfsub));

        cfsub_mat.row(i) = cfsub;

        xval2.column(i) = approx_cpp((cfsub-gx*xtest + gx*mx)*correct[i],xtest,ccut);
      } else {
        // TODO truncation due to integer division?
        xval2.column(i)=-qnorm((ccut/correct_ccut)/2.0);
      }
    }

    // Rcout << "test6" << std::endl;
  
 // https://stackoverflow.com/questions/49026407/how-can-i-do-logical-operations-on-rcppnumericmatrix-using-a-sugar-manner

    int xval2_nrow = xval2.nrow();
    int xval2_ncol = xval2.ncol();

    // The probability distribution functions in R:: deal with scalars and those in Rcpp:: with vectors, apparently the R::qnorm function has no defaults so one needs to call it verbosely like this
    double comparator = -R::qnorm(p_threshold/2.0, 0.0, 1.0, true, false);

    LogicalMatrix xval2_w(xval2_nrow, xval2_ncol);
    for(int i = 0; i < xval2_nrow; ++i) {
      for(int j = 0; j < xval2_ncol; ++j) {
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

      // Rcout << "test7a" << std::endl;

      for(int i = 0; i < yval2_length; i++) {
        yval2_mod[2+i] = yval2[i];
      }
      yval2 = yval2_mod;


      // Rcout << "test7b" << std::endl;

      //xval2=cbind(Inf,Inf,xval2,xval2[,nv],Inf)
      NumericMatrix xval2_mod(xval2_nrow, xval2_ncol+4);
      NumericVector r_posinf_vec(xval2_nrow, R_PosInf);
      xval2_mod.column(0) = r_posinf_vec;
      xval2_mod.column(1) = r_posinf_vec;

      // Rcout << "test7bi" << std::endl;
      xval2_mod.column(xval2_ncol+2) = xval2.column(nv-1);
      xval2_mod.column(xval2_ncol+3) = r_posinf_vec;

      // Rcout << "test7c" << std::endl;

      for(int i = 0; i < xval2_nrow; ++i) {
        for(int j = 0; j < xval2_ncol; ++j) {
          xval2_mod(i,j+2) = xval2(i,j);
        }
      }

      xval2 = xval2_mod;

    }

    return(List::create(Named("xval2")=xval2, Named("yval2")=yval2));

    /*
    // Rcout << "test8" << std::endl;

    NumericMatrix X;
    NumericVector Y;

    // TODO might be able to use Rcpp::pnorm for these calls

    if(scale[0]=="p") {
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
    */
}

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

  //Rprintf("mx: %.5f my: %.5f\n", mx, my);
  
    // Note that the vector constructor is mimicking rep(0, length(indices))
  NumericVector ccut = NumericVector(indices.size(),0.0);
  
  // yval2=seq(from=0,to=my,length.out=nv+1)[1:nv];
  NumericVector yval2(nv+1);

  for(int i = 0; i < yval2.size(); i++) {
    yval2[i] = my * ((double) i/(yval2.size()-1)); 
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
  NumericVector pval2 = 2.0*pnorm(-yval2);

  // xtest=seq(0,mx,length.out=nt);
  NumericVector xtest(nt,0.0);

  // nt == xtest.size()
  for(int i = 0; i < nt; i++) {
    xtest[i] = mx * ((double) i/(nt-1)); 
  }

  // ptest=2*pnorm(-xtest);
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

      /*
        TODO

        zp_index is a single scalar here, so we can probably replace this call to approx_cpp with something more explicit. Copying zp_index to a singleton vector can't be cheap.

        Is xtest always sorted in decreasing order?

       */
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

    // The probability distribution functions in R:: deal with scalars and those in Rcpp:: with vectors, apparently the R::qnorm function has no defaults so one needs to call it verbosely like this
    double comparator = -R::qnorm(p_threshold/2.0, 0.0, 1.0, true, false);

    LogicalMatrix xval2_w(xval2_nrow, xval2_ncol);
    for(int i = 0; i < xval2_nrow; ++i) {
      for(int j = 0; j < xval2_ncol; ++j) {
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

      // Rcout << "test7bi" << std::endl;
      xval2_mod.column(xval2_ncol+2) = xval2.column(nv-1);
      xval2_mod.column(xval2_ncol+3) = r_posinf_vec;

      // Rcout << "test7c" << std::endl;

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
