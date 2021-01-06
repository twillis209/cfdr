#include <Rcpp.h>
#include <Rmath.h>

using namespace Rcpp;

// From https://github.com/dmbates/ecdfExample

/*
 TODO The implementation from dmbates which I've copied below requires that we have the order of the elements in sample when they are sorted. I'm not yet sure how to get this in Rcpp, although we could write a function to get it. What about duplicates in sample, though? Should include that in the test case
*/

// [[Rcpp::export]]
NumericVector ecdf_cpp(NumericVector reference, NumericVector sample) {
  NumericVector sortedRef = clone(reference);
  NumericVector sortedSample = clone(sample);

  std::sort(sortedRef.begin(), sortedRef.end());
  std::sort(sortedSample.begin(), sortedSample.end());

  IntegerVector ord = match(sortedSample, sample);

  NumericVector estimatedQuantiles(sample.size());

  for(int i = 0, j = 0; i < sortedSample.size(); ++i) {
    while(sortedRef[j] <= sortedSample[i] && j < sortedRef.size()) ++j;
    estimatedQuantiles[ord[i]] = (j+1)/((double) reference.size());
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

// TODO maybe we should correct the 1-indexed vectors at the outset
// [[Rcpp::export]]
void vl_mode2(NumericVector p, NumericVector q, IntegerVector indices, IntegerVector fold, bool adj, NumericVector at, int nt, int nv, double p_threshold, CharacterVector scale, bool closed, bool verbose, double gx) {
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

  /*
  if (!is.null(indices)) { // set ccut. NOT equal to cfdr at the points; needs to be adjusted since an additional test point is used in determination of L
      ccut=rep(0,length(indices))
      for (i in 1:length(indices)) {
        w=which(zq[-fold] >= zq[indices[i]])
        if (length(w) >= 1) {
          cfsub= (1+(1/length(w)))*ptest/(1+(1/length(w))-ecdf(zp[-fold][w])(xtest)); cfsub=cummin(cfsub)
          ccut[i]=approx(xtest,cfsub-gx*xtest + gx*mx,zp[indices[i]],rule=2)$y
        } else ccut[i]=p[indices[i]]
      }
      //    ccut= p[indices]*sapply(indices,function(x) 
      //      (1+length(which(q[-fold] <= q[x])))/(1+length(which(p[-fold] <= p[x] & q[-fold] <= q[x]))))// set ccut
  }
  */

  if(indices != R_NilValue) {

    ccut = NumericVector(indices.size());

    for(int i = 0; i < indices.size(); i++){
      // TODO slow not to preallocate the vector w here, but we don't know ahead of time; surely it would *not* be quicker to calculate the indices first?

      // The following code simulates the R code below
      //      w=which(zq[-fold] >= zq[indices[i]])

      // warning: fold is a 1-based indexing mask and we are computing a 0-based mask
      LogicalVector negFoldMask(zq.size(), true);

      for(int j = 0; j < fold.size(); j++) {
        negFoldMask[(fold[j]-1)]=false;
      }

      NumericVector zp_negFold = zp[negFoldMask];
      NumericVector zq_negFold = zq[negFoldMask];
      NumericVector zq_indices = zq[(indices[i]-1)];

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

        //    cfsub= (1+(1/length(w)))*ptest/(1+(1/length(w))-ecdf(zp[-fold][w])(xtest));
        // cfsub=cummin(cfsub);
        // ccut[i]=approx(xtest,cfsub-gx*xtest + gx*mx,zp[indices[i]],rule=2)$y
      } else {
        ccut[i] = p[(indices[i]-1)];
      }
    }
  }

  /*
  ccut=ccut*(1+ 1e-6) // ccut + 1e-8 // prevent floating-point comparison errors
  
  out=rep(0,length(ccut))
  
    if (adj) {
      correct=cummin((1+ecdf(q[-fold][which(p[-fold]>0.5)])(pval2)*length(p[-fold]))/
          (1+ ecdf(q[-fold])(pval2)*length(p[-fold])))
      if (!is.null(indices)) correct_ccut=approx(pval2,correct,q[indices],rule=2)$y //cummin((1+ecdf(pc[which(p>0.5)])(pc)*length(p))/(1+ rank(pc)))
    } else {
      correct=rep(1,length(pval2)) // adjustment factor for pc[i]
      correct_ccut=rep(1,length(ccut))
    }
    if (!is.null(indices)) ccut=ccut*correct_ccut
    
    zp_ind=ceiling(zp[indices]*nt/mx); zp_ind=pmax(1,pmin(zp_ind,nt))
    zq_ind=ceiling(zq[indices]*nv/my); zq_ind=pmax(1,pmin(zq_ind,nv))
    
    for (i in 1:length(yval2)) {
      w=which(zq[-fold] > yval2[i])
      if (length(w) >= 1 ) {
        cfsub= (1+(1/length(w)))*ptest/(1+(1/length(w))-ecdf(zp[-fold][w])(xtest)); cfsub=cummin(cfsub)
        xval2[,i]=approx((cfsub-gx*xtest + gx*mx)*correct[i],xtest,ccut,rule=2,f=1)$y 
      }  else xval2[,i]=-qnorm((ccut/correct_ccut)/2)
    }
    
  
  xval2[which(xval2> -qnorm(p_threshold/2))]=-qnorm(p_threshold/2)
  if (closed) {
    yval2=c(Inf,0,yval2,Inf,Inf)
    xval2=cbind(Inf,Inf,xval2,xval2[,nv],Inf)
  }
  
  if (scale[1]=="p") {
    X=2*pnorm(-abs(xval2))
    Y=2*pnorm(-abs(yval2))
  } else {
    X=xval2
    Y=yval2
  }
  
  return(list(x=X,y=Y))  
*/
  //  List out(2);

  //  return out;
}

// TODO internal interpolation function, not used in mode 2
//double interpol()

// vl(p,q, indices=ind[[i]], fold=folds[[i]])
// TODO mode argument hasn't been dealt with
// TODO const declarations for arguments
// TODO noNA optimisation
// TODO all subsetting with 1-based vectors!
