#include <Rcpp.h>
#include <Rmath.h>
using namespace Rcpp;

// TODO internal interpolation function, not used in mode 2
//double interpol()

// vl(p,q, indices=ind[[i]], fold=folds[[i]])
// TODO mode argument hasn't been dealt with
// TODO const declarations for arguments
// TODO type of at 
// [[Rcpp::export]]
List vl_mode2(const NumericVector& p, const NumericVector& q, const NumericVector& indices, const NumericVector& fold, bool adj = TRUE, const NumericVector& at = R_NilValue, int nt = 5000, int nv = 1000, double p_threshold = 0, CharacterVector scale = CharacterVector({"p", "z"}), bool closed = TRUE, bool verbose = FALSE, double gx=10^-5) {

  NumericVector zp = qnorm(p/2);
  NumericVector zq = -1.0*qnorm(q/2);

  if (any(!is_finite(zp+zq))) {
    stop("P-values p,q must be in [1e-300,1]");
  }

  // maximum limits for integration
  double mx = max(NumericVector::create(10.0, max(abs(-qnorm(p/2)))));
  double my = max(NumericVector::create(10.0, max(abs(-qnorm(q/2)))));

  //gx=0 //1/1000  // use for 'tiebreaks'- if a point is on a curve with nonzero area, force the L-curve through that point
  
  if (indices==R_NilValue) {
    if (at==R_NilValue) {
      stop("One of the parameters 'indices', 'at' must be set");
    }
    NumericVector ccut=at;
    int mode=0;
  } else {
    // Note that the vector constructor is mimicking rep(0, length(indices))
    NumericVector ccut(indices.size());
  }

  
  // yval2=seq(from=0,to=my,length.out=nv+1)[1:nv];
  NumericVector yval2(nv+1);

  for(int i = 0; i < (nv+1); i++) {
    yval2[i] = my * (i/nv); 
  }
  
  yval2 = yval2[Range(0,nv-1)];

  /* 
  xval2=outer(rep(1,length(ccut)),yval2);
  pval2=2*pnorm(-yval2);
  xtest=seq(0,mx,length.out=nt);
  ptest=2*pnorm(-xtest);
  
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
  List out(2);

  return out;
}
