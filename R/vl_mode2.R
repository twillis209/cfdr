# An R equivalent to my Rcpp vl_mode2 to guide its development

Rvl_mode2=function(p,q,indices,fold,adj=TRUE,at=NULL,nt=5000, nv=1000, p_threshold=0, scale=c("p","z"), closed=TRUE,verbose=FALSE,gx=10^-5) {
  
  zp=-qnorm(p/2)
  zq=-qnorm(q/2)

  if (any(!is.finite(zp+zq))) stop("P-values p,q must be in [1e-300,1]")

  mx=max(c(10,abs(-qnorm(p/2))))
  my=max(c(10,abs(-qnorm(q/2)))) # maximum limits for integration

  #print(sprintf("mx: %.5f my: %.5f", mx, my))
  
  ccut=rep(0,length(indices))
  
  yval2=seq(0,my,length.out=nv+1)[1:nv]

  xval2=outer(rep(1,length(ccut)),yval2)

  pval2=2*pnorm(-yval2)

  xtest=seq(0,mx,length.out=nt)

  ptest=2*pnorm(-xtest)
  
  # set ccut. NOT equal to cfdr at the points; needs to be adjusted since an additional test point is used in determination of L
  for (i in 1:length(indices)) {




      w=which(zq[-fold] >= zq[indices[i]])

      # TODO remove me
      print(sprintf("length of w: %d", length(w)))
      #print(sprintf("length(zp_negFold): %d length(zq_negFold): %d", length(zp[-fold]), length(zq[-fold])))

      if (length(w) >= 1) {
        cfsub= (1+(1/length(w)))*ptest/(1+(1/length(w))-ecdf(zp[-fold][w])(xtest)); cfsub=cummin(cfsub)
        ccut[i]=approx(xtest,cfsub-gx*xtest + gx*mx,zp[indices[i]],rule=2)$y
      } else ccut[i]=p[indices[i]]
    }
    
  
  ccut=ccut*(1+ 1e-6) # ccut + 1e-8 # prevent floating-point comparison errors
  
    if (adj) {
      correct=cummin((1+ecdf(q[-fold][which(p[-fold]>0.5)])(pval2)*length(p[-fold]))/
          (1+ ecdf(q[-fold])(pval2)*length(p[-fold])))
      if (!is.null(indices)) correct_ccut=approx(pval2,correct,q[indices],rule=2)$y
    } else {
      correct=rep(1,length(pval2)) # adjustment factor for pc[i]
      correct_ccut=rep(1,length(ccut))
    }
    if (!is.null(indices)) ccut=ccut*correct_ccut
    
  zp_ind=ceiling(zp[indices]*nt/mx)
  zp_ind=pmax(1,pmin(zp_ind,nt))
  zq_ind=ceiling(zq[indices]*nv/my)
  zq_ind=pmax(1,pmin(zq_ind,nv))
    
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
  
}
