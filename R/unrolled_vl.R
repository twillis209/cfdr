                                        #vl=function(p,q,adj=TRUE,indices=NULL,at=NULL,mode=0,fold=NULL,nt=5000, nv=1000, p_threshold=0, scale=c("p","z"), closed=TRUE,verbose=FALSE,gx=10^-5) {

set.seed(1)

n=10000; n1p=100; n1pq=100; n1q=100

zp=c(rnorm(n1p,sd=3), rnorm(n1q,sd=1),
     rnorm(n1pq,sd=3), rnorm(n-n1p-n1q-n1pq,sd=1))

zq=c(rnorm(n1p,sd=1), rnorm(n1q,sd=3), 
     rnorm(n1pq,sd=3), rnorm(n-n1p-n1q-n1pq,sd=1))

p=2*pnorm(-abs(zp))
q=2*pnorm(-abs(zq))

fold_id=(1:n) %% 3

candidate_indices=which(p<0.01 | q< 0.001)

fold1=which(fold_id==0) 
fold2=which(fold_id==1) 
fold3=which(fold_id==2)

ind1=intersect(candidate_indices,fold1) 
ind2=intersect(candidate_indices,fold2) 
ind3=intersect(candidate_indices,fold3)

adj=TRUE
# default is indices=NULL
indices=ind1
at=NULL
# default is mode=0
mode=2
# default is fold=NULL
fold=fold1
nt=5000
nv=1000
p_threshold=0
scale=c("p", "z")
closed=TRUE
verbose=FALSE
gx=10^-5

# internal interpolation function
ff=function(xtest,xx,cxi) {
  if (xx[1] <= cxi) return(xtest[1]); if (xx[length(xx)] >= cxi ) return(xtest[length(xx)])
  w1=max(which(xx> cxi)); w2=1+w1 #min(which(xx< cxi))
  xtest[w1] + 
    (xtest[w2] - xtest[w1])*
    (cxi-xx[w1])/(xx[w2]-xx[w1])
}
  
  zp=-qnorm(p/2); zq=-qnorm(q/2)
  if (any(!is.finite(zp+zq))) stop("P-values p,q must be in [1e-300,1]")
  mx=max(c(10,abs(-qnorm(p/2)))); my=max(c(10,abs(-qnorm(q/2)))); # maximum limits for integration
  
  #gx=0 #1/1000  # use for 'tiebreaks'- if a point is on a curve with nonzero area, force the L-curve through that point
  
  if (is.null(indices)) {
    if (is.null(at)) stop("One of the parameters 'indices', 'at' must be set")
    ccut=at
    mode=0
  } else ccut=rep(0,length(indices))
  
  yval2=seq(0,my,length.out=nv+1)[1:nv]; xval2=outer(rep(1,length(ccut)),yval2); pval2=2*pnorm(-yval2)
  xtest=seq(0,mx,length.out=nt); ptest=2*pnorm(-xtest)
  
  
  if (!is.null(indices)) { # set ccut. NOT equal to cfdr at the points; needs to be adjusted since an additional test point is used in determination of L
    if (mode==1) {
      ccut=rep(0,length(indices))
      for (i in 1:length(indices)) {
        w=which(zq[-indices[i]] >= zq[indices[i]])
        if (length(w) >= 1) {
          cfsub= (1+(1/length(w)))*ptest/(1+(1/length(w))-ecdf(zp[-indices[i]][w])(xtest)); cfsub=cummin(cfsub)
          ccut[i]=approx(xtest,cfsub-gx*xtest + gx*mx,zp[indices[i]],rule=2)$y
        } else ccut[i]=p[indices[i]]*1/(1+ length(which(zp[-indices[i]][w]>= zp[indices[i]])))
      }
    } 
    if (mode==2) {
      ccut=rep(0,length(indices))
      for (i in 1:length(indices)) {
        w=which(zq[-fold] >= zq[indices[i]])
        if (length(w) >= 1) {
          cfsub= (1+(1/length(w)))*ptest/(1+(1/length(w))-ecdf(zp[-fold][w])(xtest)); cfsub=cummin(cfsub)
          ccut[i]=approx(xtest,cfsub-gx*xtest + gx*mx,zp[indices[i]],rule=2)$y
        } else ccut[i]=p[indices[i]]
      }
      #    ccut= p[indices]*sapply(indices,function(x) 
      #      (1+length(which(q[-fold] <= q[x])))/(1+length(which(p[-fold] <= p[x] & q[-fold] <= q[x]))))# set ccut
    }
    if (mode==0) {
      #  yval2=seq(0,my,length.out=nv+1)[1:nv]; xval2=outer(rep(1,length(indices)),yval2); pval2=2*pnorm(-yval2)
      #  xtest=seq(0,my,length.out=nt); ptest=2*pnorm(-xtest)
      
      ccut=rep(0,length(indices))
      for (i in 1:length(indices)) {
        w=which(zq >= zq[indices[i]])
        if (length(w) >= 2) {
          cfsub= (1+(1/length(w)))*ptest/(1+(1/length(w))-ecdf(zp[w])(xtest)); cfsub=cummin(cfsub)
          ccut[i]=approx(xtest,cfsub-gx*xtest + gx*mx,zp[indices[i]],rule=2)$y
        } else ccut[i]=p[indices[i]]
      }
      #    ccut=p[indices]*sapply(indices,function(x) 
      #      (1+length(which(q <= q[x])))/(1+length(which(p <= p[x] & q <= q[x]))))# set ccut
    }
  }  
  
  ccut=ccut*(1+ 1e-6) # ccut + 1e-8 # prevent floating-point comparison errors
  
  if (verbose & mode==1) print(paste0(length(ccut)," regions to calculate"))
  
  out=rep(0,length(ccut))
  
  
  
  if (mode==0) {
    
    if (adj) {
      correct=cummin((1+ecdf(q[which(p>0.5)])(pval2)*length(p))/
          (1+ ecdf(q)(pval2)*length(p)))
      if (!is.null(indices)) correct_ccut=approx(pval2,correct,q[indices],rule=2)$y #cummin((1+ecdf(pc[which(p>0.5)])(pc)*length(p))/(1+ rank(pc)))
    } else {
      correct=rep(1,length(pval2)) # adjustment factor for pc[i]
      correct_ccut=rep(1,length(ccut))
    }
    if (!is.null(indices)) ccut=ccut*correct_ccut
    
    zp_ind=ceiling(zp*nt/mx); zp_ind=pmax(1,pmin(zp_ind,nt))
    zq_ind=ceiling(zq*nv/my); zq_ind=pmax(1,pmin(zq_ind,nv))
    zq[which(zq>my)]=my; zp[which(zp>mx)]=mx; p=2*pnorm(-abs(zp)); q=2*pnorm(-abs(zq))
    
    for (i in 1:length(yval2)) {
      w=which(zq > yval2[i])
      if (length(w)>= 1) {
        cfsub= ptest/(1+(1/length(w))-ecdf(zp[w])(xtest)); cfsub=cummin(cfsub)
        xval2[,i]=approx(cfsub*correct[i]-gx*xtest + gx*mx,xtest,ccut,rule=2,method="const",f=1)$y 
      }  else xval2[,i]=-qnorm(ccut/2)
    }
    
  } 
  if (mode==1) {
    zp_ind=ceiling(zp*nt/mx); zp_ind=pmax(1,pmin(zp_ind,nt))
    zq_ind=ceiling(zq*nv/my); zq_ind=pmax(1,pmin(zq_ind,nv))
    # populate zp_ind[i],zq_ind[i] with the index of xval2,yval2 closest (least above?) zp[i],zq[i]. Only need to do at 'sub'
    
    xmat=matrix(0,nt,nv) # base
    
    pqmat=xmat
    qvec=rep(0,nv)
    
    for (i in 1:length(p)) {
      pqmat[1:zp_ind[i],1:zq_ind[i]]=1+pqmat[1:zp_ind[i],1:zq_ind[i]]
      qvec[1:zq_ind[i]]=1+qvec[1:zq_ind[i]]
    }
    # matrix [i,j] is 1+ number of SNPs with zp>xval2[i]; zq>yval2[i], only works for i in sub
    pqmat=1+pqmat
    qvec=1+qvec
    
    #qvec=1+  length(zq_ind)*(1-ecdf(zq_ind)(1:nv))
    # qvec[i] is 1+ number of SNPs with zq>yval[i], only for i in sub
    
    cf_mat=outer(ptest,qvec)/pqmat #pmat*qmat/pqmat
    # [i,j] has cFDR value at xtest[i], ytest[i]; only valid for i in sub
    cf_mat=apply(cf_mat,2,cummin) # Smooth shape of L, see paragraph in methods
    
    l_new=rep(0,length(p)) # set to L in new method
    
    
    for (i in 1:length(indices)) {
      
      if (adj) correctx=cummin((1+ecdf(q[-indices[i]][which(p[-indices[i]]>0.5)])(pval2)*length(p[-indices[i]]))/(1+ ecdf(q[-indices[i]])(pval2)*length(p[-indices[i]]))) else correctx=rep(1,length(pval2)) 
      
      ccut[i]=ccut[i]*approx(yval2,correctx,zq[indices[i]],rule=2)$y #ff(rev(correctx),rev(yval2),zq[indices[i]]) #correctx[zq_ind[indices[i]]]
      
      pqnew=pqmat
      pqnew[1:zp_ind[indices[i]],1:zq_ind[indices[i]]]=-1+pqnew[1:zp_ind[indices[i]],1:zq_ind[indices[i]]] #remove contribution of current SNP
      qnew=qvec
      qnew[1:zq_ind[indices[i]]]=-1+qnew[1:zq_ind[indices[i]]]
      cfx=apply(outer(ptest,qnew)/pqnew,2,cummin)
      cfx=t(t(cfx- gx*xtest + gx*mx)*correctx)
      cxi=ccut[i]
      xv=suppressWarnings(apply(cfx,2,function(x) ff(xtest,x,cxi))) #approx(x,xtest,cxi,rule=2)$y))
      xv[which(xv<0)]=0; 
      xv[which(!is.finite(xv))]=0
      xval2[i,]=xv;
      if (verbose) print(i)
    }
  }
  if (mode==2) {
    
    if (adj) {
      correct=cummin((1+ecdf(q[-fold][which(p[-fold]>0.5)])(pval2)*length(p[-fold]))/
          (1+ ecdf(q[-fold])(pval2)*length(p[-fold])))
      if (!is.null(indices)) correct_ccut=approx(pval2,correct,q[indices],rule=2)$y #cummin((1+ecdf(pc[which(p>0.5)])(pc)*length(p))/(1+ rank(pc)))
    } else {
      correct=rep(1,length(pval2)) # adjustment factor for pc[i]
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
