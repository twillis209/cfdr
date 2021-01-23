set.seed(1)

n=10000; n1p=100; n1pq=100; n1q=100

zp=c(rnorm(n1p,sd=3), rnorm(n1q,sd=1),
     rnorm(n1pq,sd=3), rnorm(n-n1p-n1q-n1pq,sd=1))

zq=c(rnorm(n1p,sd=1), rnorm(n1q,sd=3), 
     rnorm(n1pq,sd=3), rnorm(n-n1p-n1q-n1pq,sd=1))

p=2*pnorm(-abs(zp))
q=2*pnorm(-abs(zq))

true_q0_pars=c(1- n1q/(n-n1p-n1pq),2)
est_q0_pars=fit.2g(q[which(p>0.5)])$pars

fold_id=(1:n) %% 3

candidate_indices=which(p<0.01 | q< 0.001)

fold1=which(fold_id==0) 
fold2=which(fold_id==1) 
fold3=which(fold_id==2)

ind1=intersect(candidate_indices,fold1) 
ind2=intersect(candidate_indices,fold2) 
ind3=intersect(candidate_indices,fold3)

test_that("Test full Rcpp implementation", {
  Rcpp_vl<-vl_mode2(p=p,q=q,indices=ind1,fold=fold1)
  v1<-vl(p,q,indices=ind1,mode=2,fold=fold1)

  X<-Rcpp_vl$x
  Y<-Rcpp_vl$y
  r_X<-v1$x
  r_Y<-v1$y

  expect_equal(r_X,X)
  expect_equal(r_Y,Y)
})

test_that("Rcpp implementation of vl reproduces behaviour of R implementation of vl w.r.t. v-value output", {
  v1=vl(p,q,indices=ind1,mode=2,fold=fold1)  
  v2=vl(p,q,indices=ind2,mode=2,fold=fold2) 
  v3=vl(p,q,indices=ind3,mode=2,fold=fold3)

  v_fold=rep(1,n)

  v_fold[ind1]=il(v1,pi0_null=est_q0_pars[1],sigma_null=est_q0_pars[2],distx="norm")
  v_fold[ind2]=il(v2,pi0_null=est_q0_pars[1],sigma_null=est_q0_pars[2],distx="norm")
  v_fold[ind3]=il(v3,pi0_null=est_q0_pars[1],sigma_null=est_q0_pars[2],distx="norm")

  v1_Rcpp=vl_mode2(p,q,indices=ind1,fold=fold1)
  v2_Rcpp=vl_mode2(p,q,indices=ind2,fold=fold2)
  v3_Rcpp=vl_mode2(p,q,indices=ind3,fold=fold3)

  v_fold_rcpp = rep(1,n)
  v_fold_rcpp[ind1]=il(v1_Rcpp,pi0_null=est_q0_pars[1],sigma_null=est_q0_pars[2],distx="norm")
  v_fold_rcpp[ind2]=il(v2_Rcpp,pi0_null=est_q0_pars[1],sigma_null=est_q0_pars[2],distx="norm")
  v_fold_rcpp[ind3]=il(v3_Rcpp,pi0_null=est_q0_pars[1],sigma_null=est_q0_pars[2],distx="norm")

  expect_equal(v_fold, v_fold_rcpp)
})

test_that("Rcpp implementation of vl reproduces behaviour of R implementation of vl w.r.t. BH procedure output", {
  v1=vl(p,q,indices=ind1,mode=2,fold=fold1)  
  v2=vl(p,q,indices=ind2,mode=2,fold=fold2) 
  v3=vl(p,q,indices=ind3,mode=2,fold=fold3)

  v_fold=rep(1,n)

  v_fold[ind1]=il(v1,pi0_null=est_q0_pars[1],sigma_null=est_q0_pars[2],distx="norm")
  v_fold[ind2]=il(v2,pi0_null=est_q0_pars[1],sigma_null=est_q0_pars[2],distx="norm")
  v_fold[ind3]=il(v3,pi0_null=est_q0_pars[1],sigma_null=est_q0_pars[2],distx="norm")

  hit_fold=bh(v_fold,0.1)

  v1_Rcpp=vl_mode2(p,q,indices=ind1,fold=fold1)
  v2_Rcpp=vl_mode2(p,q,indices=ind2,fold=fold2)
  v3_Rcpp=vl_mode2(p,q,indices=ind3,fold=fold3)

  v_fold_rcpp = rep(1,n)
  v_fold_rcpp[ind1]=il(v1_Rcpp,pi0_null=est_q0_pars[1],sigma_null=est_q0_pars[2],distx="norm")
  v_fold_rcpp[ind2]=il(v2_Rcpp,pi0_null=est_q0_pars[1],sigma_null=est_q0_pars[2],distx="norm")
  v_fold_rcpp[ind3]=il(v3_Rcpp,pi0_null=est_q0_pars[1],sigma_null=est_q0_pars[2],distx="norm")

  hit_fold_rcpp=bh(v_fold_rcpp,0.1)

  expect_equal(hit_fold, hit_fold_rcpp)
})
