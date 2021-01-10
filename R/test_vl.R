library(cfdr)
library(Rcpp)
library(testthat)
library(microbenchmark)

sourceCpp('vl.cpp')
source('vl_mode2.R')

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

# 'U' for ur-implementation
U_v1<-vl(p=p,q=q,indices=ind1,fold=fold1,mode=2)

test_that("Test part 1 of Rcpp-implementation", {
  R_v1_part1<-Rvl_mode2_part1(p,q,indices=ind1,fold=fold1)
  Rcpp_v1_part1<-vl_mode2_part1(p,q,ind1,fold1)
  expect_equal(R_v1_part1["zp"], Rcpp_v1_part1["zp"])
  expect_equal(R_v1_part1["zq"], Rcpp_v1_part1["zq"])
  expect_equal(R_v1_part1["my"], Rcpp_v1_part1["my"])
  expect_equal(R_v1_part1["mx"], Rcpp_v1_part1["mx"])
  expect_equal(R_v1_part1["xval2"], Rcpp_v1_part1["xval2"])
  expect_equal(R_v1_part1["pval2"], Rcpp_v1_part1["pval2"])
  expect_equal(R_v1_part1["ptest"], Rcpp_v1_part1["ptest"])
})

test_that("Test part 2 of Rcpp-implementation", {
  R_v1_part2<-Rvl_mode2_part2(p,q,indices=ind1,fold=fold1)
  Rcpp_v1_part2<-vl_mode2_part2(p,q,ind1,fold1)
  expect_equal(R_v1_part2["ccut"], Rcpp_v1_part2["ccut"])
})

test_that("Test part 3 of Rcpp-implementation", {
  R_v1_part3<-Rvl_mode2_part3(p,q,indices=ind1,fold=fold1)
  Rcpp_v1_part3<-vl_mode2_part3(p,q,ind1,fold1)
  expect_equal(R_v1_part3[["ccut"]], Rcpp_v1_part3[["ccut"]])
  expect_equal(R_v1_part3[["correct_ccut"]], Rcpp_v1_part3[["correct_ccut"]])
  expect_equal(R_v1_part3[["correct"]], Rcpp_v1_part3[["correct"]])
})

test_that("Test part 4 of Rcpp-implementation", {
  R_v1_part4<-Rvl_mode2_part4(p,q,ind1,fold1)
  Rcpp_v1_part4<-vl_mode2_part4(p,q,ind1,fold1)
  expect_equal(R_v1_part4[["xval2"]], Rcpp_v1_part4[["xval2"]])
  expect_equal(R_v1_part4[["zp_ind"]], Rcpp_v1_part4[["zp_ind"]])
  expect_equal(R_v1_part4[["zq_ind"]], Rcpp_v1_part4[["zq_ind"]])
  expect_equal(R_v1_part4[["cfsub_mat"]], Rcpp_v1_part4[["cfsub_mat"]])
  expect_equal(R_v1_part4[["ecdf_mat"]], Rcpp_v1_part4[["ecdf_mat"]])
})

test_that("Test part 5 of Rcpp-implementation", {
  R_v1_part5<-Rvl_mode2_part5(p,q,ind1,fold1)
  Rcpp_v1_part5<-vl_mode2_part5(p,q,ind1,fold1)
  expect_equal(R_v1_part5[["xval2"]], Rcpp_v1_part5[["xval2"]])
  expect_equal(R_v1_part5[["yval2"]], Rcpp_v1_part5[["yval2"]])
})

# It seems that these actually match very closely
test_that("Test part 5.5 of Rcpp implementation", {
  Rcpp_v1_part5<-vl_mode2_part5(p,q,ind1,fold1)
  R_v1_part6<-Rvl_mode2_part6(p,q,ind1,fold1)
  X=2*pnorm(-abs(Rcpp_v1_part5$xval2))
  Y=2*pnorm(-abs(Rcpp_v1_part5$yval2))
  r_X<-R_v1_part6$x
  r_Y<-R_v1_part6$y
})

# TODO the test passes but the values of X and Y differ, look at X-r_X
# It seems that the 
test_that("Test part 6 of Rcpp-implementation", {
  R_v1_part6<-Rvl_mode2_part6(p,q,ind1,fold1)
  Rcpp_v1_part6<-vl_mode2_part6(p,q,ind1,fold1)
  X<-Rcpp_v1_part6$x
  Y<-Rcpp_v1_part6$y
  r_X<-R_v1_part6$x
  r_Y<-R_v1_part6$y
  expect_equal(r_X,X)
  expect_equal(r_Y,Y)
})

test_that("Rcpp implementation reproduces behaviour of R implementation w.r.t. BH procedure output", {
  v1=vl(p,q,indices=ind1,mode=2,fold=fold1);  
  v2=vl(p,q,indices=ind2,mode=2,fold=fold2); 
  v3=vl(p,q,indices=ind3,mode=2,fold=fold3); 
  v_fold=rep(1,n)
  v_fold[ind1]=il(v1,pi0_null=est_q0_pars[1],sigma_null=est_q0_pars[2],distx="norm")
  v_fold[ind2]=il(v2,pi0_null=est_q0_pars[1],sigma_null=est_q0_pars[2],distx="norm")
  v_fold[ind3]=il(v3,pi0_null=est_q0_pars[1],sigma_null=est_q0_pars[2],distx="norm")
  hit_fold=bh(v_fold,0.1)

  v1_Rcpp=vl_mode2_part6(p,q,indices=ind1,fold=fold1)
  v2_Rcpp=vl_mode2_part6(p,q,indices=ind2,fold=fold2)
  v3_Rcpp=vl_mode2_part6(p,q,indices=ind3,fold=fold3)

  v_fold[ind1]=il(v1_Rcpp,pi0_null=est_q0_pars[1],sigma_null=est_q0_pars[2],distx="norm")
  v_fold[ind2]=il(v2_Rcpp,pi0_null=est_q0_pars[1],sigma_null=est_q0_pars[2],distx="norm")
  v_fold[ind3]=il(v3_Rcpp,pi0_null=est_q0_pars[1],sigma_null=est_q0_pars[2],distx="norm")
  hit_fold_rcpp=bh(v_fold,0.1)
  expect_equal(hit_fold, hit_fold_rcpp)
})

mbm<-microbenchmark(
  "R"={
    v1=vl(p,q,indices=ind1,mode=2,fold=fold1);  
    v2=vl(p,q,indices=ind2,mode=2,fold=fold2); 
    v3=vl(p,q,indices=ind3,mode=2,fold=fold3); 
  },
  "Rcpp"={
    v1_Rcpp=vl_mode2_part6(p,q,indices=ind1,fold=fold1)
    v2_Rcpp=vl_mode2_part6(p,q,indices=ind2,fold=fold2)
    v3_Rcpp=vl_mode2_part6(p,q,indices=ind3,fold=fold3)
  },
  times=100
)
