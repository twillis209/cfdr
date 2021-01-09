library(cfdr)
library(Rcpp)
library(testthat)
library(microbenchmark)

sourceCpp('vl.cpp')
source('vl_mode2.R')
                                        # The idea here is to repeat the analysis from the vignette but cache the result of the call to vl so that we have something against which to validate the behaviour of the Rcpp implementation of this function

# Contains objects v1, v2, and v3
#load('vl_testData.RData')

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


R_v1<-Rvl_mode2(p,q,indices=ind1,fold=fold1)

Rcpp_v1<-vl_mode2(p=p,q=q,indices=ind1,fold=fold1)

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

# pval2 max difference is 1.11e-16
                                        # q[indices(-1)] is equal, too
# pval2 is sorted from largest to smallest, but the values are extremely small
test_that("Test part 3 of Rcpp-implementation", {
  R_v1_part3<-Rvl_mode2_part3(p,q,indices=ind1,fold=fold1)
  Rcpp_v1_part3<-vl_mode2_part3(p,q,ind1,fold1)
  expect_equal(R_v1_part3[["ccut"]], Rcpp_v1_part3[["ccut"]])
  expect_equal(R_v1_part3[["correct_ccut"]], Rcpp_v1_part3[["correct_ccut"]])
  expect_equal(R_v1_part3[["correct"]], Rcpp_v1_part3[["correct"]])
})


#test_that("Test R-based implementation of vl", {
#  test_v1=vl(p,q,indices=ind1,mode=2,fold=fold1)   
#  test_v2=vl(p,q,indices=ind2,mode=2,fold=fold2) 
#  test_v3=vl(p,q,indices=ind3,mode=2,fold=fold3)
#  expect_equal(v1, test_v1)
#  expect_equal(v2, test_v2)
#  expect_equal(v3, test_v3)
#})

