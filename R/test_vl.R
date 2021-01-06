library(testthat)
library(cfdr)
library(Rcpp)

sourceCpp('vl.cpp')
                                        # The idea here is to repeat the analysis from the vignette but cache the result of the call to vl so that we have something against which to validate the behaviour of the Rcpp implementation of this function

# Contains objects v1, v2, and v3
load('vl_testData.RData')

                                        # Not sure how we should load this code
test_that("Test ecdf_cpp on some trivial cases", {
  testSample<-seq(0.1, 1.0, length.out=10)
  expect_equal(ecdf_cpp(testSample, testSample), testSample)
})


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

test_that("Test R-based implementation of vl", {
  test_v1=vl(p,q,indices=ind1,mode=2,fold=fold1)   
  test_v2=vl(p,q,indices=ind2,mode=2,fold=fold2) 
  test_v3=vl(p,q,indices=ind3,mode=2,fold=fold3)
  expect_equal(v1, test_v1)
  expect_equal(v2, test_v2)
  expect_equal(v3, test_v3)
})
