library(testthat)
library(Rcpp)
library(microbenchmark)

sourceCpp('vl.cpp')

test_that("Test ecdf_cpp on some trivial cases", {
  testSample<-seq(0.1, 1.0, length.out=10)
  shufSample<-sample(testSample)
  expect_equal(ecdf_cpp(testSample, testSample), testSample)
  expect_equal(ecdf_cpp(shufSample, shufSample), shufSample)
  expect_equal(ecdf_cpp(testSample, c(0.1)), c(0.1))
})

test_that('Test ecdf_cpp with duplicates in the sample', {
  testSample<-c(0.1,0.2,0.2,0.3,0.4)
  expect_equal(ecdf_cpp(testSample, testSample), c(0.2,0.6,0.6,0.8,1.0))
})

test_that('Test ecdf_cpp with standard normal samples', {
  ref<-rnorm(1e3)
  sam<-rnorm(1e2)
  expect_equal(ecdf_cpp(ref, sam), ecdf(ref)(sam))
})

set.seed(42)

ref<-rnorm(1e3)
sam<-rnorm(1e2)

mbm<-microbenchmark(
  "stats::ecdf"={
    ecdf(ref)(sam)
  },
  "ecdf_cpp"={
    ecdf_cpp(ref, sam)
  },
  times=20
  )
