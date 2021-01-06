library(testthat)
library(Rcpp)
library(microbenchmark)

sourceCpp('vl.cpp')

test_that("Test approx_cpp on some trivial cases", {
  x<-1:10
  y<-rnorm(10)
  expect_equal(approx_cpp(x,y,x), approx(x,y,x,rule=2)$y)
})

test_that("Handles xout value smaller than min x in accordance with rule 2", {
  x<-1:10
  expect_equal(approx_cpp(x,x,c(0)), c(1));
})

test_that("Handles xout value larger than max x in accordance with rule 2", {
  x<-1:10
  expect_equal(approx_cpp(x,x,c(11)), c(10));
})

test_that("Handles typical use case", {
  x<-1:10
  y<-rnorm(10)
  xout<-c(0,seq(1,10,length.out=51),11)
  expect_equal(approx_cpp(x,y,xout), approx(x,y,xout,rule=2)$y)
})

