set.seed(42)

n=1e3; n1p=1e2; n1pq=1e2; n1q=1e2

zp=c(rnorm(n1p,sd=3), rnorm(n1q,sd=1), rnorm(n1pq,sd=3), rnorm(n-n1p-n1q-n1pq,sd=1))

zq=c(rnorm(n1p,sd=1), rnorm(n1q,sd=3), rnorm(n1pq,sd=3), rnorm(n-n1p-n1q-n1pq,sd=1))

p=2*pnorm(-abs(zp)); q=2*pnorm(-abs(zq))

est_q0_pars=fit.2g(q[which(p>0.5)])$pars

test_that("fit.2g.parallel estimates mixture parameters close to those estimated by the serial implementation when run with one core", {
  expect_equal(est_q0_pars, fit.2g.parallel(q[p>0.5], ncores=1)$pars)
})

test_that("fit.2g.parallel estimates mixture parameters close to those estimated by the serial implementation when run with two cores", {
  skip_if(parallel::detectCores()==1, "Need two cores for test")
  expect_equal(est_q0_pars, fit.2g.parallel(P=q[p>0.5], ncores=2)$pars)
})
