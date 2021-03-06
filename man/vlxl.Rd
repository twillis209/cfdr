% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/functions.R
\name{vlxl}
\alias{vlxl}
\title{vlxl}
\usage{
vlxl(p, q, pars, indices = NULL, at = NULL, p_threshold = 0,
  nt = 5000, nv = 1000, scale = c("p", "z"), closed = T, up = T)
}
\arguments{
\item{p}{principal p-values}

\item{q}{conditional p-values}

\item{pars}{fitted parameters governing distribution of (P,Q). A seven element-vector: (pi0,pi1,pi2,s1,s2,t1,t2) respectively. See function description.}

\item{indices}{compute v(L) at indices of points. Overrides parameter at if set.}

\item{at}{cfdr cutoff/cutoffs. Defaults to null}

\item{p_threshold}{if H0 is to be rejected automatically whenever p<p_threshold, include this in all regions L}

\item{nt}{number of test points in x-direction, default 5000}

\item{nv}{resolution for constructing L-curves, default 1000}

\item{scale}{return curves on the p- or z- plane. Y values are equally spaced on the z-plane.}

\item{closed}{determines whether curves are closed polygons encircling regions L (closed=T), or lines indicating the rightmost border of regions L}
}
\value{
list containing elements x, y. Assuming n curves are calculated and closed=T (where n=length(indices) or length(at)), x is a matrix of dimension n x (4+nv), y ix a vector of length (4+nv).
}
\description{
Return co-ordinates of L-regions for cFDR method using local four-groups method. Automatically includes an estimate of Pr(H0|Q<q).
}
\details{
Assumes (P,Q) follows a bivariate mixture-Gaussian distribution with four components, each centred at
the origin and with covariance matrices I2, (s1^2,0; 0,1), (1,0; 0,t1^2), (s2^2,0; 0,t2^2). Components
have weights pi0,pi1,pi2 and (1-pi0-pi1-pi2) respectively. Function fit.4g fits parameters to data.

This function does not have a 'mode' option (as for function vl) since the mapping from
[0,1]^2 -> [0,1] defined by L-curves depends only on 'pars' in this case and not on observed p,q.
}
\examples{
# Generate standardised simulated dataset
set.seed(1); n=10000; n1p=100; n1pq=100; n1q=100
zp=c(rnorm(n1p,sd=3), rnorm(n1q,sd=1),rnorm(n1pq,sd=3), rnorm(n-n1p-n1q-n1pq,sd=1))
zq=c(rnorm(n1p,sd=1), rnorm(n1q,sd=3),rnorm(n1pq,sd=3), rnorm(n-n1p-n1q-n1pq,sd=1))
p=2*pnorm(-abs(zp)); q=2*pnorm(-abs(zq))
fold_id=(1:n) \%\% 3

# estimate parameters of underying dataset
fit_pars=fit.4g(cbind(zp,zq))$pars

# estimate parameters of underying dataset, removing fold 1
fit_pars_fold23=fit.4g(cbind(zp[which(fold_id!=1)],zq[which(fold_id!=1)]))$pars


# points to generate L-regions for
example_indices=c(4262, 268,83,8203)

par(mfrow=c(1,2))

v=vlx(p,q,pars=fit_pars,indices=example_indices); 
plot(p,q,cex=0.5,col="gray",xlim=c(0,0.001),ylim=c(0,1), main="All points in"); 
for (i in 1:length(example_indices)) lines(v$x[i,],v$y); 
for (i in 1:length(example_indices)) points(p[example_indices[i]],q[example_indices[i]],pch=16,col="blue")

v2=vlx(p,q,pars=fit_pars_fold23,indices=example_indices); 
plot(p,q,cex=0.5,col="gray",xlim=c(0,0.001),ylim=c(0,1), main="Fold removed"); 
for (i in 1:length(example_indices)) lines(v2$x[i,],v2$y); 
for (i in 1:length(example_indices)) points(p[example_indices[i]],q[example_indices[i]],pch=16,col="blue")


}
\author{
James Liley
}
