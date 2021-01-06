library(cfdr)

set.seed(1)                                           # ensure reproducibility
n=10000; n1p=100; n1pq=100; n1q=100                   # parameters
zp=c(rnorm(n1p,sd=3), rnorm(n1q,sd=1),
     rnorm(n1pq,sd=3), rnorm(n-n1p-n1q-n1pq,sd=1))    # simulate z-scores corresponding to p
zq=c(rnorm(n1p,sd=1), rnorm(n1q,sd=3), 
     rnorm(n1pq,sd=3), rnorm(n-n1p-n1q-n1pq,sd=1))    # simulate z-scores corresponding to q
p=2*pnorm(-abs(zp))
q=2*pnorm(-abs(zq))              # convert to p-values

fold_id=(1:n) %% 3

candidate_indices=which(p<0.01 | q< 0.001)

fold1=which(fold_id==0) 
fold2=which(fold_id==1) 
fold3=which(fold_id==2)

ind1=intersect(candidate_indices,fold1) 
ind2=intersect(candidate_indices,fold2) 
ind3=intersect(candidate_indices,fold3) 

v1=vl(p,q,indices=ind1,mode=2,fold=fold1);  # parameter fold specifies indices to leave out 
v2=vl(p,q,indices=ind2,mode=2,fold=fold2); 
v3=vl(p,q,indices=ind3,mode=2,fold=fold3);

save(v1, v2, v3, file='vl_testData.RData')
