source('functions.R')
library(data.table)
library(parallel)
library(optimParallel)
library(tictoc)

setDTthreads(threads=4)

dat<-fread('../../pidProject/data/joins/joinedData/pid/aster_pid_flipped.tsv.gz', sep='\t', header=T, select=c('P.2', 'P.1', 'CHR38.1', 'BP38.1', 'SNPID.1'))

set.seed(27)

dat<-na.omit(dat, cols=c('P.2', 'P.1', 'CHR38.1'))

folds<-mclapply(unique(dat[['CHR38.1']]), function(x) which(dat[['CHR38.1']]==x), mc.cores=4)

p<-dat[['P.2']]
names(p)<-1:length(p)
q<-dat[['P.1']]

# The vl function requires p-values in the interval [1e-300, 1]
p[which(p<1e-300)]<-1e-300
q[which(q<1e-300)]<-1e-300

# 109 p-values below 1e-7
candidate_indices<-which(p<1e-7)

# Taken mostly from the tictoc manual with some sugar added 
timeMe<-function(noOfReplicates, fn, ...) {
tic.clearlog()
for(i in 1:noOfReplicates) {
  tic(i)
  fn(...)
  toc(log=TRUE, quiet=TRUE)
}
log.lst<-tic.log(format=F)
tic.clearlog()

return(unlist(lapply(log.lst, function(x) x$toc-x$tic)))
}

# Correctness of parallel implementation
fit.2g.serial<-fit.2g(q[which(p>0.5)])
fit.2g.parallel_1<-fit.2g.parallel(q[which(p>0.5)], ncores=1)
fit.2g.parallel_2<-fit.2g.parallel(q[which(p>0.5)], ncores=2)
fit.2g.parallel_4<-fit.2g.parallel(q[which(p>0.5)], ncores=4)
fit.2g.parallel_8<-fit.2g.parallel(q[which(p>0.5)], ncores=8)

# Timing
serialTimes<-timeMe(5, fit.2g, q[which(p>0.5)])
parallelOneCoresTimes<-timeMe(5, fit.2g.parallel, q[which(p>0.5)], ncores=1)
parallelTwoCoresTimes<-timeMe(5, fit.2g.parallel, q[which(p>0.5)], ncores=2)
parallelFourCoresTimes<-timeMe(5, fit.2g.parallel, q[which(p>0.5)], ncores=4)
parallelEightCoresTimes<-timeMe(5, fit.2g.parallel, q[which(p>0.5)], ncores=8)
