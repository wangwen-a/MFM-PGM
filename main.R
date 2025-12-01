library(foreach)
library(doParallel)
library(doRNG)

t1 = Sys.time()

source('ogm_dahl.R')
clnum<-detectCores() 

#cl <- makeCluster(getOption("cl.cores", clnum-1))
cl <- makeCluster(getOption("cl.cores", 100),outfile="Log.txt")

registerDoParallel(cl)

registerDoRNG()

sim1 = 100

ogm.out_rea <- foreach(xun1 = 1:sim1, .options.RNG = 10) %dopar% {
  result <- try({
    ogm(xun1 = xun1)
  }, silent = TRUE)
  
  
  if (inherits(result, "try-error")) {
    return(NULL)
  } else {
    return(result)
  }
}

stopCluster(cl)

print(ogm.out_rea)
ogm.out <- Filter(Negate(is.null), ogm.out_rea)
lens=length(ogm.out);lens
result=(matrix(data.frame(ogm.out),byrow =TRUE,ncol=4))

true_clu_prob=mean(unlist(result[,1]));true_clu_prob
ARI =mean(unlist(result[,3]));ARI
RMSE =sqrt(mean(unlist(result[,4])));RMSE
  

t2 = Sys.time()
print(t2-t1)

