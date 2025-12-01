ogm=function(xun1){
rm(list = ls(all = TRUE))
start_timeww <- Sys.time()
library(Rcpp)
library(RcppArmadillo)
library(mvtnorm)
library("MASS")

sourceCpp("ogm_MFM.cpp")


#####
#Data generation
#####

#K=2
if(TRUE){
  #''Graphical Models for Ordinal Data'', 2015 by Jian Guo et al.
  library("MASS") 
  set.seed(xun1) 
  n <- 100 
  p=10
  {
    K=3
    
    mean_vectors <- replicate(p, sample(c(0, -1, 1), K, replace = FALSE))
    m1=mean_vectors[1,]
    m2=mean_vectors[2,]
    m3=mean_vectors[3,]
    
  }
  
  true_graph <- vector("list", 2)
  S2 <- diag(1,p,p)
  for(i in 1:p){
    for(j in 1:p){
      if(j==i+1){
        S2[i,j]=0.5
        S2[j,i]=0.5
      }
    }
  }
  adj1=ifelse(S2>0,1,0)
  diag(adj1)=0
  true_graph[[1]]=adj1
  cov_matrix <- solve(S2)
  
  
  variances <- diag(cov_matrix)
  # Scale the covariance matrix so that variances are equal to 1
  scaling_matrix <- diag(1 / sqrt(variances))
  cov_matrix_scaled <- scaling_matrix %*% cov_matrix %*% scaling_matrix
  
  clus2 <- mvrnorm(n = n, mu = m2, Sigma = cov_matrix_scaled)
  
  S3 <- diag(1,p,p)
  adj2=ifelse(S3>0,1,0)
  diag(adj2)=0
  true_graph[[2]]=adj2
  clus3 <- mvrnorm(n = n, mu = m3, Sigma = S3)
  
  Z <- rbind( clus2, clus3)
  
  
  #label
  c=c(rep(1,n),rep(2,n))
  n=nrow(Z)
  p=ncol(Z)
  
  KK=4#level of ordinal data is K-1
  Theta=matrix(0,nrow=p,ncol=KK)
  for(j in 1:p){
    
    for(k in 2:(KK-1)){
      Theta[j,k]=runif(1,qnorm((k-1-0.5)/(KK-1)),qnorm((k-1+0.5)/(KK-1)))
    }
    Theta[j,1]=-Inf
    Theta[j,KK]=Inf
    
  }
  ###get ordinal data
  X=matrix(0,nrow=n,ncol=p)
  for(i in 1:n){
    for(j in 1:p){
      X[i,j]=sum(Z[i,j]>=Theta[j,])
    }
  }
  View(X)
  
  
}
#K=3
if(FALSE){
  #''Graphical Models for Ordinal Data'', 2015 by Jian Guo et al.
  library("MASS") 
  set.seed(xun1+90) 
  n <- 100 
  p=15
  {
    K=3
    
    mean_vectors <- replicate(p, sample(c(0, -1, 1), K, replace = FALSE))
    m1=mean_vectors[1,]
    m2=mean_vectors[2,]
    m3=mean_vectors[3,]
  }
  
  true_graph <- vector("list", 3)
  
  S1 <- diag(1,p,p)
  for(i in 1:p){
    for(j in 1:p){
      if(j==i+1&j%%2==0){
        S1[i,j]=0.5
        S1[j,i]=0.5
      }
      if(j==i+2&j%%2==0){
        S1[i,j]=0.25
        S1[j,i]=0.25
      }
    }
  }
  adj=ifelse(S1>0,1,0)
  diag(adj)=0
  true_graph[[1]]=adj
  cov_matrix <- solve(S1)
  
 
  variances <- diag(cov_matrix)
  # Scale the covariance matrix so that variances are equal to 1
  scaling_matrix <- diag(1 / sqrt(variances))
  cov_matrix_scaled <- scaling_matrix %*% cov_matrix %*% scaling_matrix

  clus1 <- mvrnorm(n = n, mu = m1, Sigma = cov_matrix_scaled)
  
  
  
  S2 <- diag(1,p,p)
  for(i in 1:p){
    for(j in 1:p){
      if(j==i+1){
        S2[i,j]=0.5
        S2[j,i]=0.5
      }
    }
  }
  adj1=ifelse(S2>0,1,0)
  diag(adj1)=0
  true_graph[[2]]=adj1
  cov_matrix <- solve(S2)
  
  
  variances <- diag(cov_matrix)
 
  scaling_matrix <- diag(1 / sqrt(variances))
  cov_matrix_scaled <- scaling_matrix %*% cov_matrix %*% scaling_matrix
  
  clus2 <- mvrnorm(n = n, mu = m2, Sigma = cov_matrix_scaled)
  
  S3 <- diag(1,p,p)
  adj2=ifelse(S3>0,1,0)
  diag(adj2)=0
  true_graph[[3]]=adj2
  
  clus3 <- mvrnorm(n = n, mu = m3, Sigma = S3)
  
  Z <- rbind(clus1, clus2, clus3)
  
  
  #label
  c=c(rep(1,n),rep(2,n),rep(3,n))
  
  n=nrow(Z)
  p=ncol(Z)
  
  KK=4#level of ordinal data is K-1
  Theta=matrix(0,nrow=p,ncol=KK)
  for(j in 1:p){
    
    for(k in 2:(KK-1)){
      Theta[j,k]=runif(1,qnorm((k-1-0.5)/(KK-1)),qnorm((k-1+0.5)/(KK-1)))
    }
    Theta[j,1]=-Inf
    Theta[j,KK]=Inf
    
  }
  ###get ordinal data
  X=matrix(0,nrow=n,ncol=p)
  for(i in 1:n){
    for(j in 1:p){
      X[i,j]=sum(Z[i,j]>=Theta[j,])
    }
  }
  #View(X)
  
  
}

#unbalanced data
if(FALSE){
  #''Graphical Models for Ordinal Data'', 2015 by Jian Guo et al.
  library("MASS")
  set.seed(xun1) 
  n <- 100
  p=15
  {
    K=3
    
    mean_vectors <- replicate(p, sample(c(0, -1, 1), K, replace = FALSE))
    m1=mean_vectors[1,]
    m2=mean_vectors[2,]
    m3=mean_vectors[3,]
  }
  
  true_graph <- vector("list", 3)
  
  S1 <- diag(1,p,p)
  for(i in 1:p){
    for(j in 1:p){
      if(j==i+1){
        S1[i,j]=0.5
        S1[j,i]=0.5
      }
      if(j==i+2){
        S1[i,j]=0.25
        S1[j,i]=0.25
      }
    }
  }
  adj=ifelse(S1>0,1,0)
  diag(adj)=0
  true_graph[[1]]=adj
  cov_matrix <- solve(S1)
  
  
  variances <- diag(cov_matrix)
  
  scaling_matrix <- diag(1 / sqrt(variances))
  cov_matrix_scaled <- scaling_matrix %*% cov_matrix %*% scaling_matrix
  
  clus1 <- mvrnorm(n = n, mu = m1, Sigma = cov_matrix_scaled)
  
  
  
  S2 <- diag(1,p,p)
  for(i in 1:p){
    for(j in 1:p){
      if(j==i+1){
        S2[i,j]=0.5
        S2[j,i]=0.5
      }
    }
  }
  adj1=ifelse(S2>0,1,0)
  diag(adj1)=0
  true_graph[[2]]=adj1
  cov_matrix <- solve(S2)
  
  
  variances <- diag(cov_matrix)
  
  scaling_matrix <- diag(1 / sqrt(variances))
  cov_matrix_scaled <- scaling_matrix %*% cov_matrix %*% scaling_matrix
  
  clus2 <- mvrnorm(n = n, mu = m2, Sigma = cov_matrix_scaled)
  
  S3 <- diag(1,p,p)
  adj2=ifelse(S3>0,1,0)
  diag(adj2)=0
  true_graph[[3]]=adj2
  
  clus3 <- mvrnorm(n = n, mu = m3, Sigma = S3)
  
  Z <- rbind(clus1, clus2, clus3)
  
  
  
  c=c(rep(1,n),rep(2,n),rep(3,n))
  
  n=nrow(Z)
  p=ncol(Z)
  
  KK=4#level of ordinal data is K-1
  Theta=matrix(0,nrow=p,ncol=KK)
  for(j in 1:p){
    for(k in 2:(KK-1)){
      Theta[j,k]=runif(1,qnorm((k-1-0.5)/(KK-1)),qnorm((k-1+0.5)/(KK-1)))
    }
    Theta[j,1]=-Inf
    Theta[j,KK]=Inf
    
  }
  ###get ordinal data
  X=matrix(0,nrow=n,ncol=p)
  for(i in 1:n){
    for(j in 1:p){
      X[i,j]=sum(Z[i,j]>=Theta[j,])
    }
  }
  #View(X)
  
  
}



##K=5&different level
if(FALSE){
  #''Graphical Models for Ordinal Data'', 2015 by Jian Guo et al.
  library("MASS") 
  set.seed(xun1+70) 
  n <- 200 
  p=10
  {
    K=5
    mean_vectors <- replicate(p, sample(c(0, -1, 1), K, replace = TRUE))
    m1=mean_vectors[1,]
    m2=mean_vectors[2,]
    m3=mean_vectors[3,]
    m4=mean_vectors[4,]
    m5=mean_vectors[5,]
  }
  
  true_graph <- vector("list", 5)
  
  S1 <- diag(1,p,p)
  for(i in 1:p){
    for(j in 1:p){
      if(j==i+1){
        S1[i,j]=0.5
        S1[j,i]=0.5
      }
      if(j==i+2){
        S1[i,j]=0.25
        S1[j,i]=0.25
      }
    }
  }
  adj=ifelse(S1>0,1,0)
  diag(adj)=0
  true_graph[[1]]=adj
  cov_matrix <- solve(S1)
  
  
  variances <- diag(cov_matrix)
  # Scale the covariance matrix so that variances are equal to 1
  scaling_matrix <- diag(1 / sqrt(variances))
  cov_matrix_scaled <- scaling_matrix %*% cov_matrix %*% scaling_matrix
  
  clus1 <- mvrnorm(n = n, mu = m1, Sigma = cov_matrix_scaled)
  
  
  
  S2 <- diag(1,p,p)
  for(i in 1:p){
    for(j in 1:p){
      if(j==i+1){
        S2[i,j]=0.5
        S2[j,i]=0.5
      }
    }
  }
  adj1=ifelse(S2>0,1,0)
  diag(adj1)=0
  true_graph[[2]]=adj1
  cov_matrix <- solve(S2)
  
  
  variances <- diag(cov_matrix)
  scaling_matrix <- diag(1 / sqrt(variances))
  cov_matrix_scaled <- scaling_matrix %*% cov_matrix %*% scaling_matrix
  
  clus2 <- mvrnorm(n = n, mu = m2, Sigma = cov_matrix_scaled)
  
  S3 <- diag(1,p,p)
  adj2=ifelse(S3>0,1,0)
  diag(adj2)=0
  true_graph[[3]]=adj2
  
  clus3 <- mvrnorm(n = n, mu = m3, Sigma = S3)
  
  S4 <- diag(1,p,p)
  for(i in 1:p){
    for(j in 1:p){
      if(j==i+1&j%%2==0){
        S4[i,j]=0.5
        S4[j,i]=0.5
      }
      if(j==i+2&j%%2==0){
        S4[i,j]=0.25
        S4[j,i]=0.25
      }
    }
  }
  adj4=ifelse(S4>0,1,0)
  diag(adj4)=0
  true_graph[[4]]=adj4
  cov_matrix <- solve(S4)
  
  
  variances <- diag(cov_matrix)
  # Scale the covariance matrix so that variances are equal to 1
  scaling_matrix <- diag(1 / sqrt(variances))
  cov_matrix_scaled <- scaling_matrix %*% cov_matrix %*% scaling_matrix
  
  clus4 <- mvrnorm(n = n, mu = m4, Sigma = cov_matrix_scaled)
  
  S5 <- diag(1,p,p)
  for(i in 1:p){
    for(j in 1:p){
      if(j==i+1&j%%2==1){
        S5[i,j]=0.5
        S5[j,i]=0.5
      }
      if(j==i+2&j%%2==1){
        S5[i,j]=0.25
        S5[j,i]=0.25
      }
    }
  }
  adj5=ifelse(S5>0,1,0)
  diag(adj5)=0
  true_graph[[5]]=adj5
  cov_matrix <- solve(S5)
  
  
  variances <- diag(cov_matrix)
  
  scaling_matrix <- diag(1 / sqrt(variances))
  cov_matrix_scaled <- scaling_matrix %*% cov_matrix %*% scaling_matrix
  
  clus5 <- mvrnorm(n = n, mu = m5, Sigma = cov_matrix_scaled)
  
  
  Z <- rbind(clus1, clus2, clus3,clus4,clus5)
  
  c=c(rep(1,n),rep(2,n),rep(3,n),rep(4,n),rep(5,n))
  
  n=nrow(Z)
  p=ncol(Z)
  
  KK=5+1#level of ordinal data is K-1
  Theta=matrix(0,nrow=p,ncol=KK)
  for(j in 1:p){
    
    for(k in 2:(KK-1)){
      Theta[j,k]=runif(1,qnorm((k-1-0.5)/(KK-1)),qnorm((k-1+0.5)/(KK-1)))
    }
    Theta[j,1]=-Inf
    Theta[j,KK]=Inf
    
  }
  ###get ordinal data
  X=matrix(0,nrow=n,ncol=p)
  for(i in 1:n){
    for(j in 1:p){
      X[i,j]=sum(Z[i,j]>=Theta[j,])
    }
  }
  #View(X)
  
  
}




##Get simple parameters about the data
data=X#ordinal data
n=nrow(data)
p=ncol(data)
KK=nlevels(as.factor(data[,1]))+1


## hyper-parameter values setting 
hyper <- list()
# DP prior
hyper$diri <- 1
hyper$poiss <- 1
# VN
hyper$VN <- fn_calc_VN(n, hyper$diri, hyper$poiss)

hyper$a=0.01#for the prior on the mean
hyper$g.prior=0.2
#####
#Initialization
#####

#ini.K <- 2#
ini.K <-length(unique(c))
ini.indicator <- c
#ini.indicator <- c(sample(1:ini.K, size = ini.K, replace = FALSE),
#                   sample(1:ini.K, size = n- ini.K, replace = TRUE))
ini.sample <- NULL
ini.sample$K <- as.integer(ini.K)#the number of cluster
ini.sample$indicator <- as.integer(ini.indicator)#indicator
#ini.sample$N.k <- table(ini.z)
ini.sample$N.k <- rep(NA, ini.sample$K)
for(k in 1:ini.sample$K) ini.sample$N.k[k] <- sum(ini.sample$indicator == k)



ini.sample$Theta <-matrix(0,nrow=p,ncol=KK)
for(j in 1:p){
  for(k in 1:KK){
    if(k==1){
      ini.sample$Theta[j,k]=-Inf
    }else if(k==KK){
      ini.sample$Theta[j,k]=Inf
    }else{
      ini.sample$Theta[j,k]=qnorm(sum(data[,j]<k)/n)
    }
  }
}
ini.sample$graphs <- vector("list", ini.sample$K)
ini.sample$prec <- vector("list", ini.sample$K)
ini.sample$mu <- vector("list",ini.sample$K)


library(BDgraph)
for (k in 1:ini.sample$K) {
  G <- matrix(rbinom(p * p, size = 1, prob = hyper$g.prior), nrow = p, ncol = p)
  G[lower.tri(G)] <- t(G)[lower.tri(G)]
  diag(G) <- 0
  ini.sample$graphs[[k]] <- G
  ini.sample$prec[[k]] <-rgwish( n = 1, adj = G, b = 3, D = diag( p ) )
  ini.sample$mu[[k]]=mvrnorm(n = 1, mu = rep(0,p), Sigma = diag( p ))
}


ini.sample$latent=matrix(0,nrow=n,ncol=p)



library(TruncatedNormal)
timing <- system.time({if(TRUE){
  
  for(i in 1:n){
    lower=rep(0,p)
    upper=rep(0,p)
    for(j in 1:p){
      lower[j]=max(ini.sample$Theta[j,X[i,j]],-100)
      upper[j]=min(ini.sample$Theta[j,X[i,j]+1],100)
    }
    
   
    ini.sample$latent[i,]=rtmvnorm(n=1, mu=ini.sample$mu[[ini.indicator[i]]],
                                   sigma=solve(ini.sample$prec[[ini.indicator[i]]]), lb=lower, ub=upper)
  }
}})


cur.sample=ini.sample

n.MCMC.iter <- 2000
n.burnin <- 1000
save.result<-list()
save.result$z <- array(dim = c(n.MCMC.iter, n))
save.result$K <- rep(NA, n.MCMC.iter)
save.result$N.k <- vector("list", n.MCMC.iter)
save.result$theta <- vector("list", n.MCMC.iter)
save.result$mu <- vector("list", n.MCMC.iter)
save.result$graphs <- vector("list", n.MCMC.iter)
save.result$prec <- vector("list", n.MCMC.iter)

##Gibbs
for(i.iter in 1:n.MCMC.iter){
  
  
  #Updata (\mu,\Omega,G)
  for(k in 1:cur.sample$K){
    burnin=100
    
    iter=burnin+1
    
    g.start = "empty"
    g_prior = BDgraph::get_g_prior( g.prior = hyper$g.prior, p = p )
    G       = BDgraph::get_g_start( g.start = g.start, g_prior = hyper$g_prior, p = p )
    
    
    b      = 3
    
    b_star = b + cur.sample$N.k[k]+1
    D      = diag( p )
    
    SS1=sweep(matrix(cur.sample$latent[cur.sample$indicator == k,],ncol=p), 2, as.vector(cur.sample$mu[[k]]), "-")
    SS2=t(SS1)%*%SS1
    S=SS2+hyper$a*as.vector(cur.sample$mu[[k]])%*%t(as.vector(cur.sample$mu[[k]]))
    Ds     = D + S
    Ts     = chol( solve( Ds ) )
    threshold = 1e-8
    K_hat      = matrix( 0, p, p )
    p_links = matrix( 0, p, p )
    trace_mcmc=iter + 1000
    K_graph       = BDgraph::get_K_start( G = G, g.start = g.start, Ts = Ts,
                                          b_star = b_star, threshold = threshold )
    
    result=.C( "ggm_bdmcmc_ma", as.integer(iter), as.integer(burnin), 
               G = as.integer(G), as.double(g_prior), as.double(Ts), 
               K = as.double(K_graph), as.integer(p), as.double(threshold), 
               K_hat = as.double(K_hat), p_links = as.double(p_links),
               as.integer(b), as.integer(b_star), as.double(Ds), 
               as.integer(trace_mcmc), PACKAGE = "BDgraph" )
    
    cur.sample$graphs[[k]]=matrix( result $ G    , p, p)
    cur.sample$prec[[k]]= matrix( result $ K    , p, p )
    
    mu_tem=(hyper$a+cur.sample$N.k[k])^{-1}*(colSums(matrix(cur.sample$latent[cur.sample$indicator==k,],ncol=p)))
    
    Prec_tem=(hyper$a+cur.sample$N.k[k])*cur.sample$prec[[k]]
    L <- chol(Prec_tem) 
    sigma <- solve(L) %*% t(solve(L)) 
    cur.sample$mu[[k]]=mvrnorm(n = 1, mu = mu_tem, Sigma = sigma)
  }
  
  
  
  
  #updata latent Z
  for(i in 1:n){
    lower=rep(0,p)
    upper=rep(0,p)
    for(j in 1:p){
      
      lower[j]=max(cur.sample$Theta[j,X[i,j]],-100)#>=-100
      upper[j]=min(cur.sample$Theta[j,X[i,j]+1],100)#<=100
      
    }
    
    
    L <- chol(cur.sample$prec[[cur.sample$indicator[i]]])
    sigma <- solve(L) %*% t(solve(L))
    
    cur.sample$latent[i,]=rtmvnorm(n=1, mu=cur.sample$mu[[cur.sample$indicator[i]]],
                                   sigma=sigma, lb=lower, ub=upper)
    latent_i_mu=cur.sample$mu[[cur.sample$indicator[i]]]
    latent_i_cov=solve(cur.sample$prec[[cur.sample$indicator[i]]])
    
    if (any(is.infinite(cur.sample$latent[i,])) || any(is.na(cur.sample$latent[i,]))) {
      
      latent_sample <- numeric(length(lower))
      for (m in 1:length(lower)) {
        if (is.na(cur.sample$latent[i,m])|is.finite(is.na(cur.sample$latent[i,m])))  {
          if(latent_i_mu[m]>=lower[m]&&latent_i_mu[m]<=upper[m]){
            latent_sample[m]=latent_i_mu[m]
          } else if(latent_i_mu[m]<lower[m]){
            latent_sample[m] <- runif(1, min = lower[m], max = min(upper[m],lower[m]+1)) 
          } else if(latent_i_mu[m]>upper[m]){
            latent_sample[m] <- runif(1, min = max(lower[m],upper[m]-1), max = upper[m])
          }
        } else{
          latent_sample[m]=cur.sample$latent[i,m]
        }
        
      }
      cur.sample$latent[i,]=latent_sample
      if (any(is.infinite(cur.sample$latent[i,])) || any(is.na(cur.sample$latent[i,]))) {
        stop("error: cur.sample$latent[i,] contains Inf or NaN values.")
      }
    }
    }
  
  
  #updata theta
  for(j in 1:p){
    for(kk in 2:(KK-1)){
      min_tem=max(max(cur.sample$latent[which(X[,j]==kk-1),j]),cur.sample$Theta[j,kk-1])
      max_tem=min(min(cur.sample$latent[which(X[,j]==kk),j]),cur.sample$Theta[j,kk+1])
      cur.sample$Theta[j,kk]=runif(1,min=min_tem,max=max_tem)
      
    }
    cur.sample$Theta[j,1]=-Inf
    cur.sample$Theta[j,KK]=Inf
  }
  
  
  
  
  
  
  
  
  # Update indicator/cluster label--the final step-"Cpp"
  #In R, I define .sample$N.k as a length K vector
  #But in Rcpp, I define N.k as a n vector which is its max dim.
  #So I fill the rest of it as 0. 
  N.k <- c(cur.sample$N.k, rep(as.integer(0), (n - cur.sample$K)))
  
  zero_graphs=replicate(n - cur.sample$K, matrix(0, nrow = p, ncol = p), simplify = FALSE)
  one_precisions=replicate(n - cur.sample$K, diag(1, nrow = p, ncol = p), simplify = FALSE)
  graphs<-c(cur.sample$graphs, zero_graphs)
  prec<-c(cur.sample$prec,one_precisions)
  
  zero_vectors <- replicate(n - cur.sample$K, rep(0, p), simplify = FALSE)
  mu<-c(cur.sample$mu,zero_vectors)
  
  
  
  fn_z_update_mod(cur.sample$latent, hyper$diri,hyper$VN, 
                  cur.sample$indicator, cur.sample$K, N.k, 
                  graphs,prec,mu, 
                  n,hyper$g.prior,hyper$a)
  
  
  cur.sample$N.k <- N.k[1:cur.sample$K]
  
  cur.sample$graphs <- graphs[1:cur.sample$K]
  cur.sample$prec<-prec[1:cur.sample$K]
  cur.sample$mu<-mu[1:cur.sample$K]
  
  
  
  
  
  
  
  save.result$z[i.iter, ] <- cur.sample$indicator
  save.result$K[i.iter] <- cur.sample$K
  save.result$N.k[[i.iter]] <- cur.sample$N.k
  save.result$theta[[i.iter]]<-cur.sample$Theta
  
  save.result$mu[[i.iter]] <- cur.sample$mu
  save.result$graphs[[i.iter]] <- cur.sample$graphs
  save.result$prec[[i.iter]]<-cur.sample$prec
  
  
  cat(" iteration:", i.iter,"\n")
  
  cat(" So far, take the time:", Sys.time()-start_timeww,"\n")
  
}








#####
## Dahl's method to summarize the samples from the MCMC
GetDahl <- function(MFMfit, n.burn.in)
{
  ################################################################
  
  ## Input: MFMfit = the result from MFM ##
  ##        n.burn.in = the number of burn-in iterations ##
  
  ## Output: 
  ##         Dahl.res = estimated output ##
  
  #################################################################
  z <- MFMfit$z[-(1:n.burn.in),]
  
  n.iter <- dim(z)[1]
  n.grid <- dim(z)[2]
  sum.matrix <- matrix(0, nrow = n.grid, ncol = n.grid)
  for(i in 1:n.iter){
    membership.matrix <- outer(z[i,], z[i,], "==")
    sum.matrix <- sum.matrix + membership.matrix
  }
  membership.avg <- sum.matrix / n.iter
  min.sq.error = n.grid^2
  for(i in 1:n.iter){
    membership.matrix <- outer(z[i,], z[i,], "==")
    sq.error <- sum((membership.matrix - membership.avg)^2)
    if(sq.error < min.sq.error){
      min.sq.error <- sq.error
      Dahl.index <- i
    }
  }
  iter.index <- n.burn.in + Dahl.index
  Dahl.res <- list(z = MFMfit$z[iter.index,], K = MFMfit$K[iter.index], 
                   N.k = MFMfit$N.k[[iter.index]], #lambda = MFMfit$lambda[[iter.index]]
                   mu=MFMfit$mu[[iter.index]],
                   theta=MFMfit$theta[[iter.index]],
                   graphs=MFMfit$graphs[[iter.index]],
                   prec=MFMfit$prec[[iter.index]]
                   )
  attr(Dahl.res, "iterIndex") <- iter.index
  attr(Dahl.res, "burnin") <- n.burn.in
  return(Dahl.res)
}
#####
#
#####
Dahl_train <- GetDahl(save.result, n.burnin)
print(Dahl_train)









end_time <- Sys.time()


execution_time <- end_time - start_timeww
print(execution_time)

#measure 1
true_clu_num=length(unique(c))==Dahl_train$K
clu_num=Dahl_train$K

#measure 2
library(mclust)
ARI=adjustedRandIndex(c, Dahl_train$z)

#measure 3：graph
score=0
n=nrow(data)
for(i in 1:n){
  score=score+sum((true_graph[[c[i]]]-Dahl_train$graphs[[Dahl_train$z[i]]])^2)
}
score=score/n

return(list(
  true_clu_num=true_clu_num,
  clu_num=clu_num,
  ARI=ARI,
  RMSE=score
))
}