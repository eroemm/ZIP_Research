## specify response vector, predictor of Poisson, and predictor of gamma ##
#Eventually like to get a formula function so we can just specify a data matrix, and write a y ~ x.1 +x.2 + ... | z.1 + z.2 + ...
EM.poisson <- function(y,B.0,G.0,tol = 1e-08 , max.iter = 1000, verb = F)
{
  n <- length(y)
  
  ## Posterior Probabilites
  z <- rep(NA,n)
  
  ##collection of coefficients through each run ##
  beta <- rep(NA,ncol(B.0)+1)
  gamma <- rep(NA,ncol(G.0)+1)
  
  ## Taking initial values to be beta's for the positive poisson likelihood, and logistic regression of zero - to no zero count as intial gamma ##
  beta <- summary(glm(subset(y, y>0) ~ . , family = poisson, data = subset(B.0,y>0)))$"coefficients"[,1]
  gamma <- summary(glm(as.numeric(y==0) ~ . , family = binomial, data = G.0))$"coefficients"[,1]

  ## Creating B and G matrix to calcuate needed quantities ##
  B <- as.matrix(cbind(rep(1,length(y)),B.0))
  G <- as.matrix(cbind(rep(1,length(y)),G.0))
  
  ## log.likelihood iterations
  
  
  k = 1
  log.like.new <- 2
  log.like <- 1
  while(  abs((log.like.new-log.like)/log.like) > tol &&  k <= max.iter    )
  {
    
    z <- sapply(1:n , function(i) ifelse(y[i] > 0,0,(1 + exp(-t(G[i,])%*%gamma - exp(t(B[i,])%*%beta)))^-1))
    
    log.like =  sum((1-z)*(y*(B%*%beta) - exp(B%*%beta))) + sum(z*(G%*%gamma) - log(1+exp(G%*%gamma))) - sum((1-z)*lfactorial(y))
    
    
    ## M step for Beta ##
    ##Lambert says you can use log-linear poisson with weights of 1 - Z(k) ##
    beta.new <- summary(glm(y ~ ., weights = 1 - z, family = poisson, data = B.0))$"coefficients"[,1]
    
    ## Complete log-likelihood associated with beta as defined in Lambert page 4 ##
    log.like.beta <- sum((1-z)*(y*(B%*%beta.new) - exp(B%*%beta.new)))
    ##M step for Gamma ## ## Lamber says we can use a weighted logistic regression ##
    gamma.new <- suppressWarnings(summary(glm(z ~ . , family = binomial, data = G.0))$"coefficients"[,1])
    
    ## Complete log-likelihood for gamma ##
    log.like.gamma <- sum(z*(G%*%gamma.new) - log(1+exp(G%*%gamma.new)))
    
    ##Overall log-likelihood as derived in Lambert
    log.like.new <- log.like.beta + log.like.gamma - sum((1-z)*lfactorial(y))
    
    beta <- beta.new
    gamma <- gamma.new
    
    
    if(verb ==T) {
      cat("iteration =", k , "log-lik =", log.like.new, "\n", "Beta.estimates",beta , "\n", "Gamma.estimates", gamma , "\n")
    }
    k = k + 1 
  }
  
  if( k == max.iter) cat("WARNING! NOT CONVERGENT!")
  
  ## Give them names ##
   
   l <- list(Beta.coefficients = beta,Gamma.coefficients = gamma,Posterior.Probabilitites = z, log.likelihood.iterations = log.like)
   return(l)
}








