


## Epa kernel ##
epa = function(t)
{
  ifelse(abs(t)<= 1, .75*(1-t^2),0)
}

## Completed log-likelihood ##
lczip = function(r,p,y,lambda)
{
  r*log(p) + (1-r)*log(1-p) + (1-r)*dpois(y,lambda = lambda, log = TRUE)
}

## log like for zip
lzip = function(p,y,lambda)
{
  ifelse(y == 0 , log(p + (1-p)*dpois(y, lambda = lambda)), log((1-p)*dpois(y,lambda = lambda)))
}



## After extracting response vector, X, Z matrices, fit model using given kernel and bandwidth ##
local.poisson.fit = function(Y,X,Z,kernel = c("Epa","Gauss","Unif"), h, tol,max.iter, delta) 
{
  
h.input <- h
  
n <- length(Y)
kx <- NCOL(X)
kz <- NCOL(Z)  
  
  
kernel = match.arg(kernel)
## calculate kernel weights
kw <- matrix(NA,nrow = n, ncol = n)
h.test = TRUE
while(h.test)
{
  kw = sapply(1:n , function(i) sapply(1:n, function(j) kw.cal(Z[j,-1],Z[i,-1], h = h, kernel = kernel)))
  
  zfit <- EM.poisson(Y,as.data.frame(X[,-1]),as.data.frame(Z[,-1]))
  
  ##initial values for 1st local steps ##
  btilde = as.matrix(zfit$Beta.coefficients)
  btilde.z = matrix(rep(btilde,n),nrow = n, ncol = kx, byrow = TRUE)
  
  ##initial values for pi's
  pi.tilde = sapply(1:n, function(i)  (1 + exp(- Z[i,]%*%zfit$Gamma.coefficients))^-1)
  
  ##local likelihood 1 ##
  k = 1
  ll.loc.prev.z = rep(1,n)
  ll.loc.new.z = rep(2,n)
  
  while(norm(as.matrix(ll.loc.new.z - ll.loc.prev.z )) > tol && k <= max.iter)
  {
    ##E step
    r = sapply(1:n,function(i) ifelse(Y[i] == 0, pi.tilde[i]*(pi.tilde[i]+ (1-pi.tilde[i])*dpois(0,exp(X[i,]%*%t(btilde.z[i,]))))^-1, 0))
    
    ll.loc.prev.z = sapply(1:n, function(i) lczip(r,pi.tilde,Y,exp(X%*%btilde.z[i,]))%*%kw[i,])
    
    ##M Step 
  
    pi.tilde = sapply(1:n,function(i) (sum(kw[i,])^-1)*(kw[i,]%*%r))
    
    pi.tilde = ifelse(pi.tilde == 0 , min(subset(pi.tilde, pi.tilde > 0 )), ifelse(pi.tilde == 1, max(subset(pi.tilde, pi.tilde < 1 )), pi.tilde))
  
    
    
    btilde.z = t(sapply(1:n, function(i) tryCatch(glm(Y ~ X[,-1] , family = 'poisson', weights = kw[i,]*(1-r))$coefficients, error = function(e){return(rep(NA,kx))})))
    
    ll.loc.new.z = sapply(1:n, function(i) lczip(r,pi.tilde,Y,exp(X%*%btilde.z[i,]))%*%kw[i,])
    
    if(any(is.na(btilde.z)) == TRUE || any(is.na(ll.loc.new.z)) == TRUE)
    {
      k = max.iter + 1
      h = h + delta
    
    }else { k = k + 1}

    

    
 
    
  }
  if(any(is.na(btilde.z)) == TRUE || any(is.na(ll.loc.new.z)) == TRUE ) h.test = TRUE else h.test = FALSE
}


##global likelihood step 2 ##
a = 1
ll.global.prev = 1
ll.global.new = 2
b.hat <- colMeans(btilde.z)
while(norm(as.matrix(ll.global.new - ll.global.prev)) >tol && a <= max.iter)
{
  r = sapply(1:n,function(i) ifelse(Y[i] == 0, pi.tilde[i]*(pi.tilde[i]+ (1-pi.tilde[i])*dpois(0,exp(X[i,]%*%b.hat)))^-1, 0))
  
  ll.global.prev = sum(lczip(r,pi.tilde,Y,exp(X%*%b.hat)))
  
  b.hat = glm(Y ~ X[,-1], family = 'poisson', weights = 1 - r)$coefficients
  
  ll.global.new = sum(lczip(r,pi.tilde,Y,exp(X%*%b.hat)))
  
  a = a + 1
  
}

## final local maximization step 3 ##
c = 1
ll.loc.3 <- rep(1,n)
ll.loc.3.new <- rep(2,n)
pi.hat <- pi.tilde


kw = sapply(1:n , function(i) sapply(1:n, function(j) kw.cal(Z[j,-1],Z[i,-1], h = h.input, kernel = kernel)))

## Issue of having estimated pi.hats equal to one. In this case, do the same thing as in step 1. ##

while(norm(as.matrix(ll.loc.3.new - ll.loc.3))> tol && c <= max.iter)
{
  r = sapply(1:n,function(i) ifelse(Y[i] == 0, pi.hat[i]*(pi.hat[i]+ (1-pi.hat[i])*dpois(0,exp(X[i,]%*%b.hat)))^-1, 0))
  
  ll.loc.3 <- sapply(1:n, function(i) lczip(r,pi.hat,Y,exp(X%*%b.hat))%*%kw[i,])
  
  pi.hat = sapply(1:n,function(i) (sum(kw[i,])^-1)*(kw[i,]%*%r))
  impute = sum(pi.hat == 0 || pi.hat == 1)
  pi.hat = ifelse(pi.hat == 0 , min(subset(pi.hat, pi.hat > 0 )), ifelse(pi.hat == 1, max(subset(pi.hat, pi.hat < 1 )), pi.hat))
  
  ll.loc.3.new = sapply(1:n, function(i) lczip(r,pi.hat,Y,exp(X%*%b.hat))%*%kw[i,]) 
  
  
}

log.prob <- lzip(pi.hat,Y,exp(X%*%b.hat))

ll <- sum(lzip(pi.hat,Y,exp(X%*%b.hat)))


all.estimates <- list(loglik = ll, log.prob = log.prob ,pi.tilde = pi.tilde, pi.hat = pi.hat , beta.hat = b.hat, post.prob = r, bandwidth.first = h, bandwidth.3rd = h.input, percent.impute = impute/n)
return(all.estimates)
}






## formula of y ~ x1 + x2 + ... | z1 + z2 + ... , data, vector of bandwidths , kernel type , tolerance, and max iterations ##
semizip <- function(formula, data,h, delta, kernel = c("Epa","Gauss","Unif") , tol = 1e-8,max.iter = 40) 
{ 
  ## Extracts the relevant information from the formula statement.
  cl <- match.call()
  if(missing(data)) data <- environment(formula)
  mf <- match.call(expand.dots = FALSE)
  m <- match(c("formula", "data","subset"), names(mf), 0)
  mf <- mf[c(1, m)]
  mf$drop.unused.levels <- TRUE
  if (length(formula[[3]]) > 1 && identical(formula[[3]][[1]], 
                                            as.name("|"))) {
    ff <- formula
    formula[[3]][1] <- call("+")
    mf$formula <- formula
    ffc <- . ~ .
    ffz <- ~.
    ffc[[2]] <- ff[[2]]
    ffc[[3]] <- ff[[3]][[2]]
    ffz[[3]] <- ff[[3]][[3]]
    ffz[[2]] <- NULL
  }
  else {
    ffz <- ffc <- ff <- formula
    ffz[[2]] <- NULL
  }
  if (inherits(try(terms(ffz), silent = TRUE), "try-error")) {
    ffz <- eval(parse(text = sprintf(paste("%s -", deparse(ffc[[2]])), 
                                     deparse(ffz))))
  }
  mf[[1]] <- as.name("model.frame")
  mf <- eval(mf, parent.frame())
  mt <- attr(mf, "terms")
  mtX <- terms(ffc, data = data)
  X <- model.matrix(mtX, mf)
  mtZ <- terms(ffz, data = data)
  mtZ <- terms(update(mtZ, ~.), data = data)
  Z <- model.matrix(mtZ, mf)
  Y <- model.response(mf, "numeric")
  n <- length(Y)
  kx <- NCOL(X)
  kz <- NCOL(Z)
  
  kernel = match.arg(kernel)
  
  
return(local.poisson.fit(Y,X,Z,kernel = kernel, h = h, tol = tol, max.iter = max.iter, delta = delta))  

}
















