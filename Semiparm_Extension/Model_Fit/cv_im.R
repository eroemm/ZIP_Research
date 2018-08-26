
## For calculating Kernel Weights at each target point, input z.i, target poiint, bandwidth, and kernel type ##
kw.cal = function(t.i,t.0, h, kernel = c("Epa","Gauss","Unif"))
{
  kernel = match.arg(kernel)
  if(kernel == "Epa") val = prod((h^-1)*epa((t.i-t.0)*(h^-1))) 
  else if(kernel == "Gauss") val = prod((h^-1)*dnorm((t.i-t.0)*(h^-1)))
  else val = prod((h^-1)*dunif( (t.i-t.0)*(h^-1) , -1 , 1))
  return(val)
}





CV = function(formula,data,h, k, kernel = c("Epa","Gauss","Unif") , delta ) ## formula, matrix of bandwidths, k -folds, kernel type##
{
  ## Extracts the relevant information from the formula statement.
  cl <- match.call()
  if(missing(data)) data <- environment(formula)
  mf <- match.call(expand.dots = FALSE)
  m <- match(c("formula", "data"), names(mf), 0)
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
  
  ## Returns testing classes ##
  classes = sample(rep(1:k,length.out = length(Y)))
  ## Vector to stor CV value for each h ##
  
  result = t(sapply(1:nrow(h) , function(i) 
  {
    ## Calculate log-like for each testing fold (within bandwidth)
    ## For each fold, calculate log-like ##
   quant.k <- t(sapply(1:k , function(j)
    {
      ## Training data ##
      Y.train <- subset(Y, classes != j)  
      Z.train <- subset(Z, classes != j)
      X.train <- subset(X, classes != j)
      
      
      ## Testing Data
      Y.test <- subset(Y,classes == j)
      Z.test <- subset(Z, classes == j)
      X.test <- subset(X, classes == j)
      
      
      ## Obtain trianing posterior probabilities and training betas
      
      model.train = local.poisson.fit(Y.train,X.train,Z.train, kernel = kernel , h = h[i,], delta = delta, max.iter = 40, tol = .0001)
      
      h.emp <- model.train$bandwidth.3rd
      
      r.train = model.train$post.prob
      b.train = model.train$beta.hat
      ## have estimates of betas and pi's ## 
      ## Extrapolate pi's to those in test set ##
      ## Test data
      
      ## Kernel weights associated with testing data observations
      n.k <- length(Y.test)
      n.train <- length(Y.train)
      kw.test = sapply(1:n.k , function(ii) sapply(1:(n.train), function(jj) kw.cal(Z.train[jj,-1],Z.test[ii,-1], h = h.emp, kernel = kernel)))
      
      ## extrapolate to testing set ##
      pi.hat.test = sapply(1:nrow(kw.test), function(l) (r.train%*%t(kw.test[l,]))*(sum(kw.test[l,])^-1))
      
      CV.k  = sum(lzip(pi.hat.test,Y.test, exp(X.test%*%b.train)))
      return(c(h.emp,CV.k))
    }))
    
    final.h <- apply(as.matrix(quant.k[,-ncol(quant.k)]), 2 , max)
    CV = sum(quant.k[,ncol(quant.k)])
    return(c(final.h,CV))
  }))
  return(result)
}

