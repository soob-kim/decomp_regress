# ------------------------------------------------------------------------------
# -------- R-Code Smooth Group Lasso for Functional Regression -----------------
# ------------------------------------------------------------------------------


require(splines)
require(grplasso)

# Multiple Functional Linear Model

grplFlinear <- function(Y, X, Tps, lambda, phi, dfs = 10,
adapt1 = NULL, adapt2 = NULL, ...){

	#### Observed functions x_{ij}(t) at possibly different support grids t
	## Y = nx1 vector of responses
	## X = list of matrices, each element corresponds to observations of one functional predictor.
	##		Thus, [[jj]] would provide the jj-th observed function with subjects as rows and the support grid as columns
	##		Note that each function (i.e. component of X) can have different support
	## Tps = list of vectors, each element is a suppport grid of the corresponding observed function.
	##		Thus, Tps[[jj]] would give the support grid of the jj-th observed function. This grid is
	##		assumed to be the same for all subjects ii, and equidistant
	## lambda = vector of penalty parameters
	## phi = penalty parameter for smoothing
	## dfs = dfs used for basis expansions of coefficient functions (can be a vector)
	## adapt1, adapt2 = adaptive weights for the coefficient functions and second derivatives, respectively.
	##     (if given, vectors with length matching the number of predictors)


	nsub = length(Y) ## number of subjects
	nfunc = length(Tps) ## number of functions per subject



	#### We use bsplines as basis functions for the corrsponding beta functions
  if (length(dfs) == 1)
    dfs = rep(dfs, nfunc) ## vector of intended df of each spline basis
  if (length(dfs) != nfunc)
    stop("length of dfs does not match number of predictors")

	B <- Psi <- Omega <- K <- iR <- eK <- list()
	delt <- rep(NA, nfunc)


	for (jj in 1:nfunc){
    spj = diff(range(Tps[[jj]]))#/(dfs[jj]-2)
    bknj = c(min(Tps[[jj]]) - spj, max(Tps[[jj]]) + spj) ## boundary knots
		B[[jj]] = bs(Tps[[jj]], df=dfs[jj], Boundary.knots=bknj) ## basis spline set up
		delt[jj] = Tps[[jj]][2] - Tps[[jj]][1] ## differences in Tps

		Psi[[jj]] = delt[jj] * t(B[[jj]]) %*% B[[jj]] ## approximate norm of bsplines assuming dense design
    if (length(adapt1) == nfunc)
      Psi[[jj]] = adapt1[jj]*Psi[[jj]]

    dBj <- matrix(NA,nrow(B[[jj]]),ncol(B[[jj]]))
    for (k in 1:ncol(B[[jj]])) ## computation of second derivatives
      {
        iS <- interpSpline(Tps[[jj]],B[[jj]][,k])
        dBj[,k] <- predict(iS, Tps[[jj]], deriv = 2)$y
      }
    Omega[[jj]] = delt[jj] * t(dBj) %*% dBj ## approximate norm of 2nd deriv of bsplines assuming dense design
    if (length(adapt2) == nfunc)
      Omega[[jj]] = adapt2[jj]*Omega[[jj]]

		K[[jj]] = Psi[[jj]] + phi * Omega[[jj]] ## K matrix
    eK[[jj]] <- eigen(K[[jj]])
    #iR[[jj]] <- t((1/sqrt(eK[[jj]]$values))*t(eK[[jj]]$vectors))
		iR[[jj]] = backsolve(chol(K[[jj]]), x = diag(ncol(K[[jj]])))  ## inverse of cholesky of K
	}


	## covariates for the linear model
	Z = 1
	for (jj in 1:nfunc)
    {
			tmp = delt[jj]*(X[[jj]]%*%B[[jj]])
			Z = cbind(Z, tmp%*%iR[[jj]])
		}


	## group lasso
	index = c(NA,rep(1:nfunc,dfs))
	grpl = grplasso(x = Z, y = Y, index = index, model = LinReg(), lambda = lambda,
  standardize = F, ...)


  ## output: intercept and fitted coefficient functions
  intercept = grpl$coef[1,]
  Coef <- list()
  index[1] = 0
  for (jj in 1:nfunc)
    {
      Coef[[jj]] <- B[[jj]]%*%iR[[jj]]%*%grpl$coef[index == jj,]
    }

  out = list("intercept" = intercept, "Coef" = Coef)
	return(out)
}



cv.grplFlinear <- function(k, Y, X, Tps, lambda, phi, dfs = 10, ...)
{
  # K-fold cross-validation function
  n <- length(Y)
  p <- length(X)
  os <- sample(n,n)
  
  cvError <- array(NA, c(k,length(phi),length(lambda)))
  nk <- floor(n/k)
  for (wk in 1:k)
  {
    cat("fold",wk,"\n")
    if (wk < k)
      inds <- os[(wk-1)*nk+(1:nk)]
    else
      inds <- os[-(1:((wk-1)*nk))]
    
    # training/test obs
    Xk <- X
    Xkt <- X
    delt = rep(NA, p)
    for (jj in 1:p)
    {
      Xk[[jj]] <- X[[jj]][-inds,]
      Xkt[[jj]] <- X[[jj]][inds,]
      delt[jj] = Tps[[jj]][2] - Tps[[jj]][1]
    }
    Yk <- Y[-inds]
    Ykt <- Y[inds]
    
    # phi walues
    for (wp in 1:length(phi))
    {
      # model fitting
      grplk <- grplFlinear(Y = Yk, X = Xk, Tps = Tps, lambda = lambda, 
                        phi = phi[wp], dfs = dfs, ...)
      
      # test set prediction
      predk <- matrix(NA,length(inds),length(lambda))
      for (ll in 1:length(lambda))
      {
        predll <- rep(grplk$intercept[ll],nrow(predk))
        for (jj in 1:p)
        {
          predll <- predll + delt[jj]*Xkt[[jj]]%*%grplk$Coef[[jj]][,ll]
        }
        # error
        cvError[wk,wp,ll] <- sum((Ykt - predll)^2)
      }
    }
  }
  return(cvError)
}





# Logistic Regression

grplFlogit <- function(Y, X, Tps, lambda, phi, dfs = 10,
adapt1 = NULL, adapt2 = NULL, ...){

	#### Observed functions x_{ij}(t) at possibly different support grids t
	## Y = nx1 vector of responses (0/1)
	## X = list of matrices, each element corresponds to observations of one functional predictor.
	##		Thus, [[jj]] would provide the jj-th observed function with subjects as rows and the support grid as columns
	##		Note that each function can have different support
	## Tps = list of vectors, each element is a suppport grid of the corresponding observed function.
	##		Thus, Tps[[jj]] would give the support grid of the jj-th observed function. This grid is
	##		assumed to be the same for all subjects ii, and equidistant
	## lambda = vector of penalty parameters
	## phi = penalty parameter for smoothing
	## dfs = dfs used for basis expansions of coefficient functions (can be a vector)


	nsub = length(Y) ## number of subjects
	nfunc = length(Tps) ## number of functions per subject



	#### We use bsplines as basis functions for the corrsponding beta functions
  if (length(dfs) == 1)
    dfs = rep(dfs, nfunc) ## vector of intended df of each spline basis
  if (length(dfs) != nfunc)
    stop("length of dfs does not match number of predictors")

	B <- Psi <- Omega <- K <- iR <- eK <- list()
	delt <- rep(NA, nfunc)


	for (jj in 1:nfunc){
    spj = diff(range(Tps[[jj]]))#/(dfs[jj]-2)
    bknj = c(min(Tps[[jj]]) - spj, max(Tps[[jj]]) + spj) ## boundary knots
		B[[jj]] = bs(Tps[[jj]], df=dfs[jj], Boundary.knots=bknj) ## basis spline set up
		delt[jj] = Tps[[jj]][2] - Tps[[jj]][1] ## differences in Tps

		Psi[[jj]] = delt[jj] * t(B[[jj]]) %*% B[[jj]] ## approximate norm of bsplines assuming dense design
    if (length(adapt1) == nfunc)
      Psi[[jj]] = adapt1[jj]*Psi[[jj]]

    dBj <- matrix(NA,nrow(B[[jj]]),ncol(B[[jj]]))
    for (k in 1:ncol(B[[jj]])) ## computation of second derivatives
      {
        iS <- interpSpline(Tps[[jj]],B[[jj]][,k])
        dBj[,k] <- predict(iS, Tps[[jj]], deriv = 2)$y
      }
    Omega[[jj]] = delt[jj] * t(dBj) %*% dBj ## approximate norm of 2nd deriv of bsplines assuming dense design
    if (length(adapt2) == nfunc)
      Omega[[jj]] = adapt2[jj]*Omega[[jj]]

		K[[jj]] = Psi[[jj]] + phi * Omega[[jj]] ## K matrix
    eK[[jj]] <- eigen(K[[jj]])
    #iR[[jj]] <- t((1/sqrt(eK[[jj]]$values))*t(eK[[jj]]$vectors))
		iR[[jj]] = backsolve(chol(K[[jj]]), x = diag(ncol(K[[jj]])))  ## inverse of cholesky of K
	}


	## covariates for the generalized linear model
	Z = 1
	for (jj in 1:nfunc)
    {
			tmp = delt[jj]*(X[[jj]]%*%B[[jj]])
			Z = cbind(Z, tmp%*%iR[[jj]])
		}


	## group lasso
	index = c(NA,rep(1:nfunc,dfs))
	grpl = grplasso(x = Z, y = Y, index = index, model = LogReg(), lambda = lambda,
  standardize = F, ...)


  ## output: intercept and fitted coefficient functions
  intercept = grpl$coef[1,]
  Coef <- list()
  index[1] = 0
  for (jj in 1:nfunc)
    {
      Coef[[jj]] <- B[[jj]]%*%iR[[jj]]%*%grpl$coef[index == jj,]
    }

  out = list("intercept" = intercept, "Coef" = Coef)
	return(out)
}



cv.grplFlogit <- function(k, Error = "logScore", Y, X, Tps, lambda, phi,
                         dfs = 10, ...)
{
  # K-fold cross-validation function
  n <- length(Y)
  p <- length(X)
  os <- sample(n,n)
  
  cvError <- array(NA, c(k,length(phi),length(lambda)))
  nk <- floor(n/k)
  for (wk in 1:k)
  {
    cat("fold",wk,"\n")
    if (wk < k)
      inds <- os[(wk-1)*nk+(1:nk)]
    else
      inds <- os[-(1:((wk-1)*nk))]
    
    # training/test obs
    Xk <- X
    Xkt <- X
    delt = rep(NA, p)
    for (jj in 1:p)
    {
      Xk[[jj]] <- X[[jj]][-inds,]
      Xkt[[jj]] <- X[[jj]][inds,]
      delt[jj] = Tps[[jj]][2] - Tps[[jj]][1]
    }
    Yk <- Y[-inds]
    Ykt <- Y[inds]
    
    # phi values
    for (wp in 1:length(phi))
    {
      # model fitting
      grplk <- grplFlogit(Y = Yk, X = Xk, Tps = Tps, lambda = lambda,
                         phi = phi[wp], dfs = dfs, ...)
      
      # test set prediction
      predk <- matrix(NA,length(inds),length(lambda))
      for (ll in 1:length(lambda))
      {
        etall <- rep(grplk$intercept[ll],nrow(predk))
        for (jj in 1:p)
        {
          etall <- etall + delt[jj]*Xkt[[jj]]%*%grplk$Coef[[jj]][,ll]
        }
        # probs
        predll <- 1/(1+exp(-etall))
        
        # error
        if (Error == "logScore")
        {
          cvError[wk,wp,ll] <- sum(as.numeric(Ykt==1)*log(1/predll)) +
            sum(as.numeric(Ykt==0)*log(1/(1-predll)))
        }
        else if (Error == "BrierScore")
          cvError[wk,wp,ll] <- sum((Ykt - predll)^2)
        else
          stop("error type not supported")
      }
    }
  }
  return(cvError)
}



