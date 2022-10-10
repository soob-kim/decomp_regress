library(quadprog)
library(fda.usc)
library(wavethresh)

fit01 <- function(x){
  #### Make y fit into [0,1]
  res <- x
  if(any(x<0)){
    res[which(x<0)] <- 0
  }
  if(any(x>1)){
    res[which(x>1)] <- 1
  }
  return(res)
}

#### Bspline Basis ####

jk_bsp_errmat <- function(X, y_, x_nb, b_nb,
                          tr = train_ind, tst = test_ind, verbose = F){
  #### X : list of length p(=number of variables), each is nobs*argvals matrix
  #### y_: response variable
  #### x_nb, b_nb : nbasis for x and beta, respectively. each is vector of length p
  nobs <- length(y_)
  p <- length(X)
  
  # Find error matrix
  error_mat <- matrix(0, nrow = length(tr), ncol = p)
  for(j in 1:p){
    npoints = dim(X[[j]])[2]
    X_fd <- fdata(mdata = X[[j]], argvals=0:(npoints-1)/(npoints-1))
    for(i in 1:length(tr)){
      jn_ind <- tr[setdiff(1:length(tr), i)] ### jackknife index
      X_fd_tmp <- X_fd[jn_ind]
      tmp_fit <- fregre.basis(
        fdataobj = X_fd_tmp, y = y_[jn_ind],
        basis.x = create.bspline.basis(rangeval = X_fd$rangeval, nbasis = x_nb[j]),
        basis.b = create.bspline.basis(rangeval = X_fd$rangeval, nbasis = b_nb[j]),
        lambda = 0)
      tmp_pred <- predict(tmp_fit, X_fd[tr[i]])
      #tmp_pred <- fit01(tmp_pred)
      error_mat[i, j] <- y_[tr[i]] - tmp_pred
    }
    if(verbose){
      cat("Model", j, "done\n")
    }
  }
  return(error_mat)
}

jk_bsp <- function(X, y_, x_nb, b_nb,
                   tr = train_ind, tst = test_ind, verbose = F){
  nobs <- length(y_)
  p <- length(X)
  
  res <- list()
  # Find error matrix
  res$error_mat <- jk_bsp_errmat(X = X, y_ = y_, x_nb = x_nb, b_nb = b_nb,
                                 tr = tr, tst = tst, verbose = verbose)
  # Find weight (QP)
  res$weight <- solve.QP(
    Dmat = t(res$error_mat) %*% res$error_mat / length(tr),
    dvec = rep(0, p),
    Amat = t(rbind(rep(1, p), diag(1, nrow = p, ncol = p))),
    bvec = c(1, rep(0, p)),
    meq = 1
  )
  
  res$y_m <- matrix(0, nrow = length(tst), ncol = p)
  for(j in 1:p){
    npoints = dim(X[[j]])[2]
    X_fd <- fdata(mdata = X[[j]], argvals=0:(npoints-1)/(npoints-1))
    X_fd_tr <- X_fd[tr]
    tmp_fit <- fregre.basis(
      fdataobj = X_fd_tr, y = y_[tr],
      basis.x = create.bspline.basis(rangeval = X_fd$rangeval, nbasis = x_nb[j]),
      basis.b = create.bspline.basis(rangeval = X_fd$rangeval, nbasis = b_nb[j]),
      lambda = 0)
    tmp_pred <- predict(tmp_fit, X_fd[tst])
    res$y_m[, j] <- tmp_pred
    #res$y_m[, j] <- fit01(tmp_pred)
  }
  res$y_jpa <- res$y_m %*% res$weight$solution
  res$rmse <- sqrt(mean((res$y_jpa - y_[tst])^2))
  res$mae <- mean(abs(res$y_jpa - y_[tst]))
  res$r2 <- 1 - (sum((y_[tst] - res$y_jpa)^2)/sum((y_[tst] - mean(y_[tst]))^2))
  return(res)
}

#### Wavelet Basis ####
decomp1d<-function(x, family="DaubLeAsymm", filter.number=4, 
                   bc="periodic", min.scale=2){
  nsample<-length(x[,1])
  N<-length(x[1,])
  wdt<-matrix(0, nrow=nsample, ncol=N)
  
  wds<-apply(x,1,wd, filter.number=filter.number, min.scale=min.scale, family=family, bc=bc)
  
  for(i in 1:nsample){
    base<-0
    for(j in ((log2(N)-1):min.scale)){
      temp<-accessD(wds[[i]], level=j, boundary=F)
      wdt[i,((base+1):(base+2^j))]<-temp
      base<-base+2^j
    }
    wdt[i,((N-2^min.scale+1):N)]<-accessC(wds[[i]], level=min.scale, boundary=F)
  }
  return(wdt)
}

###### Screening by mean Magnitude
scr_mag <- function(x, thres_q = 0.75, thres_var = 0.95, use_var = F){
  abs_mag <- abs(apply(x, 2, mean))
  x_var <- apply(x, 2, var)
  mag_order <- order(abs_mag, decreasing = T)
  var_cumsum <- cumsum(x_var[mag_order]) / sum(x_var)
  if(use_var){
    end_ord_ind <- which(var_cumsum >= thres_var)[1]
  } else {
    end_ord_ind <- round(dim(x)[2] * thres_q)
  }
  ind <- mag_order[1:end_ord_ind]
  var_explained <- var_cumsum[end_ord_ind]
  return(list(screened_x = x[, ind], ind = ind, var_explained = var_explained))
}


# jk_wv(X = oes_penfsvd$v, y_ = y$Y, q = 0.01,
#       tr = 1:50, tst = 51:55)

jk_wv_errmat <- function(X, y_, q = 0.1,
                         tr = train_ind, tst = test_ind, verbose = F){
  #### X : list of length p(=number of variables), each is nobs*argvals matrix
  #### y_: response variable
  #### x_nb, b_nb : nbasis for x and beta, respectively. each is vector of length p
  #### q: screening threshold
  nobs <- length(y_)
  p <- length(X)
  
  # Find error matrix
  error_mat <- matrix(0, nrow = length(tr), ncol = p)
  for(j in 1:p){
    tmp_x <- decomp1d(x = X[[j]], min.scale = 2)
    for(i in 1:length(tr)){
      jn_ind <- tr[setdiff(1:length(tr), i)] ### jackknife index
      tmp_scx <- scr_mag(tmp_x[jn_ind,], thres_q = q)
      tmp_fit <- tmp_fit <- lm(
        y~., data = data.frame(y = y_[jn_ind], tmp_scx$screened_x))
      tmp_pred <- predict(
        tmp_fit, newdata = data.frame(matrix(tmp_x[tr[i], tmp_scx$ind], nrow=1)))
      #tmp_pred <- fit01(tmp_pred)
      error_mat[i, j] <- y_[tr[i]] - tmp_pred
    }
    if(verbose){
      cat("Model", j, "done\n")
    }
  }
  return(error_mat)
}

jk_wv <- function(X, y_, q = 0.1,
                  tr = train_ind, tst = test_ind, verbose = F){
  #### X : list of length p(=number of variables), each is nobs*argvals matrix
  #### y_: response variable
  #### x_nb, b_nb : nbasis for x and beta, respectively. each is vector of length p
  #### q: screening threshold
  nobs <- length(y_)
  p <- length(X)
  
  res <- list()
  # Find error matrix
  res$error_mat <- jk_wv_errmat(X = X, y_ = y_, q = q,
                                tr = tr, tst = tst, verbose = verbose)
  
  # Find weight (QP)
  res$weight <- solve.QP(
    Dmat = t(res$error_mat) %*% res$error_mat / length(tr),
    dvec = rep(0, p),
    Amat = t(rbind(rep(1, p), diag(1, nrow = p, ncol = p))),
    bvec = c(1, rep(0, p)),
    meq = 1
  )
  
  res$y_m <- matrix(0, nrow = length(tst), ncol = p)
  for(j in 1:p){
    tmp_x <- decomp1d(x = X[[j]], min.scale = 2)
    tmp_scx <- scr_mag(tmp_x[tr,], thres_q = q)
    tmp_fit <- lm(y~., data = data.frame(y = y_[tr], tmp_scx$screened_x))
    tmp_pred <- predict(tmp_fit, newdata = data.frame(matrix(tmp_x[tst, tmp_scx$ind], nrow=length(tst))))
    res$y_m[, j] <- tmp_pred
    #res$y_m[, j] <- fit01(tmp_pred)
  }
  res$y_jpa <- res$y_m %*% res$weight$solution
  res$rmse <- sqrt(mean((res$y_jpa - y_[tst])^2))
  res$mae <- mean(abs(res$y_jpa - y_[tst]))
  res$r2 <- 1 - (sum((y_[tst] - res$y_jpa)^2)/sum((y_[tst] - mean(y_[tst]))^2))
  return(res)
}

jk_wvlasso_errmat <- function(X, y_, tr=train_ind, tst=test_ind, verbose=F) {
  nobs = length(y_)
  p = length(X)
  
  # find error matrix
  error_mat = matrix(0, nrow=length(tr), ncol=p)
  for (j in 1:p) {
    tmp_x = decomp1d(x = X[[j]], min.scale = 2)
    for (i in 1:length(tr)) {
      jn_ind = tr[setdiff(1:length(tr), i)] ### jackknife index
      tmp_fit = cv.glmnet(x=tmp_x[jn_ind,], y=y_[jn_ind], nfolds = 5, lambda=2^seq(-10,10, by=1))
      tmp_pred = as.numeric(predict(tmp_fit, newx=tmp_x[tr[i],], s=tmp_fit$lambda.min))
      error_mat[i,j] = y_[tr[i]] - tmp_pred
    }
    if(verbose){
      cat("Model", j, "done\n")
    }
  }
  error_mat
}

jk_wvlasso <- function(X, y_, 
                  tr = train_ind, tst = test_ind, verbose = F){
  #### X : list of length p(=number of variables), each is nobs*argvals matrix
  #### y_: response variable
  #### x_nb, b_nb : nbasis for x and beta, respectively. each is vector of length p
  #### q: screening threshold
  nobs <- length(y_)
  p <- length(X)
  
  res <- list()
  # Find error matrix
  res$error_mat <- jk_wvlasso_errmat(X = X, y_ = y_,
                                tr = tr, tst = tst, verbose = verbose)
  
  # Find weight (QP)
  res$weight <- solve.QP(
    Dmat = t(res$error_mat) %*% res$error_mat / length(tr),
    dvec = rep(0, p),
    Amat = t(rbind(rep(1, p), diag(1, nrow = p, ncol = p))),
    bvec = c(1, rep(0, p)),
    meq = 1
  )
  
  res$y_m <- matrix(0, nrow = length(tst), ncol = p)
  for(j in 1:p){
    tmp_x <- decomp1d(x = X[[j]], min.scale = 2)
    tmp_fit = cv.glmnet(x=tmp_x[tr,], y=y_[tr], nfolds = 5, lambda=2^seq(-10,10, by=1))
    tmp_pred = as.numeric(predict(tmp_fit, newx=tmp_x[tst,], s=tmp_fit$lambda.min))
    res$y_m[,j] = tmp_pred
  }
  res$y_jpa <- res$y_m %*% res$weight$solution
  res$rmse <- sqrt(mean((res$y_jpa - y_[tst])^2))
  res$mae <- mean(abs(res$y_jpa - y_[tst]))
  res$r2 <- 1 - (sum((y_[tst] - res$y_jpa)^2)/sum((y_[tst] - mean(y_[tst]))^2))
  return(res)
}

#### Merging models ####
jk_merge <- function(error_mat1, error_mat2, y_m1, y_m2, 
                     y_, tr, tst){
  res <- list()
  res$error_mat <- cbind(error_mat1, error_mat2)
  p <- ncol(res$error_mat)
  
  # Find weight (QP)
  res$weight <- solve.QP(
    Dmat = t(res$error_mat) %*% res$error_mat / length(tr),
    dvec = rep(0, p),
    Amat = t(rbind(rep(1, p), diag(1, nrow = p, ncol = p))),
    bvec = c(1, rep(0, p)),
    meq = 1
  )
  
  res$y_m <- cbind(y_m1, y_m2)
  res$y_jpa <- res$y_m %*% res$weight$solution
  res$rmse <- sqrt(mean((res$y_jpa - y_[tst])^2))
  res$mae <- mean(abs(res$y_jpa - y_[tst]))
  res$r2 <- 1 - (sum((y_[tst] - res$y_jpa)^2)/sum((y_[tst] - mean(y_[tst]))^2))
  return(res)
}

# jk_merge(error_mat1 = oes_res$PenBsp$error_mat[,1:10], error_mat2 = oes_res$PenBsp$error_mat[,30:50],
#          y_m1 = oes_res$PenBsp$y_m[,1:10], y_m2 = oes_res$PenBsp$y_m[,30:50], 
#          y_ = y$Y, tr = train_ind, tst = test_ind )
