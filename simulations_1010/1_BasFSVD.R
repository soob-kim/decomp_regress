library(wavethresh) # wavelet basis
library(fda)
library(parallel) # for faster computing
library(far) # for Gram-schmidt, but too slow!

orth_bu <- function(evalarg, rangeval = NULL, nb){
  if(is.null(rangeval)){
    rangeval <- c(evalarg[1], evalarg[length(evalarg)])
  }
  B <- getbasismatrix(
    evalarg, 
    create.bspline.basis(rangeval, nbasis = nb)
  )
  OB <- orthonormalization(B)[,1:nb]
  return(OB)
}

basis_svd <- function(data, nb_by = 5, ncmp = 5){
  #### data of size ntime * nwv(=1024)
  ntime = dim(data)[1]
  nwv = dim(data)[2]
  
  nb = round(ntime / nb_by)
  Bu = orth_bu(1:ntime, nb = nb)
  Bv = t(GenW(n = nwv)) #vanishing mmt = 10 by default
  
  X_tilde = t(Bu) %*% data %*% Bv
  
  ## rank one svd on X_tilde
  tmp <- svd(X_tilde)
  u = Bu %*% tmp$u[,1:ncmp]
  v = Bv %*% tmp$v[,1:ncmp]
  u_coef = tmp$u[,1:ncmp]
  v_coef = tmp$v[,1:ncmp]
  return(list(u = u, v = v, d = tmp$d,
              u_coef = u_coef, v_coef = v_coef))
}




