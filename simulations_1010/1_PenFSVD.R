library(parallel) # faster computing
library(assist)   # find smoother matrix

svd_reg_u2 <- function(data, s_u, i){
  X_tilde <- s_u %*% data   # smooth data once. no big difference empirically..
  tmp <- svd(X_tilde)
  return(list(u = tmp$u[,i], v = tmp$v[,i]))
}

tw_svd_smu <- function(data, ncomp = 5){
  plain_svd <- svd(data)
  res <- list(
    pen_u = matrix(0, nrow = nrow(data), ncol = ncomp),
    pen_v = matrix(0, nrow = ncol(data), ncol = ncomp),
    d = plain_svd$d)
  for(i in 1:ncomp){
    ssr_u <- assist::ssr(u ~ times,
                         data = data.frame(u = plain_svd$u[, i], 
                                           times = 1:nrow(data)),
                         scale = T, rk = cubic(times))
    s_u <- assist::hat.ssr(ssr_u) 
    tmp <- svd_reg_u2(data, s_u, i)
    res$pen_u[, i] <- tmp$u
    res$pen_v[, i] <- tmp$v
  }
  return(res)
}
