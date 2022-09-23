library(sde)

##### get normalized function ##### 
normalize = function(ft, delta) {
  norm_f = sqrt(sum(ft^2)*delta)
  res = ft/norm_f
  res
} 

########## Generate component (U, V) functions ##########
##### Time component U(t) #####
#### smooth function ####
comp_u <- function(nt=101, type="sine", lcomp=1) {
  times = seq(0,nt-1)/(nt-1)
  res = rep(0, nt)
  if (type == "sine") {
    res = sqrt(2) * sin(2*lcomp*pi*times) # normalized sine (cosine) fucntion
  } else {
    res = sqrt(2) * cos(2*lcomp*pi*times)
  }
  res
}

##### Wavelength component V(s) #####
#### Non-smooth bump function with selected bumps ####
comp_v <- function(ns=256, bumploc, signal=7, rsnr=7, is_noise=F, is_normalize=F) {
  x = seq(0,ns-1)/(ns-1)
  nbump = length(bumploc)
  
  w = c(0.005, 0.005, 0.006, 0.01, 0.01, 0.02, 0.02, 0.01, 
        0.005, 0.008, 0.005, 0.01, 0.008)
  h2 = c(4, 5, 3, 4, 5, 4.2, 2.1, 4.3, 3.1, 5.1, 4.2, 2.1, 3)
  bump = rep(0, ns)
  for (i in 1:nbump) {
    bump = bump + h2[i]*pmax(0, (1-abs((x-bumploc[i])/w[i])))^4
  }
  
  if (is_normalize) {
    bump = bump/sqrt(var(bump))*signal
  } 
  
  if (is_noise) {
    bump = bump + rnorm(ns, 0, sd=signal/rsnr)
  }
  
  list(bumps=bump, log_bumps=log(bump+5))
}


########## Generate regression coefficient (\beta^u, \beta^v) functions ##########
#### beta_u(t): cosine function ####
gaussprocess <- function(from = 0, to = 1, K = function(s, t) {min(s, t)},
                         start = 0, nt = 101) {
  
  t <- seq(0,nt-1)/(nt-1)
  Sigma <- sapply(t, function(s1) {
    sapply(t, function(s2) {
      K(s1, s2)
    })
  })
  
  path <- mvrnorm(mu = rep(0, times = nt), Sigma = Sigma)
  path <- path - path[1] + start  # Must always start at "start"
  
  return(data.frame("t" = t, "xt" = path))
}

#### beta_u(t): smooth gaussian process ####
beta_u <- function(nt=101) {
  gaussdf = gaussprocess(nt=nt, K=function(s,t) exp(-100*(s-t)^2))
  res = gaussdf$xt
  res
}

#### beta_v(s): standard brownian bridge ####
beta_v <- function(ns=256) {
  res = sde::BBridge(N=ns-1)
  res = as.numeric(res)
  res
}
