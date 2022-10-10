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

beta_u2 = function(nt=101, l=1) {
  tv = 0:(nt-1)/(nt-1)
  res = -3+8*exp(-4*(tv-0.125*l)^2)
  res
}

#### beta_v(s): standard brownian bridge ####
beta_v <- function(ns=256) {
  res = sde::BBridge(N=ns-1)
  res = as.numeric(res)
  res
}

beta_v2 <- function(ns=256, l=1) {
  tv = 0:(ns-1)/(ns-1)
  res = cos(10*pi*tv/l)
  res 
}

#############################################
########### Generate N functions ############
#############################################

gen_u <- function(nT=101, ncomp=5) {
  U_functions = matrix(0, nrow=nT, ncol=ncomp) 
  for (j in seq(1, ncomp, by=2)) {
    U_functions[,j] = comp_u(nt=nT, type="sine", lcomp=j) 
  }
  for (j in seq(2, ncomp-1, by=2)) {
    U_functions[,j] = comp_u(nt=nT, type="cosine", lcomp=j)
  }
  U_functions
}  

gen_v_13bumps <- function(nS=256, ncomp=5, bumpnum) {
  V_functions = matrix(0, nrow=ncomp, ncol=nS)
  nbump = 13
  loc_tot = seq(0.05, 0.95, length=nbump) 
  loc_ind = split(sample(1:nbump), rep(1:ncomp, bumpnum)) # random bump locations and numbers for each component 
  loc_list = list()
  for (j in 1:5) {
    loc_list[[j]] = loc_tot[sort(loc_ind[[j]])]
  }
  
  for (j in 1:ncomp) {
    bumpfun = comp_v(ns=nS, bumploc=loc_list[[j]], is_noise=F, is_normalize = F)$bumps
    V_functions[j,] = normalize(bumpfun, delta=1/(nS-1)) # normalized
  }
  V_functions
}

plain_svd <- function(data, ncmp=5) {
  svd_res = svd(data)
  res = list(
    u = svd_res$u[,1:ncmp],
    v = svd_res$v[,1:ncmp],
    d = svd_res$d
  )
  res
}

gen_xy = function(Umat, Vmat, beta_u, beta_v, N=300, d, 
                  idx_ucomp=NULL, idx_vcomp=NULL, noise_u, noise_v, SNR) {
  ncomp = ncol(Umat) 
  nT = nrow(Umat)
  nS = ncol(Vmat)
  
  sim_oes_comps = list(u = vector("list", length = ncomp),
                       v = vector("list", length = ncomp))
  for (j in 1:ncomp) {
    sim_oes_comps$u[[j]] = matrix(0, nrow=N, ncol=nT)
    sim_oes_comps$v[[j]] = matrix(0, nrow=N, ncol=nS)
  }
  
  #### X_tot: Simulated Two-way functional data (N X nT X nS) 
  X_tot = array(0, dim=c(N, nT, nS))
  d = d
  #### Add random noise to component of each functional observations. ex) U_i(t) <- U(t) + e_i(t)
  for (i in 1:N) {
    U_tot = matrix(0, nrow=nT, ncol=ncomp)
    for (j in 1:ncomp) {
      sim_oes_comps$u[[j]][i,] = Umat[,j] + rnorm(nT, 0, sqrt(noise_u))
      U_tot[,j] = sim_oes_comps$u[[j]][i,]
    }
    V_tot = matrix(0, nrow=ncomp, ncol=nS)
    for (j in 1:ncomp) {
      sim_oes_comps$v[[j]][i,] = Vmat[j,] + rnorm(nS, 0, sqrt(noise_v))
      V_tot[j,] = sim_oes_comps$v[[j]][i,]
    }
    X_tot[i,,] = U_tot %*% diag(d) %*% V_tot 
  }
  
  #### y = f(X) + epsilon, epsilon의 표준편차는 SNR로부터 결정
  if (is.null(idx_ucomp)) {
    sim_fx_time = rep(0, N)
  } else {
    sim_fx_time = rowSums(sapply(idx_ucomp, function(x) sim_oes_comps$u[[x]] %*% beta_u[x,] * (1/(nT-1))))
  }
  
  if (is.null(idx_vcomp)) {
    sim_fx_wave = rep(0, N)
  } else {
    sim_fx_wave = rowSums(sapply(idx_vcomp, function(x) sim_oes_comps$v[[x]] %*% beta_v[x,] * (1/(nS-1))))
  }

  sim_fx = sim_fx_time + sim_fx_wave
  var_noise = var(sim_fx)/SNR
  sim_error = rnorm(N, 0, sqrt(var_noise))
  sim_y = sim_fx + sim_error
  
  res = list(X=X_tot, compfun=sim_oes_comps, y=sim_y, fx=sim_fx)
  return(res)
}

get_peak_wv <- function(data, span = 15, npeaks = 10){
  nobs = dim(data)[1]
  nt = dim(data)[2]
  nwv = dim(data)[3] 
  compressed_data <- data.frame(matrix(NA, nrow = nobs, ncol = 2*npeaks))
  data_tmp <- array(dim = c(nobs, nt, npeaks))
  for(i in 1:nobs){
    for(t in 1:nt){
      peaks_allind <- ggpmisc:::find_peaks(data[i, t, ], span = span)
      data_tmp[i,t,] <- sort(data[i, t, peaks_allind], decreasing = T)[1:npeaks]
    }
    compressed_data[i,1:npeaks] <- apply(data_tmp[i,,], 2, mean, na.rm=TRUE)
    compressed_data[i,npeaks+1:npeaks] <- apply(data_tmp[i,,], 2, sd, na.rm=TRUE)
  }
  names(compressed_data) <- c(paste(rep("m", npeaks), 1:npeaks, sep = ""), 
                              paste(rep("sd", npeaks), 1:npeaks, sep = ""))
  return(compressed_data)
}




