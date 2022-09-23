##### package loading #####
library(ggplot2)
library(dplyr)
library(fda)
library(fda.usc)
library(wavethresh)
library(sde)
library(glmnet)
library(pls) 
library(ggpmisc)

##### source file loading ####
source("sim_generate.R")
source("1_PenFSVD.R")
source("1_BasFSVD.R")
source("2_jma_sim.R")

#################################################################

##### data dimension (size) settings #####
#### Generate 10 component functions (5 time, 5 wavelength) ####
N = 200 # number of two-way functional observations
nT = 101; nS = 256 # number of time points (0:100/100) / wavelength points (0:255/255)
nT_comp = nS_comp = 5 # number of components

##### Generate coefficient functions #####
sim_beta_u = matrix(0, nrow=nT_comp, ncol=nT) 
sim_beta_v = matrix(0, nrow=nS_comp, ncol=nS)
set.seed(2022)
for (j in 1:nT_comp) {
  sim_beta_u[j,] = beta_u(nt=nT)
  sim_beta_v[j,] = beta_v(ns=nS)
}

##############################################
################# Simulations ################
##############################################

##### functions for generating X #####
#### 시뮬레이션 할 때마다 조금씩 달라짐 (random noise) #####
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

gen_xy = function(Umat, Vmat, beta_u, beta_v, N=300, noise_x=0.01, SNR=2) {
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
  d = 2^(ncomp:1)
  #### Add random noise to component of each functional observations. ex) U_i(t) <- U(t) + e_i(t)
  for (i in 1:N) {
    U_tot = matrix(0, nrow=nT, ncol=ncomp)
    for (j in 1:ncomp) {
      sim_oes_comps$u[[j]][i,] = Umat[,j] + rnorm(nT, 0, sqrt(noise_x))
      U_tot[,j] = sim_oes_comps$u[[j]][i,]
    }
    V_tot = matrix(0, nrow=ncomp, ncol=nS)
    for (j in 1:ncomp) {
      sim_oes_comps$v[[j]][i,] = Vmat[j,] + rnorm(nS, 0, sqrt(noise_x))
      V_tot[j,] = sim_oes_comps$v[[j]][i,]
    }
    X_tot[i,,] = U_tot %*% diag(d) %*% V_tot 
  }
  
  #### y = f(X) + epsilon, epsilon의 표준편차는 SNR로부터 결정 
  sim_fx_time = rowSums(sapply(1:ncomp, function(x) sim_oes_comps$u[[x]] %*% sim_beta_u[x,] * (1/(nT-1))))
  sim_fx_wave = rowSums(sapply(1:ncomp, function(x) sim_oes_comps$v[[x]] %*% sim_beta_v[x,] * (1/(nS-1))))
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


##############################################################################


get_simulation <- function(R=100, N=200, beta_u, beta_v, noise_x=0.01, SNR=2) {
  fsvd_res = matrix(0, nrow=R, ncol=12) # (RMSE, MAE, R2) X 4 methods
  colnames(fsvd_res) = paste0(rep(c("RMSE", "MAE", "R2"), 4), "_", rep(1:4, each=3))
  time_avg_bsp_res = matrix(0, nrow=R, ncol=3) 
  time_avg_wv_res = matrix(0, nrow=R, ncol=3) 
  multi_pc_reg_res = matrix(0, nrow=R, ncol=3) 
  multi_pc_lasso_res = matrix(0, nrow=R, ncol=3) 
  multi_pls_res = matrix(0, nrow=R, ncol=3) 
  
  pre = Sys.time()
  for (r in 1:R) {
    #### training-test split (train : test = 2:1)####
    sim_trainind = sample(1:N, round(N*2/3))
    sim_testind = (1:N)[-sim_trainind]
    
    ##### 1. Decomposition+Regression #####
    ##### Generate (X_1, y_1), ... , (X_N, y_N) #####
    Umat = gen_u(nT=nT, ncomp=nT_comp)
    Vmat = gen_v_13bumps(nS=nS, ncomp=nS_comp, bumpnum = c(3,2,3,2,3))
    generated_XY = gen_xy(Umat=Umat, Vmat=Vmat, beta_u=beta_u, beta_v=beta_v, noise_x=noise_x, N=N, SNR=SNR)
    X_tot = generated_XY$X
    sim_y = generated_XY$y
    sim_fx = generated_XY$fx
    sim_oes_comps = generated_XY$compfun
    
    ##### PenFSVD on X_{i}(t_j, s_k) #####
    sim_tw_res = apply(X_tot, 1, tw_svd_smu, ncomp=nT_comp)
    sim_oes_penfsvd = list(u = vector("list", length = nT_comp),
                           v = vector("list", length = nS_comp))
    for (j in 1:nT_comp) {
      for (i in 1:N) {
        sim_oes_penfsvd$u[[j]] = rbind(sim_oes_penfsvd$u[[j]], sim_tw_res[[i]]$pen_u[,j])
        sim_oes_penfsvd$v[[j]] = rbind(sim_oes_penfsvd$v[[j]], sim_tw_res[[i]]$pen_v[,j])
      }
    }
    # flip 
    flipind = c(80, 80, 90, 90, 95)
    flipind_u = list()
    for (j in 1:nT_comp) {
      flipind_u[[j]] = which(sim_oes_penfsvd$u[[j]][,flipind[j]] > 0)
    }
    sim_oes_penfsvd_fliped = sim_oes_penfsvd
    for (j in 1:nT_comp) {
      for (i in flipind_u[[j]]) {
        sim_oes_penfsvd_fliped$u[[j]][i,] = (-1)*sim_oes_penfsvd_fliped$u[[j]][i,]
        sim_oes_penfsvd_fliped$v[[j]][i,] = (-1)*sim_oes_penfsvd_fliped$v[[j]][i,]
      }
    }
    sim_oes_penfsvd = sim_oes_penfsvd_fliped
    
    ##### BasFSVD on X_{i}(t_j, s_k) #####
    sim_basis_raw = apply(X_tot, 1, basis_svd, ncmp=nT_comp)
    sim_oes_basfsvd = list(u = vector("list", length = nT_comp),
                           v = vector("list", length = nS_comp))
    
    for (j in 1:nS_comp) {
      for (i in 1:N) {
        sim_oes_basfsvd$u[[j]] = rbind(sim_oes_basfsvd$u[[j]], sim_basis_raw[[i]]$u[,j])
        sim_oes_basfsvd$v[[j]] = rbind(sim_oes_basfsvd$v[[j]], sim_basis_raw[[i]]$v[,j])
      }
    }
    # flip
    flipind = c(80, 80, 90, 90, 95)
    flipind_u = list()
    for (j in 1:nT_comp) {
      flipind_u[[j]] = which(sim_oes_basfsvd$u[[j]][,flipind[j]] > 0)
    }
    sim_oes_basfsvd_fliped = sim_oes_basfsvd
    for (j in 1:nT_comp) {
      for (i in flipind_u[[j]]) {
        sim_oes_basfsvd_fliped$u[[j]][i,] = (-1)*sim_oes_basfsvd_fliped$u[[j]][i,]
        sim_oes_basfsvd_fliped$v[[j]][i,] = (-1)*sim_oes_basfsvd_fliped$v[[j]][i,]
      }
    }
    sim_oes_basfsvd = sim_oes_basfsvd_fliped

    ###############################################
    ##### Estimate weights from training data #####
    ###############################################

    #### Step 1 - Bsp(Time) + Bsp (Wavelength)
    #### PenFSVD + Bsp (Time) & Bsp (Wavelength) ####
    sim_oes_res = list()
    sim_oes_res$PenBsp = jk_bsp(
      X = c(sim_oes_penfsvd$u, sim_oes_penfsvd$v), y_=sim_y,
      x_nb = c(rep(10, nT_comp), rep(26, nS_comp)),
      b_nb = c(rep(10, nT_comp), rep(26, nS_comp)),
      tr = sim_trainind, tst = sim_testind
    )
    fsvd_res[r, 1:3] = c(sqrt(mean((sim_oes_res$PenBsp$y_jpa[,1]-sim_fx[sim_testind])^2)),
                         mean(abs(sim_oes_res$PenBsp$y_jpa[,1]-sim_fx[sim_testind])),
                         1 - (sum((sim_fx[sim_testind] - sim_oes_res$PenBsp$y_jpa[,1])^2)/sum((sim_fx[sim_testind] - mean(sim_fx[sim_testind]))^2)))
    
    #### BasFSVD + Bsp (Time) & Bsp (Wavelength) ####
    sim_oes_res$BasBsp = jk_bsp(
      X = c(sim_oes_basfsvd$u, sim_oes_basfsvd$v), y_=sim_y,
      x_nb = c(rep(10, nT_comp), rep(26, nS_comp)),
      b_nb = c(rep(10, nT_comp), rep(26, nS_comp)),
      tr = sim_trainind, tst = sim_testind
    )
    fsvd_res[r, 4:6] = c(sqrt(mean((sim_oes_res$BasBsp$y_jpa[,1]-sim_fx[sim_testind])^2)),
                         mean(abs(sim_oes_res$BasBsp$y_jpa[,1]-sim_fx[sim_testind])),
                         1 - (sum((sim_fx[sim_testind] - sim_oes_res$BasBsp$y_jpa[,1])^2)/sum((sim_fx[sim_testind] - mean(sim_fx[sim_testind]))^2)))
    
    #### Step 2 - Bsp(Time) + Wv(Wavelength) ####
    #### Pen + Wv 
    sim_oes_res$PenWv150 <- jk_wv(
      X = sim_oes_penfsvd$v, y_ = sim_y, q = 0.15,
      tr = sim_trainind, tst = sim_testind
    )
    sim_oes_res$PenMerged150 <- jk_merge(
      error_mat1 = sim_oes_res$PenBsp$error_mat[, 1:nT_comp],
      error_mat2 = sim_oes_res$PenWv150$error_mat,
      y_m1 = sim_oes_res$PenBsp$y_m[, 1:nT_comp],
      y_m2 = sim_oes_res$PenWv150$y_m,
      y_ = sim_y, tr = sim_trainind, tst = sim_testind
    )
    fsvd_res[r, 7:9] = c(sqrt(mean((sim_oes_res$PenMerged150$y_jpa[,1]-sim_fx[sim_testind])^2)),
                         mean(abs(sim_oes_res$PenMerged150$y_jpa[,1]-sim_fx[sim_testind])),
                         1 - (sum((sim_fx[sim_testind] - sim_oes_res$PenMerged150$y_jpa[,1])^2)/sum((sim_fx[sim_testind] - mean(sim_fx[sim_testind]))^2)))
    
    #### Bas + Wv
    sim_oes_res$BasWv150 <- jk_wv(
      X = sim_oes_basfsvd$v, y_ = sim_y, q = 0.15,
      tr = sim_trainind, tst = sim_testind
    )
    sim_oes_res$BasMerged150 <- jk_merge(
      error_mat1 = sim_oes_res$BasBsp$error_mat[, 1:nT_comp], 
      error_mat2 = sim_oes_res$BasWv150$error_mat,
      y_m1 = sim_oes_res$BasBsp$y_m[, 1:nT_comp], 
      y_m2 = sim_oes_res$BasWv150$y_m,
      y_ = sim_y, tr = sim_trainind, tst = sim_testind
    )
    fsvd_res[r, 10:12] = c(sqrt(mean((sim_oes_res$BasMerged150$y_jpa[,1]-sim_fx[sim_testind])^2)),
                           mean(abs(sim_oes_res$BasMerged150$y_jpa[,1]-sim_fx[sim_testind])),
                           1 - (sum((sim_fx[sim_testind] - sim_oes_res$BasMerged150$y_jpa[,1])^2)/sum((sim_fx[sim_testind] - mean(sim_fx[sim_testind]))^2)))
    
    ########## Compare with other methods ##########
    ##### 2. Time-averaged method #####
    #### B-spline ####
    X_time_avg = matrix(0, nrow=N, ncol=nS)
    for (i in 1:N) {
      X_time_avg[i,] = colMeans(X_tot[i,,])
    }
    
    X_time_avg_fd = fdata(X_time_avg, argvals = 0:(nS-1)/(nS-1))
    X_time_avg_tr = X_time_avg_fd[sim_trainind]
    sim_y_tr = sim_y[sim_trainind]
    X_time_avg_fit = fregre.basis(
      fdataobj = X_time_avg_tr,
      y = sim_y_tr,
      basis.x = create.bspline.basis(rangeval=X_time_avg_tr$rangeval, nbasis=26),
      basis.b = create.bspline.basis(rangeval=X_time_avg_tr$rangeval, nbasis=26), 
      lambda=0
    )
    y_avg_pred = as.numeric(predict(X_time_avg_fit, X_time_avg_fd[sim_testind]))
    time_avg_bsp_res[r,1:3] = c(sqrt(mean((y_avg_pred-sim_fx[sim_testind])^2)),
                               mean(abs(y_avg_pred-sim_fx[sim_testind])),
                               1 - (sum((sim_fx[sim_testind] - y_avg_pred)^2)/sum((sim_fx[sim_testind] - mean(sim_fx[sim_testind]))^2)))
    
    ### Wavelet ####
    wv_X_time_avg = decomp1d(X_time_avg)
    sim_y_tr = sim_y[sim_trainind]
    wv_scx_tr = scr_mag(wv_X_time_avg[sim_trainind,], thres_q = 0.15)
    wv_scx_fit = lm(y ~., data=data.frame(y=sim_y_tr, wv_scx_tr$screened_x))
    wv_y_pred = predict(wv_scx_fit, newdata = data.frame(matrix(wv_X_time_avg[sim_testind, wv_scx_tr$ind], nrow=length(sim_testind))))
    wv_y_pred = as.numeric(wv_y_pred)
    time_avg_wv_res[r,1:3] = c(sqrt(mean((wv_y_pred-sim_fx[sim_testind])^2)),
                              mean(abs(wv_y_pred-sim_fx[sim_testind])),
                              1 - (sum((sim_fx[sim_testind] - wv_y_pred)^2)/sum((sim_fx[sim_testind] - mean(sim_fx[sim_testind]))^2)))
    
    
    #### 3. Multivariate approach ####
    #### data compressing using mean and sd of selected peaks ####
    sim_peak_wv = get_peak_wv(X_tot, span=15, npeak=10)
    sim_pca_train = princomp(sim_peak_wv[sim_trainind,])
    sim_peakcomp_train = cbind(sim_pca_train$scores[,1:13]) # 98% 
    colnames(sim_peakcomp_train) = paste0("PC", 1:13)
    sim_peakcomp_test = predict(sim_pca_train, newdata=sim_peak_wv[sim_testind,])[,1:13]
    colnames(sim_peakcomp_test) = paste0("PC", 1:13)
    
    #### PC reg
    sim_pca_fit = lm(y ~ ., data=data.frame(y=sim_y[sim_trainind], sim_peakcomp_train))
    sim_pca_pred = predict(sim_pca_fit, data.frame(sim_peakcomp_test))
    multi_pc_reg_res[r,1:3] = c(sqrt(mean((sim_pca_pred-sim_fx[sim_testind])^2)),
                               mean(abs(sim_pca_pred-sim_fx[sim_testind])),
                               1 - (sum((sim_fx[sim_testind] - sim_pca_pred)^2)/sum((sim_fx[sim_testind] - mean(sim_fx[sim_testind]))^2)))
    
    ##### PC+LASSO
    sim_pca_lassofit = glmnet(x=sim_peakcomp_train, y=sim_y[sim_trainind], alpha=1)
    sim_cv_lasso = cv.glmnet(x=sim_peakcomp_train, y=sim_y[sim_trainind], alpha=1, grid=sim_pca_lassofit$lambda)
    sim_pca_lassopred = predict(sim_pca_lassofit, newx=sim_peakcomp_test, s=sim_cv_lasso$lambda.min)
    sim_pca_lassopred = as.numeric(sim_pca_lassopred)
    multi_pc_lasso_res[r,1:3] = c(sqrt(mean((sim_pca_lassopred-sim_fx[sim_testind])^2)),
                                 mean(abs(sim_pca_lassopred-sim_fx[sim_testind])),
                                 1 - (sum((sim_fx[sim_testind] - sim_pca_lassopred)^2)/sum((sim_fx[sim_testind] - mean(sim_fx[sim_testind]))^2)))
    ##### PLS using compressed data
    sim_cv_pls = pls::plsr(y ~ ., data=data.frame(y=sim_y[sim_trainind], sim_peak_wv[sim_trainind,]), ncomp=20, validation="CV")
    sim_pls_pred = predict(sim_cv_pls, ncomp=1, newdata=data.frame(y=sim_y[sim_testind], sim_peak_wv[sim_testind,]))
    sim_pls_pred = as.numeric(sim_pls_pred)
    multi_pls_res[r,1:3] = c(sqrt(mean((sim_pls_pred-sim_fx[sim_testind])^2)),
                            mean(abs(sim_pls_pred-sim_fx[sim_testind])),
                            1 - (sum((sim_fx[sim_testind] - sim_pls_pred)^2)/sum((sim_fx[sim_testind] - mean(sim_fx[sim_testind]))^2)))
  
    print(paste("rep", r, "fin"))
  }
  print(Sys.time()-pre)
  res = list(fsvd=fsvd_res, time_avg_bsp=time_avg_bsp_res, time_avg_wv=time_avg_wv_res, 
             multi_pc_reg=multi_pc_reg_res, multi_pc_lasso = multi_pc_lasso_res, multi_pls=multi_pls_res)
  
  res 
}

### noise_x = 0.01, SNR=10 ###
sim_res = get_simulation(R=100, N=200, beta_u=sim_beta_u, beta_v=sim_beta_v, noise_x=0.01, SNR=10)
sim_res_rmse = cbind(sim_res$fsvd[,c(1,4,7,10)], sim_res$time_avg_bsp[,1], sim_res$time_avg_wv[,1],
                     sim_res$multi_pc_reg[,1], sim_res$multi_pc_lasso[,1], sim_res$multi_pls[,1])
colnames(sim_res_rmse) = c("PenBsp", "BasBsp", "PenWv", "BasWv", "AvgBsp", "AvgWv",
                                                    "PCreg", "PClasso", "PLS")
boxplot(sim_res_rmse, cex.axis=0.7, ylab="RMSPE",
        main=expression(sigma[x]^2 == paste(0.01, ",   SNR=10")))

### noise_x = 1, SNR=10 ###
sim_res_p2 = get_simulation(R=100, N=200, beta_u=sim_beta_u, beta_v=sim_beta_v, noise_x=1, SNR=10)
sim_res_p2_rmse = cbind(sim_res_p2$fsvd[,c(1,4,7,10)], sim_res_p2$time_avg_bsp[,1], sim_res_p2$time_avg_wv[,1],
                        sim_res_p2$multi_pc_reg[,1], sim_res_p2$multi_pc_lasso[,1], sim_res_p2$multi_pls[,1])
colnames(sim_res_p2_rmse) = c("PenBsp", "BasBsp", "PenWv", "BasWv", "AvgBsp", "AvgWv",
                           "PCreg", "PClasso", "PLS")
boxplot(sim_res_p2_rmse, cex.axis=0.7, ylab="RMSPE",
        main=expression(sigma[x]^2 == paste(1, ",   SNR=10")))

### noise_x = 0.01, SNR=2###
sim_res_p3 = get_simulation(R=100, N=200, beta_u=sim_beta_u, beta_v=sim_beta_v, noise_x=0.01, SNR=2)
sim_res_p3_rmse = cbind(sim_res_p3$fsvd[,c(1,4,7,10)], sim_res_p3$time_avg_bsp[,1], sim_res_p3$time_avg_wv[,1],
                        sim_res_p3$multi_pc_reg[,1], sim_res_p3$multi_pc_lasso[,1], sim_res_p3$multi_pls[,1])
colnames(sim_res_p3_rmse) = c("PenBsp", "BasBsp", "PenWv", "BasWv", "AvgBsp", "AvgWv",
                              "PCreg", "PClasso", "PLS")
boxplot(sim_res_p3_rmse, cex.axis=0.7, ylab="RMSPE",
        main=expression(sigma[x]^2 == paste(0.01, ",   SNR=2")))


### noise_x = 1, SNR=2 ###
sim_res_p4 = get_simulation(R=100, N=200, beta_u=sim_beta_u, beta_v=sim_beta_v, noise_x=1, SNR=2)
sim_res_p4_rmse = cbind(sim_res_p4$fsvd[,c(1,4,7,10)], sim_res_p4$time_avg_bsp[,1], sim_res_p4$time_avg_wv[,1],
                        sim_res_p4$multi_pc_reg[,1], sim_res_p4$multi_pc_lasso[,1], sim_res_p4$multi_pls[,1])
colnames(sim_res_p4_rmse) = c("PenBsp", "BasBsp", "PenWv", "BasWv", "AvgBsp", "AvgWv",
                              "PCreg", "PClasso", "PLS")
boxplot(sim_res_p4_rmse, cex.axis=0.7, ylab="RMSPE",
        main=expression(sigma[x]^2 == paste(1, ",   SNR=2")))





