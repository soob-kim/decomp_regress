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
### Functions for model fitting ###
source("1_PenFSVD.R")
source("1_BasFSVD.R")
source("2_jma_sim.R")
source("grplFunct.r")
### Functions for generating dataset ###
source("sim_generate_0920.R")

##############################################################################
################### Main functions for simulation ############################
##############################################################################

get_simulation <- function(R=100, N=200, beta_u, beta_v, 
                           idx_ucomp, idx_vcomp,d, noise_u, noise_v, SNR) {
  fsvd_res = matrix(0, nrow=R, ncol=18) # (RMSE, MAE, R2) X 6 methods (FSVD1~4, Plain SVD 1~2)
  colnames(fsvd_res) = paste0(rep(c("RMSE", "MAE", "R2"), 6), "_", rep(1:6, each=3))
  time_avg_bsp_res = matrix(0, nrow=R, ncol=3) 
  time_avg_wv_res = matrix(0, nrow=R, ncol=3) 
  multi_pc_reg_res = matrix(0, nrow=R, ncol=3) 
  multi_pc_lasso_res = matrix(0, nrow=R, ncol=3) 
  multi_pls_res = matrix(0, nrow=R, ncol=3) 
  gflm_res = matrix(0, nrow=R, ncol=3)
  
  nT_comp = dim(sim_beta_u)[1]
  nS_comp = dim(sim_beta_v)[1]
  nT = dim(sim_beta_u)[2]
  nS = dim(sim_beta_v)[2]

  Tps = vector("list", length=10)
  for (j in 1:5) {
    Tps[[j]] = (0:(nT-1))/(nT-1)
  }
  for (j in 6:10) {
    Tps[[j]] = (0:(nS-1))/(nS-1)
  }
  lambda = 10^seq(2,-4,by=-0.5)
  phi = 10^seq(4,-3,by=-1)
  
  pre = Sys.time()
  for (r in 1:R) {
    #### training-test split (train : test = 2:1 or 1:2)####
    sim_trainind = sample(1:N, round(N*1/3))
    sim_testind = (1:N)[-sim_trainind]
    
    ##### 1. Decomposition+Regression #####
    ##### Generate (X_1, y_1), ... , (X_N, y_N) #####
    Umat = gen_u(nT=nT, ncomp=nT_comp)
    Vmat = gen_v_13bumps(nS=nS, ncomp=nS_comp, bumpnum = c(3,2,3,2,3))
    generated_XY = gen_xy(Umat=Umat, Vmat=Vmat, beta_u=beta_u, beta_v=beta_v, d=d, 
                          idx_ucomp = idx_ucomp, idx_vcomp = idx_vcomp, noise_u=noise_u, noise_v=noise_v, N=N, SNR=SNR)
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
    
    
    ##### palin SVD on X_{i}(t_j, s_k) #####
    sim_plain_raw = apply(X_tot, 1, plain_svd, ncmp=nT_comp)
    sim_oes_plain = list(u = vector("list", length = nT_comp),
                         v = vector("list", length = nS_comp))
    
    for (j in 1:nS_comp) {
      for (i in 1:N) {
        sim_oes_plain$u[[j]] = rbind(sim_oes_plain$u[[j]], sim_plain_raw[[i]]$u[,j])
        sim_oes_plain$v[[j]] = rbind(sim_oes_plain$v[[j]], sim_plain_raw[[i]]$v[,j])
      }
    }
    # flip
    flipind = c(80, 80, 90, 90, 95)
    flipind_u = list()
    for (j in 1:nT_comp) {
      flipind_u[[j]] = which(sim_oes_plain$u[[j]][,flipind[j]] > 0)
    }
    sim_oes_plain_fliped = sim_oes_plain
    for (j in 1:nT_comp) {
      for (i in flipind_u[[j]]) {
        sim_oes_plain_fliped$u[[j]][i,] = (-1)*sim_oes_plain_fliped$u[[j]][i,]
        sim_oes_plain_fliped$v[[j]][i,] = (-1)*sim_oes_plain_fliped$v[[j]][i,]
      }
    }
    sim_oes_plain = sim_oes_plain_fliped
    
    
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
    
    ### plain SVD + Bsp(wavelength)
    sim_oes_res$plain = jk_bsp(
      X = c(sim_oes_plain$u, sim_oes_plain$v), y_=sim_y,
      x_nb = c(rep(10, nT_comp), rep(26, nS_comp)),
      b_nb = c(rep(10, nT_comp), rep(26, nS_comp)),
      tr = sim_trainind, tst = sim_testind
    )
    fsvd_res[r, 13:15] = c(sqrt(mean((sim_oes_res$plain$y_jpa[,1]-sim_fx[sim_testind])^2)),
                           mean(abs(sim_oes_res$plain$y_jpa[,1]-sim_fx[sim_testind])),
                           1 - (sum((sim_fx[sim_testind] - sim_oes_res$plain$y_jpa[,1])^2)/sum((sim_fx[sim_testind] - mean(sim_fx[sim_testind]))^2))) 
    
    ### plain SVD + Wv(wavelength)
    sim_oes_res$plainWv150 <- jk_wv(
      X = sim_oes_plain$v, y_ = sim_y, q = 0.15,
      tr = sim_trainind, tst = sim_testind
    )
    sim_oes_res$plainMerged150 <- jk_merge(
      error_mat1 = sim_oes_res$plain$error_mat[, 1:nT_comp], 
      error_mat2 = sim_oes_res$plainWv150$error_mat,
      y_m1 = sim_oes_res$plain$y_m[, 1:nT_comp], 
      y_m2 = sim_oes_res$plainWv150$y_m,
      y_ = sim_y, tr = sim_trainind, tst = sim_testind
    )
    fsvd_res[r, 16:18] = c(sqrt(mean((sim_oes_res$plainMerged150$y_jpa[,1]-sim_fx[sim_testind])^2)),
                           mean(abs(sim_oes_res$plainMerged150$y_jpa[,1]-sim_fx[sim_testind])),
                           1 - (sum((sim_fx[sim_testind] - sim_oes_res$plainMerged150$y_jpa[,1])^2)/sum((sim_fx[sim_testind] - mean(sim_fx[sim_testind]))^2)))
    
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
    
    ##### GFLM method (Gerthesiss et al. 2013) #####
    #### This method uses component functions estimated from FSVD (Penfsvd) ####
    penX = c(sim_oes_penfsvd$u, sim_oes_penfsvd$v)
    penX_tr = vector("list", length=10)
    for (j in 1:10) {
      penX_tr[[j]] = penX[[j]][sim_trainind,]
    }
    
    ### 5 fold CV using training set ###
    cv_grpl = cv.grplFlinear(k=5, Y=sim_y[sim_trainind], X=penX_tr, Tps=Tps, lambda=lambda, phi=phi,dfs = 35)
    
    cv_grpl_error = apply(cv_grpl, c(2,3), mean)
    minidx = which.min(cv_grpl_error)
    minphi_id = ifelse(minidx %% nphi !=0, minidx %% nphi, minidx %% nphi + nphi)
    minlam_id = ifelse(minidx %% nphi !=0, minidx %/% nphi + 1, minidx %/% nphi)

    #cv_grpl_error[minphi_id , minlam_id] # same with min(cv_grpl_error)
    
    ## Fitting training set ##
    grpl_fit_tr = grplFlinear(Y=sim_y[sim_trainind], X=penX_tr, Tps=Tps, lambda=lambda[minlam_id], phi=phi[minphi_id], dfs=35)

    ### validation using test set ###
    penX_tst = list()
    for (j in 1:10) {
      penX_tst[[j]] = penX[[j]][sim_testind,]
    }
    
    gflm_pred <- rep(grpl_fit_tr$intercept,length(sim_y[sim_testind]))
    for (jj in 1:nT_comp) {
      gflm_pred <- gflm_pred + (1/(nT-1))*penX_tst[[jj]] %*% grpl_fit_tr$Coef[[jj]]
    }
    for (jj in ((nS_comp+1) : (nS_comp+nT_comp))){
      gflm_pred <- gflm_pred + (1/(nS-1))*penX_tst[[jj]]%*%grpl_fit_tr$Coef[[jj]]
    }
    gflm_res[r,1:3] = c(sqrt(mean((gflm_pred-sim_fx[sim_testind])^2)),
                             mean(abs(gflm_pred-sim_fx[sim_testind])),
                             1 - (sum((sim_fx[sim_testind] - gflm_pred)^2)/sum((sim_fx[sim_testind] - mean(sim_fx[sim_testind]))^2)))

    print(paste("rep", r, "fin"))
  }

  print(Sys.time()-pre)
  res = list(fsvd=fsvd_res, time_avg_bsp=time_avg_bsp_res, time_avg_wv=time_avg_wv_res, 
             multi_pc_reg=multi_pc_reg_res, multi_pc_lasso = multi_pc_lasso_res, multi_pls=multi_pls_res, gflm=gflm_res)
  
  res 
}


##############################################################################
############################ Perform simulations ############################# 
##############################################################################


R = 100
N = 300 # number of two-way functional observations
nT = 101; nS = 256 # number of time points (0:100/100) / wavelength points (0:255/255)
nT_comp = nS_comp = 5 # number of components

#### generating coefficient functions ####
sim_beta_u = matrix(0, nrow=nT_comp, ncol=nT) 
sim_beta_v = matrix(0, nrow=nS_comp, ncol=nS)
set.seed(2022)
for (j in 1:nT_comp) {
  sim_beta_u[j,] = beta_u(nt=nT)
  sim_beta_v[j,] = beta_v(ns=nS)
}
par(mfrow=c(2,5))
for (j in 1:5) {
  plot(0:(nT-1)/(nT-1), sim_beta_u[j,], type="l", xlab="time", ylab="")
}
for (j in 1:5) {
  plot(0:(nS-1)/(nS-1), sim_beta_v[j,], type="l", xlab="Wavelength", ylab="")
}
#par(mfrow=c(1,1))
#for (j in 1:nT_comp) {
#  sim_beta_u[j,] = beta_u2(nt=nT, l=j)
#  sim_beta_v[j,] = beta_v2(ns=nS, l=j)
#}
par(mfrow=c(2,2))


d = (nT_comp:1)^2
#par(mar = c(3.9, 2.8, 2.8, 3))

# test
#sim_res_set11 = get_simulation(R=1, N=100, beta_u=sim_beta_u, beta_v=sim_beta_v, 
#                               idx_ucomp=1:5, idx_vcomp=1:5, noise_u=0.05, noise_v=0.05, d=d, SNR=20)

###################################################################################

### noise_x = 0.05, SNR=20 ###
sim_res_set1 = get_simulation(R=100, N=300, beta_u=sim_beta_u, beta_v=sim_beta_v, 
                         idx_ucomp=1:5, idx_vcomp=1:5, noise_u=0.05, noise_v=0.05, d=d, SNR=20)
set1_mae = cbind(sim_res_set1$fsvd[,c(2,5,8,11,14,17)], sim_res_set1$time_avg_bsp[,2], sim_res_set1$time_avg_wv[,2],
                 sim_res_set1$multi_pc_reg[,2], sim_res_set1$multi_pc_lasso[,2], 
                 sim_res_set1$multi_pls[,2], sim_res_set1$gflm[,2])
#colnames(set1_mae) = c("FSVD1", "FSVD2", "FSVD3", "FSVD4", "Plain1", "Plain2", "Avg1", "Avg2",
#                           "PC1", "PC2", "PLS", "GFLM")
apply(set1_mae[1:50,], 2, mean) %>% round(4)
apply(set1_mae[1:50,], 2, sd) %>% round(4)
par(mfrow=c(2,2))
boxplot(set1_mae[1:50,], cex.axis=0.8, ylab="MAE", 
        main=expression(sigma[x]^2 == paste(0.05, ",   SNR=20")),
        col=c("blue", "blue", "blue", "blue","cyan", "cyan", "lightgrey","lightgrey","lightgrey","lightgrey","lightgrey","green"))


###################################################################################

par(mfrow=c(1,1))
sim_res_set2 = get_simulation(R=50, N=300, beta_u=sim_beta_u, beta_v=sim_beta_v, 
                         idx_ucomp=1:5, idx_vcomp=1:5, noise_u=0.05, noise_v=0.05, d=d, SNR=2)
set2_mae = cbind(sim_res_set2$fsvd[,c(2,5,8,11,14,17)], sim_res_set2$time_avg_bsp[,2], sim_res_set2$time_avg_wv[,2],
                 sim_res_set2$multi_pc_reg[,2], sim_res_set2$multi_pc_lasso[,2], 
                 sim_res_set2$multi_pls[,2], sim_res_set2$gflm[,2])
colnames(set2_mae) = c("FSVD1", "FSVD2", "FSVD3", "FSVD4", "Plain1", "Plain2", "Avg1", "Avg2",
                       "PC1", "PC2", "PLS", "GFLM")
apply(set2_mae, 2, mean) %>% round(4)
apply(set2_mae, 2, sd) %>% round(4)
boxplot(set2_mae, cex.axis=0.8, ylab="MAE",
        main=expression(sigma[x]^2 == paste(0.05, ",   SNR=2")),
        col=c("blue", "blue", "blue", "blue","cyan", "cyan", "lightgrey","lightgrey","lightgrey","lightgrey","lightgrey","green"))


###################################################################################

sim_res_set3 = get_simulation(R=50, N=300, beta_u=sim_beta_u, beta_v=sim_beta_v, 
                              idx_ucomp=1:5, idx_vcomp=1:5, noise_u=0.5, noise_v=0.5, d=d, SNR=20)

set3_mae = cbind(sim_res_set3$fsvd[,c(2,5,8,11,14,17)], sim_res_set3$time_avg_bsp[,2], sim_res_set3$time_avg_wv[,2],
                 sim_res_set3$multi_pc_reg[,2], sim_res_set3$multi_pc_lasso[,2], 
                 sim_res_set3$multi_pls[,2], sim_res_set3$gflm[,2])
colnames(set3_mae) = c("FSVD1", "FSVD2", "FSVD3", "FSVD4", "Plain1", "Plain2", "Avg1", "Avg2",
                       "PC1", "PC2", "PLS", "GFLM")
apply(set3_mae, 2, mean) %>% round(4)
apply(set3_mae, 2, sd) %>% round(4)
boxplot(set3_mae, cex.axis=0.8, ylab="MAE",
        main=expression(sigma[x]^2 == paste(0.5, ",   SNR=20")),
        col=c("blue", "blue", "blue", "blue","cyan", "cyan", "lightgrey","lightgrey","lightgrey","lightgrey","lightgrey","green"))


###################################################################################

sim_res_set4 = get_simulation(R=50, N=300, beta_u=sim_beta_u, beta_v=sim_beta_v, 
                              idx_ucomp=1:5, idx_vcomp=1:5, noise_u=0.5, noise_v=0.5, d=d, SNR=2)
set4_mae = cbind(sim_res_set4$fsvd[,c(2,5,8,11,14,17)], sim_res_set4$time_avg_bsp[,2], sim_res_set4$time_avg_wv[,2],
                 sim_res_set4$multi_pc_reg[,2], sim_res_set4$multi_pc_lasso[,2], 
                 sim_res_set4$multi_pls[,2], sim_res_set4$gflm[,2])
colnames(set4_mae) = c("FSVD1", "FSVD2", "FSVD3", "FSVD4", "Plain1", "Plain2", "Avg1", "Avg2",
                       "PC1", "PC2", "PLS", "GFLM")
apply(set4_mae, 2, mean) %>% round(4)
apply(set4_mae, 2, sd) %>% round(4)

str(set4_mae)
colnames(set4_mae) <- NULL
boxplot(set4_mae, cex.axis=0.8, ylab="MAE",
        main=expression(sigma[x]^2 == paste(0.5, ",   SNR=2")),
        col=c("blue", "blue", "blue", "blue","cyan", "cyan", "lightgrey","lightgrey","lightgrey","lightgrey","lightgrey","green"))



###################################################################################



