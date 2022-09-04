######### Step 1 - PenFSVD #########
source("1_PenFSVD.r")

Sys.time()
tw_res_SS1_3 <- mclapply(
  1:3,
  function(x){
    apply(data_log[[x]], 1, tw_svd_smu, ncomp = 5)
  },
  mc.cores = 20
)

Sys.time()
tw_res_SS4_6 <- mclapply(
  1:3,
  function(x){
    x <- x+3
    apply(data_log[[x]], 1, tw_svd_smu, ncomp = 5)
  },
  mc.cores = 20
)

Sys.time()
tw_res_SS7_9 <- mclapply(
  1:3,
  function(x){
    x <- x+6
    apply(data_log[[x]], 1, tw_svd_smu, ncomp = 5)
  },
  mc.cores = 20
)
Sys.time()

# Get full joined results
tw_res <- list(
  tw_res_SS1_3[[1]], tw_res_SS1_3[[2]], tw_res_SS1_3[[3]],
  tw_res_SS4_6[[1]], tw_res_SS4_6[[2]], tw_res_SS4_6[[3]],
  tw_res_SS7_9[[1]], tw_res_SS7_9[[2]], tw_res_SS7_9[[3]]
)


################################# Flip 
##### Sign of U, V might vary by observation,
#####   so manually find and 'flip' certain observations
##### Result: filp_res, list of length 9,
#####   each: list of length k(=3~5, see vector use_pair),
#####     each: list of length 2 (time component u, wv component v)

# use_pair <- list(1:5, 1:3, 1:5, 1:5, 1:3, 1:5, 1:5, 1:3, 1:5)
# 
# get_pair <- function(ss, pind){
#   list(
#     u = t(sapply(1:307, function(x){tw_res[[ss]][[x]]$pen_u[,pind]}, simplify = "array")),
#     v = t(sapply(1:307, function(x){tw_res[[ss]][[x]]$pen_v[,pind]}, simplify = "array"))
#   )
# }
# 
# pair_list <- vector("list", length = 9)
# for(i in 1:9){
#   for(j in 1:5){
#     pair_list[[i]][[j]] <- get_pair(i, j)
#   }
# }
# 
# what_to_flip <- list(NA, c(3,4), c(5), c(4,5), 
#                      3, c(3,5), c(3,4,5), c(3), c(4))
# point <- list(NA, c(30,20), c(20), c(50,50), c(15), 
#               c(20, 190), c(10, 30, 100), c(15), c(300))
# 
# flipper <- function(ss, pind, point, is_sign_plus = T){
#   flip_u = pair_list[[ss]][[pind]]$u
#   flip_v = pair_list[[ss]][[pind]]$v
#   
#   if(is_sign_plus){
#     flipind <- which(flip_u[, pind] <0)
#   } else {
#     flipind <- which(flip_v[, pind] >0)
#   }
#   flip_u[flipind, ] <- -flip_u[flipind, ]
#   flip_v[flipind, ] <- -flip_v[flipind, ]
#   list(u = flip_u, v = flip_v)
# }
# 
# filp_res <- pair_list
# for(i in 2:9){
#   for(j in what_to_flip[[i]]){
#     filp_res[[i]][[j]] <- flipper(i, j, point[[i]][[j]])
#   }
# }
# 
# filp_res[[7]][[5]] <- flipper(7, 5, 80)
# filp_res[[6]][[5]] <- flipper(6, 5, 170)
# filp_res[[2]][[3]] <- flipper(2, 3, 25)
# 
# filp_res[[2]][[3]]$u[c(8, 229), ] <- -filp_res[[2]][[3]]$u[c(8, 229), ]
# filp_res[[2]][[3]]$v[c(8, 229), ] <- -filp_res[[2]][[3]]$v[c(8, 229), ]
# matplot(t(filp_res[[2]][[3]]$u), type = "l")
# 
# tt <- which(filp_res[[7]][[5]]$u[, 75] >0)
# filp_res[[7]][[5]]$u[tt, ] <- -filp_res[[7]][[5]]$u[tt, ]
# filp_res[[7]][[5]]$v[tt, ] <- -filp_res[[7]][[5]]$v[tt, ]
# matplot(t(filp_res[[7]][[5]]$u), type = "l")
# 
# tt <- which(filp_res[[6]][[5]]$u[, 170] >0)
# filp_res[[6]][[5]]$u[tt, ] <- -filp_res[[6]][[5]]$u[tt, ]
# filp_res[[6]][[5]]$v[tt, ] <- -filp_res[[6]][[5]]$v[tt, ]
# matplot(t(filp_res[[6]][[5]]$u), type = "l")


# matplot(t(filp_res[[1]][[1]]$u[,1:20]), type = "l")
# matplot(t(filp_res[[1]][[2]]$u), type = "l")
# matplot(t(filp_res[[1]][[3]]$u[,100:200]), type = "l")
# matplot(t(filp_res[[1]][[4]]$u), type = "l")
# matplot(t(filp_res[[1]][[5]]$u), type = "l")

## 결국엔 수정하긴 해야함
oes_penfsvd <- list(u = vector("list", length = 45),
                    v = vector("list", length = 45))
for(j in 1:45){
  nc = ifelse(j %% 5 == 0, 5, j %% 5)
  ss = ((j - nc) %/% 5) + 1
  oes_penfsvd$u[[j]] <- filp_res[[ss]][[nc]]$u
  oes_penfsvd$v[[j]] <- filp_res[[ss]][[nc]]$v
}

par(mfrow = c(5, 9))
for(j in 1:45){
  matplot(t(oes_penfsvd$u[[j]]), type = "l", xlab = "")
}
# 9, 10, 24, 25, 39, 40
# 
# 
# par(mfrow = c(2,3))
# for(j in c(9, 10, 24, 25, 39, 40)){
#   matplot(t(oes_penfsvd$u[[j]]), type = "l", xlab = "")
# }
# 
# par(mfrow = c(1,1))
# ## 9
# matplot(t(oes_penfsvd$u[[9]]), type = "l")
# matplot( t(oes_penfsvd$v[[9]][oes_penfsvd$v[[9]][, 1] < 0,]), type = "l")
# matplot(-t(oes_penfsvd$v[[9]][oes_penfsvd$v[[9]][, 1] >= 0,]), type = "l",
#         add = F)
# tmp_ind <- which(oes_penfsvd$u[[9]][, 150] < 0)
# for(i in tmp_ind){
#   oes_penfsvd$u[[9]][i,] <- -oes_penfsvd$u[[9]][i,]
#   oes_penfsvd$v[[9]][i,] <- -oes_penfsvd$v[[9]][i,]
# }
# matplot(t(oes_penfsvd$u[[9]]), type = "l")
# matplot(t(oes_penfsvd$v[[9]][1:2,]), type = "l")
# 
# 
# ## 10
# matplot(t(oes_penfsvd$u[[10]]), type = "l")
# matplot(t(oes_penfsvd$v[[10]]), type = "l")
# matplot( t(oes_penfsvd$v[[10]][oes_penfsvd$v[[10]][, 1000] < 0,]), type = "l")
# matplot(-t(oes_penfsvd$v[[10]][oes_penfsvd$v[[10]][, 1000] >= 0,]), type = "l",
#         add = T)
# tmp_ind <- which(oes_penfsvd$v[[10]][, 1000] < 0)
# for(i in tmp_ind){
#   oes_penfsvd$u[[10]][i,] <- -oes_penfsvd$u[[10]][i,]
#   oes_penfsvd$v[[10]][i,] <- -oes_penfsvd$v[[10]][i,]
# }
# matplot(t(oes_penfsvd$u[[10]]), type = "l")
# 
# #24
# matplot(t(oes_penfsvd$u[[24]]), type = "l")
# matplot(t(oes_penfsvd$v[[24]]), type = "l")
# matplot( t(oes_penfsvd$u[[24]][oes_penfsvd$u[[24]][, 1] < 0,]), type = "l")
# matplot(-t(oes_penfsvd$u[[24]][oes_penfsvd$u[[24]][, 1] >= 0,]), type = "l",
#         add = T)
# tmp_ind <- which(oes_penfsvd$u[[24]][, 1] < 0)
# for(i in tmp_ind){
#   oes_penfsvd$u[[24]][i,] <- -oes_penfsvd$u[[24]][i,]
#   oes_penfsvd$v[[24]][i,] <- -oes_penfsvd$v[[24]][i,]
# }
# matplot(t(oes_penfsvd$u[[24]]), type = "l")
# matplot(t(oes_penfsvd$v[[24]]), type = "l")
# 
# 
# #25
# matplot(t(oes_penfsvd$u[[25]]), type = "l")
# tmp_ind <- which(oes_penfsvd$u[[25]][, 1] < 0)
# for(i in tmp_ind){
#   oes_penfsvd$u[[25]][i,] <- -oes_penfsvd$u[[25]][i,]
#   oes_penfsvd$v[[25]][i,] <- -oes_penfsvd$v[[25]][i,]
# }
# matplot(t(oes_penfsvd$u[[25]]), type = "l")
# matplot(t(oes_penfsvd$v[[25]]), type = "l")
# 
# #39
# matplot(t(oes_penfsvd$u[[39]]), type = "l")
# tmp_ind <- which(oes_penfsvd$u[[39]][, 1] < 0)
# for(i in tmp_ind){
#   oes_penfsvd$u[[39]][i,] <- -oes_penfsvd$u[[39]][i,]
#   oes_penfsvd$v[[39]][i,] <- -oes_penfsvd$v[[39]][i,]
# }
# matplot(t(oes_penfsvd$u[[39]]), type = "l")
# matplot(t(oes_penfsvd$v[[39]]), type = "l")


#40: ?

#### Step 1 - BasFSVD ####
source("1_BasFSVD.r")

detectCores()
Sys.time() #"2022-08-23 23:08:05 KST"
basis_raw <- mclapply(
  1:9,
  function(x){
    apply(data_log[[x]], 1, basis_svd, ncmp = 5)
  },
  mc.cores = 20
)
Sys.time() #"2022-08-23 23:33:33 KST"


oes_basfsvd <- list(u = vector("list", length = 45),
                    v = vector("list", length = 45))
for(j in 1:45){
  nc = ifelse(j %% 5 == 0, 5, j %% 5)
  ss = ((j - nc) %/% 5) + 1
  for(i in 1:307){
    oes_basfsvd$u[[j]] <- rbind(oes_basfsvd$u[[j]], basis_raw[[ss]][[i]]$u[,nc])
    oes_basfsvd$v[[j]] <- rbind(oes_basfsvd$v[[j]], basis_raw[[ss]][[i]]$v[,nc])
  }
}

par(mfrow = c(3, 5))
for(j in 1:45){
  matplot(t(oes_basfsvd$u[[j]]), type = "l", xlab = "", main = j)
}


plot_flip <- function(data = oes_basfsvd$u, ind){
  par(mfrow = c(1,1))
  matplot(t(data[[ind]]), type = "l")
}

flip_at <- function(data = oes_basfsvd, ind, at){
  tmp_ind <- which(data$u[[ind]][, at] < 0)
  for(i in tmp_ind){
    data$u[[ind]][i,] <- -data$u[[ind]][i,]
    data$v[[ind]][i,] <- -data$v[[ind]][i,]
  }
  matplot(t(data$u[[ind]]), type = "l")
  return(list(u = data$u[[ind]], v = data$v[[ind]]))
}

#8 
plot_flip(ind = 8)
tmp <- flip_at(ind = 8, at = 25)
oes_basfsvd$u[[8]] <- tmp$u
oes_basfsvd$v[[8]] <- tmp$v
matplot(t(oes_basfsvd$u[[8]]), type = "l")
matplot(t(oes_basfsvd$v[[8]]), type = "l")

#9 
plot_flip(ind = 9)
tmp <- flip_at(ind = 9, at = 25)
oes_basfsvd$u[[9]] <- tmp$u
oes_basfsvd$v[[9]] <- tmp$v
matplot(t(oes_basfsvd$u[[9]]), type = "l")
matplot(t(oes_basfsvd$v[[9]]), type = "l")

#10
plot_flip(ind = 10)
tmp <- flip_at(ind = 10, at = 70)
oes_basfsvd$u[[10]] <- tmp$u
oes_basfsvd$v[[10]] <- tmp$v
matplot(t(oes_basfsvd$v[[10]]), type = "l")

#13 
plot_flip(ind = 13)
tmp <- flip_at(ind = 13, at = 20)
oes_basfsvd$u[[13]] <- tmp$u
oes_basfsvd$v[[13]] <- tmp$v
matplot(t(oes_basfsvd$v[[13]]), type = "l")

#15 
plot_flip(ind = 15)
tmp <- flip_at(ind = 15, at = 20)
oes_basfsvd$u[[15]] <- tmp$u
oes_basfsvd$v[[15]] <- tmp$v
matplot(t(oes_basfsvd$v[[15]]), type = "l")

#20 
plot_flip(ind = 20)
tmp <- flip_at(ind = 20, at = 150)
oes_basfsvd$u[[20]] <- tmp$u
oes_basfsvd$v[[20]] <- tmp$v
matplot(t(oes_basfsvd$v[[20]]), type = "l")

#30 
plot_flip(ind = 30)
tmp <- flip_at(ind = 30, at = 195)
oes_basfsvd$u[[30]] <- tmp$u
oes_basfsvd$v[[30]] <- tmp$v
matplot(t(oes_basfsvd$v[[30]]), type = "l")

#34 
plot_flip(ind = 34)
tmp <- flip_at(ind = 34, at = 170)
oes_basfsvd$u[[34]] <- tmp$u
oes_basfsvd$v[[34]] <- tmp$v
matplot(t(oes_basfsvd$v[[34]]), type = "l")

#35 
plot_flip(ind = 35)
tmp <- flip_at(ind = 35, at = 170)
oes_basfsvd$u[[35]] <- tmp$u
oes_basfsvd$v[[35]] <- tmp$v
matplot(t(oes_basfsvd$v[[35]]), type = "l")

#44 
plot_flip(ind = 44)
tmp <- flip_at(ind = 44, at = 50)
oes_basfsvd$u[[44]] <- tmp$u
oes_basfsvd$v[[44]] <- tmp$v
matplot(t(oes_basfsvd$v[[44]]), type = "l")


par(mfrow = c(3, 5))
for(j in 1:45){
  matplot(t(oes_basfsvd$u[[j]]), type = "l", xlab = "", main = j)
}


#### Step 2 - Bspline ####
##########################################from here..
source("2_jma.R")
oes_res <- list()
#### Pen + Bsp
Sys.time() #2022-08-13 18:31:36 KST
oes_res$PenBsp <- jk_bsp(
  X = c(oes_penfsvd$u, oes_penfsvd$v), y_ = y$Y,
  x_nb = c(rep(200, 5), rep(200, 5), rep(100, 5),
           rep(168, 5), rep(58, 5), rep(199, 5),
           rep(173, 5), rep(51, 5), rep(200, 5),
           rep(100, 45)),
  b_nb = c(rep(50, 5), rep(80, 5), rep(20, 5),
           rep(33, 5), rep(11, 5), rep(40, 5),
           rep(34, 5), rep(10, 5), rep(63, 5),
           rep(100, 45)),
  tr = train_ind, tst = test_ind
)
Sys.time() # 2022-08-13 19:50:17 KST

#### Basis + Bsp
Sys.time() #"2022-08-24 01:21:41 KST"
oes_res$BasBsp <- jk_bsp(
  X = c(oes_basfsvd$u, oes_basfsvd$v), y_ = y$Y,
  x_nb = c(rep(200, 5), rep(200, 5), rep(100, 5),
           rep(168, 5), rep(58, 5), rep(199, 5),
           rep(173, 5), rep(51, 5), rep(200, 5),
           rep(100, 45)),
  b_nb = c(rep(50, 5), rep(80, 5), rep(20, 5),
           rep(33, 5), rep(11, 5), rep(40, 5),
           rep(34, 5), rep(10, 5), rep(63, 5),
           rep(100, 45)),
  tr = train_ind, tst = test_ind
)
Sys.time() #"2022-08-24 02:40:10 KST"

#### Step 2 - Bsp(Time) + Wv(Wvlngth) ####
#### Pen + Wv
oes_res$PenWv150 <- jk_wv(
  X = oes_penfsvd$v, y_ = y$Y, q = 0.15,
  tr = train_ind, tst = test_ind
)
oes_res$PenMerged150 <- jk_merge(
  error_mat1 = oes_res$PenBsp$error_mat[, 1:45],
  error_mat2 = oes_res$PenWv150$error_mat,
  y_m1 = oes_res$PenBsp$y_m[, 1:45],
  y_m2 = oes_res$PenWv150$y_m,
  y_ = y$Y, tr = train_ind, tst = test_ind )

#### Bas + Wv
oes_res$BasWv150 <- jk_wv(
  X = oes_basfsvd$v, y_ = y$Y, q = 0.15,
  tr = train_ind, tst = test_ind
)
oes_res$BasMerged150 <- jk_merge(
  error_mat1 = oes_res$BasBsp$error_mat[, 1:45], 
  error_mat2 = oes_res$BasWv150$error_mat,
  y_m1 = oes_res$BasBsp$y_m[, 1:45], 
  y_m2 = oes_res$BasWv150$y_m,
  y_ = y$Y, tr = train_ind, tst = test_ind )
Sys.time() # "2022-08-24 02:47:08 KST"

##................................................................
#### Comparison ####
#### Wavelength selection: Find Peaks

library(ggpmisc)
matplot(t(data_log[[1]][, 100,]), type = "l")
sum(ggpmisc:::find_peaks(data_log[[1]][1, 100,],
                         span = 7))


ggplot(data = data.frame(x=1:1024, y=data_log[[1]][1, 100,]), 
       aes(x = x, y = y)) + 
  geom_line() + 
  stat_peaks(col = "red", span = 7)


##### 1. Compression with time, selected wavelengths #####
# use only peak values with time too??IIRR as ref.
peak <- vector("list", length = 9)
get_peak_wv <- function(data, span = 5, npeaks = 80){
  nobs = dim(data)[1]
  nt = dim(data)[2]
  nwv = dim(data)[3] #1024
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
for(i in 1:9){
  peak[[i]] <- get_peak_wv(data = data_log[[i]])
}

anyNA(peak)
peak[[5]] <- peak[[5]][, apply(peak[[5]], 2, function(x) !any(is.na(x)))]
peak[[6]] <- peak[[6]][, apply(peak[[6]], 2, function(x) !any(is.na(x)))]
peak[[7]] <- peak[[7]][, apply(peak[[7]], 2, function(x) !any(is.na(x)))]
peak[[8]] <- peak[[8]][, apply(peak[[8]], 2, function(x) !any(is.na(x)))]
peak[[9]] <- peak[[9]][, apply(peak[[9]], 2, function(x) !any(is.na(x)))]


peak_pca_train <- vector("list", length = 9)
for(i in c(1:9)){
  peak_pca_train[[i]] <- princomp(peak[[i]][train_ind,])
}

### Number of PC = 5 --> variance explained > 0.98 for all processes
# summary(peak_pca_train[[1]]) #0.997138952
# summary(peak_pca_train[[2]]) #0.993343703
# summary(peak_pca_train[[3]]) #0.997673906
# summary(peak_pca_train[[4]]) #0.984180312
# summary(peak_pca_train[[5]]) #0.983259075
# summary(peak_pca_train[[6]]) #0.989314009
# summary(peak_pca_train[[7]]) #0.99092154
# summary(peak_pca_train[[8]]) #0.990003690
# summary(peak_pca_train[[9]]) #0.984430596


peak_compdata <- cbind(
  peak_pca_train[[1]]$scores[, 1:5],
  peak_pca_train[[2]]$scores[, 1:5],
  peak_pca_train[[3]]$scores[, 1:5],
  peak_pca_train[[4]]$scores[, 1:5],
  peak_pca_train[[5]]$scores[, 1:5],
  peak_pca_train[[6]]$scores[, 1:5],
  peak_pca_train[[7]]$scores[, 1:5],
  peak_pca_train[[8]]$scores[, 1:5],
  peak_pca_train[[9]]$scores[, 1:5]
)
colnames(peak_compdata) <- c(
  paste(rep("SS1PC", 5), 1:5, sep = ""),
  paste(rep("SS2PC", 5), 1:5, sep = ""),
  paste(rep("SS3PC", 5), 1:5, sep = ""),
  paste(rep("SS4PC", 5), 1:5, sep = ""),
  paste(rep("SS5PC", 5), 1:5, sep = ""),
  paste(rep("SS6PC", 5), 1:5, sep = ""),
  paste(rep("SS7PC", 5), 1:5, sep = ""),
  paste(rep("SS8PC", 5), 1:5, sep = ""),
  paste(rep("SS9PC", 5), 1:5, sep = "")
)
peak_compdata #dim: 250 * 45


### a-PCR
lm_fit <- lm(y~., data = data.frame(y = y$Y[train_ind], peak_compdata))

peak_compdata_test <- cbind(
  predict(peak_pca_train[[1]], newdata = peak[[1]][test_ind,])[, 1:5],
  predict(peak_pca_train[[2]], newdata = peak[[2]][test_ind,])[, 1:5],
  predict(peak_pca_train[[3]], newdata = peak[[3]][test_ind,])[, 1:5],
  predict(peak_pca_train[[4]], newdata = peak[[4]][test_ind,])[, 1:5],
  predict(peak_pca_train[[5]], newdata = peak[[5]][test_ind,])[, 1:5],
  predict(peak_pca_train[[6]], newdata = peak[[6]][test_ind,])[, 1:5],
  predict(peak_pca_train[[7]], newdata = peak[[7]][test_ind,])[, 1:5],
  predict(peak_pca_train[[8]], newdata = peak[[8]][test_ind,])[, 1:5],
  predict(peak_pca_train[[9]], newdata = peak[[9]][test_ind,])[, 1:5]
)
colnames(peak_compdata_test) <- c(
  paste(rep("SS1PC", 5), 1:5, sep = ""),
  paste(rep("SS2PC", 5), 1:5, sep = ""),
  paste(rep("SS3PC", 5), 1:5, sep = ""),
  paste(rep("SS4PC", 5), 1:5, sep = ""),
  paste(rep("SS5PC", 5), 1:5, sep = ""),
  paste(rep("SS6PC", 5), 1:5, sep = ""),
  paste(rep("SS7PC", 5), 1:5, sep = ""),
  paste(rep("SS8PC", 5), 1:5, sep = ""),
  paste(rep("SS9PC", 5), 1:5, sep = "")
)

tbl1 <- data.frame(rmse = rep(0,9), mae = rep(0,9), r2 = rep(0,9))
rownames(tbl1) <- c("a-PCR", "a-PCLasso", "a-PLS",
                    "b-Bsp", "b-wv",
                    "PenBsp", "BasBsp", "PenWv", "BasWv")
oes_res$a_pred <- matrix(nrow = 57, ncol = 9)
colnames(oes_res$a_pred) <- c("a-PCR", "a-PCLasso", "a-PLS",
                              "b-Bsp", "b-wv",
                              "PenBsp", "BasBsp", "PenWv", "BasWv")

oes_res$a_pred[,1] <- predict(lm_fit, data.frame(peak_compdata_test))
tbl1[1,] <- c(
  sqrt(mean((pred[,1] - y$Y[test_ind])^2)),
  mean(abs(pred[,1] - y$Y[test_ind])),
  1 - (sum((y$Y[test_ind] - pred[,1])^2)/sum((y$Y[test_ind] - mean(y$Y[test_ind]))^2))
)






tbl_jma <- matrix(NA, nrow = 9, ncol = 3)
rownames(tbl_jma) <- c(
  "a-PCR", "a-PCLasso", "a-PLS",
  "b-Bsp", "b-Wv",
  "PenBsp", "BasBsp", "PenWv", "BasWv"
)
colnames(tbl_jma) <- c("RMSE", "MAE", "R2")
tbl_jma[6,] <- c(oes_res$PenBsp$rmse, oes_res$PenBsp$mae, oes_res$PenBsp$r2)
tbl_jma[7,] <- c(oes_res$BasBsp$rmse, oes_res$BasBsp$mae, oes_res$BasBsp$r2)
tbl_jma[8,] <- c(oes_res$PenMerged150$rmse, oes_res$PenMerged150$mae, oes_res$PenMerged150$r2)
tbl_jma[9,] <- c(oes_res$BasMerged150$rmse, oes_res$BasMerged150$mae, oes_res$BasMerged150$r2)


tbl_jma[1,] <- c(tbl1[1,1], tbl1[1,2], tbl1[1,3])
tbl_jma[2,] <- c(tbl1[2,1], tbl1[2,2], tbl1[2,3])
tbl_jma[3,] <- c(tbl1[3,1], tbl1[3,2], tbl1[3,3])


Sys.time() # "2022-08-18 03:37:41 KST"
oes_res$b_bsp <- jk_bsp(
  X = time_red_data, y_ = y$Y,
  x_nb = rep(100, 45),
  b_nb = rep(100, 45),
  tr = train_ind, tst = test_ind
)
Sys.time() #"2022-08-18 03:46:26 KST"
oes_res$b_wv <- jk_wv(
  X = time_red_data, y_ = y$Y, q = 0.15,
  tr = train_ind, tst = test_ind
)
Sys.time() #"2022-08-18 03:47:50 KST"
tbl_jma[4,] <- c(oes_res$b_bsp$rmse, oes_res$b_bsp$mae, oes_res$b_bsp$r2)
tbl_jma[5,] <- c(oes_res$b_wv$rmse, oes_res$b_wv$mae, oes_res$b_wv$r2)
round(tbl_jma, 4)




#######################
oes_res$b_bsp$weight$solution %>%  round(., 4)
oes_res$b_wv$weight$solution %>%   round(., 4)
oes_res$PenBsp$weight$solution %>% round(., 4)







# par(mfrow = c(1,1))
# matplot(t(oes_penfsvd$u[[41]])[, 1:5], type = "l", col = "black")
# matplot(t(oes_penfsvd$u[[42]])[, 1:5], type = "l", col = "black")
# matplot(t(oes_penfsvd$u[[43]])[, 1:5], type = "l", col = "black")
# matplot(t(oes_penfsvd$u[[44]])[, 1:5], type = "l", col = "black")
# matplot(t(oes_penfsvd$u[[45]])[, 1:5], type = "l", col = "black")
