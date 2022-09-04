# Decomposition - U
tmp <- array(dim = c(307, 58, 2))
for(i in 1:307){
  tmp[i,,] <- svd(data_log[[5]][i,,])$u[,1:2]
}
set.seed(10101)
i <- sample(1:307, 10)
pdf(file = "plot4paper/decomp_u.pdf", width = 9, height = 5)
par(mfrow = c(2, 3), 
    oma = c(2, 2, 0, 0), # two rows of text at the outer left and bottom margin
    mar = c(1.5, 2, 2, 1), # space for one row of text at ticks and to separate plots
    mgp = c(2.2, 1, 0),    # axis label at 2 rows distance, tick labels at 1 row
    xpd = NA)            # allow content to protrude into outer margin (and beyond))
matplot(t(tmp[i,,1]), type = "l", col = "black",
        xlab = "", ylab = "U1", main = "FSVD")
matplot(t(oes_penfsvd$u[[21]][i,]), type = "l", col = "black", 
        xlab = "", ylab = "", main = "PenFSVD")
matplot(t(oes_basfsvd$u[[21]][i,]), type = "l", col = "black", 
        xlab = "", ylab = "", main = "BasisFSVD")
matplot(t(tmp[i,,2]), type = "l", col = "black",
        xlab = "time", ylab = "U2", main = "")
matplot(t(oes_penfsvd$u[[22]][i,]), type = "l", col = "black", 
        xlab = "time", ylab = "", main = "")
matplot(t(oes_basfsvd$u[[22]][i,]), type = "l", col = "black", 
        xlab = "time", ylab = "", main = "")
dev.off()


# Decomposition - V
pdf(file = "plot4paper/decomp_v.pdf", 
    width = 9, height = 5)
par(mfrow = c(2, 3), 
    oma = c(2, 2, 0, 0), # two rows of text at the outer left and bottom margin
    mar = c(1.5, 2, 2, 1), # space for one row of text at ticks and to separate plots
    mgp = c(2, 1, 0),    # axis label at 2 rows distance, tick labels at 1 row
    xpd = NA)
matplot(t(time_red_data[[9]][i,]), 
        type ="l", col = "black", lty = "dashed", 
        xlab = "wavelength", ylab = "average over time", main = "(b) Time-average")
matplot(-t(oes_penfsvd$v[[41]][i,]), 
        type = "l", col = "black", lty = "dashed", 
        xlab = "", ylab = "V1", main = "PenFSVD")
matplot(-t(oes_basfsvd$v[[41]][i,]), 
        type = "l", col = "black", lty = "dashed", 
        xlab = "", ylab = "", main = "BasisFSVD")
plot.new()
matplot(-t(oes_penfsvd$v[[42]][1:5,]), 
        type = "l", col = "black", lty = "dashed", 
        xlab = "wavelength", ylab = "V2", main = "")
matplot(-t(oes_basfsvd$v[[42]][1:5,]), 
        type = "l", col = "black", lty = "dashed", 
        xlab = "wavelength", ylab ="", main = "")
dev.off()


# Result table
round(tbl_jma, 4)



# RMSE boxplot



# Weight analysis
w <- data.frame(
  PenBsp = oes_res$PenBsp$weight$solution,
  BasBsp = oes_res$BasBsp$weight$solution,
  PenWav = oes_res$PenMerged150$weight$solution,
  BasWav = oes_res$BasMerged150$weight$solution
)
w$key <- c(paste(rep(1:9, each = 5), rep(1:5, 9), sep = "_"),
           paste(rep(1:9, each = 5), rep(1:5, 9), sep = "_"))
w$comp <- c(rep("U", 45), rep("V", 45))

library(reshape2)



pdf(file = "plot4paper/weight_heatmap.pdf", width = 5.5, height = 6)
melt(data = w, id.vars = c("key", "comp"), value.name = "weight") %>%
  ggplot(data = ., aes(x = comp, y = key, fill = weight)) +
  geom_tile(color = "grey") +
  scale_fill_gradient(low = "white", high = "red") +
  coord_equal(ratio = 0.25) +
  facet_wrap(~variable, ncol = 4) +
  labs(title="", x ="components", y = "coefficients")
dev.off()


apply(w[, 1:4], 2, function(x){sum(abs(x)>0.00000001)})
apply(w[1:45, 1:4], 2, function(x){sum(abs(x)>0.00000001)})
apply(w[46:90, 1:4], 2, function(x){sum(abs(x)>0.00000001)})


# Coef plot



# Scatterplot
pdf(file = "plot4paper/jma_pred.pdf",
    width = 9, height = 27/4)
par(mfrow = c(3,4),
    oma = c(2, 2, 0, 0), # two rows of text at the outer left and bottom margin
    mar = c(1.5, 2, 2, 1), # space for one row of text at ticks and to separate plots
    mgp = c(2, 1, 0),    # axis label at 2 rows distance, tick labels at 1 row
    xpd = F)
plot(x = oes_res$a_pred[,1], y = y$Y[test_ind],
     xlim = c(0,1), ylim = c(0,1), xlab = "", ylab = "y",
     main = "(a) PCR")
abline(a=0, b=1, col="red")
plot(x = oes_res$a_pred[,2], y = y$Y[test_ind],
     xlim = c(0,1), ylim = c(0,1), xlab = "", ylab = "",
     main = "(a) PC + Lasso")
abline(a=0, b=1, col="red")
plot(x = oes_res$a_pred[,3], y = y$Y[test_ind],
     xlim = c(0,1), ylim = c(0,1), xlab = "", ylab = "",
     main = "(a) PLS")
abline(a=0, b=1, col="red")
plot.new()

plot(x = oes_res$b_bsp$y_jpa, y = y$Y[test_ind],
     xlim = c(0,1), ylim = c(0,1), xlab = "", ylab = "y",
     main = "(b) B-spline reg")
abline(a=0, b=1, col="red")
plot(x = oes_res$b_wv$y_jpa, y = y$Y[test_ind],
     xlim = c(0,1), ylim = c(0,1), xlab = "", ylab = "",
     main = "(b) Wavelet reg")
abline(a=0, b=1, col="red")
plot.new()
plot.new()

plot(x = oes_res$PenBsp$y_jpa, y = y$Y[test_ind],
     xlim = c(0,1), ylim = c(0,1), xlab = "y_hat", ylab = "y",
     main = "PenFSVD + B-spline reg")
abline(a=0, b=1, col="red")
plot(x = oes_res$BasBsp$y_jpa, y = y$Y[test_ind],
     xlim = c(0,1), ylim = c(0,1), xlab = "y_hat", ylab = "",
     main = "BasisFSVD + B-spline reg")
abline(a=0, b=1, col="red")
plot(x = oes_res$PenMerged150$y_jpa, y = y$Y[test_ind],
     xlim = c(0,1), ylim = c(0,1), xlab = "y_hat", ylab = "",
     main = "PenFSVD + Wavelet reg")
abline(a=0, b=1, col="red")
plot(x = oes_res$BasMerged150$y_jpa, y = y$Y[test_ind],
     xlim = c(0,1), ylim = c(0,1), xlab = "y_hat", ylab = "",
     main = "BasisFSVD + Wavelet reg")
abline(a=0, b=1, col="red")
dev.off()



#####
par(mfrow = c(2,3))
plot(x = oes_res$b_bsp$y_jpa, y = y$Y[test_ind],
     ylim = c(0.2,0.8), xlab = "", ylab = "y",
     main = "(b) B-spline reg")
abline(a=0, b=1, col="red")
plot(x = oes_res$b_wv$y_jpa, y = y$Y[test_ind],
     ylim = c(0.2,0.8), xlab = "", ylab = "",
     main = "(b) Wavelet reg")
abline(a=0, b=1, col="red")

plot(x = oes_res$PenBsp$y_jpa, y = y$Y[test_ind],
     ylim = c(0.2,0.8), xlab = "y_hat", ylab = "y",
     main = "PenFSVD + B-spline reg")
abline(a=0, b=1, col="red")
plot(x = oes_res$BasBsp$y_jpa, y = y$Y[test_ind],
     ylim = c(0.2,0.8), xlab = "y_hat", ylab = "",
     main = "BasisFSVD + B-spline reg")
abline(a=0, b=1, col="red")
plot(x = oes_res$PenMerged150$y_jpa, y = y$Y[test_ind],
     ylim = c(0.2,0.8), xlab = "y_hat", ylab = "",
     main = "PenFSVD + Wavelet reg")
abline(a=0, b=1, col="red")
plot(x = oes_res$BasMerged150$y_jpa, y = y$Y[test_ind],
     ylim = c(0.2,0.8), xlab = "y_hat", ylab = "",
     main = "BasisFSVD + Wavelet reg")
abline(a=0, b=1, col="red")








