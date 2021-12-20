beta
show_pvmm_curve <- function(spect_em_pvmm_res, x, y, mix_ratio_init, mu_init, sigma_init, eta_init) {

  #trancated Pseudo-Voigt
  truncated_pv <- function(x, mu, sigma, eta) {
    (eta*dcauchy(x, mu, sqrt(2*log(2))*sigma) + (1-eta)*dnorm(x, mu, sigma)) /
      sum(eta*dcauchy(x, mu, sqrt(2*log(2))*sigma) + (1-eta)*dnorm(x, mu, sigma))
  }

  es_mu_PV         <- spect_em_pvmm_res$mu
  es_sigma_PV      <- spect_em_pvmm_res$sigma
  es_eta_PV        <- spect_em_pvmm_res$eta
  es_mix_ratio_PV  <- spect_em_pvmm_res$mix_ratio
  W_K              <- spect_em_pvmm_res$W_K
  cal_time         <- spect_em_pvmm_res$cal_time
  K                <- length(es_mu_PV)

  #Fitting curve
  cols <- colorRamp(c("#ff4b00", "#fff100", "#03af7a", "#005aff"))

  esPVMM_EM <- matrix(NA, nrow = K, ncol = length(x))
  for(k in 1:K) {
    esPVMM_EM[k, ] <- es_mix_ratio_PV[k] * truncated_pv(x = x, mu = es_mu_PV[k], sigma = es_sigma_PV[k], eta = es_eta_PV[k])
  }

  esPVMM    <- colSums(esPVMM_EM)
  esPVMM_c  <- matrix(NA, nrow = K, ncol = length(x))

  for(k in 1:K) {
    esPVMM_c[k,] <- esPVMM_EM[k,] * sum(y) / sum(esPVMM)
  }

  Est_spect_PV <- colSums(esPVMM_c)
  x_PV         <- x

  for(k in 1:K) {
    plot(esPVMM_c[k,]~x_PV, pch = 19,  xlim = c(min(x), max(x)), ylim = c(0, max(y)), col = "black", las = 1, typ = "l", lwd = 2, cex = 0.75, lty = 3, ylab = "Frequency", xlab = "x", main = "")
    par(new = TRUE)
  }

  plot(Est_spect_PV~x_PV, pch = 19,  xlim = c(min(x), max(x)), ylim = c(0, max(y)), col = "tomato3", las = 1, typ = "l", lwd = 2, cex = 0.75, ylab = "Frequency", xlab = "x", main = "Estimted PVMM")
  par(new = TRUE)
  plot(y~x_PV, pch=21,  xlim = c(min(x), max(x)), ylim = c(0, max(y)), col = "black", las = 1, typ = "p", lwd = 2, cex = 0.75, ylab = "Frequency", xlab = "x")

  plot(spect_em_pvmm_res$LL, ylab = "L", xlab = "iteration")

  if(K > 1) {
  for(k in 1:K) {
    plot(spect_em_pvmm_res$MU[,k], type = "l", col = rgb(cols((k-1) / (K-1)) / 255), xlim = c(0, length(spect_em_pvmm_res$LL)), ylim = c(min(x), max(x)), lwd = 2, ylab = "mu", xlab = "iteration", las = 1)
    par(new = TRUE)
  }
  par(new = FALSE)

  for(k in 1:K) {
    plot(spect_em_pvmm_res$SIGMA[,k], type = "l", col = rgb(cols((k-1) / (K-1)) / 255), xlim = c(0, length(spect_em_pvmm_res$LL)), ylim = c(0, max(spect_em_pvmm_res$SIGMA)), lwd = 2, ylab = "sigma", xlab = "iteration", las = 1)
    par(new = TRUE)
  }
  par(new = FALSE)

  for(k in 1:K) {
    plot(spect_em_pvmm_res$ETA[,k], type = "l", col = rgb(cols((k-1) / (K-1)) / 255), xlim = c(0, length(spect_em_pvmm_res$LL)), ylim = c(0, 1), lwd = 2, ylab = "eta", xlab = "iteration", las = 1)
    par(new = TRUE)
  }
  par(new = FALSE)

  for(k in 1:K) {
    plot(spect_em_pvmm_res$MIX_RATIO[,k], type = "l", col = rgb(cols((k-1) / (K-1)) / 255), xlim = c(0, length(spect_em_pvmm_res$LL)), ylim = c(0, max(spect_em_pvmm_res$MIX_RATIO)), lwd = 2, ylab = "mix_ratio", xlab = "iteration", las = 1)
    par(new = TRUE)
  }
  plot(c(0,0)~c(spect_em_pvmm_res$it, 0), type = "l", lty = 2, ylab = "", xlab = "", xlim = c(0, length(spect_em_pvmm_res$LL)), ylim = c(0, max(spect_em_pvmm_res$MIX_RATIO)), las = 1)
  par(new = FALSE)

  } else {
    plot(spect_em_pvmm_res$MU, type = "l", col = "red", xlim = c(0, length(spect_em_pvmm_res$LL)), ylim = c(min(x), max(x)), lwd = 2, ylab = "mu", xlab = "iteration", las = 1)
    plot(spect_em_pvmm_res$SIGMA, type = "l", col = "red", xlim = c(0, length(spect_em_pvmm_res$LL)), ylim = c(0, max(spect_em_pvmm_res$SIGMA)), lwd = 2, ylab = "sigma", xlab = "iteration", las = 1)
    plot(spect_em_pvmm_res$ETA, type = "l", col = "red", xlim = c(0, length(spect_em_pvmm_res$LL)), ylim = c(0, max(spect_em_pvmm_res$ETA)), lwd = 2, ylab = "eta", xlab = "iteration", las = 1)
    plot(spect_em_pvmm_res$MIX_RATIO, type = "l", col = "red", xlim = c(0, length(spect_em_pvmm_res$LL)), ylim = c(0, max(spect_em_pvmm_res$MIX_RATIO)),lwd = 2, ylab="mix_ratio", xlab = "iteration", las=1)
    plot(c(0,0)~c(spect_em_pvmm_res$it, 0), type = "l", lty = 2, ylab = "", xlab = "", xlim = c(0, length(spect_em_pvmm_res$LL)), ylim = c(0, max(spect_em_pvmm_res$MIX_RATIO)), las=1)
    par(new = FALSE)
  }

}

