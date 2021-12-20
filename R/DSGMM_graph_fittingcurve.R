show_dsgmm_curve <- function(spect_em_dsgmm_res, x, y, mix_ratio_init, mu_init, sigma_init, alpha_init, eta_init) {

  #trancated Doniach-Sunjic-Gauss
  truncated_dsg <- function(x, mu, sigma, alpha, eta) {
    ((eta*(((gamma(1-alpha)) / ((x-mu)^2+(sqrt(2*log(2))*sigma)^2)^((1-alpha)/2))*cos((pi*alpha/2)+(1-alpha)*atan((x-mu)/(sqrt(2*log(2))*sigma))))) + (1-eta)*dnorm(x, mu, sigma)) /
      sum( ((eta*(((gamma(1-alpha)) / ((x-mu)^2+(sqrt(2*log(2))*sigma)^2)^((1-alpha)/2))*cos((pi*alpha/2)+(1-alpha)*atan((x-mu)/(sqrt(2*log(2))*sigma))))) + (1-eta)*dnorm(x, mu, sigma)))
  }

  es_mu_DSG         <- spect_em_dsgmm_res$mu
  es_sigma_DSG      <- spect_em_dsgmm_res$sigma
  es_alpha_DSG      <- spect_em_dsgmm_res$alpha
  es_eta_DSG       <- spect_em_dsgmm_res$eta
  es_mix_ratio_DSG  <- spect_em_dsgmm_res$mix_ratio
  es_mu             <- es_mu_DSG + es_sigma_DSG / tan(pi / (2-es_alpha_DSG))
  W_K               <- spect_em_dsgmm_res$W_K
  cal_time          <- spect_em_dsgmm_res$cal_time
  K                 <- length(es_mu_DSG)

  #Fitting curve
  cols <- colorRamp(c("#ff4b00", "#fff100", "#03af7a", "#005aff"))

  esDSGMM_EM <- matrix(NA, nrow = K, ncol = length(x))
  for(k in 1:K) {
    esDSGMM_EM[k, ] <- es_mix_ratio_DSG[k] * truncated_dsg(x = x, mu = es_mu_DSG[k], sigma = es_sigma_DSG[k], alpha = es_alpha_DSG[k], eta = es_eta_DSG[k])
  }

  esDSGMM   <- colSums(esDSGMM_EM)

  esDSGMM_c   <- matrix(NA, nrow = K, ncol = length(x))
  for(k in 1:K) {
    esDSGMM_c[k,] <- esDSGMM_EM[k,] * sum(y) / sum(esDSGMM)
  }

  Est_spect_DSG   <- colSums(esDSGMM_c)
  x_DSG           <- x

  for(k in 1:K) {
    plot(esDSGMM_c[k,]~x_DSG, pch=19,  xlim = c(min(x), max(x)), ylim = c(0, max(y)), col = "black", las = 1, typ = "l", lwd = 2, cex = 0.75, lty = 3, ylab = "Frequency", xlab = "x", main = "")
    par(new = TRUE)
  }

  plot(Est_spect_DSG~x_DSG, pch=19,  xlim = c(min(x), max(x)), ylim=c(0, max(y)), col = "tomato3", las = 1, typ = "l", lwd = 2, cex = 0.75, ylab = "Frequency", xlab = "x", main = "Estimted DSGMM")
  par(new = TRUE)
  plot(y~x_DSG, pch=21, xlim = c(min(x), max(x)), ylim = c(0, max(y)), col = "black", las = 1, typ = "p", lwd = 2, cex = 0.75, ylab = "Frequency", xlab = "x")

  plot(spect_em_dsgmm_res$LL, ylab = "L", xlab = "iteration")

  if(K > 1) {
  for(k in 1:K) {
    plot(spect_em_dsgmm_res$MU[,k], type = "l", col = rgb(cols((k-1) / (K-1)) / 255), xlim = c(0, length(spect_em_dsgmm_res$LL)), ylim = c(0, 100), lwd = 2, ylab = "mu", xlab = "iteration", las = 1)
    par(new = TRUE)
  }
  par(new = FALSE)

  for(k in 1:K) {
    plot(spect_em_dsgmm_res$SIGMA[,k], type = "l", col = rgb(cols((k-1) / (K-1)) / 255), xlim = c(0, length(spect_em_dsgmm_res$LL)), ylim = c(0, max(spect_em_dsgmm_res$SIGMA)), lwd = 2, ylab = "sigma", xlab = "iteration", las = 1)
    par(new = TRUE)
  }
  par(new = FALSE)

  for(k in 1:K) {
    plot(spect_em_dsgmm_res$ALPHA[,k], type = "l", col = rgb(cols((k-1) / (K-1)) / 255), xlim = c(0, length(spect_em_dsgmm_res$LL)), ylim = c(0, 1), lwd = 2, ylab = "alpha", xlab = "iteration", las = 1)
    par(new = TRUE)
  }
  par(new = FALSE)

  for(k in 1:K) {
    plot(spect_em_dsgmm_res$ETA[,k], type = "l", col = rgb(cols((k-1) / (K-1)) / 255), xlim = c(0, length(spect_em_dsgmm_res$LL)), ylim = c(0, 1), lwd = 2, ylab = "eta", xlab = "iteration", las = 1)
    par(new = TRUE)
  }
  par(new = FALSE)

  for(k in 1:K) {
    plot(spect_em_dsgmm_res$MIX_RATIO[,k], type = "l", col = rgb(cols((k-1) / (K-1)) / 255), xlim = c(0, length(spect_em_dsgmm_res$LL)), ylim = c(0, max(spect_em_dsgmm_res$MIX_RATIO)),lwd = 2, ylab = "mix_ratio", xlab = "iteration", las = 1)
    par(new = TRUE)
  }
  plot(c(0,0)~c(spect_em_dsgmm_res$it, 0), type = "l", lty = 2, ylab = "", xlab = "", xlim = c(0, length(spect_em_dsgmm_res$LL)), ylim = c(0, max(spect_em_dsgmm_res$MIX_RATIO)), las = 1)
  par(new = FALSE)

  } else {
    plot(spect_em_dsgmm_res$MU, type = "l", col = "red", xlim = c(0, length(spect_em_dsgmm_res$LL)), ylim = c(min(x), max(x)), lwd = 2, ylab = "mu", xlab = "iteration", las = 1)
    plot(spect_em_dsgmm_res$SIGMA, type = "l", col = "red", xlim = c(0, length(spect_em_dsgmm_res$LL)), ylim = c(0, max(spect_em_dsgmm_res$SIGMA)), lwd = 2, ylab = "sigma", xlab = "iteration", las = 1)
    plot(spect_em_dsgmm_res$ALPHA, type = "l", col = "red", xlim = c(0, length(spect_em_dsgmm_res$LL)), ylim = c(0, max(spect_em_dsgmm_res$ALPHA)), lwd = 2, ylab = "beta", xlab = "iteration", las = 1)
    plot(spect_em_dsgmm_res$ETA, type = "l", col = "red", xlim = c(0, length(spect_em_dsgmm_res$LL)), ylim = c(0, max(spect_em_dsgmm_res$ETA)), lwd = 2, ylab = "eta", xlab = "iteration", las = 1)
    plot(spect_em_dsgmm_res$MIX_RATIO, type = "l", col = "red", xlim = c(0, length(spect_em_dsgmm_res$LL)), ylim = c(0, max(spect_em_dsgmm_res$MIX_RATIO)),lwd = 2, ylab="mix_ratio", xlab = "iteration", las=1)
    plot(c(0,0)~c(spect_em_dsgmm_res$it, 0), type = "l", lty = 2, ylab = "", xlab = "", xlim = c(0, length(spect_em_dsgmm_res$LL)), ylim = c(0, max(spect_em_dsgmm_res$MIX_RATIO)), las=1)
    par(new = FALSE)
  }

}


