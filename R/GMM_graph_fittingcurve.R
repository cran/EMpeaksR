
show_gmm_curve <- function(spect_em_gmm_res, x, y, mix_ratio_init, mu_init, sigma_init) {

  es_mu_G         <- spect_em_gmm_res$mu
  es_sigma_G      <- spect_em_gmm_res$sigma
  es_mix_ratio_G  <- spect_em_gmm_res$mix_ratio
  W_K             <- spect_em_gmm_res$W_K
  cal_time        <- spect_em_gmm_res$cal_time
  K               <- length(es_mu_G)

  #Fitting curve
  cols <- colorRamp(c("#ff4b00", "#fff100", "#03af7a", "#005aff"))

  esGMM_EM <- matrix(NA, nrow = K, ncol = length(x))
  for(k in 1:K) {
    esGMM_EM[k, ] <- es_mix_ratio_G[k] * dnorm(x = x, mean = es_mu_G[k], sd = es_sigma_G[k])
  }

  esGMM   <- colSums(esGMM_EM)

  esGMM_c   <- matrix(NA, nrow = K, ncol = length(x))
  for(k in 1:K) {
    esGMM_c[k,] <- esGMM_EM[k,] * sum(y) / sum(esGMM)
  }

  Est_spect_G    <- colSums(esGMM_c)
  x_G            <- x

  for(k in 1:K) {
    plot(esGMM_c[k,]~x_G, pch=19,  xlim = c(min(x), max(x)), ylim = c(0, max(y)), col = "black", las = 1, typ = "l", lwd = 2, cex = 0.75, lty = 3, ylab = "Frequency", xlab = "x", main = "")
    par(new = TRUE)
  }

  plot(Est_spect_G~x_G, pch=19,  xlim = c(min(x), max(x)), ylim = c(0, max(y)), col = "tomato3", las = 1, typ = "l", lwd = 2, cex = 0.75, ylab = "Frequency", xlab = "x", main = "Estimted GMM")
  par(new = TRUE)
  plot(y~x_G, pch=21,  xlim = c(min(x), max(x)), ylim = c(0, max(y)), col = "black", las = 1, typ = "p", lwd = 2, cex = 0.75, ylab = "Frequency", xlab = "x")

  plot(spect_em_gmm_res$LL, ylab = "L", xlab = "iteration")

if(K > 1) {
  for(k in 1:K) {
    plot(spect_em_gmm_res$MU[,k], type = "l", col = rgb(cols((k-1) / (K-1)) / 255), xlim = c(0, length(spect_em_gmm_res$LL)), ylim = c(min(x), max(x)), lwd = 2, ylab = "mu", xlab = "iteration", las = 1)
    par(new = TRUE)
  }
  par(new = FALSE)

  for(k in 1:K) {
    plot(spect_em_gmm_res$SIGMA[,k], type = "l", col = rgb(cols((k-1) / (K-1))/255), xlim = c(0, length(spect_em_gmm_res$LL)), ylim = c(0, max(spect_em_gmm_res$SIGMA)), lwd = 2, ylab = "sigma", xlab = "iteration", las = 1)
    par(new = TRUE)
  }
  par(new = FALSE)

  for(k in 1:K) {
    plot(spect_em_gmm_res$MIX_RATIO[,k], type = "l", col = rgb(cols((k-1) / (K-1))/255), xlim = c(0, length(spect_em_gmm_res$LL)), ylim = c(0, max(spect_em_gmm_res$MIX_RATIO)),lwd = 2, ylab="mix_ratio", xlab = "iteration", las=1)
    par(new = TRUE)
  }

  plot(c(0,0)~c(spect_em_gmm_res$it, 0), type = "l", lty = 2, ylab = "", xlab = "", xlim = c(0, length(spect_em_gmm_res$LL)), ylim = c(0, max(spect_em_gmm_res$MIX_RATIO)), las=1)
  par(new = FALSE)

} else {

    plot(spect_em_gmm_res$MU, type = "l", col = "red", xlim = c(0, length(spect_em_gmm_res$LL)), ylim = c(min(x), max(x)), lwd = 2, ylab = "mu", xlab = "iteration", las = 1)
    plot(spect_em_gmm_res$SIGMA, type = "l", col = "red", xlim = c(0, length(spect_em_gmm_res$LL)), ylim = c(0, max(spect_em_gmm_res$SIGMA)), lwd = 2, ylab = "sigma", xlab = "iteration", las = 1)
    plot(spect_em_gmm_res$MIX_RATIO, type = "l", col = "red", xlim = c(0, length(spect_em_gmm_res$LL)), ylim = c(0, max(spect_em_gmm_res$MIX_RATIO)),lwd = 2, ylab="mix_ratio", xlab = "iteration", las=1)
    plot(c(0,0)~c(spect_em_gmm_res$it, 0), type = "l", lty = 2, ylab = "", xlab = "", xlim = c(0, length(spect_em_gmm_res$LL)), ylim = c(0, max(spect_em_gmm_res$MIX_RATIO)), las=1)
    par(new = FALSE)
    }

}

