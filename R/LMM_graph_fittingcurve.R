
show_lmm_curve <- function(spect_em_lmm_res, x, y, mix_ratio_init, mu_init, gam_init) {

  #Normalized Lorentz
  dCauchy <- function(x, mu, gam) {
    (dcauchy(x, mu, gam)) / sum(dcauchy(x, mu, gam))
  }

  es_mu_L         <- spect_em_lmm_res$mu
  es_gam_L        <- spect_em_lmm_res$gam
  es_mix_ratio_L  <- spect_em_lmm_res$mix_ratio
  W_K             <- spect_em_lmm_res$W_K
  cal_time        <- spect_em_lmm_res$cal_time
  K               <- length(es_mu_L)

  #Plotting fitting curve
  #par(mfrow = c(1,1))

  esLMM_EM <- matrix(NA, nrow = K, ncol = length(x))
  for(k in 1:K) {
    esLMM_EM[k, ] <- es_mix_ratio_L[k] * dCauchy(x = x, mu = es_mu_L[k], gam = es_gam_L[k])
  }

  esLMM   <- colSums(esLMM_EM)

  esLMM_c   <- matrix(NA, nrow = K, ncol = length(x))
  for(k in 1:K) {
    esLMM_c[k,] <- esLMM_EM[k,] * sum(y) / sum(esLMM)
  }

  Est_spect_L    <- colSums(esLMM_c)
  Energy_L       <- x

  for(k in 1:K) {
    plot(esLMM_c[k,]~Energy_L, pch=19,  xlim = c(min(x), max(x)), ylim = c(0, max(y)), col = "black", las = 1, typ = "l", lwd = 2, cex = 0.75, lty = 3, ylab = "Frequency", xlab = "x", main = "")
    par(new = TRUE)
  }

  plot(Est_spect_L~Energy_L, pch=19,  xlim = c(min(x), max(x)), ylim = c(0, max(y)), col = "tomato3", las = 1, typ = "l", lwd = 2, cex = 0.75, ylab = "Frequency", xlab = "x", main = "Estimted LMM")
  par(new = TRUE)
  plot(y~Energy_L, pch=21,  xlim = c(min(x), max(x)), ylim = c(0, max(y)), col = "black", las = 1, typ = "p", lwd = 2, cex = 0.75, ylab = "Frequency", xlab = "x")

  #Plotting trace plot of the parameters
  cols <- colorRamp(c("#ff4b00", "#fff100", "#03af7a", "#005aff"))
  #par(mfrow = c(3,2))

  esLMM_EM <- matrix(NA, nrow = K, ncol = length(x))
  for(k in 1:K) {
    esLMM_EM[k, ] <- es_mix_ratio_L[k] * dCauchy(x = x, mu = es_mu_L[k], gam = es_gam_L[k])
  }

  esLMM   <- colSums(esLMM_EM)

  esLMM_c   <- matrix(NA, nrow = K, ncol = length(x))
  for(k in 1:K) {
    esLMM_c[k,] <- esLMM_EM[k,] * sum(y) / sum(esLMM)
  }

  Est_spect_L    <- colSums(esLMM_c)
  x_L            <- x

  for(k in 1:K) {
    plot(esLMM_c[k,]~x_L, pch=19,  xlim = c(min(x), max(x)), ylim = c(0, max(y)), col = "black", las = 1, typ = "l", lwd = 2, cex = 0.75, lty = 3, ylab = "Frequency", xlab = "x", main = "")
    par(new = TRUE)
  }

  plot(Est_spect_L~x_L, pch=19,  xlim = c(min(x), max(x)), ylim = c(0, max(y)), col = "tomato3", las = 1, typ = "l", lwd = 2, cex = 0.75, ylab = "Frequency", xlab = "x", main = "Estimted LMM")
  par(new = TRUE)
  plot(y~Energy_L, pch=21,  xlim = c(min(x), max(x)), ylim = c(0, max(y)), col = "black", las = 1, typ = "p", lwd = 2, cex = 0.75, ylab = "Frequency", xlab = "x")

  plot(spect_em_lmm_res$LL, ylab = "L", xlab = "iteration")

  if(K > 1) {
    for(k in 1:K) {
      plot(spect_em_lmm_res$MU[,k], type = "l", col = rgb(cols((k-1) / (K-1)) / 255), xlim = c(0, length(spect_em_lmm_res$LL)), ylim = c(min(x), max(x)), lwd = 2, ylab = "mu", xlab = "iteration", las = 1)
      par(new = TRUE)
    }
    par(new = FALSE)

    for(k in 1:K) {
      plot(spect_em_lmm_res$GAM[,k], type = "l", col = rgb(cols((k-1) / (K-1))/255), xlim = c(0, length(spect_em_lmm_res$LL)), ylim = c(0, max(spect_em_lmm_res$GAM)), lwd = 2, ylab = "gam", xlab = "iteration", las = 1)
      par(new = TRUE)
    }
    par(new = FALSE)

    for(k in 1:K) {
      plot(spect_em_lmm_res$MIX_RATIO[,k], type = "l", col = rgb(cols((k-1) / (K-1))/255), xlim = c(0, length(spect_em_lmm_res$LL)), ylim = c(0, max(spect_em_lmm_res$MIX_RATIO)),lwd = 2, ylab="mix_ratio", xlab = "iteration", las=1)
      par(new = TRUE)
    }

    plot(c(0,0)~c(spect_em_lmm_res$it, 0), type = "l", lty = 2, ylab = "", xlab = "", xlim = c(0, length(spect_em_lmm_res$LL)), ylim = c(0, max(spect_em_lmm_res$MIX_RATIO)), las=1)
    par(new = FALSE)

  } else {

    plot(spect_em_lmm_res$MU, type = "l", col = "red", xlim = c(0, length(spect_em_lmm_res$LL)), ylim = c(min(x), max(x)), lwd = 2, ylab = "mu", xlab = "iteration", las = 1)
    plot(spect_em_lmm_res$GAM, type = "l", col = "red", xlim = c(0, length(spect_em_lmm_res$LL)), ylim = c(0, max(spect_em_lmm_res$GAM)), lwd = 2, ylab = "gam", xlab = "iteration", las = 1)
    plot(spect_em_lmm_res$MIX_RATIO, type = "l", col = "red", xlim = c(0, length(spect_em_lmm_res$LL)), ylim = c(0, max(spect_em_lmm_res$MIX_RATIO)),lwd = 2, ylab="mix_ratio", xlab = "iteration", las=1)
    plot(c(0,0)~c(spect_em_lmm_res$it, 0), type = "l", lty = 2, ylab = "", xlab = "", xlim = c(0, length(spect_em_lmm_res$LL)), ylim = c(0, max(spect_em_lmm_res$MIX_RATIO)), las=1)
    par(new = FALSE)
  }

}




