
show_pvmm_lback_curve <-function(spect_em_pvmm_lback_res, x, y, mix_ratio_init, mu_init, sigma_init, eta_init, x_lower = min(x), x_upper = max(x)){

  #trancated pseudo-Voigt
  truncated_pv <- function(x, mu, sigma, eta) {
    (eta*dcauchy(x, mu, sqrt(2*log(2))*sigma) + (1-eta)*dnorm(x, mu, sigma)) /
      sum(eta*dcauchy(x, mu, sqrt(2*log(2))*sigma) + (1-eta)*dnorm(x, mu, sigma))
  }

  #Uniform function
  ramp_ind_uni <- function(x, x_upper) {
    Step_fact <- as.numeric(x >= x_upper)
    return(Step_fact)
  }

  es_mu_PV          <- spect_em_pvmm_lback_res$mu
  es_sigma_PV       <- spect_em_pvmm_lback_res$sigma
  es_eta_PV         <- spect_em_pvmm_lback_res$eta
  es_mix_ratio_PV   <- spect_em_pvmm_lback_res$mix_ratio
  W_K               <- spect_em_pvmm_lback_res$W_K
  cal_time          <- spect_em_pvmm_lback_res$cal_time

  K     <- length(es_mu_PV)
  M     <- K + length(x_lower) + 2
  x_N   <- x[length(x)]

  #Calculating  fitting curve
  for(i in 1:M) {
    plot(spect_em_pvmm_lback_res$W_K[i,]~x, main = (paste0("comp", i)), xlim = c(min(x), max(x)), ylab = "intensity", ylim = c(0, max(y)), pch = 19, cex = 0.25)
  }

  plot(colSums(spect_em_pvmm_lback_res$W_K[1:K,])~x, main = "comp_peaks", xlim = c(min(x), max(x)), ylab = "intensity", ylim = c(0, max(y)),  pch = 19, cex = 0.25)
  plot(colSums(spect_em_pvmm_lback_res$W_K[(K+1):M,])~x, main = "comp_background", xlim = c(min(x), max(x)), ylab = "intensity", ylim = c(0, max(y)),  pch = 19, cex = 0.25)


  Es_PVMM   <- matrix(NA, ncol = K, nrow = length(x))
  for(i in 1:K) {
    Es_PVMM[,i] <- es_mix_ratio_PV[i] * truncated_pv(x = x, mu = es_mu_PV[i], sigma = es_sigma_PV[i], eta = es_eta_PV[i])
  }


  Es_background       <- matrix(NA, ncol = (M-K), nrow = length(x))
  Es_background[,1]   <- es_mix_ratio_PV[K+1] * ((1 / (x_N - x[1])) * ramp_ind_uni(x, x[1])) / sum(((1 / (x_N - x[1])) * ramp_ind_uni(x, x[1])))
  Es_background[,2]   <- es_mix_ratio_PV[K+2] * (2*(x - x_lower)/(x_upper - x_lower)^2) / sum(2*(x - x_lower)/(x_upper - x_lower)^2)
  Es_background[,3]   <- es_mix_ratio_PV[M] * (2*(x_upper - x)/(x_upper - x_lower)^2) / sum(2*(x_upper - x)/(x_upper - x_lower)^2)
  Es_spect            <- rowSums(Es_PVMM) +  rowSums(Es_background)


  cols_K <- colorRamp(c("#cc6a70", "#de7065", "#ed8055", "#f68f46", "#f9a242", "#f9b641", "#f7cb44", "#efe350"))
  K_plot <- function(K_val, K_upper = K, K_lower = 1) {
    conb_K <- (K_val - K_lower) / (K_upper - K_lower)
    RGB_K  <- rgb(cols_K(conb_K)/255)
    return(RGB_K)
  }

  cols_B <- colorRamp(c("#187FC4", "#005293", "#001A43"))
  B_plot <- function(B_val, B_upper = (M-K), B_lower = 1) {
    conb_B <- (B_val - B_lower) / (B_upper - B_lower)
    RGB_B  <- rgb(cols_B(conb_B) / 255)
    return(RGB_B)
  }

  scale_cont   <-  sum(y)

  plot(y ~ x, las = 1, typ = "p", lwd = 2, cex = 0.75, lty = 3, ylab = "Frequency", xlab = "x", main = "", xlim = c(min(x), max(x)), ylim = c(0, max(y)), col = "black")
  par(new = TRUE)

  for(i in 1:K) {
    plot(scale_cont * Es_PVMM[,i] ~ x, las = 1, typ = "l", lwd = 2, cex = 0.75, lty = 3, ylab = "Frequency", xlab = "x", main = "", xlim = c(min(x), max(x)), ylim = c(0, max(y)), col = K_plot(i))
    par(new = TRUE)
  }

  plot(scale_cont * (Es_background[,1]) ~ x, las = 1, typ = "l", lwd = 2, cex = 0.75, lty = 3, ylab = "Frequency", xlab = "x", main = "", xlim = c(min(x), max(x)), ylim = c(0, max(y)), col = B_plot(1))
  par(new = TRUE)

  for(i in 2:(M-K)){
    plot(scale_cont * rowSums(Es_background[,1:i]) ~ x, las = 1, typ = "l", lwd = 2, cex = 0.75, lty = 3, ylab = "Frequency", xlab = "x", main = "", xlim = c(min(x), max(x)), ylim = c(0, max(y)), col = B_plot(i))
    par(new = TRUE)
  }

  plot(scale_cont * Es_spect ~ x, las = 1, typ = "l", lwd = 2, cex = 0.75, lty = 1, ylab = "Frequency", xlab = "x", main = "Estimated spect", xlim = c(min(x), max(x)), ylim = c(0, max(y)), col = "tomato3")

  #Trace plot
  plot(spect_em_pvmm_lback_res$LL, ylab = "L", xlab = "iteration")

  for(i in 1:K) {
    plot(spect_em_pvmm_lback_res$MU[,i], type = "l", col = K_plot(i), xlim = c(0, length(spect_em_pvmm_lback_res$LL)), ylim = c(min(x), max(x)), lwd = 2, ylab = "mu", xlab = "iteration", las = 1)
    par(new = TRUE)
  }
  par(new = FALSE)

  for(i in 1:K) {
    plot(spect_em_pvmm_lback_res$SIGMA[,i], type = "l", col = K_plot(i), xlim = c(0, length(spect_em_pvmm_lback_res$LL)), ylim = c(0, max(spect_em_pvmm_lback_res$SIGMA)), lwd = 2, ylab = "sigma", xlab = "iteration", las = 1)
    par(new = TRUE)
  }
  par(new = FALSE)

  for(i in 1:K) {
    plot(spect_em_pvmm_lback_res$ETA[,i], type = "l", col = K_plot(i), xlim = c(0, length(spect_em_pvmm_lback_res$LL)), ylim = c(0, max(spect_em_pvmm_lback_res$ETA)), lwd = 2, ylab = "eta", xlab = "iteration", las = 1)
    par(new = TRUE)
  }
  par(new = FALSE)

  for(i in 1:K) {
    plot(spect_em_pvmm_lback_res$MIX_RATIO[,i], type = "l", col = K_plot(i), xlim = c(0, length(spect_em_pvmm_lback_res$LL)), ylim = c(0, max(spect_em_pvmm_lback_res$MIX_RATIO)), lwd = 2, ylab = "mix_ratio", xlab = "iteration", las = 1)
    par(new = TRUE)
  }
  par(new = TRUE)

  for(i in (K+1):M) {
    plot(spect_em_pvmm_lback_res$MIX_RATIO[,i], type = "l", col = B_plot(i-K), xlim = c(0, length(spect_em_pvmm_lback_res$LL)), ylim = c(0, max(spect_em_pvmm_lback_res$MIX_RATIO)), lwd = 2, ylab = "", xlab = "", las = 1)
    par(new = TRUE)
  }
  par(new = FALSE)

}
