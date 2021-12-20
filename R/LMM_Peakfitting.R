
spect_em_lmm  <- function(x, y, mu, gam, mix_ratio, conv.cri, maxit) {
  #Normalized Lorentz
  dCauchy <- function(x, mu, gam) {
    (dcauchy(x, mu, gam)) / sum(dcauchy(x, mu, gam))
  }

  #Defining each value
  start_cal     <- Sys.time()
  messe         <-  "Not converged"
  N             <- length(x)
  LL_1          <- numeric(0)
  mix_ratio_1   <- numeric(0)
  gam_1         <- numeric(0)
  mu_1          <- numeric(0)
  n_k           <- numeric(0)
  K             <- length(mu)

  #log Likelihood
  f_k <- function(i) {
    mix_ratio[i]*dCauchy(x, mu[i], gam[i])
    }

  LL <- function(x, y, mu, gam, mix_ratio) {
    pL <- sapply(1:K,f_k)
    sum(y*log(apply(pL,1,sum)))
  }

  LL_1[1]        <- LL(x, y, mu, gam, mix_ratio)
  mu_1           <- rbind(mu_1, mu)
  gam_1          <- rbind(gam_1, gam)
  mix_ratio_1    <- rbind(mix_ratio_1, mix_ratio)

  #Q function
  Q_fun <- function(x, w_k, mu, gam, mix_ratio) {
    w_k %*% (log(mix_ratio) + log(dCauchy(x, mu, gam)))
  }

  #Starting ECM algorithm
  for(i in 1:maxit) {
    tmp <- sapply(1:K, f_k)
    den <- apply(tmp, 1, sum)
    w_k <- matrix(NA, nrow=K, ncol=N)

    for(j in 1:K) {
      w_k[j,] <- y * mix_ratio[j]*dCauchy(x, mu[j], gam[j])/den
    }

    n_k <- apply(w_k,1,sum)
    n_k[which(is.na(n_k))] <- 0

    #Hanger for mu and gam
    mu_cal     <- c()
    gam_cal    <- c()

    #UPdating mix_ratio
    mix_ratio  <- n_k/sum(y)

    #Updating mu
    for(k in 1:K) {
      opt    <- optimize(Q_fun, interval = c(min(x), max(x)), tol = 1e-10, x = x, gam = gam[k], w_k = w_k[k,], mix_ratio = mix_ratio[k], maximum = TRUE)
      mu_cal <- c(mu_cal, opt$maximum)
    }
    mu   <- mu_cal

    #Updating gam
    for(k in 1:K) {
      opt     <- optimize(Q_fun, interval = c(1e-3, 100), tol = 1e-10, x = x, mu = mu[k], w_k = w_k[k,], mix_ratio = mix_ratio[k], maximum=TRUE)
      gam_cal <- c(gam_cal, opt$maximum)
    }
    gam <- gam_cal

    #Updating each value
    LL_1[i+1]    <- LL(x, y, mu, gam, mix_ratio)
    mu_1         <- rbind(mu_1, mu)
    gam_1        <- rbind(gam_1, gam)
    mix_ratio_1  <- rbind(mix_ratio_1, mix_ratio)

    if(abs(LL_1[i+1]-LL_1[i]) < conv.cri) {
      messe = "Converged"; break
      }
    print(LL_1[i+1]-LL_1[i])
  }

  end_cal  <- Sys.time()
  cal_time <- difftime(end_cal, start_cal, units = "sec")

  #Hanger for result
  list(mu = mu, gam = gam, mix_ratio = mix_ratio, it = i, LL = LL_1,
       MU = mu_1, GAM = gam_1, MIX_RATIO = mix_ratio_1, convergence =messe, W_K = w_k, cal_time = cal_time)
}



