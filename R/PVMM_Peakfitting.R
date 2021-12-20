
spect_em_pvmm <- function(x, y,  mu, sigma, eta, mix_ratio, conv.cri, maxit) {

  #trancated Pseudo-Voigt
  truncated_pv <- function(x, mu, sigma, eta) {
    (eta*dcauchy(x, mu, sqrt(2*log(2))*sigma) + (1-eta)*dnorm(x, mu, sigma)) /
      sum(eta*dcauchy(x, mu, sqrt(2*log(2))*sigma) + (1-eta)*dnorm(x, mu, sigma))
  }

  #Defining each values
  start_cal     <- Sys.time()
  messe         <- "Not converged"
  N             <- length(x)
  LL_1          <- numeric(0)
  mix_ratio_1   <- numeric(0)
  sigma_1       <- numeric(0)
  mu_1          <- numeric(0)
  eta_1         <- numeric(0)
  n_k           <- numeric(0)
  K             <- length(mu)

  #log-Likelihood
  f_k <- function(i) {
    mix_ratio[i] * truncated_pv(x, mu[i], sigma[i], eta[i])
  }

  LL <- function(x, y, mu, sigma, eta, mix_ratio) {
    pL <- sapply(1:K,f_k)
    sum(y * log(apply(pL,1,sum)))
  }

  LL_1[1]        <- LL(x, y, mu, sigma, eta, mix_ratio)
  mu_1           <- rbind(mu_1, mu)
  sigma_1        <- rbind(sigma_1, sigma)
  eta_1          <- rbind(eta_1, eta)
  mix_ratio_1    <- rbind(mix_ratio_1, mix_ratio)

  #Q-function
  Q_fun <- function(x, w_k, mu, sigma, eta, mix_ratio) {
    w_k %*% (log(mix_ratio) + log(truncated_pv(x, mu, sigma, eta)))
  }

  #Starting ECM algorithm
  for(i in 1:maxit) {
    tmp <- sapply(1:K, f_k)
    den <- apply(tmp, 1, sum)
    w_k <- matrix(NA, nrow=K, ncol=N)

    for(j in 1:K) {
      w_k[j,] <- y * mix_ratio[j] * truncated_pv(x, mu[j], sigma[j], eta[j]) / den
    }

    n_k                    <- apply(w_k,1,sum)
    n_k[which(is.na(n_k))] <- 0

    #Hanger for each parameter
    mu_cal     <- c()
    sigma_cal  <- c()
    eta_cal    <- c()

    #Updating mix_ratio
    mix_ratio  <- n_k/sum(y)

    #Updating mu
    for(k in 1:K) {
      opt    <- optimize(Q_fun, interval = c(min(x), max(x)), tol = 1e-10, x = x, sigma = sigma[k], eta = eta[k], w_k = w_k[k,], mix_ratio = mix_ratio[k], maximum = TRUE)
      mu_cal <- c(mu_cal, opt$maximum)
    }
    mu       <- mu_cal

    #Updating sigma
    for(k in 1:K) {
      opt       <- optimize(Q_fun, interval = c(1e-3, 100), tol = 1e-10, x = x, mu = mu[k], eta = eta[k], w_k = w_k[k,], mix_ratio = mix_ratio[k], maximum = TRUE)
      sigma_cal <- c(sigma_cal, opt$maximum)
    }
    sigma       <- sigma_cal

    #Updating eta
    for(k in 1:K) {
      opt      <- optimize(Q_fun, interval = c(0, 1), tol = 1e-10, x = x, mu = mu[k], sigma = sigma[k], w_k = w_k[k,], mix_ratio = mix_ratio[k], maximum = TRUE)
      eta_cal  <- c(eta_cal, opt$maximum)
    }

    #Grid search eta
    # for(k in 1:K){
    #   eta_candi   <-  seq(0.3, 1.0, by = 0.1)
    #   Qfun_candi  <-  c()
    #
    #   for(j in 1:length(beta_candi)){
    #     Qfun_candi  <-  c(Qfun_candi, Q_fun(x = x, mu = mu[k], sigma = sigma[k], w_k = w_k[k,], mix_ratio = mix_ratio[k], eta_candi[j]))
    #   }
    #   #print(Qfun_candi)
    #   eta_cal <- c(eta_cal, eta_candi[Qfun_candi == max(Qfun_candi)])
    # }

    eta         <- eta_cal

    #Updating each value
    LL_1[i+1]    <- LL(x, y, mu, sigma, eta, mix_ratio)
    mu_1         <- rbind(mu_1, mu)
    sigma_1      <- rbind(sigma_1, sigma)
    eta_1        <- rbind(eta_1, eta)
    mix_ratio_1  <- rbind(mix_ratio_1, mix_ratio)

    #Convergence check
    if(abs(LL_1[i+1] - LL_1[i]) < conv.cri){ messe = "Converged"; break }
    print(LL_1[i+1] - LL_1[i])
  }

  end_cal        <- Sys.time()
  cal_time       <- difftime(end_cal, start_cal, units = "sec")

  #Hanger for result
  list(mu = mu, sigma = sigma, eta = eta, mix_ratio = mix_ratio, it = i, LL = LL_1,
       MU = mu_1, SIGMA = sigma_1, ETA = eta_1, MIX_RATIO = mix_ratio_1,
       convergence =messe, W_K = w_k, cal_time = cal_time)
}



