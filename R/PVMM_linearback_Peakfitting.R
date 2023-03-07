
spect_em_pvmm_lback <- function(x, y, mu, sigma, eta, mix_ratio, x_lower = min(x), x_upper = max(x), conv.cri, maxit) {

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

  NFun      <- length(x_lower)
  K         <- length(mu)
  M         <- K + length(x_lower) + 2
  x_N       <- x[length(x)]

  #log-likelihood
  f_k <- function(i) {
    mix_ratio[i] * truncated_pv(x, mu[i], sigma[i], eta[i])
  }

  LL  <- function(x, y, mu, sigma, eta, mix_ratio, x_lower, x_upper, x_N) {
    LL_back      <- matrix(NA, ncol=(M-K), nrow=N)
    LL_back[,1]  <- mix_ratio[K+1] * ((1/(x_N - x[1])) * ramp_ind_uni(x, x[1])) / sum(((1/(x_N - x[1])) * ramp_ind_uni(x, x[1])))
    LL_back[,2]  <- mix_ratio[K+2] * (2*(x - x_lower)/(x_upper - x_lower)^2) / sum(2*(x - x_lower)/(x_upper - x_lower)^2)
    LL_back[,3]  <- mix_ratio[M] * (2*(x_upper - x)/(x_upper - x_lower)^2) / sum(2*(x_upper - x)/(x_upper - x_lower)^2)

    pL <- cbind(sapply(1:(K), f_k),  LL_back)
    sum(y * log(apply(pL, 1, sum)))
  }

  LL_1[1]        <- LL(x, y, mu, sigma, eta, mix_ratio, x_lower, x_upper, x_N)
  mix_ratio_1    <- rbind(mix_ratio_1, mix_ratio)
  mu_1           <- rbind(mu_1, mu)
  sigma_1        <- rbind(sigma_1, sigma)
  eta_1          <- rbind(eta_1, eta)

  #Q-function
  Q_fun <- function(x, w_k, mu, sigma, eta, mix_ratio) {
    w_k %*% (log(mix_ratio) + log(truncated_pv(x, mu, sigma, eta)))
  }

  #Strating ECM algorithm
  for(i in 1:maxit) {
    tmp_back      <- matrix(NA, ncol=(M-K), nrow=N)
    tmp_back[,1]  <- mix_ratio[K+1] * ((1/(x_N - x[1])) * ramp_ind_uni(x, x[1])) / sum(((1/(x_N - x[1])) * ramp_ind_uni(x, x[1])))
    tmp_back[,2]  <- mix_ratio[K+2] * (2*(x - x_lower)/(x_upper - x_lower)^2) / sum(2*(x - x_lower)/(x_upper - x_lower)^2)
    tmp_back[,3]  <- mix_ratio[M] * (2*(x_upper - x)/(x_upper - x_lower)^2) / sum(2*(x_upper - x)/(x_upper - x_lower)^2)

    tmp          <- cbind(sapply(1:(K), f_k), tmp_back)
    tmp[tmp < 0] <- 0
    den          <- apply(tmp, 1, sum)

    w_k <- matrix(NA, nrow=M, ncol=N)

    for(j in 1:K) {
      w_k[j,] <- y * mix_ratio[j] * truncated_pv(x, mu[j], sigma[j], eta[j]) / den
    }

    w_k[K+1,]   <- y * mix_ratio[K+1] * ((1/(x_N - x[1])) * ramp_ind_uni(x, x[1])) / sum(((1/(x_N - x[1])) * ramp_ind_uni(x, x[1]))) / den
    w_k[K+2, ]  <- y * mix_ratio[K+2] * (2*(x - x_lower)/(x_upper - x_lower)^2) / sum(2*(x - x_lower)/(x_upper - x_lower)^2) / den
    w_k[M, ]    <- y * mix_ratio[M] * (2*(x_upper - x)/(x_upper - x_lower)^2) / sum(2*(x_upper - x)/(x_upper - x_lower)^2) /den
    w_k[w_k < 0] <- 0

    n_k                    <- apply(w_k, 1, sum)
    n_k[which(is.na(n_k))] <- 0

    #store for each parameter
    mu_cal      <- c()
    sigma_cal   <- c()
    eta_cal     <- c()

    #Updating mix_ratio
    mix_ratio  <- n_k/sum(y)
    mix_ratio[which(mix_ratio < 1e-10)] <- 1e-10

    #Updating mu
    for(k in 1:K) {
      opt    <- optimize(Q_fun, interval = c(min(x), max(x)), tol = 1e-10, x = x, sigma = sigma[k], eta = eta[k], w_k = w_k[k,], mix_ratio = mix_ratio[k], maximum = TRUE)
      mu_cal <- c(mu_cal, opt$maximum)
    }
    mu       <- mu_cal

    #Updating sigma
    for(k in 1:K) {
      opt       <- optimize(Q_fun, interval = c(1e-3, 1000), tol = 1e-10, x = x, mu = mu[k], eta = eta[k], w_k = w_k[k,], mix_ratio = mix_ratio[k], maximum = TRUE)
      sigma_cal <- c(sigma_cal, opt$maximum)
    }
    sigma       <- sigma_cal

    #Updating eta
    for(k in 1:K) {
      opt      <- optimize(Q_fun, interval = c(0, 1), tol = 1e-10, x = x, mu = mu[k], sigma = sigma[k], w_k = w_k[k,], mix_ratio = mix_ratio[k], maximum = TRUE)
      eta_cal <- c(eta_cal, opt$maximum)
    }
    eta      <- eta_cal


    LL_1[i+1]    <- LL(x, y, mu, sigma, eta, mix_ratio, x_lower, x_upper, x_N)
    mu_1         <- rbind(mu_1, mu)
    sigma_1      <- rbind(sigma_1, sigma)
    eta_1        <- rbind(eta_1, eta)
    mix_ratio_1  <- rbind(mix_ratio_1, mix_ratio)

    #Convergence check
    print(LL_1[i+1] - LL_1[i])
    if(abs(LL_1[i+1] - LL_1[i]) < conv.cri){ messe = "Convereged"; break }
  }

  #Store for result
  end_cal        <- Sys.time()
  cal_time       <- difftime(end_cal, start_cal, units = "sec")

  list(mu = mu, sigma = sigma, eta = eta, mix_ratio = mix_ratio,it = i, LL = LL_1,
       MU = mu_1, SIGMA = sigma_1, ETA = eta_1, MIX_RATIO = mix_ratio_1,
       convergence = messe, W_K = w_k, cal_time = cal_time)
}







