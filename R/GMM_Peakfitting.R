
spect_em_gmm <- function(x, y, mu, sigma, mix_ratio, conv.cri, maxit) {
  #数値の格納先の定義
  start_cal    <- Sys.time()
  messe        <- "Not converged"
  N            <- length(x)
  LL_1         <- numeric(0)
  mix_ratio_1  <- numeric(0)
  sigma_1      <- numeric(0)
  mu_1         <- numeric(0)
  n_k          <- numeric(0)
  K            <- length(mu)

  #重み付き対数尤度関数を計算する関数を定義
  f_k   <- function(i) {
    mix_ratio[i] * dnorm(x, mu[i], sigma[i])
  }

  LL    <- function(x, y, mu, sigma, mix_ratio) {
    pL  <- sapply(1:K,f_k)
    sum(y * log(apply(pL, 1, sum)))
  }

  LL_1[1]        <- LL(x, y, mu, sigma, mix_ratio)
  mu_1           <- rbind(mu_1, mu)
  sigma_1        <- rbind(sigma_1, sigma)
  mix_ratio_1    <- rbind(mix_ratio_1, mix_ratio)

  #反復計算の開始
  for(i in 1:maxit) {
    tmp <- sapply(1:K, f_k)
    den <- apply(tmp, 1, sum)
    w_k <- matrix(NA, nrow=K, ncol=N)

    #重み付きの負担率の計算
    for(j in 1:K) {
      w_k[j,] <- y * mix_ratio[j] * dnorm(x, mu[j], sigma[j]) / den
    }
    n_k <- apply(w_k, 1, sum)

    #負担率を使ってパラメータを更新
    for(j in 1:K) {
      mu[j]    <- sum((w_k[j,] * x)) / n_k[j]
      sigma[j] <- sqrt(sum(w_k[j,] * (x-mu[j])^2) / n_k[j])
    }
    mix_ratio        <- n_k/sum(y)

    #各値を更新
    LL_1[i+1]        <- LL(x, y, mu, sigma, mix_ratio)
    mu_1             <- rbind(mu_1, mu)
    sigma_1          <- rbind(sigma_1, sigma)
    mix_ratio_1      <- rbind(mix_ratio_1, mix_ratio)

    if(abs(LL_1[i+1] - LL_1[i]) < conv.cri) {
      messe = "Converged"; break
    }
    print(LL_1[i+1] - LL_1[i])
  }

  end_cal            <- Sys.time()
  cal_time           <- difftime(end_cal, start_cal, units = "sec")

  #出力する数値を格納
  list(mu = mu, sigma = sigma, mix_ratio = mix_ratio, it = i, LL = LL_1, MU = mu_1, SIGMA = sigma_1, MIX_RATIO = mix_ratio_1,
       convergence = messe, W_K = w_k, cal_time = cal_time)
}




