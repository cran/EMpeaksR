\name{spect_em_gmm}

\alias{spect_em_gmm}

\title{Spectrum adapted EM algorithm by GMM}

\usage{
spect_em_gmm(x, y, mu, sigma, mix_ratio, conv.cri, maxit)
}

\description{
Perform a peak fitting based on the spectrum adapted EM algorithm by Gaussian mixture model.
}

\arguments{
  \item{x}{measurement steps}
  \item{y}{intensity}
  \item{mu}{mean of the components}
  \item{sigma}{standard deviation of the components}
  \item{mix_ratio}{mixture ratio of the components}
  \item{conv.cri}{criterion of the convergence}
  \item{maxit}{maximum number of the iteration}
}

\details{
Peak fitting is conducted by spectrum adapted EM algorithm.
}

\value{
  \item{mu}{estimated mean of the components}
  \item{sigma}{estimated standard deviation of the components}
  \item{mix_ratio}{estimated mixture ratio of the components}
  \item{it}{number of the iteration to reach the convergence}
  \item{LL}{variation of the weighted log likelihood values}
  \item{MU}{variation of mu}
  \item{SIGMA}{variation of sigma}
  \item{MIX_RATIO}{variation of mix_ratio}
  \item{W_K}{decomposed component of the spectral data }
  \item{convergence}{message for the convergence in the calculation}
  \item{cal_time}{calculation time to complete the peak fitting. Unit is seconds}
}

\references{
Matsumura, T., Nagamura, N., Akaho, S., Nagata, K., & Ando, Y. (2019). Spectrum adapted expectation-maximization algorithm for high-throughput peak shift analysis. Science and technology of advanced materials, 20(1), 733-745.
}

\examples{
#generating the synthetic spectral data based on three component Gausian mixture model.
x               <- seq(0, 100, by = 0.5)
true_mu         <- c(35, 50, 65)
true_sigma      <- c(3, 3, 3)
true_mix_ratio  <- rep(1/3, 3)
degree          <- 4

y <- c(true_mix_ratio[1] * dnorm(x = x, mean = true_mu[1], sd = true_sigma[1])*10^degree +
       true_mix_ratio[2] * dnorm(x = x, mean = true_mu[2], sd = true_sigma[2])*10^degree +
       true_mix_ratio[3] * dnorm(x = x, mean = true_mu[3], sd = true_sigma[3])*10^degree)

plot(y~x, main = "genrated synthetic spectral data")

#Peak fitting by EMpeaksR
#Initial values
K <- 3

mix_ratio_init  <- c(0.2, 0.4, 0.4)
mu_init         <- c(20, 40, 70)
sigma_init      <- c(2, 5, 4)

#Coducting calculation
SP_EM_G_res <- spect_em_gmm(x, y, mu = mu_init, sigma = sigma_init, mix_ratio = mix_ratio_init,
                            conv.cri = 1e-2, maxit = 2000)

#Plot fitting curve and trace plot of parameters
show_gmm_curve(SP_EM_G_res, x, y, mix_ratio_init, mu_init, sigma_init)

#Showing the result of spect_em_gmm()
print(cbind(c(mu_init), c(sigma_init), c(mix_ratio_init)))
print(cbind(SP_EM_G_res$mu, SP_EM_G_res$sigma, SP_EM_G_res$mix_ratio))
print(cbind(true_mu, true_sigma, true_mix_ratio))

}









