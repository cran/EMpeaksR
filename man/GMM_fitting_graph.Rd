\name{show_gmm_curve}

\alias{show_gmm_curve}

\title{Visualization of the result of spect_em_gmm}

\usage{
show_gmm_curve(spect_em_gmm_res, x, y, mix_ratio_init, mu_init, sigma_init)
}

\description{
Visualization of the result of spect_em_gmm().
}

\arguments{
  \item{spect_em_gmm_res}{data set obtained by spect_em_gmm()}
  \item{x}{measurement steps}
  \item{y}{intensity}
  \item{mix_ratio_init}{initial values of the mixture ratio of the components}
  \item{mu_init}{initial values of the mean of the components}
  \item{sigma_init}{initial values of the standard deviation of the components}
}

\details{
Perform a visualization of fitting curve estimated by Gaussian mixture model.
}

\value{
Show the fitting curve and variation of the parameters.

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
                            conv.cri = 1e-8, maxit = 100000)

#Plot fitting curve and trace plot of parameters
show_gmm_curve(SP_EM_G_res, x, y, mix_ratio_init, mu_init, sigma_init)

#Showing the result of spect_em_gmm()
print(cbind(c(mu_init), c(sigma_init), c(mix_ratio_init)))
print(cbind(SP_EM_G_res$mu, SP_EM_G_res$sigma, SP_EM_G_res$mix_ratio))
print(cbind(true_mu, true_sigma, true_mix_ratio))

}









