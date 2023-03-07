\name{show_pvmm_curve}

\alias{show_pvmm_curve}

\title{Visualization of the result of spect_em_pvmm}

\usage{
show_pvmm_curve(spect_em_pvmm_res, x, y, mix_ratio_init, mu_init, sigma_init, eta_init)
}

\description{
Visualization of the result of spect_em_pvmm().
}

\arguments{
  \item{spect_em_pvmm_res}{data set obtained by spect_em_pvmm()}
  \item{x}{measurement steps}
  \item{y}{intensity}
  \item{mix_ratio_init}{initial values of the mixture ratio of the components}
  \item{mu_init}{initial values of the mean of the components}
  \item{sigma_init}{initial values of the standard deviation of the components}
  \item{eta_init}{initial values of the mixing ratio of Gauss and Lorentz distribution}
}

\details{
Perform a visualization of fitting curve estimated by Pseudo-Voigt mixture model.
}

\value{
Show the fitting curve and variation of the parameters.

}
\references{
Matsumura, T., Nagamura, N., Akaho, S., Nagata, K., & Ando, Y. (2019). Spectrum adapted expectation-maximization algorithm for high-throughput peak shift analysis. Science and technology of advanced materials, 20(1), 733-745.

Matsumura, T., Nagamura, N., Akaho, S., Nagata, K., & Ando, Y. (2021). Spectrum adapted expectation-conditional maximization algorithm for extending high–throughput peak separation method in XPS analysis. Science and Technology of Advanced Materials: Methods, 1(1), 45-55.
}

\examples{
#generating the synthetic spectral data based on three component Pseudo-Voigt mixture model.
x               <- seq(0, 100, by = 0.5)
true_mu         <- c(35, 50, 65)
true_sigma      <- c(3, 3, 3)
true_eta        <- c(0.3, 0.8, 0.5)
true_mix_ratio  <- rep(1/3, 3)
degree          <- 4

#Normalized Pseudo-Voigt distribution
  truncated_pv <- function(x, mu, sigma, eta) {
    (eta*dcauchy(x, mu, sqrt(2*log(2))*sigma) + (1-eta)*dnorm(x, mu, sigma)) /
      sum(eta*dcauchy(x, mu, sqrt(2*log(2))*sigma) + (1-eta)*dnorm(x, mu, sigma))
  }

y <- c(true_mix_ratio[1]*truncated_pv(x = x,
                                      mu = true_mu[1],
                                      sigma = true_sigma[1],
                                      eta = true_eta[1])*10^degree +
       true_mix_ratio[2]*truncated_pv(x = x,
                                      mu = true_mu[2],
                                      sigma = true_sigma[2],
                                      eta = true_eta[2])*10^degree +
       true_mix_ratio[3]*truncated_pv(x = x,
                                      mu = true_mu[3],
                                      sigma = true_sigma[3],
                                      eta = true_eta[3])*10^degree)

plot(y~x, main = "genrated synthetic spectral data")

#Peak fitting by EMpeaksR
#Initial values
K <- 3

mix_ratio_init  <- c(0.2, 0.4, 0.4)
mu_init         <- c(20, 40, 70)
sigma_init      <- c(2, 5, 4)
eta_init        <- c(0.5, 0.4, 0.3)

#Coducting calculation
SP_ECM_PV_res <- spect_em_pvmm(x = x,
                               y = y,
                               mu = mu_init,
                               sigma = sigma_init,
                               eta = eta_init,
                               mix_ratio = mix_ratio_init,
                               conv.cri = 1e-3,
                               maxit = 3000)

#Plot fitting curve and trace plot of parameters
show_pvmm_curve(SP_ECM_PV_res, x, y, mix_ratio_init, mu_init, sigma_init, eta_init)

#Showing the result of spect_em_pvmm()
print(cbind(c(mu_init), c(sigma_init), c(eta_init), c(mix_ratio_init)))
print(cbind(SP_ECM_PV_res$mu, SP_ECM_PV_res$sigma, SP_ECM_PV_res$eta, SP_ECM_PV_res$mix_ratio))
print(cbind(true_mu, true_sigma, true_eta, true_mix_ratio))

}







