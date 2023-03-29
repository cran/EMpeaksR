\name{show_pvmm_lback_curve}

\alias{show_pvmm_lback_curve}

\title{Visualization of the result of spect_em_pvmm_lback}

\usage{
show_pvmm_lback_curve(spect_em_pvmm_lback_res,
                      x, y,
                      mix_ratio_init,
                      mu_init,
                      sigma_init,
                      eta_init,
                      x_lower,
                      x_upper)
}

\description{
Visualization of the result of spect_em_pvmm_lback().
}

\arguments{
  \item{spect_em_pvmm_lback_res}{data set obtained by spect_em_pvmm_lback()}
  \item{x}{measurement steps}
  \item{y}{intensity}
  \item{mu_init}{initial values of the mean of the components}
  \item{sigma_init}{initial values of the standard deviation of the components}
  \item{eta_init}{initial values of the mixing ratio of Gauss and Lorentz distribution}
  \item{mix_ratio_init}{initial values of the mixture ratio of the components}
  \item{x_lower}{lower limit of the measurement steps. Default is a minimum of x}
  \item{x_upper}{upper limit of the measurement steps. Default is a maximum of x}
}

\details{
Perform a visualization of fitting curve estimated by pseudo-Voigt mixture model with a linear background.
}

\value{
Show the fitting curve and variation of the parameters.

}

\references{
Matsumura, T., Nagamura, N., Akaho, S., Nagata, K., & Ando, Y. (2019). Spectrum adapted expectation-maximization algorithm for high-throughput peak shift analysis. Science and technology of advanced materials, 20(1), 733-745.

Matsumura, T., Nagamura, N., Akaho, S., Nagata, K., & Ando, Y. (2021). Spectrum adapted expectation-conditional maximization algorithm for extending high–throughput peak separation method in XPS analysis. Science and Technology of Advanced Materials: Methods, 1(1), 45-55.

Matsumura, T., Nagamura, N., Akaho, S., Nagata, K., & Ando, Y. (2023). High-throughput XPS spectrum modeling with autonomous background subtraction for 3 d 5/2 peak mapping of SnS. Science and Technology of Advanced Materials: Methods, 3(1), 2159753.
}

\examples{
#generating the synthetic spectral data based on three component Pseudo-Voigt mixture model.
x              <- seq(0, 100, by = 0.5)
K              <- 3
true_mu        <- c(35, 50, 65)
true_sigma     <- c(3, 3, 3)
true_mix_ratio <- c(0.5/3, 0.5/3, 0.5/3, 0.5)
true_eta       <- c(0.4, 0.6, 0.1)
degree         <- 4

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
                                      eta = true_eta[3])*10^degree +
       true_mix_ratio[4]*(c(500*x + 15000) / sum(500*x + 15000))*10^degree)

plot(y~x, main = "genrated synthetic spectral data")

#Peak fitting by EMpeaksR
#Initial values
mu_init        <- c(30, 40, 60)
sigma_init     <- c(4, 4, 4)
mix_ratio_init <- rep(1/(length(mu_init)+3), length(mu_init)+3)
eta_init       <- c(1, 1, 1)

#Coducting calculation
SP_ECM_PV_LBACK_res <- spect_em_pvmm_lback(x = x,
                                           y = y,
                                           mu = mu_init,
                                           sigma = sigma_init,
                                           eta = eta_init,
                                           mix_ratio = mix_ratio_init,
                                           x_lower = min(x),
                                           x_upper = max(x),
                                           conv.cri = 1e-2,
                                           maxit = 2000)

#Plot fitting curve and trace plot of parameters
show_pvmm_lback_curve(spect_em_pvmm_lback_res = SP_ECM_PV_LBACK_res,
                      x = x,
                      y = y,
                      mix_ratio_init = mix_ratio_init,
                      mu_init = mu_init,
                      sigma_init = sigma_init,
                      eta_init = eta_init,
                      x_lower = min(x),
                      x_upper = max(x))


#Showing the result of spect_em_pvmm_lback()
print(cbind(SP_ECM_PV_LBACK_res$mu, SP_ECM_PV_LBACK_res$sigma, SP_ECM_PV_LBACK_res$eta,
            SP_ECM_PV_LBACK_res$mix_ratio[1:K]))

print(cbind(true_mu, true_sigma, true_eta, true_mix_ratio[1:K]))

}







