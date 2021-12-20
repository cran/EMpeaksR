\name{show_dsgmm_curve}

\alias{show_dsgmm_curve}

\title{Visualization of the result of spect_em_dsgmm}

\usage{
show_dsgmm_curve(spect_em_dsgmm_res,
                 x,
                 y,
                 mix_ratio_init,
                 mu_init,
                 sigma_init,
                 alpha_init,
                 eta_init)
}

\description{
Visualization of the result of spect_em_dsgmm().
}

\arguments{
  \item{spect_em_dsgmm_res}{data set obtained by spect_em_dsgmm()}
  \item{x}{measurement steps}
  \item{y}{intensity}
  \item{mix_ratio_init}{initial values of the mixture ratio of the components}
  \item{mu_init}{initial values of the mean of the components}
  \item{sigma_init}{initial values of the standard deviation of the components}
  \item{alpha_init}{initial values of the asymmetric parameter of the components}
  \item{eta_init}{initial values of the mixing ratio of Gauss and Lorentz distribution}
}

\details{
Perform a visualization of fitting curve estimated by Doniach-Sunjic-Gauss mixture model.
}

\value{
Show the fitting curve and variation of the parameters.

}
\references{
Matsumura, T., Nagamura, N., Akaho, S., Nagata, K., & Ando, Y. (2019). Spectrum adapted expectation-maximization algorithm for high-throughput peak shift analysis. Science and technology of advanced materials, 20(1), 733-745.

Matsumura, T., Nagamura, N., Akaho, S., Nagata, K., & Ando, Y. (2021). Spectrum adapted expectation-conditional maximization algorithm for extending high–throughput peak separation method in XPS analysis. Science and Technology of Advanced Materials: Methods, 1(1), 45-55.
}

\examples{
#generating the synthetic spectral data based on three component Doniach-Sunjic-Gauss mixture model.
x               <- seq(0, 100, by = 0.5)
true_mu         <- c(20, 50, 80)
true_sigma      <- c(3, 3, 3)
true_alpha      <- c(0.1, 0.3, 0.1)
true_eta        <- c(0.4, 0.6, 0.1)
true_mix_ratio  <- rep(1/3, 3)
degree          <- 4

#trancated Doniach-Sunjic-Gauss
truncated_dsg <- function(x, mu, sigma, alpha, eta) {
                 ((eta*(((gamma(1-alpha)) /
                 ((x-mu)^2+(sqrt(2*log(2))*sigma)^2)^((1-alpha)/2)) *
                 cos((pi*alpha/2)+(1-alpha)*atan((x-mu) /
                 (sqrt(2*log(2))*sigma))))) + (1-eta)*dnorm(x, mu, sigma)) /
                 sum( ((eta*(((gamma(1-alpha)) /
                 ((x-mu)^2+(sqrt(2*log(2))*sigma)^2)^((1-alpha)/2)) *
                 cos((pi*alpha/2)+(1-alpha)*atan((x-mu) /
                 (sqrt(2*log(2))*sigma))))) + (1-eta)*dnorm(x, mu, sigma)))
}

y <- c(true_mix_ratio[1]*truncated_dsg(x = x,
                                       mu = true_mu[1],
                                       sigma = true_sigma[1],
                                       alpha = true_alpha[1],
                                       eta = true_eta[1])*10^degree +
       true_mix_ratio[2]*truncated_dsg(x = x,
                                       mu = true_mu[2],
                                       sigma = true_sigma[2],
                                       alpha = true_alpha[2],
                                       eta = true_eta[2])*10^degree +
       true_mix_ratio[3]*truncated_dsg(x = x,
                                       mu = true_mu[3],
                                       sigma = true_sigma[3],
                                       alpha = true_alpha[3],
                                       eta = true_eta[3])*10^degree)

plot(y~x, main = "genrated synthetic spectral data")

#Peak fitting by EMpeaksR
#Initial values
K <- 3
mix_ratio_init <- c(0.2, 0.4, 0.4)
mu_init        <- c(20, 40, 70)
sigma_init     <- c(4, 3, 2)
alpha_init     <- c(0.3, 0.2, 0.4)
eta_init       <- c(0.5, 0.4, 0.3)

#Coducting calculation
SP_ECM_DSG_res <- spect_em_dsgmm(x = x,
                                 y = y,
                                 mu = mu_init,
                                 sigma = sigma_init,
                                 alpha = alpha_init,
                                 eta = eta_init,
                                 mix_ratio = mix_ratio_init,
                                 conv.cri = 1e-6,
                                 maxit = 100)

#Plot fitting curve and trace plot of parameters
show_dsgmm_curve(SP_ECM_DSG_res,
                 x,
                 y,
                 mix_ratio_init,
                 mu_init,
                 sigma_init,
                 alpha_init,
                 eta_init)

#Showing the result of spect_em_dsgmm()
print(cbind(c(mu_init),
            c(sigma_init),
            c(alpha_init),
            c(eta_init),
            c(mix_ratio_init)))

print(cbind(SP_ECM_DSG_res$mu,
            SP_ECM_DSG_res$sigma,
            SP_ECM_DSG_res$alpha,
            SP_ECM_DSG_res$eta,
            SP_ECM_DSG_res$mix_ratio))

print(cbind(true_mu,
            true_sigma,
            true_alpha,
            true_eta,
            true_mix_ratio))
}









