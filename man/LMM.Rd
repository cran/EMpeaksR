\name{spect_em_lmm}

\alias{spect_em_lmm}

\title{Spectrum adapted ECM algorithm by LMM}

\usage{
spect_em_lmm(x, y, mu, gam, mix_ratio, conv.cri, maxit)
}

\description{
Perform a peak fitting based on the spectrum adapted ECM algorithm by Lorentz mixture model.
}

\arguments{
  \item{x}{measurement steps}
  \item{y}{intensity}
  \item{mu}{mean of the components}
  \item{gam}{scale parameter of the components}
  \item{mix_ratio}{mixture ratio of the components}
  \item{conv.cri}{criterion of the convergence}
  \item{maxit}{maximum number of the iteration}
}

\details{
Peak fitting is conducted by spectrum adapted ECM algorithm.
}

\value{
  \item{mu}{estimated mean of the components}
  \item{gam}{estimated scale parameter of the components}
  \item{mix_ratio}{estimated mixture ratio of the components}
  \item{it}{number of the iteration to reach the convergence}
  \item{LL}{variation of the weighted log likelihood values}
  \item{MU}{variation of mu}
  \item{GAM}{variation of gam}
  \item{MIX_RATIO}{variation of mix_ratio}
  \item{W_K}{decomposed component of the spectral data }
  \item{convergence}{message for the convergence in the calculation}
  \item{cal_time}{calculation time to complete the peak fitting. Unit is seconds}
}

\references{
Matsumura, T., Nagamura, N., Akaho, S., Nagata, K., & Ando, Y. (2019). Spectrum adapted expectation-maximization algorithm for high-throughput peak shift analysis. Science and technology of advanced materials, 20(1), 733-745.

Matsumura, T., Nagamura, N., Akaho, S., Nagata, K., & Ando, Y. (2021). Spectrum adapted expectation-conditional maximization algorithm for extending high–throughput peak separation method in XPS analysis. Science and Technology of Advanced Materials: Methods, 1(1), 45-55.
}

\examples{
#generating the synthetic spectral data based on three component Lorentz mixture model.
x               <- seq(0, 100, by = 0.5)
true_mu         <- c(35, 50, 65)
true_gam        <- c(3, 3, 3)
true_mix_ratio  <- rep(1/3, 3)
degree          <- 4

#Normalized Lorentz distribution
dCauchy <- function(x, mu, gam) {
    (dcauchy(x, mu, gam)) / sum(dcauchy(x, mu, gam))
  }

y <- c(true_mix_ratio[1] * dCauchy(x = x, mu = true_mu[1], gam = true_gam[1])*10^degree +
       true_mix_ratio[2] * dCauchy(x = x, mu = true_mu[2], gam = true_gam[2])*10^degree +
       true_mix_ratio[3] * dCauchy(x = x, mu = true_mu[3], gam = true_gam[3])*10^degree)

plot(y~x, main = "genrated synthetic spectral data")

#Peak fitting by EMpeaksR
#Initial values
K <- 3

mix_ratio_init  <- c(0.2, 0.4, 0.4)
mu_init         <- c(20, 40, 70)
gam_init        <- c(2, 5, 4)

#Coducting calculation
SP_ECM_L_res <- spect_em_lmm(x, y, mu = mu_init, gam = gam_init, mix_ratio = mix_ratio_init,
                             conv.cri = 1e-6, maxit = 100000)

#Plot fitting curve and trace plot of parameters
show_lmm_curve(SP_ECM_L_res, x, y, mix_ratio_init, mu_init, gam_init)

#Showing the result of spect_em_lmm()
print(cbind(c(mu_init), c(gam_init), c(mix_ratio_init)))
print(cbind(SP_ECM_L_res$mu, SP_ECM_L_res$gam, SP_ECM_L_res$mix_ratio))
print(cbind(true_mu, true_gam, true_mix_ratio))

}









