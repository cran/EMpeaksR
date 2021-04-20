
# EMpeaksR

<!-- badges: start -->
<!-- badges: end -->

The goal of EMpeaks is to efficiently perform the peak fitting on a large number of spectral data. 
The EMpeaks can support the investigation of peak fitting with two advantages: (1) a large amount of data can be processed at high speed; and (2) stable and automatic calculation can be easily performed.

## Installation

You can install the released version of EMpeaksR from [CRAN](https://CRAN.R-project.org) with:

``` r
install.packages("EMpeaksR")
```

## Example

This is a basic example which shows you how to solve a common problem:

``` r
library(EMpeaksR)
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

```

