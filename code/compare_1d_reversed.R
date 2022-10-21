library(tidyverse)
#source('code/multi_matern_source.R')
theme_set(theme_bw())

spec_dens <- function(r, h, nu1, nu2, a1 = 1, a2 = 1, re_z, im_z) {
  val <- exp(complex(imaginary = h * r)) *
    complex(real = a1, imaginary = r)^(-nu1 - 1/2) *
    complex(real = a2, imaginary = -r)^(-nu2 - 1/2) * 
    complex(real = re_z, imaginary = sign(r) * im_z)
  Re(val)
}
spec_dens_opp <- function(r, h, nu1, nu2, a1 = 1, a2 = 1, re_z, im_z) {
  val <- exp(complex(imaginary = h * r)) *
    complex(real = a1, imaginary = -r)^(-nu1 - 1/2) *
    complex(real = a2, imaginary = r)^(-nu2 - 1/2) * 
    complex(real = re_z, imaginary = sign(r) * im_z)
  Re(val)
}
plot_seq <- seq(-5, 5, by = .01)

re_z <- .4
im_z <- .5
a1 <- 1
a2 <- 1
nu2 <- .5
limits <- 300
subdivisions <- 4000
nu_vals <- c(.3, 1.5, 2.5, .8)
full_function <- function(plot_seq, nu1, nu2, a1, a2, re_z, im_z, limits, subdivisions, norm_type = 'D'){
  sapply(plot_seq, function(z) {
    integrate(spec_dens, h = z, nu1 = nu1, nu2 = nu2, re_z = re_z, im_z = im_z,
              a1 = a1, a2 = a2,subdivisions = subdivisions,
              lower = -limits, upper = limits)$value 
  })
}
full_function_opp <- function(plot_seq, nu1, nu2, a1, a2, re_z, im_z, limits, subdivisions, norm_type = 'D'){
  sapply(plot_seq, function(z) {
    integrate(spec_dens_opp, h = z, nu1 = nu1, nu2 = nu2, re_z = re_z, im_z = im_z,
              a1 = a1, a2 = a2,subdivisions = subdivisions,
              lower = -limits, upper = limits)$value 
  })
}
test25 <- full_function(plot_seq, nu2, nu1 = nu_vals[1], a1, a2, re_z, im_z, limits, subdivisions, norm_type = 'D')
test15 <- full_function(plot_seq, nu2, nu1 = nu_vals[2], a1, a2, re_z, im_z, limits, subdivisions, norm_type = 'D')
test250 <- full_function(plot_seq, nu2, nu1 = nu_vals[3], a1, a2, re_z, im_z, limits, subdivisions, norm_type = 'D')
test8 <- full_function(plot_seq, nu2, nu1 = nu_vals[4], a1, a2, re_z, im_z, limits, subdivisions, norm_type = 'D')

df <- data.frame(plot_seq = plot_seq, 
                 value = c(test25,  test15, test250, test8),
                 nu = rep(nu_vals, each = length(plot_seq)))

ggplot(data = df, aes(x = plot_seq, y = value, color = factor(nu), linetype = factor(nu))) + 
  geom_line() + 
  labs(x = 'Lags', y = 'Cross-covariance function value',
       color = expression(nu[j]),
       linetype = expression(nu[j]))

## opposite
test25_opp <- full_function_opp(plot_seq, nu2, nu1 = nu_vals[1], a1, a2, re_z, -im_z, limits, subdivisions, norm_type = 'D')
test15_opp <- full_function_opp(plot_seq, nu2, nu1 = nu_vals[2], a1, a2, re_z, -im_z, limits, subdivisions, norm_type = 'D')
test250_opp <- full_function_opp(plot_seq, nu2, nu1 = nu_vals[3], a1, a2, re_z, -im_z, limits, subdivisions, norm_type = 'D')
test8_opp  <- full_function_opp(plot_seq, nu2, nu1 = nu_vals[4], a1, a2, re_z, -im_z, limits, subdivisions, norm_type = 'D')

df_opp <- data.frame(plot_seq = plot_seq, 
                 value = c(test25_opp, test15_opp, test250_opp, test8_opp),
                 nu = rep(nu_vals, each = length(plot_seq)))

ggplot(data = df_opp, aes(x = plot_seq, y = value, color = factor(nu), linetype = factor(nu))) + 
  geom_line() + 
  labs(x = 'Lags', y = 'Cross-covariance function value',
       color = expression(nu[j]),
       linetype = expression(nu[j]))
ggplot(data = rbind(cbind(df, type = 'original'),
                    cbind(df_opp, type = 'opposite')), aes(x = plot_seq, y = value, color = factor(nu), linetype = factor(nu))) + 
  geom_line() + 
  facet_wrap(~type, ncol = 1) + 
  labs(x = 'Lags', y = 'Cross-covariance function value',
       color = expression(nu[j]),
       linetype = expression(nu[j]))

ggplot(data = rbind(cbind(df, type = 'original'),
                    cbind(dplyr::mutate(df_opp, plot_seq = -plot_seq), type = 'opposite')), 
       aes(x = plot_seq, y = value, color = factor(nu), linetype = factor(nu))) + 
  geom_line() + 
  facet_wrap(~type, ncol = 1) + 
  labs(x = 'Lags', y = 'Cross-covariance function value',
       color = expression(nu[j]),
       linetype = expression(nu[j]))
