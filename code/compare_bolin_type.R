library(tidyverse)
source('code/multi_matern_source.R')
library(ggplot2)
theme_set(theme_bw())
n_points <- 2^15
grid_info <- create_grid_info_1d(n_points, x_max = 80)

fft_1d <- function(grid_info, nu1 = .5, nu2 = .5, a1 = 1, a2 = 1, 
                   re, im, norm_type = 'A') {
  phase_factor = grid_info[['phase_factor']]
  delta_t = grid_info[['delta_t']]
  n_points = grid_info[['n_points']]
  x_vals = grid_info[['x_vals']]
  freq_points = grid_info[['freq_points']]
  x_max = grid_info[['x_max']]
  # https://stackoverflow.com/questions/24077913/discretized-continuous-fourier-transform-with-numpy
  tv <- (a1^2 + freq_points^2)^(-nu1/2 - .5/2) *
    (a2^2 + freq_points^2)^(-nu2/2 - .5/2) * 
    complex(real = re, imaginary = im*sign(freq_points))
  ff_res <- fftwtools::fftw_c2c(data = tv, inverse = 1)
  p <- length(ff_res)/2
  ff_res_adj <- c(ff_res[(p + 1):(2*p)], ff_res[1:p]) * phase_factor
  cbind(x_vals, 'val' = 
          as.double(Re(ff_res_adj) * norm_constant(nu1, nu2, a1, a2, norm_type = norm_type)) / x_max * n_points * 2.512596 )
}



full_function <- function(grid_info, nu, a, norm_type = 'A') {
  fft_1d(grid_info, nu1 = nu, nu2 = nu, a1 = 1, a2 = a, re = 1, im = 0, norm_type = norm_type)[,2]
}

a_vals <- c(.1, .4, 1, 2)
fun_vals <- sapply(a_vals, FUN = full_function, grid_info = grid_info,
                   nu = 1.5, norm_type = 'A')
df <- data.frame(x = grid_info$x_vals, 
                 value = as.double(fun_vals),
                 a = rep(a_vals, each = n_points))

ggplot(data = df %>%
         filter(abs(x) < 8), aes(x = x, y = value, color = factor(a), linetype = factor(a))) + 
  geom_line(size = .65) + 
  labs(x = 'Lags', y = 'Cross-covariance function value', color = expression(a[j]),
       linetype = expression(a[j]))
ggsave('images/bolin_version_range.png', height = 3, width = 6)

true_func <- function(x, a1 = 1, a2) {
  pi * (a1 * exp(-a2 * abs(x)) - a2 * exp(-a1 * abs(x))) / (a1^2 - a2^2) / a1 / a2 * 
    norm_constant(3/2, 3/2, a1, a2, norm_type = 'A')
}

fun_vals <- sapply(a_vals, FUN = true_func, x = grid_info$x_vals,
                   a1 = 1)
df <- data.frame(x = grid_info$x_vals, 
                 value = as.double(fun_vals),
                 a = rep(a_vals, each = n_points))

ggplot(data = df %>%
         filter(abs(x) < 8), aes(x = x, y = value, color = factor(a), linetype = factor(a))) + 
  geom_line() + 
#  scale_y_continuous(limits = c(-1,1)) + 
  labs(x = 'Lags', y = 'Cross-covariance function value', color = expression(a),
       linetype = expression(a))



full_function <- function(grid_info, nu, a, norm_type = 'A') {
  fft_1d(grid_info, nu1 = nu, nu2 = nu, a1 = 1, a2 = a, re = 0, im = 1, norm_type = norm_type)[,2]
}

a_vals <- c(.1, .4, 1, 2)
fun_vals <- sapply(a_vals, FUN = full_function, grid_info = grid_info,
                   nu = 1.5, norm_type = 'A')
df <- data.frame(x = grid_info$x_vals, 
                 value = as.double(fun_vals),
                 a = rep(a_vals, each = n_points))

ggplot(data = df %>%
         filter(abs(x) < 8), aes(x = x, y = value, color = factor(a), linetype = factor(a))) + 
  geom_line(size = .65) + 
  labs(x = 'Lags', y = 'Cross-covariance function value', color = expression(a[j]),
       linetype = expression(a[j]))
ggsave('images/bolin_version_range_im.png', height = 3, width = 6)
