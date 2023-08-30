
values <- seq(-5, 5, length.out = 400)
values <- seq(-20, 20, length.out = 4000)

cc_hilbert <- function(value, a) {
  a_abs_value <- a * abs(value)
  -sign(value)/ pi * (exp(a_abs_value) * expint::expint_E1(a_abs_value) + 
                       exp(-a_abs_value) * expint::expint_Ei(a_abs_value))
} 
a_value = 10
result <- sapply(values, cc_hilbert, a = a_value)
plot(values, result, type = 'l')
library(tidyverse)
source('code/multi_matern_source.R')
theme_set(theme_bw())
grid_info <- create_grid_info_1d(2^14, 80)
original_values <- fft_1d(grid_info, nu2 = .5, nu1 = .5, a1 = a_value, a2 = a_value, re = 0, im = 1)[,2]
plot(grid_info$x_vals[abs(grid_info$x_vals)<5], original_values[abs(grid_info$x_vals)<5], type = 'l')
lines(values,result, type = 'l', col = 2)
lines(values, -2/pi/a_value/values, type = 'l', col = 3) #asymptotic expression


cc_hilbert <- function(value, a1, a2) {
  if (value <= 0) {
    -sign(value)/ pi * (exp(a1 * abs(value)) * expint::expint_E1(a1 * abs(value)) + 
                         exp(-a2 * abs(value)) * expint::expint_Ei(a2 * abs(value))) *
      2 * sqrt(a1 *a2) /(a1 + a2)
  } else {
    -sign(value)/ pi * (exp(a2 * abs(value)) * expint::expint_E1(a2 * abs(value)) + 
                         exp(-a1 * abs(value)) * expint::expint_Ei(a1 * abs(value))) *
      2 * sqrt(a1 *a2) /(a1 + a2)
  }
} 
a1_value <- .3
a2_value <- 3
result <- sapply(values, cc_hilbert, a1 = a1_value, a2 = a2_value)
original_values <- fft_1d(grid_info, nu2 = .5, nu1 = .5, a1 = a1_value, a2 = a2_value, re = 0, im = 1)[,2]
plot(grid_info$x_vals[abs(grid_info$x_vals)<20], original_values[abs(grid_info$x_vals)<20], type = 'l')
lines(values,result, type = 'l', col = 2)
lines(values,-2*sqrt(a1_value*a2_value)/(a1_value + a2_value)/pi*(1/a1_value/values + 1/a2_value/values), type = 'l', col = 3) #asymptotic expression

plot(grid_info$x_vals[abs(grid_info$x_vals)<5], original_values[abs(grid_info$x_vals)<5], type = 'l')
lines(values,result, type = 'l', col = 2)


# \nu = 3/2
cc_hilbert <- function(value, a) {
  a_abs_value <- a * abs(value)
  R_value <- -sign(value)/ pi * (exp(a_abs_value) * expint::expint_E1(a_abs_value) + 
                     exp(-a_abs_value) * expint::expint_Ei(a_abs_value))
  (a_abs_value + 1) * R_value - value * a * 2 * exp(a_abs_value)/pi * expint::expint_Ei(-a_abs_value)
} 
a_value <- 4
result <- sapply(values, cc_hilbert, a = a_value)
original_values <- fft_1d(grid_info, nu2 = 1.5, nu1 = 1.5, a1 = a_value, a2 = a_value, re = 0, im = 1)[,2]
plot(grid_info$x_vals[abs(grid_info$x_vals)<5], original_values[abs(grid_info$x_vals)<5], type = 'l',
     ylim = c(-1,1))
lines(values,result, type = 'l', col = 2)
lines(values,-2 * gamma(2)/sqrt(pi)/gamma(3/2)/a_value/values, type = 'l', col = 3) #asymptotic expression
lines(values,-4/pi/a_value/values, type = 'l', col = 4) #asymptotic expression


