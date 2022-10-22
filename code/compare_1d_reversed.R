library(tidyverse)
source('code/multi_matern_source.R')
theme_set(theme_bw())
grid_info <- create_grid_info_1d(2^14, 80)

re_z <- .4
im_z <- .5
a1 <- 1
a2 <- 1
nu2 <- .5
nu_vals <- c(.3, 1.5, 2.5, .8)
test25 <- fft_1d(grid_info, nu2 = nu2, nu1 = nu_vals[1], a1 = a1, a2 = a2, re = re_z, im = im_z)[,2]
test15 <- fft_1d(grid_info, nu2 = nu2, nu1 = nu_vals[2], a1 = a1, a2 = a2, re = re_z, im = im_z)[,2]
test250 <- fft_1d(grid_info, nu2 = nu2, nu1 = nu_vals[3], a1 = a1, a2 = a2, re = re_z, im = im_z)[,2]
test8 <- fft_1d(grid_info, nu2 = nu2, nu1 = nu_vals[4], a1 = a1, a2 = a2, re = re_z, im = im_z)[,2]
df <- data.frame(x = grid_info$x_vals, 
                 value = c(test25,  test15, test250, test8),
                 nu = rep(nu_vals, each = length(grid_info$x_vals)))
ggplot(data = filter(df, abs(x) < 5), aes(x = x, y = value, color = factor(nu), linetype = factor(nu))) + 
  geom_line() + 
  labs(x = 'Lags', y = 'Cross-covariance function value',
       color = expression(nu[j]),
       linetype = expression(nu[j]))

## opposite
test25_opp <- fft_1d(grid_info, nu1 = nu2, nu2 = nu_vals[1], a2 = a1, a1 = a2, re = re_z, im = -im_z)[,2]
test15_opp <- fft_1d(grid_info, nu1 = nu2, nu2 = nu_vals[2], a2 = a1, a1 = a2, re = re_z, im = -im_z)[,2]
test250_opp <- fft_1d(grid_info, nu1 = nu2, nu2 = nu_vals[3], a2 = a1, a1 = a2, re = re_z, im = -im_z)[,2]
test8_opp <- fft_1d(grid_info, nu1 = nu2, nu2 = nu_vals[4], a2 = a1, a1 = a2, re = re_z, im = -im_z)[,2]
df_opp <- data.frame(x = grid_info$x_vals,
                 value = c(test25_opp, test15_opp, test250_opp, test8_opp),
                 nu = rep(nu_vals, each = length(grid_info$x_vals)))

ggplot(data = filter(df_opp, abs(x) < 5), aes(x = x, y = value, color = factor(nu), linetype = factor(nu))) + 
  geom_line() + 
  labs(x = 'Lags', y = 'Cross-covariance function value',
       color = expression(nu[j]),
       linetype = expression(nu[j]))
ggplot(data = rbind(cbind(df, type = 'original'),
                    cbind(df_opp, type = 'opposite')) %>%
         filter(abs(x) < 5), aes(x = x, y = value, color = factor(nu), linetype = factor(nu))) + 
  geom_line() + 
  facet_wrap(~type, ncol = 1) + 
  labs(x = 'Lags', y = 'Cross-covariance function value',
       color = expression(nu[j]),
       linetype = expression(nu[j]))

ggplot(data = rbind(cbind(df, type = 'original'),
                    cbind(dplyr::mutate(df_opp, x = -x), type = 'opposite')) %>%
         filter(abs(x) < 5), 
       aes(x = x, y = value, color = factor(nu), linetype = factor(nu))) + 
  geom_line() + 
  facet_wrap(~type, ncol = 1) + 
  labs(x = 'Lags', y = 'Cross-covariance function value',
       color = expression(nu[j]),
       linetype = expression(nu[j]))
