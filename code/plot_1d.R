library(tidyverse)
source('code/multi_matern_source.R')
library(ggplot2)
theme_set(theme_bw())
n_points <- 2^13
grid_info <- create_grid_info_1d(n_points, x_max = 60)

full_function <- function(grid_info, nu, a, norm_type = 'A') {
  fft_1d(grid_info, nu1 = nu, nu2 = nu, a1 = a, a2 = a, re = 0, im = 1, norm_type = norm_type)[,2]
}
nu_vals <- c(.25, .75, 1.75, 2.25, 2.75)
fun_vals <- sapply(nu_vals, full_function, grid_info = grid_info, a = 1)
df <- data.frame(x = grid_info[['x_vals']], 
                 value = as.double(fun_vals),
                 nu = rep(nu_vals, each = n_points))

ggplot(data = df %>%
         filter(abs(x) < 8), aes(x = x, y = value, color = factor(nu), linetype = factor(nu))) + 
  geom_line() + 
  labs(x = 'Lags', y = 'Cross-covariance function value', color = expression(nu),
       linetype = expression(nu))
ggsave('images/example_fun.png', height = 3, width = 6)

a_vals <- c(.2, .4, .7, 1, 1.4)
fun_vals <- sapply(a_vals, FUN = full_function, grid_info = grid_info,
                  nu = .75, norm_type = 'A')
df <- data.frame(x = grid_info$x_vals, 
                 value = as.double(fun_vals),
                 a = rep(a_vals, each = n_points))

ggplot(data = df %>%
         filter(abs(x) < 8), aes(x = x, y = value, color = factor(a), linetype = factor(a))) + 
  geom_line() + 
  labs(x = 'Lags', y = 'Cross-covariance function value', color = expression(a),
       linetype = expression(a))
ggsave('images/example_fun_range.png', height = 3, width = 6)

width <- 6
nu_vals <- c(.2, .5, .8, 1.1, 1.8, 2.8)
whitt_vals <- sapply(nu_vals, fft_1d, grid_info = grid_info, re = 1, im = 0,
                     nu2 = .5, a1 = 1, a2 = 1, norm_type = 'A')[-c(1:n_points),]


df <- data.frame(x = grid_info[['x_vals']], 
                 value = as.double(whitt_vals),
                 nu = rep(nu_vals, each = length(grid_info$x_vals)))

ggplot(data = filter(df, abs(x) < 5), aes(x = x, y = value, color = factor(nu), linetype = factor(nu))) + 
  geom_line() + 
  labs(x = 'Lags', y = 'Cross-covariance function', color = expression(nu[j]),
       linetype = expression(nu[j]))
ggsave('images/example_fun_whittaker_norm_A.png', height = 3, width = width)

whitt_vals <- sapply(nu_vals, fft_1d, grid_info = grid_info, re = 1, im = 0,
                     nu2 = .5, a1 = 1, a2 = 1, norm_type = 'B')[-c(1:n_points),]
df <- data.frame(x = grid_info[['x_vals']], 
                 value = as.double(whitt_vals),
                 nu = rep(nu_vals, each = n_points))
ggplot(data = filter(df, abs(x) < 5), aes(x = x, y = value, color = factor(nu), linetype = factor(nu))) + 
  geom_line() + 
  labs(x = 'Lags', y = 'Cross-covariance function', color = expression(nu[j]),
       linetype = expression(nu[j]))
ggsave('images/example_fun_whittaker_norm_B.png', height = 3, width = width)

avals <- c(.2, .5, 1, 2)
whitt_vals <- sapply(avals, fft_1d, grid_info = grid_info, re = 1, im = 0,
                     nu1 = 1.2, nu2 = 1.2, a2 = 1, norm_type = 'A')[-c(1:n_points),]
df <- data.frame(x = grid_info[['x_vals']], 
                 value = as.double(whitt_vals),
                 a = rep(avals, each = length( grid_info$x_vals)))
ggplot(data = filter(df, abs(x) < 5), aes(x = x, y = value, color = factor(a), linetype = factor(a))) + 
  geom_line() + 
  labs(x = 'Lags', y = 'Cross-covariance function', color = expression(a[j]),
       linetype = expression(a[j]))
ggsave('images/example_fun_whittaker_a_norm_A.png', height = 3, width = width)

whitt_vals <- sapply(avals, fft_1d, grid_info = grid_info, re = 1, im = 0,
                     nu1 = 1.2, nu2 = 1.2, a2 = 1, norm_type = 'B')[-c(1:n_points),]
df <- data.frame(x = grid_info$x_vals, 
                 value = as.double(whitt_vals),
                 a = rep(avals, each = length(grid_info$x_vals)))
ggplot(data = filter(df, abs(x) < 5), aes(x = x, y = value, color = factor(a), linetype = factor(a))) + 
  geom_line() + 
  labs(x = 'Lags', y = 'Cross-covariance function', color = expression(a[j]),
       linetype = expression(a[j]))
ggsave('images/example_fun_whittaker_a_norm_B.png', height = 3, width = width)

cov_val_lag <- fft_1d(grid_info, nu1 = .75, nu2 = .75, a1 = 1, a2 = 1, re = 1, im = 0, norm_type = 'A')[,2]
cov_val_lag2 <- fft_1d(grid_info, nu1 = .75, nu2 = .75, a1 = 1, a2 = 1, re = 0, im = 1, norm_type = 'A')[,2]
cov_val_lag3 <- fft_1d(grid_info, nu1 = .75, nu2 = .75, a1 = 1, a2 = 1, re = sqrt(2)/sqrt(3), im =  1/sqrt(3), norm_type = 'A')[,2]
cov_val_lag4 <- fft_1d(grid_info, nu1 = .75, nu2 = .75, a1 = 1, a2 = 1, re = 1/sqrt(3), im = sqrt(2)/sqrt(3), norm_type = 'A')[,2]

df <- data.frame(x = grid_info$x_vals, 
                 value = c(cov_val_lag, cov_val_lag2, cov_val_lag3, cov_val_lag4),
                 nu = factor(rep(c('1', 'i', '(2 + 1i)/sqrt(3)', '(1 + 2i)/sqrt(3)'), 
                          each = length(grid_info$x_vals)),
                          levels = c('1', '(2 + 1i)/sqrt(3)', '(1 + 2i)/sqrt(3)', 'i')))

ggplot(data = filter(df, abs(x) < 5), aes(x = x, y = value, color = factor(nu), linetype = factor(nu))) + 
  geom_line() + 
  labs(x = 'Lags', y = 'Cross-covariance function value', color = expression(Sigma[jk]),
       linetype = expression(Sigma[jk])) +
  scale_color_discrete(labels = expression(1, frac(sqrt(2) + i,sqrt(3)),frac(1 + sqrt(2)*i,sqrt(3)),i)) +
  scale_linetype_discrete(labels = expression(1, frac(sqrt(2) + i, sqrt(3)),frac(1 + sqrt(2)*i,sqrt(3)), i))
ggsave('images/example_combination_function.png', height = 3, width = width)

re_z <- .4
im_z <- .5
a1 <- 1
a2 <- 1
nu2 <- .5
nu_vals <- c(.3, 1.5, 2.5, .8)
cov_val_lag <- fft_1d(grid_info, nu1 = nu_vals[1], nu2 = nu2, a1 = a1, a2 = a2, re = re_z, im = im_z, norm_type = 'A')[,2]
cov_val_lag2 <- fft_1d(grid_info, nu1 = nu_vals[2], nu2 = nu2, a1 = a1, a2 = a2, re = re_z, im = im_z, norm_type = 'A')[,2]
cov_val_lag3 <- fft_1d(grid_info, nu1 = nu_vals[3], nu2 = nu2, a1 = a1, a2 = a2, re = re_z, im = im_z, norm_type = 'A')[,2]
cov_val_lag4 <- fft_1d(grid_info, nu1 = nu_vals[4], nu2 = nu2, a1 = a1, a2 = a2, re = re_z, im = im_z, norm_type = 'A')[,2]

df <- data.frame(x = grid_info[['x_vals']], 
                 value = c(cov_val_lag,  cov_val_lag2, cov_val_lag3, cov_val_lag4),
                 nu = rep(nu_vals, each = length(grid_info[['x_vals']])))

ggplot(data = filter(df, abs(x) < 5), aes(x = x, y = value, color = factor(nu), linetype = factor(nu))) + 
  geom_line() + 
  labs(x = 'Lags', y = 'Cross-covariance function value',
       color = expression(nu[j]),
       linetype = expression(nu[j]))
ggsave('images/comb_varying_nu.png', height = 3, width = width)

nu2 <- 1.25
nu_vals <- c(.25, .75,  1.25, 2.25, 3.25)
test_vals <- sapply(nu_vals, function(x) {
  fft_1d(grid_info, nu1 = x, nu2 = nu2, a1 = a1, a2 = a2, re = 0, im = 1, norm_type = 'A')[,2]
})

df <- data.frame(x = grid_info[['x_vals']], 
                 value = as.double(test_vals),
                 nu = rep(nu_vals, each = length(grid_info[['x_vals']])))

ggplot(data = filter(df, abs(x) < 5), aes(x = x, y = value, color = factor(nu), linetype = factor(nu))) + 
  geom_line() + 
  labs(x = 'Lags', y = 'Cross-covariance function value', color = expression(nu[j]),
       linetype = expression(nu[j]))
ggsave('images/im_varying_nu.png', height = 3, width = width)

a_vals <- c(.2, .5, 1, 1.5, 2)
test_vals <- sapply(a_vals, function(x) {
  fft_1d(grid_info, nu1 = 1.25, nu2 = 1.25, a1 = x, a2 = 1, re = 0, im = 1, norm_type = 'A')[,2]
})

df <- data.frame(x = grid_info[['x_vals']], 
                 value = as.double(test_vals),
                 nu = rep(a_vals, each = length(grid_info[['x_vals']])))

ggplot(data = filter(df, abs(x) < 5), aes(x = x, y = value, color = factor(nu), linetype = factor(nu))) + 
  geom_line() + 
  labs(x = 'Lags', y = 'Cross-covariance function value', color = expression(a[j]),
       linetype = expression(a[j]))
ggsave('images/im_varying_a.png', height = 3, width = width)

