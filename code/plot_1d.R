library(tidyverse)
source('code/multi_matern_source.R')
theme_set(theme_bw())
plot_cc <- function(plot_seq =  seq(-5, 5, by = .05), cons, nu) {
  params <- par()$mfrow
  cov_val_lag <- sapply(1:length(plot_seq), function(x) cross_cov(0, plot_seq[x], nu = nu, z_ij = cons))
  plot(plot_seq, cov_val_lag, type = 'l', xlab = 'lag', main = 'Cross Covariance')
}

plot_function <- function(h,nu, a = 1) {
  if(h == 0) {
    return(0)
  }
  sign(h) * (abs(h)/a)^nu*
    (besselI(a*abs(h), nu = nu) - struve(a*abs(h), -nu))
}

norm_constant <- function(nu_1, nu_2, a_1 = 1, a_2 = 1, norm_type = 'A') {
  if (norm_type == 'A') {
    (a_1 + a_2)^(nu_1 + nu_2)  /2/pi/gamma(nu_1 + nu_2) * gamma(nu_1 + 1/2) * gamma(nu_2 + 1/2)
  } else if (norm_type == 'B') {
    (2*a_1)^(nu_1) * (2*a_2)^(nu_2) /2/pi/sqrt(gamma(2*nu_1)*gamma(2 * nu_2)) * 
      gamma(nu_1 + 1/2) * gamma(nu_2 + 1/2)
  } else if (norm_type == 'C') {
    a_1^(nu_1 + 1/2) * a_2^(nu_2 + 1/2) /2/pi
  } else if (norm_type == 'D') {
    (a_1)^(nu_1) * (a_2)^(nu_2)*
      sqrt(gamma(nu_1 + 1/2)) * sqrt(gamma(nu_2 + 1/2))/pi^(1/2)/sqrt(gamma(nu_1)*gamma(nu_2))
  }
}


full_function <- function(plot_seq, nu, a, norm_type = 'A') {
  -sapply(1:length(plot_seq), function(x) plot_function(h = plot_seq[x], nu= nu, a = a))* 
           2^(-nu) * pi^(3/2)/cos(nu * pi)/gamma(nu + .5)*
   norm_constant(nu_1 = nu, nu_2 = nu, a_1 = a, a_2 = a, norm_type = norm_type) 
}
nu_vals <- c(.25, .75, 1.75, 2.25, 2.75)
plot_seq <- seq(-8, 8, by = .05)
fun_vals <- sapply(nu_vals, full_function, plot_seq = plot_seq,
                   a = 1, norm_type = 'D')
df <- data.frame(plot_seq, 
                 value = as.double(fun_vals),
                 nu = rep(nu_vals, each = length(plot_seq)))

ggplot(data = df, aes(x = plot_seq, y = value, color = factor(nu), linetype = factor(nu))) + 
  geom_line() + 
  labs(x = 'Lags', y = 'Cross-covariance function value', color = expression(nu),
       linetype = expression(nu))
ggsave('images/example_fun.png', height = 3, width = 6)

a_vals <- c(.2, .4, .7, 1, 1.4)
plot_seq <- seq(-8, 8, by = .05)
fun_vals <- sapply(a_vals, full_function, plot_seq = plot_seq,
                  nu = .75, norm_type = 'D')
df <- data.frame(plot_seq, 
                 value = as.double(fun_vals),
                 a = rep(a_vals, each = length(plot_seq)))

ggplot(data = df, aes(x = plot_seq, y = value, color = factor(a), linetype = factor(a))) + 
  geom_line() + 
  labs(x = 'Lags', y = 'Cross-covariance function value', color = expression(a),
       linetype = expression(a))
ggsave('images/example_fun_range.png', height = 3, width = 6)

width <- 6
plot_seq <- seq(-5, 5, by = .01)
nu_vals <- c(.2, .5, .8, 1.1, 1.8, 2.8)
full_function <- function(plot_seq, nu1, nu2, a1, a2, norm_type = 'A') {
  sapply(1:length(plot_seq), function(x) whitt_version(h = plot_seq[x], nu1 = nu1, nu2 = nu2,
                                                       c2 = 1,c11 = 1, c12 = 1, a1 = a1, a2 = a2)*
           norm_constant(nu_1 = nu1, nu_2 = nu2, a_1 = a1, a_2 = a2, norm_type = norm_type))[2,]
}

whitt_vals <- sapply(nu_vals, full_function, plot_seq = plot_seq, nu1 = .5, a1 = 1, a2 = 1, norm_type = 'A')


df <- data.frame(plot_seq, 
                 value = as.double(whitt_vals),
                 nu = rep(nu_vals, each = length(plot_seq)))

library(ggplot2)
theme_set(theme_bw())
ggplot(data = df, aes(x = plot_seq, y = value, color = factor(nu), linetype = factor(nu))) + 
  geom_line() + 
  labs(x = 'Lags', y = 'Cross-covariance function', color = expression(nu[k]),
       linetype = expression(nu[k]))#+
  #theme(legend.position = 'bottom')
ggsave('images/example_fun_whittaker_norm_A.png', height = 3, width = width)

whitt_vals <- sapply(nu_vals, full_function, plot_seq = plot_seq, 
                     nu1 = .5, a1 = 1, a2 = 1, norm_type = 'D')
df <- data.frame(plot_seq, 
                 value = as.double(whitt_vals),
                 nu = rep(nu_vals, each = length(plot_seq)))
ggplot(data = df, aes(x = plot_seq, y = value, color = factor(nu), linetype = factor(nu))) + 
  geom_line() + 
  labs(x = 'Lags', y = 'Cross-covariance function', color = expression(nu[k]),
       linetype = expression(nu[k]))#+
  #theme(legend.position = 'bottom')
ggsave('images/example_fun_whittaker_norm_B.png', height = 3, width = width)

avals <- c(.2, .5, 1, 2)
whitt_vals <- sapply(avals, full_function, plot_seq = plot_seq, 
                     nu1 = 1.2, nu2 = 1.2, a1 = 1, norm_type = 'A')
df <- data.frame(plot_seq, 
                 value = as.double(whitt_vals),
                 a = rep(avals, each = length(plot_seq)))
ggplot(data = df, aes(x = plot_seq, y = value, color = factor(a), linetype = factor(a))) + 
  geom_line() + 
  labs(x = 'Lags', y = 'Cross-covariance function', color = expression(a[k]),
       linetype = expression(a[k]))#+
  #theme(legend.position = 'bottom')
ggsave('images/example_fun_whittaker_a_norm_A.png', height = 3, width = width)

whitt_vals <- sapply(avals, full_function, plot_seq = plot_seq, 
                     nu1 = 1.2, nu2 = 1.2, a1 = 1, norm_type = 'D')
df <- data.frame(plot_seq, 
                 value = as.double(whitt_vals),
                 a = rep(avals, each = length(plot_seq)))
ggplot(data = df, aes(x = plot_seq, y = value, color = factor(a), linetype = factor(a))) + 
  geom_line() + 
  labs(x = 'Lags', y = 'Cross-covariance function', color = expression(a[k]),
       linetype = expression(a[k]))#+
 # theme(legend.position = 'bottom')
ggsave('images/example_fun_whittaker_a_norm_B.png', height = 3, width = width)


full_function <- function(plot_seq, nu, a, realp, imp, norm_type = 'A') {
  -imp * sapply(1:length(plot_seq), function(x) plot_function(h = plot_seq[x], nu= nu, a = a))* 
    2^(-nu) * pi^(3/2)/cos(nu * pi)/gamma(nu + .5)*
    norm_constant(nu_1 = nu, nu_2 = nu, a_1 = a, a_2 = a, norm_type = norm_type)  +
    realp *  sapply(1:length(plot_seq), function(x) whitt_version(h = plot_seq[x], nu1 = nu, nu2 = nu,
                                                                c2 = 1,c11 = 1, c12 = 1, a1 = a, a2 = a)*
                    norm_constant(nu_1 = nu, nu_2 = nu, a_1 = a, a_2 = a, norm_type = norm_type))[2,]
}

cov_val_lag <- full_function(plot_seq = plot_seq, nu = .75, a = 1, realp = 1, imp = 0, norm_type = 'D')
cov_val_lag2 <- full_function(plot_seq = plot_seq, nu = .75, a = 1, realp = 0, imp = 1, norm_type = 'D')
cov_val_lag3 <- full_function(plot_seq = plot_seq, nu = .75, a = 1, realp = sqrt(2)/sqrt(3), imp = 1/sqrt(3), norm_type = 'D')
cov_val_lag4 <- full_function(plot_seq = plot_seq, nu = .75, a = 1, realp = 1/sqrt(3), imp = sqrt(2)/sqrt(3), norm_type = 'D')


df <- data.frame(plot_seq, 
                 value = c(cov_val_lag, cov_val_lag2, cov_val_lag3, cov_val_lag4),
                 nu = factor(rep(c('1', 'i', '(2 + 1i)/sqrt(3)', '(1 + 2i)/sqrt(3)'), 
                          each = length(plot_seq)),
                          levels = c('1', '(2 + 1i)/sqrt(3)', '(1 + 2i)/sqrt(3)', 'i')))

ggplot(data = df, aes(x = plot_seq, y = value, color = factor(nu), linetype = factor(nu))) + 
  geom_line() + 
  labs(x = 'Lags', y = 'Cross-covariance function value', color = expression(Sigma[jk]),
       linetype = expression(Sigma[jk]))+
  scale_color_discrete(labels = expression(1, frac(sqrt(2) + i,sqrt(3)),frac(1 + sqrt(2)*i,sqrt(3)),i))+
  scale_linetype_discrete(labels = expression(1, frac(sqrt(2) + i, sqrt(3)),frac(1 + sqrt(2)*i,sqrt(3)), i))
ggsave('images/example_combination_function.png', height = 3, width = width)

spec_dens <- function(r, h, nu1, nu2, a1 = 1, a2 = 1, re_z, im_z) {
  val <- exp(complex(imaginary = -h * r)) *
    complex(real = a1, imaginary = r)^(-nu1 - 1/2) *
    complex(real = a2, imaginary = -r)^(-nu2 - 1/2) * 
    complex(real = re_z, imaginary = -sign(r) * im_z)
  Re(val)
}
re_z <- .4
im_z <- .5
a1 <- 1
a2 <- 1
nu1 <- .5
limits <- 300
subdivisions <- 4000
nu_vals <- c(.3, 1.5, 2.5, .8)
full_function <- function(plot_seq, nu1, nu2, a1, a2, re_z, im_z, limits, subdivisions, norm_type = 'D'){
  sapply(plot_seq, function(z) {
    integrate(spec_dens, h = z, nu1 = nu1, nu2 = nu2, re_z = re_z, im_z = im_z,
              a1 = a1, a2 = a2,subdivisions = subdivisions,
              lower = -limits, upper = limits)$value*
      norm_constant(
        nu_1 = nu1, nu_2 = nu2, 
        a_1 = a1, a_2 = a2, norm_type = 'D')
  })
}
test25 <- full_function(plot_seq, nu1, nu2 = nu_vals[1], a1, a2, re_z, im_z, limits, subdivisions, norm_type = 'D')
test15 <- full_function(plot_seq, nu1, nu2 = nu_vals[2], a1, a2, re_z, im_z, limits, subdivisions, norm_type = 'D')
test250 <- full_function(plot_seq, nu1, nu2 = nu_vals[3], a1, a2, re_z, im_z, limits, subdivisions, norm_type = 'D')
test8 <- full_function(plot_seq, nu1, nu2 = nu_vals[4], a1, a2, re_z, im_z, limits, subdivisions, norm_type = 'D')

df <- data.frame(plot_seq = plot_seq, 
                 value = c(test25,  test15, test250, test8),
                 nu = rep(nu_vals, each = length(plot_seq)))

ggplot(data = df, aes(x = plot_seq, y = value, color = factor(nu), linetype = factor(nu))) + 
  geom_line() + 
  labs(x = 'Lags', y = 'Cross-covariance function value',
       color = expression(nu[k]),
       linetype = expression(nu[k]))
ggsave('images/comb_varying_nu.png', height = 3, width = width)

nu1 <- 1.25
nu_vals <- c(.25, .75,  1.25, 2.25, 3.25)
rm(re_z);rm(im_z)
test_vals <- sapply(nu_vals, function(x) {
  full_function(plot_seq, nu1, nu2 = x, a1, a2, re_z = 0, im_z = 1, 
                limits = 40, subdivisions = 1000, norm_type = 'D')
})

df <- data.frame(plot_seq, 
                 value = as.double(test_vals),
                 nu = rep(nu_vals, each = length(plot_seq)))

ggplot(data = df, aes(x = plot_seq, y = value, color = factor(nu), linetype = factor(nu))) + 
  geom_line() + 
  labs(x = 'Lags', y = 'Cross-covariance function value', color = expression(nu[k]),
       linetype = expression(nu[k]))
ggsave('images/im_varying_nu.png', height = 3, width = width)

a_vals <- c(.2, .5, 1, 1.5, 2)
test_vals <- sapply(a_vals, function(x) {
  full_function(plot_seq, nu1 = 1.25, nu2 = 1.25, a1 = 1, a2 = x, re_z = 0, im_z = 1, 
                limits = 50, subdivisions = 1000, norm_type = 'D')
})

df <- data.frame(plot_seq, 
                 value = as.double(test_vals),
                 nu = rep(a_vals, each = length(plot_seq)))

ggplot(data = df, aes(x = plot_seq, y = value, color = factor(nu), linetype = factor(nu))) + 
  geom_line() + 
  labs(x = 'Lags', y = 'Cross-covariance function value', color = expression(a[k]),
       linetype = expression(a[k]))
ggsave('images/im_varying_a.png', height = 3, width = width)

