
source('code/multi_matern_source.R')
A <- complex(real = rnorm(1), 
             imaginary = rnorm(1))

plot_cc <- function(plot_seq =  seq(-5, 5, by = .05), cons, nu) {
  params <- par()$mfrow
  cov_val_lag <- sapply(1:length(plot_seq), function(x) cross_cov(0, plot_seq[x], nu = nu, z_ij = cons))
  plot(plot_seq, cov_val_lag, type = 'l', xlab = 'lag', main = 'Cross Covariance')
}

plot_function <- function(s,t,nu) {
  if(s-t == 0) {
    return(0)
  }
 sign(t-s) * (abs(t-s))^nu*
    (besselI(abs(t-s), nu = nu) - struve(abs(t-s), -nu))
}

norm_constant <- function(d, x, nu_1, nu_2, a_1 = 1, a_2 = 1) {
  a_1^(nu_1) * a_2^(nu_2)
}

plot_seq <- seq(-8, 8, by = .05)
nu25 <- sapply(1:length(plot_seq), function(x) plot_function(0, plot_seq[x], nu = .25) *
                 norm_constant(d = 1, x = plot_seq[x], nu_1 = .25, nu_2 = .25))
nu75 <- sapply(1:length(plot_seq), function(x) plot_function(0, plot_seq[x], nu = .75) *
                 norm_constant(d =1, x = plot_seq[x], nu_1 = .75, nu_2 = .75))
nu125 <- sapply(1:length(plot_seq), function(x) plot_function(0, plot_seq[x], nu = 1.25)*
                  norm_constant(d =1, x = plot_seq[x], nu_1 = 1.25, nu_2 = 1.25))
nu01 <- sapply(1:length(plot_seq), function(x) plot_function(0, plot_seq[x], nu = .1)*
                 norm_constant(d =1, x = plot_seq[x], nu_1 = .1, nu_2 = .1))
nu175 <- sapply(1:length(plot_seq), function(x) plot_function(0, plot_seq[x], nu = 1.75)*
                  norm_constant(d =1, x = plot_seq[x], nu_1 = 1.75, nu_2 = 1.75))
nu225 <- sapply(1:length(plot_seq), function(x) plot_function(0, plot_seq[x], nu = 2.25) *
                  norm_constant(d =1, x = plot_seq[x], nu_1 = 2.25, nu_2 = 2.25))
nu275 <- sapply(1:length(plot_seq), function(x) plot_function(0, plot_seq[x], nu = 2.75)*
                  norm_constant(d =1, x = plot_seq[x], nu_1 = 2.75, nu_2 = 2.75))

df <- data.frame(plot_seq, 
                 value = c(nu01, nu25, nu75, nu125, nu175, nu225, nu275),
                 nu = rep(c(.01, .25, .75, 1.25, 1.75, 2.25, 2.75), each = length(plot_seq)))

library(ggplot2)
theme_set(theme_bw())
ggplot(data = df, aes(x = plot_seq, y = value, color = factor(nu), linetype = factor(nu))) + 
  geom_line() + 
  labs(x = 'Lags', y = 'Cross-covariance function value', color = expression(nu),
       linetype = expression(nu))+
  theme(legend.position = 'bottom')
ggsave('images/example_fun.png', height = 4, width = 4)

plot_seq <- seq(-8, 8, by = .01)
nu11 <- sapply(1:length(plot_seq), function(x) whitt_version(h = plot_seq[x], nu1 = .5, nu2 = .5, c2 = 1,c11 = 1, c12 = 1)*
                 norm_constant(d = 1, x = plot_seq[x], nu_1 = .5, nu_2 = .5, a_1 = 1, a_2 = 1))
nu22 <- sapply(1:length(plot_seq), function(x) whitt_version(h = plot_seq[x], nu1 = .5, nu2 = .8,c2 = 1,c11 = 1, c12 = 1)*
                 norm_constant(d = 1, x = plot_seq[x], nu_1 = .5, nu_2 = .8, a_1 = 1, a_2 = 1))
nu33 <- sapply(1:length(plot_seq), function(x) whitt_version(h = plot_seq[x], nu1 = .5, nu2 = 1.1,c2 = 1,c11 = 1, c12 = 1)*
                 norm_constant(d = 1, x = plot_seq[x], nu_1 = .5, nu_2 = 1.1, a_1 = 1, a_2 = 1))
nu44 <- sapply(1:length(plot_seq), function(x) whitt_version(h = plot_seq[x], nu1 = .5, nu2 = 1.8,c2 = 1,c11 = 1, c12 = 1)*
                 norm_constant(d = 1, x = plot_seq[x], nu_1 = .5, nu_2 = 1.8, a_1 = 1, a_2 = 1))
nu55 <- sapply(1:length(plot_seq), function(x) whitt_version(h = plot_seq[x], nu1 = .5, nu2 = 2.8,c2 = 1,c11 = 1, c12 = 1)*
                 norm_constant(d = 1, x = plot_seq[x], nu_1 = .5, nu_2 = 2.8, a_1 = 1, a_2 = 1))
nu66 <- sapply(1:length(plot_seq), function(x) whitt_version(h = plot_seq[x], nu1 = .5, nu2 = 3.8,c2 = 1,c11 = 1, c12 = 1)*
                 norm_constant(d = 1, x = plot_seq[x], nu_1 = .5, nu_2 = 3.8, a_1 = 1, a_2 = 1))

df <- data.frame(plot_seq, 
                 value = c(nu11[2,], nu22[2,], nu33[2,], nu44[2,], nu55[2,], nu66[2,]),
                 nu = rep(c(.5, .8, 1.1, 1.8, 2.8, 3.8), each = length(plot_seq)))

library(ggplot2)
theme_set(theme_bw())
ggplot(data = df, aes(x = plot_seq, y = value, color = factor(nu), linetype = factor(nu))) + 
  geom_line() + 
  labs(x = 'Lags', y = 'Cross-covariance function value', color = expression(nu[2]),
       linetype = expression(nu[2]))+
  theme(legend.position = 'bottom')
ggsave('images/example_fun_whittaker.png', height = 4, width = 4)

nu11 <- sapply(1:length(plot_seq), function(x) whitt_version(h = plot_seq[x], nu1 = .5, nu2 = .8,c2 = 1,c11 = 1, c12 = 1,
                                                             a1 = 1, a2 = .2)*
                 norm_constant(d = 1, x = plot_seq[x], nu_1 = .5, nu_2 = .8, a_1 = 1, a_2 = .2))
nu22 <- sapply(1:length(plot_seq), function(x) whitt_version(h = plot_seq[x], nu1 = .5, nu2 = .8,c2 = 1,c11 = 1, c12 = 1,
                                                             a1 = 1, a2 = .5)*
                 norm_constant(d = 1, x = plot_seq[x], nu_1 = .5, nu_2 = .8, a_1 = 1, a_2 = .5))
nu33 <- sapply(1:length(plot_seq), function(x) whitt_version(h = plot_seq[x], nu1 = .5, nu2 = .8,c2 = 1,c11 = 1, c12 = 1,
                                                             a1 = 1, a2 = 1)*
                 norm_constant(d = 1, x = plot_seq[x], nu_1 = .5, nu_2 = .8, a_1 = 1, a_2 = 1))
nu44 <- sapply(1:length(plot_seq), function(x) whitt_version(h = plot_seq[x], nu1 = .5, nu2 = .8,c2 = 1,c11 = 1, c12 = 1,
                                                             a1 = 1, a2 = 2)*
                 norm_constant(d = 1, x = plot_seq[x], nu_1 = .5, nu_2 = .8, a_1 = 1, a_2 = 2))
df <- data.frame(plot_seq, 
                 value = c(nu11[2,], nu22[2,], nu33[2,], nu44[2,]),
                 nu = rep(c('0.2', '0.5','1.0', '2.0'), each = length(plot_seq)))

ggplot(data = df, aes(x = plot_seq, y = value, color = factor(nu), linetype = factor(nu))) + 
  geom_line() + 
  labs(x = 'Lags', y = 'Cross-covariance function value', color = expression(a[2]),
       linetype = expression(a[2]))+
  theme(legend.position = 'bottom')
ggsave('images/example_fun_whittaker_scale.png', height = 4, width = 4)

cov_val_lag <- sapply(1:length(plot_seq), function(x) 
  cross_cov(0, plot_seq[x], nu = .75, z_ij = complex(real = 1, imaginary = 0), a = 1) *
  norm_constant(d = 1, x = plot_seq[x], nu_1 = .75, nu_2 = .75, a_1 = 1, a_2 = 1))
cov_val_lag2 <- sapply(1:length(plot_seq), function(x) 
  cross_cov(0, plot_seq[x], nu = .75, z_ij = complex(real = 0, imaginary = 1), a = 1)*
    norm_constant(d = 1, x = plot_seq[x], nu_1 = .75, nu_2 = .75, a_1 = 1, a_2 = 1))
cov_val_lag3 <- sapply(1:length(plot_seq), function(x) 
  cross_cov(0, plot_seq[x], nu = .75, z_ij = complex(real = sqrt(2)/sqrt(3), imaginary = 1/sqrt(3)), a = 1)*
    norm_constant(d = 1, x = plot_seq[x], nu_1 = .75, nu_2 = .75, a_1 = 1, a_2 = 1))
cov_val_lag4 <- sapply(1:length(plot_seq), function(x) 
  cross_cov(0, plot_seq[x], nu = .75, z_ij = complex(real = 1/sqrt(3), imaginary = sqrt(2)/sqrt(3)), a = 1)*
    norm_constant(d = 1, x = plot_seq[x], nu_1 = .75, nu_2 = .75, a_1 = 1, a_2 = 1))


df <- data.frame(plot_seq, 
                 value = c(cov_val_lag, cov_val_lag2, cov_val_lag3, cov_val_lag4),
                 nu = factor(rep(c('1', 'i', '(2 + 1i)/sqrt(3)', '(1 + 2i)/sqrt(3)'), 
                          each = length(plot_seq)),
                          levels = c('1', '(2 + 1i)/sqrt(3)', '(1 + 2i)/sqrt(3)', 'i')))

ggplot(data = df, aes(x = plot_seq, y = value, color = factor(nu), linetype = factor(nu))) + 
  geom_line() + 
  labs(x = 'Lags', y = 'Cross-covariance function value', color = expression(Sigma[12]),
       linetype = expression(Sigma[12]))+
  theme(legend.position = 'bottom') + 
  scale_color_discrete(labels = expression(1, frac(sqrt(2) + i,sqrt(3)),frac(1 + sqrt(2)*i,sqrt(3)),i))+
  scale_linetype_discrete(labels = expression(1, frac(sqrt(2) + i, sqrt(3)),frac(1 + sqrt(2)*i,sqrt(3)), i))
ggsave('images/example_combination_function.png', height = 4, width = 7)



spec_dens <- function(x, h, nu1, nu2, a1 = 1, a2 = 1, re_z, im_z) {
  val1 <- complex(imaginary = 1) * (im_z * cos(h * x) * 
        (a1 + complex(imaginary= x))^(-nu1 - 1/2) * 
        (a2 - complex(imaginary= x))^(-nu2 - 1/2)) -
    complex(imaginary = 1)*(im_z * cos(h * x) * 
        (a2 + complex(imaginary= x))^(-nu2 - 1/2) * 
        (a1 - complex(imaginary= x))^(-nu1 - 1/2))
  val2 <-  -(im_z * sin(h * x) * 
               (a1 + complex(imaginary= x))^(-nu1 - 1/2) * 
               (a2 - complex(imaginary= x))^(-nu2 - 1/2)) +
    -(im_z * sin(h * x) * 
        (a2 + complex(imaginary= x))^(-nu2 - 1/2) * 
        (a1 - complex(imaginary= x))^(-nu1 - 1/2))
  Re(val1 + val2)
}
nu1_val <- 1.25
test25 <- sapply(plot_seq, function(z) {
  integrate(spec_dens, h = z, nu1 = nu1_val, nu2 = .25, re_z = 0, im_z = 1,
            a1 = 1, a2 = 1,
            lower = 0, upper = 100)$value*
    norm_constant(d = 1, x = plot_seq[x], nu_1 = nu1_val, nu_2 = .25, a_1 = 1, a_2 = 1)
})
test75 <- sapply(plot_seq, function(z) {
  integrate(spec_dens, h = z, nu1 = nu1_val, nu2 = .75, re_z = 0, im_z = 1,
            a1 = 1, a2 = 1,
            lower = 0, upper = 100)$value*
    norm_constant(d = 1, x = plot_seq[x], nu_1 = nu1_val, nu_2 = .75, a_1 = 1, a_2 = 1)
})
test125 <- sapply(plot_seq, function(z) {
  integrate(spec_dens, h = z, nu1 = nu1_val, nu2 = 1.25, re_z = 0, im_z = 1,
            a1 = 1, a2 = 1,
            lower = 0, upper = 100)$value*
    norm_constant(d = 1, x = plot_seq[x], nu_1 = nu1_val, nu_2 = 1.25, a_1 = 1, a_2 = 1)
})

test225 <- sapply(plot_seq, function(z) {
  integrate(spec_dens, h = z, nu1 = nu1_val, nu2 = 2.25, re_z = 0, im_z = 1,
            a1 = 1, a2 = 1,
            lower = 0, upper = 100)$value*
    norm_constant(d = 1, x = plot_seq[x], nu_1 = nu1_val, nu_2 = 2.25, a_1 = 1, a_2 = 1)
})

test325 <- sapply(plot_seq, function(z) {
  integrate(spec_dens, h = z, nu1 = nu1_val, nu2 = 3.25, re_z = 0, im_z = 1,
            a1 = 1, a2 = 1,
            lower = 0, upper = 100)$value*
    norm_constant(d = 1, x = plot_seq[x], nu_1 = nu1_val, nu_2 = 3.25, a_1 = 1, a_2 = 1)
})

df <- data.frame(plot_seq, 
                 value = c(test25, test75, test125, test225, test325),
                 nu = rep(c(.25, .75, 1.25, 2.25, 3.25), each = length(plot_seq)))

ggplot(data = df, aes(x = plot_seq, y = value, color = factor(nu), linetype = factor(nu))) + 
  geom_line() + 
  labs(x = 'Lags', y = 'Cross-covariance function value', color = expression(nu[2]),
       linetype = expression(nu[2]))+
  theme(legend.position = 'bottom')
ggsave('images/im_varying_nu.png', height = 4, width = 4)

nu1_val <- 1.25
test25 <- sapply(plot_seq, function(z) {
  integrate(spec_dens, h = z, nu1 = nu1_val, nu2 = nu1_val, re_z = 0, im_z = 1,
            a1 = 1, a2 = .6,
            lower = 0, upper = 100)$value*
    norm_constant(d = 1, x = plot_seq[x], nu_1 = nu1_val, nu_2 = nu1_val, a_1 = 1, a_2 = .6)
})
test75 <- sapply(plot_seq, function(z) {
  integrate(spec_dens, h = z, nu1 = nu1_val, nu2 = nu1_val, re_z = 0, im_z = 1,
            a1 = 1, a2 = .8,
            lower = 0, upper = 100)$value*
    norm_constant(d = 1, x = plot_seq[x], nu_1 = nu1_val, nu_2 = nu1_val, a_1 = 1, a_2 = .8)
})
test125 <- sapply(plot_seq, function(z) {
  integrate(spec_dens, h = z, nu1 = nu1_val, nu2 = nu1_val, re_z = 0, im_z = 1,
            a1 = 1, a2 = 1,
            lower = 0, upper = 100)$value*
    norm_constant(d = 1, x = plot_seq[x], nu_1 = nu1_val, nu_2 = nu1_val, a_1 = 1, a_2 = 1)
})

test225 <- sapply(plot_seq, function(z) {
  integrate(spec_dens, h = z, nu1 = nu1_val, nu2 = nu1_val, re_z = 0, im_z = 1,
            a1 = 1, a2 = 1.2,
            lower = 0, upper = 100)$value*
    norm_constant(d = 1, x = plot_seq[x], nu_1 = nu1_val, nu_2 = nu1_val, a_1 = 1, a_2 = 1.2)
})

test325 <- sapply(plot_seq, function(z) {
  integrate(spec_dens, h = z, nu1 = nu1_val, nu2 = nu1_val, re_z = 0, im_z = 1,
            a1 = 1, a2 = 1.4,
            lower = 0, upper = 100)$value*
    norm_constant(d = 1, x = plot_seq[x], nu_1 = nu1_val, nu_2 = nu1_val, a_1 = 1, a_2 = 1.4)
})

df <- data.frame(plot_seq, 
                 value = c(test25, test75, test125, test225, test325),
                 nu = rep(c(.6, .8, 1, 1.2, 1.4), each = length(plot_seq)))

ggplot(data = df, aes(x = plot_seq, y = value, color = factor(nu), linetype = factor(nu))) + 
  geom_line() + 
  labs(x = 'Lags', y = 'Cross-covariance function value', color = expression(a[2]),
       linetype = expression(a[2]))+
  theme(legend.position = 'bottom')
ggsave('images/im_varying_a.png', height = 4, width = 4)



