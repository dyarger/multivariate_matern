spec_dens <- function(r, h, nu1, nu2, a1 = 1, a2 = 1, re_z, im_z) {
  val <- exp(complex(imaginary = h * r)) *
    complex(real = a1, imaginary = r)^(-nu1 - 1/2) *
    complex(real = a2, imaginary = -r)^(-nu2 - 1/2) * 
    complex(real = re_z, imaginary = sign(r) * im_z)
  Re(val)
}
re_z <- 1
im_z <- 0
a1 <- 1
a2 <- 1
nu1 <- .5
nu2 <- 0.8
plot_seq <- seq(-8, 8, by = .05)

norm_constant <- function(d, x, nu_1, nu_2, a_1 = 1, a_2 = 1) {
  a_1^(nu_1 + d/2) * a_2^(nu_2 + d/2)  / 2/pi
}

test25 <- sapply(plot_seq, function(z) {
  integrate(spec_dens, h = z, nu1 = nu1, nu2 = nu2, re_z = re_z, im_z = im_z,
            a1 = a1, a2 = a2,subdivisions = 1000,
            lower = -200, upper = 200)$value*
    norm_constant(d = 1, x = plot_seq[x], 
                  nu_1 = nu1, nu_2 = nu2, 
                  a_1 = a1, a_2 = a2)
})

testa <- sapply(plot_seq, function(z) {
  integrate(spec_dens, h = z, nu1 = nu1, nu2 = nu1, re_z = re_z, im_z = im_z,
            a1 = a1, a2 = a2,subdivisions = 1000,
            lower = -200, upper = 200)$value*
    norm_constant(d = 1, x = plot_seq[x], 
                  nu_1 = nu1, nu_2 = nu1, 
                  a_1 = a1, a_2 = a2)
})
testb <- sapply(plot_seq, function(z) {
  integrate(spec_dens, h = z, nu1 = nu2, nu2 = nu2, re_z = re_z, im_z = im_z,
            a1 = a1, a2 = a2,subdivisions = 1000,
            lower = -200, upper = 200)$value*
    norm_constant(d = 1, x = plot_seq[x], 
                  nu_1 = nu2, nu_2 = nu2, 
                  a_1 = a1, a_2 = a2)
})

plot(testa)
plot(testb)
plot(test25)
plot(test25/sqrt(testa) * sqrt(testb))
plot(test25/testa)
plot(test25/testb)

nu1 <- nu2 <- .6

test25 <- sapply(plot_seq, function(z) {
  integrate(spec_dens, h = z, nu1 = nu1, nu2 = nu2, re_z = 0, im_z = 1,
            a1 = a1, a2 = a2,subdivisions = 1000,
            lower = -200, upper = 200)$value*
    norm_constant(d = 1, x = plot_seq[x], 
                  nu_1 = nu1, nu_2 = nu2, 
                  a_1 = a1, a_2 = a2)
})

testa <- sapply(plot_seq, function(z) {
  integrate(spec_dens, h = z, nu1 = nu1, nu2 = nu1, re_z = 1, im_z = 0,
            a1 = a1, a2 = a2,subdivisions = 1000,
            lower = -200, upper = 200)$value*
    norm_constant(d = 1, x = plot_seq[x], 
                  nu_1 = nu1, nu_2 = nu1, 
                  a_1 = a1, a_2 = a2)
})
plot(test25/testb )
plot(test25/testb[plot_seq == 0] )

plot(test25)
testb <- sapply(plot_seq, function(z) {
  integrate(spec_dens, h = z, nu1 = nu2, nu2 = nu2, re_z = re_z, im_z = im_z,
            a1 = a1, a2 = a2,subdivisions = 1000,
            lower = -200, upper = 200)$value*
    norm_constant(d = 1, x = plot_seq[x], 
                  nu_1 = nu2, nu_2 = nu2, 
                  a_1 = a1, a_2 = a2)
})




library(tidyverse)
library(fields)
source('code/multi_matern_source.R')

plot_seq <- seq(-5, 5, by = .05)
nu25 <- Matern(abs(plot_seq), range = 1, smoothness = .5, phi = 1)
nu5 <- Matern(abs(plot_seq), range = 1, smoothness = 1.5, phi = 1)
nu1 <- Matern(abs(plot_seq), range = 1, smoothness = 2.5, phi = 1)


df <- data.frame(plot_seq, 
                 value = c(nu25, nu5, nu1),
                 nu = rep(c(.5, 1.5, 2.5), each = length(plot_seq)))

library(ggplot2)
theme_set(theme_bw())
ggplot(data = df , aes(x = plot_seq, y = value, color = factor(nu), linetype = factor(nu))) + 
  geom_line() + 
  labs(x = 'Lags', y = 'Covariance function value', color = expression(nu),
       linetype = expression(nu))+
  theme(legend.position = 'bottom')
ggsave('images/Matern.png', height = 4, width = 4)

ggplot(data = df %>% filter(plot_seq >= 0) , aes(x = plot_seq, y = value, color = factor(nu), linetype = factor(nu))) + 
  geom_line() + 
  labs(x = 'Lags', y = 'Covariance function value', color = expression(nu),
       linetype = expression(nu))+
  theme(legend.position = 'bottom')
ggsave('images/Matern_one_side.png', height = 4, width = 4)


