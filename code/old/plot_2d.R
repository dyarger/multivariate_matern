
library(ggplot2)
library(scales)
theme_set(theme_bw())
spatial_integrate <- function(h,d, a_1, a_2, nu_1, nu_2) {
  spec_dens_single <- function(x, h, d, a_1, a_2, nu_1, nu_2) {
    Re(besselJ(h*x, d/2 - 1) * x^(-(d/2 - 1)) * (a_1 + complex(imaginary = x))^(-nu_1-d/2)*
      (a_2 - complex(imaginary = x))^(-nu_2-d/2)) * x^(d-1)
  }
  test_integrate <- integrate(lower = .0001, upper = 1000, 
                              f = spec_dens_single,
                              h = h, d = d, a_1 = a_1, a_2 = a_2, nu_1 = nu_1, 
                              nu_2 = nu_2,subdivisions = 2000)
  (2 * pi)^(d/2) *h^(-d/2 + 1) * test_integrate[['value']] #/    (a_1^(nu_1) * a_2^(nu_2))
}
a_1 <- a_2 <- 1
d <- 2


lag_seq <- seq(0, 5, by = .005)
res <- sapply(lag_seq, spatial_integrate, d = d, a_1 = a_1, a_2 = a_2, nu_1 = .5, nu_2 = .5)
res2 <- sapply(lag_seq, spatial_integrate, d = d, a_1 = a_1, a_2 = a_2, nu_1 = 1, nu_2 = .5)
res3 <- sapply(lag_seq, spatial_integrate, d = d, a_1 = a_1, a_2 = a_2, nu_1 = 1, nu_2 = 1)
res4 <- sapply(lag_seq, spatial_integrate, d = d, a_1 = a_1, a_2 = a_2, nu_1 = .75, nu_2 = .75)
df <- data.frame(lag_seq, 
                 value = c(res, res2, res3, res4),
                 nu1 = rep(c(.5, .5, 1, .75), each = length(lag_seq)),
                 nu2 = rep(c(.5, 1, 1, .75), each = length(lag_seq)), 
                 nu_label = rep(c('nu[j]=0.2 nu[k]=0.2', 'nu[j]=0.2~nu[k]=0.5', 
                                  'nu[j]=0.5~nu[k]=0.5',  'nu[j]=0.5~nu[k]=0.5'), each = length(lag_seq)))

ggplot(data = df, aes(x = lag_seq, y = value, 
                      color = paste0(nu1, ', ', nu2), 
                      linetype = paste0(nu1, ', ', nu2))) + 
  geom_line() + 
  labs(x = 'Lags', y = 'Cross-covariance function value', #color = 'Parameters',
       linetype = expression(paste('Values of ', nu[j], ', ', nu[k])),
       color = expression(paste('Values of ', nu[j], ', ', nu[k])))+
  theme(legend.position = 'bottom') + 
  guides(linetype=  guide_legend(nrow = 2),
         color=  guide_legend(nrow = 2))
ggsave('images/symmetric_spat.png', height = 4, width = 4)

lag_seq <- seq(0, 5, by = .005)
res <- sapply(lag_seq, spatial_integrate, d = d, a_1 = a_1, a_2 = a_2, nu_1 = .5, nu_2 = .25)
res2 <- sapply(lag_seq, spatial_integrate, d = d, a_1 = a_1, a_2 = a_2, nu_1 = .5, nu_2 = .75)
res3 <- sapply(lag_seq, spatial_integrate, d = d, a_1 = a_1, a_2 = a_2, nu_1 = .5, nu_2 = 1.25)
res4 <- sapply(lag_seq, spatial_integrate, d = d, a_1 = a_1, a_2 = a_2, nu_1 = .5, nu_2 = 1.75)
df <- data.frame(lag_seq, 
                 value = c(res, res2, res3, res4),
                 nu1 = rep(c(.5, .5, .5, .5), each = length(lag_seq)),
                 nu2 = rep(c(.25, .75, 1.25, 1.75), each = length(lag_seq)), 
                 nu_label = rep(c('nu[2]=0.25', 'nu[2]=0.75', 
                                  'nu[2]=1.25',  'nu[2]=1.75'), each = length(lag_seq)))

ggplot(data = df, aes(x = lag_seq, y = value, 
                      color = as.character(nu2), 
                      linetype = as.character( nu2))) + 
  geom_line() + 
  labs(x = 'Lags', y = 'Cross-covariance function value', #color = 'Parameters',
       linetype = expression(paste('Values of ', nu[k])),
       color = expression(paste('Values of ', nu[k])))+
  theme(legend.position = 'bottom') + 
  guides(linetype=  guide_legend(nrow = 2),
         color=  guide_legend(nrow = 2))
ggsave('images/symmetric_spat2.png', height = 4, width = 4)

lag_seq <- seq(0, 5, by = .005)
res <- sapply(lag_seq, spatial_integrate, d = d, a_1 = .2, a_2 = a_2, nu_1 = .5, nu_2 = .5)
res2 <- sapply(lag_seq, spatial_integrate, d = d, a_1 = .5, a_2 = a_2, nu_1 = .5, nu_2 = .5)
res3 <- sapply(lag_seq, spatial_integrate, d = d, a_1 = 1, a_2 = a_2, nu_1 = .5, nu_2 = .5)
res4 <- sapply(lag_seq, spatial_integrate, d = d, a_1 = 1.5, a_2 = a_2, nu_1 = .5, nu_2 = .5)
res5 <- sapply(lag_seq, spatial_integrate, d = d, a_1 = 1.8, a_2 = a_2, nu_1 = .5, nu_2 = .5)
df <- data.frame(lag_seq, 
                 value = c(res, res2, res3, res4, res5),
                 nu1 = rep(c(.5, .5, .5, .5, .5), each = length(lag_seq)),
                 nu2 = rep(c(.2, .5, 1, 1.5, 1.8), each = length(lag_seq)), 
                 nu_label = rep(c('nu[2]=0.25', 'nu[2]=0.75', 
                                  'nu[2]=1.25',  'nu[2]=1.75', ''), each = length(lag_seq)))

ggplot(data = df, aes(x = lag_seq, y = value, 
                      color = as.character(nu2), 
                      linetype = as.character( nu2))) + 
  geom_line() + 
  labs(x = 'Lags', y = 'Cross-covariance function value', #color = 'Parameters',
       linetype = expression(paste('Values of ', a[k])),
       color = expression(paste('Values of ', a[k])))+
  theme(legend.position = 'bottom') + 
  guides(linetype=  guide_legend(nrow = 2),
         color=  guide_legend(nrow = 2))
ggsave('images/symmetric_spat_a.png', height = 4, width = 4)


# test imaginary part?

spatial_integrate <- function(h,d, a_1, a_2, nu_1, nu_2) {
  spec_dens_single <- function(x, h, d, a_1, a_2, nu_1, nu_2) {
    (besselJ(h*x, d/2 - 1) * x^(-(d/2 - 1)) * (a_1 + complex(imaginary = x))^(-nu_1-d/2)*
         (a_2 - complex(imaginary = x))^(-nu_2-d/2)) * x^(d-1)
  }
  test_integrate <- integrate(lower = .0001, upper = 1000, 
                              f = spec_dens_single,
                              h = h, d = d, a_1 = a_1, a_2 = a_2, nu_1 = nu_1, 
                              nu_2 = nu_2,subdivisions = 2000)
  (2 * pi)^(d/2) *h^(-d/2 + 1) * test_integrate[['value']] #/    (a_1^(nu_1) * a_2^(nu_2))
}
a_1 <- a_2 <- 1
d <- 2
(spec_dens_single(x = 1, h = 1, d = 2, a_1 = 1, a_2 = 1, nu_1 = .4, nu_2 = .78))


lag_seq <- seq(0, 5, by = .005)
res <- sapply(lag_seq, spatial_integrate, d = d, a_1 = a_1, a_2 = a_2, nu_1 = .5, nu_2 = .5)



spatial_integrate_d2 <- function(h, a_1, a_2, nu_1, nu_2, e_1, e_2, Delta, approx_seq) {
  spec_dens_theta <- function(x, h, e_1, e_2, Delta) {
    f_theta <- function(theta, x, h, d, e_1, e_2, Delta) {
      theta_x <- cos(theta);theta_y <- sin(theta)
      Re(exp(complex(imaginary = h[1]*x^e_1 * theta_x + h[2]*x^e_2 * theta_y)) *
           Delta(theta_x, theta_y, 1, 2))
    }
    f_theta_im <- function(theta, x, h, d, e_1, e_2, Delta) {
      theta_x <- cos(theta);theta_y <- sin(theta)
      Im(exp(complex(imaginary = h[1]*x^e_1 * theta_x + h[2]*x^e_2 * theta_y)) * 
           Delta(theta_x, theta_y, 1, 2))
    }
    Re_part <- integrate(lower = 0, upper = 2 * pi, f = f_theta, stop.on.error = F,
                         x = x, h = h, e_1 = e_1, e_2 = e_2, Delta = Delta)$value
    Im_part <- integrate(lower = 0, upper = 2 * pi, f = f_theta_im, stop.on.error = F,
                         x = x, h = h, e_1 = e_1, e_2 = e_2, Delta = Delta)$value
    complex(real = Re_part, imaginary = Im_part)
  }
  spec_dens <- function(x, h, a_1, a_2, nu_1, nu_2, e_1, e_2, Delta) {
    Re(spec_dens_theta(x, h, e_1, e_2, Delta) * 
         (a_1 + complex(imaginary = x))^(-nu_1-2/2)*
         (a_2 - complex(imaginary = x))^(-nu_2-2/2)) * x^(d-1)
  }
  approx_seq_lag <- c(approx_seq[1], approx_seq[2:length(approx_seq)] - 
    approx_seq[1:(length(approx_seq)-1)])
  test_final <- sapply(approx_seq,  FUN = spec_dens,
                       h = h, a_1 = a_1, a_2 = a_2, nu_1 = nu_1, 
                       nu_2 = nu_2, e_1 = e_1, e_2 = e_2, Delta = Delta)
  sum(test_final * approx_seq_lag)
}

e_1 <- e_2 <- 1
nu_2 <- nu_1 <- .5

lag_seq <- seq(-1, 1, by = .05)
grid <- expand.grid(lag_seq, 0)

approx_seq <- exp(seq(-10, 7, length.out = 300))
Delta <- function(theta_x, theta_y, entry_1, entry_2) {
  1
}
# make sure it matches with Matern
res <- sapply(1:nrow(grid), function(x, a_1, a_2, nu_1, nu_2, e_1, e_2,
                                     Delta, approx_seq) {
  h = as.double(grid[x,])
  spatial_integrate_d2(h = h, a_1, a_2, nu_1 = nu_1, nu_2 = nu_2, e_1 = e_1, e_2 = e_2, Delta = Delta, approx_seq = approx_seq)
}, a_1 = a_1, a_2 = a_2, nu_1 = nu_1, nu_2 = nu_2, e_1 = e_1, e_2 = e_2, Delta = Delta,approx_seq = approx_seq
)
plot(res)
plot(fields::Matern(abs(lag_seq), nu = .5))
plot(res, fields::Matern(abs(lag_seq), nu = .5))
abline(lm(data.frame(x = res, y =  fields::Matern(abs(lag_seq), nu = .5)), formula = y~x))

# now try more complicated things
Delta <- function(theta_x, theta_y, entry_1, entry_2) {
  complex(imaginary = sign(theta_x) * .3)
}

lag_seq <- seq(-3, 3, by = .025)
grid <- expand.grid(lag_seq, lag_seq)
res <- lapply(1:nrow(grid), function(x, a_1, a_2, nu_1, nu_2, e_1, e_2,
                                     Delta, approx_seq) {
  h = as.double(grid[x,])
  if (x %% 100 == 0) {
    print(x)
    print(h)
  }
  spatial_integrate_d2(h = h, a_1, a_2, nu_1 = nu_1, nu_2 = nu_2, e_1 = e_1, e_2 = e_2, Delta = Delta, approx_seq = approx_seq)
}, a_1 = a_1, a_2 = a_2, nu_1 = nu_1, nu_2 = nu_2, e_1 = e_1, e_2 = e_2, Delta = Delta,approx_seq = approx_seq
)
res_mat <- matrix( nrow = length(lag_seq), ncol = length(lag_seq), unlist(res))
fields::image.plot(res_mat)

grid_vals <- data.frame(grid, res_vals = as.double(unlist(res)))
save(grid_vals, approx_seq, file = paste0('results/asymmetric_2d_', length(approx_seq), '.RData'))
ggplot(data = grid_vals) + 
  geom_tile(aes(x = Var1, y = Var2, fill = res_vals)) + 
  scale_fill_gradient2() + 
  labs(x = expression(t[1]), y = expression(t[2]),
       fill = 'Cross\nCovariance') + 
  theme(legend.position = 'bottom')
ggsave(paste0('images/asymmetric_2d_', length(approx_seq), '.png'), height = 4, width = 4)

Delta <- function(theta_x, theta_y, entry_1, entry_2) {
  1
}

res <- lapply(1:nrow(grid), function(x, a_1, a_2, nu_1, nu_2, e_1, e_2,
                                       Delta, approx_seq) {
  h = as.double(grid[x,])
  if (x %% 100 == 0) {
    print(x)
    print(h)
  }
  spatial_integrate_d2(h = h, a_1, a_2, nu_1 = nu_1, nu_2 = nu_2, e_1 = e_1, e_2 = e_2, Delta = Delta, approx_seq = approx_seq)
}, a_1 = a_1, a_2 = a_2, nu_1 = nu_1, nu_2 = nu_2, e_1 = 1/4, e_2 = e_2, Delta = Delta,approx_seq = approx_seq
)
res_mat <- matrix( nrow = length(lag_seq), ncol = length(lag_seq), unlist(res))
fields::image.plot(res_mat)

grid_vals <- data.frame(grid, res_vals = as.double(unlist(res)))
ggplot(data = grid_vals) + 
  geom_tile(aes(x = Var1, y = Var2, fill = res_vals)) + 
  scale_fill_gradientn(colors = rainbow(5, rev = T)) + 
  labs(x = expression(t[1]), y = expression(t[2]),
       fill = 'Cross\nCovariance') + 
  theme(legend.position = 'bottom')
ggsave(paste0('images/asymmetric_2d_E_', length(approx_seq), '.png'), height = 4, width = 4)
save(grid_vals, approx_seq, file = paste0('results/asymmetric_2d_E_', length(approx_seq), '.RData'))

Delta <- function(theta_x, theta_y, entry_1, entry_2) {
  sign(theta_x)*sign(theta_y) * .3
}

res <- lapply(1:nrow(grid), function(x, a_1, a_2, nu_1, nu_2, e_1, e_2,
                                       Delta, approx_seq) {
  h = as.double(grid[x,])
  if (x %% 100 == 0) {
    print(x)
  }
  spatial_integrate_d2(h = h, a_1, a_2, nu_1 = nu_1, nu_2 = nu_2, e_1 = e_1, e_2 = e_2, Delta = Delta, approx_seq = approx_seq)
}, a_1 = a_1, a_2 = a_2, nu_1 = nu_1, nu_2 = nu_2, e_1 = e_1, e_2 = e_2, Delta = Delta,approx_seq = approx_seq
)
res_mat <- matrix( nrow = length(lag_seq), ncol = length(lag_seq), unlist(res))
fields::image.plot(res_mat)

grid_vals <- data.frame(grid, res_vals = as.double(unlist(res)))
ggplot(data = grid_vals) + 
  geom_tile(aes(x = Var1, y = Var2, fill = res_vals)) + 
  scale_fill_gradient2(midpoint = 0) + 
  labs(x = expression(t[1]), y = expression(t[2]),
       fill = 'Cross\nCovariance') + 
  theme(legend.position = 'bottom')
ggsave(paste0('images/asymmetric_2d_quadrant_', length(approx_seq), '.png'), height = 4, width = 4)
save(grid_vals, approx_seq, file = paste0('results/asymmetric_2d_quadrant_', length(approx_seq), '.RData'))
