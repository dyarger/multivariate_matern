# informal comparison of the difference computing methods
source('code/multi_matern_source.R')
library(fields)

# d = 1

# exponential covariance
grid_info_1d <- create_grid_info_1d(2^14, 10)
nu <- .5
a <- 1
df <- fft_1d(nu1 = nu, nu2 = nu, a1 = a, a2 = a, Sigma_re = 1, Sigma_im = 0, grid_info = grid_info_1d)
ii <- abs(df[,1]) < 5
x <- df[ii,1]
fft_fun <- df[ii,2]
exp_fun <- exp(-abs(df[ii,1]))
int_fun <- sapply(df[ii,1], function(y) integrate(spec_dens, h = y, nu1 = nu, nu2 = nu, a1 = 1, a2 = 1, re_z = 1, im_z = 0,
                                                  lower = -10^2, upper = 10^2)$value *
                    norm_constant(nu1 = nu, nu2 = nu, a1 = 1, a2 = 1, d = 1))
plot(x, fft_fun, type = 'l')
lines(x, exp_fun, col = 2)
lines(x, int_fun, col = 3)
mean((fft_fun - exp_fun)^2)
mean((int_fun - exp_fun)^2)
mean((fft_fun - int_fun)^2)
x_seq <- seq(-3, 3, length.out = 100)
library(fAsianOptions)
grid_info_1d_smaller <- create_grid_info_1d(2^12, 10)
grid_info_1d_smallest <- create_grid_info_1d(2^10, 10)

nu <- .5
a <- 1
exp_vals <-     sapply(x_seq, function(h) {1 * exp(-a * abs(h))})

matern_vals <- sapply(x_seq, function(h) {2^(1 - nu)/gamma(nu) * (a * abs(h))^nu * 
      besselK(a * abs(h), nu)})
mean((matern_vals - exp_vals)^2)

df <- fft_1d(nu1 = nu, nu2 = nu, a1 = a, a2 = a,
             Sigma_re = 1, Sigma_im = 0, grid_info = grid_info_1d_smaller)
fft_interp_vals_smaller <- approx(xout = x_seq, y = df[,2], x = df[,1], method = 'linear')$y
mean((fft_interp_vals_smaller - exp_vals)^2)

df <- fft_1d(nu1 = nu, nu2 = nu, a1 = a, a2 = a,
             Sigma_re = 1, Sigma_im = 0, grid_info = grid_info_1d)
fft_interp_vals <- approx(xout = x_seq, y = df[,2], x = df[,1], method = 'linear')$y
mean((fft_interp_vals - exp_vals)^2)

df <- fft_1d(nu1 = nu, nu2 = nu, a1 = a, a2 = a,
             Sigma_re = 1, Sigma_im = 0, grid_info = grid_info_1d_smallest)
fft_interp_vals_smallest <- approx(xout = x_seq, y = df[,2], x = df[,1], method = 'linear')$y
mean((fft_interp_vals_smallest - exp_vals)^2)

fAsianOptions_vals <-     Re(sapply(x_seq, function(h) {
  a1 <- a
  a2 <- a
  nu1 = nu
  nu2 <- nu + .000000001
  Sigma12 <- 1
  if (h < 0) {
    Re(2 * pi * Sigma12 * 1/gamma(nu2 + 1/2) * abs(h)^(nu1/2 + nu2/2 - 1/2) * (a1 + a2)^(-nu1/2 - nu2/2 - 1/2) *
      exp((a1 - a2)/2 * -h) * 
      fAsianOptions::whittakerW(x = (a1 + a2)*abs(h),  kappa = -nu1/2 + nu2/2, mu = (nu1 + nu2)/2 , ip = 0)) *
      norm_constant(nu1, nu2, a1 = a1, a2 = a2, d = 1, norm_type = 'A')
  } else {
    Re(
      2 * pi * Sigma12 * 1/gamma(nu1 + 1/2) * abs(h)^(nu1/2 + nu2/2 - 1/2) *
        (a1 + a2)^(-nu1/2 - nu2/2 - 1/2) * exp((a1 - a2)/2 * -h) *  
        whittakerW(x = (a1 + a2)*abs(h), kappa = nu1/2 - nu2/2,mu = (nu1 + nu2)/2 , ip = 0)
    ) * norm_constant(nu1, nu2, a1 = a1, a2 = a2, d = 1, norm_type = 'A')
  }
}))
mean((fAsianOptions_vals - exp_vals)^2)

integrate_vals <-     Re(sapply(x_seq, function(h) {
  sapply(x_seq, function(y) integrate(spec_dens, h = y, nu1 = nu, nu2 = nu, a1 = a, a2 = a, re_z = 1, im_z = 0,
                                      lower = -10^2, upper = 10^2)$value *
           norm_constant(nu1 = nu, nu2 = nu, a1 = a, a2 = a, d = 1))
}))
mean((integrate_vals - exp_vals)^2)



rbenchmark::benchmark(
  fft_only = {
    df <- fft_1d(nu1 = nu, nu2 = nu, a1 = a, a2 = a,
                 Sigma_re = 1, Sigma_im = 0, grid_info = grid_info_1d)
  },
  fft_interp = {
    df <- fft_1d(nu1 = nu, nu2 = nu, a1 = a, a2 = a,
                 Sigma_re = 1, Sigma_im = 0, grid_info = grid_info_1d)
    yout <- approx(xout = x_seq, y = df[,2], x = df[,1], method = 'linear')$y
    yout
  }, 
  fft_only_smaller = {
    df <- fft_1d(nu1 = nu, nu2 = nu, a1 = a, a2 = a,
                 Sigma_re = 1, Sigma_im = 0, grid_info = grid_info_1d_smaller)
  },
  fft_interp_smaller = {
    df <- fft_1d(nu1 = nu, nu2 = nu, a1 = a, a2 = a,
                 Sigma_re = 1, Sigma_im = 0, grid_info = grid_info_1d_smaller)
    yout <- approx(xout = x_seq, y = df[,2], x = df[,1], method = 'linear')$y
    yout
  }, 
  fft_interp_smallest = {
    df <- fft_1d(nu1 = nu, nu2 = nu, a1 = a, a2 = a,
                 Sigma_re = 1, Sigma_im = 0, grid_info = grid_info_1d_smallest)
    yout <- approx(xout = x_seq, y = df[,2], x = df[,1], method = 'linear')$y
    yout
  }, 
  integrate = {
    sapply(x_seq, function(y) integrate(spec_dens, h = y, nu1 = nu, nu2 = nu, a1 = a, a2 = a, re_z = 1, im_z = 0,
                                        lower = -10^2, upper = 10^2)$value *
             norm_constant(nu1 = nu, nu2 = nu, a1 = a, a2 = a, d = 1))
  }, 
  fAsianOptions = {
    a1 <- a
    a2 <- a
    nu1 = nu
    nu2 <- nu + .000000001
    Sigma12 <- 1
    Re(sapply(x_seq, function(h) {
      2 * pi * Sigma12 * 1/gamma(nu1 + 1/2) * abs(h)^(nu1/2 + nu2/2 - 1/2) *
        (a1 + a2)^(-nu1/2 - nu2/2 - 1/2) * exp((a1 - a2)/2 * -h) *  
        whittakerW(x = (a1 + a2)*abs(h), kappa = nu1/2 - nu2/2,mu = (nu1 + nu2)/2 , ip = 10)
    }))
  },
  fAsianOptions_ip0 = {
    a1 <- a
    a2 <- a
    nu1 = nu
    nu2 <- nu + .000000001
    Sigma12 <- 1
    Re(sapply(x_seq, function(h) {
      2 * pi * Sigma12 * 1/gamma(nu1 + 1/2) * abs(h)^(nu1/2 + nu2/2 - 1/2) *
        (a1 + a2)^(-nu1/2 - nu2/2 - 1/2) * exp((a1 - a2)/2 * -h) *  
        whittakerW(x = (a1 + a2)*abs(h), kappa = nu1/2 - nu2/2,mu = (nu1 + nu2)/2 , ip = 0)
    }))
  }, 
  exponential = {
    sapply(x_seq, function(h) {1 * exp(-a * abs(h))})
  },
  matern = {
    sapply(x_seq, function(h) {2^(1 - nu)/gamma(nu) * (a * abs(h))^nu * 
        besselK(a * abs(h), nu)})
  },
  replications = 500, order = 'elapsed')
### other various plots


# Matern covariance
nu <- .88
df <- fft_1d(nu1 = nu, nu2 = nu, a1 = a, a2 = a, Sigma_re = 1, Sigma_im = 0, grid_info = grid_info_1d)
fields_fun <- fields::Matern(d = abs(x), smoothness = nu)
our_matern_fun <- sapply(x, function(y) re_part_d1(h = abs(y), nu1 = nu,nu2 = nu, a1 = 1, a2 = 1, Sigma12 = 1))
int_fun <- sapply(df[ii,1], function(y) integrate(spec_dens, h = y, nu1 = nu, nu2 = nu, a1 = 1, a2 = 1, re_z = 1, im_z = 0,
                                                  lower = -10^2, upper = 10^2)$value *
                    norm_constant(nu1 = nu, nu2 = nu, a1 = 1, a2 = 1, d = 1))
fft_fun <- df[ii,2]
plot(x, fft_fun, type = 'l')
lines(x, fields_fun, col = 2)
lines(x, our_matern_fun, col = 3)
lines(x, int_fun, col = 4)
mean((fft_fun - fields_fun)^2)
mean((fft_fun - our_matern_fun)^2)
mean((fields_fun - our_matern_fun)^2)
mean((int_fun - our_matern_fun)^2)

# Whittaker function
a = 1
df <- fft_1d(nu1 = nu + .6, nu2 = nu, a1 = a - .2, a2 = a, Sigma_re = 1, Sigma_im = 0, grid_info = grid_info_1d)
fft_fun <- df[ii,2]
whitt_fun <- sapply(x, function(y) re_part_d1(y, nu1 = nu + .6, nu2 = nu, a1 = a - .2,
                                                     a2 = a, Sigma12 = 1))
int_fun <- sapply(df[ii,1], function(y) integrate(spec_dens, h = y, nu1 = nu + .6, nu2 = nu, a1 = a - .2, a2 = a, re_z = 1, im_z = 0,
                                                  lower = -10^2, upper = 10^2)$value *
                    norm_constant(nu1 = nu + .6, nu2 = nu, a1 = a - .2, a2 = a, d = 1))
plot(x, fft_fun, type = 'l')
lines(x, whitt_fun, col = 2)
lines(x, int_fun, col = 3)
mean((fft_fun - whitt_fun)^2)
mean((int_fun - whitt_fun)^2)
mean((fft_fun - int_fun)^2)

grid_info_1d <- create_grid_info_1d(2^14, 20)
library(fAsianOptions)
rbenchmark::benchmark(
  fft_only = {
    df <- fft_1d(nu1 = nu + .6, nu2 = nu, a1 = a - .2, a2 = a,
               Sigma_re = 1, Sigma_im = 0, grid_info = grid_info_1d)
    },
  fft_interp = {
    df <- fft_1d(nu1 = nu + .6, nu2 = nu, a1 = a - .2, a2 = a,
                 Sigma_re = 1, Sigma_im = 0, grid_info = grid_info_1d)
    yout <- approx(xout = x_seq, y = df[,2], x = df[,1], method = 'linear')$y
    yout
  }, 
  integrate = {
    sapply(x_seq, function(y) integrate(spec_dens, h = y, nu1 = nu + .6, nu2 = nu, a1 = a - .2, a2 = a, re_z = 1, im_z = 0,
                                           lower = -10^2, upper = 10^2)$value *
             norm_constant(nu1 = nu + .6, nu2 = nu, a1 = a - .2, a2 = a, d = 1))
  }, 
  fAsianOptions = {
    a1 <- a - .2
    a2 <- a
    nu1 = nu + .6
    nu2 <- nu
    Sigma12 <- 1
    Re(sapply(x_seq, function(h) {
      2 * pi * Sigma12 * 1/gamma(nu1 + 1/2) * abs(h)^(nu1/2 + nu2/2 - 1/2) *
        (a1 + a2)^(-nu1/2 - nu2/2 - 1/2) * exp((a1 - a2)/2 * -h) *  
      whittakerW(x = (a1 + a2)*abs(h), kappa = nu1/2 - nu2/2,mu = (nu1 + nu2)/2 , ip = 10)
  }))
    },
  fAsianOptions_ip0 = {
    a1 <- a - .2
    a2 <- a
    nu1 = nu + .6
    nu2 <- nu
    Sigma12 <- 1
    Re(sapply(x_seq, function(h) {
      2 * pi * Sigma12 * 1/gamma(nu1 + 1/2) * abs(h)^(nu1/2 + nu2/2 - 1/2) *
        (a1 + a2)^(-nu1/2 - nu2/2 - 1/2) * exp((a1 - a2)/2 * -h) *  
        whittakerW(x = (a1 + a2)*abs(h), kappa = nu1/2 - nu2/2,mu = (nu1 + nu2)/2 , ip = 0)
    }))
  },replications = 100)



#struve function
df <- fft_1d(nu1 = nu, nu2 = nu, a1 = a, a2 = a, Sigma_re = 0, Sigma_im = 1, grid_info = grid_info_1d)
fft_fun <- df[ii,2]
struve_fun <- sapply(x, function(y) -im_part_d1(y,nu = nu, 
                                                     a = a))
int_fun <- sapply(df[ii,1], function(y) integrate(spec_dens, h = y, nu1 = nu, nu2 = nu, a1 = a, a2 = a, re_z = 0, im_z = 1,
                                                  lower = -10^2, upper = 10^2)$value *
                    norm_constant(nu1 = nu, nu2 = nu, a1 = a, a2 = a, d = 1))
plot(x, fft_fun, type = 'l')
lines(x, struve_fun, col = 2)
lines(x, int_fun, col = 3)
mean((fft_fun - struve_fun)^2)
mean((int_fun - struve_fun)^2)
mean((fft_fun - int_fun)^2)

# imaginary part, general 
df <- fft_1d(nu1 = nu + 1.6, nu2 = nu, a1 = a + .2, a2 = a, Sigma_re = 0, Sigma_im = 1, grid_info = grid_info_1d)
fft_fun <- df[ii,2]
int_fun <- sapply(df[ii,1], function(y) integrate(spec_dens, h = y, nu1 = nu + 1.6, nu2 = nu, a1 = a + .2, a2 = a, re_z = 0, im_z = 1,
                                                  lower = -10^2, upper = 10^2)$value *
                    norm_constant(nu1 = nu + 1.6, nu2 = nu, a1 = a + .2, a2 = a, d = 1))
plot(x, fft_fun, type = 'l')
lines(x, int_fun, col = 2)
mean((fft_fun - int_fun)^2)


# full 1d cross-covariance
df <- fft_1d(nu1 = nu + 1.6, nu2 = nu, a1 = a + .2, a2 = a, Sigma_re = .2, Sigma_im = .2, grid_info = grid_info_1d)
fft_fun <- df[ii,2]
int_fun <- sapply(df[ii,1], function(y) integrate(spec_dens, h = y, nu1 = nu + 1.6, nu2 = nu, a1 = a + .2, a2 = a, re_z = .2, im_z = .2,
                                                  lower = -10^2, upper = 10^2)$value *
                    norm_constant(nu1 = nu + 1.6, nu2 = nu, a1 = a + .2, a2 = a, d = 1))
plot(x, fft_fun, type = 'l')
lines(x, int_fun, col = 2)
mean((fft_fun - int_fun)^2)


############ 2d ########
library(dplyr)
# set up integration challenge in polar coordinates
# we're not going to do this because it takes forever, but a reminder that it exists
angle_grid <- seq(0, 2 * pi, length.out = 100)
r_grid <- exp(seq(-9, 8, length.out = 800))
r_lag <- data.frame('r' = r_grid,
                    'r_lag' = c(r_grid[2] - r_grid[1], (r_grid - dplyr::lag(r_grid))[-1]))
angle_lag <- data.frame('angle' = angle_grid,
                        'angle_lag' = c(angle_grid[2] - angle_grid[1], (angle_grid - dplyr::lag(angle_grid))[-1]))
approx_grid <- expand.grid('angle' = angle_grid,
                           'r' = r_grid) %>%
  left_join(r_lag) %>% left_join(angle_lag)

approx_grid$x <- approx_grid$r * cos(approx_grid$angle)
approx_grid$y <- approx_grid$r * sin(approx_grid$angle)
approx_grid$theta_x <- cos(approx_grid$angle)
approx_grid$theta_y <- sin(approx_grid$angle)

# Matern
grid_info_2d <- create_grid_info_2d(n_points = 2^9, x_max = 10)
df <- fft_2d(nu1 = 1.5, nu2 = 1.5, a1 = 1, a2 = 1, Psi = sign, Mu = function(x) 1, grid_info = grid_info_2d)
ii <- df[['Var2']] == min(abs(df[['Var2']]))
min(abs(df[['Var2']]))
x <- df[['Var1']][ii]
fft_fun <- df[['val']][ii]
fields_fun <- Matern(abs(x), nu = 1.5)

# only integrate at a few places
n_compare <- 500
indexes_use <- sample(1:nrow(df), size = n_compare, replace = FALSE)
int_fun <- sapply(indexes_use, function(x, a1, a2, nu1, nu2, Mu, Psi,  approx_grid,
                                       d = 2) {
  h = as.double(df[x,1:2])
  spatial_integrate_d2(h = h, a1, a2, nu1 = nu1, nu2 = nu2, Mu = Mu, Psi = Psi, approx_grid = approx_grid)
}, a1 = 1, a2 = 1, nu1 = 1.5, nu2 = 1.5,
Mu = function(theta_x, theta_y, i1,i2) 1, Psi = function(theta_x, theta_y) sign(theta_x), approx_grid = approx_grid, d = 2
)


plot(x, fft_fun, type = 'l')
lines(x, col = 2, fields_fun)
mean((fft_fun - fields_fun)^2)
mean((df[['val']][indexes_use] - int_fun)^2)
plot(df[['val']][indexes_use], int_fun)


fft_mat <- matrix(df[['val']], nrow = sqrt(length(ii)))
fields_mat <- matrix(Matern(sqrt(df[['Var1']]^2 + df[['Var2']]^2), nu = 1.5) , nrow = sqrt(length(ii)))
image.plot(fft_mat)
image.plot(fields_mat)
image.plot(fft_mat - fields_mat)
mean((fft_mat - fields_mat)^2)

# with imaginary directional measure
df <- fft_2d(nu1 = 1.1, nu2 = 1.5, a1 = 1.1, a2 = 1, Psi = sign, 
             Mu = function(y) complex(real = .2, imaginary = .2 * sign(y)), 
             grid_info = grid_info_2d)
# only integrate at a few places
n_compare <- 500
indexes_use <- sample(1:nrow(df), size = n_compare, replace = FALSE)
int_fun <- sapply(indexes_use, function(x, a1, a2, nu1, nu2, Mu, Psi,  approx_grid,
                                       d = 2) {
  h = as.double(df[x,2:1])
  spatial_integrate_d2(h = h, a1, a2, nu1 = nu1, nu2 = nu2, Mu = Mu, Psi = Psi, approx_grid = approx_grid)
}, a1 = 1.1, a2 = 1, nu1 = 1.1, nu2 = 1.5,
Mu = function(theta_x, theta_y, i1, i2) complex(real = .2, imaginary = .2 * sign(theta_x)), Psi = function(theta_x, theta_y) sign(theta_x), approx_grid = approx_grid, d = 2
)
plot(df[['val']][indexes_use], int_fun)

