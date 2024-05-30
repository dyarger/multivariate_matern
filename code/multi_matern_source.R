
## Covariance and plotting

# the normalization constants as written in the paper
norm_constant <- function(nu1, nu2, a1 = 1, a2 = 1, d = 1, norm_type = 'A') {
  if (norm_type == 'B') {
    (a1 + a2)^(nu1 + nu2)  / 2 / pi / gamma(nu1 + nu2) * gamma(nu1 + d/2) * gamma(nu2 + d/2)
  } else if (norm_type == 'C') {
    a1^(nu1 + d/2) * a2^(nu2 + d/2) / 2 / pi
  } else if (norm_type == 'A') {
    (a1)^(nu1) * (a2)^(nu2) *
      sqrt(gamma(nu1 + d/2)) * sqrt(gamma(nu2 + d/2))/pi^(d/2)/sqrt(gamma(nu1)*gamma(nu2))
  }
}

# naive implementation of modified struve function L
struve_funL <- function(z, nu_eval) {
  if (nu_eval == -1/2) {
    return(sqrt(2/(pi * z)) * sinh(z))
  } else if (nu_eval == -3/2) {
    return(sqrt(2/pi) * (z * cosh(z) - sinh(z))/(z^(3/2)))
  }
  k <- 0:400
  (z/2)^(nu_eval + 1) * sum((z/2)^(2*k) / (gamma(k + 3/2) * gamma( k  + nu_eval + 3/2)))
}

# combine real and imaginary parts
cross_cov_same_params <- function(s, t, nu, Sigma, a) {
  Re(Sigma) * fields::Matern(abs(t - s), nu, a) - 
    Im(Sigma) * im_part_d1(t - s, nu, a)
}

im_part_d1 <- function(h, nu, a = 1) {
  if (h == 0) {
    return(0)
  }
  sign(h) * pi / 2 / cos(pi * nu) * 2^(1 - nu) / gamma(nu) *
    (a * abs(h))^nu  *
    (besselI(a*abs(h), nu = nu) - struve_funL(a*abs(h), -nu))
  # sign(h) * (abs(h)/a)^nu *
  #   (besselI(a*abs(h), nu = nu) - struve_funL(a*abs(h), -nu)) * 
  #   2^(-nu) * pi^(3/2)/cos(nu * pi)/gamma(nu + .5) *
  #   norm_constant(nu, nu, a1 = a, a2 = a, d = 1, norm_type = 'A')
}

######### whittaker function where nu1, nu2, a1, and a2 are allowed to vary ########
# only with real directional measure, bivariate case
re_part_d1 <- function(h, nu1, nu2, Sigma12, a1 = 1, a2 = 1) {
  # cross covariance function
  if (nu1 == nu2 & a1 == a2) { # if these are the same, have matern
    p12 <- Sigma12 * (2*sqrt(pi)/(gamma(nu1 + .5)))*((abs(h))^nu1) * (2*a1)^(-nu1) * 
      besselK(x = a1 * abs(h), nu = nu1)
    if (h == 0) {
      p12 <- Sigma12 * (2*sqrt(pi)/(gamma(nu1 + .5)))*((abs(1e-10))^nu1) * (2*a1)^(-nu1) * 
        besselK(x = a1 * abs(1e-10), nu = nu1)
    }
  } else if (nu1 + nu2 == round(nu1 + nu2)) {  #weird instability in whittaker function
    p12 <- 0
  } else if (h == 0) { # manually do lag of 0
    p12 <- 2 * pi*Sigma12 * 1/gamma(nu2 + 1/2) * abs(1e-10)^(nu1/2 + nu2/2 - 1/2) * 
      fAsianOptions::whittakerW(x = (a1 + a2)*abs(1e-10),kappa = -nu1/2 + nu2/2,
                                mu = -(nu1 + nu2)/2 , ip = 10) *
      (a1 + a2)^(-nu1/2 - nu2/2 - 1/2)*exp((a1 - a2)/2 * 1e-10) 
  } else if (h < 0) { # do different things for h < 0 vs h > 0
    p12 <- 2 * pi * Sigma12 * 1/gamma(nu2 + 1/2) * abs(h)^(nu1/2 + nu2/2 - 1/2) * (a1 + a2)^(-nu1/2 - nu2/2 - 1/2) *
      exp((a1 - a2)/2 * -h) * 
      fAsianOptions::whittakerW(x = (a1 + a2)*abs(h),  kappa = -nu1/2 + nu2/2, mu = (nu1 + nu2)/2 , ip = 10)
  } else {
    p12 <- 2 * pi * Sigma12 * 1/gamma(nu1 + 1/2) * abs(h)^(nu1/2 + nu2/2 - 1/2) * (a1 + a2)^(-nu1/2 - nu2/2 - 1/2) *
      exp((a1 - a2)/2 * -h) *  
      fAsianOptions::whittakerW(x = (a1 + a2)*abs(h), kappa = nu1/2 - nu2/2,mu = (nu1 + nu2)/2 , ip = 10)
  }
  Re(p12) * norm_constant(nu1, nu2, a1 = a1, a2 = a2, d = 1, norm_type = 'A')
}


### FFT work ####
# instead of exact functional forms, we implement d=1 and d=2 computations of the fft
library(fftwtools)


# 1d versions
# define fourier transform grid so we don't have to recompute thing
create_grid_info_1d <- function(n_points, x_max) {
  delta_t <-  pi / x_max
  x_vals <- 1:n_points * (2 * pi) / n_points / delta_t
  freq_points <- seq(-delta_t * n_points/2 + delta_t/2, 
                     delta_t * n_points/2 - delta_t/2, 
                     length.out = n_points)
  phase_factor <- 1/(sqrt(2*pi))  * 
    exp(complex(imaginary = freq_points[1] * 2 * pi * 
                  (1:length(freq_points)) / (delta_t * length(freq_points))))
  x_vals <- x_vals - x_max - abs(abs(x_vals[length(x_vals)] - 2*x_max) - abs(x_vals[1]))
  list('freq_points' = freq_points,
       'delta_t' = delta_t, 'n_points' = n_points, 
       'x_vals' = x_vals, 'x_max' = x_max,
       'phase_factor' = phase_factor)
}
# given a grid, compute covariance on grid with given parameters
fft_1d <- function(grid_info, nu1 = .5, nu2 = .5, a1 = 1, a2 = 1, 
                   Sigma_re, Sigma_im, norm_type = 'A') {
  # https://stackoverflow.com/questions/24077913/discretized-continuous-fourier-transform-with-numpy
  tv <- complex(real = a1, imaginary = grid_info[['freq_points']])^(-nu1 - .5) *
    complex(real = a2, imaginary = -grid_info[['freq_points']])^(-nu2 - .5) * 
    complex(real = Sigma_re, imaginary = Sigma_im*sign(grid_info[['freq_points']]))
  ff_res <- fftwtools::fftw_c2c(data = tv, inverse = 1)
  p <- length(ff_res)/2
  ff_res_adj <- c(ff_res[(p + 1):(2*p)], ff_res[1:p]) * grid_info[['phase_factor']]
  cbind(grid_info[['x_vals']], 'val' = 
          as.double(Re(ff_res_adj) * norm_constant(nu1, nu2, a1, a2, norm_type = norm_type)) / 
          grid_info[['x_max']] * grid_info[['n_points']] * sqrt(2 * pi) )
}

# given a grid, compute covariance on grid with given parameters
fft_1d_alt <- function(grid_info, nu1 = .5, nu2 = .5, a1 = 1, a2 = 1, 
                       Sigma_re, Sigma_im, norm_type = 'A') {
  # https://stackoverflow.com/questions/24077913/discretized-continuous-fourier-transform-with-numpy
  tv <- (a1^2 + grid_info[['freq_points']]^2)^(-nu1/2 - .5/2) *
    (a2^2 + grid_info[['freq_points']]^2)^(-nu2/2 - .5/2) * 
    complex(real = Sigma_re, imaginary = Sigma_im*sign(grid_info[['freq_points']]))
  ff_res <- fftwtools::fftw_c2c(data = tv, inverse = 1)
  p <- length(ff_res)/2
  ff_res_adj <- c(ff_res[(p + 1):(2*p)], ff_res[1:p]) * grid_info[['phase_factor']]
  cbind(grid_info[['x_vals']], 'val' = 
          as.double(Re(ff_res_adj) * norm_constant(nu1, nu2, a1, a2, norm_type = norm_type)) / 
          grid_info[['x_max']] * grid_info[['n_points']] * sqrt(2 * pi) )
}


# 2d versions
# define Fourier transform grid 
create_grid_info_2d <- function(n_points, x_max) {
  delta_t <-  pi / x_max
  x_vals <- 1:n_points * (2 * pi) / n_points / delta_t
  freq_points <- seq(-delta_t * n_points/2 + delta_t/2, 
                     delta_t * n_points/2 - delta_t/2, 
                     length.out = n_points)
  freq_grid <- as.data.frame(expand.grid('x' = freq_points, 'y' = freq_points))
  freq_grid$r <- sqrt(freq_grid[['x']]^2 + freq_grid[['y']]^2)
  freq_grid$theta <- atan2(freq_grid[['y']], freq_grid[['x']])
  
  phase_factor <- 1/(2*pi)  * 
    exp(complex(imaginary = rowSums(cbind(freq_grid[['x']][1], freq_grid[['y']][1]) * 2 * pi *
                                      (expand.grid((1:length(freq_points)) / (delta_t * length(freq_points)), 
                                                   (1:length(freq_points)) / (delta_t * length(freq_points))) ))))
  phase_factor_mat <- matrix(nrow = length(freq_points), ncol = length(freq_points),
                             phase_factor)
  x_vals <- x_vals - x_max - abs(abs(x_vals[length(x_vals)] - 2*x_max) - abs(x_vals[1]))
  x_vals_eg <- expand.grid(x_vals, x_vals)
  list('freq_grid' = freq_grid, 'freq_points' = freq_points,
       'delta_t' = delta_t, 'n_points' = n_points, 
       'x_vals' = x_vals, 'x_max' = x_max,
       'phase_factor_mat' = phase_factor_mat,
       'x_vals_eg' = x_vals_eg)
}
# given a grid and model parameters, computes value of Matern cross-covariance on grid
fft_2d <- function(grid_info, nu1 = .5, nu2 = .5, a1 = 1, a2 = 1, 
                   Psi, 
                   Mu, d = 2) {
  freq_grid <- grid_info[['freq_grid']]
  # https://stackoverflow.com/questions/24077913/discretized-continuous-fourier-transform-with-numpy
  tv <- complex(real = a1, imaginary = Psi(freq_grid[['theta']])*freq_grid[['r']])^(-nu1 - d/2) *
    complex(real = a2, imaginary = -Psi(freq_grid[['theta']])*freq_grid[['r']])^(-nu2 - d/2) * 
    Mu(freq_grid[['theta']])
  tv_mat <- matrix(nrow = length(grid_info[['freq_points']]), ncol = length(grid_info[['freq_points']]), tv)
  ff_res <- fftwtools::fftw_c2c_2d(data = tv_mat, inverse = 1)/grid_info[['n_points']]^2
  p <- ncol(ff_res)/2
  ff_res_adj <- cbind(rbind(ff_res[(p + 1):(2*p),(p + 1):(2*p)], ff_res[1:p,(p + 1):(2*p)]),
                      rbind(ff_res[(p + 1):(2*p),1:p], ff_res[1:p,1:p])) * -grid_info[['phase_factor_mat']]
  #x_vals <- grid_info[['x_vals']] - (grid_info[['x_vals']][length(grid_info[['x_vals']])] - grid_info[['x_vals']][1]) / 2
  cbind(grid_info[['x_vals_eg']], 'val' = (length(grid_info[['x_vals']]))^(2) *
          as.double(Re(ff_res_adj) * norm_constant(nu1, nu2, a1, a2, d = d)) * 2 / pi / 
          grid_info[['x_max']]^2 * pi^4)
}


#### Finally, a more rudimentary (bad) way of computing things
# approximate integrals

# d=1
spec_dens <- function(r, h, nu1, nu2, a1 = 1, a2 = 1, re_z, im_z) {
  val <- exp(complex(imaginary = h * r)) *
    complex(real = a1, imaginary = r)^(-nu1 - 1/2) *
    complex(real = a2, imaginary = -r)^(-nu2 - 1/2) * 
    complex(real = re_z, imaginary = sign(r) * im_z)
  Re(val)
}

# d=2 
spatial_integrate_d2 <- function(h, a1, a2, nu1, nu2, Mu, Psi= Psi, approx_grid,
                                 d = 2) {
  theta_x <-  approx_grid[['theta_x']]
  theta_y <-  approx_grid[['theta_y']]
  r <-  approx_grid[['r']]
  Psi_val <- Psi(theta_x = theta_x, theta_y = theta_y)
  complex_r <- complex(imaginary = r)
  values <- exp(complex_r*(h[1] * theta_x + h[2] * theta_y)) *
    Mu(theta_x, theta_y, 1, 2) *
    (a1 + Psi_val * complex_r)^(-nu1 - d/2) *
    (a2 - Psi_val * complex_r)^(-nu2 - d/2) * r^(d - 1) 
  Re(sum(values*approx_grid[['angle_lag']]*approx_grid[['r_lag']]))
}
