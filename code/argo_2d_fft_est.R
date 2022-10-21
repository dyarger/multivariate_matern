
library(dplyr)
dat300 <- read.table('bolin_code/article_code/Application/Argo/data300.txt',header=T)
dat1500 <- read.table('bolin_code/article_code/Application/Argo/data1500.txt',header=T)
dat300 <- dat300[dat300$year==2015,]
dat1500 <- dat1500[dat1500$year==2015,]

data_all <- rbind(cbind(dat300, level = 300), cbind(dat1500, level = 1500))

ggplot(data =data_all,
       aes(x = lon, y = lat, color = Y))+
  facet_grid(level~year) + 
  geom_point() + 
  scale_color_gradient2()
library(fields)
library(Rcpp)
library(tidyverse)
library(fftw)
library(fftwtools)

norm_constant <- function(nu_1, nu_2, a_1 = 1, a_2 = 1, d = 2) {
  (a_1)^(nu_1) * (a_2)^(nu_2) *
    sqrt(gamma(nu_1 + d/2)) * sqrt(gamma(nu_2 + d/2))/pi^(d/2)/sqrt(gamma(nu_1)*gamma(nu_2))
}
# given a grid and model parameters, computes value of Matern cross-covariance on grid
fft_2d <- function(grid_info, nu1 = .5, nu2 = .5, a1 = 1, a2 = 1, 
                   Psi, 
                   Delta, d = 2) {
  freq_grid <- grid_info[['freq_grid']]
  freq_points <- grid_info[['freq_points']]
  delta_t <- grid_info[['delta_t']]
  n_points <- grid_info[['n_points']]
  x_vals <- grid_info[['x_vals']]
  x_max <- grid_info[['x_max']]
  x_vals_eg <- grid_info[['x_vals_eg']]
  phase_factor_mat <- grid_info[['phase_factor_mat']]
  # https://stackoverflow.com/questions/24077913/discretized-continuous-fourier-transform-with-numpy
  tv <- complex(real = a1, imaginary = Psi(freq_grid[['theta']])*freq_grid[['r']])^(-nu1 - d/2) *
    complex(real = a2, imaginary = -Psi(freq_grid[['theta']])*freq_grid[['r']])^(-nu2 - d/2) * 
    Delta(freq_grid[['theta']])
  tv_mat <- matrix(nrow = length(freq_points), ncol = length(freq_points), tv)
  #ff_res <- armafft(tv_mat)
  ff_res <- fftwtools::fftw_c2c_2d(data = tv_mat, inverse = 1)/n_points^2
  p <- ncol(ff_res)/2
  ff_res_adj <- cbind(rbind(ff_res[(p + 1):(2*p),(p + 1):(2*p)], ff_res[1:p,(p + 1):(2*p)]),
                      rbind(ff_res[(p + 1):(2*p),1:p], ff_res[1:p,1:p])) * -phase_factor_mat
  x_vals <- x_vals - (x_vals[length(x_vals)] - x_vals[1]) / 2
  cbind(x_vals_eg, 'val' = (length(x_vals))^(2) *
          as.double(Re(ff_res_adj) * norm_constant(nu1, nu2, a1, a2)) * 2 / pi / x_max^2 / 0.01026171)
}
Delta <- function(x) {
  1
}
Psi <- function(x) {
  sign(x)
}

# define Fourier transform grid 
create_grid_info <- function(n_points, x_max) {
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
  #x_vals <- x_vals - (x_vals[length(x_vals)] - x_vals[1]) / 2
  x_vals <- x_vals - x_max - abs(abs(x_vals[length(x_vals)] - 2*x_max) - abs(x_vals[1]))
  x_vals_eg <- expand.grid(x_vals, x_vals)
  list('freq_grid' = freq_grid, 'freq_points' = freq_points,
       'delta_t' = delta_t, 'n_points' = n_points, 
       'x_vals' = x_vals, 'x_max' = x_max,
       'phase_factor_mat' = phase_factor_mat,
       'x_vals_eg' = x_vals_eg)
}


n1 <- nrow(dat300)
n2 <- nrow(dat1500)

dist1 <- fields::rdist.earth(x1 = dat300[, c('lon', 'lat')],
                             x2 = dat300[, c('lon', 'lat')], miles = F)
diag(dist1) <- 0
dist2 <- fields::rdist.earth(x1 = dat1500[, c('lon', 'lat')],
                             x2 = dat1500[, c('lon', 'lat')], miles = F)
diag(dist2) <- 0
dist12 <- fields::rdist.earth(x1 = dat300[, c('lon', 'lat')],
                             x2 = dat1500[, c('lon', 'lat')], miles = F)
diag(dist12) <- 0
dist_all <- rbind(cbind(dist1, dist12), cbind(t(dist12), dist2))

response <- c(dat300$Y, dat1500$Y)

# create long/lat distance matrices
dist_one <- dist
dist_tens <- array(NA, c(dim(dist_all), 2))
reference_location <- colMeans(rbind(dat300[, c('lon', 'lat')], 
                                     dat1500[, c('lon', 'lat')]))

for (i in 1:dim(dist_tens)[1]) {
  for (j in 1:dim(dist_tens)[1]) {
    dist_long_prelim <- fields::rdist.earth.vec(cbind(data_all$lon[i],
                                                      reference_location['lat']), 
                                                cbind(data_all$lon[j],
                                                      reference_location['lat']),
                                                miles = F)
    if (data_all$lon[i] < data_all$lon[j]) {
      dist_long_prelim <- -dist_long_prelim
    }
    
    dist_lat_prelim <- fields::rdist.earth.vec(cbind(reference_location['lon'],
                                                     data_all$lat[i]), 
                                               cbind(reference_location['lon'],
                                                     data_all$lat[j]),
                                               miles = F)
    if (data_all$lat[i] < data_all$lat[j]) {
      dist_lat_prelim <- -dist_lat_prelim
    }
    
    dist_tens[i,j,1] <- dist_long_prelim
    dist_tens[i,j,2] <- dist_lat_prelim
  }
}
dist_tens_1 <- dist_tens[1:nrow(dat300), 1:nrow(dat300),]
dist_tens_2 <- dist_tens[-(1:nrow(dat300)), -(1:nrow(dat300)),]
dist_tens_12 <- dist_tens[1:nrow(dat300), -(1:nrow(dat300)),]
dist_tens_mat_all <- cbind(as.vector(dist_tens[,,1]), as.vector(dist_tens[,,2]))
dist_tens_mat_1 <- cbind(as.vector(dist_tens_1[,,1]), as.vector(dist_tens_1[,,2]))
dist_tens_mat_2 <- cbind(as.vector(dist_tens_2[,,1]), as.vector(dist_tens_2[,,2]))
dist_tens_mat_12 <- cbind(as.vector(dist_tens_12[,,1]), as.vector(dist_tens_12[,,2]))

# given distances, compute fft on grid, then interpolate onto distances
construct_matrix <- function(nu1, nu2, a1, a2, 
                             grid_info,
                             Psi, Delta, dist_tens_mat, n1, n2) {
  C_val <- fft_2d(grid_info = grid_info, nu1 = nu1, nu2 = nu2, 
                  a1 = a1, a2 = a2,  Psi = Psi, Delta = Delta)
  
  tmat <- matrix(C_val[['val']], nrow = sqrt(nrow(C_val)))
  long1 <- grid_info[['x_vals']]
  lat1 <- grid_info[['x_vals']]
  
  interp_points <- fields::interp.surface(obj = list(x = long1, y = lat1, z = tmat),
                                          loc = dist_tens_mat)
  matrix(interp_points, nrow = n1, ncol = n2)
}

# construct bivariate Matern covariance matrix
construct_entire_matrix <- function(nu1, nu2, a1, a2, 
                                    grid_info,
                                    Psi_list, Delta_list, dist_tens_mat_1, 
                                    dist_tens_mat_2, dist_tens_mat_12, nugget1, nugget2,
                                    n1, n2) {
  C1 <- construct_matrix(nu1 = nu1, nu2 = nu1, a1 = a1, a2 = a1,
                         Psi = Psi_list[[1]], Delta = Delta_list[[1]], 
                         grid_info = grid_info, dist_tens_mat_1, 
                         n1 = n1, n2 = n1)
  diag(C1) <- diag(C1) + nugget1
  C12 <- construct_matrix(nu1 = nu1, nu2 = nu2, a1 = a1, a2 = a2,
                          Psi = Psi_list[[3]], Delta = Delta_list[[3]], 
                          grid_info = grid_info, dist_tens_mat_12, 
                          n1 = n1, n2 = n2)
  C2 <- construct_matrix(nu1 = nu2, nu2 = nu2, a1 = a2, a2 = a2,
                         Psi = Psi_list[[2]], Delta = Delta_list[[2]],
                         grid_info = grid_info, dist_tens_mat_2, 
                         n1 = n2, n2 = n2)
  diag(C2) <- diag(C2) + nugget2
  rbind(cbind(C1, C12), cbind(t(C12), C2))
}
grid_info <- create_grid_info(2^10, x_max = 3000)

nu1_init <- nu2_init <- .5
nu_lower <- 0.001; nu_higher = 7
a1_init <- a2_init <- 10^-3
a_lower <- 10^-6
Sigma11_init <- var(dat300$Y)
Sigma22_init <- var(dat1500$Y)
cor_re_init <- 0
cor_im_init <- 0
nugget1_init <- var(dat1500$Y)/25
nugget2_init <- var(dat1500$Y)/25
var_lower <- 10^-9

# log likelihood as function of parameters
ll_fun <- function(par, dist_tens_mat_1, dist_tens_mat_2, dist_tens_mat_12, response, grid_info,
                   n1, n2) {
  print(exp(par[1:8]))
  nu1 <- exp(par[1]); nu2 <- exp(par[2])
  a1 <- exp(par[3]); a2 <- exp(par[4])
  nugget1 <- exp(par[5]); nugget2 <- exp(par[6])
  Sigma11 <- exp(par[7]); Sigma22 <- exp(par[8])
  Sigma12 <- par[9]*sqrt(Sigma11)*sqrt(Sigma22)
  cov_mat <- construct_entire_matrix(nu1 = nu1, nu2 = nu2, a1 = a1, a2 = a2, 
                                     Delta_list = list(function(x) Sigma11, function(x) Sigma22, 
                                                       function(x) Sigma12),
                                     Psi_list = replicate(3, {function(x) {
                                       sign(x)
                                     }}),
                                     dist_tens_mat_1 = dist_tens_mat_1,
                                     dist_tens_mat_2 = dist_tens_mat_2, 
                                     dist_tens_mat_12 = dist_tens_mat_12, grid_info = grid_info,
                                     nugget1 = nugget1, nugget2 = nugget2,
                                     n1 = n1, n2 = n2)
  chol_mat <- base::chol(cov_mat)
  ll_val <- -nrow(cov_mat)/2 * log(2 * pi) - 1/2 * sum(backsolve(chol_mat, response, transpose = T)^2) - 
    sum(log(diag(chol_mat)))
  print(ll_val)
  - ll_val
}


#### with psi fixed at sign(theta)
test_optim <- optim(
  par = c(log(c(nu1_init, nu2_init, a1_init, a2_init, nugget1_init, nugget2_init, 
                Sigma11_init, Sigma22_init)), cor_re_init),
  fn = ll_fun, dist_tens_mat_1 = dist_tens_mat_1,
  dist_tens_mat_2 = dist_tens_mat_2, 
  dist_tens_mat_12 = dist_tens_mat_12, response = response, n1 = n1, n2 = n2,
  lower = c(log(c(nu_lower, nu_lower, a_lower, a_lower, var_lower, var_lower, var_lower, var_lower)), -1),
  upper = c(log(c(nu_higher, nu_higher, NA, NA, NA, NA, NA, NA)), 1),
  method = 'L-BFGS-B',
  grid_info = grid_info,
  control = list(parscale = c(rep(1, 8), .1))
)

test_optim
par <- test_optim$par
(nu1 <- exp(par[1]))
(nu2 <- exp(par[2]))
(a1 <- exp(par[3]))
(a2 <- exp(par[4]))

(nugget1 <-  exp(par[5]))
(nugget2 <-  exp(par[6]))
(Sigma11 <- exp(par[7]))
(Sigma22 <- exp(par[8]))
(Sigma12re <- par[9]*sqrt(Sigma11)*sqrt(Sigma22))
par[9]

df_fix <- fft_2d(nu1 = nu1, nu2 = nu2, a1 = a1, a2 = a2, 
                 d = 2,
                 Delta = function(x) Sigma12re, Psi = function(x) {
                   sign(x)
                 }, grid_info = grid_info)

ggplot(data = df_fix %>%
         filter(abs(Var1) < 750, abs(Var2) < 750), aes(x = Var2, y = Var1, fill = val)) +
  geom_raster() + 
  scale_fill_gradient2() +
  labs(x = 'Dimension 1', y = 'Dimension 2',
       fill = 'Cross-\ncovariance') +
  coord_equal() + 
  theme(legend.position = 'left',legend.key.height = unit(.8, "cm")) 


ll_fun_psi_real <- function(par, dist_tens_mat_1, dist_tens_mat_2, dist_tens_mat_12,
                            response, grid_info, n1, n2) {
  print(exp(par[1:8]))
  nu1 <- exp(par[1]); nu2 <- exp(par[2])
  a1 <- exp(par[3]); a2 <- exp(par[4])
  nugget1 <- exp(par[5]); nugget2 <- exp(par[6])
  Sigma11 <- exp(par[7]); Sigma22 <- exp(par[8])
  Sigma12 <- par[9]*sqrt(Sigma11)*sqrt(Sigma22)
  theta_star <- atan2(par[11], par[10])
  Psi_fun <- function(theta) {
    if (theta_star < 0) {
      ifelse(theta > theta_star & theta < theta_star + pi, 1, -1)
    } else {
      ifelse(theta > theta_star | theta < theta_star - pi, 1, -1)
    }
  }
  cov_mat <- construct_entire_matrix(nu1 = nu1, nu2 = nu2, a1 = a1, a2 = a2, 
                                     Delta_list = list(function(x) Sigma11, function(x) Sigma22, 
                                                       function(x) Sigma12),
                                     Psi_list = replicate(3, Psi_fun),
                                     dist_tens_mat_1 = dist_tens_mat_1, 
                                     dist_tens_mat_2 = dist_tens_mat_2, 
                                     dist_tens_mat_12 = dist_tens_mat_12, grid_info = grid_info,
                                     nugget1 = nugget1, nugget2 = nugget2,
                                     n1 = n1, n2 = n2)
  chol_mat <- base::chol(cov_mat)
  ll_val <- -nrow(cov_mat)/2 * log(2 * pi) - 1/2 * sum(backsolve(chol_mat, response, transpose = T)^2) - 
    sum(log(diag(chol_mat)))
  print(ll_val)
  - ll_val
}
test_optim_real <- optim(
  par = c(test_optim$par, 0.1, 0.1),
  fn = ll_fun_psi_real, dist_tens_mat_1 = dist_tens_mat_1, 
  dist_tens_mat_12 = dist_tens_mat_12, n1 = n1, n2 = n2,
  dist_tens_mat_2 = dist_tens_mat_2,response = response, 
  lower = c(log(c(nu_lower, nu_lower, a_lower, a_lower, var_lower, var_lower, var_lower, var_lower)), -1, NA, NA),
  upper = c(log(c(nu_higher, nu_higher, NA, NA, NA, NA, NA, NA)), 1, NA, NA),
  method = 'L-BFGS-B',
  grid_info = grid_info,
  control = list(parscale = c(rep(1, 8), .1, .1, .1))
)
test_optim_real
par <- test_optim_real$par
(nu1 <- exp(par[1]))
(nu2 <- exp(par[2]))
(a1 <- exp(par[3]))
(a2 <- exp(par[4]))

(nugget1 <-  exp(par[5]))
(nugget2 <-  exp(par[6]))
(Sigma11 <- exp(par[7]))
(Sigma22 <- exp(par[8]))
(Sigma12re <- par[9]*sqrt(Sigma11)*sqrt(Sigma22))
par[9]
atan2(par[11],par[10])

theta_star <- atan2(par[11], par[10])
Psi_fun <- function(theta) {
  if (theta_star < 0) {
    ifelse(theta > theta_star & theta < theta_star + pi, 1, -1)
  } else {
    ifelse(theta > theta_star | theta < theta_star - pi, 1, -1)
  }
}

df_re <- fft_2d(nu1 = nu1, nu2 = nu2, a1 = a1, a2 = a2, 
                d = 2,
                Delta = function(x) Sigma12re, Psi = Psi_fun, grid_info = grid_info)

ggplot(data = df_re %>%
         filter(abs(Var1) < 750, abs(Var2) < 750), aes(x = Var2, y = Var1, fill = val)) +
  geom_raster() + 
  scale_fill_gradient2() +
  labs(x = 'Dimension 1', y = 'Dimension 2',
       fill = 'Cross-\ncovariance') +
  coord_equal() + 
  theme(legend.position = 'left',legend.key.height = unit(.8, "cm")) 





# with an imaginary entry
ll_fun_psi_im <- function(par, response, grid_info,
                          dist_tens_mat_1 = dist_tens_mat_1, 
                          dist_tens_mat_12 = dist_tens_mat_12, n1 = n1, n2 = n2,
                          dist_tens_mat_2 = dist_tens_mat_2) {
  print(exp(par[1:8]))
  print((par[9:10]))
  nu1 <- exp(par[1]); nu2 <- exp(par[2])
  a1 <- exp(par[3]); a2 <- exp(par[4])
  nugget1 <- exp(par[5]); nugget2 <- exp(par[6])
  Sigma11 <- exp(par[7]); Sigma22 <- exp(par[8])
  Sigma12 <- par[9]*sqrt(Sigma11)*sqrt(Sigma22)
  Sigma12im <- par[10]*sqrt(Sigma11)*sqrt(Sigma22)
  theta_star <- atan2(par[12], par[11])
  Psi_fun <- function(theta) {
    if (theta_star < 0) {
      ifelse(theta > theta_star & theta < theta_star + pi, 1, -1)
    } else {
      ifelse(theta > theta_star | theta < theta_star - pi, 1, -1)
    }
  }
  theta_star2 <- atan2(par[14], par[13])
  Psi_fun2 <- function(theta) {
    if (theta_star < 0) {
      ifelse(theta > theta_star2 & theta < theta_star2 + pi, 1, -1)
    } else {
      ifelse(theta > theta_star2 | theta < theta_star2 - pi, 1, -1)
    }
  }
  cov_mat <- construct_entire_matrix(nu1 = nu1, nu2 = nu2, a1 = a1, a2 = a2, 
                                     Delta_list = list(function(x) Sigma11, 
                                                       function(x) Sigma22, 
                                                       function(x) complex(real = Sigma12,
                                                                           imaginary = Sigma12im*Psi_fun2(x))),
                                     Psi_list = replicate(3, Psi_fun),
                                     grid_info = grid_info,
                                     nugget1 = nugget1, nugget2 = nugget2,
                                     dist_tens_mat_1 = dist_tens_mat_1, 
                                     dist_tens_mat_12 = dist_tens_mat_12, n1 = n1, n2 = n2,
                                     dist_tens_mat_2 = dist_tens_mat_2)
  chol_mat <- base::chol(cov_mat)
  ll_val <- -nrow(cov_mat)/2 * log(2 * pi) - 1/2 * sum(backsolve(chol_mat, response, transpose = T)^2) - 
    sum(log(diag(chol_mat)))
  print(ll_val)
  - ll_val
}

test_optim_im <- optim(
  par = c(log(c(nu1_init, nu2_init, a1_init, a2_init, nugget1_init, nugget2_init, 
                Sigma11_init, Sigma22_init)), cor_re_init, .2, .1, .1, .1, .1),
  #par = c(test_optim_real$par[1:9], 0, test_optim_real$par[10:11]),
  fn = ll_fun_psi_im,  dist_tens_mat_1 = dist_tens_mat_1, 
  dist_tens_mat_12 = dist_tens_mat_12, n1 = n1, n2 = n2,
  dist_tens_mat_2 = dist_tens_mat_2, response = response, 
  lower = c(log(c(nu_lower, nu_lower, a_lower, a_lower, var_lower, var_lower, var_lower, var_lower)), -1, -1, NA, NA, NA, NA),
  upper = c(log(c(nu_higher, nu_higher, NA, NA, NA, NA, NA, NA)), 1, 1, NA, NA, NA, NA),
  method = 'L-BFGS-B',
  grid_info = grid_info,
  control = list(parscale = c(rep(1, 8), .1, .1, .1, .1, .1, .1))
)
test_optim_im
par <- test_optim_im$par
(nu1 <- exp(par[1]))
(nu2 <- exp(par[2]))
(a1 <- exp(par[3]))
(a2 <- exp(par[4]))

(nugget1 <-  exp(par[5]))
(nugget2 <-  exp(par[6]))
(Sigma11 <- exp(par[7]))
(Sigma22 <- exp(par[8]))
(Sigma12re <- par[9]*sqrt(Sigma11)*sqrt(Sigma22))
(Sigma12im <- par[10]*sqrt(Sigma11)*sqrt(Sigma22))
par[9]
par[10]
atan2(par[12],par[11])

theta_star <- atan2(par[12], par[11])
theta_star2 <- atan2(par[14], par[13])
Psi_fun <- function(theta) {
  if (theta_star < 0) {
    ifelse(theta > theta_star & theta < theta_star + pi, 1, -1)
  } else {
    ifelse(theta > theta_star | theta < theta_star - pi, 1, -1)
  }
}
Psi_fun2 <- function(theta) {
  if (theta_star < 0) {
    ifelse(theta > theta_star2 & theta < theta_star2 + pi, 1, -1)
  } else {
    ifelse(theta > theta_star2 | theta < theta_star2 - pi, 1, -1)
  }
}

grid_info_fine <- create_grid_info(2^11, x_max = 3000)


df_im <- fft_2d(nu1 = nu1, nu2 = nu2, a1 = a1, a2 = a2, 
                d = 2,
                Delta = function(x) complex(real = Sigma12re, imaginary = Sigma12im*Psi_fun2(x)),
                Psi = Psi_fun, grid_info = grid_info_fine)

ggplot(data = df_im %>%
         filter(abs(Var1) < 300, abs(Var2) < 300), aes(x = Var2, y = Var1, fill = val)) +
  geom_raster() + 
  scale_fill_gradientn(colors = rev(rainbow(10))) +
  labs(x = 'Dimension 1', y = 'Dimension 2',
       fill = 'Cross-\ncovariance') +
  coord_equal() + 
  theme(legend.position = 'left', legend.key.height = unit(.8, "cm")) 


# single correlation function
# independent Matern
ll_fun_single <- function(par, response, grid_info,
                          dist_tens_mat_1 = dist_tens_mat_1, 
                          dist_tens_mat_12 = dist_tens_mat_12, n1 = n1, n2 = n2,
                          dist_tens_mat_2 = dist_tens_mat_2) {
  #print(exp(par[1:7]))
  nu1 <- exp(par[1]); a1 <- exp(par[2])
  nugget1 <- exp(par[3]); nugget2 <- exp(par[4])
  Sigma11 <- exp(par[5]); Sigma22 <- exp(par[6])
  Sigma12 <- par[7]*sqrt(Sigma11)*sqrt(Sigma22)
  cov_mat <- construct_entire_matrix(nu1 = nu1, nu2 = nu1, a1 = a1, a2 = a1, 
                                     Delta_list = list(function(x) Sigma11, 
                                                       function(x) Sigma22, 
                                                       function(x) Sigma12),
                                     Psi_list = replicate(3, function(x) sign(x)),
                                     dist_tens_mat_1 = dist_tens_mat_1, dist_tens_mat_2 = dist_tens_mat_2,
                                     dist_tens_mat_12 = dist_tens_mat_12, n1 = n1, 
                                     n2 = n2, grid_info = grid_info,
                                     nugget1 = nugget1, nugget2 = nugget2)
  chol_mat <- base::chol(cov_mat)
  ll_val <- -nrow(cov_mat)/2 * log(2 * pi) - 1/2 * sum(backsolve(chol_mat, response, transpose = T)^2) - 
    sum(log(diag(chol_mat)))
  #print(ll_val)
  
  - ll_val
}

test_optim_single <- optim(
  par = c(log(c(nu1_init, a1_init, nugget1_init,nugget2_init ,
                Sigma11_init, Sigma22_init)), cor_re_init),
  fn = ll_fun_single,  dist_tens_mat_1 = dist_tens_mat_1, dist_tens_mat_2 = dist_tens_mat_2,
  dist_tens_mat_12 = dist_tens_mat_12, n1 = n1, 
  n2 = n2, response = response, 
  lower = c(log(c(nu_lower, a_lower, var_lower, var_lower, var_lower, var_lower)), -1),
  upper = c(log(c(nu_higher, NA, NA, NA, NA, NA)), 1),
  method = 'L-BFGS-B',
  grid_info = grid_info,
  control = list(parscale = c(rep(1, 6), .1))
)
test_optim_single
par <- test_optim_single$par
(nu1 <- exp(par[1]))
(a1 <- exp(par[2]))

(nugget1 <-  exp(par[3]))
(nugget2 <-  exp(par[4]))
(Sigma11 <- exp(par[5]))
(Sigma22 <- exp(par[6]))
par[7]

df_single <- fft_2d(nu1 = nu1, nu2 = nu1, a1 = a1, a2 = a1, Delta = 
                      function(x) par[7] * sqrt(Sigma11) * sqrt(Sigma22), Psi = function(x) sign(x), grid_info = grid_info)

# bivariate Matern
construct_matrix_mm <- function(nu1, nu2, nu12, a1, a2, a12,
                                grid_info,
                                Psi_list, Delta_list, dist_tens_mat_1, dist_tens_mat_2,
                                dist_tens_mat_12, n1, n2, nugget1, nugget2) {
  C1 <- construct_matrix(nu1 = nu1, nu2 = nu1, a1 = a1, a2 = a1,
                         Psi = Psi_list[[1]], Delta = Delta_list[[1]], 
                         grid_info = grid_info, dist_tens_mat_1, n1 = n1, n2 = n1)
  diag(C1) <- diag(C1) + nugget1
  C12 <- construct_matrix(nu1 = nu12, nu2 = nu12, a1 = a12, a2 = a12,
                          Psi = Psi_list[[3]], Delta = Delta_list[[3]], 
                          grid_info = grid_info, dist_tens_mat_12, 
                          n1 = n1, n2 = n2)
  C2 <- construct_matrix(nu1 = nu2, nu2 = nu2, a1 = a2, a2 = a2,
                         Psi = Psi_list[[2]], Delta = Delta_list[[2]],
                         grid_info = grid_info, dist_tens_mat_2, 
                         n1 = n2, n2 = n2)
  diag(C2) <- diag(C2) + nugget2
  rbind(cbind(C1, C12), cbind(t(C12), C2))
}

# log likelihood as function of parameters
ll_fun_mm <- function(par, response, grid_info,
                      dist_tens_mat_1, dist_tens_mat_2,
                      dist_tens_mat_12, n1, n2) {
  #print(exp(par[1:10]))
  nu1 <- exp(par[1]); nu2 <- exp(par[2]); nu12 <- exp(par[3])
  a1 <- exp(par[4]); a2 <- exp(par[5]); a12 <- exp(par[6])
  nugget1 <- exp(par[7]); nugget2 <- exp(par[8])
  Sigma11 <- exp(par[9]); Sigma22 <- exp(par[10])
  Sigma12 <- par[11]*sqrt(Sigma11)*sqrt(Sigma22)
  cov_mat <- construct_matrix_mm(nu1 = nu1, nu2 = nu2, a1 = a1, a2 = a2, 
                                 nu12 = nu12, a12 = a12,
                                 Delta_list = list(function(x) Sigma11, function(x) Sigma22, 
                                                   function(x) Sigma12),
                                 Psi_list = replicate(3, {function(x) {
                                   sign(x)
                                 }}),
                                 dist_tens_mat_1 = dist_tens_mat_1, dist_tens_mat_2 = dist_tens_mat_2,
                                 dist_tens_mat_12 = dist_tens_mat_12, n1 = n1, 
                                 n2 = n2,  grid_info = grid_info,
                                 nugget1 = nugget1, nugget2 = nugget2)
  if (min(eigen(cov_mat)$values) < .Machine$double.eps) {
    #print(NA)
    return(8000)
  }
  chol_mat <- base::chol(cov_mat)
  ll_val <- -nrow(cov_mat)/2 * log(2 * pi) - 1/2 * sum(backsolve(chol_mat, response, transpose = T)^2) - 
    sum(log(diag(chol_mat)))
  #print(ll_val)
  
  - ll_val
}

test_optim_mm <- optim(
  par = c(log(c(nu1_init, nu2_init, (nu1_init + nu2_init)/2, 
                a1_init, a2_init,  (a1_init + a1_init)/2,
                nugget1_init, nugget2_init,
                Sigma11_init, Sigma22_init)), cor_re_init),
  fn = ll_fun_mm, dist_tens_mat_1 = dist_tens_mat_1, dist_tens_mat_2 = dist_tens_mat_2,
  dist_tens_mat_12 = dist_tens_mat_12, n1 = n1, 
  n2 = n2, response = response, 
  lower = c(log(c(nu_lower,nu_lower, nu_lower,  a_lower,a_lower, a_lower,
                  var_lower, var_lower, var_lower, var_lower)), -1),
  upper = c(log(c(nu_higher, nu_higher, nu_higher, 1, 1, 1, NA, NA, NA, NA)), 1),
  method = 'L-BFGS-B',
  grid_info = grid_info,
  control = list(parscale = c(rep(1, 10), .1))
)

test_optim_mm
par <- test_optim_mm$par
(nu1 <- exp(par[1]))
(nu2 <- exp(par[2]))
(nu12 <- exp(par[3]))
(a1 <- exp(par[4]))
(a2 <- exp(par[5]))
(a12 <- exp(par[6]))

(nugget1 <-  exp(par[7]))
(nugget2 <-  exp(par[8]))
(Sigma11 <- exp(par[9]))
(Sigma22 <- exp(par[10]))
par[11]

df_mm <- fft_2d(nu1 = nu12, nu2 = nu12, a1 = a12, a2 = a12, Delta = 
                  function(x) par[11] * sqrt(Sigma11) * sqrt(Sigma22), 
                Psi = function(x) sign(x), grid_info = grid_info)


# independent Matern
construct_entire_matrix_independent <- function(nu1, nu2, a1, a2, 
                                                grid_info,
                                                Psi_list, Delta_list, dist_tens_mat_1, dist_tens_mat_2,
                                                dist_tens_mat_12, n1, n2, nugget1, nugget2) {
  C1 <- construct_matrix(nu1 = nu1, nu2 = nu1, a1 = a1, a2 = a1,
                         Psi = Psi_list[[1]], Delta = Delta_list[[1]], 
                         grid_info = grid_info, dist_tens_mat_1, n1 = n1, n2 = n1)
  diag(C1) <- diag(C1) + nugget1
  C2 <- construct_matrix(nu1 = nu2, nu2 = nu2, a1 = a2, a2 = a2,
                         Psi = Psi_list[[2]], Delta = Delta_list[[2]],
                         grid_info = grid_info, dist_tens_mat_2, n1 = n2, n2 = n2)
  diag(C2) <- diag(C2) + nugget2
  as.matrix(Matrix::bdiag(C1, C2))
}

ll_fun_ind <- function(par, response, grid_info, dist_tens_mat_1, dist_tens_mat_2,
                       dist_tens_mat_12, n1, n2) {
  #print(exp(par[1:8]))
  nu1 <- exp(par[1]); nu2 <- exp(par[2])
  a1 <- exp(par[3]); a2 <- exp(par[4])
  nugget1 <- exp(par[5]); nugget2 <- exp(par[6])
  Sigma11 <- exp(par[7]); Sigma22 <- exp(par[8])
  cov_mat <- construct_entire_matrix_independent(nu1 = nu1, nu2 = nu2, a1 = a1, a2 = a2, 
                                                 Delta_list = list(function(x) Sigma11, 
                                                                   function(x) Sigma22, 
                                                                   function(x) 0),
                                                 Psi_list = replicate(3, function(x) sign(x)),
                                                 dist_tens_mat_1=dist_tens_mat_1, dist_tens_mat_2=dist_tens_mat_2,
                                                 dist_tens_mat_12=dist_tens_mat_12, n1 =n1, n2=n2, grid_info = grid_info,
                                                 nugget1 = nugget1, nugget2 = nugget2)
  chol_mat <- base::chol(cov_mat)
  ll_val <- -nrow(cov_mat)/2 * log(2 * pi) - 1/2 * sum(backsolve(chol_mat, response, transpose = T)^2) - 
    sum(log(diag(chol_mat)))
  #print(ll_val)
  - ll_val
}

test_optim_ind <- optim(
  par = c(log(c(nu1_init, nu2_init, a1_init, a2_init, nugget1_init, nugget2_init, 
                Sigma11_init, Sigma22_init))),
  fn = ll_fun_ind, dist_tens_mat_2=dist_tens_mat_2, dist_tens_mat_1=dist_tens_mat_1,
  dist_tens_mat_12=dist_tens_mat_12, n1 =n1, n2=n2, response = response, 
  lower = c(log(c(nu_lower, nu_lower, a_lower, a_lower, var_lower, var_lower, var_lower, var_lower))),
  upper = c(log(c(nu_higher, nu_higher, NA, NA, NA, NA, NA, NA))),
  method = 'L-BFGS-B',
  grid_info = grid_info,
  control = list(parscale = rep(1, 8))
)
test_optim_ind
par <- test_optim_ind$par
(nu1 <- exp(par[1]))
(nu2 <- exp(par[2]))
(a1 <- exp(par[3]))
(a2 <- exp(par[4]))

(nugget1 <-  exp(par[5]))
(nugget2 <-  exp(par[6]))
(Sigma11 <- exp(par[7]))
(Sigma22 <- exp(par[8]))




df_all <- rbind(cbind(df_fix, type = 'theta_star_fixed'),
                cbind(df_re, type = 'real'),
                cbind(df_im, type = 'imaginary'),
                #cbind(df_im_only, type = 'imaginary_only'),
                cbind(df_mm, type = 'multi_matern'),
                cbind(df_single, type = 'single')
                )

# df_params <- data.frame(type = c('theta_star_fixed', 'real', 'imaginary', 'multi_matern', 'single'),
#                         var1 = exp(c(test_optim$par[7], test_optim_real$par[7],
#                                      test_optim_im$par[7], test_optim_mm$par[9],
#                                      test_optim_single$par[5])),
#                         var2 = exp(c(test_optim$par[8], test_optim_real$par[8],
#                                      test_optim_im$par[8], test_optim_mm$par[10],
#                                      test_optim_single$par[6])))

save(test_optim, test_optim_im, test_optim_single, test_optim_mm, 
     test_optim_real,dist_tens_mat_all,
     dist_tens_mat_1,
     dist_tens_mat_2,
     dist_tens_mat_12,
     df_all,
     file = 'results/argo_parameters.RData')

#load('pres_temp_spectral_results_11.RData')



labels <- data.frame(type = c('real', 'imaginary', 'multi_matern', 'single', 'theta_star_fixed'),
                     label = factor(c('TPMM w/ real directional measure',
                                      'TPMM w/ complex directional measure',
                                      'MM of Gnieting et al. (2010)',
                                      'Single covariance function',
                                      'TPMM w/ phi fixed'),
                                    levels = c('Single covariance function',
                                               'MM of Gnieting et al. (2010)',
                                               'TPMM w/ phi fixed', 
                                               'TPMM w/ real directional measure',
                                               'TPMM w/ complex directional measure')))

ggplot(data = df_all %>%
         filter(abs(Var1) < 300, abs(Var2) < 300, type != 'imaginary_only') %>%
         left_join(labels), aes(x = Var2, y = Var1, fill = val)) +
  geom_raster() + 
  facet_wrap(~label) + 
  scale_fill_gradientn(colors = rev(rainbow(10))) +
  #scale_fill_gradient2() +
  labs(x = 'Zonal distance (km)', y = 'Meridional distance (km)',
       fill = 'Cross-\ncovariance') +
  coord_equal() + 
  theme(legend.position = 'left',legend.key.height = unit(.8, "cm")) 
ggsave('images/argo_cov_fun_comparison_data.png', height = 6, width = 9)