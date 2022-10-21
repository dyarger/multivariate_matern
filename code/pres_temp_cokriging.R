load('results/pres_temp_actual_spectral_results_10_final.RData')
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

grid_info <- create_grid_info(n_points = 2^10, x_max = 10)


# test it out
df <- fft_2d(nu1 = 1.5, nu2 = 1.5, a1 = 1, a2 = 1, Psi = Psi, Delta = Delta, grid_info = grid_info)
library(fields)
var1_vals <- df[['Var1']][df[['Var2']] == min(abs(df[['Var2']]))]
plot(var1_vals, Re(df[['val']][df[['Var2']] == min(abs(df[['Var2']]))]), type = 'l')
lines(var1_vals, col = 2, Matern(abs(var1_vals), nu = 1.5))
plot(var1_vals, Re(df[['val']][df[['Var2']] == min(abs(df[['Var2']]))]) - 
       Matern(abs(var1_vals), nu = 1.5), type = 'l')
df_mat <- matrix(df[['val']], nrow = sqrt(length(df[['val']])))
image.plot(df_mat)


df <- fft_2d(nu1 = 1.5, nu2 = .2, a1 = 1, a2 = 1, Delta = Delta, Psi = Psi, grid_info = grid_info)
var1_vals <- df$Var1[df[['Var2']] == min(abs(df[['Var2']]))]
plot(var1_vals, Re(df[['val']][df[['Var2']] == min(abs(df[['Var2']]))]), type = 'l')
df_mat <- matrix(df[['val']], nrow = sqrt(length(df[['val']])))
image.plot(df_mat)


Delta <- function(x) {
  .97 * complex(imaginary = sign(x))
}

df <- fft_2d(nu1 = 1.5, nu2 = 1.5, a1 = 1, a2 = 1, Delta = Delta, Psi = Psi, grid_info = grid_info)

ggplot(data = df %>%
         filter(abs(Var1) < 5, abs(Var2) < 5), aes(x = Var2, y = Var1, fill = val)) +
  geom_raster() + 
  scale_fill_gradient2() +
  labs(x = 'Dimension 1', y = 'Dimension 2',
       fill = 'Cross-\ncovariance') +
  coord_equal() + 
  theme(legend.position = 'left',legend.key.height = unit(.8, "cm")) 

# load in data
library(R.matlab)
library(fields)
weather <- R.matlab::readMat('bolin_code/article_code/Application/TempPress/weather_data.mat')

weather <- as.data.frame(weather)
colnames(weather) <- c('lat', 'long', 'pres', 'temp')
n <- nrow(weather)
ggplot(data = weather, aes(x = long, y = lat, color = temp)) +
  geom_point()
ggplot(data = weather, aes(x = long, y = lat, color = pres)) +
  geom_point()

dist <- fields::rdist.earth(x1 = weather[, c('long', 'lat')],
                            x2 = weather[, c('long', 'lat')], miles = F)
diag(dist) <- 0
response <- unlist(weather[, c('pres', 'temp')])
response <- response -
  rep(colMeans(weather[, c('pres', 'temp')]), each = nrow(dist))

# create long/lat distance matrices
dist_one <- dist
dist_tens <- array(NA, c(dim(dist), 2))
reference_location <- colMeans(weather[, c('long', 'lat')])
for (i in 1:dim(dist_tens)[1]) {
  for (j in 1:dim(dist_tens)[1]) {
    dist_long_prelim <- fields::rdist.earth.vec(cbind(weather$long[i],
                                                      reference_location['lat']), 
                                                cbind(weather$long[j],
                                                      reference_location['lat']),
                                                miles = F)
    if (weather$long[i] < weather$long[j]) {
      dist_long_prelim <- -dist_long_prelim
    }
    
    dist_lat_prelim <- fields::rdist.earth.vec(cbind(reference_location['long'],
                                                     weather$lat[i]), 
                                               cbind(reference_location['long'],
                                                     weather$lat[j]),
                                               miles = F)
    if (weather$lat[i] < weather$lat[j]) {
      dist_lat_prelim <- -dist_lat_prelim
    }
    
    dist_tens[i,j,1] <- dist_long_prelim
    dist_tens[i,j,2] <- dist_lat_prelim
  }
}
dist_tens_mat <- cbind(as.vector(dist_tens[,,1]), as.vector(dist_tens[,,2]))



# given distances, compute fft on grid, then interpolate onto distances
construct_matrix <- function(nu1, nu2, a1, a2, 
                             grid_info,
                             Psi, Delta, dist_tens_mat) {
  C_val <- fft_2d(grid_info = grid_info, nu1 = nu1, nu2 = nu2, 
                  a1 = a1, a2 = a2,  Psi = Psi, Delta = Delta)
  
  tmat <- matrix(C_val[['val']], nrow = sqrt(nrow(C_val)))
  long1 <- grid_info[['x_vals']]
  lat1 <- grid_info[['x_vals']]
  
  interp_points <- fields::interp.surface(obj = list(x = long1, y = lat1, z = tmat),
                                          loc = dist_tens_mat)
  matrix(interp_points, nrow = sqrt(length(interp_points)))
}

# construct bivariate Matern covariance matrix
construct_entire_matrix <- function(nu1, nu2, a1, a2, 
                                    grid_info,
                                    Psi_list, Delta_list, dist_tens_mat, nugget1, nugget2) {
  C1 <- construct_matrix(nu1 = nu1, nu2 = nu1, a1 = a1, a2 = a1,
                         Psi = Psi_list[[1]], Delta = Delta_list[[1]], 
                         grid_info = grid_info, dist_tens_mat)
  diag(C1) <- diag(C1) + nugget1
  C12 <- construct_matrix(nu1 = nu1, nu2 = nu2, a1 = a1, a2 = a2,
                          Psi = Psi_list[[3]], Delta = Delta_list[[3]], 
                          grid_info = grid_info, dist_tens_mat)
  C2 <- construct_matrix(nu1 = nu2, nu2 = nu2, a1 = a2, a2 = a2,
                         Psi = Psi_list[[2]], Delta = Delta_list[[2]],
                         grid_info = grid_info, dist_tens_mat)
  diag(C2) <- diag(C2) + nugget2
  rbind(cbind(C1, C12), cbind(t(C12), C2))
}

grid_info <- create_grid_info(2^10, x_max = 2500)

# test the matrix creation
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
test_mat_real <- construct_entire_matrix(nu1 = nu1, nu2 = nu2, a1 = a1, a2 = a2, 
                                    Delta_list = list(function(x) Sigma11, function(x) Sigma22, 
                                                      function(x) {Sigma12re}),
                                    Psi_list = replicate(3, Psi_fun),
                                    grid_info = grid_info,
                                    dist_tens_mat = dist_tens_mat, 
                                    nugget1 = nugget1, nugget2 = nugget2)

response1 <- response[1:n]
response2 <- response[-c(1:n)]
sqrt(mean(response1^2))
sqrt(mean(response2^2))

do_cokriging <- function(cov_mat, folds, n, response1, response2) {
  rmse_matrix <- matrix(ncol = 7, nrow = max(folds))
  pred_values <- list()
  for (i in 1:nfolds) {
    r_have <- which(folds != i)
    r_ind <- dplyr::setdiff(1:n, r_have)
    
    response1_have <- response1[r_have]
    response1_out <- response1[r_ind]
    
    cov_mat_11 <- cov_mat[c(r_have, (n + 1):(2*n)), c(r_have, (n + 1):(2*n))]
    cov_mat_12 <- cov_mat[c(r_have, (n + 1):(2*n)), r_ind]
    cov_mat_22 <- cov_mat[r_ind, r_ind]
    
    cov_mat_11_single <- cov_mat[r_have, r_have]
    cov_mat_12_single <- cov_mat[r_have, r_ind]
    
    cov_mat_11_none <- cov_mat[(n + 1):(2*n), (n + 1):(2*n)]
    cov_mat_12_none <- cov_mat[(n + 1):(2*n), r_ind]
    
    pred <- t(cov_mat_12) %*% solve(cov_mat_11, c(response1_have, response2))
    var_est <- cov_mat_22 - t(cov_mat_12) %*% solve(cov_mat_11,cov_mat_12)
    
    pred_single <- t(cov_mat_12_single) %*% solve(cov_mat_11_single, response1_have)
    var_est_single <- cov_mat_22 - t(cov_mat_12_single) %*% solve(cov_mat_11_single,cov_mat_12_single)
    
    pred_none <- t(cov_mat_12_none) %*% solve(cov_mat_11_none, response2)
    var_est_none <- cov_mat_22 - t(cov_mat_12_none) %*% solve(cov_mat_11_none,cov_mat_12_none)
    
    rmse_matrix[i,] <- c((mean(response1_out^2)),
                         (mean((response1_out - pred)^2)),
                         (mean((response1_out - pred_single)^2)),
                         (mean((response1_out - pred_none)^2)),
                         (mean(diag(var_est))),
                         (mean(diag(var_est_single))), 
                         (mean(diag(var_est_none))))
    pred_values[[i]] <- list(pred, pred_single, pred_none, var_est, var_est_single, 
                             var_est_none, response1_out)
  }
  return(list(rmse_matrix, pred_values))
}

# test the matrix creation
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
  if (theta_star2 < 0) {
    ifelse(theta > theta_star2 & theta < theta_star2 + pi, 1, -1)
  } else {
    ifelse(theta > theta_star2 | theta < theta_star2 - pi, 1, -1)
  }
}
test_mat_im <- construct_entire_matrix(nu1 = nu1, nu2 = nu2, a1 = a1, a2 = a2, 
                                         Delta_list = list(function(x) Sigma11, function(x) Sigma22, 
                                                           function(x) {complex(real = Sigma12re, 
                                                                                imaginary = Sigma12im * Psi_fun2(x))}),
                                         Psi_list = replicate(3, Psi_fun),
                                         grid_info = grid_info,
                                         dist_tens_mat = dist_tens_mat, 
                                         nugget1 = nugget1, nugget2 = nugget2)


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
test_mat_fix <- construct_entire_matrix(nu1 = nu1, nu2 = nu2, a1 = a1, a2 = a2, 
                                       Delta_list = list(function(x) Sigma11, function(x) Sigma22, 
                                                         function(x) Sigma12re),
                                       Psi_list = replicate(3, function(x) sign(x)),
                                       grid_info = grid_info,
                                       dist_tens_mat = dist_tens_mat, 
                                       nugget1 = nugget1, nugget2 = nugget2)

par <- test_optim_single$par
nu1 <- exp(par[1]); a1 <- exp(par[2])
nugget1 <- exp(par[3]); nugget2 <- exp(par[4])
Sigma11 <- exp(par[5]); Sigma22 <- exp(par[6])
Sigma12 <- par[7]*sqrt(Sigma11)*sqrt(Sigma22)
test_mat_single <- construct_entire_matrix(nu1 = nu1, nu2 = nu1, a1 = a1, a2 = a1, 
                                   Delta_list = list(function(x) Sigma11, 
                                                     function(x) Sigma22, 
                                                     function(x) Sigma12),
                                   Psi_list = replicate(3, function(x) sign(x)),
                                   dist_tens_mat = dist_tens_mat, grid_info = grid_info,
                                   nugget1 = nugget1, nugget2 = nugget2)


construct_matrix_mm <- function(nu1, nu2, nu12, a1, a2, a12,
                                grid_info,
                                Psi_list, Delta_list, dist_tens_mat, nugget1, nugget2) {
  C1 <- construct_matrix(nu1 = nu1, nu2 = nu1, a1 = a1, a2 = a1,
                         Psi = Psi_list[[1]], Delta = Delta_list[[1]], 
                         grid_info = grid_info, dist_tens_mat)
  diag(C1) <- diag(C1) + nugget1
  C12 <- construct_matrix(nu1 = nu12, nu2 = nu12, a1 = a12, a2 = a12,
                          Psi = Psi_list[[3]], Delta = Delta_list[[3]], 
                          grid_info = grid_info, dist_tens_mat)
  C2 <- construct_matrix(nu1 = nu2, nu2 = nu2, a1 = a2, a2 = a2,
                         Psi = Psi_list[[2]], Delta = Delta_list[[2]],
                         grid_info = grid_info, dist_tens_mat)
  diag(C2) <- diag(C2) + nugget2
  rbind(cbind(C1, C12), cbind(t(C12), C2))
}
par <- test_optim_mm$par
nu1 <- exp(par[1]); nu2 <- exp(par[2]); nu12 <- exp(par[3])
a1 <- exp(par[4]); a2 <- exp(par[5]); a12 <- exp(par[6])
nugget1 <- exp(par[7]); nugget2 <- exp(par[8])
Sigma11 <- exp(par[9]); Sigma22 <- exp(par[10])
Sigma12 <- par[11]*sqrt(Sigma11)*sqrt(Sigma22)
test_mat_mm <- construct_matrix_mm(nu1 = nu1, nu2 = nu2, a1 = a1, a2 = a2, 
                               nu12 = nu12, a12 = a12,
                               Delta_list = list(function(x) Sigma11, function(x) Sigma22, 
                                                 function(x) Sigma12),
                               Psi_list = replicate(3, {function(x) {
                                 sign(x)
                               }}),
                               dist_tens_mat = dist_tens_mat, grid_info = grid_info,
                               nugget1 = nugget1, nugget2 = nugget2)


nfolds <- 5
set.seed(8)
folds <- c(sample(1:nfolds,nfolds, replace = F), sample(1:nfolds, n - nfolds, replace = T))
cokriging_real <- do_cokriging(test_mat_real, folds, n, response1, response2)
cokriging_im <- do_cokriging(test_mat_im, folds, n, response1, response2)
cokriging_fix <- do_cokriging(test_mat_fix, folds, n, response1, response2)
cokriging_single <- do_cokriging(test_mat_single, folds, n, response1, response2)
cokriging_mm <- do_cokriging(test_mat_mm, folds, n, response1, response2)
folds_n <- as.vector(table(folds))
apply(cokriging_single[[1]], 2, function(x) sqrt(sum(x * folds_n) / sum(folds_n)))
apply(cokriging_mm[[1]], 2, function(x) sqrt(sum(x * folds_n) / sum(folds_n)))
apply(cokriging_fix[[1]], 2, function(x) sqrt(sum(x * folds_n) / sum(folds_n)))
apply(cokriging_real[[1]], 2, function(x) sqrt(sum(x * folds_n) / sum(folds_n)))
apply(cokriging_im[[1]], 2, function(x) sqrt(sum(x * folds_n) / sum(folds_n)))


nfolds <- n
set.seed(10)
folds <- c(sample(1:nfolds, nfolds, replace = F), sample(1:nfolds, n - nfolds, replace = T))
cokriging_real <- do_cokriging(test_mat_real, folds, n, response1, response2)
cokriging_im <- do_cokriging(test_mat_im, folds, n, response1, response2)
cokriging_fix <- do_cokriging(test_mat_fix, folds, n, response1, response2)
cokriging_single <- do_cokriging(test_mat_single, folds, n, response1, response2)
cokriging_mm <- do_cokriging(test_mat_mm, folds, n, response1, response2)
folds_n <- rep(1, n)
apply(cokriging_single[[1]], 2, function(x) sqrt(sum(x * folds_n) / sum(folds_n)))
apply(cokriging_mm[[1]], 2, function(x) sqrt(sum(x * folds_n) / sum(folds_n)))
apply(cokriging_fix[[1]], 2, function(x) sqrt(sum(x * folds_n) / sum(folds_n)))
apply(cokriging_real[[1]], 2, function(x) sqrt(sum(x * folds_n) / sum(folds_n)))
apply(cokriging_im[[1]], 2, function(x) sqrt(sum(x * folds_n) / sum(folds_n)))

switch_variables <- function(test_mat, n) {
  test_mat2 = test_mat
  test_mat2[1:n, 1:n] = test_mat[(n + 1):(2*n), (n + 1):(2*n)]
  test_mat2[(n + 1):(2*n), (n + 1):(2*n)] = test_mat[1:n, 1:n]
  test_mat2[(n + 1):(2*n), 1:n] = t(test_mat2[(n + 1):(2*n), 1:n])
  test_mat2[1:n, (n + 1):(2*n)] = t(test_mat2[1:n, (n + 1):(2*n)])
  test_mat2
}

nfolds <- 157
set.seed(12)
folds <- c(sample(1:nfolds,nfolds, replace = F), sample(1:nfolds, n - nfolds, replace = T))
folds_n <- rep(1, n)
cokriging_real <- do_cokriging(switch_variables(test_mat_real, n), folds, n, response2, response1)
cokriging_im <- do_cokriging(switch_variables(test_mat_im, n), folds, n, response2, response1)
cokriging_fix <- do_cokriging(switch_variables(test_mat_fix, n), folds, n, response2, response1)
cokriging_single <- do_cokriging(switch_variables(test_mat_single, n), folds, n, response2, response1)
cokriging_mm <- do_cokriging(switch_variables(test_mat_mm, n), folds, n, response2, response1)
apply(cokriging_single[[1]], 2, function(x) sqrt(sum(x * folds_n) / sum(folds_n)))
apply(cokriging_mm[[1]], 2, function(x) sqrt(sum(x * folds_n) / sum(folds_n)))
apply(cokriging_fix[[1]], 2, function(x) sqrt(sum(x * folds_n) / sum(folds_n)))
apply(cokriging_real[[1]], 2, function(x) sqrt(sum(x * folds_n) / sum(folds_n)))
apply(cokriging_im[[1]], 2, function(x) sqrt(sum(x * folds_n) / sum(folds_n)))

nfolds <- 5
set.seed(12)
folds <- c(sample(1:nfolds,nfolds, replace = F), sample(1:nfolds, n - nfolds, replace = T))
folds_n <- as.vector(table(folds))
cokriging_real <- do_cokriging(switch_variables(test_mat_real, n), folds, n, response2, response1)
cokriging_im <- do_cokriging(switch_variables(test_mat_im, n), folds, n, response2, response1)
cokriging_fix <- do_cokriging(switch_variables(test_mat_fix, n), folds, n, response2, response1)
cokriging_single <- do_cokriging(switch_variables(test_mat_single, n), folds, n, response2, response1)
cokriging_mm <- do_cokriging(switch_variables(test_mat_mm, n), folds, n, response2, response1)
apply(cokriging_single[[1]], 2, function(x) sqrt(sum(x * folds_n) / sum(folds_n)))
apply(cokriging_mm[[1]], 2, function(x) sqrt(sum(x * folds_n) / sum(folds_n)))
apply(cokriging_fix[[1]], 2, function(x) sqrt(sum(x * folds_n) / sum(folds_n)))
apply(cokriging_real[[1]], 2, function(x) sqrt(sum(x * folds_n) / sum(folds_n)))
apply(cokriging_im[[1]], 2, function(x) sqrt(sum(x * folds_n) / sum(folds_n)))






