load('results/pres_temp_actual_spectral_results_10_final.RData')
library(tidyverse)
library(fftw)
library(fftwtools)
library(R.matlab)
library(fields)
source('code/multi_matern_source.R')

# load in data
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
response1 <- response[1:n]
response2 <- response[-c(1:n)]
sqrt(mean(response1^2))
sqrt(mean(response2^2))

# create joint matrices for each model

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
construct_bivariate_matrix <- function(nu1, nu2, a1, a2, 
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
grid_info <- create_grid_info_2d(2^10, x_max = 2500)



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
test_mat_real <- construct_bivariate_matrix(nu1 = nu1, nu2 = nu2, a1 = a1, a2 = a2, 
                                    Delta_list = list(function(x) Sigma11, function(x) Sigma22, 
                                                      function(x) {Sigma12re}),
                                    Psi_list = replicate(3, Psi_fun),
                                    grid_info = grid_info,
                                    dist_tens_mat = dist_tens_mat, 
                                    nugget1 = nugget1, nugget2 = nugget2)


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
test_mat_im <- construct_bivariate_matrix(nu1 = nu1, nu2 = nu2, a1 = a1, a2 = a2, 
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
test_mat_fix <- construct_bivariate_matrix(nu1 = nu1, nu2 = nu2, a1 = a1, a2 = a2, 
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
test_mat_single <- construct_bivariate_matrix(nu1 = nu1, nu2 = nu1, a1 = a1, a2 = a1, 
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






