
source('code/multi_matern_source.R')
source('code/simulation/estimation_source_1d.R')
#grid_info <- create_grid_info_1d(2^12, 25)
grid_info <- create_grid_info_1d(2^14, 10)
library(fields)

test <- fft_1d(grid_info, re = 1, im = 0)
sim_grid <- as.data.frame(expand.grid(n = 300, gen = c('A', 'B', 'C'), 
                                      Sigma12_re = c(0, 0.4), Sigma12_im = c(0, .4),
                                      nu1 = .5, nu2 = .75, nu12 = (.5 + .75)/2, a1 = 8, a2 = 12, a12 = 8, Sigma11 = 1, Sigma22 = 1,
                                      nugget1 = .1, nugget2 = .1, sim = 1:200))
sim_grid <- sim_grid[!((sim_grid$Sigma12_im == 0) & (sim_grid$Sigma12_re == 0)),]

print(Sys.getenv('THISJOBVALUE'))
job = as.integer(Sys.getenv('THISJOBVALUE'))
params <- sim_grid[job,]
set.seed(params[['sim']])
n <- params[['n']]

pred_loc <- seq(0, 1, length.out = 100)
loc <- runif(params[['n']])
loc <- sort(loc)

dists <- create_distance_matrices(loc, pred_loc)
dist_all <- dists[[1]]
dist_pred <- dists[[2]]
dist_pred_data <- dists[[3]]
dist <- dists[[4]]

if (params[['gen']] == 'C') {
  test_mat <- construct_bivariate_matrix_mm(nu1 = params[['nu1']], nu2 = params[['nu2']], nu12 = params[['nu12']],
                                            a1 = params[['a1']], a2 = params[['a2']], a12 = params[['a12']],
                                            Sigma11 = params[['Sigma11']], Sigma22 = params[['Sigma22']], 
                                            Sigma12 =  complex(real = params[['Sigma12_re']],
                                                               imaginary = params[['Sigma12_im']]),
                                            dist_mat = dist_all, grid_info = grid_info,
                                            nugget1 = params[['nugget1']], nugget2 = params[['nugget2']])
  
} else if (params[['gen']] == 'B') {
  test_mat <- construct_bivariate_matrix_alt(nu1 = params[['nu1']], nu2 = params[['nu2']], 
                                             a1 = params[['a1']], a2 = params[['a2']], 
                                             Sigma11 = params[['Sigma11']], Sigma22 = params[['Sigma22']], 
                                             Sigma12 =  complex(real = params[['Sigma12_re']],
                                                                imaginary = params[['Sigma12_im']]),
                                             dist_mat = dist_all, grid_info = grid_info,
                                             nugget1 = params[['nugget1']], nugget2 = params[['nugget2']])
} else {
  test_mat <- construct_bivariate_matrix_orig(nu1 = params[['nu1']], nu2 = params[['nu2']], 
                                              a1 = params[['a1']], a2 = params[['a2']], 
                                              Sigma11 = params[['Sigma11']], Sigma22 = params[['Sigma22']], 
                                              Sigma12 =  complex(real = params[['Sigma12_re']],
                                                                 imaginary = params[['Sigma12_im']]),
                                              dist_mat = dist_all, grid_info = grid_info,
                                              nugget1 = params[['nugget1']], nugget2 = params[['nugget2']])
}
response_all <- as.vector(mvtnorm::rmvnorm(n = 1, sigma = test_mat))
indexes_data_dist <- 1:n
indexes_pred_dist <- dplyr::setdiff(1:(length(response_all)/2), indexes_data_dist)

indexes_data1 <- 1:n
indexes_data2 <-  (n + length(pred_loc) + 1):(2*n + length(pred_loc))
indexes_data <- c(indexes_data1, indexes_data2)
indexes_pred <- dplyr::setdiff(1:length(response_all), indexes_data)
response <- response_all[indexes_data]
response_pred <- response_all[indexes_pred]

dist <- dist_all[indexes_data_dist, indexes_data_dist]

test_optim <- optim(
  par = c(log(c(1.5, 1.5, 2, 2, .2, .2, var(response[1:n]), var(response[-c(1:n)]))), .01, .01),
  fn = ll_fun, dist_one = dist, response = response, 
  lower = c(log(c(.01, .01, 10^-3, 10^-3, .005, .005, .01, .01)), -1, -1),
  upper = c(log(c(5, 5, 40, 40, 1.5, 1.5, 20, 20)),1, 1),
  method = 'L-BFGS-B',
  grid_info = grid_info,
  control = list(parscale = c(rep(1, 8), .1, .1))
)
par <- test_optim$par
nu1 <- exp(par[1])
nu2 <- exp(par[2])
a1 <- exp(par[3])
a2 <- exp(par[4])
nugget1 <-  exp(par[5])
nugget2 <-  exp(par[6])
Sigma11 <- exp(par[7])
Sigma22 <- exp(par[8])
Sigma12re <- par[9]*sqrt(Sigma11)*sqrt(Sigma22)
Sigma12im <- par[10]*sqrt(Sigma11)*sqrt(Sigma22)
cov_mat <- construct_bivariate_matrix_orig(nu1 = nu1, nu2 = nu2, a1 = a1, a2 = a2, 
                                           Sigma11 = Sigma11, Sigma22 = Sigma22, 
                                           Sigma12 = complex(real = Sigma12re, imaginary = Sigma12im),
                                           dist_mat = dist_all, grid_info = grid_info,
                                           nugget1 = nugget1, nugget2 = nugget2)
pred <- cbind(response_pred,
              cov_mat[indexes_pred, indexes_data] %*% solve(cov_mat[indexes_data, indexes_data], response),
              cov_mat[indexes_pred, indexes_data1] %*% solve(cov_mat[indexes_data1, indexes_data1], response_all[indexes_data1]),
              cov_mat[indexes_pred, indexes_data2] %*% solve(cov_mat[indexes_data2, indexes_data2], response_all[indexes_data2]),
              diag(cov_mat[indexes_pred, indexes_pred] - 
                     cov_mat[indexes_pred, indexes_data] %*% solve(cov_mat[indexes_data, indexes_data]) %*% cov_mat[indexes_data, indexes_pred]),
              diag(cov_mat[indexes_pred, indexes_pred] - 
                     cov_mat[indexes_pred, indexes_data1] %*% solve(cov_mat[indexes_data1, indexes_data1]) %*% cov_mat[indexes_data1, indexes_pred]), 
              diag(cov_mat[indexes_pred, indexes_pred] - 
                     cov_mat[indexes_pred, indexes_data2] %*% solve(cov_mat[indexes_data2, indexes_data2]) %*% cov_mat[indexes_data2, indexes_pred]))

test_optim_real <- optim(
  par = c(log(c(1.5, 1.5, 2, 2, .2, .2, var(response[1:n]), var(response[-c(1:n)]))), .01),
  fn = ll_fun_re, dist_one = dist, response = response, 
  lower = c(log(c(.01, .01, 10^-3, 10^-3, .005, .005, .01, .01)), -1),
  upper = c(log(c(5, 5, 40, 40, 1.5, 1.5, 20, 20)),1),
  method = 'L-BFGS-B',
  grid_info = grid_info,
  control = list(parscale = c(rep(1, 8), .1))
)

par <- test_optim_real$par
nu1 <- exp(par[1])
nu2 <- exp(par[2])
a1 <- exp(par[3])
a2 <- exp(par[4])
nugget1 <-  exp(par[5])
nugget2 <-  exp(par[6])
Sigma11 <- exp(par[7])
Sigma22 <- exp(par[8])
Sigma12re <- par[9]*sqrt(Sigma11)*sqrt(Sigma22)
cov_mat <- construct_bivariate_matrix_orig(nu1 = nu1, nu2 = nu2, a1 = a1, a2 = a2, 
                                           Sigma11 = Sigma11, Sigma22 = Sigma22, 
                                           Sigma12 = complex(real = Sigma12re, imaginary = 0),
                                           dist_mat = dist_all, grid_info = grid_info,
                                           nugget1 = nugget1, nugget2 = nugget2)


pred_re <-  cbind(response_pred,
                  cov_mat[indexes_pred, indexes_data] %*% solve(cov_mat[indexes_data, indexes_data], response),
                  cov_mat[indexes_pred, indexes_data1] %*% solve(cov_mat[indexes_data1, indexes_data1], response_all[indexes_data1]),
                  cov_mat[indexes_pred, indexes_data2] %*% solve(cov_mat[indexes_data2, indexes_data2], response_all[indexes_data2]),
                  diag(cov_mat[indexes_pred, indexes_pred] - 
                         cov_mat[indexes_pred, indexes_data] %*% solve(cov_mat[indexes_data, indexes_data]) %*% cov_mat[indexes_data, indexes_pred]),
                  diag(cov_mat[indexes_pred, indexes_pred] - 
                         cov_mat[indexes_pred, indexes_data1] %*% solve(cov_mat[indexes_data1, indexes_data1]) %*% cov_mat[indexes_data1, indexes_pred]), 
                  diag(cov_mat[indexes_pred, indexes_pred] - 
                         cov_mat[indexes_pred, indexes_data2] %*% solve(cov_mat[indexes_data2, indexes_data2]) %*% cov_mat[indexes_data2, indexes_pred]))


test_optim_alt <- optim(
  par = c(log(c(1.5, 1.5, 2, 2, .2, .2, var(response[1:n]), var(response[-c(1:n)]))), .01, .01),
  fn = ll_fun_alt, dist_one = dist, response = response, 
  lower = c(log(c(.01, .01, 10^-3, 10^-3, .005, .005, .01, .01)), -1, -1),
  upper = c(log(c(5, 5, 40, 40, 1.5, 1.5, 20, 20)),1, 1),
  method = 'L-BFGS-B',
  grid_info = grid_info,
  control = list(parscale = c(rep(1, 8), .1, .1))
)

par <- test_optim_alt$par
nu1 <- exp(par[1])
nu2 <- exp(par[2])
a1 <- exp(par[3])
a2 <- exp(par[4])
nugget1 <-  exp(par[5])
nugget2 <-  exp(par[6])
Sigma11 <- exp(par[7])
Sigma22 <- exp(par[8])
Sigma12re <- par[9]*sqrt(Sigma11)*sqrt(Sigma22)
Sigma12im <- par[10]*sqrt(Sigma11)*sqrt(Sigma22)
cov_mat <- construct_bivariate_matrix_alt(nu1 = nu1, nu2 = nu2, a1 = a1, a2 = a2, 
                                          Sigma11 = Sigma11, Sigma22 = Sigma22, 
                                          Sigma12 = complex(real = Sigma12re, imaginary = Sigma12im),
                                          dist_mat = dist_all, grid_info = grid_info,
                                          nugget1 = nugget1, nugget2 = nugget2)
pred_alt <- cbind(response_pred,
                  cov_mat[indexes_pred, indexes_data] %*% solve(cov_mat[indexes_data, indexes_data], response),
                  cov_mat[indexes_pred, indexes_data1] %*% solve(cov_mat[indexes_data1, indexes_data1], response_all[indexes_data1]),
                  cov_mat[indexes_pred, indexes_data2] %*% solve(cov_mat[indexes_data2, indexes_data2], response_all[indexes_data2]),
                  diag(cov_mat[indexes_pred, indexes_pred] - 
                         cov_mat[indexes_pred, indexes_data] %*% solve(cov_mat[indexes_data, indexes_data]) %*% cov_mat[indexes_data, indexes_pred]),
                  diag(cov_mat[indexes_pred, indexes_pred] - 
                         cov_mat[indexes_pred, indexes_data1] %*% solve(cov_mat[indexes_data1, indexes_data1]) %*% cov_mat[indexes_data1, indexes_pred]), 
                  diag(cov_mat[indexes_pred, indexes_pred] - 
                         cov_mat[indexes_pred, indexes_data2] %*% solve(cov_mat[indexes_data2, indexes_data2]) %*% cov_mat[indexes_data2, indexes_pred]))


test_optim_real_alt <- optim(
  par = c(log(c(1.5, 1.5, 2, 2, .2, .2, var(response[1:n]), var(response[-c(1:n)]))), .01),
  fn = ll_fun_alt_re, dist_one = dist, response = response, 
  lower = c(log(c(.01, .01, 10^-3, 10^-3, .005, .005, .01, .01)), -1),
  upper = c(log(c(5, 5, 40, 40, 1.5, 1.5, 20, 20)),1),
  method = 'L-BFGS-B',
  grid_info = grid_info,
  control = list(parscale = c(rep(1, 8), .1))
)

par <- test_optim_real_alt$par
nu1 <- exp(par[1])
nu2 <- exp(par[2])
a1 <- exp(par[3])
a2 <- exp(par[4])
nugget1 <-  exp(par[5])
nugget2 <-  exp(par[6])
Sigma11 <- exp(par[7])
Sigma22 <- exp(par[8])
Sigma12re <- par[9]*sqrt(Sigma11)*sqrt(Sigma22)
cov_mat <- construct_bivariate_matrix_alt(nu1 = nu1, nu2 = nu2, a1 = a1, a2 = a2, 
                                          Sigma11 = Sigma11, Sigma22 = Sigma22, 
                                          Sigma12 = complex(real = Sigma12re, imaginary = 0),
                                          dist_mat = dist_all, grid_info = grid_info,
                                          nugget1 = nugget1, nugget2 = nugget2)


pred_re_alt <-  cbind(response_pred,
                      cov_mat[indexes_pred, indexes_data] %*% solve(cov_mat[indexes_data, indexes_data], response),
                      cov_mat[indexes_pred, indexes_data1] %*% solve(cov_mat[indexes_data1, indexes_data1], response_all[indexes_data1]),
                      cov_mat[indexes_pred, indexes_data2] %*% solve(cov_mat[indexes_data2, indexes_data2], response_all[indexes_data2]),
                      diag(cov_mat[indexes_pred, indexes_pred] - 
                             cov_mat[indexes_pred, indexes_data] %*% solve(cov_mat[indexes_data, indexes_data]) %*% cov_mat[indexes_data, indexes_pred]),
                      diag(cov_mat[indexes_pred, indexes_pred] - 
                             cov_mat[indexes_pred, indexes_data1] %*% solve(cov_mat[indexes_data1, indexes_data1]) %*% cov_mat[indexes_data1, indexes_pred]), 
                      diag(cov_mat[indexes_pred, indexes_pred] - 
                             cov_mat[indexes_pred, indexes_data2] %*% solve(cov_mat[indexes_data2, indexes_data2]) %*% cov_mat[indexes_data2, indexes_pred]))

test_optim_mm <- optim(
  par = c(log(c(.25, .25, .25, 10, 10, 10, .2, .2, var(response[1:n]), var(response[-c(1:n)]))), .01, .01),
  fn = ll_fun_mm, dist_one = dist, response = response, 
  lower = c(log(c(.2, .2, .2, 5, 5, 5, .005, .005, .01, .01)), -1,-1),
  upper = c(log(c(3, 3, 3, 15, 15, 15, 1.5, 1.5, 20, 20)),1, 1),
  # lower = c(log(c(.01, .01, .01, 2, 2,2, .005, .005, .01, .01)), -1, -1),
  # upper = c(log(c(1.5, 1.5,1.5,  15, 15, 15, 1.5, 1.5, 20, 20)),1, 1),
  method = 'L-BFGS-B',
  grid_info = grid_info,
  control = list(parscale = c(rep(1, 10), .1, .1))
)

par <- test_optim_mm$par
nu1 <- exp(par[1])
nu2 <- exp(par[2])
nu12 <- exp(par[3])
a1 <- exp(par[4])
a2 <- exp(par[5])
a12 <- exp(par[6])
nugget1 <- exp(par[7])
nugget2 <- exp(par[8])
Sigma11 <- exp(par[9])
Sigma22 <- exp(par[10])
Sigma12re <- par[11]*sqrt(Sigma11)*sqrt(Sigma22)
Sigma12im <- par[12]*sqrt(Sigma11)*sqrt(Sigma22)
cov_mat <- construct_bivariate_matrix_mm(nu1 = nu1, nu2 = nu2, nu12 = nu12, a1 = a1, a2 = a2, 
                                         a12 = a12,
                                         Sigma11 = Sigma11, Sigma22 = Sigma22, 
                                         Sigma12 = complex(real = Sigma12re, imaginary = Sigma12im),
                                         dist_mat = dist_all, grid_info = grid_info,
                                         nugget1 = nugget1, nugget2 = nugget2)

pred_mm <-  cbind(response_pred,
                  cov_mat[indexes_pred, indexes_data] %*% solve(cov_mat[indexes_data, indexes_data], response),
                  cov_mat[indexes_pred, indexes_data1] %*% solve(cov_mat[indexes_data1, indexes_data1], response_all[indexes_data1]),
                  cov_mat[indexes_pred, indexes_data2] %*% solve(cov_mat[indexes_data2, indexes_data2], response_all[indexes_data2]),
                  diag(cov_mat[indexes_pred, indexes_pred] - 
                         cov_mat[indexes_pred, indexes_data] %*% solve(cov_mat[indexes_data, indexes_data]) %*% cov_mat[indexes_data, indexes_pred]),
                  diag(cov_mat[indexes_pred, indexes_pred] - 
                         cov_mat[indexes_pred, indexes_data1] %*% solve(cov_mat[indexes_data1, indexes_data1]) %*% cov_mat[indexes_data1, indexes_pred]), 
                  diag(cov_mat[indexes_pred, indexes_pred] - 
                         cov_mat[indexes_pred, indexes_data2] %*% solve(cov_mat[indexes_data2, indexes_data2]) %*% cov_mat[indexes_data2, indexes_pred]))



test_optim_mm_re <- optim(
  par = c(log(c(1.5, 1.5, 1.5, 10, 10, 10, .2, .2, var(response[1:n]), var(response[-c(1:n)]))), .01),
  fn = ll_fun_mm_re, dist_one = dist, response = response, 
  # lower = c(log(c(.01, .01, .01, 10^-3, 10^-3,10^-3, .005, .005, .01, .01)), -1),
  # upper = c(log(c(1.75, 1.75,1.75,  20, 20, 20, 1.5, 1.5, 20, 20)),1),
  lower = c(log(c(.01, .01, .01, 5, 5,5, .005, .005, .01, .01)), -1),
  upper = c(log(c(5, 5, 5, 15, 15, 15, 1.5, 1.5, 20, 20)),1),
  method = 'L-BFGS-B',
  grid_info = grid_info,
  control = list(parscale = c(rep(1, 10), .1))
)

par <- test_optim_mm_re$par
nu1 <- exp(par[1])
nu2 <- exp(par[2])
nu12 <- exp(par[3])
a1 <- exp(par[4])
a2 <- exp(par[5])
a12 <- exp(par[6])
nugget1 <- exp(par[7])
nugget2 <- exp(par[8])
Sigma11 <- exp(par[9])
Sigma22 <- exp(par[10])
Sigma12re <- par[11]*sqrt(Sigma11)*sqrt(Sigma22)
cov_mat <- construct_bivariate_matrix_mm(nu1 = nu1, nu2 = nu2, nu12 = nu12, a1 = a1, a2 = a2, 
                                         a12 = a12,
                                         Sigma11 = Sigma11, Sigma22 = Sigma22, 
                                         Sigma12 = complex(real = Sigma12re, imaginary = 0),
                                         dist_mat = dist_all, grid_info = grid_info,
                                         nugget1 = nugget1, nugget2 = nugget2)

pred_mm_re <-  cbind(response_pred,
                     cov_mat[indexes_pred, indexes_data] %*% solve(cov_mat[indexes_data, indexes_data], response),
                     cov_mat[indexes_pred, indexes_data1] %*% solve(cov_mat[indexes_data1, indexes_data1], response_all[indexes_data1]),
                     cov_mat[indexes_pred, indexes_data2] %*% solve(cov_mat[indexes_data2, indexes_data2], response_all[indexes_data2]),
                     diag(cov_mat[indexes_pred, indexes_pred] - 
                            cov_mat[indexes_pred, indexes_data] %*% solve(cov_mat[indexes_data, indexes_data]) %*% cov_mat[indexes_data, indexes_pred]),
                     diag(cov_mat[indexes_pred, indexes_pred] - 
                            cov_mat[indexes_pred, indexes_data1] %*% solve(cov_mat[indexes_data1, indexes_data1]) %*% cov_mat[indexes_data1, indexes_pred]), 
                     diag(cov_mat[indexes_pred, indexes_pred] - 
                            cov_mat[indexes_pred, indexes_data2] %*% solve(cov_mat[indexes_data2, indexes_data2]) %*% cov_mat[indexes_data2, indexes_pred]))


test_optim_mm_pars <- optim(
  par = c(log(c(.25, .25, 2, .2, .2, var(response[1:n]), var(response[-c(1:n)]))), .01, .01),
  fn = ll_fun_mm_pars, dist_one = dist, response = response, 
  lower = c(log(c(.01, .01, 10^-3, .005, .005, .01, .01)), -1,-1),
  upper = c(log(c(5, 5, 40, 1.5, 1.5, 20, 20)),1, 1),
  # lower = c(log(c(.01, .01, .01, 2, 2,2, .005, .005, .01, .01)), -1, -1),
  # upper = c(log(c(1.5, 1.5,1.5,  15, 15, 15, 1.5, 1.5, 20, 20)),1, 1),
  method = 'L-BFGS-B',
  grid_info = grid_info,
  control = list(parscale = c(rep(1, 7), .1, .1))
)

par <- test_optim_mm_pars$par
nu1 <- exp(par[1])
nu2 <- exp(par[2])
nu12 <- (nu1 + nu2)/2
a1 <- a2 <- a12 <- exp(par[3])
nugget1 <- exp(par[4])
nugget2 <- exp(par[5])
Sigma11 <- exp(par[6])
Sigma22 <- exp(par[7])
Sigma12re <- par[8]*sqrt(Sigma11)*sqrt(Sigma22)
Sigma12im <- par[9]*sqrt(Sigma11)*sqrt(Sigma22)
cov_mat <- construct_bivariate_matrix_mm(nu1 = nu1, nu2 = nu2, nu12 = nu12, a1 = a1, a2 = a2, 
                                         a12 = a12,
                                         Sigma11 = Sigma11, Sigma22 = Sigma22, 
                                         Sigma12 = complex(real = Sigma12re, imaginary = Sigma12im),
                                         dist_mat = dist_all, grid_info = grid_info,
                                         nugget1 = nugget1, nugget2 = nugget2)

pred_mm_pars <-  cbind(response_pred,
                       cov_mat[indexes_pred, indexes_data] %*% solve(cov_mat[indexes_data, indexes_data], response),
                       cov_mat[indexes_pred, indexes_data1] %*% solve(cov_mat[indexes_data1, indexes_data1], response_all[indexes_data1]),
                       cov_mat[indexes_pred, indexes_data2] %*% solve(cov_mat[indexes_data2, indexes_data2], response_all[indexes_data2]),
                       diag(cov_mat[indexes_pred, indexes_pred] - 
                              cov_mat[indexes_pred, indexes_data] %*% solve(cov_mat[indexes_data, indexes_data]) %*% cov_mat[indexes_data, indexes_pred]),
                       diag(cov_mat[indexes_pred, indexes_pred] - 
                              cov_mat[indexes_pred, indexes_data1] %*% solve(cov_mat[indexes_data1, indexes_data1]) %*% cov_mat[indexes_data1, indexes_pred]), 
                       diag(cov_mat[indexes_pred, indexes_pred] - 
                              cov_mat[indexes_pred, indexes_data2] %*% solve(cov_mat[indexes_data2, indexes_data2]) %*% cov_mat[indexes_data2, indexes_pred]))



test_optim_mm_re_pars <- optim(
  par = c(log(c(1.5, 1.5, 2, .2, .2, var(response[1:n]), var(response[-c(1:n)]))), .01),
  fn = ll_fun_mm_re_pars, dist_one = dist, response = response, 
  lower = c(log(c(.01, .01, 10^-3, .005, .005, .01, .01)), -1),
  upper = c(log(c(5, 5, 40, 1.5, 1.5, 20, 20)),1),
  method = 'L-BFGS-B',
  grid_info = grid_info,
  control = list(parscale = c(rep(1, 7), .1))
)

par <- test_optim_mm_re_pars$par
nu1 <- exp(par[1])
nu2 <- exp(par[2])
nu12 <- (nu1 + nu2)/2
a1 <- a2 <- a12 <- exp(par[3])
nugget1 <- exp(par[4])
nugget2 <- exp(par[5])
Sigma11 <- exp(par[6])
Sigma22 <- exp(par[7])
Sigma12re <- par[8]*sqrt(Sigma11)*sqrt(Sigma22)
cov_mat <- construct_bivariate_matrix_mm(nu1 = nu1, nu2 = nu2, nu12 = nu12, a1 = a1, a2 = a2, 
                                         a12 = a12,
                                         Sigma11 = Sigma11, Sigma22 = Sigma22, 
                                         Sigma12 = complex(real = Sigma12re, imaginary = 0),
                                         dist_mat = dist_all, grid_info = grid_info,
                                         nugget1 = nugget1, nugget2 = nugget2)

pred_mm_re_pars <-  cbind(response_pred,
                          cov_mat[indexes_pred, indexes_data] %*% solve(cov_mat[indexes_data, indexes_data], response),
                          cov_mat[indexes_pred, indexes_data1] %*% solve(cov_mat[indexes_data1, indexes_data1], response_all[indexes_data1]),
                          cov_mat[indexes_pred, indexes_data2] %*% solve(cov_mat[indexes_data2, indexes_data2], response_all[indexes_data2]),
                          diag(cov_mat[indexes_pred, indexes_pred] - 
                                 cov_mat[indexes_pred, indexes_data] %*% solve(cov_mat[indexes_data, indexes_data]) %*% cov_mat[indexes_data, indexes_pred]),
                          diag(cov_mat[indexes_pred, indexes_pred] - 
                                 cov_mat[indexes_pred, indexes_data1] %*% solve(cov_mat[indexes_data1, indexes_data1]) %*% cov_mat[indexes_data1, indexes_pred]), 
                          diag(cov_mat[indexes_pred, indexes_pred] - 
                                 cov_mat[indexes_pred, indexes_data2] %*% solve(cov_mat[indexes_data2, indexes_data2]) %*% cov_mat[indexes_data2, indexes_pred]))

save(test_optim, test_optim_real, test_optim_alt, test_optim_real_alt, test_optim_mm, test_optim_mm_re, 
     test_optim_mm_pars, test_optim_mm_re_pars,
     pred, pred_re, pred_alt, pred_re_alt, pred_mm, pred_mm_re, pred_mm_pars, pred_mm_re_pars, 
     file = paste0('results/sim_results_1d/', paste(unlist(params), collapse = '_'), '.RData'))

q(save = 'no')
