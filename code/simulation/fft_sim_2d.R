
source('code/multi_matern_source.R')
source('code/simulation/estimation_source_2d.R')
grid_info <- create_grid_info_2d(2^8, 4)
library(fields)

Mu <- function(theta_x, theta_y, entry_1, entry_2) {
  complex(imaginary = sign(theta_x))
}
Psi <- function(theta_x, theta_y) {
  sign(theta_x)
}
test <- fft_2d(grid_info, d = 2, nu1 = .5, nu2 = .5, Mu = Mu, Psi = Psi, a1 = 1, a2 = 1)
library(dplyr)
library(ggplot2)
ggplot(data = test %>% filter(abs(Var1) < 1, abs(Var2) < 1), aes(x = Var1, y = Var2, fill = val)) +
  geom_raster() + 
  scale_fill_gradient2()

sim_grid <- as.data.frame(expand.grid(n = 300, gen = c('A'), 
                                      Sigma12_re = c(0, 0.4), Sigma12_im = c(0, .4),
                                      angle = c(0, 1),
                                      nu1 = .5, nu2 = .75, nu12 = (.5 + .75)/2, a1 = 8, a2 = 12, a12 = 8, Sigma11 = 1, Sigma22 = 1,
                                      nugget1 = .01, nugget2 = .01, sim = 1:200))
sim_grid <- sim_grid[!((sim_grid$Sigma12_im == 0) & (sim_grid$Sigma12_re == 0)),]

print(Sys.getenv('THISJOBVALUE'))
job = as.integer(Sys.getenv('THISJOBVALUE'))
params <- sim_grid[job,]
set.seed(params[['sim']])
n <- params[['n']]

pred_loc <- expand.grid(seq(0, 1, length.out = 100),
                        seq(0, 1, length.out = 100))
loc <- cbind(runif(params[['n']]), runif(params[['n']]))
loc <- loc[order(loc[,1]),]
dist_mats <- create_distance_matrices(loc, pred_loc)
dist1 <- dist_mats[,,1]
dist2 <- dist_mats[,,2]


dist_vals <- cbind(as.vector(dist_mats[,,1]), as.vector(dist_mats[,,2]))

if (params[['angle']] == 0) {
  Mu <- function(theta_x, theta_y, entry_1, entry_2) {
    complex(real = params[['Sigma12_re']], imaginary = sign(theta_x) * params[['Sigma12_im']])
  }
  Psi <- function(theta_x, theta_y) {
    sign(theta_x)
  }
} else {
  Mu <- function(theta_x, theta_y, entry_1, entry_2) {
    complex(real = params[['Sigma12_re']], imaginary = sign(cos(theta_x + pi/4)) * params[['Sigma12_im']])
  }
  Psi <- function(theta_x, theta_y) {
    sign(cos(theta_x + pi/4))
  }
}

test_mat <- construct_bivariate_matrix_orig(nu1 = params[['nu1']], nu2 = params[['nu2']], 
                                            a1 = params[['a1']], a2 = params[['a2']], 
                                            Sigma11 = params[['Sigma11']], Sigma22 = params[['Sigma22']], 
                                            Psi = Psi, Mu = Mu, dist_vals = dist_vals,
                                            grid_info = grid_info,
                                            nugget1 = params[['nugget1']], nugget2 = params[['nugget2']])
response_all <- as.vector(mvtnorm::rmvnorm(n = 1, sigma = test_mat))

response <- response_all

nu1_init <- nu2_init <- .5
nu_lower <- 0.001; nu_higher = 7
a1_init <- a2_init <- 10
a_lower <- 10^-6; a_upper = 20
Sigma11_init <- var(response[1:n])
Sigma22_init <- var(response[-c(1:n)])
cor_re_init <- 0
cor_im_init <- 0
nugget1_init <- var(response[1:n])/25
nugget2_init <- var(response[-c(1:n)])/25
var_lower <- 10^-3



#### with psi fixed at sign(theta)
test_optim_im_psi_fixed <- optim(
  par = c(log(c(nu1_init, nu2_init, a1_init, a2_init, nugget1_init, nugget2_init, 
                Sigma11_init, Sigma22_init)), cor_re_init, cor_im_init),
  fn = ll_fun, dist_vals = dist_vals, response = response,# n1 = n1, n2 = n2,
  lower = c(log(c(nu_lower, nu_lower, a_lower, a_lower, var_lower, var_lower, var_lower, var_lower)), -.999, -.999),
  upper = c(log(c(nu_higher, nu_higher, a_upper, a_upper, NA, NA, NA, NA)), .999,.999),
  method = 'L-BFGS-B',
  grid_info = grid_info,
  control = list(parscale = c(rep(1, 8), .1, .1))
)

test_optim_re_psi_fixed <- optim(
  par = c(log(c(nu1_init, nu2_init, a1_init, a2_init, nugget1_init, nugget2_init, 
                Sigma11_init, Sigma22_init)), cor_re_init),
  fn = ll_fun_re, dist_vals = dist_vals, response = response,# n1 = n1, n2 = n2,
  lower = c(log(c(nu_lower, nu_lower, a_lower, a_lower, var_lower, var_lower, var_lower, var_lower)), -.999),
  upper = c(log(c(nu_higher, nu_higher, a_upper, a_upper, NA, NA, NA, NA)), .999),
  method = 'L-BFGS-B',
  grid_info = grid_info,
  control = list(parscale = c(rep(1, 8), .1))
)

test_optim_re <- optim(
  par = c(c(log(c(nu1_init, nu2_init, a1_init, a2_init, nugget1_init, nugget2_init, 
                  Sigma11_init, Sigma22_init)), cor_re_init), 0.1, 0.1),
  fn = ll_fun_psi_real, dist_vals = dist_vals, response = response, 
  lower = c(log(c(nu_lower, nu_lower, a_lower, a_lower, var_lower, var_lower, var_lower, var_lower)), -.999, NA, NA),
  upper = c(log(c(nu_higher, nu_higher, a_upper, a_upper, NA, NA, NA, NA)), .999, NA, NA),
  method = 'L-BFGS-B',
  grid_info = grid_info,
  control = list(parscale = c(rep(1, 8), .1, .1, .1))
)

test_optim_im <- optim(
  par = c(log(c(nu1_init, nu2_init, a1_init, a2_init, nugget1_init, nugget2_init, 
                Sigma11_init, Sigma22_init)), cor_re_init, cor_im_init, .1, .1, .1, .1),
  fn = ll_fun_psi_im,  dist_vals = dist_vals, response = response, 
  lower = c(log(c(nu_lower, nu_lower, a_lower, a_lower, var_lower, var_lower, var_lower, var_lower)), -.999, -.999, NA, NA, NA, NA),
  upper = c(log(c(nu_higher, nu_higher, a_upper, a_upper, NA, NA, NA, NA)), .999, .999, NA, NA, NA, NA),
  method = 'L-BFGS-B',
  grid_info = grid_info,
  control = list(parscale = c(rep(1, 8), .1, .1, .1, .1, .1, .1))
)

save(test_optim_im, test_optim_re,test_optim_re_psi_fixed,test_optim_im_psi_fixed,
     file = paste0('results/sim_results_2d/', paste(unlist(params), collapse = '_'), '.RData'))

q(save = 'no')
