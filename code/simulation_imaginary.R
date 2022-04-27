
library(parallel)
source('code/multi_matern_source.R')
n_vals <- c(100, 150, 200, 250, 300, 350, 400)
n_simu <- 200

# try with imaginary part
# implement marginal first, then joint
set.seed(60)
seeds <- sample(1:(10^6), (length(n_vals) * n_simu))
simu_info <- data.frame(expand.grid('simu' = 1:n_simu,
                                    'n' = n_vals),
                        seeds)

likelihood <- function(theta, response, locs) {
  if (abs(theta[3]) > sqrt(exp(theta[2])) * sqrt(exp(theta[1]))) {
    return(1000)
  }
  cov_mat <- imaginary_covariance_matrix_lags(locs, 
                                              nu1 = exp(theta[4]), 
                                              nu2 = exp(theta[4]),
                                              c11 = exp(theta[1]), 
                                              c22 = exp(theta[2]), 
                                              c12 = complex(imaginary = theta[3]), a1 = 1,
                                              a2 = 1)
  c_chol <- base::chol(cov_mat)
  v2 <- .Internal(backsolve(r = c_chol, x = response, 
                            k = ncol(c_chol), upper.tri = T, transpose = T))
  quad_form <- sum(v2^2)
  det_val <-  2* sum(log(diag(c_chol)))
  quad_form+ det_val
}

likelihood_single_var <- function(theta, response, locs, type) {
  cov_mat <- whittaker_covariance_matrix_lags(locs, 
                                              nu1 = exp(theta[2]), 
                                              nu2 = .949,
                                              c11 = exp(theta[1]), 
                                              c2 = 1.22, 
                                              c12 = .01, a1 = 1,
                                              a2 = 1)
  
  n_obs <- length(locs)
  if (type ==1) {
    c_chol <- base::chol(cov_mat[1:n_obs, 1:n_obs])
    v2 <- .Internal(backsolve(r = c_chol, x = response[1:n_obs], 
                              k = n_obs, upper.tri = T, transpose = T))
  } else if (type == 2) {
    c_chol <- base::chol(cov_mat[1:n_obs, 1:n_obs])
    v2 <- .Internal(backsolve(r = c_chol, x = response[-c(1:n_obs)], 
                              k = n_obs, upper.tri = T, transpose = T))
  }
  quad_form <- sum(v2^2)
  det_val <-  2* sum(log(diag(c_chol)))
  quad_form+ det_val
}

likelihood_joint_only <- function(theta, response, locs, estimated_params) {
  if (abs(theta) > sqrt(exp(estimated_params[2])) * sqrt(exp(estimated_params[1]))) {
    return(1000)
  }
  cov_mat <- imaginary_covariance_matrix_lags(locs, 
                                              nu1 = exp(estimated_params[3]), 
                                              nu2 = exp(estimated_params[3]),
                                              c11 = exp(estimated_params[1]), 
                                              c2 = exp(estimated_params[2]), 
                                              c12 = complex(imaginary = theta), a1 = 1,
                                              a2 = 1)
  n_obs <- length(locs)
  c_chol <- base::chol(cov_mat)
  v2 <- .Internal(backsolve(r = c_chol, x = response, 
                            k = n_obs, upper.tri = T, transpose = T))
  quad_form <- sum(v2^2)
  det_val <-  2* sum(log(diag(c_chol)))
  quad_form+ det_val
}


likelihood_nu_only <- function(theta, response, locs, estimated_params) {
  cov_mat <- imaginary_covariance_matrix_lags(locs, 
                                              nu1 = exp(theta[1]), 
                                              nu2 = exp(theta[1]),
                                              c11 = estimated_params[1], 
                                              c2 = estimated_params[2], 
                                              c12 = estimated_params[3], a1 = 1,
                                              a2 = 1)
  n_obs <- length(locs)
  c_chol <- base::chol(cov_mat)
  v2 <- .Internal(backsolve(r = c_chol, x = response, 
                            k = n_obs, upper.tri = T, transpose = T))
  quad_form <- sum(v2^2)
  det_val <-  2* sum(log(diag(c_chol)))
  quad_form+ det_val
}


estimate_from_simulation <- function(array_id, simulations, locs, simu_info) {
  print(array_id)
  n <- simu_info[['n']][array_id]
  seed <- simu_info[['seeds']][array_id]
  set.seed(seed)
  locs <- (1:n)/n
  
  cov_mat <- imaginary_covariance_matrix_lags(locs, nu1 = .8, nu2 = .8, c11 = 1, 
                                              c22 = 1, c12 = complex(imaginary = .8), 
                                              a1 = 1, a2 = 1)
  cov_mat_chol <- base::chol(cov_mat)
  simulations <- simulate_manually(cholesky = cov_mat_chol, n_simu = 1, locs)
  
  response <- c(simulations[['var1']], simulations[['var2']])
  
  joint_optim <- optim(par = c(log(1.3),log(.7), .2, log(.51)), fn = likelihood,
                       response = response, method = 'L-BFGS-B', 
                       hessian = T, lower = c(-1, -1, -2,-2),
                       upper = c(1, 1, 2,.25), locs = locs)
  single_optim1 <- optim(par = c(log(1.3),log(.51)), fn = likelihood_single_var, 
                         response = response, type = 1,
                         locs = locs, method = 'L-BFGS-B', hessian = T, 
                         lower = c(-1, -2), upper = c(1,.25))
  single_optim2 <- optim(par = c(log(1.3),log(.51)), fn = likelihood_single_var, 
                         response = response, type = 2,
                         locs = locs, method = 'L-BFGS-B', hessian = T, 
                         lower = c(-1, -2), upper = c(1,.25))
  
  joint_optim_init <- optim(par = c(single_optim1$par[1],single_optim2$par[1], .2,
                                    mean(c(single_optim1$par[2],single_optim2$par[2]))), fn = likelihood,
                            response = response, method = 'L-BFGS-B', 
                            hessian = T, lower = c(-1, -1, -2,-2, -2),
                            upper = c(1, 1, 2,.25, .25), locs = locs)
  
  joint_optim_seq <- optim(par = .2, fn = likelihood_joint_only, estimated_params = 
                             c(single_optim1$par[1],single_optim2$par[1],
                               mean(c(single_optim1$par[2],single_optim2$par[2]))),
                           response = response, method = 'L-BFGS-B', 
                           hessian = T, lower = -2,
                           upper =2, locs = locs)
  
  AA_star_mom <- crossprod(cbind(simulations[['var1']], simulations[['var2']]))/n
  
  joint_optim_mom <- optim(par = c(log(.51)), fn = likelihood_nu_only, estimated_params = 
                             c(AA_star_mom[1,1], AA_star_mom[2,2], AA_star_mom[1,2]),
                           response = response, method = 'L-BFGS-B', 
                           hessian = T, lower = c(-2),
                           upper = c(.25), locs = locs)
  
  list(joint_optim, single_optim1, single_optim2, joint_optim_init, joint_optim_seq, 
       AA_star_mom, joint_optim_mom)
}
simu_results <- mclapply(1:nrow(simu_info), estimate_from_simulation, 
                         simulations = simulations, locs = locs,
                         simu_info = simu_info,
                         mc.cores = 12, mc.preschedule = F)

save(simu_results, simu_info, file = 'results/simu_results_imaginary.RData')
