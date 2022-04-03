
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
  if ((abs(theta[3]) + abs(theta[4])) > sqrt(exp(theta[2])) * sqrt(exp(theta[1]))) {
    return(10^6)
  }
  cov_mat <- imaginary_covariance_matrix_lags(locs, 
                                              nu1 = exp(theta[5]), 
                                              nu2 = exp(theta[5]),
                                              c11 = exp(theta[1]), 
                                              c22 = exp(theta[2]), 
                                              c12 = complex(real = theta[3], 
                                                            imaginary = theta[4]), a1 = 1,
                                              a2 = 1)
  c_chol <- base::chol(cov_mat)
  v2 <- .Internal(backsolve(r = c_chol, x = response, 
                            k = ncol(c_chol), upper.tri = T, transpose = T))
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
                                              c22 = 1, c12 = complex(real = .3, imaginary = .3), 
                                              a1 = 1, a2 = 1)
  cov_mat_chol <- base::chol(cov_mat)
  simulations <- simulate_manually(cholesky = cov_mat_chol, n_simu = 1, locs)
  
  response <- c(simulations[['var1']], simulations[['var2']])
  
  joint_optim <- optim(par = c(log(1.3),log(.7), 0, 0, log(.51)), fn = likelihood,
                       response = response, method = 'L-BFGS-B', 
                       hessian = T, lower = c(-1, -1, -2, -2,-2),
                       upper = c(1, 1, 2, 2, .25), locs = locs, 
                       control = list(parscale = c(1, 1, 1/5, 1/5, 1)))
  
  list(joint_optim)
}
simu_results <- mclapply(1:nrow(simu_info), estimate_from_simulation, 
                         simulations = simulations, locs = locs,
                         simu_info = simu_info,
                         mc.cores = 12, mc.preschedule = F)

save(simu_results, simu_info, file = 'results/simu_results_both.RData')
