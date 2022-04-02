
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
  cov_mat <- whittaker_covariance_matrix_lags(locs, 
                                         nu1 = exp(theta[4]), 
                                         nu2 = exp(theta[5]),
                                         c11 = exp(theta[1]), 
                                         c2 = exp(theta[2]), 
                                         c12 = theta[3], a1 = 1,
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
  cov_mat <- whittaker_covariance_matrix_lags(locs, 
                                              nu1 = exp(estimated_params[3]), 
                                              nu2 = exp(estimated_params[4]),
                                              c11 = exp(estimated_params[1]), 
                                              c2 = exp(estimated_params[2]), 
                                              c12 = theta, a1 = 1,
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
  cov_mat <- whittaker_covariance_matrix_lags(locs, 
                                              nu1 = exp(theta[1]), 
                                              nu2 = exp(theta[2]),
                                              c11 = exp(estimated_params[1]), 
                                              c2 = exp(estimated_params[2]), 
                                              c12 = exp(estimated_params[3]), a1 = 1,
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
  
  cov_mat <- whittaker_covariance_matrix_lags(locs, nu1 = .4, nu2 = .8, c11 = 1, 
                                              c2 = 1, c12 = .8, a1 = 1, a2 = 1)
  cov_mat_chol <- base::chol(cov_mat)
  simulations <- simulate_manually(cholesky = cov_mat_chol, n_simu = 1, locs)
  
  response <- c(simulations[['var1']], simulations[['var2']])
  
  joint_optim <- optim(par = c(log(1.3),log(.7), .2, log(.51),log(.51)), fn = likelihood,
                       response = response, method = 'L-BFGS-B', 
                       hessian = T, lower = c(-1, -1, -2,-2, -2),
                       upper = c(1, 1, 2,.25, .25), locs = locs)
  single_optim1 <- optim(par = c(log(1.3),log(.51)), fn = likelihood_single_var, 
                        response = response, type = 1,
                        locs = locs, method = 'L-BFGS-B', hessian = T, 
                        lower = c(-1, -2), upper = c(1,.25))
  single_optim2 <- optim(par = c(log(1.3),log(.51)), fn = likelihood_single_var, 
                        response = response, type = 2,
                        locs = locs, method = 'L-BFGS-B', hessian = T, 
                        lower = c(-1, -2), upper = c(1,.25))
  
  joint_optim_init <- optim(par = c(single_optim1$par[1],single_optim2$par[1], .2,
                               single_optim1$par[2],single_optim2$par[2]), fn = likelihood,
                       response = response, method = 'L-BFGS-B', 
                       hessian = T, lower = c(-1, -1, -2,-2, -2),
                       upper = c(1, 1, 2,.25, .25), locs = locs)
  
  joint_optim_seq <- optim(par = .2, fn = likelihood_joint_only, estimated_params = 
                             c(single_optim1$par[1],single_optim2$par[1],
                               single_optim1$par[2],single_optim2$par[2]),
                            response = response, method = 'L-BFGS-B', 
                            hessian = T, lower = -2,
                            upper =2, locs = locs)
  
  AA_star_mom <- crossprod(cbind(simulations[['var1']], simulations[['var2']]))/n
  
  joint_optim_mom <- optim(par = c(log(.51),log(.51)), fn = likelihood_nu_only, estimated_params = 
                             c(AA_star_mom[1,1], AA_star_mom[2,2], AA_star_mom[1,2]),
                           response = response, method = 'L-BFGS-B', 
                           hessian = T, lower = c(-2,-2),
                           upper = c(.25,.25), locs = locs)
  
  list(joint_optim, single_optim1, single_optim2, joint_optim_init, joint_optim_seq, 
       AA_star_mom, joint_optim_mom)
}
simu_results <- mclapply(1:nrow(simu_info), estimate_from_simulation, 
                         simulations = simulations, locs = locs,
                         simu_info = simu_info,
                         mc.cores = 12, mc.preschedule = F)

save(simu_results, simu_info, file = 'results/simu_results.RData')

q()


# simulation_long <- simulations %>%
#   pivot_longer(cols = starts_with('var'))
# 
# ggplot(data = simulation_long %>% filter(simulation ==1), 
#        aes(x = t, y = value, color = name))+
#   geom_line() + geom_point()

load('results/simu_results.RData')

simu_info <- simu_info[sapply(simu_results, length)!=1,]
simu_results <- simu_results[sapply(simu_results, length)!=1]
df_list <- list()
theta_vals <- matrix(nrow = length(simu_results), ncol = 5)
theta_vals_no_trans <- theta_vals
hessian_diag <- theta_vals
theta_single_vals <- matrix(nrow = length(simu_results), ncol = 2)
hessian_single_diag <- theta_single_vals
theta_single2_vals <- theta_single_vals
for (i in 1:length(simu_results)) {
  
  
  df1 <- data.frame(simu = i, parameter = c('Sigma11', 'Sigma22', 'Sigma12',
                                           'nu1', 'nu2'),
                   type = 'joint',
                   value = c(exp(simu_results[[i]][[1]]$par[1:2]),
                              simu_results[[i]][[1]]$par[3],
                              exp(simu_results[[i]][[1]]$par[4:5])),
                   n = simu_info[['n']][i])
  
  df2 <- data.frame(simu = i, parameter = c('Sigma11', 'Sigma22', 'Sigma12',
                                           'nu1', 'nu2'),
                   type = 'single',
                   value = c(exp(simu_results[[i]][[2]]$par[1]),
                              exp(simu_results[[i]][[3]]$par[1]),
                              simu_results[[i]][[5]]$par,
                              exp(simu_results[[i]][[2]]$par[2]),
                              exp(simu_results[[i]][[3]]$par[2])),
                   n = simu_info[['n']][i])
  
  df3 <- data.frame(simu = i, parameter = c('Sigma11', 'Sigma22', 'Sigma12',
                                           'nu1', 'nu2'),
                   type = 'init',
                   value = c(exp(simu_results[[i]][[4]]$par[1]),
                              exp(simu_results[[i]][[4]]$par[2]),
                              simu_results[[i]][[4]]$par[3],
                              exp(simu_results[[i]][[4]]$par[4]),
                              exp(simu_results[[i]][[4]]$par[5])),
                   n = simu_info[['n']][i])
  
  
  theta_single_vals[i,] <- exp(simu_results[[i]][[2]]$par)
  theta_single2_vals[i,] <- exp(simu_results[[i]][[3]]$par)
 # hessian_diag[i,] <- sqrt(diag(solve(simu_results[[i]][[1]]$hessian)))
 # hessian_single_diag[i,] <- sqrt(diag(solve(simu_results[[i]][[2]]$hessian)))
  df_list[[i]] <- rbind(df1, df2, df3)
}

true_values <- data.frame(parameter = c('Sigma11', 'Sigma22', 'Sigma12',
                                        'nu1', 'nu2'), true_value = c(1, 1, .8, .4, .8))
plot_labels <- data.frame(parameter = c('Sigma11', 'Sigma22', 'Sigma12',
                                        'nu1', 'nu2'), 
                          plot_labels = factor(c('Sigma[11]', 'Sigma[22]', 'Sigma[12]','nu[1]', 'nu[2]'),
                                               levels = c('Sigma[11]', 'Sigma[22]', 'Sigma[12]','nu[1]', 'nu[2]')))
strat_labels <- data.frame(type = c('init', 'joint', 'single'), 
                          strat_lab = c('c)', 'a)', 'b)'))


library(dplyr)
df <- dplyr::bind_rows(df_list) %>%
  left_join(true_values) %>% left_join(plot_labels) %>% left_join(strat_labels)

error_summary <- df %>%
  group_by(plot_labels, n, strat_lab) %>%
  summarise(rmse = sqrt(mean((value - true_value)^2)),
            bias = (mean(value) - mean(true_value))^2,
            var = var(value)) %>%
  as.data.frame()
library(ggplot2)
ggplot(data = error_summary, aes(x = n, y = rmse, color = strat_lab, linetype = strat_lab))+
  geom_point()+
  geom_line() + 
  facet_wrap(~plot_labels, scales = 'free_y', labeller = label_parsed) + 
  labs(x = 'Sample Size', y = 'Root-Mean-Squared-Error of Parameters',
       color = 'Estimation\nStrategy', linetype = 'Estimation\nStrategy')+
  theme_bw() + 
  theme(legend.position = 'bottom') 
ggsave('simulation_plot_real.png', height=  3.5, width = 6)

ggplot(data = error_summary, aes(x = n, y = bias, color = strat_lab, linetype = strat_lab))+
  geom_point()+
  geom_line() + 
  facet_wrap(~plot_labels, scales = 'free_y', labeller = label_parsed) + 
  labs(x = 'Sample Size', y = 'Root-Mean-Squared-Error of Parameters',
       color = 'Estimation\nStrategy', linetype = 'Estimation\nStrategy')+
  theme(legend.position = 'bottom')

ggplot(data = error_summary, aes(x = n, y = var, color = strat_lab, linetype = strat_lab))+
  geom_point()+
  geom_line() + 
  facet_wrap(~plot_labels, scales = 'free_y', labeller = label_parsed) + 
  labs(x = 'Sample Size', y = 'Root-Mean-Squared-Error of Parameters',
       color = 'Estimation\nStrategy', linetype = 'Estimation\nStrategy')+
  theme(legend.position = 'bottom')

error_summary


colMeans(hessian_diag)
colMeans(hessian_single_diag)





sapply(split.data.frame(abs(theta_vals_no_trans -  matrix(vec_params_no_trans, nrow = nrow(theta_vals),
                 ncol = ncol(theta_vals), byrow = T)) <  2 * hessian_diag, simu_info$n),
       colMeans)
colMeans(abs(theta_vals_no_trans[,c(1,4)] -  matrix(vec_params_no_trans[c(1,4)], nrow = nrow(hessian_single_diag),
                                           ncol = ncol(hessian_single_diag), byrow = T)) <  2 *hessian_single_diag)
sapply(split.data.frame(hessian_diag, simu_info$n), colMeans)
sapply(split.data.frame(hessian_single_diag, simu_info$n), colMeans)
