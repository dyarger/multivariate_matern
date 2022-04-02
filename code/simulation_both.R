
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
  # single_optim1 <- optim(par = c(log(1.3),log(.51)), fn = likelihood_single_var, 
  #                        response = response, type = 1,
  #                        locs = locs, method = 'L-BFGS-B', hessian = T, 
  #                        lower = c(-1, -2), upper = c(1,.25))
  # single_optim2 <- optim(par = c(log(1.3),log(.51)), fn = likelihood_single_var, 
  #                        response = response, type = 2,
  #                        locs = locs, method = 'L-BFGS-B', hessian = T, 
  #                        lower = c(-1, -2), upper = c(1,.25))
  # 
  # joint_optim_init <- optim(par = c(single_optim1$par[1],single_optim2$par[1], 0,0,
  #                                   mean(c(single_optim1$par[2],single_optim2$par[2]))), fn = likelihood,
  #                           response = response, method = 'L-BFGS-B', 
  #                           hessian = T, lower = c(-1, -1, -2,-2, -2),
  #                           upper = c(1, 1, 2,2,.25), locs = locs)
  # 
  # joint_optim_seq <- optim(par = c(0,0), fn = likelihood_joint_only, estimated_params = 
  #                            c(single_optim1$par[1],single_optim2$par[1],
  #                              mean(c(single_optim1$par[2],single_optim2$par[2]))),
  #                          response = response, method = 'L-BFGS-B', 
  #                          hessian = T, lower = c(-2, -2),
  #                          upper =c(2, 2), locs = locs)
  # 
  # AA_star_mom <- crossprod(cbind(simulations[['var1']], simulations[['var2']]))/n
  # 
  # joint_optim_mom <- optim(par = c(log(.51)), fn = likelihood_nu_only, estimated_params = 
  #                            c(AA_star_mom[1,1], AA_star_mom[2,2], AA_star_mom[1,2]),
  #                          response = response, method = 'L-BFGS-B', 
  #                          hessian = T, lower = c(-2),
  #                          upper = c(.25), locs = locs)
  
  list(joint_optim)
}
simu_results <- mclapply(1:nrow(simu_info), estimate_from_simulation, 
                         simulations = simulations, locs = locs,
                         simu_info = simu_info,
                         mc.cores = 12, mc.preschedule = F)

save(simu_results, simu_info, file = 'results/simu_results_both.RData')

q()


# simulation_long <- simulations %>%
#   pivot_longer(cols = starts_with('var'))
# 
# ggplot(data = simulation_long %>% filter(simulation ==1), 
#        aes(x = t, y = value, color = name))+
#   geom_line() + geom_point()

load('results/simu_results_both.RData')

#simu_info <- simu_info[sapply(simu_results, length)!=1,]
#simu_results <- simu_results[sapply(simu_results, length)!=1]
df_list <- list()
theta_vals <- matrix(nrow = length(simu_results), ncol = 5)
theta_vals_no_trans <- theta_vals
hessian_diag <- theta_vals
theta_single_vals <- matrix(nrow = length(simu_results), ncol = 2)
hessian_single_diag <- theta_single_vals
theta_single2_vals <- theta_single_vals
for (i in 1:length(simu_results)) {
  
  
  df1 <- data.frame(simu = i, parameter = c('Sigma11', 'Sigma22', 'Re_Sigma12',
                                            'Im_Sigma12',
                                            'nu'),
                    type = 'joint',
                    value = c(exp(simu_results[[i]][[1]]$par[1:2]),
                              simu_results[[i]][[1]]$par[3:4],
                              exp(simu_results[[i]][[1]]$par[5])),
                    n = simu_info[['n']][i])
  
 
  
  hessian_diag[i,] <- sqrt(diag(solve(simu_results[[i]][[1]]$hessian)))
  df_list[[i]] <- df1
}

true_values <- data.frame(parameter = c('Sigma11', 'Sigma22', 'Re_Sigma12',
                                        'Im_Sigma12', 'nu'), true_value = c(1, 1, .3, .3, .8))
plot_labels <- data.frame(parameter = c('Sigma11', 'Sigma22', 'Re_Sigma12',
                                        'Im_Sigma12', 'nu'), 
                          plot_labels = factor(c('Sigma[11]', 'Sigma[22]', 'ReSigma[12]','ImSigma[12]', 'nu'),
                                               levels = c('Sigma[11]', 'Sigma[22]', 'ReSigma[12]','ImSigma[12]', 'nu')))
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
ggsave('simulation_plot_both.png', height=  3.5, width = 6)



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

df %>%
  group_by(plot_labels, n, strat_lab) %>%
  summarise(rmse = sqrt(mean((value - true_value)^2)),
            bias = (mean(value) - mean(true_value)),
            var = sqrt(mean((value - mean(value))^2))) %>%
  filter(plot_labels == 'ReSigma[12]')
df %>%
  group_by(plot_labels, n, strat_lab) %>%
  summarise(rmse = sqrt(mean((value - true_value)^2)),
            bias = (mean(value) - mean(true_value)),
            var = sqrt(mean((value - mean(value))^2))) %>%
  filter(plot_labels == 'ImSigma[12]')

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
