library(tidyverse)
##### simulation with both imaginary and real ####
load('results/simu_Results_both.RData')

df_list <- list()
for (i in 1:length(simu_results)) {
  
  
  df1 <- data.frame(simu = i, parameter = c('Sigma11', 'Sigma22', 'Re_Sigma12',
                                            'Im_Sigma12',
                                            'nu'),
                    type = 'joint',
                    value = c(exp(simu_results[[i]][[1]]$par[1:2]),
                              simu_results[[i]][[1]]$par[3:4],
                              exp(simu_results[[i]][[1]]$par[5])),
                    n = simu_info[['n']][i])
  
  
  
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

df <- dplyr::bind_rows(df_list) %>%
  left_join(true_values) %>% left_join(plot_labels) %>% left_join(strat_labels)

error_summary <- df %>%
  group_by(plot_labels, n, strat_lab) %>%
  summarise(rmse = sqrt(mean((value - true_value)^2)),
            bias = (mean(value) - mean(true_value))^2,
            var = var(value)) %>%
  as.data.frame()

error_summary <- df %>%
  group_by(plot_labels, n, strat_lab) %>%
  summarise(rmse = sqrt(mean((value - true_value)^2)),
            bias = (mean(value) - mean(true_value)),
            var = sd(value)) %>%
  as.data.frame()

error_summary



###### from simulation with real entries #####

load('results/simu_results_real.RData')
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


df <- dplyr::bind_rows(df_list) %>%
  left_join(true_values) %>% left_join(plot_labels) %>% left_join(strat_labels)

error_summary <- df %>%
  group_by(plot_labels, n, strat_lab) %>%
  summarise(rmse = sqrt(mean((value - true_value)^2)),
            bias = (mean(value) - mean(true_value))^2,
            var = var(value)) %>%
  as.data.frame()
ggplot(data = error_summary, aes(x = n, y = rmse, color = strat_lab, linetype = strat_lab))+
  geom_point()+
  geom_line() + 
  facet_wrap(~plot_labels, scales = 'free_y', labeller = label_parsed) + 
  labs(x = 'Sample Size', y = 'Root-Mean-Squared-Error of Parameters',
       color = 'Estimation\nStrategy', linetype = 'Estimation\nStrategy')+
  theme_bw() + 
  theme(legend.position = 'bottom') 
ggsave('simulation_plot_real.png', height = 3.5, width = 6)


#### from simulation with imaginary entries ######

load('results/simu_results_imaginary.RData')
df_list <- list()
for (i in 1:length(simu_results)) {
  
  
  df1 <- data.frame(simu = i, parameter = c('Sigma11', 'Sigma22', 'Sigma12',
                                            'nu1'),
                    type = 'joint',
                    value = c(exp(simu_results[[i]][[1]]$par[1:2]),
                              simu_results[[i]][[1]]$par[3],
                              exp(simu_results[[i]][[1]]$par[4])),
                    n = simu_info[['n']][i])
  
  df2 <- data.frame(simu = i, parameter = c('Sigma11', 'Sigma22', 'Sigma12',
                                            'nu1'),
                    type = 'single',
                    value = c(exp(simu_results[[i]][[2]]$par[1]),
                              exp(simu_results[[i]][[3]]$par[1]),
                              simu_results[[i]][[5]]$par,
                              mean(c(exp(simu_results[[i]][[2]]$par[2]),
                                     exp(simu_results[[i]][[3]]$par[2])))),
                    n = simu_info[['n']][i])
  
  df3 <- data.frame(simu = i, parameter = c('Sigma11', 'Sigma22', 'Sigma12',
                                            'nu1'),
                    type = 'init',
                    value = c(exp(simu_results[[i]][[4]]$par[1]),
                              exp(simu_results[[i]][[4]]$par[2]),
                              simu_results[[i]][[4]]$par[3],
                              exp(simu_results[[i]][[4]]$par[4])),
                    n = simu_info[['n']][i])
  df_list[[i]] <- rbind(df1, df2, df3)
}

true_values <- data.frame(parameter = c('Sigma11', 'Sigma22', 'Sigma12',
                                        'nu1', 'nu2'), true_value = c(1, 1, .8, .8, .8))
plot_labels <- data.frame(parameter = c('Sigma11', 'Sigma22', 'Sigma12',
                                        'nu1', 'nu2'), 
                          plot_labels = factor(c('Sigma[11]', 'Sigma[22]', 'Im(Sigma[12])','nu[1]', 'nu[2]'),
                                               levels = c('Sigma[11]', 'Sigma[22]', 'Im(Sigma[12])','nu[1]', 'nu[2]')))
strat_labels <- data.frame(type = c('init', 'joint', 'single'), 
                           strat_lab = c('c)', 'a)', 'b)'))

df <- dplyr::bind_rows(df_list) %>%
  left_join(true_values) %>% left_join(plot_labels) %>% left_join(strat_labels)

error_summary <- df %>%
  group_by(plot_labels, n, strat_lab) %>%
  summarise(rmse = sqrt(mean((value - true_value)^2)),
            bias = (mean(value) - mean(true_value))^2,
            var = var(value)) %>%
  as.data.frame()
ggplot(data = error_summary, aes(x = n, y = rmse, color = strat_lab, linetype = strat_lab))+
  geom_point()+
  geom_line() + 
  facet_wrap(~plot_labels, scales = 'free_y', labeller = label_parsed) + 
  labs(x = 'Sample Size', y = 'Root-Mean-Squared-Error of Parameters',
       color = 'Estimation\nStrategy', linetype = 'Estimation\nStrategy')+
  theme_bw() + 
  theme(legend.position = 'bottom') 
ggsave('simulation_plot_imaginary.png', height=  3.5, width = 6)

