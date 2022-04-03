library(tidyverse)
theme_set(theme_bw())
source('code/multi_matern_source.R')
AA_star <- matrix(nrow = 2, ncol = 2, complex(real = c(1,0,0,1),
                                              imaginary = c(0, -.95, .95, 0)))
nu <- .8

#### Plot covariance
# note - plot does not show reversability for nu = 1/2, 3/2, ...
plot_cov(AA_star = AA_star, nu = nu)
test_seq <- seq(-2, 2, by = .005)
grid <- expand.grid(test_seq, test_seq)
cov_val1 <- sapply(1:nrow(grid), function(x) {
  cross_cov(grid$Var1[x], grid$Var2[x], nu = nu, z_ij = AA_star[1,1], a = 1)
})
cov_val2 <- sapply(1:nrow(grid), function(x) {
  cross_cov(grid$Var1[x], grid$Var2[x], nu = nu, z_ij = AA_star[2,2], a = 1)
})
cov_val12 <- sapply(1:nrow(grid), function(x) {
  cross_cov(grid$Var1[x], grid$Var2[x], nu = nu, z_ij = AA_star[1,2], a = 1)
})

cov_mat1 <- matrix(cov_val1, nrow = length(test_seq))
cov_mat2 <- matrix(cov_val2, nrow = length(test_seq))
cov_mat12 <- matrix(cov_val12, nrow = length(test_seq))
cov_mat_all <- rbind(cbind(cov_mat1, cov_mat12),
                     cbind(t(cov_mat12), cov_mat2))

cov_mat_chol <- chol(cov_mat_all)
n_simu <- 3
set.seed(25)
simulated <- t(cov_mat_chol) %*% matrix(nrow = nrow(cov_mat_all),
                                        ncol = n_simu,
                                        rnorm(nrow(cov_mat_all) *n_simu))
s1 <- simulated[1:length(test_seq),]
s2 <- simulated[-(1:length(test_seq)),]

simulation <- data.frame(var1 = as.double(s1), var2 = as.double(s2),
                         t = test_seq, 
                         simulation = rep(1:n_simu, each = length(test_seq)))
simulation_long <- simulation %>%
  tidyr::pivot_longer(cols = starts_with('var'))
ggplot(data = simulation_long, aes(x = t, y = value, group = name, color = name))+
  geom_line()+
  facet_wrap(~simulation, ncol = 1) + 
  labs(x = 't', y = 'Value', color = 'Variable') +
  theme(strip.background = element_blank(),
        strip.text.x = element_blank())+
  scale_color_discrete(labels = c('1', '2'))
ggsave('images/example_simulation_asymm.png', height = 5, width = 4)

# try whittaker 
test_seq <- seq(-2, 2, by = .004)
grid <- expand.grid(test_seq, test_seq)
nu1 <- .4
nu2 <- .8
cov_val <- sapply(1:nrow(grid), function(x) {
  whitt_version(grid$Var1[x]-grid$Var2[x], nu1= nu1, nu2 = nu2, c11 = 1, c12 = .95, c2 = 1)
})

cov_mat1 <- matrix(cov_val[1,], nrow = length(test_seq))
cov_mat2 <- matrix(cov_val[3,], nrow = length(test_seq))
cov_mat12 <- matrix(cov_val[2,], nrow = length(test_seq))
cov_mat_all <- rbind(cbind(cov_mat1, cov_mat12),
                     cbind(t(cov_mat12), cov_mat2))

cov_mat_chol <- chol(cov_mat_all)

n_simu <- 3
set.seed(25)
simulated <- t(cov_mat_chol) %*% matrix(nrow = nrow(cov_mat_all),
                                        ncol = n_simu,
                                        rnorm(nrow(cov_mat_all) *n_simu))
s1 <- simulated[1:length(test_seq),]
s2 <- simulated[-(1:length(test_seq)),]

simulation <- data.frame(var1 = as.double(s1), var2 = as.double(s2),
                              t = test_seq, 
                              simulation = rep(1:n_simu, each = length(test_seq)))
simulation_long <- simulation %>%
  pivot_longer(cols = starts_with('var'))
ggplot(data = simulation_long, aes(x = t, y = value, group = name, color = name))+
  geom_line()+
  facet_wrap(~simulation, ncol = 1) + 
  labs(x = 't', y = 'Value', color = 'Variable') +
  theme(strip.background = element_blank(),
        strip.text.x = element_blank())+
  scale_color_discrete(labels = c('1', '2'))
ggsave('images/example_simulation.png', height = 5, width = 4)

#### Simulate in 1 dimension of time ###
# t_eval - points to evaluate
# N - number of basis function approximation
# I'm not convinced this works
N <- 15000
data <- sim_bivariate(AA_star = AA_star, nu = nu, t_eval = seq(-5, 5, by = .01),
                      N = N)
par(mfrow = c(2,1))
plot(data[,1], data[,2], type = 'l', main = 'Simulation, Variable 1')
plot(data[,1], data[,3], type = 'l', main = 'Simulation, Variable 2')
