library(tidyverse)
theme_set(theme_bw())
source('code/multi_matern_source.R')
# AA_star <- matrix(nrow = 2, ncol = 2, complex(real = c(1,0,0,1),
#                                               imaginary = c(0, -.95, .95, 0)))
AA_star <- matrix(nrow = 2, ncol = 2, complex(real = c(1,0,0,1),
                                              imaginary = c(0, -.95, .95, 0)))
nu <- .8

#### Plot covariance
# gap <- .005
# range <- 2
gap <- .01
range <- 2
test_seq <- seq(-range, range, by = gap)
#test_seq <- seq(-2, 2, by = .05)
grid <- expand.grid(test_seq, test_seq)
grid$lag <- round((grid$Var1 - grid$Var2)/gap)*gap
grid_unique <- grid[!duplicated(grid$lag),]
cov_val1 <- sapply(1:nrow(grid_unique), function(x) {
  full_cross_cov_single(h = grid_unique$lag[x], nu = nu, a = 1, 
                        realp = Re(AA_star[1,1]), imp = Im(AA_star[1,1]))
})
cov_val2 <- sapply(1:nrow(grid_unique), function(x) {
  full_cross_cov_single(h = grid_unique$lag[x], nu = nu, a = 1, 
                        realp = Re(AA_star[2,2]), imp = Im(AA_star[2,2]))
})
cov_val12 <- sapply(1:nrow(grid_unique), function(x) {
  full_cross_cov_single(h =  grid_unique$lag[x], nu = nu, a = 1, 
                        realp = Re(AA_star[1,2]), imp = Im(AA_star[1,2]))
})
grid_unique$cov1 <- cov_val1
grid_unique$cov2 <- cov_val2
grid_unique$cov12 <- cov_val12

grid_all <- full_join(grid, grid_unique[,c('lag','cov1', 'cov2', 'cov12')], 
                      by = c('lag'))
grid_all <- grid_all[order(grid_all$Var2, grid_all$Var1),]
cov_mat1 <- matrix(grid_all$cov1, nrow = length(test_seq))
cov_mat2 <- matrix(grid_all$cov2, nrow = length(test_seq))
cov_mat12 <- matrix(grid_all$cov12, nrow = length(test_seq))
cov_mat_all <- rbind(cbind(cov_mat1, cov_mat12),
                     cbind(t(cov_mat12), cov_mat2))
# cov_mat_all <- rbind(cbind(cov_mat1, t(cov_mat12)),
#                      cbind(cov_mat12, cov_mat2))

cov_mat_chol <- chol(cov_mat_all)
n_simu <- 3
set.seed(104)
simulated <- t(cov_mat_chol) %*% matrix(nrow = nrow(cov_mat_all),
                                        ncol = n_simu,
                                        rnorm(nrow(cov_mat_all) *n_simu))
s1 <- simulated[1:length(test_seq),]
s2 <- simulated[-(1:length(test_seq)),]
cov(simulated)
cor(simulated)

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
  scale_color_discrete(labels = c('j', 'k'))
ggsave('images/example_simulation_asymm.png', height = 5, width = 4)

cc_summary <- grid_all %>%
  left_join(data.frame('Var1' = test_seq, p1l1 = s1[,1], p2l1 = s2[,1]), by = 'Var1') %>%
  left_join(data.frame('Var2' = test_seq, p1l2 = s1[,1], p2l2 = s2[,1]), by = 'Var2') %>%
  group_by(lag) %>%
  summarise(cc_val = mean(p1l1 * p2l2),
            vval1 = mean(p1l1 * p1l2),
            vval2 = mean(p2l1*p2l2)) %>%
  filter(abs(lag) < 2)
plot(cc_summary$lag, cc_summary$vval1, cex = .2)
points(cc_summary$lag, sapply(cc_summary$lag, function(x) {
  full_cross_cov_single(h = x, nu = nu, a = 1, 
                        realp = Re(AA_star[1,1]), imp = Im(AA_star[1,1]))
}),col = 2, cex = .2)
plot(cc_summary$lag, cc_summary$vval2, cex = .2)
points(cc_summary$lag, sapply(cc_summary$lag, function(x) {
  full_cross_cov_single(h = x, nu = nu, a = 1, 
                        realp = Re(AA_star[2,2]), imp = Im(AA_star[2,2]))
}),col = 2, cex = .2)
plot(cc_summary$lag, cc_summary$cc_val, cex = .2, ylim = c(-.4, .4))
points(cc_summary$lag, sapply(cc_summary$lag, function(x) {
  full_cross_cov_single(h = x, nu = nu, a = 1, 
                        realp = Re(AA_star[1,2]), imp = Im(AA_star[1,2]))
}),col = 2, cex = .2)
# try whittaker 
nu1 <- .4
nu2 <- .8
AA_star <- matrix(nrow = 2, ncol = 2, complex(real = c(1,.95,.95,1),
                                              imaginary = c(0, 0, 0, 0))) 
cov_val1 <- sapply(1:nrow(grid_unique), function(x) {
  whitt_only_single(h = grid_unique$lag[x], nu1 = nu1, nu2 = nu1, a1 = 1, a2 =1, 
                    realp = Re(AA_star[1,1]), imp = Im(AA_star[1,1]), which_val = 1)
})
cov_val2 <- sapply(1:nrow(grid_unique), function(x) {
  whitt_only_single(h = grid_unique$lag[x], nu1 = nu2, nu2 = nu2, a1 = 1, a2 =1, 
                    realp = Re(AA_star[2,2]), imp = Im(AA_star[2,2]), which_val = 3)
})
cov_val12 <- sapply(1:nrow(grid_unique), function(x) {
  whitt_only_single(h =  grid_unique$lag[x], nu1 = nu1, nu2 = nu2, a1 = 1, a2 =1, 
                    realp = Re(AA_star[1,2]), imp = Im(AA_star[1,2]), which_val = 2)
})
grid_unique$cov1 <- cov_val1
grid_unique$cov2 <- cov_val2
grid_unique$cov12 <- cov_val12

grid_all <- merge(grid, grid_unique[,c('lag','cov1', 'cov2', 'cov12')], 
                  by = c('lag'), all = T)
grid_all <- grid_all[order(grid_all$Var2, grid_all$Var1),]
cov_mat1 <- matrix(grid_all$cov1, nrow = length(test_seq))
cov_mat2 <- matrix(grid_all$cov2, nrow = length(test_seq))
cov_mat12 <- matrix(grid_all$cov12, nrow = length(test_seq))

cov_mat_all <- rbind(cbind(cov_mat1, cov_mat12),
                     cbind(t(cov_mat12), cov_mat2))

# cov_mat_all <- rbind(cbind(cov_mat1, t(cov_mat12)),
#                      cbind(cov_mat12, cov_mat2))

cov_mat_chol <- chol(cov_mat_all)

n_simu <- 3
set.seed(25)
simulated <- t(cov_mat_chol) %*% matrix(nrow = nrow(cov_mat_all),
                                        ncol = n_simu,
                                        rnorm(nrow(cov_mat_all) *n_simu))
s1 <- simulated[1:length(test_seq),]
s2 <- simulated[-(1:length(test_seq)),]
cov(simulated)
cor(simulated)

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
  scale_color_discrete(labels = c('j', 'k'))
ggsave('images/example_simulation.png', height = 5, width = 4)


cc_summary <- grid_all %>%
  left_join(data.frame('Var1' = test_seq, p1l1 = s1[,1], p2l1 = s2[,1]), by = 'Var1') %>%
  left_join(data.frame('Var2' = test_seq, p1l2 = s1[,1], p2l2 = s2[,1]), by = 'Var2') %>%
  group_by(lag) %>%
  summarise(cc_val = mean(p1l1 * p2l2),
            vval1 = mean(p1l1 * p1l2),
            vval2 = mean(p2l1 * p2l2)) %>%
  filter(abs(lag) < 1)

plot(cc_summary$lag, cc_summary$vval1, cex = .2)
points(cc_summary$lag, sapply(cc_summary$lag, function(x) {
  whitt_only_single(x, nu1 = nu1, nu2 = nu2, a1 = 1, a2 = 1, realp = Re(AA_star[1,1]),
                    imp = Im(AA_star[1,1]), norm_type = 'D', which_val = 1)
}),col = 2, cex = .2)
plot(cc_summary$lag, cc_summary$vval2, cex = .2)
points(cc_summary$lag, sapply(cc_summary$lag, function(x) {
  whitt_only_single(x, nu1 = nu1, nu2 = nu2, a1 = 1, a2 = 1, realp = Re(AA_star[2,2]),
                    imp = Im(AA_star[2,2]), norm_type = 'D', which_val = 1)
}),col = 2, cex = .2)
plot(cc_summary$lag, cc_summary$cc_val, cex = .2, ylim = c(-.1, 1))
points(cc_summary$lag, sapply(cc_summary$lag, function(x) {
  whitt_only_single(x, nu1 = nu1, nu2 = nu2, a1 = 1, a2 = 1, realp = Re(AA_star[1,2]),
                    imp = Im(AA_star[1,2]), norm_type = 'D', which_val = 2)
}),col = 2, cex = .2)

#### Simulate in 1 dimension of time ###
# t_eval - points to evaluate
# N - number of basis function approximation
N <- 15000
data <- sim_bivariate(AA_star = AA_star, nu = nu, t_eval = test_seq,
                      N = N)
#par(mfrow = c(2,1))
plot(data[,1], data[,2], type = 'l', main = 'Simulation, Variable 1')
lines(col = 2,data[,1], data[,3], type = 'l', main = 'Simulation, Variable 2')
AA_star <- matrix(nrow = 2, ncol = 2, complex(real = c(1,0,0,1),
                                              imaginary = c(0, -.95, .95, 0)))
nu <- 1.2

data <- sim_bivariate(AA_star = AA_star, nu = nu, t_eval = test_seq,
                      N = N)
#par(mfrow = c(2,1))
plot(data[,1], data[,2], type = 'l', main = 'Simulation, Variable 1')
lines(col = 2,data[,1], data[,3], type = 'l', main = 'Simulation, Variable 2')
