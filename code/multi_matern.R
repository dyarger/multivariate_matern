
source('code/multi_matern_source.R')
A <- matrix(nrow = 2, ncol = 2,
            complex(real = rnorm(4), 
                    imaginary = rnorm(4)))
# A <- matrix(nrow = 2, ncol = 2,
#             complex(real = c(1,0,0,1)))
# A <- matrix(nrow = 2, ncol = 2,
#             complex(real = c(1,0,0,1)))
A <- matrix(nrow = 2, ncol = 2,
            complex(real = c(1.1,1,1,1),
                    imaginary = c(.4,.2, .1, 0)))
A <- matrix(nrow = 2, ncol = 2,
            complex(real = c(-1,1,1,0),
                    imaginary = c(1,2,-1,1)))
(A %*% t(Conj(A)))
AA_star <- A %*% t(Conj(A))
AA_star <- matrix(nrow = 2, ncol = 2, complex(real = c(1,0,0,1),
                                              imaginary = c(0, -.95, .95, 0)))
# AA_star <- matrix(nrow = 2, ncol = 2, complex(real = c(1,.9,.9,1),
#                                               imaginary = c(0, 0, 0, 0)))
# AA_star <- matrix(nrow = 2, ncol = 2,
#                   complex(real = c(1,0,0,1),
#                           imaginary = c(0, 10000, 10000,0)))
nu <- .8

#### Plot covariance
# note - plot does not show reversability for nu = 1/2, 3/2, ...
plot_cov(AA_star = AA_star, nu = nu)

test_seq <- seq(-2, 2, by = .005)
grid <- expand.grid(test_seq, test_seq)
cov_val1 <- sapply(1:nrow(grid), function(x) {
  cross_cov(grid$Var1[x], grid$Var2[x], nu = nu, z_ij = AA_star[1,1])
})
cov_val2 <- sapply(1:nrow(grid), function(x) {
  cross_cov(grid$Var1[x], grid$Var2[x], nu = nu, z_ij = AA_star[2,2])
})
cov_val12 <- sapply(1:nrow(grid), function(x) {
  cross_cov(grid$Var1[x], grid$Var2[x], nu = nu, z_ij = AA_star[1,2])
})

cov_mat1 <- matrix(cov_val1, nrow = length(test_seq))
cov_mat2 <- matrix(cov_val2, nrow = length(test_seq))
cov_mat12 <- matrix(cov_val12, nrow = length(test_seq))
cov_mat_all <- rbind(cbind(cov_mat1, cov_mat12),
                     cbind(t(cov_mat12), cov_mat2))

cov_mat_chol <- chol(cov_mat_all)


n_simu <- 3
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
ggsave('example_simulation_asymm.png', height = 5, width = 4)



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
simulated <- t(cov_mat_chol) %*% matrix(nrow = nrow(cov_mat_all),
                                        ncol = n_simu,
                                        rnorm(nrow(cov_mat_all) *n_simu))
s1 <- simulated[1:length(test_seq),]
s2 <- simulated[-(1:length(test_seq)),]
library(tidyverse)
theme_set(theme_bw())
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
ggsave('example_simulation.png', height = 5, width = 4)


ggplot(data = simulation, aes(x = var1, y = var2))+
  geom_path()+
  facet_wrap(~simulation, ncol = 1) + 
  theme(strip.background = element_blank(),
        strip.text.x = element_blank())
  
plot(test_seq, s1, type = 'l', ylim = range(c(s1, s2)))
lines(test_seq, s2, type = 'l')
plot(s1, s2, type = 'l')


#### Simulate in 1 dimension of time ###
# t_eval - points to evaluate
# N - number of basis function approximation
N <- 15000
data <- sim_bivariate(AA_star = AA_star, nu = nu, t_eval = seq(-5, 5, by = .01),
                      N = N)
par(mfrow = c(2,1))
plot(data[,1], data[,2], type = 'l', main = 'Simulation, Variable 1')
plot(data[,1], data[,3], type = 'l', main = 'Simulation, Variable 2')

# #### Compare with matern generated by RandomFields package - limited to non-reversible
# library(RandomFields)
# x <- seq(-5,5, by =.01)
# sigmaone <- cross_cov(0, 0, nu = nu, z_ij = AA_star[1,1])
# sigmatwo <- cross_cov(0, 0, nu = nu, z_ij = AA_star[2,2])
# rho <- cross_cov(0, 0, nu = nu, z_ij = AA_star[1,2])/ sqrt(sigmaone * sigmatwo)
# model <- RMbiwm(nu = c(nu,nu,nu), cdiag = Re(c(sigmaone, sigmatwo)), 
#                 rhored  = Re(rho))# consider only non-reversible part
# z <- as.data.frame(RFsimulate(model=model, x,spConform=FALSE))
# plot(x,z[,1], type = 'l')
# plot(x,z[,2], type = 'l')

#### Turning bands method for multiple dimensions in space
# n_d <- 30 # number of directions
# directions <- runif(n_d, 0, pi)
# directions <- seq(from = 0, to = pi,length.out =  n_d)
# vec_dir <- as.matrix(cbind(cos(directions), sin(directions)))
# plot(vec_dir, main = 'Directions for Turning Band method')
# 
# t_eval <- as.matrix(expand.grid(seq(-5, 5, length.out = 125),seq(-5, 5, length.out = 125)))
# data <- list()
# N <- 400
# for (r in 1:n_d) {
#   data[[r]] <- directions[r]^20* sim_bivariate(AA_star = AA_star, nu = nu, t_eval = t_eval %*% vec_dir[r,],
#                            N = N)
# }
# field1 <- apply(sapply(data, function(x) {return(x[,2])}),1, sum)/sqrt(n_d)
# field2 <- apply(sapply(data, function(x) {return(x[,3])}),1, sum)/sqrt(n_d)
# library(fields)
# par(mfrow = c(2,1))
# image.plot(matrix(field1, nrow =sqrt(nrow(t_eval))), axes=F)
# image.plot(matrix(field2, nrow =sqrt(nrow(t_eval))), axes=F)


# profvis::profvis({t_eval <- as.matrix(expand.grid(seq(-5, 5, length.out = 250),seq(-5, 5, length.out = 250)))
# data <- list()
# N <- 100
# for (r in 1:n_d) {
#   data[[r]] <- sim_bivariate(AA_star = AA_star, nu = nu, t_eval = t_eval %*% vec_dir[r,],
#                              N = N)
# }})
