

source('code/multi_matern_source.R')
AA_star <- matrix(nrow = 2, ncol = 2, complex(real = c(1,0,0,1),
                                              #imaginary = c(0, -.5,.5, 0)),
                                              imaginary = c(0,-.5,.5, 0))
                  )

nu <- .7

#### Plot covariance
# note - plot does not show reversability for nu = 1/2, 3/2, ...
plot_cov(AA_star = AA_star, nu = nu)

#### Simulate in 1 dimension of time ###
# t_eval - points to evaluate
# N - number of basis function approximation
N <- 15000
gap <- .005
n_samples <- 1200
data <- sim_bivariate(AA_star = AA_star, nu = nu, t_eval = seq(0, 2*gap*n_samples, gap),
                      N = N)
range_vals <- data[,2:3]- matrix(nrow = nrow(data), ncol = 2, byrow =2, as.double(data[1,2:3]))
plot(data[,1], data[,2] - data[1,2], type = 'l', ylim =range(range_vals))
lines(data[,1]- .5, data[,3]- data[1,3],col = 2)

plot(data[,2], data[,3], type = 'l')



theta <- c(Re(AA_star[1,1]), Re(AA_star[2,2]), Im(AA_star[1,2]))
lag_mat <- toeplitz(seq(0, 2*gap*n_samples, gap))
lag_mat[lower.tri(lag_mat)] <- -lag_mat[lower.tri(lag_mat)]
lag_vec <- seq(-2* gap*n_samples,  2*gap*n_samples, gap)


# theta <- c(1,1,.8)
# at <- matrix(nrow = 2, ncol = 2, complex(real = c(1,0,0,1), 
#                                          imaginary = c(0, .8, .8, 0)))
# corpcor::is.positive.definite(at)
likelihood <- function(theta, data) {
  print(theta)
  if (theta[1] < 0 | theta[2] < 0 | 
      abs(theta[3]) > sqrt(theta[1]) * sqrt(theta[2])) {
    return(1000)
  }
  cov_val_lag <- sapply(1:length(lag_vec), function(x) {
    c(cross_cov(0, lag_vec[x], nu = nu, z_ij = theta[1]),
      cross_cov(0, lag_vec[x], nu = nu, z_ij = complex(imaginary = -theta[3])),
      cross_cov(0, lag_vec[x], nu = nu, z_ij = theta[2])
    )
  })
  
  cov_mat_lag_1 <- sapply(as.vector(lag_mat), FUN = function(x) { cov_val_lag[1,abs(lag_vec - x) < .1*gap]} )
  cov_mat_lag_1 <- matrix(cov_mat_lag_1, nrow = sqrt(length(cov_mat_lag_1)), byrow = T)

  cov_mat_lag_2 <- sapply(as.vector(lag_mat), FUN = function(x) { cov_val_lag[3,abs(lag_vec - x) < .1*gap]} )
  cov_mat_lag_2 <- matrix(cov_mat_lag_2, nrow = sqrt(length(cov_mat_lag_2)), byrow = T)
  
  cov_mat_lag_12 <- sapply(as.vector(lag_mat), FUN = function(x) { cov_val_lag[2,abs(lag_vec - x) < .1*gap]} )
  cov_mat_lag_12 <- matrix(cov_mat_lag_12, nrow = sqrt(length(cov_mat_lag_12)), byrow = T)
  c_mat <- cbind(rbind(cov_mat_lag_1, cov_mat_lag_12),
                 rbind(t(cov_mat_lag_12), cov_mat_lag_2))
  c_chol <- chol(c_mat)
  response <- c(data$f1, data$f2)
  v2 <- .Internal(backsolve(r = c_chol, x = response, 
                            k = ncol(c_chol), upper.tri = T, transpose = T))
  quad_form <- sum(v2^2)
  #quad_form <- sum(solve(c_chol, response)^2)
  det_val <-  2* sum(log(diag(c_chol)))
  #sum(log(prod(diag(c_chol))))
  
  print( quad_form+ det_val)
  quad_form+ det_val
}



test <- optim(par = c(1,1,-.1), fn = likelihood, data = data, method = 'L-BFGS-B', 
              hessian = T, lower = c(0.01, 0.01, -Inf))
test$par  
fisherI <- solve(test$hessian)
sqrt(diag(fisherI))

plot(c(1,2,3), test$par, ylim = c(-1.3, 1.3))
points(c(1,2,3), test$par + 2*sqrt(diag(fisherI)), col = 2, pch = 5)
points(c(1,2,3), test$par - 2*sqrt(diag(fisherI)), col = 2, pch = 5)
points(c(1,2,3), c(1,1,.2), cex = .5)



# optimize using each of them independently
# theta = nu, 1, 2, 3, 
likelihood_w_smoothness <- function(theta, data) {
  print(theta)
  if (theta[1] < 0 |  theta[2] < 0 | theta[3] < 0 | 
      abs(theta[4]) > sqrt(theta[2]) * sqrt(theta[3])) {
    return(1000)
  }
  cov_val_lag <- sapply(1:length(lag_vec), function(x) {
    c(cross_cov(0, lag_vec[x], nu = theta[1], z_ij = theta[2]),
      cross_cov(0, lag_vec[x], nu = theta[1], z_ij = complex(imaginary = -theta[4])),
      cross_cov(0, lag_vec[x], nu = theta[1], z_ij = theta[3])
    )
  })
  
  cov_mat_lag_1 <- sapply(as.vector(lag_mat), FUN = function(x) { cov_val_lag[1,abs(lag_vec - x) < .1*gap]} )
  cov_mat_lag_1 <- matrix(cov_mat_lag_1, nrow = sqrt(length(cov_mat_lag_1)), byrow = T)
  
  cov_mat_lag_2 <- sapply(as.vector(lag_mat), FUN = function(x) { cov_val_lag[3,abs(lag_vec - x) < .1*gap]} )
  cov_mat_lag_2 <- matrix(cov_mat_lag_2, nrow = sqrt(length(cov_mat_lag_2)), byrow = T)
  
  cov_mat_lag_12 <- sapply(as.vector(lag_mat), FUN = function(x) { cov_val_lag[2,abs(lag_vec - x) < .1*gap]} )
  cov_mat_lag_12 <- matrix(cov_mat_lag_12, nrow = sqrt(length(cov_mat_lag_12)), byrow = T)
  c_mat <- cbind(rbind(cov_mat_lag_1, cov_mat_lag_12),
                 rbind(t(cov_mat_lag_12), cov_mat_lag_2))
  c_chol <- chol(c_mat)
  response <- c(data$f1, data$f2)
  v2 <- .Internal(backsolve(r = c_chol, x = response, 
                            k = ncol(c_chol), upper.tri = T, transpose = T))
  quad_form <- sum(v2^2)
  #quad_form <- sum(solve(c_chol, response)^2)
  det_val <-  2* sum(log(diag(c_chol)))
  #sum(log(prod(diag(c_chol))))
  
  print( quad_form+ det_val)
  quad_form+ det_val
}
test <- optim(par = c(1,1,1,-.1), fn = likelihood_w_smoothness, data = data, method = 'L-BFGS-B', 
              hessian = T, lower = c(0.01, 0.01, 0.01, -Inf))
fisherI <- solve(test$hessian)
plot(c(1,2,3,4), test$par, ylim = c(-0, 2.5))
points(c(1,2,3,4), test$par + 2*sqrt(diag(fisherI)), col = 2, pch = 5)
points(c(1,2,3,4), test$par - 2*sqrt(diag(fisherI)), col = 2, pch = 5)
points(c(1,2,3,4), c(.7, 1,1,.2), cex = .5)

# likelihood with both real and imaginary parts 
AA_star_both <- matrix(complex(real = c(1,.1,.1, 1), c(0, .1, .1, 0)), nrow = 2, ncol = 2)

data <- sim_bivariate(AA_star = AA_star_both, nu = nu, t_eval = seq(0, 2*gap*n_samples, gap),
                      N = N)
likelihood_w_smoothness_real <- function(theta, data) {
  print(theta)
  if (theta[1] < 0 |  theta[2] < 0 | theta[3] < 0 | 
      abs(theta[4]) + abs(theta[5]) > sqrt(theta[2]) * sqrt(theta[3])) {
    return(1000)
  }
  cov_val_lag <- sapply(1:length(lag_vec), function(x) {
    c(cross_cov(0, lag_vec[x], nu = theta[1], z_ij = theta[2]),
      cross_cov(0, lag_vec[x], nu = theta[1], z_ij = complex(real = theta[4], 
                                                             imaginary = -theta[5])),
      cross_cov(0, lag_vec[x], nu = theta[1], z_ij = theta[3])
    )
  })
  
  cov_mat_lag_1 <- sapply(as.vector(lag_mat), FUN = function(x) { cov_val_lag[1,abs(lag_vec - x) < .1*gap]} )
  cov_mat_lag_1 <- matrix(cov_mat_lag_1, nrow = sqrt(length(cov_mat_lag_1)), byrow = T)
  
  cov_mat_lag_2 <- sapply(as.vector(lag_mat), FUN = function(x) { cov_val_lag[3,abs(lag_vec - x) < .1*gap]} )
  cov_mat_lag_2 <- matrix(cov_mat_lag_2, nrow = sqrt(length(cov_mat_lag_2)), byrow = T)
  
  cov_mat_lag_12 <- sapply(as.vector(lag_mat), FUN = function(x) { cov_val_lag[2,abs(lag_vec - x) < .1*gap]} )
  cov_mat_lag_12 <- matrix(cov_mat_lag_12, nrow = sqrt(length(cov_mat_lag_12)), byrow = T)
  c_mat <- cbind(rbind(cov_mat_lag_1, cov_mat_lag_12),
                 rbind(t(cov_mat_lag_12), cov_mat_lag_2))
  c_chol <- chol(c_mat)
  response <- c(data$f1, data$f2)
  v2 <- .Internal(backsolve(r = c_chol, x = response, 
                            k = ncol(c_chol), upper.tri = T, transpose = T))
  quad_form <- sum(v2^2)
  det_val <-  2* sum(log(diag(c_chol)))

  print( quad_form+ det_val)
  quad_form+ det_val
}
test <- optim(par = c(1,1,1,-.1,-.1), fn = likelihood_w_smoothness, data = data, method = 'L-BFGS-B', 
              hessian = T, lower = c(0.01, 0.01, 0.01, -Inf, -Inf))
fisherI <- solve(test$hessian)
plot(c(1,2,3,4), test$par, ylim = c(-0, 2.5))
points(c(1,2,3,4), test$par + 2*sqrt(diag(fisherI)), col = 2, pch = 5)
points(c(1,2,3,4), test$par - 2*sqrt(diag(fisherI)), col = 2, pch = 5)
points(c(1,2,3,4), c(.7, 1,1,.2), cex = .5)

# make a function that plots the empirical cross-correlation

##### LIKELIHOOD with scale parameters, both real and imaginary parts

# check likelihood for one of the parameters
likelihood_single_change <- function(theta, data) {
  #print(theta)
  if (1 < 0 | 1 < 0 | 
      abs(theta) > sqrt(1) * sqrt(1)/3) {
    return(1000)
  }
  cov_val_lag <- sapply(1:length(lag_vec), function(x) {
    c(cross_cov(0, lag_vec[x], nu = nu, z_ij = 1),
      cross_cov(0, lag_vec[x], nu = nu, z_ij = complex(imaginary = -theta)),
      cross_cov(0, lag_vec[x], nu = nu, z_ij = 1)
    )
  })
  
  cov_mat_lag_1 <- sapply(as.vector(lag_mat), FUN = function(x) { cov_val_lag[1,abs(lag_vec - x) < .1*gap]} )
  cov_mat_lag_1 <- matrix(cov_mat_lag_1, nrow = sqrt(length(cov_mat_lag_1)), byrow = T)
  
  cov_mat_lag_2 <- sapply(as.vector(lag_mat), FUN = function(x) { cov_val_lag[3,abs(lag_vec - x) < .1*gap]} )
  cov_mat_lag_2 <- matrix(cov_mat_lag_2, nrow = sqrt(length(cov_mat_lag_2)), byrow = T)
  
  cov_mat_lag_12 <- sapply(as.vector(lag_mat), FUN = function(x) { cov_val_lag[2,abs(lag_vec - x) < .1*gap]} )
  cov_mat_lag_12 <- matrix(cov_mat_lag_12, nrow = sqrt(length(cov_mat_lag_12)), byrow = T)
  c_mat <- cbind(rbind(cov_mat_lag_1, cov_mat_lag_12),
                 rbind(t(cov_mat_lag_12), cov_mat_lag_2))
  c_chol <- chol(c_mat)
  response <- c(data$f1, data$f2)
  v2 <- .Internal(backsolve(r = c_chol, x = response, 
                            k = ncol(c_chol), upper.tri = T, transpose = T))
  quad_form <- sum(v2^2)
  #quad_form <- sum(solve(c_chol, response)^2)
  det_val <-  2* sum(log(diag(c_chol)))
  #sum(log(prod(diag(c_chol))))
  
  print( quad_form+ det_val)
  quad_form+ det_val
}
one_dim <- sapply(seq(-.9, .9, by = .01), likelihood_single, data = data)
plot(seq(-.9, .9, by = .01)[one_dim < 1000], one_dim[one_dim < 1000], type = 'l')


