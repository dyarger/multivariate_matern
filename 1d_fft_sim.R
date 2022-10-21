
norm_constant <- function(nu_1, nu_2, a_1 = 1, a_2 = 1) {
  (a_1)^(nu_1) * (a_2)^(nu_2) *
    sqrt(gamma(nu_1 + 1/2)) * sqrt(gamma(nu_2 + 1/2))/pi^(1/2)/sqrt(gamma(nu_1)*gamma(nu_2))
}

create_grid_info <- function(n_points, x_max) {
  delta_t <-  pi / x_max
  x_vals <- 1:n_points * (2 * pi) / n_points / delta_t
  freq_points <- seq(-delta_t * n_points/2 + delta_t/2, 
                     delta_t * n_points/2 - delta_t/2, 
                     length.out = n_points)
  phase_factor <- 1/(sqrt(2*pi))  * 
    exp(complex(imaginary = freq_points[1] * 2 * pi * 
                  (1:length(freq_points)) / (delta_t * length(freq_points))))
  x_vals <- x_vals - x_max - abs(abs(x_vals[length(x_vals)] - 2*x_max) - abs(x_vals[1]))
  list('freq_points' = freq_points,
       'delta_t' = delta_t, 'n_points' = n_points, 
       'x_vals' = x_vals, 'x_max' = x_max,
       'phase_factor' = phase_factor)
}

grid_info <- create_grid_info(2^10, 15)



fft_1d <- function(grid_info, nu1 = .5, nu2 = .5, a1 = 1, a2 = 1, 
                   re, im) {
  phase_factor = grid_info[['phase_factor']]
  delta_t = grid_info[['delta_t']]
  n_points = grid_info[['n_points']]
  x_vals = grid_info[['x_vals']]
  freq_points = grid_info[['freq_points']]
  x_max = grid_info[['x_max']]
  # https://stackoverflow.com/questions/24077913/discretized-continuous-fourier-transform-with-numpy
  tv <- complex(real = a1, imaginary = freq_points)^(-nu1 - .5) *
    complex(real = a2, imaginary = -freq_points)^(-nu2 - .5) * 
    complex(real = re, imaginary = im*sign(freq_points))
  ff_res <- fftwtools::fftw_c2c(data = tv, inverse = 1)
  p <- length(ff_res)/2
  ff_res_adj <- c(ff_res[(p + 1):(2*p)], ff_res[1:p]) * phase_factor
  cbind(x_vals, 'val' = 
          as.double(Re(ff_res_adj) * norm_constant(nu1, nu2, a1, a2)) /x_max * n_points * 2.512596 )
}
nu_test <- .5
df <- fft_1d(nu1 = nu_test, nu2 = nu_test, a1 = 1, a2 = 1, re = 1, im = 0, grid_info = grid_info)
plot(df, type = 'l')
lines(df[,1], exp(-abs(df[,1])), col = 2)
mean((abs(df[,2]) - exp(-abs(df[,1])))^2)
mean(exp(-abs(df[,1]))/ abs(df[,2]))

df[df[,1] == 0,2] / exp(-abs(df[df[,1] == 0,1]))
exp(-abs(df[df[,1] == 0,1]))/ df[df[,1] == 0,2] 

library(fields)
df <- fft_1d(nu1 = 1.5, nu2 = 1.5, a1 = 2, a2 = 2, re = 1, im = 0,  grid_info = grid_info)
plot(df, type = 'l')
lines(df[,1], Matern(abs(df[,1]), range = 1/2, smoothness = 1.5), col = 2)
mean((abs(df[,2]) - Matern(abs(df[,1]), range = 1/2, smoothness = 1.5))^2)



df <- fft_1d(nu1 = 1, nu2 = 2, a1 = .8, a2 = .8, re = 1, im = 0, grid_info = grid_info)
plot(df, type = 'l')
df <- fft_1d(nu1 = 1, nu2 = 2, a1 = .8, a2 = .8, re = 0, im = 1,  grid_info = grid_info)
plot(df, type = 'l')


library(R.matlab)
library(ggplot2)
library(fields)
n <- 1000
loc <- runif(n)

dist <- fields::rdist(loc)
diag(dist) <- 0

for (i in 1:nrow(dist)) {
  for (j in 1:nrow(dist)) {
    if (i == j) {
      
    } else if (loc[i] < loc[j]) {
      dist[i,j] = -dist[i,j]
    }
  }
}

construct_matrix <- function(nu1, nu2, a1, a2, 
                             grid_info,
                             re, im, dist_mat) {
  C_val <- fft_1d(grid_info = grid_info, nu1 = nu1, nu2 = nu2, 
                  a1 = a1, a2 = a2, re = re, im = im)
  
  yout <- approx(xout = as.vector(dist_mat), y = C_val[,2], x = C_val[,1], method = 'linear')$y
  matrix(yout, nrow = (nrow(dist_mat)), ncol = (nrow(dist_mat)))
}
construct_entire_matrix <- function(nu1, nu2, a1, a2, 
                                    grid_info = grid_info,
                                    Sigma11, Sigma12, Sigma22, dist_mat, nugget1, nugget2) {
  C1 <- construct_matrix(nu1 = nu1, nu2 = nu1, a1 = a1, a2 = a1,
                         re = Re(Sigma11), im = Im(Sigma11), grid_info = grid_info, dist_mat)
  diag(C1) <- diag(C1) + nugget1
  C12 <- construct_matrix(nu1 = nu1, nu2 = nu2, a1 = a1, a2 = a2,
                          re = Re(Sigma12), im = Im(Sigma12), grid_info = grid_info, dist_mat)
  C2 <- construct_matrix(nu1 = nu2, nu2 = nu2, a1 = a2, a2 = a2,
                         re = Re(Sigma22), im = Im(Sigma22), grid_info = grid_info, dist_mat)
  diag(C2) <- diag(C2) + nugget2
  rbind(cbind(C1, C12), cbind(t(C12), C2))
}

test_mat <- construct_entire_matrix(nu1 = 1.2, nu2 = 1.5, a1 = 4, a2 = 4, Sigma11 = 1, Sigma22 = 1, 
                                    Sigma12 = complex(real = .4, imaginary = .4),
                                    dist_mat = dist, grid_info = grid_info,
                                    nugget1 = .1, nugget2 = .1)
image.plot(test_mat)
plot(diag(test_mat))
summary(Re(eigen(test_mat)$values))
solve(test_mat)
image.plot(test_mat[1:50, 1:50])

response <- as.vector(mvtnorm::rmvnorm(n = 1, sigma = test_mat))

ll_fun <- function(par, dist_one, grid_info, response) {
  print(exp(par[1:8]))
  print((par[9:10]))
  if (sum(abs(par[9:10])) > 1) {
    return(8000)
  }
  
  nu1 <- exp(par[1])
  nu2 <- exp(par[2])
  a1 <- exp(par[3])
  a2 <- exp(par[4])
  nugget1 <-  exp(par[5])
  nugget2 <-  exp(par[6])
  Sigma11 <- exp(par[7])
  Sigma22 <- exp(par[8])
  Sigma12re <- par[9]*sqrt(Sigma11)*sqrt(Sigma22)
  Sigma12im <- par[10]*sqrt(Sigma11)*sqrt(Sigma22)
  cov_mat <- construct_entire_matrix(nu1 = nu1, nu2 = nu2, a1 = a1, a2 = a2, 
                                     Sigma11 = Sigma11, Sigma22 = Sigma22, 
                                     Sigma12 = complex(real = Sigma12re, imaginary = Sigma12im),
                                     dist_mat = dist_one, grid_info = grid_info,
                                     nugget1 = nugget1, nugget2 = nugget2)
  inv_all <- solve(cov_mat)
  ll_val <- (- nrow(cov_mat)/2 * log(2 * pi) - 1/2 * as.double(t(response) %*% inv_all %*% response) - 
               1/2 * as.vector(determinant(cov_mat)$modulus))
  print(ll_val)
  
  - ll_val
}

test_optim <- optim(
  par = c(log(c(1.5, 1.5, 2, 2, .2, .2, var(response[1:n]), var(response[-c(1:n)]))), .01, .01),
  fn = ll_fun, dist_one = dist, response = response, 
  lower = c(log(c(.01, .01, 10^-6, 10^-6, NA, NA, NA, NA)), -1, -1),
  upper = c(log(c(3, 3, NA, NA, NA, NA, NA, NA)),1, 1),
  method = 'L-BFGS-B',
  grid_info = grid_info,
  control = list(parscale = c(rep(1, 8), .1, .1))
)

par <- test_optim$par
(nu1 <- exp(par[1]))
(nu2 <- exp(par[2]))
(a1 <- exp(par[3]))
(a2 <- exp(par[4]))

(nugget1 <-  exp(par[5]))
(nugget2 <-  exp(par[6]))
(Sigma11 <- exp(par[7]))
(Sigma22 <- exp(par[8]))
(Sigma12re <- par[9]*sqrt(Sigma11)*sqrt(Sigma22))
(Sigma12im <- par[10]*sqrt(Sigma11)*sqrt(Sigma22))
par[9]
par[10]
