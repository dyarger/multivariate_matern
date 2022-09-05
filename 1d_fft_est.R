
norm_constant <- function(nu_1, nu_2, a_1 = 1, a_2 = 1) {
  (a_1)^(nu_1) * (a_2)^(nu_2) *
    sqrt(gamma(nu_1 + 1/2)) * sqrt(gamma(nu_2 + 1/2))/pi^(1/2)/sqrt(gamma(nu_1)*gamma(nu_2))
}


fft_1d <- function(n_points = 2^10, x_max = 10, nu1 = .5, nu2 = .5, a1 = 1, a2 = 1, 
                   re, im) {
  delta_t <-  pi / x_max
  x_vals <- 1:n_points * (2 * pi) / n_points / delta_t
  freq_points <- seq(-delta_t * n_points/2 + delta_t/2, delta_t * n_points/2 - 
                       delta_t/2, length.out = n_points)
  # https://stackoverflow.com/questions/24077913/discretized-continuous-fourier-transform-with-numpy
  tv <- complex(real = a1, imaginary = freq_points)^(-nu1 - .5) *
    complex(real = a2, imaginary = -freq_points)^(-nu2 - .5) * 
    complex(real = re, imaginary = im*sign(freq_points))
  phase_factor <- 1/(sqrt(2*pi))  * 
    exp(complex(imaginary = freq_points[1] * 2 * pi * 
                  (1:length(freq_points)) / (delta_t * length(freq_points))))
  ff_res <- fft(tv, inverse = T)*length(tv) #* phase_factor
  #ff_res_adj <- rev((c(ff_res[freq_points >= 0], ff_res[freq_points < 0])) * phase_factor)
  ff_res_adj <- (c(ff_res[freq_points >= 0], ff_res[freq_points < 0])) * phase_factor
  x_vals <- x_vals - (x_vals[length(x_vals)] - x_vals[1])/2
  cbind(x_vals, Re(ff_res_adj) * norm_constant(nu1, nu2, a1, a2) * (5/2/x_max))
}
nu_test <- .5
df <- fft_1d(nu1 = nu_test, nu2 = nu_test, a1 = 1, a2 = 1, re = 1, im = 0, n_points = 2^16,x_max = 20)
plot(df, type = 'l')
lines(df[,1], exp(-abs(df[,1])), col = 2)
mean((abs(df[,2]) - exp(-abs(df[,1])))^2)

library(fields)
df <- fft_1d(nu1 = 1.5, nu2 = 1.5, a1 = 2, a2 = 2, re = 1, im = 0, n_points = 2^16,x_max = 5)
plot(df, type = 'l')
lines(df[,1], Matern(abs(df[,1]), range = 1/2, smoothness = 1.5), col = 2)
mean((abs(df[,2]) - Matern(abs(df[,1]), range = 1/2, smoothness = 1.5))^2)



df <- fft_1d(nu1 = 1, nu2 = 2, a1 = .8, a2 = .8, re = 1, im = 0, n_points = 2^16,x_max = 10)
plot(df, type = 'l')
df <- fft_1d(nu1 = 1, nu2 = 2, a1 = .8, a2 = .8, re = 0, im = 1, n_points = 2^16,x_max = 20)
plot(df, type = 'l')


library(R.matlab)
library(ggplot2)
library(fields)
weather <- R.matlab::readMat('bolin_code/article_code/Application/TempPress/weather_data.mat')

weather <- as.data.frame(weather)
colnames(weather) <- c('lat', 'long', 'pres', 'temp')
n <- nrow(weather)
ggplot(data = weather, aes(x = long, y = lat, color = temp)) +
  geom_point()
ggplot(data = weather, aes(x = long, y = lat, color = pres)) +
  geom_point()

dist <- fields::rdist.earth(x1 = weather[, c('long', 'lat')],
                            x2 = weather[, c('long', 'lat')], miles = F)
diag(dist) <- 0
response <- unlist(weather[, c('pres', 'temp')])
response <- response -
  rep(colMeans(weather[, c('pres', 'temp')]), each = nrow(dist))

dist_one <- dist
# dist_one <- matrix(weather[['long']], nrow = n, ncol = n) - 
#   matrix(weather[['long']], nrow = n, ncol = n, byrow = T)
# image.plot(dist_one)
# diag(dist_one)
# dist_pos <- abs(dist_one)

construct_matrix <- function(nu1, nu2, a1, a2, 
                             n_points = 2^10, 
                             x_max = 10,
                             re, im, dist_mat) {
  C_val <- fft_1d(n_points = n_points, x_max = x_max, nu1 = nu1, nu2 = nu2, 
               a1 = a1, a2 = a2, re = re, im = im)
  
  yout <- approx(xout = as.vector(dist_mat), y = C_val[,2], x = C_val[,1], method = 'linear')$y
  matrix(yout, nrow = (nrow(dist_mat)), ncol = (nrow(dist_mat)))
}
construct_entire_matrix <- function(nu1, nu2, a1, a2, 
                                    n_points = 2^10, 
                                    x_max = 10,
                                    Sigma11, Sigma12, Sigma22, dist_mat, nugget1, nugget2) {
  C1 <- construct_matrix(nu1 = nu1, nu2 = nu1, a1 = a1, a2 = a1,
                         re = Re(Sigma11), im = Im(Sigma11), n_points = n_points,
                         x_max = x_max, dist_mat)
  diag(C1) <- diag(C1) + nugget1
  C12 <- construct_matrix(nu1 = nu1, nu2 = nu2, a1 = a1, a2 = a2,
                         re = Re(Sigma12), im = Im(Sigma12), n_points = n_points,
                         x_max = x_max, dist_mat)
  C2 <- construct_matrix(nu1 = nu2, nu2 = nu2, a1 = a2, a2 = a2,
                         re = Re(Sigma22), im = Im(Sigma22), n_points = n_points,
                         x_max = x_max, dist_mat)
  diag(C2) <- diag(C2) + nugget2
  rbind(cbind(C1, C12), cbind(t(C12), C2))
  #as.matrix(Matrix::bdiag(C1, C2))
  
  # C1 <- Sigma11 * Matern(a1 * abs(dist_mat), smoothness = nu1)
  # diag(C1) <- diag(C1) + nugget1
  # C2 <- Sigma22 * Matern(a2 * abs(dist_mat), smoothness = nu2)
  # diag(C2) <- diag(C2) + nugget2
  # as.matrix(Matrix::bdiag(C1, C2))
  
  #rbind(cbind(C1, t(C12)), cbind((C12), C2))
}

test_mat <- construct_entire_matrix(nu1 = .4, nu2 = .7, a1 = 1, a2 = 1, Sigma11 = 1, Sigma22 = 1, 
                                    Sigma12 = complex(real = .4, imaginary = .4),
                                    dist_mat = dist_one, n_points = 2^20, x_max = 2000, 
                                    nugget1 = .01, nugget2 = .01)
image.plot(test_mat)
plot(diag(test_mat))
summary(Re(eigen(dist_one)$values))
solve(test_mat)
image.plot(test_mat[1:50, 1:50])

#y <- weather[, c('pres', 'temp')]
ll_fun <- function(par, dist_one, n_points = 2^20, x_max = 30, response) {
  print(exp(par[1:8]))
  #if (abs(par[9]) + abs(par[10]) > 1) {
  #  return(Inf)
  #}
  nu1 <- exp(par[1])
  nu2 <- exp(par[2])
  a1 <- exp(par[3])
  a2 <- exp(par[4])
  nugget1 <-  exp(par[5])
  nugget2 <-  exp(par[6])
  Sigma11 <- exp(par[7])
  Sigma22 <- exp(par[8])
  Sigma12re <- par[9]*sqrt(Sigma11)*sqrt(Sigma22)
  #Sigma12im <- par[10]*sqrt(Sigma11)*sqrt(Sigma22)
  Sigma12im <- 0
  cov_mat <- construct_entire_matrix(nu1 = nu1, nu2 = nu2, a1 = a1, a2 = a2, 
                                     Sigma11 = Sigma11, Sigma22 = Sigma22, 
                          Sigma12 = complex(real = Sigma12re, imaginary = Sigma12im),
                          dist_mat = dist_one, n_points = n_points, x_max = x_max, 
                          nugget1 = nugget1, nugget2 = nugget2)
  #cov_mat <- Matrix::symmpart(cov_mat)
  inv_all <- solve(cov_mat)
  ll_val <- (- nrow(cov_mat)/2 * log(2 * pi) - 1/2 * as.double(t(response) %*% inv_all %*% response) - 
      1/2 * as.vector(determinant(cov_mat)$modulus))
  #ll_val <- mvtnorm::dmvnorm(x = response, sigma = cov_mat, log = T)
  print(ll_val)
  
  - ll_val
}

test_optim <- optim(#par = c(log(c(.5, .5, .1, .1, 50, 1, var(weather$pres), var(weather$temp))), c(0, 0)),
                    par = c(log(c(1.5, 1.5, .1, .1, 50, 1, var(weather$pres), var(weather$temp))), .01),
                    fn = ll_fun, dist_one = dist_one, response = response, 
                    lower = log(c(.01, .01, 10^-6, 10^-6, NA, NA, NA, NA), -1),
                    upper = log(c(5, 5, NA, NA, NA, NA, NA, NA),1),
                    method = 'L-BFGS-B',
                    n_points = 2^15, x_max = 2200,
                    control = list(parscale = c(rep(1, 8), .1))
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
