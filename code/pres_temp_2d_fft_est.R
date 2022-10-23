
library(Rcpp)
library(tidyverse)
library(fftw)
library(fftwtools)
library(R.matlab)
library(fields)

source('code/multi_matern_source.R')
# load in data
weather <- R.matlab::readMat('data/weather_data.mat')

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

# create long/lat distance matrices
dist_one <- dist
dist_tens <- array(NA, c(dim(dist), 2))
reference_location <- colMeans(weather[, c('long', 'lat')])
for (i in 1:dim(dist_tens)[1]) {
  for (j in 1:dim(dist_tens)[1]) {
    dist_long_prelim <- fields::rdist.earth.vec(cbind(weather$long[i],
                                                      reference_location['lat']), 
                                                cbind(weather$long[j],
                                                      reference_location['lat']),
                                                miles = F)
    if (weather$long[i] < weather$long[j]) {
      dist_long_prelim <- -dist_long_prelim
    }
    
    dist_lat_prelim <- fields::rdist.earth.vec(cbind(reference_location['long'],
                                                     weather$lat[i]), 
                                                cbind(reference_location['long'],
                                                      weather$lat[j]),
                                                miles = F)
    if (weather$lat[i] < weather$lat[j]) {
      dist_lat_prelim <- -dist_lat_prelim
    }
    
    dist_tens[i,j,1] <- dist_long_prelim
    dist_tens[i,j,2] <- dist_lat_prelim
  }
}
dist_tens_mat <- cbind(as.vector(dist_tens[,,1]), as.vector(dist_tens[,,2]))

# given distances, compute fft on grid, then interpolate onto distances
construct_matrix <- function(nu1, nu2, a1, a2, 
                             grid_info,
                             Psi, Delta, dist_tens_mat) {
  C_val <- fft_2d(grid_info = grid_info, nu1 = nu1, nu2 = nu2, 
                  a1 = a1, a2 = a2,  Psi = Psi, Delta = Delta)
  
  tmat <- matrix(C_val[['val']], nrow = sqrt(nrow(C_val)))
  long1 <- grid_info[['x_vals']]
  lat1 <- grid_info[['x_vals']]
  
  interp_points <- fields::interp.surface(obj = list(x = long1, y = lat1, z = tmat),
                                          loc = dist_tens_mat)
  matrix(interp_points, nrow = sqrt(length(interp_points)))
}

# construct bivariate Matern covariance matrix
construct_bivariate_matrix <- function(nu1, nu2, a1, a2, 
                                    grid_info,
                                    Psi_list, Delta_list, dist_tens_mat, nugget1, nugget2) {
  C1 <- construct_matrix(nu1 = nu1, nu2 = nu1, a1 = a1, a2 = a1,
                         Psi = Psi_list[[1]], Delta = Delta_list[[1]], 
                         grid_info = grid_info, dist_tens_mat)
  diag(C1) <- diag(C1) + nugget1
  C12 <- construct_matrix(nu1 = nu1, nu2 = nu2, a1 = a1, a2 = a2,
                          Psi = Psi_list[[3]], Delta = Delta_list[[3]], 
                          grid_info = grid_info, dist_tens_mat)
  C2 <- construct_matrix(nu1 = nu2, nu2 = nu2, a1 = a2, a2 = a2,
                         Psi = Psi_list[[2]], Delta = Delta_list[[2]],
                         grid_info = grid_info, dist_tens_mat)
  diag(C2) <- diag(C2) + nugget2
  rbind(cbind(C1, C12), cbind(t(C12), C2))
}

grid_info <- create_grid_info_2d(2^10, x_max = 2500)

# test the matrix creation
test_mat <- construct_bivariate_matrix(nu1 = .5, nu2 = .5, a1 = 10^-3, a2 = 10^-3, 
                                    Delta_list = list(function(x) 1, function(x) 1, 
                                                 function(x) {
                                                   .97 * complex(imaginary = sign(x))
                                                 }),
                                    Psi_list = replicate(3, function(x) {sign(x)}),
                                    grid_info = grid_info,
                                    dist_tens_mat = dist_tens_mat, 
                                    nugget1 = .01, nugget2 = .01)
test_mat_true <- Matern(d = abs(dist), range = 10^3, smoothness = .5)
summary(as.vector(abs(test_mat_true - test_mat[1:157, 1:157])))
image.plot(test_mat)
plot(diag(test_mat))
summary(Re(eigen(test_mat)$values))
solve(test_mat)
image.plot(test_mat[1:50, 1:50])
isSymmetric(test_mat)


nu1_init <- nu2_init <- .5
nu_lower <- 0.001; nu_higher = 7
a1_init <- a2_init <- 10^-3
a_lower <- 10^-6
Sigma11_init <- var(weather[['pres']])
Sigma22_init <- var(weather[['temp']])
cor_re_init <- cor(weather[['temp']], weather[['pres']])
cor_im_init <- cor(weather[['temp']], weather[['pres']])/10
nugget1_init <- var(weather[['pres']])/25
nugget2_init <- var(weather[['temp']])/25
var_lower <- 10^-9

# log likelihood as function of parameters
ll_fun <- function(par, dist_tens_mat, response, grid_info) {
  #print(exp(par[1:8]))
  nu1 <- exp(par[1]); nu2 <- exp(par[2])
  a1 <- exp(par[3]); a2 <- exp(par[4])
  nugget1 <- exp(par[5]); nugget2 <- exp(par[6])
  Sigma11 <- exp(par[7]); Sigma22 <- exp(par[8])
  Sigma12 <- par[9]*sqrt(Sigma11)*sqrt(Sigma22)
  cov_mat <- construct_bivariate_matrix(nu1 = nu1, nu2 = nu2, a1 = a1, a2 = a2, 
                                     Delta_list = list(function(x) Sigma11, function(x) Sigma22, 
                                                       function(x) Sigma12),
                                     Psi_list = replicate(3, {function(x) {
                                       sign(x)
                                     }}),
                                     dist_tens_mat = dist_tens_mat, grid_info = grid_info,
                                     nugget1 = nugget1, nugget2 = nugget2)
  chol_mat <- base::chol(cov_mat)
  ll_val <- -nrow(cov_mat)/2 * log(2 * pi) - 1/2 * sum(backsolve(chol_mat, response, transpose = T)^2) - 
      sum(log(diag(chol_mat)))
  #print(ll_val)
  - ll_val
}


#### with psi fixed at sign(theta)
test_optim <- optim(
  par = c(log(c(nu1_init, nu2_init, a1_init, a2_init, nugget1_init, nugget2_init, 
                Sigma11_init, Sigma22_init)), cor_re_init),
  fn = ll_fun, dist_tens_mat = dist_tens_mat, response = response, 
  lower = c(log(c(nu_lower, nu_lower, a_lower, a_lower, var_lower, var_lower, var_lower, var_lower)), -1),
  upper = c(log(c(nu_higher, nu_higher, NA, NA, NA, NA, NA, NA)), 1),
  method = 'L-BFGS-B',
  grid_info = grid_info,
  control = list(parscale = c(rep(1, 8), .1))
)

test_optim
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
par[9]

df_fix <- fft_2d(nu1 = nu1, nu2 = nu2, a1 = a1, a2 = a2, 
                d = 2,
                Delta = function(x) Sigma12re, Psi = function(x) {
                  sign(x)
                }, grid_info = grid_info)

ll_fun_psi_real <- function(par, dist_tens_mat, response, grid_info) {
  #print(exp(par[1:8]))
  nu1 <- exp(par[1]); nu2 <- exp(par[2])
  a1 <- exp(par[3]); a2 <- exp(par[4])
  nugget1 <- exp(par[5]); nugget2 <- exp(par[6])
  Sigma11 <- exp(par[7]); Sigma22 <- exp(par[8])
  Sigma12 <- par[9]*sqrt(Sigma11)*sqrt(Sigma22)
  theta_star <- atan2(par[11], par[10])
  Psi_fun <- function(theta) {
    if (theta_star < 0) {
      ifelse(theta > theta_star & theta < theta_star + pi, 1, -1)
    } else {
      ifelse(theta > theta_star | theta < theta_star - pi, 1, -1)
    }
  }
  cov_mat <- construct_bivariate_matrix(nu1 = nu1, nu2 = nu2, a1 = a1, a2 = a2, 
                                     Delta_list = list(function(x) Sigma11, function(x) Sigma22, 
                                                       function(x) Sigma12),
                                     Psi_list = replicate(3, Psi_fun),
                                     dist_tens_mat = dist_tens_mat, grid_info = grid_info,
                                     nugget1 = nugget1, nugget2 = nugget2)
  chol_mat <- base::chol(cov_mat)
  ll_val <- -nrow(cov_mat)/2 * log(2 * pi) - 1/2 * sum(backsolve(chol_mat, response, transpose = T)^2) - 
    sum(log(diag(chol_mat)))
  #print(ll_val)
  - ll_val
}
test_optim_real <- optim(
  par = c(test_optim$par, 0.1, 0.1),
  fn = ll_fun_psi_real, dist_tens_mat = dist_tens_mat, response = response, 
  lower = c(log(c(nu_lower, nu_lower, a_lower, a_lower, var_lower, var_lower, var_lower, var_lower)), -1, NA, NA),
  upper = c(log(c(nu_higher, nu_higher, NA, NA, NA, NA, NA, NA)), 1, NA, NA),
  method = 'L-BFGS-B',
  grid_info = grid_info,
  control = list(parscale = c(rep(1, 8), .1, .1, .1))
)
test_optim_real
par <- test_optim_real$par
(nu1 <- exp(par[1]))
(nu2 <- exp(par[2]))
(a1 <- exp(par[3]))
(a2 <- exp(par[4]))

(nugget1 <-  exp(par[5]))
(nugget2 <-  exp(par[6]))
(Sigma11 <- exp(par[7]))
(Sigma22 <- exp(par[8]))
(Sigma12re <- par[9]*sqrt(Sigma11)*sqrt(Sigma22))
par[9]
atan2(par[11],par[10])

theta_star <- atan2(par[11], par[10])
Psi_fun <- function(theta) {
  if (theta_star < 0) {
    ifelse(theta > theta_star & theta < theta_star + pi, 1, -1)
  } else {
    ifelse(theta > theta_star | theta < theta_star - pi, 1, -1)
  }
}

df_re <- fft_2d(nu1 = nu1, nu2 = nu2, a1 = a1, a2 = a2, 
                d = 2,
                Delta = function(x) Sigma12re, Psi = Psi_fun, grid_info = grid_info)

# with an imaginary entry
ll_fun_psi_im <- function(par, dist_tens_mat, response, grid_info) {
  #print(exp(par[1:8]))
  #print((par[9:10]))
  nu1 <- exp(par[1]); nu2 <- exp(par[2])
  a1 <- exp(par[3]); a2 <- exp(par[4])
  nugget1 <- exp(par[5]); nugget2 <- exp(par[6])
  Sigma11 <- exp(par[7]); Sigma22 <- exp(par[8])
  Sigma12 <- par[9]*sqrt(Sigma11)*sqrt(Sigma22)
  Sigma12im <- par[10]*sqrt(Sigma11)*sqrt(Sigma22)
  theta_star <- atan2(par[12], par[11])
  theta_star2 <- atan2(par[14], par[13])
  Psi_fun <- function(theta) {
    if (theta_star < 0) {
      ifelse(theta > theta_star & theta < theta_star + pi, 1, -1)
    } else {
      ifelse(theta > theta_star | theta < theta_star - pi, 1, -1)
    }
  }
  Psi_fun2 <- function(theta) {
    if (theta_star2 < 0) {
      ifelse(theta > theta_star2 & theta < theta_star2 + pi, 1, -1)
    } else {
      ifelse(theta > theta_star2 | theta < theta_star2 - pi, 1, -1)
    }
  }
  cov_mat <- construct_bivariate_matrix(nu1 = nu1, nu2 = nu2, a1 = a1, a2 = a2, 
                                     Delta_list = list(function(x) Sigma11, 
                                                       function(x) Sigma22, 
                                                       function(x) complex(real = Sigma12,
                                                                           imaginary = Sigma12im*Psi_fun2(x))),
                                     Psi_list = replicate(3, Psi_fun),
                                     dist_tens_mat = dist_tens_mat, grid_info = grid_info,
                                     nugget1 = nugget1, nugget2 = nugget2)
  chol_mat <- base::chol(cov_mat)
  ll_val <- -nrow(cov_mat)/2 * log(2 * pi) - 1/2 * sum(backsolve(chol_mat, response, transpose = T)^2) - 
    sum(log(diag(chol_mat)))
  #print(ll_val)
  - ll_val
}

test_optim_im <- optim(
  par = c(test_optim_real$par[1:9],0, test_optim_real$par[10:11],.1,.1),
  fn = ll_fun_psi_im, dist_tens_mat = dist_tens_mat, response = response, 
  lower = c(log(c(nu_lower, nu_lower, a_lower, a_lower, var_lower, var_lower, var_lower, var_lower)), -1, -1, NA, NA, NA, NA),
  upper = c(log(c(nu_higher, nu_higher, NA, NA, NA, NA, NA, NA)), 1, 1, NA, NA, NA, NA),
  method = 'L-BFGS-B',
  grid_info = grid_info,
  control = list(parscale = c(rep(1, 8), .1, .1, .1, .1, .1, .1))
)
test_optim_im
par <- test_optim_im$par
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
atan2(par[12],par[11])

theta_star <- atan2(par[12], par[11])
Psi_fun <- function(theta) {
  if (theta_star < 0) {
    ifelse(theta > theta_star & theta < theta_star + pi, 1, -1)
  } else {
    ifelse(theta > theta_star | theta < theta_star - pi, 1, -1)
  }
}
theta_star2 <- atan2(par[14], par[13])
Psi_fun2 <- function(theta) {
  if (theta_star2 < 0) {
    ifelse(theta > theta_star2 & theta < theta_star2 + pi, 1, -1)
  } else {
    ifelse(theta > theta_star2 | theta < theta_star2 - pi, 1, -1)
  }
}


df_im <- fft_2d(nu1 = nu1, nu2 = nu2, a1 = a1, a2 = a2, 
                      d = 2,
                    Delta = function(x) complex(real = Sigma12re,
                                                imaginary = Sigma12im*Psi_fun2(x)), Psi = Psi_fun, grid_info = grid_info)

# only an imaginary entry, not presented
ll_fun_psi_im_only <- function(par, dist_tens_mat, response, grid_info) {
  #print(exp(par[1:8]))
 # print((par[9:10]))
  nu1 <- exp(par[1]); nu2 <- exp(par[2])
  a1 <- exp(par[3]); a2 <- exp(par[4])
  nugget1 <- exp(par[5]); nugget2 <- exp(par[6])
  Sigma11 <- exp(par[7]); Sigma22 <- exp(par[8])
  Sigma12im <- par[9]*sqrt(Sigma11)*sqrt(Sigma22)
  theta_star <- atan2(par[11], par[10])
  theta_star2 <- atan2(par[13], par[12])
  Psi_fun <- function(theta) {
    if (theta_star < 0) {
      ifelse(theta > theta_star & theta < theta_star + pi, 1, -1)
    } else {
      ifelse(theta > theta_star | theta < theta_star - pi, 1, -1)
    }
  }
  Psi_fun2 <- function(theta) {
    if (theta_star2 < 0) {
      ifelse(theta > theta_star2 & theta < theta_star2 + pi, 1, -1)
    } else {
      ifelse(theta > theta_star2 | theta < theta_star2 - pi, 1, -1)
    }
  }
  cov_mat <- construct_bivariate_matrix(nu1 = nu1, nu2 = nu2, a1 = a1, a2 = a2, 
                                     Delta_list = list(function(x) Sigma11, 
                                                       function(x) Sigma22, 
                                                       function(x) complex(real = 0,
                                                                           imaginary = Sigma12im*Psi_fun2(x))),
                                     Psi_list = replicate(3, Psi_fun),
                                     dist_tens_mat = dist_tens_mat, grid_info = grid_info,
                                     nugget1 = nugget1, nugget2 = nugget2)
  chol_mat <- base::chol(cov_mat)
  ll_val <- -nrow(cov_mat)/2 * log(2 * pi) - 1/2 * sum(backsolve(chol_mat, response, transpose = T)^2) - 
    sum(log(diag(chol_mat)))
  #print(ll_val)
  - ll_val
}

test_optim_im_only <- optim(
  par = c(test_optim_real$par[1:8],0, test_optim_real$par[10:11],.1,.1),
  fn = ll_fun_psi_im_only, dist_tens_mat = dist_tens_mat, response = response, 
  lower = c(log(c(nu_lower, nu_lower, a_lower, a_lower, var_lower, var_lower, var_lower, var_lower)), -1, NA, NA, NA, NA),
  upper = c(log(c(nu_higher, nu_higher, NA, NA, NA, NA, NA, NA)), 1, NA, NA, NA, NA),
  method = 'L-BFGS-B',
  grid_info = grid_info,
  control = list(parscale = c(rep(1, 8), .1, .1, .1, .1, .1))
)
test_optim_im_only
par <- test_optim_im_only$par
(nu1 <- exp(par[1]))
(nu2 <- exp(par[2]))
(a1 <- exp(par[3]))
(a2 <- exp(par[4]))

(nugget1 <-  exp(par[5]))
(nugget2 <-  exp(par[6]))
(Sigma11 <- exp(par[7]))
(Sigma22 <- exp(par[8]))
(Sigma12im <- par[9]*sqrt(Sigma11)*sqrt(Sigma22))
par[9]
atan2(par[11],par[10])

theta_star <- atan2(par[11], par[10])
Psi_fun <- function(theta) {
  if (theta_star < 0) {
    ifelse(theta > theta_star & theta < theta_star + pi, 1, -1)
  } else {
    ifelse(theta > theta_star | theta < theta_star - pi, 1, -1)
  }
}
theta_star2 <- atan2(par[13], par[12])
Psi_fun2 <- function(theta) {
  if (theta_star2 < 0) {
    ifelse(theta > theta_star2 & theta < theta_star2 + pi, 1, -1)
  } else {
    ifelse(theta > theta_star2 | theta < theta_star2 - pi, 1, -1)
  }
}


df_im_only <- fft_2d(nu1 = nu1, nu2 = nu2, a1 = a1, a2 = a2, 
                d = 2,
                Delta = function(x) complex(real = 0,
                                            imaginary = Sigma12im*Psi_fun2(x)), Psi = Psi_fun, grid_info = grid_info)


# independent Matern
construct_bivariate_matrix_independent <- function(nu1, nu2, a1, a2, 
                                    grid_info,
                                    Psi_list, Delta_list, dist_tens_mat, nugget1, nugget2) {
  C1 <- construct_matrix(nu1 = nu1, nu2 = nu1, a1 = a1, a2 = a1,
                         Psi = Psi_list[[1]], Delta = Delta_list[[1]], 
                         grid_info = grid_info, dist_tens_mat)
  diag(C1) <- diag(C1) + nugget1
  C2 <- construct_matrix(nu1 = nu2, nu2 = nu2, a1 = a2, a2 = a2,
                         Psi = Psi_list[[2]], Delta = Delta_list[[2]],
                         grid_info = grid_info, dist_tens_mat)
  diag(C2) <- diag(C2) + nugget2
  as.matrix(Matrix::bdiag(C1, C2))
}

ll_fun_ind <- function(par, dist_tens_mat, response, grid_info) {
  #print(exp(par[1:8]))
  nu1 <- exp(par[1]); nu2 <- exp(par[2])
  a1 <- exp(par[3]); a2 <- exp(par[4])
  nugget1 <- exp(par[5]); nugget2 <- exp(par[6])
  Sigma11 <- exp(par[7]); Sigma22 <- exp(par[8])
  cov_mat <- construct_bivariate_matrix_independent(nu1 = nu1, nu2 = nu2, a1 = a1, a2 = a2, 
                                     Delta_list = list(function(x) Sigma11, 
                                                       function(x) Sigma22, 
                                                       function(x) 0),
                                     Psi_list = replicate(3, function(x) sign(x)),
                                     dist_tens_mat = dist_tens_mat, grid_info = grid_info,
                                     nugget1 = nugget1, nugget2 = nugget2)
  chol_mat <- base::chol(cov_mat)
  ll_val <- -nrow(cov_mat)/2 * log(2 * pi) - 1/2 * sum(backsolve(chol_mat, response, transpose = T)^2) - 
    sum(log(diag(chol_mat)))
  #print(ll_val)
  - ll_val
}

test_optim_ind <- optim(
  par = c(log(c(nu1_init, nu2_init, a1_init, a2_init, nugget1_init, nugget2_init, 
                Sigma11_init, Sigma22_init))),
  fn = ll_fun_ind, dist_tens_mat = dist_tens_mat, response = response, 
  lower = c(log(c(nu_lower, nu_lower, a_lower, a_lower, var_lower, var_lower, var_lower, var_lower))),
  upper = c(log(c(nu_higher, nu_higher, NA, NA, NA, NA, NA, NA))),
  method = 'L-BFGS-B',
  grid_info = grid_info,
  control = list(parscale = rep(1, 8))
)
test_optim_ind
par <- test_optim_ind$par
(nu1 <- exp(par[1]))
(nu2 <- exp(par[2]))
(a1 <- exp(par[3]))
(a2 <- exp(par[4]))

(nugget1 <-  exp(par[5]))
(nugget2 <-  exp(par[6]))
(Sigma11 <- exp(par[7]))
(Sigma22 <- exp(par[8]))




# single correlation function
# independent Matern
ll_fun_single <- function(par, dist_tens_mat, response, grid_info) {
  #print(exp(par[1:7]))
  nu1 <- exp(par[1]); a1 <- exp(par[2])
  nugget1 <- exp(par[3]); nugget2 <- exp(par[4])
  Sigma11 <- exp(par[5]); Sigma22 <- exp(par[6])
  Sigma12 <- par[7]*sqrt(Sigma11)*sqrt(Sigma22)
  cov_mat <- construct_bivariate_matrix(nu1 = nu1, nu2 = nu1, a1 = a1, a2 = a1, 
                                     Delta_list = list(function(x) Sigma11, 
                                                       function(x) Sigma22, 
                                                       function(x) Sigma12),
                                     Psi_list = replicate(3, function(x) sign(x)),
                                     dist_tens_mat = dist_tens_mat, grid_info = grid_info,
                                     nugget1 = nugget1, nugget2 = nugget2)
  chol_mat <- base::chol(cov_mat)
  ll_val <- -nrow(cov_mat)/2 * log(2 * pi) - 1/2 * sum(backsolve(chol_mat, response, transpose = T)^2) - 
    sum(log(diag(chol_mat)))
  #print(ll_val)
  
  - ll_val
}

test_optim_single <- optim(
  par = c(log(c(nu1_init, a1_init, nugget1_init,nugget2_init ,
                Sigma11_init, Sigma22_init)), cor_re_init),
  fn = ll_fun_single, dist_tens_mat = dist_tens_mat, response = response, 
  lower = c(log(c(nu_lower, a_lower, var_lower, var_lower, var_lower, var_lower)), -1),
  upper = c(log(c(nu_higher, NA, NA, NA, NA, NA)), 1),
  method = 'L-BFGS-B',
  grid_info = grid_info,
  control = list(parscale = c(rep(1, 6), .1))
)
test_optim_single
par <- test_optim_single$par
(nu1 <- exp(par[1]))
(a1 <- exp(par[2]))

(nugget1 <-  exp(par[3]))
(nugget2 <-  exp(par[4]))
(Sigma11 <- exp(par[5]))
(Sigma22 <- exp(par[6]))
par[7]

df_single <- fft_2d(nu1 = nu1, nu2 = nu1, a1 = a1, a2 = a1, Delta = 
                  function(x) par[7] * sqrt(Sigma11) * sqrt(Sigma22), Psi = function(x) sign(x), grid_info = grid_info)


# bivariate Matern
construct_matrix_mm <- function(nu1, nu2, nu12, a1, a2, a12,
                                    grid_info,
                                    Psi_list, Delta_list, dist_tens_mat, nugget1, nugget2) {
  C1 <- construct_matrix(nu1 = nu1, nu2 = nu1, a1 = a1, a2 = a1,
                         Psi = Psi_list[[1]], Delta = Delta_list[[1]], 
                         grid_info = grid_info, dist_tens_mat)
  diag(C1) <- diag(C1) + nugget1
  C12 <- construct_matrix(nu1 = nu12, nu2 = nu12, a1 = a12, a2 = a12,
                          Psi = Psi_list[[3]], Delta = Delta_list[[3]], 
                          grid_info = grid_info, dist_tens_mat)
  C2 <- construct_matrix(nu1 = nu2, nu2 = nu2, a1 = a2, a2 = a2,
                         Psi = Psi_list[[2]], Delta = Delta_list[[2]],
                         grid_info = grid_info, dist_tens_mat)
  diag(C2) <- diag(C2) + nugget2
  rbind(cbind(C1, C12), cbind(t(C12), C2))
}

# log likelihood as function of parameters
ll_fun_mm <- function(par, dist_tens_mat, response, grid_info) {
  #print(exp(par[1:10]))
  nu1 <- exp(par[1]); nu2 <- exp(par[2]); nu12 <- exp(par[3])
  a1 <- exp(par[4]); a2 <- exp(par[5]); a12 <- exp(par[6])
  nugget1 <- exp(par[7]); nugget2 <- exp(par[8])
  Sigma11 <- exp(par[9]); Sigma22 <- exp(par[10])
  Sigma12 <- par[11]*sqrt(Sigma11)*sqrt(Sigma22)
  cov_mat <- construct_matrix_mm(nu1 = nu1, nu2 = nu2, a1 = a1, a2 = a2, 
                                 nu12 = nu12, a12 = a12,
                                     Delta_list = list(function(x) Sigma11, function(x) Sigma22, 
                                                       function(x) Sigma12),
                                     Psi_list = replicate(3, {function(x) {
                                       sign(x)
                                     }}),
                                     dist_tens_mat = dist_tens_mat, grid_info = grid_info,
                                     nugget1 = nugget1, nugget2 = nugget2)
  if (min(eigen(cov_mat)$values) < .Machine$double.eps) {
    #print(NA)
    return(8000)
  }
  chol_mat <- base::chol(cov_mat)
  ll_val <- -nrow(cov_mat)/2 * log(2 * pi) - 1/2 * sum(backsolve(chol_mat, response, transpose = T)^2) - 
    sum(log(diag(chol_mat)))
  #print(ll_val)
  
  - ll_val
}

test_optim_mm <- optim(
  par = c(log(c(nu1_init, nu2_init, (nu1_init + nu2_init)/2, 
                a1_init, a2_init,  (a1_init + a1_init)/2,
                nugget1_init, nugget2_init,
                Sigma11_init, Sigma22_init)), cor_re_init),
  fn = ll_fun_mm, dist_tens_mat = dist_tens_mat, response = response, 
  lower = c(log(c(nu_lower,nu_lower, nu_lower,  a_lower,a_lower, a_lower,
                  var_lower, var_lower, var_lower, var_lower)), -1),
  upper = c(log(c(nu_higher, nu_higher, nu_higher, 1, 1, 1, NA, NA, NA, NA)), 1),
  method = 'L-BFGS-B',
  grid_info = grid_info,
  control = list(parscale = c(rep(1, 10), .1))
)

test_optim_mm
par <- test_optim_mm$par
(nu1 <- exp(par[1]))
(nu2 <- exp(par[2]))
(nu12 <- exp(par[3]))
(a1 <- exp(par[4]))
(a2 <- exp(par[5]))
(a12 <- exp(par[6]))

(nugget1 <-  exp(par[7]))
(nugget2 <-  exp(par[8]))
(Sigma11 <- exp(par[9]))
(Sigma22 <- exp(par[10]))
par[11]

df_mm <- fft_2d(nu1 = nu12, nu2 = nu12, a1 = a12, a2 = a12, Delta = 
               function(x) par[11] * sqrt(Sigma11) * sqrt(Sigma22), 
               Psi = function(x) sign(x), grid_info = grid_info)

df_all <- rbind(cbind(df_fix, type = 'theta_star_fixed'),
                cbind(df_re, type = 'real'),
                cbind(df_im, type = 'imaginary'),
                cbind(df_im_only, type = 'imaginary_only'),
                cbind(df_mm, type = 'multi_matern'),
                cbind(df_single, type = 'single'))

df_params <- data.frame(type = c('theta_star_fixed', 'real', 'imaginary', 'multi_matern', 'single'),
                        var1 = exp(c(test_optim$par[7], test_optim_real$par[7],
                                     test_optim_im$par[7], test_optim_mm$par[9],
                                     test_optim_single$par[5])),
                        var2 = exp(c(test_optim$par[8], test_optim_real$par[8],
                                     test_optim_im$par[8], test_optim_mm$par[10],
                                     test_optim_single$par[6])))

labels <- data.frame(type = c('real', 'imaginary', 'multi_matern', 'single', 'theta_star_fixed'),
                     label = factor(c('TPMM w/ real directional measure',
                               'TPMM w/ complex directional measure',
                               'MM of Gnieting et al. (2010)',
                               'Single covariance function',
                               'TPMM w/ phi fixed'),
                               levels = c('Single covariance function',
                                          'MM of Gnieting et al. (2010)',
                                          'TPMM w/ phi fixed', 
                                          'TPMM w/ real directional measure',
                                          'TPMM w/ complex directional measure')))

ggplot(data = df_all %>%
         filter(abs(Var1) < 500, abs(Var2) < 500, type != 'imaginary_only') %>%
         left_join(labels), aes(x = Var2, y = Var1, fill = val)) +
  geom_raster() + 
  facet_wrap(~label) + 
  scale_fill_gradientn(colors = rev(rainbow(10))) +
  labs(x = 'Zonal distance (km)', y = 'Meridional distance (km)',
       fill = 'Cross-\ncovariance') +
  coord_equal() + 
  theme(legend.position = 'left',legend.key.height = unit(.8, "cm")) 
ggsave('images/cov_fun_comparison_data.png', height = 6, width = 9)

save(test_optim, test_optim_im, test_optim_ind, test_optim_single, test_optim_mm, test_optim_real,
     df_all,dist_tens_mat, df_params,
     file = 'results/pres_temp_actual_spectral_results_10_final.RData')
