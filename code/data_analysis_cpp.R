# cpp only

library(R.matlab)
library(ggplot2)
library(fields)
weather <- R.matlab::readMat('bolin_code/article_code/Application/TempPress/weather_data.mat')

weather <- as.data.frame(weather)
colnames(weather) <- c('lat', 'long', 'pres', 'temp')

ggplot(data = weather, aes(x = long, y = lat, color = temp))+
  geom_point()
ggplot(data = weather, aes(x = long, y = lat, color = pres))+
  geom_point()

dist <- fields::rdist.earth(x1 = weather[, c('long', 'lat')],
                            x2 = weather[, c('long', 'lat')], miles = F)
diag(dist) <- 0
response <- unlist(weather[, c('pres', 'temp')])
response <- response - 
  rep(colMeans(weather[, c('pres', 'temp')]), each = nrow(dist))


library(Rcpp)
cppFunction(depends = c('RcppArmadillo', 'RcppGSL'), 
            'arma::mat spec_dens_single_cpp(arma::mat fvals, arma::mat foutvals, double xmax, double nuval, std::size_t size) {
double* fvals_mem = fvals.memptr();
double* foutvals_mem = foutvals.memptr();
const gsl_dht *t = gsl_dht_new(size, nuval, xmax);
gsl_dht_apply(t, fvals_mem, foutvals_mem);
return foutvals;
}         ', includes = '#include <gsl/gsl_dht.h>')

cppFunction(depends = c('RcppArmadillo', 'RcppGSL'), 
            'arma::mat hankel_get_info_x(arma::mat fvals, double xmax, double nuval, std::size_t size) {
const gsl_dht *t = gsl_dht_new(size, nuval, xmax);
for (int i = 0; i < size; i++) {
  fvals[i] = gsl_dht_x_sample(t, i);
}
return fvals;
}         ', includes = '#include <gsl/gsl_dht.h>')

cppFunction(depends = c('RcppArmadillo', 'RcppGSL'), 
            'arma::mat hankel_get_info_k(arma::mat fvals, double xmax, double nuval, std::size_t size) {
const gsl_dht *t = gsl_dht_new(size, nuval, xmax);
for (int i = 0; i < size; i++) {
  fvals[i] = gsl_dht_k_sample(t, i);
}
return fvals;
}         ', includes = '#include <gsl/gsl_dht.h>')


spec_dens_cpp <- function(x, d, a_1, a_2, nu_1, nu_2) {
  y <- complex(imaginary = x, length.out = length(x))
  Re((a_1 + y)^(-nu_1-d/2)*(a_2 - y)^(-nu_2-d/2)) *
    x^(d/2 - 1)
}

ml_val <- function(par, dist, response, x_vals, k_vals, d = 2) {
  #theta: a_1, a_2, nu_1, nu_2, sigma1, sigma2, sigma12, nugget1, nuggets
  print(format(exp(par), scientific = F))
  print(par[7])
  a_1 <- exp(par[1]);a_2 <- exp(par[2])
  nu_1 <- exp(par[3]);nu_2 <- exp(par[4])
  sigma_1 <- exp(par[5]);sigma_2 <- exp(par[6])
  sigma_12 <- par[7]*sqrt(sigma_1) * sqrt(sigma_2)
  nugget_var_1 <- exp(par[8]);nugget_var_2 <- exp(par[9])
  
  
  # approximate, then interpolate
  test_input <- spec_dens_cpp(x_vals, d = d, 
                              a_1 = a_1, a_2 = a_2, nu_1 = nu_1, nu_2 = nu_2)
  cc_vals <- spec_dens_single_cpp(matrix(ncol = 1, nrow= approx_n, test_input),
                               matrix(ncol = 1, nrow= approx_n, 0),
                               x_max,
                               nu_val,
                               approx_n)
  cc_vals_new <- (2 * pi)^(d/2-1) *k_vals^(-d/2 + 1) * as.double(cc_vals)* sqrt(nu_1) * sqrt(nu_2)*2* 
    a_1^nu_1 * a_2^nu_2 
  values <- approx(xout = dist[upper.tri(dist, diag = T)], x = k_vals, 
                   y = cc_vals_new, method = 'linear', rule = 2)$y
  
  
  cov12 <- dist
  cov12[upper.tri(dist, diag = T)] <- values
  cov12[lower.tri(dist, diag = F)] <- t(cov12)[lower.tri(dist, diag = F)]
  
  cov1 <- Matern(dist/a_1, smoothness = nu_1)
  
  cov11 <- dist
  cov11[upper.tri(dist, diag = T)] <- Matern(dist[upper.tri(dist, diag = T)] * a_1, 
                                             smoothness = nu_1)
  cov11[lower.tri(dist, diag = F)] <- t(cov11)[lower.tri(dist, diag = F)]
  
  
  cov22 <- dist
  cov22[upper.tri(dist, diag = T)] <- Matern(dist[upper.tri(dist, diag = T)] * a_2, 
                                             smoothness = nu_2)
  cov22[lower.tri(dist, diag = F)] <- t(cov22)[lower.tri(dist, diag = F)]
  
  cov_all <- rbind(cbind(sigma_1 * cov11+ diag(nrow(dist), x = nugget_var_1), sigma_12 * cov12),
                   cbind(sigma_12 * t(cov12), sigma_2  * cov22+ diag(nrow(dist), x = nugget_var_2)))
  inv_all <- solve(cov_all)
  
  print(  format(- nrow(cov_all)/2 * log(2 * pi) - 1/2 * as.double(t(response) %*% inv_all %*% response) - 
                   1/2 * as.vector(determinant(cov_all)$modulus), scientific = F))
  -(- nrow(cov_all)/2 * log(2 * pi) - 1/2 * as.double(t(response) %*% inv_all %*% response) - 
      1/2 * as.vector(determinant(cov_all)$modulus))
}

d <- 2
approx_n <- 400
x_max <- 2
nu_val <- (d-2)/2
summary(as.double(dist))
x_vals  <- hankel_get_info_x(matrix(ncol = 1, nrow= approx_n, 1),
                             x_max,
                             nu_val,
                             approx_n)
k_vals <- hankel_get_info_k(matrix(ncol = 1, nrow= approx_n, 1),
                            x_max,
                            nu_val,
                            approx_n)

ml_product <- function(par, dist, response) {
  #theta: a, nu, sigma1, sigma2, sigma12
  a_1 <- exp(par[1])
  nu_1 <- exp(par[2])
  sigma_1 <- exp(par[3]);sigma_2 <- exp(par[4])
  sigma_12 <- par[5]*sqrt(sigma_1) * sqrt(sigma_2)
  nugget_var_1 <- exp(par[6]);nugget_var_2 <- exp(par[7])
  
  cov1 <- Matern(dist*a_1, smoothness = nu_1)
  
  cov_all <- rbind(cbind(sigma_1 * cov1+ diag(nrow(dist), x = nugget_var_1), sigma_12 * cov1),
                   cbind(sigma_12 * t(cov1), sigma_2  * cov1+ diag(nrow(dist), x = nugget_var_2)))
  inv_all <- solve(cov_all)
  
  -(- nrow(cov_all)/2 * log(2 * pi) - 1/2 * as.double(t(response) %*% inv_all %*% response) - 
      1/2 * as.vector(determinant(cov_all)$modulus))
}

ml_pars <- function(par, dist, response) {
  #theta: a, nu, sigma1, sigma2, sigma12
  a_1 <- exp(par[1])
  nu_1 <- exp(par[2]);nu_2 <- exp(par[3])
  nu_12 = (nu_1 + nu_2)/2
  sigma_1 <- exp(par[4]);sigma_2 <- exp(par[5])
  sigma_12 <- par[6]*sqrt(sigma_1) * sqrt(sigma_2)
  nugget_var_1 <- exp(par[7]);nugget_var_2 <- exp(par[8])
  
  cov1 <- Matern(dist*a_1, smoothness = nu_1)
  cov2 <- Matern(dist*a_1, smoothness = nu_2)
  cov12 <- Matern(dist*a_1, smoothness = nu_12)
  
  cov_all <- rbind(cbind(sigma_1 * cov1 + diag(nrow(dist), x = nugget_var_1), sigma_12 * cov12),
                   cbind(sigma_12 * t(cov12), sigma_2  * cov2 + diag(nrow(dist), x = nugget_var_2)))
  inv_all <- solve(cov_all)
  
  -(- nrow(cov_all)/2 * log(2 * pi) - 1/2 * as.double(t(response) %*% inv_all %*% response) - 
      1/2 * as.vector(determinant(cov_all)$modulus))
}


ml_full_bivariate <- function(par, dist, response) {
  print(exp(par))
  #theta: a, nu, sigma1, sigma2, sigma12
  a_1 <- exp(par[1]); a_2 <- exp(par[2]); a_12 <- exp(par[3])
  nu_1 <- exp(par[4]);nu_2 <- exp(par[5]); nu_12 =  exp(par[6])
  if (nu_12 < .5 *(nu_1 + nu_2)) {
    return(10^8)
  }
  sigma_1 <- exp(par[7]);sigma_2 <- exp(par[8])
  sigma_12 <- par[9]*sqrt(sigma_1) * sqrt(sigma_2)
  nugget_var_1 <- exp(par[10]);nugget_var_2 <- exp(par[11])
  
  cov1 <- Matern(dist*a_1, smoothness = nu_1)
  cov2 <- Matern(dist*a_2, smoothness = nu_2)
  cov12 <- Matern(dist*a_12, smoothness = nu_12)
  
  cov_all <- rbind(cbind(sigma_1 * cov1 + diag(nrow(dist), x = nugget_var_1), sigma_12 * cov12),
                   cbind(sigma_12 * t(cov12), sigma_2  * cov2 + diag(nrow(dist), x = nugget_var_2)))
  inv_all <- solve(cov_all)
  
  v_val <- -(- nrow(cov_all)/2 * log(2 * pi) - 1/2 * as.double(t(response) %*% inv_all %*% response) - 
      1/2 * as.vector(determinant(cov_all)$modulus))
  print(v_val)
  v_val
}


ml_ind <- function(par, dist, response) {
  a_1 <- exp(par[1]);a_2 <- exp(par[2])
  nu_1 <- exp(par[3]);nu_2 <- exp(par[4])
  sigma_1 <- exp(par[5]);sigma_2 <- exp(par[6])
  nugget_var_1 <- exp(par[7]);nugget_var_2 <- exp(par[8])
  
  cov1 <- Matern(dist*a_1, smoothness = nu_1)
  cov2 <- Matern(dist*a_2, smoothness = nu_2)
  
  cov_all <- rbind(cbind(sigma_1 * cov1 + diag(nrow(dist), x = nugget_var_1), 0 * cov1),
                   cbind(0 *cov1, sigma_2*cov2 + diag(nrow(dist), x = nugget_var_2)))
  inv_all <- solve(cov_all)
  -(- nrow(cov_all)/2 * log(2 * pi) - 1/2 * as.double(t(response) %*% inv_all %*% response) - 
      1/2 * as.vector(determinant(cov_all)$modulus))
}

nu_lower <- .02
nu_upper <- 3
a_lower <- 10^-7
a_upper <- 1
d <- 2

var(weather)
cor(weather)
a <- Sys.time()
test_optim <- optim(par = c(log(c(1/100, 1/100, .5, .5, 37974.705,7.388775)), -.4692,
                            log(c(67^2, .05))),
                    lower = c(log(c(10^(-7), 10^(-7), .02, .02, NA, NA)), -.96, NA, NA),
                    upper = c(log(c(1, 1, 2, 2, NA, NA)), .96, NA, NA),
                    method = 'L-BFGS-B',
                    x_vals = x_vals, k_vals = k_vals, 
                    d = d,
                    fn = ml_val, dist = dist, response = response)
b <- Sys.time()
b-a

test_optim$value # 1263.651
exp(test_optim$par)
test_optim$par
# 1.001108e-02 9.941713e-03 1.259391e+00 5.857067e-01 5.294730e+04 7.208793e+00 4.857589e-01
# 4.411825e+03 2.056416e-06
#11.58803 minutes


test_optim_product <- optim(par = c(log(c(1/100, .5, 37974.705,7.388775)), -.4692, log(c(67^2, .05))),
                            lower = c(log(c(10^(-7), .02, NA, NA)), -.96, NA, NA),
                            upper = c(log(c(1, 2, NA, NA)), .96, NA, NA),
                            method = 'L-BFGS-B',
                            fn = ml_product, dist = dist, response = response)
-test_optim_product$value
exp(test_optim_product$par)
test_optim_product$par
# -test_optim_product$value
# [1] -1268.251
# > exp(test_optim_product$par)
# [1] 6.713058e-03 5.439571e-01 3.860455e+04 9.059972e+00 6.690639e-01 2.463819e+03 5.343909e-06


test_optim_pars <- optim(par = c(log(c(1/100, .5, .5, 37974.705,7.388775)), -.4692,
                                 log(c(67^2, .05))),
                         lower = c(log(c(10^(-7), .02, .02, NA, NA)), -.96, NA, NA),
                         upper = c(log(c(1,2, 2, NA, NA)), .96, NA, NA),
                         method = 'L-BFGS-B',
                         fn = ml_pars, dist = dist, response = response)
-test_optim_pars$value
exp(test_optim_pars$par)
test_optim_pars$par
# -test_optim_pars$value
# [1] -1263.651
# > exp(test_optim_pars$par)
# [1] 1.176902e-02 1.539742e+00 6.118016e-01 5.146804e+04 6.668356e+00 6.027517e-01 4.786829e+03
# [8] 1.502564e-06

test_optim_ind <- optim(par = c(log(c(1/100,1/100, .5, .5, 37974.705,7.388775)),
                                log(c(67^2, .05))),
                        lower = c(log(c(10^(-7), 10^(-7),.02, .02, NA, NA, NA, NA))),
                        upper = c(log(c(1, 1, 2, 2, NA, NA, NA, NA))),
                        method = 'L-BFGS-B',
                        fn = ml_ind, dist = dist, response = response)
-test_optim_ind$value
exp(test_optim_ind$par)

# -test_optim_ind$value
# [1] -1274.429
# > exp(test_optim_ind$par)
# [1] 1.474127e-02 1.092680e-02 2.000000e+00 5.960592e-01 5.190552e+04 6.817124e+00 4.759082e+03
# [8] 6.553402e-07


test_optim_full <- optim(par = c(log(c(1/100, 1/100, 1/100, .5, .5, 1, 37974.705,7.388775)),
                                 -.4692, log(c(67^2, .05))),
                        lower = c(log(c(10^(-7), 10^(-7), 10^(-7),.02, .02,.02, NA, NA, NA, NA, NA))),
                        upper = c(log(c(1, 1, 1, 2, 2, 2, NA, NA, NA, NA, NA))),
                        method = 'L-BFGS-B',
                        fn = ml_full_bivariate, dist = dist, response = response)
test_optim_full$value
exp(test_optim_full$par)

save(test_optim_ind, test_optim_pars, test_optim_product, test_optim,
     test_optim_full, 
     file = 'results/data_analysis_optim.RData')



