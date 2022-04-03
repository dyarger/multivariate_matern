
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

ml_val <- function(par, dist, response) {
  #theta: a_1, a_2, nu_1, nu_2, sigma1, sigma2, sigma12, nugget1, nugget2
  print(format(exp(par), scientific = F))
  print(par[7])
  a_1 <- exp(par[1]);a_2 <- exp(par[2])
  nu_1 <- exp(par[3]);nu_2 <- exp(par[4])
  sigma_1 <- exp(par[5]);sigma_2 <- exp(par[6])
  sigma_12 <- par[7]*sqrt(sigma_1) * sqrt(sigma_2)
  nugget_var_1 <- exp(par[8]);nugget_var_2 <- exp(par[9])
  
  cov12 <- dist
  cov12[upper.tri(dist, diag = T)] <- sapply(dist[upper.tri(dist, diag = T)], spatial_test,
                                             d = 2, a_1 = a_1, a_2 =a_2,  nu_1 = nu_1, nu_2 = nu_2,
                                             limits = c(.00000001, 2), subdivisions = 100)
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

spatial_test <- function(h,d, a_1, a_2, nu_1, nu_2, limits = c(.00000001, .5), subdivisions = 300) {
  spec_dens_single <- function(x, h, d, a_1, a_2, nu_1, nu_2) {
    y <- complex(imaginary = x, length.out = length(x))
    besselJ(h*x, d/2 - 1) *
      Re((a_1 + y)^(-nu_1-d/2)*(a_2 - y)^(-nu_2-d/2)) *
      x^(d/2)
  }
  if (h == 0 & nu_1 == nu_2 & a_1 == a_2) {
    return(1)
  }
    
  # test <- sapply(seq(limits[1], limits[2], by = .01),
  #                FUN = spec_dens_single,
  #                h = h, d = d, a_1 = a_1, a_2 = a_2, nu_1 = nu_1, 
  #                nu_2 = nu_2)
  test_integrate <- integrate(lower = limits[1], upper = limits[2], 
                              f = spec_dens_single, stop.on.error = F,
                              h = h, d = d, a_1 = a_1, a_2 = a_2, nu_1 = nu_1, 
                              nu_2 = nu_2, subdivisions = subdivisions)
  (2 * pi)^(d/2-1) *h^(-d/2 + 1) * test_integrate[['value']]* sqrt(nu_1) * sqrt(nu_2)*2* 
    a_1^nu_1 * a_2^nu_2 
}

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

var(weather)
cor(weather)
test_optim <- optim(par = c(log(c(1/100, 1/100, .5, .5, 37974.705,7.388775)), -.4692,
                            log(c(67^2, .05))),
                   lower = c(log(c(10^(-7), 10^(-7), .02, .02, NA, NA)), -.96, NA, NA),
                   upper = c(log(c(1, 1, 2, 2, NA, NA)), .96, NA, NA),
                   method = 'L-BFGS-B',
                    fn = ml_val, dist = dist, response = response)
-test_optim$value
exp(test_optim$par)


test_optim_product <- optim(par = c(log(c(1/100, .5, 37974.705,7.388775)), -.4692, log(c(67^2, .05))),
                    lower = c(log(c(10^(-7), .02, NA, NA)), -.96, NA, NA),
                    upper = c(log(c(1, 2, NA, NA)), .96, NA, NA),
                    method = 'L-BFGS-B',
                    fn = ml_product, dist = dist, response = response)
-test_optim_product$value
exp(test_optim_product$par)

test_optim_pars <- optim(par = c(log(c(1/100,1/100, .5, .5, 37974.705,7.388775)), -.4692,
                                 log(c(67^2, .05))),
                            lower = c(log(c(10^(-7),10^(-7), .02, .02, NA, NA)), -.96, NA, NA),
                            upper = c(log(c(1, 2, 2, NA, NA)), .96, NA, NA),
                            method = 'L-BFGS-B',
                            fn = ml_pars, dist = dist, response = response)
-test_optim_pars$value
exp(test_optim_pars$par)

test_optim_ind <- optim(par = c(log(c(1/100,1/100, .5, .5, 37974.705,7.388775)),
                                log(c(67^2, .05))),
                         lower = c(log(c(10^(-7), 10^(-7),.02, .02, NA, NA, NA, NA))),
                         upper = c(log(c(1, 1, 2, 2, NA, NA, NA, NA))),
                         method = 'L-BFGS-B',
                         fn = ml_ind, dist = dist, response = response)
-test_optim_ind$value
exp(test_optim_ind$par)

save(test_optim_ind, test_optim_pars, test_optim_product, test_optim,
     file = 'data_analysis_optim.RData')

# -test_optim_product$value
# [1] -1268.251
# > exp(test_optim_product$par)
# [1] 6.713058e-03 5.439571e-01 3.860455e+04 9.059972e+00 6.690639e-01 2.463819e+03 5.343909e-06

# -test_optim_pars$value
# [1] -1263.651
# > exp(test_optim_pars$par)
# [1] 1.176902e-02 1.539742e+00 6.118016e-01 5.146804e+04 6.668356e+00 6.027517e-01 4.786829e+03
# [8] 1.502564e-06

# -test_optim_ind$value
# [1] -1274.429
# > exp(test_optim_ind$par)
# [1] 1.474127e-02 1.092680e-02 2.000000e+00 5.960592e-01 5.190552e+04 6.817124e+00 4.759082e+03
# [8] 6.553402e-07

lag_seq <- seq(0, 5, by = .05)
nu <- .5
res <- sapply(lag_seq, spatial_test, d = 2, a_1 = 1.5, a_2 = 1.5, nu_1 = nu, nu_2 = nu,
              limits = c(.0000000001, 1000))
plot(res, fields::Matern(lag_seq, smoothness = nu), type  = 'l')
plot(lag_seq, res, type  = 'l')
lines(lag_seq, fields::Matern(lag_seq, smoothness = nu), type  = 'l', col = 2)


approx_n <- 800
x_max <- 300
nu_val <- 0
x_vals  <- hankel_get_info_x(matrix(ncol = 1, nrow= approx_n, 1),
                             x_max,
                             nu_val,
                             approx_n)
k_vals <- hankel_get_info_k(matrix(ncol = 1, nrow= approx_n, 1),
                            x_max,
                            nu_val,
                            approx_n)
plot(x_vals, k_vals, type= 'l')
a_1 <- a_2 <- 2
nu_1 = .2; nu_2 = .2
d = 2

spec_dens_cpp <- function(x, d, a_1, a_2, nu_1, nu_2) {
  y <- complex(imaginary = x, length.out = length(x))
  Re((a_1 + y)^(-nu_1-d/2)*(a_2 - y)^(-nu_2-d/2)) *
    x^(d/2 - 1)
}
test_input <- spec_dens_cpp(x_vals, d =d , 
                               a_1 = a_1, a_2 = a_2, nu_1 = nu_1, nu_2 = nu_2)
test_integrate <- k_vals
for (iter_k in 1:length(k_vals)) {
  test_integrate[iter_k] <- integrate(lower = 10^-12, upper = 500, 
                              f = spec_dens_single_cpp, stop.on.error = F,
                              h = k_vals[iter_k], d = d, a_1 = a_1, a_2 = a_2, nu_1 = nu_1, 
                              nu_2 = nu_2, subdivisions = 500)$value
}
test <- spec_dens_single_cpp(matrix(ncol = 1, nrow= approx_n, test_input),
                             matrix(ncol = 1, nrow= approx_n, 0),
                             x_max,
                             nu_val,
                             approx_n)

test_matern <- (2 * pi)^(d/2-1) *k_vals^(-d/2 + 1) * test* sqrt(nu_1) * sqrt(nu_2)*2* 
  a_1^nu_1 * a_2^nu_2 
test_integralmatern <- (2 * pi)^(d/2-1) *k_vals^(-d/2 + 1) * test_integrate* sqrt(nu_1) * sqrt(nu_2)*2* 
  a_1^nu_1 * a_2^nu_2 

plot(k_vals, test_matern, type = 'l')
lines(k_vals, test_integralmatern, col = 2)
lines(k_vals, Matern(k_vals, smoothness = .2), col = 2)
all.equal(test_integralmatern, Matern(k_vals, smoothness = .6))
all.equal(test_matern, Matern(k_vals, smoothness = .6))
all.equal(test_matern, test_integralmatern)

plot(test_integrate, test, type = 'l')
abline(a = 0, b = 1)
all.equal(test_integrate, test)
abs(test_integrate - test)

# timing


approx_n <- 400
x_max <- 10
nu_val <- 0
x_vals  <- hankel_get_info_x(matrix(ncol = 1, nrow= approx_n, 1),
                             x_max,
                             nu_val,
                             approx_n)
k_vals <- hankel_get_info_k(matrix(ncol = 1, nrow= approx_n, 1),
                            x_max,
                            nu_val,
                            approx_n)
a_1 <- a_2 <- 1
nu_1 = .6; nu_2 = .6
d = 2
rbenchmark::benchmark(one = {

spec_dens_cpp <- function(x, d, a_1, a_2, nu_1, nu_2) {
  y <- complex(imaginary = x, length.out = length(x))
  Re((a_1 + y)^(-nu_1-d/2)*(a_2 - y)^(-nu_2-d/2)) *
    x^(d/2 - 1)
}
test_input <- spec_dens_cpp(x_vals, d =d , 
                            a_1 = a_1, a_2 = a_2, nu_1 = nu_1, nu_2 = nu_2)
test <- spec_dens_single_cpp(matrix(ncol = 1, nrow= approx_n, test_input),
                             matrix(ncol = 1, nrow= approx_n, 0),
                             x_max,
                             nu_val,
                             approx_n)},
                      two = {test_integrate <- k_vals
                      for (iter_k in 1:length(k_vals)) {
                        test_integrate[iter_k] <- integrate(lower = 10^-12, upper = 500, 
                                                            f = spec_dens_single, stop.on.error = F,
                                                            h = k_vals[iter_k], d = d, a_1 = a_1, a_2 = a_2, nu_1 = nu_1, 
                                                            nu_2 = nu_2, subdivisions = 500)$value
                      }},
                      replications = 5)



test <- spec_dens_single_cpp(matrix(ncol = 1, nrow= approx_n, test_input),
                     matrix(ncol = 1, nrow= approx_n, 0),
                     x_max,
                     nu_val,
                     approx_n)
plot(k_vals, test, type = 'l')



test
spec_dens_single_cpp(test,
                     matrix(ncol = 1, nrow= 2, c(50,50)),
                     1,
                     0,
                     2)
gsl_dht_x_sample
# use this to figure out where we need to evaluate the function

gsl_dht_k_sample
# use this to figure out where the result is defined



# double *finvals[size];
# for (int i =0; i < size; i++) {
#   finvals[i] = &fvals[i];
#   foutvals[i] = 0.0;
# }
spec_dens_single_cpp(c(4,4), 4, 0, 2)
int arr[5] = { 1, 2, 3, 4, 5 };
int *ptr = arr;

# library(Rcpp)
# 
# cppFunction('std::vector<double> spec_dens_single_cpp(std::vector<double> x, double h, double d, double a1, double a2, double nu1, double nu2) {
# std::vector<double> z = x;
# for(int i = 0; i < x.size(); i++) {
#   std::complex<double> y (0, x[i]);
#   z[i] = R::bessel_j(h *x[i], d/2 - 1) * std::real(pow(a1 + y, -nu1 - d/2) * pow(a2 - y, -nu2 - d/2)) * pow(x[i],d/2);
# }
# return z;
# }         ')
# #return R::bessel_j(x[0], al);
# spec_dens_single_cpp(c(1,1, 6),2, 2, 2, 2, 2, 2)
# spec_dens_single(c(1,1, 6),2, 2, 2, 2, 2, 2)
# rbenchmark::benchmark(one = {spec_dens_single(c(1,1, 6),2, 2, 2, 2, 2, 2)},
#                       two = {spec_dens_single_cpp(c(1,1, 6),2, 2, 2, 2, 2, 2)},
#                       replications = 20000)
# library(devtools)
# 
# devtools::install_bitbucket("davidbolin/ngme",ref="default")
