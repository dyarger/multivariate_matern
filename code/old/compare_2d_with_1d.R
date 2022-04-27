load(file = paste0('results/asymmetric_2d_300.RData'))


head(grid_vals)
library(tidyverse)
ggplot(data = grid_vals) + 
  geom_tile(aes(x = Var1, y = Var2, fill = res_vals)) + 
  scale_fill_gradient2() + 
  labs(x = expression(t[1]), y = expression(t[2]),
       fill = 'Cross\nCovariance') + 
  theme(legend.position = 'bottom')
grid_vals1d <- grid_vals %>%
  filter(Var2 == 0)

source('code/multi_matern_source.R')

plot_function <- function(s,t,nu) {
  if(s-t == 0) {
    return(0)
  }
  sign(t-s) * (abs(t-s))^nu*
    (besselI(abs(t-s), nu = nu) - struve(abs(t-s), -nu))
}

plot_seq <- grid_vals1d$Var1
nu25 <- sapply(1:length(plot_seq), function(x) plot_function(0, plot_seq[x], nu = .499))
plot(grid_vals1d$Var1, grid_vals1d$res_vals, type = 'l')
lines(plot_seq, -nu25 * 730, col = 2)
plot(grid_vals1d$res_vals, nu25)
abline(lm(nu25~grid_vals1d$res_vals))
summary(lm(nu25~grid_vals1d$res_vals))

I_incomplete <- function(nu, h, c) {
  f <- function(theta) {
    cosh(h * cos(theta)) * sin(theta)^(2 * nu)
  }
  int_val <- integrate(f = f, lower = 0, upper = c)
  int_val[['value']] * 
    2 * h^nu/ (2^nu * gamma(nu + 1/2) * gamma(1/2))
}
L_incomplete <- function(nu_prime, h, c, tol = .00001) {
  f <- function(theta) {
    sinh(h * cos(theta)) * sin(theta)^(2 * nu_prime)
  }
  # test <- seq(0, c, length.out = 100)
  # plot(test, f(test))
  # f(0)
  # sinh(h * cos(0)) * sin(0)^(2 * nu_prime)
  int_val <- integrate(f = f, lower = 0 + tol, upper = c, subdivisions = 10000)
  int_val[['value']] * 
    2 * h^nu_prime/ (2^nu_prime * gamma(nu_prime + 1/2) * gamma(1/2))
}
besselI(1, nu = nu)
I_incomplete(nu = nu, h = 1, c = pi/2)

struve(1, nu = -nu)
L_incomplete(nu = nu, h = 1, c = pi/2, tol = 10^-8)

nu <- 0.499
grid_vals_mat <- as.matrix(grid_vals)
test_run <- rep(NA, nrow(grid_vals_mat))
for (i in 1:nrow(grid_vals)) {
  psi <- atan2(grid_vals_mat[i,2], grid_vals_mat[i,1])
  psi_abs <- ifelse(abs(psi) > pi/2, pi - abs(psi), abs(psi))
  h_val <- sqrt(grid_vals_mat[i,2]^2 + grid_vals_mat[i,1]^2)
  I_val <- I_incomplete(nu, abs(h_val), psi_abs)
  L_val <- L_incomplete(-nu, abs(h_val), psi_abs)
  #test_run[i] <- sign(psi) * h_val^nu * (I_val - L_val)
  #test_run[i] <- sign(psi) * h_val^nu * (I_val - L_val)
  test_run[i] <- sign(psi) * h_val^nu * 
    (I_val - L_val)
  test_run[i] <- psi_abs
  if (i %% 1000 == 0) {
    print(i)
  }
}

ggplot(data = cbind(grid_vals, abs(test_run))) + 
  geom_tile(aes(x = Var1, y = Var2, fill = test_run)) + 
  scale_fill_gradient2() + 
  labs(x = expression(t[1]), y = expression(t[2]),
       fill = 'Cross\nCovariance') + 
  theme(legend.position = 'bottom')
grid_vals1d_test <- cbind(grid_vals, vals = abs(test_run)) %>%
  filter(abs(Var2 -  0.1) < .001)
plot(grid_vals1d_test$Var1, grid_vals1d_test$vals, type = 'l')
grid_vals1d_test <- cbind(grid_vals, vals = abs(test_run)) %>%
  filter(abs(Var1 -  0.1) < .001)
plot(grid_vals1d_test$Var2, grid_vals1d_test$vals, type = 'l')



# alternate plan H function

H_incomplete <- function(nu, h, c) {
  f <- function(theta) {
    sin(h * cos(theta)) * sin(theta)^(2 * nu)
  }
  int_val <- integrate(f = f, lower = 0, upper = c)
  int_val[['value']] * 
    2 * h^nu/ (2^nu * gamma(nu + 1/2) * gamma(1/2))
}
struveH <- function(z, nu_eval) {
  if (nu_eval == -1/2) {
    return(sqrt(2/(pi * z)) * sinh(z))
  } else if (nu_eval == -3/2) {
    return(sqrt(2/pi) * (z * cosh(z) - sinh(z))/(z^(3/2)))
  }
  k <- 0:200
  (z/2)^(nu_eval + 1) *
    sum((z/2)^(2*k) *(-1)^k /(gamma(k + 3/2) * gamma( k  + nu_eval + 3/2)))
}

eval_function <- function(nu, h, c, d = 2, approx_seq = seq(0, 40, length.out = 500)) {
  full_function <- function(x, nu, h, c, d = 2) {
    H_incomplete((d-2)/2, h * x, c) * 
      (1 + x^2)^(-nu - d/2) * x^(d/2)
  }
  rr_val <- rr <- approx_seq
  for (i in 1:length(rr_val)){
    rr_val[i] <- full_function(rr[i], nu, h, d, c = c)
  }
  sum(rr_val/(rr[2] - rr[1]), na.rm = T)* 2 * h^(-(d-2)/2)/ 
    (2^nu * gamma(nu + 1/2) * gamma(1/2))
  # int_val <- integrate(f = full_function, lower = 0 , upper = 20, subdivisions = 10000,
  #                      c = c, h = h, nu = nu, d = d)
  # int_val[['value']] * 
  #   2 * h^((-d+2)/2)/ (2^nu * gamma(nu + 1/2) * gamma(1/2))
}
besselI(1, nu = nu)
I_incomplete(nu = nu, h = 1, c = pi/2)
struveH(1, nu_eval = .3)
H_incomplete(.3, 1, pi/2)
struve(1, nu = -nu)
L_incomplete(nu = nu, h = 1, c = pi/2, tol = 10^-8)

nu <- 0.499
grid_vals_mat <- as.matrix(grid_vals)
test_run <- rep(NA, nrow(grid_vals_mat))
for (i in 1:nrow(grid_vals)) {
  psi <- atan2(grid_vals_mat[i,2], grid_vals_mat[i,1])
  psi_abs <- ifelse(abs(psi) > pi/2, pi - abs(psi), abs(psi))
  h_val <- sqrt(grid_vals_mat[i,2]^2 + grid_vals_mat[i,1]^2)
  # I_val <- I_incomplete(nu, abs(h_val), psi_abs)
  # L_val <- L_incomplete(-nu, abs(h_val), psi_abs)
  #test_run[i] <- sign(psi) * h_val^nu * (I_val - L_val)
  #test_run[i] <- sign(psi) * h_val^nu * (I_val - L_val)
  # test_run[i] <- sign(psi) * h_val^nu * 
  #   (I_val - L_val)
  test_run[i] <- sign(psi) * eval_function(nu, h_val, c = psi_abs, d = 2)
  if (i %% 1000 == 0) {
    print(i)
  }
}

ggplot(data = cbind(grid_vals, test1 = test_run)) + 
  geom_tile(aes(x = -Var2, y = Var1, fill = test1)) + 
  scale_fill_gradient2() + 
  labs(x = expression(t[1]), y = expression(t[2]),
       fill = 'Cross\nCovariance') + 
  theme(legend.position = 'bottom')
grid_vals1d_test_one <- cbind(grid_vals, vals = (test_run)) %>%
  filter(abs(Var2 -  0) < .001)
plot(grid_vals1d_test_one$Var1, grid_vals1d_test_one$vals, type = 'l')
grid_vals1d_test <- cbind(grid_vals, vals = (test_run)) %>%
  filter(abs(Var1 -  0) < .001)
plot(grid_vals1d_test$Var2, grid_vals1d_test$vals, type = 'l')
lines(grid_vals1d_test_one$Var1, -grid_vals1d_test_one$res_vals*15, type = 'l')
plot(grid_vals1d_test_one$Var1, grid_vals1d_test_one$res_vals, type = 'l')


eval_function_simple <- function(nu, h, d, tol = .00001, approx_seq = 
                                   seq(0, 20, length.out = 200)) {
  full_function <- function(x, nu, h, d) {
    struveH(h*x, (d-2)/2) * 
      (1 + x^2)^(-nu - d/2) * x^(d/2)
  }
  rr_val <- rr <- approx_seq
  for (i in 1:length(rr_val)){
    rr_val[i] <- full_function(rr[i], nu, h, d)
  }
  sum(rr_val/(rr[2] - rr[1]), na.rm = T)* 2 * h^(-(d-2)/2)/ 
    (2^nu * gamma(nu + 1/2) * gamma(1/2))
}

tt <- seq(-2, 2, by = .02)
tt_val <- tt
for (i in 1:length(tt)) {
  tt_val[i] <- sign(tt[i]) * eval_function_simple(nu, abs(tt[i]), d = 2)
}
plot(tt, tt_val, type = 'l')


plot(grid_vals1d_test_one$Var1, grid_vals1d_test_one$res_vals, type = 'l')

test_vals_one <- rep(0, nrow(grid_vals1d_test_one))
for (i in 1:nrow(grid_vals1d_test_one)) {
  test_vals_one[i] <- eval_function_simple(nu, grid_vals1d_test_one$Var1[i], d = 2)
}
plot(grid_vals1d_test_one$Var1, grid_vals1d_test_one$res_vals, type = 'l')
lines(grid_vals1d_test_one$Var1, -test_vals_one/53, type = 'l', col = 2)

tt <- seq(-10, 10, by = .1)
tt_val <- tt
for (i in 1:length(tt)) {
  tt_val[i] <- eval_function_simple(nu, tt[i], d = 2)
}
plot(tt, tt_val, type = 'l')



