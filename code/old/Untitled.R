
struve_version <- function(h,nu, d) {
  if(h == 0) {
    return(0)
  }
  sign(h)*(abs(h))^nu * pi * 2^(-nu-d/2) / (gamma(nu + d/2) * cos(pi * nu))*
    (besselI(abs(h), nu = nu) - RandomFieldsUtils::struveL(abs(h), -nu))
}

# using incomplete bessel/struve functions

alpha_grid <- seq(0, pi/2, length.out = 200)
alpha_delta <- alpha_grid[2] - alpha_grid[1]

A_nu = function(nu) {
  2^nu * gamma(nu + 1/2) *gamma(1/2)
}
f_plus <- function(nu, alpha, z, alpha_grid, alpha_delta) {
  alpha_chosen <- alpha_grid[alpha_grid < alpha]
  f_val <- exp(z * cos(alpha_chosen)) * sin(alpha_chosen)^(2 * nu)
  f_val <- f_val[f_val < Inf]
  z^nu/A_nu(nu) *sum(f_val * alpha_delta)
}
f_minus <- function(nu, alpha, z, alpha_grid, alpha_delta) {
  alpha_chosen <- alpha_grid[alpha_grid < alpha]
  f_val <- exp(-z * cos(alpha_chosen)) * sin(alpha_chosen)^(2 * nu)
  f_val <- f_val[f_val < Inf]
  z^nu/A_nu(nu) *sum(f_val * alpha_delta)
}
I_nu <- function(nu, alpha, z, alpha_grid, alpha_delta) {
  f_plus(nu, alpha, z, alpha_grid, alpha_delta) +
    f_minus(nu, alpha, z, alpha_grid, alpha_delta)
}
L_nu <- function(nu, alpha, z, alpha_grid, alpha_delta) {
  alpha_chosen <- alpha_grid[alpha_grid < alpha]
  f_val <- sinh(z * cos(alpha_chosen))* sin(alpha_chosen)^(2 * nu)
  f_val <- f_val[f_val < Inf]
  z^nu/A_nu(nu) * sum(f_val * alpha_delta)*2
}
H_nu <- function(nu, alpha, z, alpha_grid, alpha_delta) {
  # alpha_chosen <- alpha_grid[alpha_grid < alpha]
  # f_val <- sin(z * cos(alpha_chosen))* sin(alpha_chosen)^(2 * nu)
  # 2 * z^nu/A_nu(nu) * sum(f_val[f_val < Inf] * alpha_delta)
  m_seq <- 0:1000
  in_beta <- pbeta(sin(alpha)^2, shape1=  nu+1/2, shape2 = m_seq + 1) * 
    beta(nu + 1/2, m_seq + 1)
  test <- (-1)^m_seq  *(in_beta) * (z)^(2*m_seq + 1)/ factorial(2*m_seq + 1)
  test <- ifelse(is.nan(test), 0, 
                 ifelse(test == Inf, 0, test))
  z^nu/A_nu(nu) * sum(test)
}
whole_int_h <- function(nu, alpha, z, alpha_grid, alpha_delta,d) {
  r_delta <- .1
  r_seq <- seq(0, 1.3*z, length.out = 300)
  int_vals <- sapply(1:length(r_seq), function(x) {
    H_nu(nu = d/2-1, alpha = alpha, z = r_seq[x] * z, alpha_grid, alpha_delta)*
      (1 + r_seq[x]^2)^(-nu - 1/2) * r_seq[x]^(d/2)
  }) 
  sum(int_vals * r_delta)
}

ImL <- function(nu, alpha, z, alpha_grid, alpha_delta){
  delta <- .01
  r_seq <- seq(0, 200, by = delta)
  # alpha_chosen <- alpha_grid[alpha_grid < alpha]
  # f_val <- sin(alpha_chosen)^(2 * nu)
  # #f_val <- f_val[f_val < Inf]
  # 2 * (1/2 * z)^nu/sqrt(pi)/gamma(nu + 1/2) * sum(f_val * alpha_delta)
  f_val <- sin(r_seq * z) *(1 + r_seq^2)^(-nu - 1/2)
  2 * (z/2)^(-nu)/sqrt(pi)/gamma(nu + 1/2) * sum(f_val * delta)
}

test <- sapply(alpha_grid, function(x) {
  I_nu(.25, x, 4, alpha_grid, alpha_delta)
})
plot(alpha_grid, test, type = 'l')
points(pi/2, y=besselI(4, .25), col = 2)

test <- sapply(alpha_grid, function(x) {
  L_nu(-.25, x, 3, alpha_grid, alpha_delta)
})
plot(alpha_grid, test, type = 'l')
points(pi/2, y=RandomFieldsUtils::struveL(x = 3,-.25), col = 2)

x_seq <- seq(-5, 5, by = .01)
x_grid <- as.matrix(expand.grid(x_seq, x_seq))
d=2
nu = .25
test <- sapply(x_seq, function(x) {
  # a <- pi * 2^(-nu-d/2) / (gamma(nu + d/2) * cos(pi * nu)) *
  #   sign(x)* (abs(x)^nu)*
 ((abs(x))^nu)*sign(x)*pi * 2^(-nu-d/2) / (gamma(nu + d/2) * cos(pi * nu))*
  ImL(nu, pi/2, abs(x), alpha_grid, alpha_delta)
    # (I_nu(nu, pi/2, abs(x), alpha_grid, alpha_delta)-
    #    L_nu(-nu, pi/2, abs(x), alpha_grid, alpha_delta))
  #L_nu(nu, pi/2, abs(x), alpha_grid, alpha_delta)
  #H_nu(nu, pi/2, abs(x), alpha_grid, alpha_delta)
})
test2 <- sapply(x_seq, function(x) {
  #RandomFieldsUtils::struveH(x, nu = nu)
  struve_version(x, nu, d=d)
})
test3 <- sapply(x_seq, function(x) {
  sign(x)*whole_int_h(nu, pi/2, abs(x), alpha_grid, alpha_delta, d = d)
})
plot(x_seq, test, type = 'l')
lines(x_seq, y=test2, col = 2)
lines(x_seq, y=test3, col = 3)

library(fields)
a <- c(1,0)
a_perp <- c(0,1)
d <- length(a)
nu = .25

res <- sapply(1:nrow(x_grid), function(y) {
  a_norm <- sqrt(sum(a_perp^2))
  h_norm <- sqrt(sum(x_grid[y,]^2))
  theta <- acos(sum(a_perp * x_grid[y,])/a_norm/h_norm)
  if (is.nan(theta)) {theta = pi/2}
  if (theta > pi/2) { theta = pi - theta}
  test <- (h_norm^nu) * (I_nu(nu, theta, h_norm, alpha_grid, alpha_delta) -
      L_nu(-nu, theta, h_norm, alpha_grid, alpha_delta))
  # test <- (I_nu(nu, theta, h_norm, alpha_grid, alpha_delta) - 
  #                          L_nu(-nu, theta, h_norm, alpha_grid, alpha_delta))
  if (x_grid[y,1] < 0) {
    -test
  } else {
    test
  }
  #theta
  # (I_nu(nu, theta, h_norm, alpha_grid, alpha_delta) - 
  #   L_nu(nu, theta, h_norm, alpha_grid, alpha_delta))
  #theta
})
image.plot(matrix(res, nrow = length(x_seq)), main = 'nu = 0.75')


fm2 <- f_minus(nu = -nu, alpha = gamma, z = h_norm, alpha_grid, alpha_delta)
RandomFieldsUtils::struveL(4, .75)
L_nu(nu = .75, alpha = pi/2, z = 4, alpha_grid, alpha_delta)
besselI(4, .75)
I_nu(nu = .75, alpha = pi/2, z = 4, alpha_grid, alpha_delta)

RandomFieldsUtils::struveL(8, .75)
L_nu(nu = .75, alpha = pi/2, z = 8, alpha_grid, alpha_delta)
besselI(8, .75)
I_nu(nu = .75, alpha = pi/2, z = 8, alpha_grid, alpha_delta)

RandomFieldsUtils::struveL(8, .75)
L_nu(nu = .75, alpha = pi/4, z = 8, alpha_grid, alpha_delta)
besselI(8, .75)
I_nu(nu = .75, alpha = pi/4, z = 8, alpha_grid, alpha_delta)


struve_intern(x = 4, nu = -5.25, factor_Sign = 1, expscaled = F)
struve_intern <- function(x,nu, factor_Sign, expscaled) { 
  if ((x == 0.0) && (nu>-1.0)) return(0.0);
  if (x <= 0.0) return(RF_NA)
  dummy = 0.0
  logx = 2.0 * log(0.5 * x)
  x1 = 1.5
  x2 = nu + 1.5
  value = 1.0
  fsign = factor_Sign
  epsilon=1e-20
  exp_dummy = exp(dummy)
  
  while(exp_dummy > abs(value) * epsilon) {
    dummy = dummy + logx - log(x1) - log(abs(x2))
    exp_dummy = exp(dummy);
    value = value +  (1 - 2 * (x2 < 0))  * fsign * exp_dummy;
    print(paste(value, x1,x2, fsign));
    x1 = x1 +  1.0;
    x2 = x2 + 1.0;
    fsign = factor_Sign * fsign; 
  }
  
  x1 = 1.5;
  x2 = nu + 1.5;
  if (x2 > 0.0) { 
    dummy = (nu + 1.0) * 0.5 * logx - gamma(x1) - gamma(x2);
    if (expscaled) dummy = dummy - x;
    value = value * exp(dummy);
  } else {
    value = value * (0.5 * x)^(nu + 1.0) / (gamma(x1) * gamma(x2));
    if (expscaled) value =  value * exp(-x);
  }
  value
}

