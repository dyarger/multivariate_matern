
z <- 1

omega <- function(a, b, z, n) {
  cons <- gamma(b)/gamma(a)/gamma(b-a)
  vals <- sapply(0:n, function(k) {
    (b_half(a = b-a, b = a + k) -
      b_half(b = b-a, a = a + k)) * z^k/factorial(k)})
  #vals2 <- sapply(0:n, function(k) {b_half(b = b-a, a = a + k)})
  cons * sum(vals)
}
# Re(b) > Re(a) > 0
b_half <- function(a, b) {
  pbeta(1/2,a,b) * beta(a,b)
}

Omega_a_c_0 <- function(a,c) { 
  if (a > c | a < 0) {
    stop('need c > a > 0')
  }
  gamma(c)/gamma(a)/gamma(c-a) * (b_half(c-a, a) - b_half(a, c-a))
}
Omega_a_c_0(.2, .3)
Omega_a_c_0(.5, .3)
Omega_a_c_0(a,c)
Omega_a_c_0_full(a,c)
Omega_a_c_0_full <- function(a,c) {
  vv <- find_valid_values(a,c)
  vv_value <- Omega_a_c_0(vv[1], vv[2])
  #print(vv)
  #print(vv_value)
  vv_difference <- vv - c(a,c)
  #print(vv_difference)
  
  if (vv_difference[1] > vv_difference[2]) {
    stop('uh oh')
  }
  if (abs(vv_difference[1] == vv_difference[2]) < 10^-5) {
    
  } else {
    #print(paste('subtracting both', vv_difference[1], 'times'))
    for (i in 1:(vv_difference[1])) {
      #print(i)
      vv_difference[1] <- vv_difference[1] - 1
      vv_difference[2] <- vv_difference[2] - 1
      #print(vv_difference)
      vv_value <- Omega_a_c_0_subtract_one_each(vv[1], vv[2], vv_value)
      vv[1] <- vv[1] - 1
      vv[2] <- vv[2] - 1
      #print(vv)
      #print(vv_value)
    }
  }
  
  if (vv_difference[2] > 0) {
    #print(paste('subtracting c', vv_difference[2], 'times'))
    for (i in 1:vv_difference[2]) {
      #print(i)
      vv_difference[2] <- vv_difference[2] - 1
      #print(vv_difference)
      vv_value <- Omega_a_c_0_subtract_one_c(vv[1], vv[2], vv_value)
      vv[2] <- vv[2] - 1
      # print(vv)
      # print(vv_value)
    }
  }
  vv_value
}
Omega_a_c_0_full(-.3, -2.1)

test_seq <- seq()

Omega_a_c_0_subtract_one_each <- function(a, c, current_value) { # 4.60
  current_value-2^(2-c)* gamma(c)/gamma(a)/gamma(c-a)/(c-1)
}
Omega_a_c_0_subtract_one_c <- function(a, c, current_value) { # 4.58
  current_value+2^(2-c) *gamma(c)/gamma(a)/gamma(c-a)/(c-1)
}

find_valid_values <- function(a,c) {
  a_good_value <- a + ceiling(abs(a)) 
  c_good_value <- c + ceiling(abs(c))
  while(a_good_value - c_good_value >= -10^-5) {
    c_good_value <- c_good_value + 1
  }
  c(a_good_value, c_good_value)
}
find_valid_values(-4.2, -2.2)
find_valid_values(-2.1, -2.2)
find_valid_values(-5.2, -2.2)
find_valid_values(-2.2, -5.2)
find_valid_values(-2.1, -5.2)
find_valid_values(-2.1, -5.1)
find_valid_values(-2.2, -5.1)


Omega <- function(a,b,z, n) { # 4.108
  sum(sapply(0:n, function(x) 
    gamma(a + x)/gamma(a) *gamma(b)/gamma(b+x)*
      Omega_a_c_0_full(a + x, b + x) * z^x/factorial(x)))
}

nu1 <- .501
nu2 <- nu1#+ 1
a <- 1/2 - nu1
b = 1 - nu1 - nu2
test_seq <- val_seq <- seq(0, 5, by = .01)
for (i in 1:length(test_seq)) {
  val_seq[i] <- exp(-test_seq[i]) *
    (2*test_seq[i])^(1/2-nu1/2 - nu2/2)
    Omega(a = a, b = b, z =  2 * test_seq[i], n  =40)
}

i_seq <- fAsianOptions::whittakerM(x = 2 * test_seq, kappa = nu1/2 - nu2/2, mu = nu1/2 + nu2/2)

full_cov <- (2 * test_seq)^(nu1/2 + nu2/2-1/2)/gamma(1/2 + nu1) *
  (val_seq/gamma(1 - nu1 - nu2)/gamma(1/2 + nu2) - 
     i_seq/gamma(1 + nu1 + nu2)/gamma(1/2 - nu1))

plot(test_seq, full_cov, type = 'l')

plot(i_seq, val_seq, type = 'l')
plot(test_seq, val_seq, type = 'l')

plot(i_seq, 2^(nu1 + nu2 + 1/2) * sqrt(test_seq) *besselI(nu = nu1/2 + nu2/2, x = test_seq)*
       gamma(1 + nu1/2 + nu2/2),
     type = 'l')
plot(2^(nu1 + nu2 + 1/2) * sqrt(test_seq) *
       sapply(test_seq, function(x) struve(nu_eval = -nu1/2- nu2/2, z = x))*
       gamma(1 + nu1/2 + nu2/2),
     2^(nu1 + nu2 + 1/2) * sqrt(test_seq) *besselI(nu = nu1/2 + nu2/2, x = test_seq)*
       gamma(1 + nu1/2 + nu2/2), type = 'l')
plot(test_seq, 2^(nu1 + nu2 + 1/2) * sqrt(test_seq) *
       sapply(test_seq, function(x) struve(nu_eval = -nu1/2- nu2/2, z = x))*
       gamma(1 + nu1/2 + nu2/2)-
       2^(nu1 + nu2 + 1/2) * sqrt(test_seq) *besselI(nu = nu1/2 + nu2/2, x = test_seq)*
       gamma(1 + nu1/2 + nu2/2), type = 'l')
abline(b = 1, a = 0)
source('code/multi_matern_source.R')
plot(val_seq, 2^(nu1 + nu2 + 1/2) * sqrt(test_seq) *
       sapply(test_seq, function(x) struve(nu_eval = -nu1/2- nu2/2, z = x))*
       gamma(1 + nu1/2 + nu2/2),
     type = 'l')
abline(b = 1, a = 0)

Omega_a_c_0_full(2, 4)

a <- -.9
b <- -.8
Omega_a_c_0(a,b)
Omega_a_c_0_full(a, b)
Omega_a_c_0_subtract_one_c(a = a, c = b+1, Omega_a_c_0(a,b + 1))
Omega_a_c_0_subtract_one_each(a = a+1, c = b+1, Omega_a_c_0(a + 1,b + 1))
vv_value <- Omega_a_c_0_subtract_one_c(vv[1], vv[2], vv_value)

Omega_a_c_0(1,b) - 2^(2-b) + 1



test <- Omega(a = a, b = b,z =  1, n  =10)
sapply(0:n, function(x) 
  Omega_a_c_0(a + x, b + x))
Omega(1, 1.2, .5, 90)

A_fun <- function(kappa, mu, z, n) {
  omega_results <- sapply()
  
  exp(-1/2 * z) *  z^(1/2 + mu) * Omega(1/2 - kappa + mu, 1 + 2 * mu, z= z, n = n)
}

test <- seq(0, 5, by = .1)
test_vals <- sapply(test, omega, a = 3, b = 2.2, n = 100)

Omega_all_inside <- function(a,b,z, n) { # 4.108
  sum(sapply(0:n, function(x) {
    # if (x+(b - 1)/2  < 0) {
    #   return(0)
    # } 
    gamma(a + x)/gamma(a) *gamma(b)/gamma(b+x)*
      Omega_a_c_0_full(a + x, b + x) * z^(x)/factorial(x)
  }
    ))
}
Omega_bar <- function(a,c,z, n) {
  z^(1-c) * Omega_all_inside(a-c+1, 2-c, z, n)
}
struve2 <- function(z, nu_eval, n_approx = 200) {
  if (nu_eval == -1/2) {
    return(sqrt(2/(pi * z)) * sinh(z))
  } else if (nu_eval == -3/2) {
    return(sqrt(2/pi) * (z * cosh(z) - sinh(z))/(z^(3/2)))
  }
  k <- 0:n_approx
  (z/2)^(nu_eval + 1) *sum((z/2)^(2*k) /(gamma(k + 3/2) * gamma( k  + nu_eval + 3/2)))
}
L_by_omega <- function(nu, z, n_approx = 70) {
  exp(-z) * (z/2)^(nu) * Omega_all_inside(a = 1/2 + nu, b = 1 + 2*nu, z = 2*z, n_approx)/
    gamma(nu+1)
}

# test_seq <- val_seq <- seq(0, 5, by = .01)
# for (i in 1:length(test_seq)) {
#   val_seq[i] <- 2^(-2*nu)/gamma(nu+1) *
#     exp(-test_seq[i]) *
#     (2*test_seq[i])^(nu)*
#     Omega(a = 1/2+nu, b = 1 + 2*nu, z =  2 * test_seq[i], n  =80)
# }
nu <- -1.2
test_seq <- val_seq <- seq(0, 5, by = .01)
for (i in 1:length(test_seq)) {
  val_seq[i] <- L_by_omega(nu, test_seq[i],n_approx = 80)
    # 2^(-2*nu)/gamma(nu+1) *
    # exp(-z/2)* z^(-b/2+1/2)*
    # #exp(-z/2)* z^(b/2-1/2)*
    # Omega_bar(a = 1/2+nu, c = 1 + 2*nu, z =  2 * test_seq[i], n  =80)
  #  Omega_all_inside(a = 1/2+nu, b = 1 + 2*nu, z =  2 * test_seq[i], n  =80)
}
plot(test_seq, val_seq, type = 'l')
lines(test_seq, sapply(test_seq, function(x) struve2(z = x, nu_eval = nu, n_approx = 100)),
      col = 2)
plot(val_seq, sapply(test_seq, function(x) struve2(z = x, nu_eval = nu, n_approx = 100)))
plot(test_seq, sapply(test_seq, function(x) struve(z = x, nu_eval = nu)))
lines(test_seq, sapply(test_seq, function(x) struve(z = x, nu_eval = nu)))
plot(val_seq, sapply(test_seq, function(x) struve(z = x, nu_eval = nu)))
abline(a =0, b = 1)
#L_nu(z) = A_{0, nu}(2z) 2^(-2nu)/gamma(nu + 1) /(2z)^(1/2)
#L_nu(z) = 2^(-2nu)/gamma(nu + 1) /(2z)^(1/2) * 
              # * e^(-z) * (2z)^(1/2 + nu) * Omega(1/2 + nu, 1 + 2nu; 2z)

i_seq <- fAsianOptions::whittakerM(x = 2 * test_seq, kappa = nu1/2 - nu2/2, mu = nu1/2 + nu2/2)

full_cov <- (2 * test_seq)^(nu1/2 + nu2/2-1/2)/gamma(1/2 + nu1) *
  (val_seq/gamma(1 - nu1 - nu2)/gamma(1/2 + nu2) - 
     i_seq/gamma(1 + nu1 + nu2)/gamma(1/2 - nu1))

plot(test_seq, full_cov, type = 'l')

plot(i_seq, val_seq, type = 'l')
plot(test_seq, val_seq, type = 'l')

plot(i_seq, 2^(nu1 + nu2 + 1/2) * sqrt(test_seq) *besselI(nu = nu1/2 + nu2/2, x = test_seq)*
       gamma(1 + nu1/2 + nu2/2),
     type = 'l')
plot(2^(nu1 + nu2 + 1/2) * sqrt(test_seq) *
       sapply(test_seq, function(x) struve(nu_eval = -nu1/2- nu2/2, z = x))*
       gamma(1 + nu1/2 + nu2/2),
     2^(nu1 + nu2 + 1/2) * sqrt(test_seq) *besselI(nu = nu1/2 + nu2/2, x = test_seq)*
       gamma(1 + nu1/2 + nu2/2), type = 'l')


