
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

Omega <- function(a,b,z, n) {
  sum(sapply(0:n, function(x) 
    gamma(a + x)/gamma(a) *gamma(b)/gamma(b+x)*Omega_a_c_0(a + x, b + x)))
  
}
sapply(0:n, function(x) 
  Omega_a_c_0(a + x, b + x))
Omega(1, 1.2, .5, 90)

A_fun <- function(kappa, mu, z, n) {
  omega_results <- sapply()
  
  exp(-1/2 * z) *  z^(1/2 + mu)
}

test <- seq(0, 5, by = .1)
test_vals <- sapply(test, omega, a = 3, b = 2.2, n = 100)
