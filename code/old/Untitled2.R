

# fractional brownian motion

cov_fbm <- function(s,t,h) {
  .5 * (abs(t)^(2*h) + abs(s)^(2*h) - abs(t-s)^(2*h))
}


test_seq <- seq(-5, 5, .1)
test_grid <- expand.grid(test_seq, test_seq)
h = .7
results <- sapply(1:nrow(test_grid), function(index) {cov_fbm(test_grid[index,1],
                                                              test_grid[index,2],
                                                              h = h)})
fields::image.plot(matrix(results, nrow = length(test_seq)))


cov_ofbm <- function(s,t,h1,h2, rho, eta) {
  .5 * ((rho + eta * sign(s))  * abs(s)^(h1 + h2) +
          (rho - eta * sign(t)) * abs(t)^(h1 + h2) - 
          (rho - eta * sign(t-s)) * abs(t-s)^(h1 + h2))
}
h1 = .2
h2 = .9
rho = .5
eta = .8
results <- sapply(1:nrow(test_grid), function(index) {cov_ofbm(test_grid[index,1],
                                                               test_grid[index,2],
                                                               h1 = h1, h2 = h2,
                                                               rho = rho, eta = eta)})
fields::image.plot(matrix(results, nrow = length(test_seq)))

plot(abs(test_grid[,1] - test_grid[,2]),results)


library(fAsianOptions)
nu1 = 1.2
nu2 = 2.2
h <- seq(-5, 5, by = .011)
test <- sapply(h, function(x) {
  if (x > 0) {
    1/gamma(nu2 + 1/2) * abs(x)^(nu1/2 + nu2/2 - 1/2)* fAsianOptions::whittakerW(x = 2*abs(x), 
                                                                                 kappa = -nu1/2 + nu2/2,
                                                            mu = - (nu1+nu2)/2 , ip = 0)
  } else {
    1/gamma(nu1 + 1/2) * abs(x)^(nu1/2 + nu2/2 - 1/2)* fAsianOptions::whittakerW(x = 2*abs(x),
                                                            kappa = nu1/2 - nu2/2,
                                                            mu = - (nu1+nu2)/2 , ip = 0)
  }
})
plot(h, Re(test), type = 'l')

A <- function(kappa, mu, z, list_vals) {
  exp(-1/2 * z) * z^(1/2 + mu) * Omega(1/2 - kappa + mu, 1 + 2 * mu, z, list_vals[[1]],
                                                list_vals[[2]], list_vals[[3]])
}
nu1 <- .3
nu2 <- .3
nu_minus <- nu1/2-nu2/2
nu_plus <- nu1/2+nu2/2
n_series <- 500
list_vals <- make_b_vals(kappa = nu_minus, mu = nu_plus, n_series = n_series)
test <- sapply(seq(0,10, by = .1), function(y) y^(nu_plus- 1/2)* A(z = y, 0, .2,
                                                                   list_vals = list_vals))
test2 <- sapply(seq(0,10, by = .1), function(y) y^(nu_plus- 1/2) * whittakerM(x = y, kappa = nu_minus, mu =- nu_plus))

plot(seq(0, 10, by = .1), test)
lines(seq(0, 10, by = .05),4 * ( besselI(x = seq(0, 10, by = .05), nu = nu1+nu2) -
        RandomFieldsUtils::struveL(x = seq(0, 10, by = .05), nu = -( nu1+nu2))))
plot(seq(0, 10, by = .1), test2)

make_b_vals <- function(kappa, mu, n_series = n_series) {
  a = 1/2 - kappa - mu
  c = 1+2 * mu
  seq1 <- c-a
  seq2 <- a + 0:(n_series-1)
  b1 <- sapply(1:(n_series), function(x) b_half(seq1, seq2[x]))
  b2 <- sapply(1:(n_series), function(x) b_half(seq2[x], seq1))
  bvals <- cbind(b1,b2)
  gamma_val <- gamma(c)/gamma(a)/gamma(c-a)
  factorials <- sapply(0:n_series, factorial)
  return(list(bvals, gamma_val, factorials))
}

b_half <- function(a, b) {
  pbeta(1/2,a,b) * beta(a,b)
}

Omega <- function(a, c, z, bvals, gamma_val,factorials) {
  series_val <- sapply(0:(nrow(bvals)-1), function(n) {
    z^n/factorials[n+1]
  })
  gamma_val * sum(series_val * (bvals[,1] - bvals[,2]))
}

Omega_a_c_0 <- function(a,c) { 
  if (a > c | a < 0) {
    stop('need c > a > 0')
  }
  gamma(c)/gamma(a)/gamma(c-a) * (b_half(c-a, a) - b_half(a, c-a))
}
Omega_a_c_0(.2, .3)
Omega_a_c_0(.5, .3)

Omega_a_c_0_a <- function(a,c) {
  a_prime <- a
  a_final <- a - floor(a)
  iter <- a_final - a
  a_seq <- a_final:a
  if (c < a_final) {
    stop('need c bigger')
  }
  
  i= 1
  cons <- 2^(2-c) * gamma(c)/gamma(a_seq[i])/gamma(c-a_seq[i]+1)
  omega_init <- Omega_a_c_0(a_final, c) - cons
  omega_final <- omega_init- cons
  for (i in 2:iter) {
    cons <- 2^(2-c) * gamma(c)/gamma(a_seq[i])/gamma(c-a_seq[i]+1)
    omega_final <- omega_final- cons
    #print(omega_final)
  }
  omega_final
}

# c > a
# 1 + 2mu > 1/2 - kappa - mu
# 1/2 + mu > -kappa
# mu > -kappa - 1/2
# nu_+ > - nu_- - 1/2
# nu_+ + nu_- > - 1/2

# c > a > 0
# 1/2 - kappa - mu > 0
# kappa + mu < 1/2 
# nu_+ + nu_- < 1/2

# |nu_+ + nu_-| < 1/2
# |nu1/2 - nu2/2 + nu1/2 + nu2/2| < 1/2
# |nu1| < 1/2


