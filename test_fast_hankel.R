
besselJ(h*x, d/2 - 1) *
  Re((a_1 + y)^(-nu_1-d/2)*(a_2 - y)^(-nu_2-d/2)) *
  x^(d/2)
f_fun <- function(x, a_1, a_2, nu_1, nu_2, d = 2) {
  y <- complex(imaginary = x)
  Re((a_1 + y)^(-nu_1-d/2)*(a_2 - y)^(-nu_2-d/2)) *
    x^(d/2-1)
}
h <- 0
a_1 = 1
a_2 = 1
nu_1 = .6
nu_2 = .6

alpha <- 2
rho0 <- 2
r0 <- 2
N <- 100
n <- 0:(N-1)
m <- 0
rho = rho0 * exp(alpha * (0:(N-1)))
r = r0 * exp(alpha * (0:(N-1)))
fn = r * f_fun(r, a_1 = 1, a_2 = 1, nu_1 = .6, nu_2 = .6)
jnm <- 2 * pi * alpha





