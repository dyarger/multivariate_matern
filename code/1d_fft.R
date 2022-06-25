
norm_constant <- function(nu_1, nu_2, a_1 = 1, a_2 = 1) {
  (a_1)^(nu_1) * (a_2)^(nu_2)*
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
  x_vals <- x_vals - (x_vals[length(x_vals)] - x_vals[1]) /2
  cbind(x_vals, Re(ff_res_adj) * norm_constant(nu1, nu2, a1, a2) * (5/2/x_max))
}
nu_test <- .5
df <- fft_1d(nu1 =nu_test, nu2 =nu_test, a1 = 1, a2 = 1, re = 1, im = 0, n_points = 2^16,x_max = 10)
plot(df, type = 'l')
lines(df[,1], fields::Matern(abs(df[,1]), nu = nu_test), col = 2)

a_test <- 2
df <- fft_1d(nu1 =nu_test, nu2 =nu_test, a1 = a_test, a2 = a_test, re = 1, im = 0, 
             n_points = 2^15,x_max = 10)
plot(df, type = 'l')
lines(df[,1], fields::Matern(abs(df[,1]), nu = nu_test, range = 1/a_test), col = 2)


df <- fft_1d(nu1 = .8, nu2 =.8, a1 = 1, a2 = 1, re = 0, im = .5, n_points = 2^16,x_max = 50)
plot(df, type = 'l')
points(0,0)


df <- fft_1d(nu1 =2.5, nu2 =.8, a1 = 1, a2 = 1, re = 1, im = 0, n_points = 2^16,x_max = 15)
plot(df, type = 'l')
points(0,0)
