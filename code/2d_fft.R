
norm_constant <- function(nu_1, nu_2, a_1 = 1, a_2 = 1, d = 2) {
  (a_1)^(nu_1) * (a_2)^(nu_2)*
    sqrt(gamma(nu_1 + d/2)) * sqrt(gamma(nu_2 + d/2))/pi^(d/2)/sqrt(gamma(nu_1)*gamma(nu_2))
}
library(Rcpp)
cppFunction("arma::cx_mat armafft(arma::cx_mat x) { return ifft2(x); }",
                         depends="RcppArmadillo")
fft_2d <- function(n_points = 2^10, x_max = 10, nu1 = .5, nu2 = .5, a1 = 1, a2 = 1, 
                   Psi, 
                   Delta, d = 2) {
  delta_t <-  pi / x_max
  x_vals <- 1:n_points * (2 * pi) / n_points / delta_t
  freq_points <- seq(-delta_t * n_points/2 + delta_t/2, delta_t * n_points/2 - 
                       delta_t/2, length.out = n_points)
  freq_grid <- as.data.frame(expand.grid('x' = freq_points, 'y' = freq_points))
  freq_grid$r <- sqrt(freq_grid[['x']]^2 + freq_grid[['y']]^2)
  freq_grid$theta <- atan2(freq_grid[['y']], freq_grid[['x']])
  # https://stackoverflow.com/questions/24077913/discretized-continuous-fourier-transform-with-numpy
  tv <- complex(real = a1, imaginary = Psi(freq_grid[['theta']])*freq_grid[['r']])^(-nu1 - d/2) *
    complex(real = a2, imaginary = -Psi(freq_grid[['theta']])*freq_grid[['r']])^(-nu2 - d/2) * 
    Delta(freq_grid[['theta']])
  tv_mat <- matrix(nrow = length(freq_points), ncol = length(freq_points), tv)
  phase_factor <- 1/(2*pi)  * 
    exp(complex(imaginary = rowSums(cbind(freq_grid[['x']][1], freq_grid[['y']][1]) * 2 * pi *
                  (expand.grid((1:length(freq_points)) /(delta_t * length(freq_points)), 
                               (1:length(freq_points)) /(delta_t * length(freq_points))) ))))
  phase_factor_mat <- matrix(nrow = length(freq_points), ncol = length(freq_points),
                             phase_factor)
  ff_res <- armafft(tv_mat)
  p <- ncol(ff_res)/2
  ff_res_adj <- cbind(rbind(ff_res[(p+1):(2*p),(p+1):(2*p)], ff_res[1:p,(p+1):(2*p)]),
                      rbind(ff_res[(p+1):(2*p),1:p], ff_res[1:p,1:p])) * -phase_factor_mat
  x_vals <- x_vals - (x_vals[length(x_vals)] - x_vals[1]) /2
  cbind(expand.grid(x_vals, x_vals), 'val' = (length(x_vals))^(2)*
          as.double(Re(ff_res_adj) * norm_constant(nu1, nu2, a1, a2)) * 2/pi/x_max^2/ 0.01026171)
}
Delta <- function(x) {
  1
}
Psi <- function(x) {
  sign(x)
}

df <- fft_2d(nu1 =1.2, nu2 =1.2, a1 = 1, a2 = 1, Psi = Psi, Delta = Delta, n_points = 2^10,
             x_max = 10)

var1_vals <- df[['Var1']][df[['Var2']] == min(abs(df[['Var2']]))]
plot(var1_vals, Re(df[['val']][df[['Var2']] == min(abs(df[['Var2']]))]), type = 'l')
lines(var1_vals, col = 2, Matern(abs(var1_vals), nu = 1.2))
df_mat <- matrix(df[['val']], nrow= sqrt(length(df[['val']])))

df <- fft_2d(nu1 =1.5, nu2 =1.2, a1 = 1, a2 = 1, Delta = Delta, Psi = Psi, n_points = 2^9)
var1_vals <- df$Var1[df[['Var2']] == min(abs(df[['Var2']]))]
plot(var1_vals, Re(df[['val']][df[['Var2']] == min(abs(df[['Var2']]))]), type = 'l')
df_mat <- matrix(df[['val']], nrow= sqrt(length(df[['val']])))
image.plot(df_mat)

Delta <- function(x) {
  .97 * complex(imaginary = sign(x))
}
df <- fft_2d(nu1 =1.5, nu2 =1.5, a1 = 1, a2 = 1, Delta = Delta, Psi = Psi, n_points = 2^9)

ggplot(data = df %>%
         filter(abs(Var1) < 5, abs(Var2) < 5), aes(x = Var2, y = Var1, fill = val))+
  geom_raster() + 
  scale_fill_gradient2() +
  labs(x = 'Dimension 1', y = 'Dimension 2',
       fill = 'Cross-\ncovariance') +
  coord_equal() + 
  theme(legend.position = 'left',legend.key.height = unit(.8, "cm")) 

df <- fft_2d(nu1 =.5, nu2 = 1.5, a1 = 1.5, a2 = 1, Delta = Delta, Psi = Psi, n_points = 2^9)

ggplot(data = df %>%
         filter(abs(Var1) < 5, abs(Var2) < 5), aes(x = Var2, y = Var1, fill = val))+
  geom_raster() + 
  scale_fill_gradient2() +
  labs(x = 'Dimension 1', y = 'Dimension 2',
       fill = 'Cross-\ncovariance') +
  coord_equal() + 
  theme(legend.position = 'left',legend.key.height = unit(.8, "cm")) 

