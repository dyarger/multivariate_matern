
library(ggplot2)
library(scales)
theme_set(theme_bw())
library(fields)
library(dplyr)

spatial_integrate_d2 <- function(h, a_1, a_2, nu_1, nu_2, Delta, Psi= Psi, approx_grid,
                                 d = 2) {
  theta_x <-  approx_grid[['theta_x']]
  theta_y <-  approx_grid[['theta_y']]
  r <-  approx_grid[['r']]
  Psi_val <- Psi(theta_x = theta_x, theta_y = theta_y)
  complex_r <- complex(imaginary = r)
  values <- exp(complex_r*(h[1] * theta_x + h[2] * theta_y)) *
    Delta(theta_x, theta_y, 1, 2) *
    (a_1 + Psi_val * complex_r)^(-nu_1-d/2)*
    (a_2 - Psi_val * complex_r)^(-nu_2-d/2) * r^(d-1) 
  Re(sum(values*approx_grid[['angle_lag']]*approx_grid[['r_lag']]
         # * approx_grid[['r']]
  ))
}

norm_constant <- function(d,nu_1, nu_2, a_1 = 1, a_2 = 1, norm_type = 'D') {
  if (norm_type == 'A') {
    (a_1 + a_2)^(nu_1 + nu_2)  /2/pi/gamma(nu_1 + nu_2) * gamma(nu_1 + d/2) * gamma(nu_2 + d/2)
  } else if (norm_type == 'B') {
    (2*a_1)^(nu_1) * (2*a_2)^(nu_2) /2/pi/sqrt(gamma(2*nu_1)*gamma(2 * nu_2)) * 
      gamma(nu_1 + d/2) * gamma(nu_2 + d/2)
  } else if (norm_type == 'C') {
    a_1^(nu_1 + 1/2) * a_2^(nu_2 + 1/2) /2/pi
  } else if (norm_type == 'D') {
    (a_1)^(nu_1) * (a_2)^(nu_2)*
      sqrt(gamma(nu_1 + d/2)) * sqrt(gamma(nu_2 + d/2))/pi^(d/2)/sqrt(gamma(nu_1)*gamma(nu_2))
  }
}

# lag_seq <- seq(-1, 1, by = .1)
# lag_seq <- seq(-7, 7, by = .04)
lag_seq <- seq(-4, 4, by = .2)
grid <- expand.grid(lag_seq, lag_seq)

#angle_grid <- seq(0, 2 * pi, length.out = 40)
#r_grid <- exp(seq(-9, 10, length.out = 40))

angle_grid <- seq(0, 2 * pi, length.out = 200)
#r_grid <- exp(seq(-9, 10, length.out = 300))
r_grid <- exp(seq(-9, 8, length.out = 800))
r_grid <- exp(seq(-12, 12, length.out = 1000))
r_lag <- data.frame('r' = r_grid,
                    'r_lag' = c(r_grid[2] - r_grid[1], (r_grid - dplyr::lag(r_grid))[-1]))
angle_lag <- data.frame('angle' = angle_grid,
                        'angle_lag' = c(angle_grid[2] - angle_grid[1], (angle_grid - dplyr::lag(angle_grid))[-1]))
approx_grid <- expand.grid('angle' = angle_grid,
                           'r' = r_grid) %>%
  left_join(r_lag) %>% left_join(angle_lag)

approx_grid$x <- approx_grid$r * cos(approx_grid$angle)
approx_grid$y <- approx_grid$r * sin(approx_grid$angle)
approx_grid$theta_x <- cos(approx_grid$angle)
approx_grid$theta_y <- sin(approx_grid$angle)

Delta <- function(theta_x, theta_y, entry_1, entry_2) {
  1
}
Psi <- function(theta_x, theta_y) {
  sign(theta_x)
}
nu_1 <- .7
nu_2 <- 1.7
res <- sapply(1:nrow(grid), function(x, a_1, a_2, nu_1, nu_2, Delta, Psi,  approx_grid, 
                                     d = 2) {
  h = as.double(grid[x,])
  spatial_integrate_d2(h = h, a_1, a_2, nu_1 = nu_1, nu_2 = nu_2, Delta = Delta, Psi = Psi, approx_grid = approx_grid)
}, a_1 = 1, a_2 = 1, nu_1 = nu_1, nu_2 = nu_2, 
Delta = Delta, Psi = Psi, approx_grid = approx_grid, d = 2
)
mat <- matrix(unlist(res) * norm_constant(d =2, nu_1 = nu_1, nu_2 = nu_2, a_1 = 1, a_2 = 1, norm_type = 'D'),
              length(lag_seq), length(lag_seq))
image.plot(mat)

df <- data.frame(grid, res = unlist(res) * norm_constant(d =2, nu_1 = nu_1, nu_2 = nu_2, a_1 = 1, a_2 = 1, norm_type = 'D'))

ggplot(data = df, aes(x = Var1, y = Var2, fill = res))+
  geom_raster() + 
  scale_fill_gradientn(
    colors = rainbow(10)) +
  labs(x = expression(t[1]), y = expression(t[2]),
       fill = 'Cross Covariance') +
  theme(legend.position = 'top',legend.key.width = unit(.8, "cm")) 

source('code/multi_matern_source.R')
plot(df$Var1[df$Var2 == 0], df$res[df$Var2 == 0])
points(col = 2,-df$Var1[df$Var2 == 0], sapply(df$Var1[df$Var2 == 0], function(x) {whitt_only_single(x,
                                                nu1 = nu_1, nu2 = nu_2, a1 = 1, a2 = 1, realp = 1, 
                                                imp = 0, norm_type = 'D', which_val = 2)}))

plot(df$Var2[df$Var1 == 0], df$res[df$Var1 == 0])
points(col = 2, df$Var2[df$Var1 == 0], Matern(abs(df$Var2)[df$Var1 == 0], nu = nu_1/2 + nu_2/2))


library(tidyverse)
source('code/multi_matern_source.R')
theme_set(theme_bw())
plot_cc <- function(plot_seq =  seq(-5, 5, by = .05), cons, nu) {
  params <- par()$mfrow
  cov_val_lag <- sapply(1:length(plot_seq), function(x) cross_cov(0, plot_seq[x], nu = nu, z_ij = cons))
  plot(plot_seq, cov_val_lag, type = 'l', xlab = 'lag', main = 'Cross Covariance')
}

plot_function <- function(h,nu, a = 1) {
  if(h == 0) {
    return(0)
  }
  sign(h) * (abs(h)/a)^nu*
    (besselI(a*abs(h), nu = nu) - struve(a*abs(h), -nu, n_approx = 100))
}

norm_constant <- function(nu_1, nu_2, a_1 = 1, a_2 = 1, norm_type = 'A') {
  if (norm_type == 'A') {
    (a_1 + a_2)^(nu_1 + nu_2)  /2/pi/gamma(nu_1 + nu_2) * gamma(nu_1 + 1/2) * gamma(nu_2 + 1/2)
  } else if (norm_type == 'B') {
    (2*a_1)^(nu_1) * (2*a_2)^(nu_2) /2/pi/sqrt(gamma(2*nu_1)*gamma(2 * nu_2)) * 
      gamma(nu_1 + 1/2) * gamma(nu_2 + 1/2)
  } else if (norm_type == 'C') {
    a_1^(nu_1 + 1/2) * a_2^(nu_2 + 1/2) /2/pi
  } else if (norm_type == 'D') {
    (a_1)^(nu_1) * (a_2)^(nu_2)*
      sqrt(gamma(nu_1 + 1/2)) * sqrt(gamma(nu_2 + 1/2))/pi^(1/2)/sqrt(gamma(nu_1)*gamma(nu_2))
  }
}

struve <- function(z, nu_eval, n_approx = 200) {
  if (nu_eval == -1/2) {
    return(sqrt(2/(pi * z)) * sinh(z))
  } else if (nu_eval == -3/2) {
    return(sqrt(2/pi) * (z * cosh(z) - sinh(z))/(z^(3/2)))
  }
  k <- 0:n_approx
  (z/2)^(nu_eval + 1) *sum((z/2)^(2*k) /(gamma(k + 3/2) * gamma( k  + nu_eval + 3/2)))
}


full_function <- function(plot_seq, nu, a, norm_type = 'A') {
  -sapply(1:length(plot_seq), function(x) plot_function(h = plot_seq[x], nu= nu, a = a))* 
    2^(-nu) * pi^(3/2)/cos(nu * pi)/gamma(nu + .5)*
    norm_constant(nu_1 = nu, nu_2 = nu, a_1 = a, a_2 = a, norm_type = norm_type) 
}
nu_vals <- c(1, 4, 6, 8, 10, 50, 90)
plot_seq <- seq(-20, 20, by = .05)
fun_vals <- sapply(nu_vals, full_function, plot_seq = plot_seq,
                   a = 1, norm_type = 'D')
df <- data.frame(plot_seq, 
                 value = as.double(fun_vals),
                 nu = rep(nu_vals, each = length(plot_seq)))

ggplot(data = df, aes(x = plot_seq, y = value, color = factor(nu), linetype = factor(nu))) + 
  geom_line() + 
  labs(x = 'Lags', y = 'Cross-covariance function value', color = expression(nu),
       linetype = expression(nu))

