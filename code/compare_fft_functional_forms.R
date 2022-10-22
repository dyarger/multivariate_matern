

source('code/multi_matern_source.R')
grid_info <- create_grid_info_1d(2^12, 25)

nu_test <- .5
df <- fft_1d(nu1 = .5, nu2 = .5, a1 = 1, a2 = 1, re = 1, im = 0, grid_info = grid_info)
plot(df, type = 'l')
lines(df[,1], exp(-abs(df[,1])), col = 2)
plot(df[,2],exp(-abs(df[,1])), type = 'l')
abline(col = 2, a = 0, b = 1)
mean((abs(df[,2]) - exp(-abs(df[,1])))^2)
mean(exp(-abs(df[,1])) / abs(df[,2]))

library(fields)
df <- fft_1d(nu1 = 1.5, nu2 = 1.5, a1 = 2, a2 = 2, re = 1, im = 0,  grid_info = grid_info)
plot(df, type = 'l')
lines(df[,1], Matern(abs(df[,1]), range = 1/2, smoothness = 1.5), col = 2)
mean((abs(df[,2]) - Matern(abs(df[,1]), range = 1/2, smoothness = 1.5))^2)
plot(df[,2], Matern(abs(df[,1]), range = 1/2, smoothness = 1.5), type = 'l')
abline(col = 2, a = 0, b = 1)


# 2d
# test it out
grid_info <- create_grid_info_2d(n_points = 2^9, x_max = 10)
df <- fft_2d(nu1 = 1.5, nu2 = 1.5, a1 = 1, a2 = 1, Psi = Psi, Delta = Delta, grid_info = grid_info)
library(fields)
var1_vals <- df[['Var1']][df[['Var2']] == min(abs(df[['Var2']]))]
plot(var1_vals, Re(df[['val']][df[['Var2']] == min(abs(df[['Var2']]))]), type = 'l')
lines(var1_vals, col = 2, Matern(abs(var1_vals), nu = 1.5))
plot(var1_vals, Re(df[['val']][df[['Var2']] == min(abs(df[['Var2']]))]) - 
       Matern(abs(var1_vals), nu = 1.5), type = 'l')
df_mat <- matrix(df[['val']], nrow = sqrt(length(df[['val']])))
image.plot(df_mat)

# compare fft, spectral density, and functional forms

lag_seq <- seq(-3, 3, by = .02)
grid <- expand.grid(lag_seq, lag_seq)
angle_grid <- seq(0, 2 * pi, length.out = 100)
r_grid <- exp(seq(-9, 8, length.out = 800))
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
# res <- sapply(1:nrow(grid), function(x, a_1, a_2, nu_1, nu_2, Delta, Psi,  approx_grid, 
#                                      d = 2) {
#   h = as.double(grid[x,])
#   spatial_integrate_d2(h = h, a_1, a_2, nu_1 = nu_1, nu_2 = nu_2, Delta = Delta, Psi = Psi, approx_grid = approx_grid)
# }, a_1 = 1, a_2 = 1, nu_1 = 0.8, nu_2 = 1.2, 
# Delta = Delta, Psi = Psi, approx_grid = approx_grid, d = 2
# )
# mat <- matrix(unlist(res) * norm_constant(d =2, nu_1 = .8, nu_2 = 1.2, a_1 = 1, a_2 = 1, norm_type = 'D'),
#               length(lag_seq), length(lag_seq))
# image.plot(mat)
# df <- data.frame(grid_info$x_vals_eg, res)
# df <- data.frame(grid, res = unlist(res) * norm_constant(d =2, nu_1 = .8, nu_2 = 1.2, a_1 = 1, a_2 = 1, norm_type = 'D'))
