
library(ggplot2)
library(scales)
theme_set(theme_bw())
library(fields)
library(dplyr)
source('code/multi_matern_source.R')
grid_info_1d <- create_grid_info_1d(2^14, x_max = 80) 
df_1d <- fft_1d(nu1 = .5, nu2 = .5, a1 = 1, a2 = 1, re = 0, im = 1, grid_info = grid_info_1d)

plot(df_1d[abs(df_1d[,1]) < 10,], type = 'l')

grid_info_2d <- create_grid_info_2d(2^12, x_max = 50) 
df_2d <- fft_2d(nu1 = .5, nu2 = .5, a1 = 1, a2 = 1, d = 2, Psi = sign, 
                Delta =  function(x) complex(real = 0,
                                             imaginary = sign(x)), grid_info = grid_info_2d)

min_abs_val <- min(abs(df_2d$Var2))

plot(df_1d[abs(df_1d[,1]) < 10,], type = 'l')
lines(df_2d$Var2[df_2d$Var1 == min_abs_val], df_2d$val[df_2d$Var1 == min_abs_val], col = 2)

ggplot(data = filter(df_2d, abs(Var1) < 10, abs(Var2) < 10), aes(x = Var2, y = Var1, fill = val)) +
  geom_raster() + 
  scale_fill_gradient2() +
  labs(x = expression(t[1]), y = expression(t[2]),
       fill = 'Cross Covariance') +
  theme(legend.position = 'top',legend.key.width = unit(.8, "cm")) 

df_1d <- fft_1d(nu1 = 1.5, nu2 = .5, a1 = 1, a2 = 1, re = 0, im = 1, grid_info = grid_info_1d)
df_2d <- fft_2d(nu1 = 1.5, nu2 = .5, a1 = 1, a2 = 1, d = 2, Psi = sign, 
                Delta =  function(x) complex(real = 0,
                                             imaginary = sign(x)), grid_info = grid_info_2d)
plot(df_1d[abs(df_1d[,1]) < 10,], type = 'l')
lines(df_2d$Var2[df_2d$Var1 == min_abs_val], df_2d$val[df_2d$Var1 == min_abs_val], col = 2)
# does not hold when nu1\neq nu2
# which explains why the real directional measures do not match up except when Matern

df_1d <- fft_1d(nu1 = .5, nu2 = .5, a1 = 1.2, a2 = 1, re = 0, im = 1, grid_info = grid_info_1d)
df_2d <- fft_2d(nu1 = .5, nu2 = .5, a1 = 1.2, a2 = 1, d = 2, Psi = sign, 
                Delta =  function(x) complex(real = 0,
                                             imaginary = sign(x)), grid_info = grid_info_2d)
plot(df_1d[abs(df_1d[,1]) < 10,], type = 'l')
lines(df_2d$Var2[df_2d$Var1 == min_abs_val], df_2d$val[df_2d$Var1 == min_abs_val], col = 2)

# compare real directional measure - Matern
df_1d <- fft_1d(nu1 = .5, nu2 = .5, a1 = 1, a2 = 1, re = 1, im = 0, grid_info = grid_info_1d)
df_2d <- fft_2d(nu1 = .5, nu2 = .5, a1 = 1, a2 = 1, d = 2, Psi = sign, 
                Delta =  function(x) 1, grid_info = grid_info_2d)
plot(df_1d[abs(df_1d[,1]) < 10,], type = 'l')
lines(df_2d$Var2[df_2d$Var1 == min_abs_val], df_2d$val[df_2d$Var1 == min_abs_val], col = 2)

# compare real directional measure - mismatched nu
df_1d <- fft_1d(nu1 = 1.5, nu2 = .5, a1 = 1, a2 = 1, re = 1, im = 0, grid_info = grid_info_1d)
df_2d <- fft_2d(nu1 = 1.5, nu2 = .5, a1 = 1, a2 = 1, d = 2, Psi = sign, 
                Delta =  function(x) 1, grid_info = grid_info_2d)
plot(df_1d[abs(df_1d[,1]) < 10,], type = 'l')
lines(df_2d$Var2[df_2d$Var1 == min_abs_val], df_2d$val[df_2d$Var1 == min_abs_val], col = 2)

