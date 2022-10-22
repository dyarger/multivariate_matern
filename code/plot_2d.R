
library(ggplot2)
library(scales)
theme_set(theme_bw())
library(fields)
library(dplyr)

source('code/multi_matern_source.R')
Delta <- function(theta_x, theta_y, entry_1, entry_2) {
  1
}
Psi <- function(theta_x, theta_y) {
  sign(theta_x)
}
grid_info <- create_grid_info_2d(n_points = 2^12, x_max = 45)

df <- fft_2d(grid_info, nu1 = .5, nu2 = 1.4, a1 = 1, a2 = 1, Psi = Psi, d = 2, Delta = Delta)
ggplot(data = filter(df, abs(Var1) < 5, abs(Var2) < 5), aes(x = Var2, y = Var1, fill = val)) +
  geom_raster() + 
  scale_fill_gradientn(colors = rev(rainbow(10))) +
  labs(x = 'Dimension 1', y = 'Dimension 2',
       fill = 'Cross-\ncovariance') +
  coord_equal() + 
  theme(legend.position = 'left',legend.key.height = unit(.8, "cm")) 
ggsave('images/cc_fun_2d_7.png', height = 4, width = 5.1, dpi = 150)

df <- fft_2d(grid_info, nu1 = 1, nu2 = 1, a1 = 1.5, a2 = .75, Psi = Psi, d = 2, Delta = Delta)
ggplot(data = filter(df, abs(Var1) < 5, abs(Var2) < 5), aes(x = Var2, y = Var1, fill = val)) +
  geom_raster() + 
  scale_fill_gradientn(colors = rev(rainbow(10))) +
  labs(x = 'Dimension 1', y = 'Dimension 2',
       fill = 'Cross-\ncovariance') +
  coord_equal() + 
  theme(legend.position = 'left',legend.key.height = unit(.8, "cm")) 
ggsave('images/cc_fun_2d_9.png', height = 4, width = 5.1, dpi = 150)

df <- fft_2d(grid_info, nu1 = .5, nu2 = 1.5, a1 = 1.5, a2 = .75, Psi = Psi, d = 2, Delta = Delta)
ggplot(data = filter(df, abs(Var1) < 5, abs(Var2) < 5), aes(x = Var2, y = Var1, fill = val)) +
  geom_raster() + 
  scale_fill_gradientn(
                       colors = rev(rainbow(10))) +
  labs(x = 'Dimension 1', y = 'Dimension 2',
       fill = 'Cross-\ncovariance') +
  coord_equal() + 
  theme(legend.position = 'left',legend.key.height = unit(.8, "cm")) 
ggsave('images/cc_fun_2d_11.png', height = 4, width = 5.1, dpi = 150)

Delta <- function(theta_x, theta_y, entry_1, entry_2) {
  complex(imaginary = sign(theta_x) * .97)
}
df <- fft_2d(grid_info, nu1 = 1.5, nu2 = 1.5, a1 = 1, a2 = 1, Psi = Psi, d = 2, Delta = Delta)
ggplot(data = filter(df, abs(Var1) < 5, abs(Var2) < 5), aes(x = Var2, y = Var1, fill = val)) +
  geom_raster() + 
  scale_fill_gradient2() +
  labs(x = 'Dimension 1', y = 'Dimension 2',
       fill = 'Cross-\ncovariance') +
  coord_equal() + 
  theme(legend.position = 'left',legend.key.height = unit(.8, "cm")) 
ggsave('images/cc_fun_2d_1.png', height = 4, width = 5.1, dpi = 150)
plot(df$Var1[df$Var2 == 0], df$val[df$Var2 == 0], type = 'l')
plot(df$Var2[df$Var1 == 0], df$val[df$Var1 == 0], type = 'l')

df <- fft_2d(grid_info, nu1 = .4, nu2 = 2.5, a1 = 1, a2 = 1, Psi = Psi, d = 2, Delta = Delta)
ggplot(data = filter(df, abs(Var1) < 5, abs(Var2) < 5), aes(x = Var2, y = Var1, fill = val)) +
  geom_raster() + 
  scale_fill_gradient2() +
  labs(x = 'Dimension 1', y = 'Dimension 2',
       fill = 'Cross-\ncovariance') +
  coord_equal() + 
  theme(legend.position = 'left',legend.key.height = unit(.8, "cm")) 
ggsave('images/cc_fun_2d_3.png', height = 4, width = 5.1, dpi = 150)

Delta <- function(theta) {
  complex(real = ifelse(theta < -pi/2, 1, 
                        ifelse(theta < 0, -1, 
                               ifelse(theta < pi/2, 1, -1))) * .97)
}
df <- fft_2d(grid_info, nu1 = 1.5, nu2 = 1.5, a1 = 1, a2 = 1, Psi = Psi, d = 2, Delta = Delta)
ggplot(data = filter(df, abs(Var1) < 5, abs(Var2) < 5), aes(x = Var2, y = Var1, fill = val)) +
  geom_raster() + 
  scale_fill_gradient2() +
  labs(x = 'Dimension 1', y = 'Dimension 2',
       fill = 'Cross-\ncovariance') +
  coord_equal() + 
  theme(legend.position = 'left', legend.key.height = unit(.8, "cm")) 
ggsave('images/cc_fun_2d_5.png', height = 4, width = 5.1, dpi = 150)