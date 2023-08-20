
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

labels <- labs(x = 'Dimension 1', y = 'Dimension 2',
               fill = 'Cross-\ncovariance')
theme_use <- theme(legend.position = 'left',legend.key.height = unit(.8, "cm")) 
color_scale <- scale_fill_gradientn(colors = rev(rainbow(10)))
color_scale2 <- scale_fill_gradientn(colors = rev(rainbow(10)))

df <- fft_2d(grid_info, nu1 = .5, nu2 = 1.5, a1 = 1, a2 = 1, Psi = Psi, d = 2, Delta = Delta)
ggplot(data = filter(df, abs(Var1) < 5, abs(Var2) < 5), aes(x = Var2, y = Var1, fill = val)) +
  geom_raster() + coord_equal() + 
  color_scale + labels + theme_use
ggsave('images/cc_fun_2d_7.png', height = 4, width = 5.1, dpi = 150)

df <- fft_2d(grid_info, nu1 = 1, nu2 = 1, a1 = 1.2, a2 = .8, Psi = Psi, d = 2, Delta = Delta)
ggplot(data = filter(df, abs(Var1) < 5, abs(Var2) < 5), aes(x = Var2, y = Var1, fill = val)) +
  geom_raster() + coord_equal() + 
  color_scale + labels + theme_use
ggsave('images/cc_fun_2d_9.png', height = 4, width = 5.1, dpi = 150)

df <- fft_2d(grid_info, nu1 = .5, nu2 = 1.5, a1 = 1.1, a2 = .9, Psi = Psi, d = 2, Delta = Delta)
ggplot(data = filter(df, abs(Var1) < 5, abs(Var2) < 5), aes(x = Var2, y = Var1, fill = val)) +
  geom_raster() + coord_equal() + 
  color_scale + labels + theme_use
ggsave('images/cc_fun_2d_11.png', height = 4, width = 5.1, dpi = 150)

Delta <- function(theta_x, theta_y, entry_1, entry_2) {
  complex(imaginary = sign(theta_x) * .97)
}
df <- fft_2d(grid_info, nu1 = 1.5, nu2 = 1.5, a1 = 1, a2 = 1, Psi = Psi, d = 2, Delta = Delta)
ggplot(data = filter(df, abs(Var1) < 5, abs(Var2) < 5), aes(x = Var2, y = Var1, fill = val)) +
  geom_raster() + coord_equal() + 
  scale_fill_gradient2() + labels + theme_use
ggsave('images/cc_fun_2d_1.png', height = 4, width = 5.1, dpi = 150)
plot(df$Var1[df$Var2 == 0], df$val[df$Var2 == 0], type = 'l')
plot(df$Var2[df$Var1 == 0], df$val[df$Var1 == 0], type = 'l')

df <- fft_2d(grid_info, nu1 = .4, nu2 = 2.5, a1 = 1, a2 = 1, Psi = Psi, d = 2, Delta = Delta)
ggplot(data = filter(df, abs(Var1) < 5, abs(Var2) < 5), aes(x = Var2, y = Var1, fill = val)) +
  geom_raster() + coord_equal() + 
  scale_fill_gradient2() + labels + theme_use
ggsave('images/cc_fun_2d_3.png', height = 4, width = 5.1, dpi = 150)

Delta <- function(theta) {
  complex(real = ifelse(theta < -pi/2, 1, 
                        ifelse(theta < 0, -1, 
                               ifelse(theta < pi/2, 1, -1))) * .97)
}
df <- fft_2d(grid_info, nu1 = 1.5, nu2 = 1.5, a1 = 1, a2 = 1, Psi = Psi, d = 2, Delta = Delta)
ggplot(data = filter(df, abs(Var1) < 5, abs(Var2) < 5), aes(x = Var2, y = Var1, fill = val)) +
  geom_raster() + coord_equal() + 
  scale_fill_gradient2() + labels + theme_use 
ggsave('images/cc_fun_2d_5.png', height = 4, width = 5.1, dpi = 150)

Psi <- function(theta_x, theta_y) {
  sign(theta_x)
}
Delta <- function(theta) {
  complex(real = ifelse(theta < -pi/2, 1, 
                        ifelse(theta < 0, 0.25, 
                               ifelse(theta < pi/2, 1, 0.25))))
}
df1 <- fft_2d(grid_info, nu1 = 1.5, nu2 = 1.5, a1 = 1, a2 = 1, Psi = Psi, d = 2, Delta = Delta) %>%
  filter(abs(Var1) < 5, abs(Var2) < 5)
ggplot(data = filter(df1, abs(Var1) < 5, abs(Var2) < 5), aes(x = Var2, y = Var1, fill = val)) +
  geom_raster() + coord_equal() + 
  color_scale2 + labels + theme_use 
ggsave('images/cc_fun_2d_aniso.png', height = 4, width = 5.1, dpi = 150)

Delta <- function(theta) {
  complex(real = ifelse(theta < -pi/2, 1, 
                        ifelse(theta < 0, 1, 
                               ifelse(theta < pi/2, 1, 1))))
}
df1 <- fft_2d(grid_info, nu1 = 1.5, nu2 = 1.5, a1 = 1, a2 = 1, Psi = Psi, d = 2, Delta = Delta) %>%
  filter(abs(Var1) < 5, abs(Var2) < 5)
ggplot(data = filter(df1, abs(Var1) < 5, abs(Var2) < 5), aes(x = Var2, y = Var1, fill = val)) +
  geom_raster() + coord_equal() + 
  color_scale2 + labels + theme_use 
ggsave('images/cc_fun_2d_aniso_compare.png', height = 4, width = 5.1, dpi = 150)




