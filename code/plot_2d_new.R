
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
  if (norm_type == 'B') {
    (a_1 + a_2)^(nu_1 + nu_2)  /2/pi/gamma(nu_1 + nu_2) * gamma(nu_1 + d/2) * gamma(nu_2 + d/2)
  } else if (norm_type == 'A') {
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
#lag_seq <- seq(-2, 2, by = .01)
lag_seq <- seq(-3, 3, by = .02)
grid <- expand.grid(lag_seq, lag_seq)

#angle_grid <- seq(0, 2 * pi, length.out = 40)
#r_grid <- exp(seq(-9, 10, length.out = 40))

angle_grid <- seq(0, 2 * pi, length.out = 100)
#r_grid <- exp(seq(-9, 10, length.out = 300))
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
# }, a_1 = 1, a_2 = 1, nu_1 = 1.5, nu_2 = 1.5,
# Delta = Delta, Psi = Psi, approx_grid = approx_grid, d = 2
# )
# mat <- matrix(unlist(res)*norm_constant(d =2, nu_1 = 1.5, nu_2 = 1.5, a_1 = 1, a_2 = 1, norm_type = 'D')
#                 , length(lag_seq), length(lag_seq))
# plot(lag_seq, mat[11,], ylim = c(.4,1))
# points(col = 2, lag_seq[lag_seq >= 0], Matern(lag_seq[lag_seq >=0], nu = 1.5))
# mat[11,]
# image.plot(mat)
ggplot(data = data.frame(approx_grid),
       aes(x = x, y = y, color = theta_y))+
  geom_point()+scale_color_gradient2()

# approx_seq <- exp(seq(-10, 7, length.out = 300))
# approx_seq <- seq(-100, 100, length.out = 300)
Delta <- function(theta_x, theta_y, entry_1, entry_2) {
  1
}
Psi <- function(theta_x, theta_y) {
  sign(theta_x)
}
# res <- sapply(1:nrow(grid), function(x, a_1, a_2, nu_1, nu_2, Delta, Psi,  approx_grid, 
#                                      d = 2) {
#   h = as.double(grid[x,])
#   spatial_integrate_d2(h = h, a_1, a_2, nu_1 = nu_1, nu_2 = nu_2, Delta = Delta, Psi = Psi, approx_grid = approx_grid)
# }, a_1 = 1, a_2 = 1, nu_1 = .5, nu_2 = .5, 
# Delta = Delta, Psi = Psi, approx_grid = approx_grid, d = 2
# )
# mat <- matrix(unlist(res), length(lag_seq), length(lag_seq))
# image.plot(mat)

res <- sapply(1:nrow(grid), function(x, a_1, a_2, nu_1, nu_2, Delta, Psi,  approx_grid, 
                                     d = 2) {
  h = as.double(grid[x,])
  spatial_integrate_d2(h = h, a_1, a_2, nu_1 = nu_1, nu_2 = nu_2, Delta = Delta, Psi = Psi, approx_grid = approx_grid)
}, a_1 = 1, a_2 = 1, nu_1 = 0.8, nu_2 = 1.2, 
Delta = Delta, Psi = Psi, approx_grid = approx_grid, d = 2
)
mat <- matrix(unlist(res) * norm_constant(d =2, nu_1 = .8, nu_2 = 1.2, a_1 = 1, a_2 = 1, norm_type = 'D'),
              length(lag_seq), length(lag_seq))
image.plot(mat)

df <- data.frame(grid, res = unlist(res) * norm_constant(d =2, nu_1 = .8, nu_2 = 1.2, a_1 = 1, a_2 = 1, norm_type = 'D'))

ggplot(data = df, aes(x = Var1, y = Var2, fill = res))+
  geom_raster() + 
  scale_fill_gradientn(colors = rev(rainbow(10))) +
  labs(x = 'Dimension 1', y = 'Dimension 2',
       fill = 'Cross-\ncovariance') +
  coord_equal() + 
  theme(legend.position = 'left',legend.key.height = unit(.8, "cm")) 
ggsave('images/cc_fun_2d_7.png', height = 4, width = 5.1, dpi = 150)

res <- sapply(1:nrow(grid), function(x, a_1, a_2, nu_1, nu_2, Delta, Psi,  approx_grid, 
                                     d = 2) {
  h = as.double(grid[x,])
  spatial_integrate_d2(h = h, a_1, a_2, nu_1 = nu_1, nu_2 = nu_2, Delta = Delta, Psi = Psi, approx_grid = approx_grid)
}, a_1 = .8, a_2 = 1.2, nu_1 = 1, nu_2 = 1, 
Delta = Delta, Psi = Psi, approx_grid = approx_grid, d = 2
)
mat <- matrix(unlist(res)* norm_constant(d =2, nu_1 = 1, nu_2 = 1, a_1 = .8, a_2 = 1.2, norm_type = 'D'), length(lag_seq), length(lag_seq))
image.plot(mat)

df <- data.frame(grid, res = unlist(res)* norm_constant(d =2, nu_1 = 1, nu_2 = 1, a_1 = .8, a_2 = 1.2, norm_type = 'D'))
ggplot(data = df, aes(x = Var1, y = Var2, fill = res))+
  geom_raster() + 
  scale_fill_gradientn(colors = rev(rainbow(10))) +
  #scale_fill_viridis_b(breaks = round(seq(min(unlist(res)), max(unlist(res)), length.out = 21), 2)) +
  labs(x = 'Dimension 1', y = 'Dimension 2',
       fill = 'Cross-\ncovariance') +
  coord_equal() + 
  theme(legend.position = 'left',legend.key.height = unit(.8, "cm")) 
ggsave('images/cc_fun_2d_9.png', height = 4, width = 5.1, dpi = 150)

res <- sapply(1:nrow(grid), function(x, a_1, a_2, nu_1, nu_2, Delta, Psi,  approx_grid, 
                                     d = 2) {
  h = as.double(grid[x,])
  spatial_integrate_d2(h = h, a_1, a_2, nu_1 = nu_1, nu_2 = nu_2, Delta = Delta, Psi = Psi, approx_grid = approx_grid)
}, a_1 = .8, a_2 = 1.2, nu_1 = .8, nu_2 = 1.2, 
Delta = Delta, Psi = Psi, approx_grid = approx_grid, d = 2
)
mat <- matrix(unlist(res) * norm_constant(d =2, nu_1 = .8, nu_2 = 1.2, a_1 = .8, a_2 = 1.2, norm_type = 'D'), length(lag_seq), length(lag_seq))
image.plot(mat)

df <- data.frame(grid, res = unlist(res)* norm_constant(d =2, nu_1 = .8, nu_2 = 1.2, a_1 = .8, a_2 = 1.2, norm_type = 'D'))
ggplot(data = df, aes(x = Var1, y = Var2, fill = res))+
  geom_raster() + 
  scale_fill_gradientn(#breaks = round(seq(min(unlist(res)), max(unlist(res)), length.out = 21), 2),
                       colors = rev(rainbow(10))) +
  #scale_fill_viridis_b(breaks = round(seq(min(unlist(res)), max(unlist(res)), length.out = 21), 2)) +
  labs(x = 'Dimension 1', y = 'Dimension 2',
       fill = 'Cross-\ncovariance') +
  coord_equal() + 
  theme(legend.position = 'left',legend.key.height = unit(.8, "cm")) 
  # guides(fill = guide_legend(keywidth = .5, direction = 'vertical',
  #                            ncol = 21, label.position = 'bottom', byrow = T,
  #                            label.theme = element_text(angle = 45, size = 6)))
ggsave('images/cc_fun_2d_11.png', height = 4, width = 5.1, dpi = 150)


lag_seq <- seq(-5, 5, by = .05)
grid <- expand.grid(lag_seq, lag_seq)

angle_grid <- seq(0, 2 * pi, length.out = 200)
r_grid <- exp(seq(-9, 10, length.out = 700))
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
  complex(imaginary = sign(theta_x) * .97)
}

res <- sapply(1:nrow(grid), function(x, a_1, a_2, nu_1, nu_2, Delta, Psi,  approx_grid, 
                                     d = 2) {
  h = as.double(grid[x,])
  spatial_integrate_d2(h = h, a_1, a_2, nu_1 = nu_1, nu_2 = nu_2, Delta = Delta, Psi = Psi, approx_grid = approx_grid)
}, a_1 = 1, a_2 = 1, nu_1 = 1.5, nu_2 = 1.5, 
Delta = Delta, Psi = Psi, approx_grid = approx_grid, d = 2
)
mat <- matrix(unlist(res)* norm_constant(d =2, nu_1 = 1.5, nu_2 = 1.5, a_1 = 1, a_2 = 1, norm_type = 'D'), length(lag_seq), length(lag_seq))
image.plot(mat)

df <- data.frame(grid, res = unlist(res)* norm_constant(d =2, nu_1 = 1.5, nu_2 = 1.5, a_1 = 1, a_2 = 1, norm_type = 'D'))
ggplot(data = df, aes(x = Var1, y = Var2, fill = res))+
  geom_raster() + 
  scale_fill_gradient2() +
  labs(x = 'Dimension 1', y = 'Dimension 2',
       fill = 'Cross-\ncovariance') +
  coord_equal() + 
  theme(legend.position = 'left',legend.key.height = unit(.8, "cm")) 
ggsave('images/cc_fun_2d_1.png', height = 4, width = 5.1, dpi = 150)
plot(df$Var1[df$Var2 == 0], df$res[df$Var2 == 0])


res <- sapply(1:nrow(grid), function(x, a_1, a_2, nu_1, nu_2, Delta, Psi,  approx_grid, 
                                     d = 2) {
  h = as.double(grid[x,])
  spatial_integrate_d2(h = h, a_1, a_2, nu_1 = nu_1, nu_2 = nu_2, Delta = Delta, Psi = Psi, approx_grid = approx_grid)
}, a_1 = 1, a_2 = 1, nu_1 = .4, nu_2 = 2.5, 
Delta = Delta, Psi = Psi, approx_grid = approx_grid, d = 2
)
mat <- matrix(unlist(res)* norm_constant(d =2, nu_1 = .4, nu_2 = 2.5, a_1 = 1, a_2 = 1, norm_type = 'D'), length(lag_seq), length(lag_seq))
image.plot(mat)

df <- data.frame(grid, res = unlist(res)* norm_constant(d =2, nu_1 = .4, nu_2 = 2.5, a_1 = 1, a_2 = 1, norm_type = 'D'))
ggplot(data = df, aes(x = Var1, y = Var2, fill = res))+
  geom_raster() + 
  scale_fill_gradient2() +
  labs(x = 'Dimension 1', y = 'Dimension 2',
       fill = 'Cross-\ncovariance') +
  coord_equal() + 
  theme(legend.position = 'left',legend.key.height = unit(.8, "cm")) 
ggsave('images/cc_fun_2d_3.png', height = 4, width = 5.1, dpi = 150)

Psi <- function(theta_x, theta_y) {
  sign(theta_x * .3 + theta_y * .7)
}

res <- sapply(1:nrow(grid), function(x, a_1, a_2, nu_1, nu_2, Delta, Psi,  approx_grid, 
                                     d = 2) {
  h = as.double(grid[x,])
  spatial_integrate_d2(h = h, a_1, a_2, nu_1 = nu_1, nu_2 = nu_2, Delta = Delta, Psi = Psi, approx_grid = approx_grid)
}, a_1 = 1, a_2 = 1, nu_1 = .4, nu_2 = 2.5, 
Delta = Delta, Psi = Psi, approx_grid = approx_grid, d = 2
)
mat <- matrix(unlist(res)* norm_constant(d =2, nu_1 = .4, nu_2 = 2.5, a_1 = 1, a_2 = 1, norm_type = 'D'), length(lag_seq), length(lag_seq))
image.plot(mat)

df <- data.frame(grid, res = unlist(res)* norm_constant(d =2, nu_1 = .4, nu_2 = 2.5, a_1 = 1, a_2 = 1, norm_type = 'D'))
ggplot(data = df, aes(x = Var1, y = Var2, fill = res))+
  geom_raster() + 
  scale_fill_gradient2() +
  labs(x = 'Dimension 1', y = 'Dimension 2',
       fill = 'Cross-\ncovariance') +
  coord_equal() + 
  theme(legend.position = 'left',legend.key.height = unit(.8, "cm")) 
ggsave('images/cc_fun_2d_3_psi.png', height = 4, width = 5.1, dpi = 150)


Delta <- function(theta_x, theta_y, entry_1, entry_2) {
  complex(imaginary = sign(theta_x) * sign(theta_y) * .97)
}

res <- sapply(1:nrow(grid), function(x, a_1, a_2, nu_1, nu_2, Delta, Psi,  approx_grid, 
                                     d = 2) {
  h = as.double(grid[x,])
  spatial_integrate_d2(h = h, a_1, a_2, nu_1 = nu_1, nu_2 = nu_2, Delta = Delta, Psi = Psi, approx_grid = approx_grid)
}, a_1 = 1, a_2 = 1, nu_1 = 1.5, nu_2 = 1.5, 
Delta = Delta, Psi = Psi, approx_grid = approx_grid, d = 2
)
mat <- matrix(unlist(res)* norm_constant(d =2, nu_1 = 1.5, nu_2 = 1.5, a_1 = 1, a_2 = 1, norm_type = 'D'), length(lag_seq), length(lag_seq))
image.plot(mat)

df <- data.frame(grid, res = unlist(res)* norm_constant(d =2, nu_1 = 1.5, nu_2 = 1.5, a_1 = 1, a_2 = 1, norm_type = 'D'))
ggplot(data = df, aes(x = Var1, y = Var2, fill = res))+
  geom_raster() + 
  scale_fill_gradient2() +
  labs(x = 'Dimension 1', y = 'Dimension 2',
       fill = 'Cross-\ncovariance') +
  coord_equal() + 
  theme(legend.position = 'left',legend.key.height = unit(.8, "cm")) 
  # theme(legend.position = 'top',legend.key.width = unit(.8, "cm"), 
  #       legend.text = element_text(angle = 45) )
ggsave('images/cc_fun_2d_5.png', height = 4, width = 5.1, dpi = 150)







  









# Compare with incomplete function representation
J_incomplete <- function(nu, h, c) {
  f <- function(theta) {
    cos(h * cos(theta)) * sin(theta)^(2 * nu)
  }
  int_val <- integrate(f = f, lower = 0, upper = c, subdivisions = 500)
  int_val[['value']] * 
    2 * h^nu/ (2^nu * gamma(nu + 1/2) * gamma(1/2))
}

H_incomplete <- function(nu, h, c) {
  f <- function(theta) {
    sin(h * cos(theta)) * sin(theta)^(2 * nu)
  }
  int_val <- integrate(f = f, lower = 0, upper = c, subdivisions = 500)
  int_val[['value']] * 
    2 * h^nu/ (2^nu * gamma(nu + 1/2) * gamma(1/2))
}
spatial_integrate_d2_bessel <- function(h, a_1, a_2, nu_1, nu_2, approx_grid,
                                 d = 2) {
  r <-  approx_grid[['r']]
#  Psi_val <- Psi(theta_x = theta_x, theta_y = theta_y)
  # J_values <- sapply(1:nrow(approx_grid), function(ind) {J_incomplete(nu = d/2-1,
  #                                                                     h = r[ind], 
  #                                                                     c = abs(atan2(h[2], h[1])))})
  norm_h <- sqrt(h[1]^2 + h[2]^2)
  # print(ifelse(abs(atan2(h[2], h[1])) > pi/2,
  #              abs(atan2(h[2], h[1])) - pi/2,
  #              abs(atan2(h[2], h[1]))))
  J_values <- sapply(1:length(r), function(ind) {besselJ(nu = d/2-1, x = abs(r)[ind] *norm_h )})
  J_values1 <- sapply(1:length(r), function(ind) {
    J_incomplete(nu = d/2-1, h = abs(r)[ind] *norm_h,
                 c = pi/2)})
  J_values2 <- sapply(1:length(r), function(ind) {
    J_incomplete(nu = d/2-1, h = abs(r)[ind] *norm_h,
                 c = 0)})
  
  H_values <- sapply(1:length(r), function(ind) {H_incomplete(nu = d/2-1, 
                                                              h = abs(r)[ind] *norm_h,
                                                              c = pi/2)})
  H_values1 <- sapply(1:length(r), function(ind) {H_incomplete(nu = d/2-1, 
                                                               h = abs(r)[ind] *norm_h,
                                                               c = 0)})
  H_values2 <- sapply(1:length(r), function(ind) {H_incomplete(nu = d/2-1, 
                                                               h = abs(r)[ind] *norm_h,
                                                               c = pi/2)})
  values <- norm_h^(-d/2 + 1) *
    # complex(real = 2*J_values - J_values1 - J_values2,
    #         imaginary = 2*H_values - H_values1 - H_values2)*
    complex(real = J_values,
            imaginary = H_values)*
    (a_1 + complex(imaginary = r))^(-nu_1-d/2)*
    (a_2 -  complex(imaginary = r))^(-nu_2-d/2) * abs(r)^(d/2) 
  Re(sum(values*approx_grid[['r_lag']]))
}
r_grid2 <-  exp(seq(-8, 8, length.out = 150))
approx_grid2 <- data.frame(r = sort(c(-r_grid2,r_grid2)))
approx_grid2[['r_lag']] <- c(approx_grid2$r[2] - approx_grid2$r[1],
                             (approx_grid2$r - lag(approx_grid2$r))[-1])
rm(a_1)
rm(a_2)
rm(nu_1)
rm(nu_2)
lag_seq <- seq(-1, 1, by = .1)
#lag_seq <- seq(-1, 1, by = .08)
grid <- expand.grid(lag_seq, lag_seq)


res2 <- sapply(1:nrow(grid), function(x, a_1, a_2, nu_1, nu_2,  approx_grid, 
                                     d = 2) {
  h = as.double(grid[x,])
  spatial_integrate_d2_bessel(h = h, a_1, a_2, nu_1 = nu_1, nu_2 = nu_2, 
                              approx_grid = approx_grid)/2
}, a_1 = 1, a_2 = 1, nu_1 = .8, nu_2 = .8, 
approx_grid = approx_grid2, d = 2
)
mat2 <- matrix(unlist(res2), length(lag_seq), length(lag_seq))
image.plot(mat2)
plot(lag_seq, mat2[nrow(mat2)/2,], type = 'l')
lines(lag_seq, Matern(abs(lag_seq), smoothness = .8)/1.6, col = 2)


approx_grid2$r
test_r <- seq(0, 20, by =.01)
rho_values <- seq(0, pi/2, length.out = 10)
J_values <- matrix(nrow = length(test_r), ncol = 10)
for (i in 1:length(rho_values)) {
  J_values[,i] <- sapply(1:length(test_r), function(ind) {
    J_incomplete(nu = 0, h = abs(test_r)[ind],
                 c = rho_values[i])})
}
plot(test_r,  J_values[,1], type = 'l')
lines(test_r,  J_values[,2], type = 'l')
lines(test_r,  J_values[,3], type = 'l')
lines(test_r,  J_values[,4], type = 'l')
lines(test_r,  J_values[,5], type = 'l')
lines(test_r,  J_values[,6], type = 'l')
lines(test_r,  J_values[,7], type = 'l')
lines(test_r,  J_values[,8], type = 'l')
lines(test_r,  J_values[,9], type = 'l')
lines(test_r,  J_values[,10], type = 'l')


# now try more complicated things
Delta <- function(theta_x, theta_y, entry_1, entry_2) {
  complex(imaginary = sign(theta_x) * .3)
}

#lag_seq <- seq(-3, 3, by = .025)
grid <- expand.grid(lag_seq, lag_seq)
res <- lapply(1:nrow(grid), function(x, a_1, a_2, nu_1, nu_2, e_1, e_2,
                                     Delta, approx_seq) {
  h = as.double(grid[x,])
  if (x %% 100 == 0) {
    print(x)
    print(h)
  }
  spatial_integrate_d2(h = h, a_1, a_2, nu_1 = nu_1, nu_2 = nu_2, e_1 = e_1, e_2 = e_2, Delta = Delta, approx_seq = approx_seq)
}, a_1 = a_1, a_2 = a_2, nu_1 = nu_1, nu_2 = nu_2, e_1 = e_1, e_2 = e_2, Delta = Delta,approx_seq = approx_seq
)
res_mat <- matrix( nrow = length(lag_seq), ncol = length(lag_seq), unlist(res))
fields::image.plot(res_mat)

grid_vals <- data.frame(grid, res_vals = as.double(unlist(res)))
#save(grid_vals, approx_seq, file = paste0('results/asymmetric_2d_', length(approx_seq), '.RData'))
ggplot(data = grid_vals) + 
  geom_tile(aes(x = Var1, y = Var2, fill = res_vals)) + 
  scale_fill_gradient2() + 
  labs(x = expression(t[1]), y = expression(t[2]),
       fill = 'Cross\nCovariance') + 
  theme(legend.position = 'bottom')
ggsave(paste0('images/asymmetric_2d_', length(approx_seq), '.png'), height = 4, width = 4)

Delta <- function(theta_x, theta_y, entry_1, entry_2) {
  1
}

lag_seq <- seq(-3, 3, by = .05)
grid3 <- expand.grid(lag_seq, lag_seq, lag_seq)
set.seed(40)
theta_star <- rnorm(3)
theta_star <- theta_star / sqrt(sum(theta_star^2))

grid3_mat <- as.matrix(grid3)
is_int <- sapply(1:nrow(grid3), function(x) {
  sign(sum(grid3_mat[x,] * theta_star))
})

h_val <- c(2, .4, 4)
inner_prod_h <- sapply(1:nrow(grid3), function(x) {
  sum(grid3_mat[x,] * h_val)
})



grid3_info <- cbind(grid3, is_int, inner_prod_h)
a <- ggplot(data = grid3_info %>% filter(Var1 %in% c(-3, 0, 3)), 
       aes(x = Var2, y =Var3, fill = factor(is_int))) + 
  geom_raster() +facet_wrap(~Var1, ncol = 1) + theme(legend.position = 'none')
a
b <- ggplot(data = grid3_info %>%filter(Var1 %in% c(-3, 0, 3)), 
       aes(x = Var2, y =Var3, fill = inner_prod_h)) + 
  geom_raster() + scale_fill_gradient2()+facet_wrap(~Var1, ncol = 1)
b
library(patchwork)
a + b
s
res <- lapply(1:nrow(grid), function(x, a_1, a_2, nu_1, nu_2, e_1, e_2,
                                     Delta, approx_seq) {
  h = as.double(grid[x,])
  if (x %% 100 == 0) {
    print(x)
    print(h)
  }
  spatial_integrate_d2(h = h, a_1, a_2, nu_1 = nu_1, nu_2 = nu_2, e_1 = e_1, e_2 = e_2, Delta = Delta, approx_seq = approx_seq)
}, a_1 = a_1, a_2 = a_2, nu_1 = nu_1, nu_2 = nu_2, e_1 = 1/4, e_2 = e_2, Delta = Delta,approx_seq = approx_seq
)
res_mat <- matrix( nrow = length(lag_seq), ncol = length(lag_seq), unlist(res))
fields::image.plot(res_mat)

grid_vals <- data.frame(grid, res_vals = as.double(unlist(res)))
ggplot(data = grid_vals) + 
  geom_tile(aes(x = Var1, y = Var2, fill = res_vals)) + 
  scale_fill_gradientn(colors = rainbow(5, rev = T)) + 
  labs(x = expression(t[1]), y = expression(t[2]),
       fill = 'Cross\nCovariance') + 
  theme(legend.position = 'bottom')
ggsave(paste0('images/asymmetric_2d_E_', length(approx_seq), '.png'), height = 4, width = 4)
save(grid_vals, approx_seq, file = paste0('results/asymmetric_2d_E_', length(approx_seq), '.RData'))

Delta <- function(theta_x, theta_y, entry_1, entry_2) {
  sign(theta_x)*sign(theta_y) * .3
}

res <- lapply(1:nrow(grid), function(x, a_1, a_2, nu_1, nu_2, e_1, e_2,
                                     Delta, approx_seq) {
  h = as.double(grid[x,])
  if (x %% 100 == 0) {
    print(x)
  }
  spatial_integrate_d2(h = h, a_1, a_2, nu_1 = nu_1, nu_2 = nu_2, e_1 = e_1, e_2 = e_2, Delta = Delta, approx_seq = approx_seq)
}, a_1 = a_1, a_2 = a_2, nu_1 = nu_1, nu_2 = nu_2, e_1 = e_1, e_2 = e_2, Delta = Delta,approx_seq = approx_seq
)
res_mat <- matrix( nrow = length(lag_seq), ncol = length(lag_seq), unlist(res))
fields::image.plot(res_mat)

grid_vals <- data.frame(grid, res_vals = as.double(unlist(res)))
ggplot(data = grid_vals) + 
  geom_tile(aes(x = Var1, y = Var2, fill = res_vals)) + 
  scale_fill_gradient2(midpoint = 0) + 
  labs(x = expression(t[1]), y = expression(t[2]),
       fill = 'Cross\nCovariance') + 
  theme(legend.position = 'bottom')
ggsave(paste0('images/asymmetric_2d_quadrant_', length(approx_seq), '.png'), height = 4, width = 4)
save(grid_vals, approx_seq, file = paste0('results/asymmetric_2d_quadrant_', length(approx_seq), '.RData'))
