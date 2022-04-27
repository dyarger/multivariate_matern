
library(ggplot2)
library(scales)
theme_set(theme_bw())
spatial_integrate <- function(h,d, a_1, a_2, nu_1, nu_2) {
  spec_dens_single <- function(x, h, d, a_1, a_2, nu_1, nu_2) {
    Re(besselJ(h*x, d/2 - 1) * x^(-(d/2 - 1)) * (a_1 + complex(imaginary = x))^(-nu_1-d/2)*
         (a_2 - complex(imaginary = x))^(-nu_2-d/2)) * x^(d-1)
  }
  test_integrate <- integrate(lower = .0001, upper = 1000, 
                              f = spec_dens_single,
                              h = h, d = d, a_1 = a_1, a_2 = a_2, nu_1 = nu_1, 
                              nu_2 = nu_2,subdivisions = 2000)
  (2 * pi)^(d/2) *h^(-d/2 + 1) * test_integrate[['value']]
}
a_1 <- a_2 <- 1
d <- 2


lag_seq <- seq(0, 5, by = .005)
res <- sapply(lag_seq, spatial_integrate, d = d, a_1 = a_1, a_2 = a_2, nu_1 = .5, nu_2 = .5)
res2 <- sapply(lag_seq, spatial_integrate, d = d, a_1 = a_1, a_2 = a_2, nu_1 = 1, nu_2 = .5)
res3 <- sapply(lag_seq, spatial_integrate, d = d, a_1 = a_1, a_2 = a_2, nu_1 = 1, nu_2 = 1)
res4 <- sapply(lag_seq, spatial_integrate, d = d, a_1 = a_1, a_2 = a_2, nu_1 = .75, nu_2 = .75)
df <- data.frame(lag_seq, 
                 value = c(res, res2, res3, res4),
                 nu1 = rep(c(.5, .5, 1, .75), each = length(lag_seq)),
                 nu2 = rep(c(.5, 1, 1, .75), each = length(lag_seq)), 
                 nu_label = rep(c('nu[1]=0.2 nu[2]=0.2', 'nu[1]=0.2~nu[2]=0.5', 
                                  'nu[1]=0.5~nu[2]=0.5',  'nu[1]=0.5~nu[2]=0.5'), each = length(lag_seq)))


spatial_integrate_d2 <- function(h, a_1, a_2, nu_1, nu_2, Delta, approx_seq) {
  spec_dens_theta <- function(x, h, Delta) {
    f_theta <- function(theta, x, h, d, Delta) {
      theta_x <- cos(theta);theta_y <- sin(theta)
      Re(exp(complex(imaginary = h[1]*x* theta_x + h[2]*x * theta_y)) *
           Delta(theta_x, theta_y, 1, 2))
    }
    f_theta_im <- function(theta, x, h, d, Delta) {
      theta_x <- cos(theta);theta_y <- sin(theta)
      Im(exp(complex(imaginary = h[1]*x * theta_x + h[2]*x* theta_y)) * 
           Delta(theta_x, theta_y, 1, 2))
    }
    Re_part <- integrate(lower = 0, upper = 2 * pi, f = f_theta, stop.on.error = F,
                         x = x, h = h, Delta = Delta)$value
    Im_part <- integrate(lower = 0, upper = 2 * pi, f = f_theta_im, stop.on.error = F,
                         x = x, h = h, Delta = Delta)$value
    complex(real = Re_part, imaginary = Im_part)
  }
  spec_dens <- function(x, h, a_1, a_2, nu_1, nu_2, Delta) {
    Re(spec_dens_theta(x, h, Delta) * 
         (a_1 + complex(imaginary = x))^(-nu_1-2/2)*
         (a_2 - complex(imaginary = x))^(-nu_2-2/2)) * x^(d-1)
  }
  approx_seq_lag <- c(approx_seq[1], approx_seq[2:length(approx_seq)] - 
                        approx_seq[1:(length(approx_seq)-1)])
  test_final <- sapply(approx_seq,  FUN = spec_dens,
                       h = h, a_1 = a_1, a_2 = a_2, nu_1 = nu_1, 
                       nu_2 = nu_2, Delta = Delta)
  sum(test_final * approx_seq_lag)
}

nu_2 <- nu_1 <- .5

lag_seq <- seq(-1, 1, by = .1)
grid <- expand.grid(lag_seq, 0)

approx_seq <- exp(seq(-10, 7, length.out = 300))
Delta <- function(theta_x, theta_y, entry_1, entry_2) {
  1
}
# make sure it matches with Matern
res <- sapply(1:nrow(grid), function(x, a_1, a_2, nu_1, nu_2,
                                     Delta, approx_seq) {
  h = as.double(grid[x,])
  spatial_integrate_d2(h = h, a_1, a_2, nu_1 = nu_1, nu_2 = nu_2,  Delta = Delta, approx_seq = approx_seq)
}, a_1 = a_1, a_2 = a_2, nu_1 = nu_1, nu_2 = nu_2, Delta = Delta,approx_seq = approx_seq
)
plot(lag_seq, res)
plot(fields::Matern(abs(lag_seq), nu = .5))
plot(res, fields::Matern(abs(lag_seq), nu = .5))
abline(lm(data.frame(x = res, y =  fields::Matern(abs(lag_seq), nu = .5)), formula = y~x))

# now try more complicated things
# Delta <- function(theta_x, theta_y, entry_1, entry_2) {
#   complex(imaginary = sign(theta_x) * .3)
# }

#lag_seq <- seq(-3, 3, by = .025)
lag_seq <- seq(-3, 3, by = .4)
grid <- expand.grid(lag_seq, lag_seq)
nu_1 <- .5
nu_2 <- 1.5
res <- lapply(1:nrow(grid), function(x, a_1, a_2, nu_1, nu_2, 
                                     Delta, approx_seq) {
  h = as.double(grid[x,])
  if (x %% 100 == 0) {
    print(x)
    print(h)
  }
  spatial_integrate_d2(h = h, a_1, a_2, nu_1 = nu_1, nu_2 = nu_2,Delta = Delta, approx_seq = approx_seq)
}, a_1 = a_1, a_2 = a_2, nu_1 = nu_1, nu_2 = nu_2, Delta = Delta,approx_seq = approx_seq
)
res_mat <- matrix( nrow = length(lag_seq), ncol = length(lag_seq), unlist(res))
fields::image.plot(res_mat)
fields::image.plot(Matrix::symmpart(res_mat))
fields::image.plot(Matrix::skewpart(res_mat))


grid_vals <- data.frame(grid, res_vals = as.double(unlist(res)))
#save(grid_vals, approx_seq, file = paste0('results/asymmetric_2d_', length(approx_seq), '.RData'))
ggplot(data = grid_vals) + 
  geom_tile(aes(x = Var1, y = Var2, fill = res_vals)) + 
  scale_fill_gradient2() + 
  labs(x = expression(t[1]), y = expression(t[2]),
       fill = 'Cross\nCovariance') + 
  theme(legend.position = 'bottom')
#ggsave(paste0('images/asymmetric_2d_', length(approx_seq), '.png'), height = 4, width = 4)
