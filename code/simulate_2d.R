library(fields)
library(ggplot2)
library(tidyverse)
theme_set(theme_bw()+ theme(text = element_text(size = 16), legend.position = 'bottom'))
n_samples <- 2000
dim_grid <- 501

x <- seq(0, 25, length.out = dim_grid)
grid <- matrix(ncol = 2, unlist(expand.grid(x, x)))

# proposal
g_new <- function(h, a = 1, nu = .5, d = 2) {
  1/((a^2 + h^2)^(nu + d/2)) * a^(2 * nu) * gamma(nu + d/2) / gamma(nu) / pi^(d/2)
}

# spectral density
h_fun2 <- function(h, a1, a2, nu1, nu2, d = 2) {
  if (h[1] > 0) {
    complex(imaginary = 0.97) * (a1 + complex(imaginary = sqrt(sum(h^2))))^(-nu1 - d/2) *
      (a2 + complex(imaginary = -sqrt(sum(h^2))))^(-nu2 - d/2) *
      a1^(nu1) * a2^(nu2) * sqrt(gamma(nu1 + d/2)*gamma(nu2 + d/2)) / sqrt(gamma(nu1)*gamma(nu2)) / pi^(d/2)
  } else {
    complex(imaginary = -0.97) * (a1 + complex(imaginary = -sqrt(sum(h^2))))^(-nu1 - d/2) *
      (a2 + complex(imaginary = sqrt(sum(h^2))))^(-nu2 - d/2) *
      a1^(nu1) * a2^(nu2) * sqrt(gamma(nu1 + d/2)*gamma(nu2 + d/2)) / sqrt(gamma(nu1)*gamma(nu2)) / pi^(d/2)
  }
}
a <- 1

# implementation of Emery et al 2016
simu_bivariate <- function(grid, L = 1000, nu1 = 1.5, nu2 = 1.5, a1 = 1, a2 = 1, g_new, cc_fun) {
  sim_field <- matrix(nrow = nrow(grid), ncol = 2, 0)
  for (p in 1:2) {
    for (l in 1:L) {
      phi <- runif(n = 1, 0, 2 * pi)
      sim_vals <- rnorm(n = 2)/sqrt(rgamma(n = 1, shape = .5))
      test_g <- g_new(sqrt(sum(sim_vals^2)), a = a1, nu = nu1)
      test_g12 <- cc_fun(sim_vals, a1 = a1, a2 = a2, nu1 = nu1, nu2 = nu2)
      test_g22 <- g_new(sqrt(sum(sim_vals^2)), a = a2, nu = nu2)
      
      right_mat <- matrix(2*c(test_g, test_g12, Conj(test_g12), test_g22), nrow = 2, ncol = 2) /
        g_new(sqrt(sum(sim_vals^2)), nu = .5)
      eig_mat <- eigen(right_mat)
      sqrt_mat <- eig_mat$vectors %*% diag(sqrt(eig_mat$values))

      sqrt_mat_re <- Re(sqrt_mat)
      sqrt_mat_im <- -Im(sqrt_mat)
      
      ap <- sqrt_mat_re[,p]
      bp <- sqrt_mat_im[,p]
      # do vectorized version
      ap_mat <- matrix(nrow = nrow(sim_field),ncol = 2, byrow = T, ap)
      sim_vals_mat <- matrix(nrow = nrow(sim_field),ncol = 2, byrow = T, sim_vals)
      bp_mat <- matrix(nrow = nrow(sim_field),ncol = 2, byrow = T, bp)
      sim_field <- sim_field + 1/sqrt(L) * ap_mat * cos(rowSums(grid * sim_vals_mat) + phi) +
        1/sqrt(L) * bp_mat * sin(rowSums(grid * sim_vals_mat) + phi)  
    }
  }
  sim_field
}

set.seed(40)
nu1 <- 1.5; nu2 <- 1.5
sim_f <- simu_bivariate(grid, L = n_samples, nu1 = nu1, nu2 = nu2, 
                        g_new = g_new, cc_fun = h_fun2)
var(sim_f)
cor(sim_f)
var(cbind(sim_f[,1], sim_f[,2]))
image.plot(matrix(sim_f[,1], nrow = length(x)), axes = F)
image.plot(matrix(sim_f[,2] , nrow = length(x)), axes = F)

ggplot(data = cbind(as.data.frame(grid), sim1 = sim_f[,1], sim2 = sim_f[,2]) %>%
         tidyr::pivot_longer(cols = starts_with('sim')) %>%
         left_join(data.frame(name = c('sim1', 'sim2'), label = c('Process 1', 'Process 2'))) %>%
         mutate(`Dimension 1` = V1, `Dimension 2` = V2),
       aes(x = `Dimension 1`, y = `Dimension 2`, fill = value)) + 
  facet_wrap(~label, ncol = 2, labeller = ) + 
  geom_raster() +
  coord_equal() +
  scale_fill_gradientn(colors = rev(rainbow(10))) + 
  labs(fill = 'Simulated\nvalue') +
  theme(strip.background = element_blank(),
        strip.text.x = element_blank(),
        legend.key.height = unit(.8, "cm"),
        legend.key.width = unit(.8, "cm")
  )
ggsave(filename = 'images/simu_1_ggplot.png', height = 5.1, width = 8, dpi = 150)

set.seed(45)
sim_f <- simu_bivariate(grid, L = n_samples, nu1 = .4, nu2 = 2.5, 
                        g_new = g_new, cc_fun = h_fun2)

ggplot(data = cbind(as.data.frame(grid), sim1 = sim_f[,1], sim2 = sim_f[,2]) %>%
         tidyr::pivot_longer(cols = starts_with('sim')) %>%
         left_join(data.frame(name = c('sim1', 'sim2'), label = c('Process 1', 'Process 2'))) %>%
         mutate(`Dimension 1` = V1, `Dimension 2` = V2),
       aes(x = `Dimension 1`, y = `Dimension 2`, fill = value)) + 
  facet_wrap(~label, ncol = 2) + 
  geom_raster() +
  coord_equal() +
  scale_fill_gradientn(colors = rev(rainbow(10)))  + 
  labs(fill = 'Simulated\nvalue') +
  theme(strip.background = element_blank(),
        strip.text.x = element_blank(),
        legend.key.height = unit(.8, "cm"),
        legend.key.width = unit(.8, "cm")
  )
ggsave(filename = 'images/simu_3_ggplot.png', height = 5.1, width = 8, dpi = 150)

h_fun3 <- function(h, a1, a2, nu1, nu2, d = 2) {
  if (h[1] > 0) {
    if (h[2] > 0) {
      complex(real = 0.97)/(a1 + complex(imaginary = sqrt(sum(h^2))))^(nu1 + d/2) / 
        (a2 + complex(imaginary = -sqrt(sum(h^2))))^(nu2 + d/2) *
        a1^(nu1) * a2^(nu2) * sqrt(gamma(nu1 + d/2)*gamma(nu2 + d/2)) / sqrt(gamma(nu1)*gamma(nu2)) / pi^(d/2)
    } else {
      complex(real = -0.97)/(a1 + complex(imaginary = sqrt(sum(h^2))))^(nu1 + d/2) / 
        (a2 + complex(imaginary = -sqrt(sum(h^2))))^(nu2 + d/2) *
        a1^(nu1) * a2^(nu2) * sqrt(gamma(nu1 + d/2)*gamma(nu2 + d/2)) / sqrt(gamma(nu1)*gamma(nu2)) / pi^(d/2)
    }
  } else {
    if (h[2] > 0) {
      complex(real = -0.97)/(a1 + complex(imaginary = -sqrt(sum(h^2))))^(nu1 + d/2) / 
        (a2 + complex(imaginary = sqrt(sum(h^2))))^(nu2 + d/2) *
        a1^(nu1) * a2^(nu2) * sqrt(gamma(nu1 + d/2)*gamma(nu2 + d/2)) / sqrt(gamma(nu1)*gamma(nu2)) / pi^(d/2)
    } else {
      complex(real = 0.97)/(a1 + complex(imaginary = -sqrt(sum(h^2))))^(nu1 + d/2) / 
        (a2 + complex(imaginary = sqrt(sum(h^2))))^(nu2 + d/2) *
        a1^(nu1) * a2^(nu2) * sqrt(gamma(nu1 + d/2)*gamma(nu2 + d/2)) / sqrt(gamma(nu1)*gamma(nu2)) / pi^(d/2)
    }
  }
}

set.seed(22)

sim_f <- simu_bivariate(grid, L = n_samples, nu1 = 1.5, nu2 = 1.5,
                        g_new = g_new, cc_fun = h_fun3)
ggplot(data = cbind(as.data.frame(grid), sim1 = sim_f[,1], sim2 = sim_f[,2]) %>%
         tidyr::pivot_longer(cols = starts_with('sim')) %>%
         left_join(data.frame(name = c('sim1', 'sim2'), label = c('Process 1', 'Process 2'))) %>%
         mutate(`Dimension 1` = V1, `Dimension 2` = V2),
       aes(x = `Dimension 1`, y = `Dimension 2`, fill = value)) + 
  facet_wrap(~label, ncol = 2) + 
  geom_raster() +
  coord_equal() +
  scale_fill_gradientn(colors = rev(rainbow(10)))  + 
  labs(fill = 'Simulated\nvalue') +
  theme(strip.background = element_blank(),
        strip.text.x = element_blank(),
        legend.key.height = unit(.8, "cm"),
        legend.key.width = unit(.8, "cm")
  )
ggsave(filename = 'images/simu_5_ggplot.png', height = 5.1, width = 8, dpi = 150)



h_fun4 <- function(h, a1, a2, nu1, nu2, d = 2) {
  if (h[1] > 0) {
    complex(real = 0.85)/(1 + complex(imaginary = a1 * sqrt(sum(h^2))))^(nu1 + d/2) / 
      (1 + complex(imaginary = -a2 * sqrt(sum(h^2))))^(nu2 + d/2) *
      a1^(nu1) * a2^(nu2) * sqrt(gamma(nu1 + d/2)*gamma(nu2 + d/2)) / sqrt(gamma(nu1)*gamma(nu2)) / pi^(d/2)
  } else {
    complex(real = 0.85)/(a1 + complex(imaginary = -sqrt(sum(h^2))))^(nu1 + d/2) / 
      (a2 + complex(imaginary = sqrt(sum(h^2))))^(nu2 + d/2) *
      a1^(nu1) * a2^(nu2) * sqrt(gamma(nu1 + d/2)*gamma(nu2 + d/2)) / sqrt(gamma(nu1)*gamma(nu2)) / pi^(d/2)
  }
}
set.seed(22^2)
sim_f <- simu_bivariate(grid, L = n_samples,
                        nu1 = .5, nu2 = 1.5, g_new = g_new, cc_fun = h_fun4)

ggplot(data = cbind(as.data.frame(grid), sim1 = sim_f[,1], sim2 = sim_f[,2]) %>%
         tidyr::pivot_longer(cols = starts_with('sim')) %>%
         left_join(data.frame(name = c('sim1', 'sim2'), label = c('Process 1', 'Process 2'))) %>%
         mutate(`Dimension 1` = V1, `Dimension 2` = V2),
       aes(x = `Dimension 1`, y = `Dimension 2`, fill = value)) + 
  facet_wrap(~label, ncol = 2) + 
  geom_raster() +
  coord_equal() +
  scale_fill_gradientn(colors = rev(rainbow(10)))  + 
  labs(fill = 'Simulated\nvalue') +
  theme(strip.background = element_blank(),
        strip.text.x = element_blank(),
        legend.key.height = unit(.8, "cm"),
        legend.key.width = unit(.8, "cm")
  )
ggsave(filename = 'images/simu_7_ggplot.png', height = 5.1, width = 8, dpi = 150)



set.seed(23^2)
sim_f <- simu_bivariate(grid, L = n_samples,
                        nu1 = 1, nu2 = 1, a1 = 1.2, a2 = .8,
                        g_new = g_new, cc_fun = h_fun4)

ggplot(data = cbind(as.data.frame(grid), sim1 = sim_f[,1], sim2 = sim_f[,2]) %>%
         tidyr::pivot_longer(cols = starts_with('sim')) %>%
         left_join(data.frame(name = c('sim1', 'sim2'), label = c('Process 1', 'Process 2'))) %>%
         mutate(`Dimension 1` = V1, `Dimension 2` = V2),
       aes(x = `Dimension 1`, y = `Dimension 2`, fill = value)) + 
  facet_wrap(~label, ncol = 2, labeller = ) + 
  geom_raster() +
  coord_equal() +
  scale_fill_gradientn(colors = rev(rainbow(10)))  + 
  labs(fill = 'Simulated\nvalue') +
  theme(strip.background = element_blank(),
        strip.text.x = element_blank(),
        legend.key.height = unit(.8, "cm"),
        legend.key.width = unit(.8, "cm")
  )
ggsave(filename = 'images/simu_9_ggplot.png', height = 5.1, width = 8, dpi = 150)



set.seed(24^2)
sim_f <- simu_bivariate(grid, L = n_samples, nu1 = .5, nu2 = 1.5, a1 = 1.1, a2 = .9,
                        g_new = g_new, cc_fun = h_fun4)

ggplot(data = cbind(as.data.frame(grid), sim1 = sim_f[,1], sim2 = sim_f[,2]) %>%
         tidyr::pivot_longer(cols = starts_with('sim')) %>%
         left_join(data.frame(name = c('sim1', 'sim2'), label = c('Process 1', 'Process 2'))) %>%
         mutate(`Dimension 1` = V1, `Dimension 2` = V2),
       aes(x = `Dimension 1`, y = `Dimension 2`, fill = value)) + 
  facet_wrap(~label, ncol = 2, labeller = ) + 
  geom_raster() +
  coord_equal() +
  scale_fill_gradientn(colors = rev(rainbow(10)))  + 
  labs(fill = 'Simulated\nvalue') +
  theme(strip.background = element_blank(),
    strip.text.x = element_blank(),
    legend.key.height = unit(.8, "cm"),
    legend.key.width = unit(.8, "cm")
  )
ggsave(filename = 'images/simu_11_ggplot.png', height = 5.1, width = 8, dpi = 150)




