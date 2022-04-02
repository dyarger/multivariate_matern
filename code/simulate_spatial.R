
x <- seq(0, 25, length.out = 301)

grid <- matrix(ncol = 2, unlist(expand.grid(x, x)))
library(fields)
g <- function(h, a, nu, d = 2) {
  a^d * gamma(nu + d/2)/gamma(nu)/pi^(d/2)/(1 + a^2 * h^2)^(nu + d/2)
}

g_new <- function(h, a, nu, d = 2) {
  1/(1 + a^2 * h^2)^(nu + d/2)
}

h_fun <- function(h, a1, a2, nu1, nu2, d = 2) {
  1/(1 + complex(imaginary = a1 * h))^(nu1 + d/2)/ 
    (1 + complex(imaginary = -a2 * h))^(nu2 + d/2)
}

h_fun2 <- function(h, a1, a2, nu1, nu2, d = 2) {
  if (h[1] > 0) {
    1/(1 + complex(imaginary = a1 * sqrt(sum(h^2))))^(nu1 + d/2)/ 
      (1 + complex(imaginary = -a2 * sqrt(sum(h^2))))^(nu2 + d/2)
  } else {
    -1/(1 + complex(imaginary = a1 * sqrt(sum(h^2))))^(nu1 + d/2)/ 
      (1 + complex(imaginary = -a2 * sqrt(sum(h^2))))^(nu2 + d/2)
  }
}

h_fun2 <- function(h, a1, a2, nu1, nu2, d = 2) {
  if (h[1] > 0) {
    complex(imaginary = 0.97)/(1 + complex(imaginary = a1 * sqrt(sum(h^2))))^(nu1 + d/2)/ 
      (1 + complex(imaginary = -a2 * sqrt(sum(h^2))))^(nu2 + d/2)
  } else {
    complex(imaginary = -0.97)/(1 + complex(imaginary = a1 * sqrt(sum(h^2))))^(nu1 + d/2)/ 
      (1 + complex(imaginary = -a2 * sqrt(sum(h^2))))^(nu2 + d/2)
  }
}


nu1 <- 1.5
nu2 <- 1.5
a <- 1

simu_bivariate <- function(grid, L = 1000, nu1 = 1.5, nu2 = 1.5, g_new, cc_fun) {
  sim_field <- matrix(nrow = nrow(grid), ncol = 2, 0)
  for (p in 1:2) {
    for(l in 1:L) {
      phi <- runif(n = 1, 0, 2 * pi)
      sim_vals <- rnorm(2)/sqrt(rgamma(n = 1,shape = .5))
      test_g <- g_new(sqrt(sum(sim_vals^2)), a = 1, nu = nu1)
      #test_g12 <- .9*h_fun(sqrt(sum(sim_vals^2)), a1 = 1, a2 = 1, nu1 = nu1, nu2 = nu2)
      test_g12 <- cc_fun(sim_vals, a1 = 1, a2 = 1, nu1 = nu1, nu2 = nu2)
      test_g22 <- g_new(sqrt(sum(sim_vals^2)), a = 1, nu = nu2)
      
      right_mat <- matrix(2*c(test_g, test_g12, test_g12, test_g22), nrow = 2, ncol = 2) /
        g(sqrt(sum(sim_vals^2)), a = 1, nu = .5)
      eig_mat <- eigen(right_mat, symmetric = T)
      sqrt_mat <- eig_mat$vectors %*% diag(sqrt(eig_mat$values))
      sqrt_mat_re <- Re(sqrt_mat)
      sqrt_mat_im <- Im(sqrt_mat)
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
sim_f <- simu_bivariate(grid, L = 1000, nu1 = 1.5, nu2 = 1.5, g_new = g_new, cc_fun = h_fun2)
image.plot(matrix(sim_f[,1], nrow = length(x)), axes = F)
image.plot(matrix(sim_f[,2], nrow = length(x)), axes = F)

png('simu_spatial1.png', height = 550, width = 550)
image.plot(matrix(sim_f[,1], nrow = length(x)), axes = F)
dev.off()

png('simu_spatial2.png', height = 550, width = 550)
image.plot(matrix(sim_f[,2], nrow = length(x)), axes = F)
dev.off()

set.seed(45)
sim_f <- simu_bivariate(grid, L = 1000, nu1 = .4, nu2 = 2.5, g_new = g_new, cc_fun = h_fun2)
image.plot(matrix(sim_f[,1], nrow = length(x)), axes = F)
image.plot(matrix(sim_f[,2], nrow = length(x)), axes = F)

png('simu_spatial3.png', height = 550, width = 550)
image.plot(matrix(sim_f[,1], nrow = length(x)), axes = F)
dev.off()

png('simu_spatial4.png', height = 550, width = 550)
image.plot(matrix(sim_f[,2], nrow = length(x)), axes = F)
dev.off()


h_fun3 <- function(h, a1, a2, nu1, nu2, d = 2) {
  if (h[1] > 0) {
    if (h[2] > 0) {
      complex(real = 0.97)/(1 + complex(imaginary = a1 * sqrt(sum(h^2))))^(nu1 + d/2)/ 
        (1 + complex(imaginary = -a2 * sqrt(sum(h^2))))^(nu2 + d/2)
    } else {
      complex(real = -0.97)/(1 + complex(imaginary = a1 * sqrt(sum(h^2))))^(nu1 + d/2)/ 
        (1 + complex(imaginary = -a2 * sqrt(sum(h^2))))^(nu2 + d/2)
    }
    
  } else {
    if (h[2] > 0) {
      complex(real = -0.97)/(1 + complex(imaginary = a1 * sqrt(sum(h^2))))^(nu1 + d/2)/ 
        (1 + complex(imaginary = -a2 * sqrt(sum(h^2))))^(nu2 + d/2)
    } else {
      complex(real = 0.97)/(1 + complex(imaginary = a1 * sqrt(sum(h^2))))^(nu1 + d/2)/ 
        (1 + complex(imaginary = -a2 * sqrt(sum(h^2))))^(nu2 + d/2)
    }
  }
}

set.seed(22)

sim_f <- simu_bivariate(grid, L = 1000, nu1 = 1.5, nu2 = 1.5, g_new = g_new, cc_fun = h_fun3)
image.plot(matrix(sim_f[,1], nrow = length(x)), axes = F)
image.plot(matrix(sim_f[,2], nrow = length(x)), axes = F)

png('simu_spatial5.png', height = 550, width = 550)
image.plot(matrix(sim_f[,1], nrow = length(x)), axes = F)
dev.off()

png('simu_spatial6.png', height = 550, width = 550)
image.plot(matrix(sim_f[,2], nrow = length(x)), axes = F)
dev.off()



