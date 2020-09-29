# numerically approximate integral
struve_version <- function(h,nu, d) {
  if(h == 0) {
    return(0)
  }
  sign(h)*(abs(h))^nu * pi * 2^(-nu-d/2) / (gamma(nu + d/2) * cos(pi * nu))*
    (besselI(abs(h), nu = nu) - RandomFieldsUtils::struveL(abs(h), -nu))
}
library(fields)

library(RandomFieldsUtils)
struve_H_sign <- function(x, nu){
  pos <- x > 0 
  nu_parity <- (nu %% 2) == 0
  if (pos) {
    return(struveH(x,nu))
  }
  if (!pos & !nu_parity) {
    return(struveH(-x,nu))
  }
  if (!pos & nu_parity) {
    return(-struveH(-x,nu))
  }
}

rbenchmark::benchmark(theirs = {struveH(5, nu = .6)},
                      #  ours = {struveH_mine(5, nu_eval = .6)},
                      combo = {struve_H_sign(5, nu = .6)}, replications = 10000)


### approximate the integral on a grid
whole_fun <- function(u, nu, r, d, h, gamma) {
  sin(abs(gamma) *(u - r *h ))/(u-r*h) * struve_H_sign(x = u, nu =d/2 -1)/(u^(d/2-1)) *
    r^(d-1) * (1 + r^2)^(-nu - d/2)
}

straight_up <- function(x, nu, omega, a) {
  d <- length(x)
  if (sum(omega * a) > 0) {
    return(0)
  }
  sin(sum(x * omega)) *(1 + sum(omega^2))^(-nu - d/2)
}
delta_omega <- .1
omega_seq <- seq(-5, 5, by = delta_omega)
omega_grid <- as.matrix(expand.grid(omega_seq, omega_seq))
x_seq <- seq(-5, 5, by = .12)
x_grid <- as.matrix(expand.grid(x_seq, x_seq))
x_results_test <- sapply(1:nrow(x_grid), function(y) {
  vals <- sapply(1:nrow(omega_grid), function(z) {
    straight_up(x = x_grid[y,],omega = omega_grid[z,],
                 nu = .75, a = c(1,0))
 })
  sum(vals * delta_omega * delta_omega)
})

png('images/nu0_75.png')
image.plot(matrix(x_results_test, nrow = length(x_seq)), main = 'nu = 0.75')
dev.off()

matern <- Matern(apply(x_grid, 1, function(x) sqrt(sum(x^2))), smoothness = 1.25)
png('images/nu1_25_both.png')
image.plot(matrix(matern*.25+ x_results, nrow = length(x_seq)), main = 'nu = 1.25, both symmetric and asymmetic parts')
dev.off()

### approximate using simplified version

full_fun <- function(index, h, gamma, nu, d) {
  whole_fun(u = ur_eval[index,1], r = ur_eval[index,2],
            h = h, gamma = gamma, nu = nu, d = d)
}

delta_r <- .25
max_r <- 5*2
r_seq <- seq(0 + delta_r, max_r, by = delta_r)

delta_u <- .25
max_u <- 5*3
u_seq <- seq( - max_u, max_u, by = delta_u)

# test integration grid
ur_eval <- as.matrix(expand.grid(u = u_seq, r =r_seq))
test <- sapply(1:nrow(ur_eval), function(x) {
  whole_fun(u  = ur_eval[x, 1], nu = .75, r = ur_eval[x,2], d = 2, h = .4214, gamma = .5)
})
test_mat <- matrix(test, nrow = length(u_seq))
image.plot(as.matrix(test_mat))

cov_seq <- seq(-5, 5, by = .12)
cov_grid <- expand.grid(x = cov_seq, y = cov_seq)

d = 2
nu = .75
a <- proc.time(); print(nrow(cov_grid)); print(nrow(ur_eval))
test_x_new <- sapply(1:nrow(cov_grid), function(z) {
  h = sqrt(sum(as.double(cov_grid[z,])^2))
  gamma = cos(atan2(as.double(cov_grid[z,2]), as.double(cov_grid[z,1])))
  fun_adj_vals <- sapply(1:nrow(ur_eval), full_fun, h = h , 
                         gamma = gamma, nu = nu, d = d)
  h <- ifelse(cov_grid[z,2] <= 0, -h, h)
  fun_adj <- sum(fun_adj_vals*delta_r*delta_u, na.rm = T)
  p2 <- fun_adj/pi
  p2
})

### FOR ODD d, Convert from H to J by formula for 1/2 integer orders
#H_{d/2 - 1} = (-1)^(1/2 - d/2) * J_{1-d/2}(z)
b <- proc.time();print(b-a)
p1_vals <- sapply(1:nrow(cov_grid), function(z) {
  h = sqrt(sum(as.double(cov_grid[z,])^2))
  gamma = cos(atan2(as.double(cov_grid[z,2]), as.double(cov_grid[z,1])))
  p1 <- struve_version(h = h, nu = nu, d = d)
})
h_vals <- sapply(1:nrow(cov_grid), function(z) {
  h = sqrt(sum(as.double(cov_grid[z,])^2))
  h <- ifelse(cov_grid[z,2] <= 0, -h, h)
})
image.plot(matrix(p1_vals, nrow = length(cov_seq)))
image.plot(matrix(test_x_new, nrow = length(cov_seq)))
image.plot(t(matrix(ifelse(h_vals < 0, p1_vals-test_x_new,
                         -(p1_vals-test_x_new)) * 2, nrow = length(cov_seq))))

test_struve <- sapply(1:length(cov_seq), function(x) {
  struve_version(h = cov_seq[x], nu =nu, d = d)
  })
plot(cov_seq, test_struve, type = 'l')

############# crap ###############
## in 3 D
delta_omega <- 1
omega_seq <- seq(-5, 5, by = delta_omega)
omega_grid <- as.matrix(expand.grid(omega_seq, omega_seq,omega_seq))

x_seq <- seq(-5, 5, by = .3)
x_grid <- as.matrix(expand.grid(x_seq, x_seq, x_seq))

x_results3 <- sapply(1:nrow(x_grid), function(y) {
  vals <- sapply(1:nrow(omega_grid), function(z) {
    straight_up(x = x_grid[y,],omega = omega_grid[z,],
                nu = .75, a = c(1,0))
  })
  sum(vals * delta_omega * delta_omega)
})

image.plot(matrix(x_results, nrow = length(x_seq)))


gamma_vals1 <- gamma(0:1000 + 3/2)

struveH_mine <- function(z, nu_eval) {
  k_max <- 200
  k <- 0:k_max
  if (nu_eval == 0) {
    sum( (-1)^k *(z/2)^(2*k+nu_eval+ 1) /(gamma_vals1[1:(k_max+1)]^2))
  } else {
    sum( (-1)^k *(z/2)^(2*k+nu_eval+ 1) /(gamma_vals1[1:(k_max+1)] * gamma( k  + nu_eval + 3/2)))
  }
}

struve <- function(z, nu_eval) {
  if (nu_eval == -1/2) {
    return(sqrt(2/(pi * z)) * sinh(z))
  } else if (nu_eval == -3/2) {
    return(sqrt(2/pi) * (z * cosh(z) - sinh(z))/(z^(3/2)))
  }
  k <- 0:200
  (z/2)^(nu_eval + 1) *sum((z/2)^(2*k) /(gamma(k + 3/2) * gamma( k  + nu_eval + 3/2)))
}


# sign(h)*(abs(h))^nu * pi^(3/2) * 2^(-d/2)* gamma(nu + 1/2) / gamma(nu + d/2)/
#   cos(pi * nu)*
#   (besselI(abs(h), nu = nu) - RandomFieldsUtils::struveL(abs(h), -nu))
# sign(h)*(abs(h))^nu * pi^(3/2) * 2^(-nu-1)* gamma(d/2 + 1/2) / gamma(nu + d/2)/
#   cos(pi * nu)*
#   (besselI(abs(h), nu = nu) - RandomFieldsUtils::struveL(abs(h), -nu))


straight_up2 <- function(x, nu, omega, a) {
  d <- length(x)
  if (sum(omega * a) >= 0) {
    return(0)
  }
  h_norm <- sqrt(sum(x^2))
  theta <- acos(sum(x * omega)/h_norm/sqrt(sum(omega^2)))
  r <- sqrt(sum(omega^2))
  sin(r * h_norm *cos(theta)) * sin(theta)^(d-2)*(1 + r^2)^(-nu - d/2)
}

straight_up_from_theta <- function(theta, r, nu, h_norm,d) {
  sin(r * h_norm *cos(theta)) * sin(theta)^(d-2)*(1 + r^2)^(-nu - d/2) *r^(d-1)
}

straight_up_from_theta <- function(theta, r, nu, h_norm,d) {
  sin(r * h_norm *cos(theta)) * sin(theta)^(d-2)*(1 + r^2)^(-nu - d/2) *r^(d-1)
}



a = c(0,1)

delta_r <- .1
delta_theta <- .1
x_results3 <- sapply(1:nrow(x_grid), function(y) {
  h = x_grid[y,]
  h_norm = sqrt(sum(h^2))
  a_norm = sqrt(sum(a^2))
  c = acos(sum(h * a)/h_norm/a_norm)
  if (is.nan(c)) {
    return(0)
  }
  theta_grid <- seq(-pi + c, c, by = delta_theta)
  r_grid <- seq(.01, 5 * sqrt(2), by = delta_r)
  all_grid <- expand.grid(theta_grid, r_grid)
  vals <- sapply(1:nrow(all_grid), function(z) {
    straight_up3_from_theta(theta = all_grid$Var1[z], r = all_grid$Var2[z],
                            nu =.75, h_norm = h_norm, d = 2) #3a = c(0,1))
  })
  sign(-sum(c(a[2], a[1]) * h))*sum(vals * delta_theta *delta_r)
  #print(y)
})
image.plot(matrix(x_results3, nrow = length(x_seq)))
image.plot(matrix(x_results, nrow = length(x_seq)))
