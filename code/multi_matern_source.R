
## Covariance and plotting

matern_cov <- function(s,t, nu, a) {
  if (t-s == 0 & a == 1) {
    return(gamma(nu)/gamma(nu + 1/2) * sqrt(pi))
  } else if(t-s == 0) {
    t <- -s + 10^-18
  }
  (2*pi^(1/2) * abs(t-s)^nu)/( (2*a)^(nu) * gamma(nu + 1/2)) * besselK(a*abs(t-s), nu = nu) 
}
struve_version <- function(s,t,nu, a = 1) {
  if(s-t == 0) {
    return(0)
  }
  sign(t-s)*(abs(t-s)/a)^nu * pi^(1/2) * 2^(-nu)* gamma(-nu + 1/2)*
    (besselI(a * abs(t-s), nu = nu) - struve(a* abs(t-s), -nu))
}
plot_function <- function(s,t,nu) {
  if(s-t == 0) {
    return(0)
  }
  sign(t-s)*(abs(t-s))^nu*
    (besselI(abs(t-s), nu = nu) - struve(abs(t-s), -nu))
}

cross_cov <- function(s,t, nu, z_ij, a) {
  Re(z_ij) * matern_cov(s,t,nu, a) - 
    Im(z_ij) * struve_version(s,t,nu, a)
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
plot_cov <- function(plot_seq =  seq(-5, 5, by = .2), AA_star, nu) {
  params <- par()$mfrow
  cov_val_lag <- sapply(1:length(plot_seq), function(x) {
    c(cross_cov(0, plot_seq[x], nu = nu, z_ij = AA_star[1,1], a = 1),
      cross_cov(0, plot_seq[x], nu = nu, z_ij = AA_star[1,2], a = 1),
      #cross_cov(0, plot_seq[x], nu = nu, z_ij = AA_star[2,1]),
      cross_cov(0, plot_seq[x], nu = nu, z_ij = AA_star[2,2], a = 1)
      )
  })
  par(mfrow = c(2,2))
  plot(plot_seq, cov_val_lag[1,], type = 'l', xlab = 'lag', main = 'Covariance 1')
  plot(plot_seq, cov_val_lag[2,], type = 'l', xlab = 'lag', main = 'Cross Cov')
  plot(-plot_seq, cov_val_lag[2,], type = 'l', xlab = 'lag', main = 'Cross Cov')
  plot(plot_seq, cov_val_lag[3,], type = 'l', xlab = 'lag', main = 'Covariance 2')
  par(mfrow = params)
}

# simulation

H <- function(omega,i,j, AA_star, nu) {
  AA_star_conj <- Conj(AA_star)
  if (omega >=0) {
    H11 <- sqrt( (1 + omega^2)^(-nu -1/2) *  AA_star[1,1])
    if (i == 1 & j == 1) {
      return(H11)
    } else if (i == 2 & j == 1) {
      H21 <- (1 + omega^2)^(-nu -1/2) * AA_star[2,1]/ H11
      return(H21)
    } else if (i == 2 & j == 2) {
      H21 <- (1 + omega^2)^(-nu -1/2) *AA_star[2,1]/ H11
      return(sqrt((1 + omega^2)^(-nu -1/2) * AA_star[2,2] -
                    abs(H21)^2))
    }
  } else {
    H11 <- sqrt( (1 + omega^2)^(-nu -1/2) *  AA_star_conj[1,1])
    if (i == 1 & j == 1) {
      return(H11)
    } else if (i == 2 & j == 1) {
      H21 <- (1 + omega^2)^(-nu -1/2) * AA_star_conj[2,1] /H11
      return(H21)
    } else if (i == 2 & j == 2) {
      H21 <- (1 + omega^2)^(-nu -1/2) * AA_star_conj[2,1] /H11
      return(sqrt((1 + omega^2)^(-nu -1/2) * AA_star_conj[2,2] -
                    abs(H21)^2))
    }
  }
}

gamma_fun <- function(omega, i,j, AA_star, nu) {
  abs(H(omega,i,j,AA_star, nu))/abs(H(omega, j,j,AA_star, nu))
}

theta_fun <- function(omega, i,j, AA_star, nu) {
  atan2(Im(H(omega, i,j,AA_star, nu)),Re(H(omega,i,j,AA_star,nu)))
}

sigmaj2 <- function(j,AA_star,nu) {
  if (j == 1) {
    return(sqrt(pi) * gamma(nu) / gamma(nu + 1/2) * (abs(AA_star)[1,1]))
  } else if (j == 2) {
    val1 <- abs(AA_star[2,2] - abs(AA_star[2,1]^2 / AA_star[1,1])  )
    val2 <- abs(AA_star[2,2] - abs(Conj(AA_star)[2,1]^2 / Conj(AA_star)[1,1])  )
    sqrt(pi) * gamma(nu)/2 / gamma(nu + 1/2) * (val1 + val2)
  }
}


gj <- function(omega, j, AA_star, nu, sigma_vec) {
  abs(H(omega, j,j,AA_star,nu))^2/sigma_vec[j]
}

generate_samples <- function(n_samples, nu) {
  (1/sqrt(2*nu)) * rt(n_samples, df = 2*nu)
}

sim_bivariate <- function(AA_star, nu, t_eval = seq(-5, 5, by = .02), N) {
  phi <- matrix(ncol = 2, nrow = N, runif(n = 2 * N, min = 0, max = 2*pi))
  sigma_vec <- c(sigmaj2(1,AA_star,nu), sigmaj2(2,AA_star,nu))
  
  omega_ij <- matrix(ncol = 2, nrow = N, generate_samples(2*N,nu = nu))
  if (sum(is.na(omega_ij)) > 0) {
    print('increase n_samples, few samples accepted')
    stop()
  }
  f1 <- sqrt(sigma_vec[1] * 2/N) * sapply(t_eval, compute_marginal, omega = omega_ij[,1], phi = phi[,1])
  cross_vec <- sapply(omega_ij[,1], gamma_fun, i = 2, j = 1, AA_star =AA_star, nu =nu)
  theta_eval <- sapply(omega_ij[,1], theta_fun, i = 2, j= 1, AA_star= AA_star, nu = nu)
  f2 <- sqrt(sigma_vec[1] * 2/N) * sapply(t_eval, compute_marginal_cross, omega = omega_ij[,1], phi = phi[,1],
                                          c_vec = cross_vec,
                                          theta = theta_eval)+
    sqrt(sigma_vec[2] * 2/N) * sapply(t_eval, compute_marginal, omega = omega_ij[,2], phi = phi[,2])
  return(data.frame(t_eval, f1, f2))
}

library(Rcpp)
cppFunction('double compute_marginal(double eval_point, NumericVector omega, NumericVector phi) {
              return sum(cos(eval_point * omega + phi));
            }')
cppFunction('double compute_marginal_cross(double eval_point, NumericVector omega, NumericVector phi, NumericVector c_vec, NumericVector theta) {
              return sum(c_vec * cos(eval_point * omega + phi + theta));
            }')

whitt_version <- function(h,nu1, nu2,c11, c12, c2, a1 = 1, a2 = 1) {
  p11 <- c11* (2*sqrt(pi)/(gamma(nu1 + .5)))*((abs(h))^nu1) * (2*a1)^(-nu1) * 
    besselK(x = a1 * abs(h), nu = nu1)
  p22 <- c2* (2*sqrt(pi)/(gamma(nu2 + .5)))*((abs(h))^nu2) * (2*a2)^(-nu2) * 
    besselK(x = a2 * abs(h), nu = nu2)
  if (nu1 == nu2 & a1 == a2) {
#    p12 <- c12 * (abs(h))^nu2 * besselK(x = abs(h), nu = nu2)
    p12 <- c12* (2*sqrt(pi)/(gamma(nu1 + .5)))*((abs( h))^nu1) * (2*a1)^(-nu1) * 
      besselK(x = a1 * abs(h), nu = nu1)
    if (h == 0) {
      p12 <- c12* (2*sqrt(pi)/(gamma(nu1 + .5)))*((abs( 1e-10))^nu1) * (2*a1)^(-nu1) * 
        besselK(x = a1 * abs( 1e-10), nu = nu1)
      p11 <- c11* (2*sqrt(pi)/(gamma(nu1 + .5)))*((abs( 1e-10))^nu1) * (2*a1)^(-nu1) * 
        besselK(x = a1 * abs( 1e-10), nu = nu1)
      p22 <- c2* (2*sqrt(pi)/(gamma(nu2 + .5)))*((abs( 1e-10))^nu2) * (2*a2)^(-nu2) * 
        besselK(x = a2 * abs( 1e-10), nu = nu2)
    }
  } else if (nu1 + nu2 == round(nu1 + nu2)) { 
    p12 <- 0
  } else if(h == 0) {
    p11 <- c11* (2*sqrt(pi)/(gamma(nu1 + .5)))*((abs( 1e-10))^nu1) * (2*a1)^(-nu1) * 
      besselK(x = a1 * abs( 1e-10), nu = nu1)
   # p11 <- c11 * 1e-10^nu1 * (2*a1)^(-nu1) * sqrt(pi)* (2/( gamma(nu1 + .5)))* besselK(x = 1e-10, nu = nu1)
    #p22 <- c2 * 1e-10^nu2 * (2*a2)^(-nu2) *sqrt(pi)* (2/(gamma(nu2 + .5))) * besselK(x = 1e-10, nu = nu2)
    p22 <- c2* (2*sqrt(pi)/(gamma(nu2 + .5)))*((abs( 1e-10))^nu2) * (2*a2)^(-nu2) * 
      besselK(x = a2 * abs( 1e-10), nu = nu2)
    p12 <- 2* pi*c12 * 1/gamma(nu2 + 1/2) * abs(1e-10)^(nu1/2 + nu2/2 - 1/2)* 
      fAsianOptions::whittakerW(x = (a1 + a2)*abs(1e-10),kappa = -nu1/2 + nu2/2,
                                mu = - (nu1+nu2)/2 , ip = 10) *
      (a1 + a2)^(-nu1/2 - nu2/2 - 1/2)*exp((a1 - a2)/2 *1e-10) 
  } else if (h < 0) {
    p12 <- 2* pi *c12 * 1/gamma(nu1 + 1/2) * abs(h)^(nu1/2 + nu2/2 - 1/2)* (a1 + a2)^(-nu1/2 - nu2/2 - 1/2)*
      exp((a1 - a2)/2 * h) * 
      fAsianOptions::whittakerW(x = (a1 + a2)*abs(h),  kappa = nu1/2 - nu2/2, mu = - (nu1+nu2)/2 , ip = 10)
  } else {
    p12 <- 2* pi *c12 * 1/gamma(nu2 + 1/2) * abs(h)^(nu1/2 + nu2/2 - 1/2)* (a1 + a2)^(-nu1/2 - nu2/2 - 1/2)*
      exp((a1 - a2)/2 * h) *  
      fAsianOptions::whittakerW(x = (a1 + a2)*abs(h), kappa = -nu1/2 + nu2/2,mu = - (nu1+nu2)/2 , ip = 10)
  }
  return(c(p11, Re(p12),p22))
}

whittaker_covariance_matrix <- function(locs, nu1, nu2, c11, c12, c2, a1,a2) {
  grid <- expand.grid(locs, locs)
  cov_val <- sapply(1:nrow(grid), function(x) {
    whitt_version(grid[['Var1']][x]-grid[['Var2']][x], nu1= nu1, nu2 = nu2, c11 = c11, c12 = c12, c2 = c2)
  })
  cov_mat1 <- matrix(cov_val[1,], nrow = length(locs))
  cov_mat2 <- matrix(cov_val[3,], nrow = length(locs))
  cov_mat12 <- matrix(cov_val[2,], nrow = length(locs))
  cov_mat_all <- rbind(cbind(cov_mat1, cov_mat12),
                       cbind(t(cov_mat12), cov_mat2))
}
whittaker_covariance_matrix_lags <- function(locs, nu1, nu2, c11, c12, c2, a1,a2) {
  lags1 <- locs - locs[1]
  lags2 <- locs - locs[length(locs)]
  cov_val1 <- sapply(1:length(lags1), function(x) {
    whitt_version(lags1[x], nu1= nu1, nu2 = nu2, c11 = c11, c12 = c12, c2 = c2)
  })
  cov_val2 <- sapply(1:length(lags2), function(x) {
    whitt_version(lags2[x], nu1= nu1, nu2 = nu2, c11 = c11, c12 = c12, c2 = c2)
  })
  cov_mat1 <- toeplitz(cov_val1[1,])
  cov_mat2 <- toeplitz(cov_val1[3,])
  cov_mat12 <- toeplitz(rev(cov_val2[2,]))
  lt <- lower.tri(cov_mat12)
  cov_mat12[lt] <- toeplitz(cov_val1[2,])[lt]
  cov_mat_all <- rbind(cbind(cov_mat1, cov_mat12),
                       cbind(t(cov_mat12), cov_mat2))
}

imaginary_covariance_matrix_lags <- function(locs, nu1, nu2, c11, c12, c22, a1, a2) {
  lags1 <- locs - locs[1]
  lags2 <- locs - locs[length(locs)]
  
  cov_val11 <- sapply(1:length(lags1), function(x) {
    cross_cov(s = 0, t = lags1[x], nu = nu1, z_ij = c11, a  = a1)
  })
  cov_val22 <- sapply(1:length(lags1), function(x) {
    cross_cov(s = 0, t = lags1[x], nu = nu1, z_ij = c22, a  = a1)
  })
  cross_cov1 <- sapply(1:length(lags1), function(x) {
    cross_cov(s = 0, t = lags1[x], nu = nu1, z_ij = c12, a  = a1)
  })
  cross_cov2 <- sapply(1:length(lags1), function(x) {
    cross_cov(s = 0, t = lags2[x], nu = nu1, z_ij = c12, a  = a1)
  })
  cov_mat1 <- toeplitz(cov_val11)
  cov_mat2 <- toeplitz(cov_val22)
  cov_mat12 <- toeplitz(rev(cross_cov2))
  lt <- lower.tri(cov_mat12)
  cov_mat12[lt] <- toeplitz(cross_cov1)[lt]
  cov_mat_all <- rbind(cbind(cov_mat1, cov_mat12),
                       cbind(t(cov_mat12), cov_mat2))
}

simulate_manually <- function(cholesky, n_simu, locs) {
  simulated <- t(cholesky) %*% matrix(nrow = nrow(cholesky),
                                          ncol = n_simu,
                                          rnorm(nrow(cholesky) *n_simu))
  s1 <- simulated[1:length(locs),]
  s2 <- simulated[-(1:length(locs)),]
  simulation <- data.frame(var1 = as.double(s1), var2 = as.double(s2),
                           t = locs, 
                           simulation = rep(1:n_simu, each = length(locs)))
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


full_cross_cov_single <- function(h, nu, a, realp, imp, norm_type = 'D') {
  -imp * plot_function( h, nu= nu, a = a)* 
    2^(-nu) * pi^(3/2)/cos(nu * pi)/gamma(nu + .5)*
    norm_constant(nu_1 = nu, nu_2 = nu, a_1 = a, a_2 = a, norm_type = norm_type)  +
    realp *  whitt_version( h, nu1 = nu, nu2 = nu,c2 = 1,c11 = 1, c12 = 1, a1 = a, a2 = a)[2]*
                      norm_constant(nu_1 = nu, nu_2 = nu, a_1 = a, a_2 = a, norm_type = norm_type)
}

whitt_only_single <- function(h, nu1, nu2, a1, a2, realp, imp, norm_type = 'D', which_val = 2) {
  realp *  whitt_version( h, nu1 = nu1, nu2 = nu2,c2 = 1,c11 = 1, c12 = 1, a1 = a1, a2 = a2)[which_val]*
    norm_constant(nu_1 = nu1, nu_2 = nu2, a_1 = a1, a_2 = a2, norm_type = norm_type)
}

plot_function <- function(h,nu, a = 1) {
  if(h == 0) {
    return(0)
  }
  sign(h) * (abs(h)/a)^nu*
    (besselI(a*abs(h), nu = nu) - struve(a*abs(h), -nu))
}

