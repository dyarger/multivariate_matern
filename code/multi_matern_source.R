
## Covariance and plotting

matern_cov <- function(s,t, nu) {
  if (t-s == 0) {
    return(gamma(nu)/gamma(nu + 1/2) * sqrt(pi))
  }
  (pi^(1/2) * abs(t-s)^nu)/(gamma(nu + 1/2) * 2^(nu-1)) * besselK(abs(t-s), nu = nu) 
}
struve_version <- function(s,t,nu) {
  if(s-t == 0) {
    return(0)
  }
  sign(t-s)*(abs(t-s))^nu * pi^(1/2) * 2^(nu - 1)* gamma(-nu + 1/2)*
    (besselI(abs(t-s), nu = nu) - struve(abs(t-s), -nu))
}
plot_function <- function(s,t,nu) {
  if(s-t == 0) {
    return(0)
  }
  #(t-s)^nu*
    sign(t-s)*(abs(t-s))^nu*
    (besselI(abs(t-s), nu = nu) - struve(abs(t-s), -nu))
}

cross_cov <- function(s,t, nu, z_ij) {
  # if (nu == 1/2 | nu == 3/2) {
  #   return(Re(z_ij) * matern_cov(s,t,nu))
  # }
  Re(z_ij) * matern_cov(s,t,nu) - 
    Im(z_ij) * struve_version(s,t,nu)
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
    c(cross_cov(0, plot_seq[x], nu = nu, z_ij = AA_star[1,1]),
      cross_cov(0, plot_seq[x], nu = nu, z_ij = AA_star[1,2]),
      #cross_cov(0, plot_seq[x], nu = nu, z_ij = AA_star[2,1]),
      cross_cov(0, plot_seq[x], nu = nu, z_ij = AA_star[2,2])
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

# generate_samples1 <- function(n_samples, j, AA_star, nu, sigma_vec, constant) {
#   n <- rt(n_samples, df = nu)
#   u <- runif(n_samples)
#   vals <- sapply(n, gj, j= j , AA_star = AA_star, nu =nu, sigma_vec = sigma_vec)
#   n[ u < vals/ ( constant * dt(n, df = nu)) ]
# }

generate_samples <- function(n_samples, nu) {
  (1/sqrt(2*nu)) * rt(n_samples, df = 2*nu)
}

#generate_samples <- generate_samples1

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

#library(Rcpp)
# cppFunction('double compute_marginal(double eval_point, NumericVector omega, NumericVector phi) {
#               return sum(cos(eval_point * omega + phi));
#             }')
# cppFunction('double compute_marginal_cross(double eval_point, NumericVector omega, NumericVector phi, NumericVector c_vec, NumericVector theta) {
#               return sum(c_vec * cos(eval_point * omega + phi + theta));
#             }')
# cppFunction('double f1_eval(double eval_point, NumericVector omega, NumericVector phi) {
#               return sqrt();
#             }')
# eval_point = 1; omega <- c(0,0,0); phi = c(1,2,1)
# compute_marginal(eval_point = eval_point, omega = omega, phi = phi)
# compute_marginal_cross(eval_point = eval_point, omega = omega, phi = phi, c_vec = c(1,1,1))
# sum(cos(eval_point * omega + phi ))
# f2 <- vector(mode = 'double', length = length(t_eval))
# for (q in 1:length(t_eval)) {
#   f2[q] <- sqrt(sigma_vec[1] * 2/N) *sum(
#     cross_vec *
#       cos(t_eval[q] *omega_ij[,1] + phi[,1]+theta_eval))+
#     sqrt(sigma_vec[2] * 2/N) *sum(
#       cos(t_eval[q] *omega_ij[,2] + phi[,2]))
# }
# f1 <- vector(mode = 'double', length = length(t_eval))
# for (q in 1:length(t_eval)) {
#   f1[q] <- sqrt(sigma_vec[1] * 2/N) *
#     sum(cos(t_eval[q] *omega_ij[,1] + phi[,1]))
# }