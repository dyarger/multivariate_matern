create_distance_matrices <- function(loc, pred_loc) {
  # create long/lat distance matrices
  
  dist <- fields::rdist.earth(x1 = loc, x2 = loc, miles = F)
  diag(dist) <- 0

  dist_tens <- array(NA, c(dim(dist), 2))
  reference_location <- cbind(.5, .5)
  for (i in 1:dim(dist_tens)[1]) {
    for (j in 1:dim(dist_tens)[1]) {
      dist_long_prelim <- fields::rdist.vec(cbind(loc[i,1],
                                                        reference_location[,2]), 
                                                  cbind(loc[j,1],
                                                        reference_location[,2]))
      if (loc[i,1] < loc[j,1]) {
        dist_long_prelim <- -dist_long_prelim
      }
      
      dist_lat_prelim <- fields::rdist.vec(cbind(reference_location[,1],
                                                       loc[i,2]), 
                                                 cbind(reference_location[,1],
                                                       loc[j,2]))
      if (loc[i,2] < loc[j,2]) {
        dist_lat_prelim <- -dist_lat_prelim
      }
      
      dist_tens[i,j,1] <- dist_long_prelim
      dist_tens[i,j,2] <- dist_lat_prelim
    }
  }
  dist_tens
}

fft_2d_alt <- function(grid_info, nu1 = .5, nu2 = .5, a1 = 1, a2 = 1, 
         Psi, 
         Mu, d = 2) {
  freq_grid <- grid_info[['freq_grid']]
  freq_points <- grid_info[['freq_points']]
  delta_t <- grid_info[['delta_t']]
  n_points <- grid_info[['n_points']]
  x_vals <- grid_info[['x_vals']]
  x_max <- grid_info[['x_max']]
  x_vals_eg <- grid_info[['x_vals_eg']]
  phase_factor_mat <- grid_info[['phase_factor_mat']]
  # https://stackoverflow.com/questions/24077913/discretized-continuous-fourier-transform-with-numpy
  tv <- ( a1^2 + freq_grid[['r']]^2)^(-nu1/2 - d/4) *
    ( a2^2 + freq_grid[['r']]^2)^(-nu2/2 - d/4) * 
    Mu(freq_grid[['theta']])
  tv_mat <- matrix(nrow = length(freq_points), ncol = length(freq_points), tv)
  ff_res <- fftwtools::fftw_c2c_2d(data = tv_mat, inverse = 1)/n_points^2
  p <- ncol(ff_res)/2
  ff_res_adj <- cbind(rbind(ff_res[(p + 1):(2*p),(p + 1):(2*p)], ff_res[1:p,(p + 1):(2*p)]),
                      rbind(ff_res[(p + 1):(2*p),1:p], ff_res[1:p,1:p])) * -phase_factor_mat
  x_vals <- x_vals - (x_vals[length(x_vals)] - x_vals[1]) / 2
  cbind(x_vals_eg, 'val' = (length(x_vals))^(2) *
          as.double(Re(ff_res_adj) * norm_constant(nu1, nu2, a1, a2, d = d)) * 2 / pi / x_max^2 * pi^4)
}


construct_matrix_alt <- function(nu1, nu2, a1, a2, 
                             grid_info,
                             Psi, Mu, dist_vals) {
  C_val <- fft_2d_alt(grid_info = grid_info, nu1 = nu1, nu2 = nu2, 
                  a1 = a1, a2 = a2,  Psi = Psi, Mu = Mu)
  
  tmat <- matrix(C_val[['val']], nrow = sqrt(nrow(C_val)))
  long1 <- grid_info[['x_vals']]
  lat1 <- grid_info[['x_vals']]
  
  interp_points <- fields::interp.surface(obj = list(x = long1, y = lat1, z = tmat),
                                          loc = dist_vals)
  matrix(interp_points, nrow = sqrt(length(dist_vals)/2), ncol = sqrt(length(dist_vals)/2))
}

construct_matrix_orig <- function(nu1, nu2, a1, a2, 
                                 grid_info,
                                 Psi, Mu, dist_vals) {
  C_val <- fft_2d(grid_info = grid_info, nu1 = nu1, nu2 = nu2, 
                      a1 = a1, a2 = a2,  Psi = Psi, Mu = Mu)
  
  tmat <- matrix(C_val[['val']], nrow = sqrt(nrow(C_val)))
  long1 <- grid_info[['x_vals']]
  lat1 <- grid_info[['x_vals']]
  
  interp_points <- fields::interp.surface(obj = list(x = long1, y = lat1, z = tmat),
                                          loc = dist_vals)
  matrix(interp_points, nrow = sqrt(length(dist_vals)/2), ncol = sqrt(length(dist_vals)/2))
}
construct_bivariate_matrix_mm <- function(nu1, nu2, nu12, a1, a2, a12,
                                          grid_info,
                                          Sigma11, Sigma22,Psi, Mu, nugget1, nugget2,
                                          dist_vals) {
  Mu1 <- function(theta_x, theta_y, entry_1, entry_2) {
    Sigma11 
  }
  Mu2 <- function(theta_x, theta_y, entry_1, entry_2) {
    Sigma22 
  }
  
  C1 <- construct_matrix_alt(nu1 = nu1, nu2 = nu1, a1 = a1, a2 = a1,
                             Psi = Psi, Mu = Mu1, grid_info = grid_info, dist_vals = dist_vals)
  diag(C1) <- diag(C1) + nugget1
  C12 <- construct_matrix_alt(nu1 = nu12, nu2 = nu12, a1 = a12, a2 = a12,
                              Psi = Psi, Mu = Mu, grid_info = grid_info, dist_vals = dist_vals)
  C2 <- construct_matrix_alt(nu1 = nu2, nu2 = nu2, a1 = a2, a2 = a2,
                             Psi = Psi, Mu = Mu2, grid_info = grid_info, 
                             dist_vals = dist_vals)
  diag(C2) <- diag(C2) + nugget2
  rbind(cbind(C1, C12), cbind(t(C12), C2))
}
construct_bivariate_matrix_alt <- function(nu1, nu2, a1, a2,
                                          grid_info,
                                          Sigma11, Sigma22,Psi, Mu, nugget1, nugget2,
                                          dist_vals) {
  Mu1 <- function(theta_x, theta_y, entry_1, entry_2) {
    Sigma11 
  }
  Mu2 <- function(theta_x, theta_y, entry_1, entry_2) {
    Sigma22
  }
  
  C1 <- construct_matrix_alt(nu1 = nu1, nu2 = nu1, a1 = a1, a2 = a1,
                             Psi = Psi, Mu = Mu1, grid_info = grid_info, dist_vals = dist_vals)
  diag(C1) <- diag(C1) + nugget1
  C12 <- construct_matrix_alt(nu1 = nu1, nu2 = nu2, a1 = a1, a2 = a2,
                              Psi = Psi, Mu = Mu, grid_info = grid_info, dist_vals = dist_vals)
  C2 <- construct_matrix_alt(nu1 = nu2, nu2 = nu2, a1 = a2, a2 = a2,
                             Psi = Psi, Mu = Mu2, grid_info = grid_info, 
                             dist_vals = dist_vals)
  diag(C2) <- diag(C2) + nugget2
  rbind(cbind(C1, C12), cbind(t(C12), C2))
}

construct_bivariate_matrix_orig <- function(nu1, nu2, a1, a2,
                                           grid_info,
                                           Sigma11, Sigma22,Psi, Mu, nugget1, nugget2,
                                           dist_vals) {
  Mu1 <- function(theta_x, theta_y, entry_1, entry_2) {
    Sigma11 
  }
  Mu2 <- function(theta_x, theta_y, entry_1, entry_2) {
    Sigma22 
  }
  
  C1 <- construct_matrix_orig(nu1 = nu1, nu2 = nu1, a1 = a1, a2 = a1,
                             Psi = Psi, Mu = Mu1, grid_info = grid_info, dist_vals = dist_vals)
  diag(C1) <- diag(C1) + nugget1
  C12 <- construct_matrix_orig(nu1 = nu1, nu2 = nu2, a1 = a1, a2 = a2,
                              Psi = Psi, Mu = Mu, grid_info = grid_info, dist_vals = dist_vals)
  C2 <- construct_matrix_orig(nu1 = nu2, nu2 = nu2, a1 = a2, a2 = a2,
                             Psi = Psi, Mu = Mu2, grid_info = grid_info, 
                             dist_vals = dist_vals)
  diag(C2) <- diag(C2) + nugget2
  rbind(cbind(C1, C12), cbind(t(C12), C2))
}
# log likelihood as function of parameters
ll_fun <- function(par, dist_vals, response, grid_info, n1, n2) {
  print(exp(par))
  nu1 <- exp(par[1]); nu2 <- exp(par[2])
  a1 <- exp(par[3]); a2 <- exp(par[4])
  nugget1 <- exp(par[5]); nugget2 <- exp(par[6])
  Sigma11 <- exp(par[7]); Sigma22 <- exp(par[8])
  Sigma12_re <- par[9]*sqrt(Sigma11)*sqrt(Sigma22)
  Sigma12_im <- par[10]*sqrt(Sigma11)*sqrt(Sigma22)
  if (sum(abs(par[9:10])) > 1) {
    return(80000)
  }
  Psi <- function(x) {
    sign(x)
  }
  Mu <- function(x) complex(real = Sigma12_re, imaginary = Sigma12_im * sign(x))
  cov_mat <- construct_bivariate_matrix_orig(nu1 = nu1, nu2 = nu2, 
                                             a1 = a1, a2 = a2, 
                                             Sigma11 = Sigma11, Sigma22 =Sigma22, 
                                             Psi = Psi, Mu = Mu, dist_vals = dist_vals,
                                             grid_info = grid_info,
                                             nugget1 = nugget1, nugget2 = nugget2)
  chol_mat <- base::chol(cov_mat)
  ll_val <- -nrow(cov_mat)/2 * log(2 * pi) - 1/2 * sum(backsolve(chol_mat, response, transpose = T)^2) - 
    sum(log(diag(chol_mat)))
  print(ll_val)
  - ll_val
}

# log likelihood as function of parameters
ll_fun_re <- function(par, dist_vals, response, grid_info, n1, n2) {
  print(exp(par))
  nu1 <- exp(par[1]); nu2 <- exp(par[2])
  a1 <- exp(par[3]); a2 <- exp(par[4])
  nugget1 <- exp(par[5]); nugget2 <- exp(par[6])
  Sigma11 <- exp(par[7]); Sigma22 <- exp(par[8])
  Sigma12 <- par[9]*sqrt(Sigma11)*sqrt(Sigma22)
  Psi <- function(x) {
    sign(x)
  }
  Mu <- function(x) Sigma12
  cov_mat <- construct_bivariate_matrix_orig(nu1 = nu1, nu2 = nu2, 
                                             a1 = a1, a2 = a2, 
                                             Sigma11 = Sigma11, Sigma22 =Sigma22, 
                                             Psi = Psi, Mu = Mu, dist_vals = dist_vals,
                                             grid_info = grid_info,
                                             nugget1 = nugget1, nugget2 = nugget2)
  chol_mat <- base::chol(cov_mat)
  ll_val <- -nrow(cov_mat)/2 * log(2 * pi) - 1/2 * sum(backsolve(chol_mat, response, transpose = T)^2) - 
    sum(log(diag(chol_mat)))
  print(ll_val)
  - ll_val
}


ll_fun_psi_real <- function(par, dist_vals,
                            response, grid_info, n1, n2) {
  print(exp(par))
  nu1 <- exp(par[1]); nu2 <- exp(par[2])
  a1 <- exp(par[3]); a2 <- exp(par[4])
  nugget1 <- exp(par[5]); nugget2 <- exp(par[6])
  Sigma11 <- exp(par[7]); Sigma22 <- exp(par[8])
  Sigma12 <- par[9]*sqrt(Sigma11)*sqrt(Sigma22)
  theta_star <- atan2(par[11], par[10])
  Psi <- function(theta) {
    if (theta_star < 0) {
      ifelse(theta > theta_star & theta < theta_star + pi, 1, -1)
    } else {
      ifelse(theta > theta_star | theta < theta_star - pi, 1, -1)
    }
  }
  Mu <- function(x) Sigma12
  cov_mat <- construct_bivariate_matrix_orig(nu1 = nu1, nu2 = nu2, 
                                             a1 = a1, a2 = a2, 
                                             Sigma11 = Sigma11, Sigma22 =Sigma22, 
                                             Psi = Psi, Mu = Mu, dist_vals = dist_vals,
                                             grid_info = grid_info,
                                             nugget1 = nugget1, nugget2 = nugget2)
  chol_mat <- base::chol(cov_mat)
  ll_val <- -nrow(cov_mat)/2 * log(2 * pi) - 1/2 * sum(backsolve(chol_mat, response, transpose = T)^2) - 
    sum(log(diag(chol_mat)))
  print(ll_val)
  - ll_val
}


# with an imaginary entry
ll_fun_psi_im <- function(par, response, grid_info,
                          dist_vals, n1, n2) {
  print(exp(par[1:8]))
  print((par[9:10]))
  nu1 <- exp(par[1]); nu2 <- exp(par[2])
  a1 <- exp(par[3]); a2 <- exp(par[4])
  nugget1 <- exp(par[5]); nugget2 <- exp(par[6])
  Sigma11 <- exp(par[7]); Sigma22 <- exp(par[8])
  Sigma12 <- par[9]*sqrt(Sigma11)*sqrt(Sigma22)
  Sigma12im <- par[10]*sqrt(Sigma11)*sqrt(Sigma22)
  if (sum(abs(par[9:10])) > 1) {
    return(80000)
  }
  theta_star <- atan2(par[12], par[11])
  Psi_fun <- function(theta) {
    if (theta_star < 0) {
      ifelse(theta > theta_star & theta < theta_star + pi, 1, -1)
    } else {
      ifelse(theta > theta_star | theta < theta_star - pi, 1, -1)
    }
  }
  theta_star2 <- atan2(par[14], par[13])
  Psi_fun2 <- function(theta) {
    if (theta_star < 0) {
      ifelse(theta > theta_star2 & theta < theta_star2 + pi, 1, -1)
    } else {
      ifelse(theta > theta_star2 | theta < theta_star2 - pi, 1, -1)
    }
  }
  Mu <- function(x) {complex(real = Sigma12,
                   imaginary = Sigma12im*Psi_fun2(x))}
  cov_mat <- construct_bivariate_matrix_orig(nu1 = nu1, nu2 = nu2, a1 = a1, a2 = a2, 
                                        Mu = Mu,
                                        Psi = Psi_fun, Sigma11 = Sigma11, 
                                        Sigma22 = Sigma22,
                                        grid_info = grid_info,
                                        nugget1 = nugget1, nugget2 = nugget2,
                                        dist_vals = dist_vals)
  chol_mat <- base::chol(cov_mat)
  ll_val <- -nrow(cov_mat)/2 * log(2 * pi) - 1/2 * sum(backsolve(chol_mat, response, transpose = T)^2) - 
    sum(log(diag(chol_mat)))
  print(ll_val)
  - ll_val
}


ll_fun_alt <- function(par, dist_one, grid_info, response) {
  if (sum(abs(par[9:10])) > 1) {
    return(8000)
  }
  
  nu1 <- exp(par[1])
  nu2 <- exp(par[2])
  a1 <- exp(par[3])
  a2 <- exp(par[4])
  nugget1 <-  exp(par[5])
  nugget2 <-  exp(par[6])
  Sigma11 <- exp(par[7])
  Sigma22 <- exp(par[8])
  Sigma12re <- par[9]*sqrt(Sigma11)*sqrt(Sigma22)
  Sigma12im <- par[10]*sqrt(Sigma11)*sqrt(Sigma22)
  cov_mat <- construct_bivariate_matrix_alt(nu1 = nu1, nu2 = nu2, a1 = a1, a2 = a2, 
                                            Sigma11 = Sigma11, Sigma22 = Sigma22, 
                                            Sigma12 = complex(real = Sigma12re, imaginary = Sigma12im),
                                            dist_mat = dist_one, grid_info = grid_info,
                                            nugget1 = nugget1, nugget2 = nugget2)
  inv_all <- solve(cov_mat)
  ll_val <- (-nrow(cov_mat)/2 * log(2 * pi) - 1/2 * as.double(t(response) %*% inv_all %*% response) - 
               1/2 * as.vector(determinant(cov_mat)$modulus))
  - ll_val
}

ll_fun_alt_re <- function(par, dist_one, grid_info, response) {
  if (sum(abs(par[9])) > 1) {
    return(8000)
  }
  
  nu1 <- exp(par[1])
  nu2 <- exp(par[2])
  a1 <- exp(par[3])
  a2 <- exp(par[4])
  nugget1 <-  exp(par[5])
  nugget2 <-  exp(par[6])
  Sigma11 <- exp(par[7])
  Sigma22 <- exp(par[8])
  Sigma12re <- par[9]*sqrt(Sigma11)*sqrt(Sigma22)
  cov_mat <- construct_bivariate_matrix_alt(nu1 = nu1, nu2 = nu2, a1 = a1, a2 = a2, 
                                            Sigma11 = Sigma11, Sigma22 = Sigma22, 
                                            Sigma12 = complex(real = Sigma12re, imaginary = 0),
                                            dist_mat = dist_one, grid_info = grid_info,
                                            nugget1 = nugget1, nugget2 = nugget2)
  inv_all <- solve(cov_mat)
  ll_val <- (-nrow(cov_mat)/2 * log(2 * pi) - 1/2 * as.double(t(response) %*% inv_all %*% response) - 
               1/2 * as.vector(determinant(cov_mat)$modulus))
  - ll_val
}

ll_fun_mm <- function(par, dist_one, grid_info, response) {
  if (sum(abs(par[11:12])) > 1) {
    return(8000)
  }
  print(exp(par))
  #print(par)
  
  nu1 <- exp(par[1])
  nu2 <- exp(par[2])
  nu12 <- exp(par[3])
  a1 <- exp(par[4])
  a2 <- exp(par[5])
  a12 <- exp(par[6])
  nugget1 <- exp(par[7])
  nugget2 <- exp(par[8])
  Sigma11 <- exp(par[9])
  Sigma22 <- exp(par[10])
  Sigma12re <- par[11]*sqrt(Sigma11)*sqrt(Sigma22)
  Sigma12im <- par[12]*sqrt(Sigma11)*sqrt(Sigma22)
  cov_mat <- construct_bivariate_matrix_mm(nu1 = nu1, nu2 = nu2, nu12 = nu12, a1 = a1, a2 = a2, 
                                           a12 = a12,
                                           Sigma11 = Sigma11, Sigma22 = Sigma22, 
                                           Sigma12 = complex(real = Sigma12re, imaginary = Sigma12im),
                                           dist_mat = dist_one, grid_info = grid_info,
                                           nugget1 = nugget1, nugget2 = nugget2)
  if ((sum(is.na(cov_mat)) > 0) | (sum(cov_mat == Inf) > 0) | (sum(cov_mat == -Inf) > 0) |
      (sum(is.nan(cov_mat)) > 0)) {
    print('issue mat')
    return(8000)
  }
  eig_values <- Re(eigen(cov_mat)$values)
  if (sum(is.na(eig_values), na.rm = T) > 0) {
    print('issue eig')
    return(8000)
  }
  if (min(eig_values, na.rm = T) < min(c(nugget1, nugget2))) {
    print('issue eig')
    return(8000)
  }
  inv_all <- tryCatch(expr = {solve(cov_mat)}, error = function(e) {return(NULL)})
  if (is.null(inv_all)) {
    print('issue')
    return(8000)
  }
  ll_val <- (-nrow(cov_mat)/2 * log(2 * pi) - 1/2 * as.double(t(response) %*% inv_all %*% response) - 
               1/2 * as.vector(determinant(cov_mat)$modulus))
  if (is.na(ll_val)) {
    print('issue2')
    return(8000)
  }
  if (abs(ll_val) == Inf) {
    print('issue3')
    return(8000)
  }
  print(ll_val)
  - ll_val
}


ll_fun_mm_re <- function(par, dist_one, grid_info, response) {
  if (sum(abs(par[11])) > 1) {
    return(8000)
  }
  print(exp(par))
  
  nu1 <- exp(par[1])
  nu2 <- exp(par[2])
  nu12 <- exp(par[3])
  a1 <- exp(par[4])
  a2 <- exp(par[5])
  a12 <- exp(par[6])
  nugget1 <- exp(par[7])
  nugget2 <- exp(par[8])
  Sigma11 <- exp(par[9])
  Sigma22 <- exp(par[10])
  Sigma12re <- par[11]*sqrt(Sigma11)*sqrt(Sigma22)
  cov_mat <- construct_bivariate_matrix_mm(nu1 = nu1, nu2 = nu2, nu12 = nu12, a1 = a1, a2 = a2, 
                                           a12 = a12,
                                           Sigma11 = Sigma11, Sigma22 = Sigma22, 
                                           Sigma12 = complex(real = Sigma12re, imaginary = 0),
                                           dist_mat = dist_one, grid_info = grid_info,
                                           nugget1 = nugget1, nugget2 = nugget2)
  if ((sum(is.na(cov_mat)) > 0) | (sum(cov_mat == Inf) > 0) | (sum(cov_mat == -Inf) > 0) |
      (sum(is.nan(cov_mat)) > 0)) {
    #print('issue mat')
    return(8000)
  }
  eig_values <- Re(eigen(cov_mat)$values)
  if (sum(is.na(eig_values), na.rm = T) > 0) {
    #print('issue eig')
    return(8000)
  }
  if (min(eig_values, na.rm = T) < min(c(nugget1, nugget2))) {
    #print('issue eig')
    return(8000)
  }
  inv_all <- tryCatch(expr = {solve(cov_mat)}, error = function(e) {return(NULL)})
  if (is.null(inv_all)) {
    #print('issue')
    return(8000)
  }
  ll_val <- (-nrow(cov_mat)/2 * log(2 * pi) - 1/2 * as.double(t(response) %*% inv_all %*% response) - 
               1/2 * as.vector(determinant(cov_mat)$modulus))
  if (is.na(ll_val)) {
    #print('issue2')
    return(8000)
  }
  if (abs(ll_val) == Inf) {
    #print('issue3')
    return(8000)
  }
  print(ll_val)
  - ll_val
}


ll_fun_mm_pars <- function(par, dist_one, grid_info, response) {
  if (sum(abs(par[8:9])) > 1) {
    return(8000)
  }
  print(exp(par))
  #print(par)
  
  nu1 <- exp(par[1])
  nu2 <- exp(par[2])
  nu12 <- (nu1 + nu2)/2
  a1 <- a2 <- a12 <- exp(par[3])
  nugget1 <- exp(par[4])
  nugget2 <- exp(par[5])
  Sigma11 <- exp(par[6])
  Sigma22 <- exp(par[7])
  Sigma12re <- par[8]*sqrt(Sigma11)*sqrt(Sigma22)
  Sigma12im <- par[9]*sqrt(Sigma11)*sqrt(Sigma22)
  cov_mat <- construct_bivariate_matrix_mm(nu1 = nu1, nu2 = nu2, nu12 = nu12, a1 = a1, a2 = a2, 
                                           a12 = a12,
                                           Sigma11 = Sigma11, Sigma22 = Sigma22, 
                                           Sigma12 = complex(real = Sigma12re, imaginary = Sigma12im),
                                           dist_mat = dist_one, grid_info = grid_info,
                                           nugget1 = nugget1, nugget2 = nugget2)
  if ((sum(is.na(cov_mat)) > 0) | (sum(cov_mat == Inf) > 0) | (sum(cov_mat == -Inf) > 0) |
      (sum(is.nan(cov_mat)) > 0)) {
    #print('issue mat')
    return(8000)
  }
  eig_values <- Re(eigen(cov_mat)$values)
  if (sum(is.na(eig_values), na.rm = T) > 0) {
    #print('issue eig')
    return(8000)
  }
  if (min(eig_values, na.rm = T) < min(c(nugget1, nugget2))) {
    #print('issue eig')
    return(8000)
  }
  inv_all <- tryCatch(expr = {solve(cov_mat)}, error = function(e) {return(NULL)})
  if (is.null(inv_all)) {
    #print('issue')
    return(8000)
  }
  ll_val <- (-nrow(cov_mat)/2 * log(2 * pi) - 1/2 * as.double(t(response) %*% inv_all %*% response) - 
               1/2 * as.vector(determinant(cov_mat)$modulus))
  if (is.na(ll_val)) {
    #print('issue2')
    return(8000)
  }
  if (abs(ll_val) == Inf) {
    #print('issue3')
    return(8000)
  }
  print(ll_val)
  - ll_val
}


ll_fun_mm_re_pars <- function(par, dist_one, grid_info, response) {
  if (sum(abs(par[8])) > 1) {
    return(8000)
  }
  print(exp(par))
  
  nu1 <- exp(par[1])
  nu2 <- exp(par[2])
  nu12 <- (nu1 + nu2)/2
  a1 <- a2 <- a12 <- exp(par[3])
  nugget1 <- exp(par[4])
  nugget2 <- exp(par[5])
  Sigma11 <- exp(par[6])
  Sigma22 <- exp(par[7])
  Sigma12re <- par[8]*sqrt(Sigma11)*sqrt(Sigma22)
  cov_mat <- construct_bivariate_matrix_mm(nu1 = nu1, nu2 = nu2, nu12 = nu12, a1 = a1, a2 = a2, 
                                           a12 = a12,
                                           Sigma11 = Sigma11, Sigma22 = Sigma22, 
                                           Sigma12 = complex(real = Sigma12re, imaginary = 0),
                                           dist_mat = dist_one, grid_info = grid_info,
                                           nugget1 = nugget1, nugget2 = nugget2)
  if ((sum(is.na(cov_mat)) > 0) | (sum(cov_mat == Inf) > 0) | (sum(cov_mat == -Inf) > 0) |
      (sum(is.nan(cov_mat)) > 0)) {
    #print('issue mat')
    return(8000)
  }
  eig_values <- Re(eigen(cov_mat)$values)
  if (sum(is.na(eig_values), na.rm = T) > 0) {
    #print('issue eig')
    return(8000)
  }
  if (min(eig_values, na.rm = T) < min(c(nugget1, nugget2))) {
    #print('issue eig')
    return(8000)
  }
  inv_all <- tryCatch(expr = {solve(cov_mat)}, error = function(e) {return(NULL)})
  if (is.null(inv_all)) {
    #print('issue')
    return(8000)
  }
  ll_val <- (-nrow(cov_mat)/2 * log(2 * pi) - 1/2 * as.double(t(response) %*% inv_all %*% response) - 
               1/2 * as.vector(determinant(cov_mat)$modulus))
  if (is.na(ll_val)) {
    #print('issue2')
    return(8000)
  }
  if (abs(ll_val) == Inf) {
    #print('issue3')
    return(8000)
  }
  print(ll_val)
  - ll_val
}

