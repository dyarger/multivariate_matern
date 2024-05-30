create_distance_matrices <- function(loc, pred_loc) {
  dist <- fields::rdist(loc)
  dist_pred <- fields::rdist(pred_loc)
  dist_pred_data <- fields::rdist(pred_loc, loc)
  diag(dist) <- 0
  diag(dist_pred) <- 0
  
  for (i in 1:nrow(dist_pred)) {
    for (j in 1:nrow(dist_pred)) {
      if (i == j) {
        
      } else if (pred_loc[i] < pred_loc[j]) {
        dist_pred[i,j] = -dist_pred[i,j]
      }
    }
  }
  
  for (i in 1:nrow(dist_pred_data)) {
    for (j in 1:ncol(dist_pred_data)) {
      if (pred_loc[i] < loc[j]) {
        dist_pred_data[i,j] = -dist_pred_data[i,j]
      }
    }
  }
  
  for (i in 1:nrow(dist)) {
    for (j in 1:nrow(dist)) {
      if (i == j) {
        
      } else if (loc[i] < loc[j]) {
        dist[i,j] = -dist[i,j]
      }
    }
  }
  
  dist_all <- rbind(cbind(dist, t(dist_pred_data)),
                    cbind(dist_pred_data, dist_pred))
  dist_all <- rdist(c(loc, pred_loc))
  for (i in 1:nrow(dist_all)) {
    for (j in 1:nrow(dist_all)) {
      if (i == j) {
        
      } else if (c(loc, pred_loc)[i] < c(loc, pred_loc)[j]) {
        dist_all[i,j] = -dist_all[i,j]
      }
    }
  }
  list(dist_all, dist_pred, dist_pred_data, dist)
}


construct_matrix_alt <- function(nu1, nu2, a1, a2, 
                             grid_info,
                             Sigma_re, Sigma_im, dist_mat) {
  C_val <- fft_1d_alt(grid_info = grid_info, nu1 = nu1, nu2 = nu2, 
                      a1 = a1, a2 = a2, Sigma_re = Sigma_re, Sigma_im = Sigma_im)
  yout <- approx(xout = as.vector(dist_mat), y = C_val[,2], x = C_val[,1], method = 'linear')$y
  matrix(yout, nrow = (nrow(dist_mat)), ncol = (nrow(dist_mat)))
}
construct_bivariate_matrix_mm <- function(nu1, nu2, nu12, a1, a2, a12,
                                       grid_info = grid_info,
                                       Sigma11, Sigma12, Sigma22, dist_mat, nugget1, nugget2) {
  C1 <- construct_matrix_alt(nu1 = nu1, nu2 = nu1, a1 = a1, a2 = a1,
                         Sigma_re = Re(Sigma11), Sigma_im = Im(Sigma11), grid_info = grid_info, dist_mat)
  diag(C1) <- diag(C1) + nugget1
  C12 <- construct_matrix_alt(nu1 = nu12, nu2 = nu12, a1 = a12, a2 = a12,
                              Sigma_re = Re(Sigma12), Sigma_im = Im(Sigma12), grid_info = grid_info, dist_mat)
  C2 <- construct_matrix_alt(nu1 = nu2, nu2 = nu2, a1 = a2, a2 = a2,
                             Sigma_re = Re(Sigma22), Sigma_im = Im(Sigma22), grid_info = grid_info, dist_mat)
  diag(C2) <- diag(C2) + nugget2
  rbind(cbind(C1, C12), cbind(t(C12), C2))
}
construct_bivariate_matrix_alt <- function(nu1, nu2, a1, a2, 
                                       grid_info = grid_info,
                                       Sigma11, Sigma12, Sigma22, dist_mat, nugget1, nugget2) {
  C1 <- construct_matrix_alt(nu1 = nu1, nu2 = nu1, a1 = a1, a2 = a1,
                             Sigma_re = Re(Sigma11), Sigma_im = Im(Sigma11), grid_info = grid_info, dist_mat)
  diag(C1) <- diag(C1) + nugget1
  C12 <- construct_matrix_alt(nu1 = nu1, nu2 = nu2, a1 = a1, a2 = a2,
                              Sigma_re = Re(Sigma12), Sigma_im = Im(Sigma12), grid_info = grid_info, dist_mat)
  C2 <- construct_matrix_alt(nu1 = nu2, nu2 = nu2, a1 = a2, a2 = a2,
                             Sigma_re = Re(Sigma22), Sigma_im = Im(Sigma22), grid_info = grid_info, dist_mat)
  diag(C2) <- diag(C2) + nugget2
  rbind(cbind(C1, C12), cbind(t(C12), C2))
}

construct_matrix_orig <- function(nu1, nu2, a1, a2, 
                             grid_info,
                             Sigma_re, Sigma_im, dist_mat) {
  C_val <- fft_1d(grid_info = grid_info, nu1 = nu1, nu2 = nu2, 
                  a1 = a1, a2 = a2, Sigma_re = Sigma_re, Sigma_im = Sigma_im)
  yout <- approx(xout = as.vector(dist_mat), y = C_val[,2], x = C_val[,1], method = 'linear')$y
  matrix(yout, nrow = (nrow(dist_mat)), ncol = (nrow(dist_mat)))
}
construct_bivariate_matrix_orig <- function(nu1, nu2, a1, a2, 
                                       grid_info = grid_info,
                                       Sigma11, Sigma12, Sigma22, dist_mat, nugget1, nugget2) {
  C1 <- construct_matrix_orig(nu1 = nu1, nu2 = nu1, a1 = a1, a2 = a1,
                              Sigma_re = Re(Sigma11), Sigma_im = Im(Sigma11), grid_info = grid_info, dist_mat)
  diag(C1) <- diag(C1) + nugget1
  C12 <- construct_matrix_orig(nu1 = nu1, nu2 = nu2, a1 = a1, a2 = a2,
                               Sigma_re = Re(Sigma12), Sigma_im = Im(Sigma12), grid_info = grid_info, dist_mat)
  C2 <- construct_matrix_orig(nu1 = nu2, nu2 = nu2, a1 = a2, a2 = a2,
                              Sigma_re = Re(Sigma22), Sigma_im = Im(Sigma22), grid_info = grid_info, dist_mat)
  diag(C2) <- diag(C2) + nugget2
  rbind(cbind(C1, C12), cbind(t(C12), C2))
}

ll_fun <- function(par, dist_one, grid_info, response) {
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
  cov_mat <- construct_bivariate_matrix_orig(nu1 = nu1, nu2 = nu2, a1 = a1, a2 = a2, 
                                             Sigma11 = Sigma11, Sigma22 = Sigma22, 
                                             Sigma12 = complex(real = Sigma12re, imaginary = Sigma12im),
                                             dist_mat = dist_one, grid_info = grid_info,
                                             nugget1 = nugget1, nugget2 = nugget2)
  inv_all <- solve(cov_mat)
  ll_val <- (-nrow(cov_mat)/2 * log(2 * pi) - 1/2 * as.double(t(response) %*% inv_all %*% response) - 
               1/2 * as.vector(determinant(cov_mat)$modulus))
  - ll_val
}

ll_fun_re <- function(par, dist_one, grid_info, response) {
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
  cov_mat <- construct_bivariate_matrix_orig(nu1 = nu1, nu2 = nu2, a1 = a1, a2 = a2, 
                                             Sigma11 = Sigma11, Sigma22 = Sigma22, 
                                             Sigma12 = complex(real = Sigma12re, imaginary = 0),
                                             dist_mat = dist_one, grid_info = grid_info,
                                             nugget1 = nugget1, nugget2 = nugget2)
  inv_all <- solve(cov_mat)
  ll_val <- (-nrow(cov_mat)/2 * log(2 * pi) - 1/2 * as.double(t(response) %*% inv_all %*% response) - 
               1/2 * as.vector(determinant(cov_mat)$modulus))
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
  if (min(eig_values, na.rm = T) < .0001) {
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
  if (min(eig_values, na.rm = T) < .0001) {
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

