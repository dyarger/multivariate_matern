
library(ncdf4)
source('code/multi_matern_source.R')
file_folder <- '~/Downloads/SOCCOM_HRQC_MLR_netcdf_20220322/'
#files <- list.files(file_folder, pattern = '.nc')[100]
files <- "5906030_HRQC.nc"
data_list <- list()
for (j in 1:length(files)) {
  file <- nc_open(paste0(file_folder, files[j]))
  
  # core info
  lat <- ncvar_get(file, 'Lat')
  lon <- ncvar_get(file, 'Lon')
  lon_qc <- ncvar_get(file, 'Lat_QFA')
  temp <- ncvar_get(file, 'Temperature')
  psal <- ncvar_get(file, 'Salinity')
  pressure <- ncvar_get(file, 'Pressure')
  pressure_QC <- ncvar_get(file, 'Pressure_QFA')
  temp_QC <- ncvar_get(file, 'Temperature_QFA')
  psal_QC <- ncvar_get(file, 'Salinity_QFA')
  day <- ncvar_get(file, 'mon_day_yr')
  time <- ncvar_get(file, 'hh_mm')
  
  # BGC info
  oxy <- ncvar_get(file, 'Oxygen')
  oxy_sat <- ncvar_get(file, 'OxygenSat') # Calculation assumes atmospheric pressure= 1013.25 mbar
  oxy_QC <- ncvar_get(file, 'Oxygen_QFA')
  oxy_sat_QC <- ncvar_get(file, 'OxygenSat_QFA')
  pdens <- ncvar_get(file, 'Sigma_theta')
  pdens_QC <- ncvar_get(file, 'Sigma_theta_QFA')
  depth <- ncvar_get(file, 'Depth')
  nitrate <- ncvar_get(file, 'Nitrate')
  nitrate_QC <- ncvar_get(file, 'Nitrate_QFA')
  pH <- ncvar_get(file, 'pHinsitu')
  pH_QC <- ncvar_get(file, 'pHinsitu_QFA')
  chl <- ncvar_get(file, 'Chl_a')
  chl_QC <- ncvar_get(file, 'Chl_a_QFA')
  poc <- ncvar_get(file, 'POC')
  poc_QC <- ncvar_get(file, 'POC_QFA')
  
  nprof <- ncol(pressure)
  npres <- nrow(pressure)
  
  if (length(dim(pressure)) == 1) {
    repeated_vars <- data.frame('float' =as.numeric(ncvar_get(file, 'Cruise')),
                                'profile' =  as.numeric(ncvar_get(file, 'Station'))-1,
                                'latitude' =                    rep(lat, each = npres),
                                'longitude'=                rep(lon, each = npres),
                                'longitude_QC' =                 rep(lon_qc, each = npres),
                                'day'=                rep(day, each = npres),
                                'time' = rep(time, each = npres))
  } else {
    repeated_vars <- data.frame('float' = rep(as.numeric(ncvar_get(file, 'Cruise')), times = npres * nprof),
                                'profile' =  rep(as.numeric(ncvar_get(file, 'Station')-1), each = npres),
                                'latitude' =                    rep(lat, each = npres),
                                'longitude'=                rep(lon, each = npres),
                                'longitude_QC' =                 rep(lon_qc, each = npres),
                                'day'=                rep(day, each = npres),
                                'time' = rep(time, each = npres))
  }
  
  mat_test <- data.frame(
    repeated_vars,
    # core
    'pressure' = as.double(pressure),
    'temp' = as.double(temp),
    'psal' = as.double(psal),
    'pdens' = as.double(pdens), 
    #BGC
    'oxy' = as.double(oxy),
    'oxy_sat' = as.double(oxy_sat),
    'nitrate' = as.double(nitrate),
    'pH' = as.double(pH),
    'chl' = as.double(chl),
    'poc' = as.double(poc),
    # Quality flags
    'pressure_QC' = as.double(pressure_QC), 'temp_QC' = as.double(temp_QC),
    'psal_QC' = as.double(psal_QC), 'oxy_QC' = as.double(oxy_QC),
    'pdens_QC' = as.double(pdens_QC), 
    'oxy_sat_QC' = as.double(oxy_sat_QC),
    'nitrate_QC' = as.double(nitrate_QC),
    'pH_QC' = as.double(pH_QC),
    'chl_QC' = as.double(chl),
    'poc_QC' = as.double(poc_QC))
  data_list[[j]] <- mat_test
  print(file[['dim']][['N_PROF']][['len']])
  nc_close(file)
}

df <- data_list[[1]]  %>%
  filter(!is.na(pressure), !is.na(temp), !is.na(psal), 
         pressure < 950) %>%
  mutate(date = as.Date(day, format = '%m/%d/%Y'))
ggplot(data = df, aes(x = pressure, y = temp, color = profile)) +
  geom_point()
ggplot(data = df, aes(x = pressure, y = psal, color = profile)) +
  geom_point()
p_vec_max <- c(15, 105, 205, 305)
p_vec_min <- p_vec_max - 10
df_single_level <- df %>%
  mutate(p_group1 = findInterval(pressure, vec = p_vec_max),
         p_group2 = findInterval(pressure, p_vec_min)) %>%
  filter(p_group1 != p_group2) %>%
  group_by(profile, p_group1) %>%
  summarize(temp = mean(temp), psal = mean(psal), 
            date = date[1]) %>%
  ungroup() %>%
  group_by(p_group1) %>%
  mutate(temp_0 = temp - mean(temp),
         psal_0 = psal - mean(psal)) %>%
  ungroup()
ggplot(data = df_single_level) +
  geom_boxplot(aes(x = p_group1, y = temp_0, group = p_group1))

df_wide <- pivot_wider(df_single_level, id_cols = c(profile, date), names_from =  p_group1, values_from = psal_0)
head(df_wide)
cov(df_wide %>% dplyr::select(-profile, -date))
cor(df_wide %>% dplyr::select(-profile, -date))

nobs <- nrow(df_wide)
response <- as.double(unlist(df_wide %>% dplyr::select(-profile, -date)))
dist_matrix <- matrix(df_wide[['date']], nrow = nobs, ncol = nobs, byrow = F)  -
  matrix(df_wide[['date']], nrow = nobs, ncol = nobs, byrow = T) 
dist_upper_tri <- dist_matrix[upper.tri(dist_matrix, diag = T)]

likelihood <- function(theta, response, dist_matrix, dist_upper_tri) {
  K <- length(response) / nrow(dist_matrix)
  print(exp(theta))
  nu_vals <- exp(theta[1:K]);theta <- theta[-c(1:K)]
  range_vals <- exp(theta[1:K]);theta <- theta[-c(1:K)]
  nugget_vals <- exp(theta[1:K]);theta <- theta[-c(1:K)]
  var_diag <- exp(theta[1:K]);theta <- theta[-c(1:K)]
  var_offdiag <- theta
  Sigma <- matrix(nrow = K, ncol = K, 0)
  diag(Sigma) <- var_diag
  Sigma[upper.tri(Sigma)] <- var_offdiag
  lt <- lower.tri(dist_matrix, diag = T)
  ut <- upper.tri(dist_matrix, diag = T)
  if (tail(eigen(Sigma, only.values = T)[['values']], 1) <= 0.0001) {
    return(10^6)
  }
  
  list_mat <- list()
  for (i in 1:K) {
    list_mat[[i]] <- list()
    for (j in 1:K) {
      if (j < i) {
        next
      } else if (i == j) {
        cov_val <- sapply(dist_upper_tri, function(x) {
          whitt_version(x, nu1= nu_vals[i], nu2 = nu_vals[j], 
                        c11 = Sigma[i,i], c12 = Sigma[i,j], c2 = Sigma[j,j], 
                        a1 = range_vals[i], a2 = range_vals[j])
        })
        cov_mat <- dist_matrix
        cov_mat[ut] <- cov_val[1,]
        cov_mat[lt] <- t(cov_mat)[lt]
        diag(cov_mat) <- diag(cov_mat) + nugget_vals[i]
        list_mat[[i]][[j]] <- cov_mat
      } else {
        cov_val <- sapply(dist_matrix, function(x) {
          whitt_version(x, nu1= nu_vals[i], nu2 = nu_vals[j], 
                        c11 = Sigma[i,i], c12 = Sigma[i,j], c2 = Sigma[j,j], 
                        a1 = range_vals[i], a2 = range_vals[j])
        })
        cov_mat <- dist_matrix
        cov_mat[matrix(T, nrow = nrow(dist_matrix),
                         ncol = ncol(dist_matrix))] <- cov_val[2,]
        list_mat[[i]][[j]] <- cov_mat
      }
    }
  }
  test_list_mat <- lapply(list_mat, function(x) {
    test <- matrix(nrow = 110, ncol = 0)
    for (i in 1:length(x)) {
      if (length(x[[i]]) > 0) {
        test <- cbind(test, x[[i]])
      } else {
        test <- cbind(test,  matrix(nrow = 110, ncol = 110))
      }
    }
    return(test)
  })
  final_mat <- test_list_mat[[1]]
  for (i in 2:length(test_list_mat)) {
    final_mat <- rbind(final_mat, test_list_mat[[i]])
  }
  
  c_chol <- base::chol(cov_mat)
  v2 <- .Internal(backsolve(r = c_chol, x = response, 
                            k = ncol(c_chol), upper.tri = T, transpose = T))
  quad_form <- sum(v2^2)
  det_val <-  2* sum(log(diag(c_chol)))
  print(quad_form+ det_val)
  quad_form+ det_val
}


K <- length(response) / nrow(dist_matrix)
theta_nu <- runif(K)
theta_range <- rep(.005, K)
theta_nugget <- rep(.00001, K)
theta_diag <- rep(.2, K)
theta_offdiag <- rep(.02, (K * K  - K)/2)
theta <- c(theta_nu, theta_range, theta_nugget, theta_diag, theta_offdiag)
rm(theta)
init_params <- c(rep(.55, K),
                 rep(.005, K),
                 rep(.0001, K),
                 rep(.2, K),
                 exp(rep(.02, K * (K-1)/2)))
max_params  <- c(rep(2.5, K),
                 rep(.03, K),
                 rep(.02, K),
                 rep(3, K),
                 exp(rep(3, K * (K-1)/2)))
min_params  <- c(rep(.02, K),
                 rep(.000001, K),
                 rep(.00000001, K),
                 rep(.0001, K),
                 exp(rep(-2, K * (K-1)/2)))


joint_optim <- optim(par = log(init_params), fn = likelihood,
                     response = response, method = 'L-BFGS-B', 
                     hessian = T, lower = log(min_params),
                     upper = log(max_params), 
                    # control = list(parscale = c(1, 1, 1/50, 1/50, 1, 1, 1, 1)),
                     dist_matrix = dist_matrix, dist_upper_tri = dist_upper_tri)

off_diag <- sqrt(1 - .2^2)
re_off_diag <- .2
a1 <- 1.1
a2 <- 1.05
A <- matrix(nrow = 2, complex(real = c(a1,re_off_diag, re_off_diag, a2), imaginary = c(0, off_diag, -off_diag, 0)))
eigen(A)
corpcor::is.positive.definite(A)
tail(eigen(A, only.values = T)$values, 1) >0
abs(A)
abs(complex(real = re_off_diag, imaginary = off_diag))
sqrt(off_diag^2 + re_off_diag^2)
A
#  233.1487
Sigma = matrix(complex(real = c(exp(joint_optim[['par']][1]),
                                joint_optim[['par']][3],
                                joint_optim[['par']][3],
                                exp(joint_optim[['par']][2])),
                       imaginary = c(0, joint_optim[['par']][4],
                                     -joint_optim[['par']][4],
                                     0)), nrow = 2)
abs(Sigma)
abs(Sigma)[1,2]/sqrt(abs(Sigma)[1,1] * abs(Sigma)[2,2])
Re(Sigma[1,2])/abs(Sigma)[1,2]
Im(Sigma[1,2])/abs(Sigma)[1,2]
nu <- exp(joint_optim[['par']])[5]
a <- exp(joint_optim[['par']])[6]
nugget1 <- exp(joint_optim[['par']])[7]
nugget2 <- exp(joint_optim[['par']])[8]


cov_mat <- imaginary_covariance_matrix(dist_matrix, dist_upper_tri,
                                       nu1 = nu, 
                                       nu2 = nu,
                                       c11 = (Sigma[1,1]), 
                                       c22 = (Sigma[2,2]), 
                                       c12 = complex(real = Re(Sigma[1,2]), 
                                                     imaginary = Im(Sigma[1,2])), 
                                       a1 = a, a2 = a, nugget1 = nugget1,
                                       nugget2 = nugget2)
image.plot(cov_mat)
indexes_temp <- 1:nobs
indexes_psal <- setdiff(1:(2 * nobs), indexes_temp)
cov_mat_temp <- cov_mat[1:nobs, 1:nobs]
cov_mat_psal <- cov_mat[-c(1:nobs), -c(1:nobs)]
cov_mat_cc <- cov_mat[-c(1:nobs), 1:nobs]
pred_psal <- cov_mat_cc %*% solve(cov_mat_temp, response1)
plot(pred_psal, response2)
sum(response2^2)
sum((response2- pred_psal)^2)

pred_temp <- t(cov_mat_cc) %*% solve(cov_mat_psal, response2)
plot(pred_temp, response1)
sum(response1^2)
sum((response1- pred_temp)^2)
abline(b = 1, a = 0)

response <- c(response1, response2)
r_indexes_add <- sample(nobs + 1:nobs, 10)
indexes_use <- c(indexes_temp, r_indexes_add)
pred_using_temp <- cov_mat[,indexes_use] %*% 
  solve(cov_mat[indexes_use, indexes_use], 
        response[indexes_use])
plot(pred_using_temp[indexes_psal], response[indexes_psal])
plot(pred_using_temp[indexes_temp], response[indexes_temp])
sum( response[indexes_psal]^2)
sum(( response[indexes_psal]- pred_using_temp[indexes_psal])^2)

r_indexes_add <- sample(1:nobs, 4)
indexes_use <- c(indexes_psal, r_indexes_add)
pred_using_psal <- cov_mat[,indexes_use] %*% 
  solve(cov_mat[indexes_use, indexes_use], 
        response[indexes_use])
plot(pred_using_psal[indexes_psal], response[indexes_psal])
plot(pred_using_psal[indexes_temp], response[indexes_temp])
sum( response[indexes_psal]^2)
sum(( response[indexes_psal]- pred_using_psal[indexes_psal])^2)
sum( response[indexes_temp]^2)
sum(( response[indexes_temp]- pred_using_psal[indexes_temp])^2)

