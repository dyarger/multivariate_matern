
library(ncdf4)
source('code/multi_matern_source.R')
file_folder <- '~/Downloads/SOCCOM_HRQC_MLR_netcdf_20220322/'
files <- list.files(file_folder, pattern = '.nc')[93]#49, 93
#files <- "5904470_HRQC.nc"
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
         pressure < 950, psal_QC != 8) %>%
  mutate(date = as.Date(day, format = '%m/%d/%Y'))


df_single_level <- df %>%
  filter(pressure > 145, pressure < 150) %>%
  group_by(profile) %>%
  summarize(temp = mean(temp), psal = mean(psal), 
            date = date[1]) %>%
  ungroup()
X_val <- cbind(1, sin(julian(df_single_level$date) / 366* 2*pi), 
               cos(julian(df_single_level$date) / 366* 2*pi),
               sin(2*julian(df_single_level$date) / 366* 2*pi), 
               cos(2*julian(df_single_level$date) / 366* 2*pi))
mean_fun_temp <- solve(crossprod(X_val), crossprod(X_val, df_single_level$temp))
mean_fun_psal <- solve(crossprod(X_val), crossprod(X_val, df_single_level$psal))
df_single_level[['temp_0']] <- df_single_level$temp - as.double(X_val %*% mean_fun_temp)
df_single_level[['psal_0']] <- df_single_level$psal - as.double(X_val %*% mean_fun_psal)

nobs <- nrow(df_single_level)
response1 <- df_single_level[['temp_0']]
response2 <- df_single_level[['psal_0']]

grid <- expand.grid('V1' = 1:length(response2),
                    'V2' = 1:length(response2))


temp_val <- data.frame('V' = 1:length(response2), 'temp' = response1)
psal_val <- data.frame('V' = 1:length(response2), 'psal' = response2)

grid_temp <-grid %>%
  left_join(rename(temp_val, V1 = V), by = 'V1') %>%
  left_join(rename(temp_val, V2 = V), by = 'V2') %>%
  mutate(prod = temp.x * temp.y,
         lag_val = (V1 - V2)) %>%
  group_by(lag_val) %>%
  summarize(cov = mean(prod))
grid_psal <-grid %>%
  left_join(rename(psal_val, V1 = V), by = 'V1') %>%
  left_join(rename(psal_val, V2 = V), by = 'V2') %>%
  mutate(prod = psal.x * psal.y,
         lag_val = (V1 - V2)) %>%
  group_by(lag_val) %>%
  summarize(cov = mean(prod))

grid_cc <-grid %>%
  left_join(rename(temp_val, V1 = V), by = 'V1') %>%
  left_join(rename(psal_val, V2 = V), by = 'V2') %>%
  mutate(prod = psal * temp,
         lag_val = (V1 - V2)) %>%
  group_by(lag_val) %>%
  summarize(cov = mean(prod))

df_sum <- data.frame(grid_cc, 'temp_cov' = grid_temp[['cov']],
                     'psal_cov' = grid_psal[['cov']]) %>%
  rename(temp_psal_cov = cov) %>%
  pivot_longer(cols = ends_with('cov'))
label_df <- data.frame('name'  =  c('temp_cov', 'psal_cov', 'temp_psal_cov'),
                       'label' = factor(c('Temperature covariogram',
                                          'Salinity covariogram', 
                                          'Temperature and salinity cross-covariogram'),
                                        c('Temperature covariogram',
                                          'Temperature and salinity cross-covariogram',
                                          'Salinity covariogram')))



a <- ggplot(data = df_sum %>% filter(lag_val > -35, lag_val < 35) %>%
         left_join(label_df), aes(x = lag_val, y = value)) + 
  geom_point() + 
  geom_line() + 
  facet_wrap(~label, ncol = 2, scales = 'free_y') + 
  labs(x = 'Lag (number of profiles)',
       y = 'Covariances and cross-covariances')+
  theme_bw()
ggsave(plot = a, 'images/argo_cc.png', height = 4.5, width = 6.5)

label_df <- data.frame('name'  =  c('temp_0', 'psal_0'),
                       'label' = factor(c('Centered temperature (°C)',
                                          'Centered salinity (Practical Salinity Units)'),
                                        c('Centered temperature (°C)',
                                          'Centered salinity (Practical Salinity Units)')))

b <-ggplot(data = df_single_level %>%
         pivot_longer(ends_with('_0')) %>%
           left_join(label_df), aes(x = date, y = value)) + 
  geom_point() + 
  geom_line() + 
  facet_wrap(~label, ncol = 1, scales = 'free_y') + 
  labs(x = 'Date',
       y = 'Centered temperature and salinity values at 150 decibars')+
  theme_bw()
ggsave(plot = b, 'images/argo_plot.png', height = 4.5, width = 6.5)
library(patchwork)
a+b

dist_matrix <- matrix(df_single_level[['date']], nrow = nobs, ncol = nobs, byrow = F)  -
  matrix(df_single_level[['date']], nrow = nobs, ncol = nobs, byrow = T) 
dist_upper_tri <- dist_matrix[upper.tri(dist_matrix, diag = T)]
plot(response1)
plot(response2)

cross_cov

imag_part_struve <- function(h,nu, a) {
  if(h == 0) {
    return(0)
  }
  sign(h) * abs(h)^nu*
    (besselI(abs(h) * a, nu = nu) - struve(abs(h) * a, -nu))
}

real_part_matern <- function(h, nu, a) {
  if (h == 0) {
    h = 10^-18
  }
  abs(h)^nu*
    besselK(a*abs(h), nu = nu) 
}

imaginary_covariance_matrix <- function(dist_matrix, dist_upper_tri, nu1, nu2, c11, 
                                        c12, c22, a1, a2, nugget1, nugget2) {
  nu <- nu1
  a = a1
  cov_val11 <- c11 * sapply(dist_upper_tri, function(x) {
    real_part_matern(h = x, nu = nu1, a  = a1)
  })* 2*pi^(1/2) * a^(2 * nu + 1)/2 / pi /( (2* a)^(nu) * gamma(nu + 1/2))
  cov_val22 <- c22/c11 * cov_val11
  real_part <- sapply(dist_matrix, function(x) {
    real_part_matern(h = x, nu = nu1, a  = a1) 
  })* 2*pi^(1/2) *
    a^(2 * nu + 1)/2 / pi /( (2* a)^(nu) * gamma(nu + 1/2))
  
  im_part <- sapply(dist_matrix, function(x) {
    imag_part_struve(h = x, nu = nu1, a  = a1)
  })*a^(-nu) * a^(2 * nu + 1)/2/pi / cos(nu * pi) * 2^(-nu) * pi^(3/2)/gamma(nu + 1/2)
  
  cross_cov1 <- Re(c12) * real_part + Im(c12) * im_part
  
  lt <- lower.tri(dist_matrix, diag = T)
  ut <- upper.tri(dist_matrix, diag = T)
  cov_mat1 <- dist_matrix
  cov_mat1[ut] <- cov_val11
  cov_mat1[lt] <- t(cov_mat1)[lt]
  diag(cov_mat1) <- diag(cov_mat1) + nugget1

  cov_mat2 <- dist_matrix
  cov_mat2[ut] <- cov_val22
  cov_mat2[lt] <- t(cov_mat2)[lt]
  diag(cov_mat2) <- diag(cov_mat2) + nugget2
  
  cov_mat12 <- dist_matrix
  cov_mat12[matrix(T, nrow = nrow(dist_matrix),
                   ncol = ncol(dist_matrix))] <- cross_cov1
  cov_mat_all <- rbind(cbind(cov_mat1, cov_mat12),
                       cbind(t(cov_mat12), cov_mat2))
}

likelihood <- function(theta, response, dist_matrix, dist_upper_tri) {
  print(exp(theta))
  Sigma <- matrix(nrow = 2, ncol = 2, 
                  complex(real = c(exp(theta[1]), theta[3]*sqrt(exp(theta[1])*exp(theta[2])),
                                   theta[3]*sqrt(exp(theta[1])*exp(theta[2])), exp(theta[2])),
                          imaginary = c(0, theta[4]*sqrt(exp(theta[1])*exp(theta[2])), 
                                        -theta[4]*sqrt(exp(theta[1])*exp(theta[2])), 0)))
  if (tail(eigen(Sigma, only.values = T)[['values']], 1) <= 0.0001) {
    return(10^6)
  }
  cov_mat <- imaginary_covariance_matrix(dist_matrix, dist_upper_tri,
                                         nu1 = exp(theta[5]), 
                                         nu2 = exp(theta[5]),
                                         c11 = exp(theta[1]), 
                                         c22 = exp(theta[2]), 
                                         c12 = complex(real = theta[3], 
                                                       imaginary = theta[4])*
                                           sqrt(exp(theta[1])*exp(theta[2])), 
                                         a1 = exp(theta[6]), a2 = exp(theta[6]),
                                         nugget1 = exp(theta[7]), 
                                         nugget2 = exp(theta[8]))
  c_chol <- base::chol(cov_mat)
  v2 <- .Internal(backsolve(r = c_chol, x = response, 
                            k = ncol(c_chol), upper.tri = T, transpose = T))
  quad_form <- sum(v2^2)
  det_val <-  2* sum(log(diag(c_chol)))
  print(  quad_form+ det_val)
  quad_form+ det_val
}

response <- c(response1, response2)
cov(cbind(response1, response2))
cor(cbind(response1, response2))

joint_optim <- optim(par = c(log(var(response1)*1000),log(var(response2)*1000), 
                             cor(response1, response2), 0,
                             log(.51), log(.0001),
                             log(.06), log(.001)), fn = likelihood,
                     response = response, method = 'L-BFGS-B', 
                     hessian = T, 
                     lower = c(log(var(response1)/100000), log(var(response2)/100000), -.99,-.99, 
                               log(.05),log(.000005), log(.00000002), log(.00000002)),
                     upper = c(log(2000), log(2000), .99, .99, 
                               log(2), log(.01), log(2), log(2)),
                     control = list(parscale = c(1, 1, 1/10,1/10, 1, 1, 1, 1)),
                     dist_matrix = dist_matrix, dist_upper_tri = dist_upper_tri)

Sigma = matrix(complex(real = c(exp(joint_optim[['par']][1]),
                                joint_optim[['par']][3]*sqrt(exp(joint_optim[['par']][1])*
                                                               exp(joint_optim[['par']][2])),
                                joint_optim[['par']][3]*sqrt(exp(joint_optim[['par']][1])*
                                                               exp(joint_optim[['par']][2])),
                                exp(joint_optim[['par']][2])),
                       imaginary = c(0, joint_optim[['par']][4]*sqrt(exp(joint_optim[['par']][1])*
                                                                       exp(joint_optim[['par']][2])),
                                     -joint_optim[['par']][4]*sqrt(exp(joint_optim[['par']][1])*
                                                                     exp(joint_optim[['par']][2])),
                                     0)), nrow = 2)
#Sigma <- Re(Sigma)
abs(Sigma)
abs(Sigma)[1,2]/sqrt(abs(Sigma)[1,1] * abs(Sigma)[2,2])
Re(Sigma[1,2])/abs(Sigma)[1,2]
Im(Sigma[1,2])/abs(Sigma)[1,2]
nu <- exp(joint_optim[['par']])[5]
a <- exp(joint_optim[['par']])[6]
nugget1 <- exp(joint_optim[['par']])[7]
nugget2 <- exp(joint_optim[['par']])[8]


test_seq <- seq(-35*10, 35*10, by = 1)
real_part <- Sigma[1,1] *sapply(test_seq, function(x) {
  real_part_matern(h = x, nu = nu, a  = a) * 2*pi^(1/2) *
    a^(2 * nu + 1)/2 / pi /( (2* a)^(nu) * gamma(nu + 1/2))
})

plot(df_sum$lag_val[df_sum$name == 'temp_cov' & abs(df_sum$lag_val) < 36],
     df_sum$value[df_sum$name == 'temp_cov'& abs(df_sum$lag_val) < 36])
lines(test_seq/10, real_part, col = 2)

real_part <- Sigma[2,2] *sapply(test_seq, function(x) {
  real_part_matern(h = x, nu = nu, a  = a) * 2*pi^(1/2) *
    a^(2 * nu + 1)/2 / pi /( (2* a)^(nu) * gamma(nu + 1/2))
})

plot(df_sum$lag_val[df_sum$name == 'psal_cov' & abs(df_sum$lag_val) < 36],
     df_sum$value[df_sum$name == 'psal_cov'& abs(df_sum$lag_val) < 36])
lines(test_seq/10, real_part, col = 2)


real_part <- sapply(test_seq, function(x) {
  real_part_matern(h = x, nu = nu, a  = a) * 2*pi^(1/2) *
    a^(2 * nu + 1)/2 / pi /( (2* a)^(nu) * gamma(nu + 1/2))
})

im_part <- sapply(test_seq, function(x) {
  imag_part_struve(h = x, nu = nu, a  = a)*
    a^(-nu) * a^(2 * nu + 1)/2/pi / cos(nu * pi) * 2^(-nu) * pi^(3/2)/gamma(nu + 1/2)
})

cross_cov1 <- Re(Sigma[1,2]) * real_part + Im(Sigma[1,2]) * im_part

plot(df_sum$lag_val[df_sum$name == 'temp_psal_cov' & abs(df_sum$lag_val) < 72],
     df_sum$value[df_sum$name == 'temp_psal_cov'& abs(df_sum$lag_val) < 72])
lines(test_seq/10, cross_cov1, col = 2)


real_covariance_matrix <- function(dist_matrix, dist_upper_tri, nu1, nu2, c11, 
                                        c12, c22, a1, a2, nugget1, nugget2) {
  cov_vals <- sapply(dist_matrix, function(x) {
    whitt_version(h = x, nu1 = nu1, nu2 = nu2, c11 = c11, c12 = c12, c2 = c22,a1 = a1,
                  a2 = a2)
  })
  
  cov_mat1 <- matrix(cov_vals[1,], nrow = sqrt(ncol(cov_vals)),
                     ncol = sqrt(ncol(cov_vals)))
  diag(cov_mat1) <- diag(cov_mat1) + nugget1
  
  cov_mat2 <- matrix(cov_vals[3,], nrow = sqrt(ncol(cov_vals)),
                     ncol = sqrt(ncol(cov_vals)))
  diag(cov_mat2) <- diag(cov_mat2) + nugget2
  
  
  cov_mat12 <- matrix(cov_vals[2,], nrow = sqrt(ncol(cov_vals)),
                     ncol = sqrt(ncol(cov_vals)))
  cov_mat_all <- rbind(cbind(cov_mat1, cov_mat12),
                       cbind(t(cov_mat12), cov_mat2))
}

likelihood_real <- function(theta, response, dist_matrix, dist_upper_tri) {
  print(exp(theta))
  Sigma <- matrix(nrow = 2, ncol = 2, c(exp(theta[1]), theta[3]*sqrt(exp(theta[1])*exp(theta[2])),
                                   theta[3]*sqrt(exp(theta[1])*exp(theta[2])), exp(theta[2])))
  if (tail(eigen(Sigma, only.values = T)[['values']], 1) <= 0.0001 |
      (exp(theta[6])/exp(theta[7])> 5) | (exp(theta[7])/exp(theta[6]) > 5)) {
    return(10^6)
  }
  cov_mat <- real_covariance_matrix(dist_matrix,
                                         nu1 = exp(theta[4]), 
                                         nu2 = exp(theta[5]),
                                         c11 = exp(theta[1]), 
                                         c22 = exp(theta[2]), 
                                         c12 = theta[3]*
                                           sqrt(exp(theta[1])*exp(theta[2])), 
                                         a1 = exp(theta[6]), a2 = exp(theta[7]),
                                         nugget1 = exp(theta[8]), 
                                         nugget2 = exp(theta[9]))
  c_chol <- base::chol(cov_mat)
  v2 <- .Internal(backsolve(r = c_chol, x = response, 
                            k = ncol(c_chol), upper.tri = T, transpose = T))
  quad_form <- sum(v2^2)
  det_val <-  2* sum(log(diag(c_chol)))
  print(  quad_form+ det_val)
  quad_form+ det_val
}

joint_optim <- optim(par = c(log(var(response1)),log(var(response2)), 
                             cor(response1, response2), 
                             log(.57), log(.51), log(.001),
                             log(.001),
                             log(.001), log(.001)), fn = likelihood_real,
                     response = response, method = 'L-BFGS-B', 
                     hessian = T, 
                     lower = c(log(var(response1)/1500), log(var(response2)/1500), -1, log(.01),
                               log(.01),
                               log(.000005), log(.000005), 
                               log(.00000002),
                               log(.00000002)),
                     upper = c(log(2), log(2), 1, log(1.5), log(1.5),log(.05),log(.05),
                               log(2), log(2)),
                     control = list(parscale = c(1, 1, 1/10, 1, 1, 1, 1, 1,1)),
                     dist_matrix = dist_matrix, dist_upper_tri = dist_upper_tri)
theta <- joint_optim$par
nu1 <- exp(theta[4])
nu2 <- exp(theta[5])
a1 <- exp(theta[6])
a2 <- exp(theta[7])
nugget1 <- exp(theta[8])
nugget2 <- exp(theta[9])
s11 <- exp(theta[1])
s22 <- exp(theta[2])
s12 <- theta[3] * sqrt(s11 * s22)

cov_vals <- sapply(test_seq, function(x) {
  whitt_version(h = x, nu1 = nu1, nu2 = nu2, c11 = s11, c12 = s12, c2 = s22,a1 = a1,
            a2 = a2)
})

plot(df_sum$lag_val[df_sum$name == 'temp_psal_cov' & abs(df_sum$lag_val) < 36],
     df_sum$value[df_sum$name == 'temp_psal_cov'& abs(df_sum$lag_val) < 36],
     ylim = c(-.1, .1))
lines(-test_seq/10, cov_vals[2,], col = 2)

plot(df_sum$lag_val[df_sum$name == 'temp_cov' & abs(df_sum$lag_val) < 36],
     df_sum$value[df_sum$name == 'temp_cov'& abs(df_sum$lag_val) < 36])
lines(-test_seq/10, cov_vals[1,] + 
        ifelse(test_seq == 0, nugget1, 0), col = 2)

plot(df_sum$lag_val[df_sum$name == 'psal_cov' & abs(df_sum$lag_val) < 36],
     df_sum$value[df_sum$name == 'psal_cov'& abs(df_sum$lag_val) < 36])
lines(-test_seq/10, cov_vals[3,]+ 
        ifelse(test_seq == 0, nugget2, 0), col = 2)

# move to approximate integrals

spec_dens <- function(r, h, nu1, nu2, a1 = 1, a2 = 1, re_z, im_z) {
  val <- exp(complex(imaginary = h * r)) *
    complex(real = a1, imaginary = r)^(-nu1 - 1/2) *
    complex(real = a2, imaginary = -r)^(-nu2 - 1/2) * 
    complex(real = re_z, imaginary = sign(r) * im_z)
  Re(val)
}

spec_dens_symm <- function(r, h, nu, a, re_z) {
  val <- exp(complex(imaginary = h * r)) *
    (a^2 + r^2)^(-nu - 1/2) * re_z
  Re(val)
}

int_covariance_matrix <- function(dist_matrix, dist_upper_tri, nu1, nu2, c11, 
                                   c12, c22, a1, a2, nugget1, nugget2) {
  cov_vals12 <- sapply(dist_matrix, function(x) {
    #print(x)
    integrate(spec_dens, h = x, nu1 = nu1, nu2 = nu2, re_z = Re(c12), im_z = Im(c12),
              a1 = a1, a2 = a2,
              lower = -.2, upper = .2)$value#+
      # integrate(spec_dens, h = x, nu1 = nu1, nu2 = nu2, re_z = Re(c12), im_z = Im(c12),
      #           a1 = a1, a2 = a2,
      #           lower = -.2, upper = 10^-8)$value#*
      #norm_constant(d = 1, x = plot_seq[x], nu_1 = nu1_val, nu_2 = .25, a_1 = 1, a_2 = 1)
  })/2/pi * a1^(nu1 + 1/2) * a2^(nu2 + 1/2)
  
  # test <- sapply(seq(-.01, .01, by = .001), function(x) {
  #   #print(x)
  #   spec_dens(h = 0, r = x, nu1 = nu1, nu2 = nu2, re_z = Re(c12), im_z = Im(c12),
  #             a1 = a1, a2 = a2)#*
  #   #norm_constant(d = 1, x = plot_seq[x], nu_1 = nu1_val, nu_2 = .25, a_1 = 1, a_2 = 1)
  # })
  
  
  
  
  cov_vals1 <- sapply(dist_matrix, function(x) {
    integrate(spec_dens_symm, h = x, nu = nu1, re_z = Re(c11), 
              a = a1,
              lower = -.2, upper = .2)$value#*
    #norm_constant(d = 1, x = plot_seq[x], nu_1 = nu1_val, nu_2 = .25, a_1 = 1, a_2 = 1)
  }) * a1^(2*nu1 + 1)/2/pi
  
  cov_vals2 <- sapply(dist_matrix, function(x) {
    integrate(spec_dens_symm, h = x, nu = nu2, re_z = Re(c22), 
              a = a2,
              lower = -.2, upper = .2)$value#*
    #norm_constant(d = 1, x = plot_seq[x], nu_1 = nu1_val, nu_2 = .25, a_1 = 1, a_2 = 1)
  })* a2^(2*nu2 + 1)/2/pi
  
  cov_mat1 <- matrix(cov_vals1, nrow = (ncol(dist_matrix)),
                     ncol = (ncol(dist_matrix)))
  diag(cov_mat1) <- diag(cov_mat1) + nugget1
  
  cov_mat2 <- matrix(cov_vals2, nrow = (ncol(dist_matrix)),
                     ncol = (ncol(dist_matrix)))
  diag(cov_mat2) <- diag(cov_mat2) + nugget2
  
  
  cov_mat12 <- matrix(cov_vals12, nrow = (ncol(dist_matrix)),
                      ncol = (ncol(dist_matrix)))
  cov_mat_all <- rbind(cbind(cov_mat1, cov_mat12),
                       cbind(t(cov_mat12), cov_mat2))
  # cov_mat_all <- rbind(cbind(cov_mat1, t(cov_mat12)),
  #                      cbind(cov_mat12, cov_mat2))
}


likelihood_int <- function(theta, response, dist_matrix, dist_upper_tri) {
  print(exp(theta))
  Sigma <- matrix(nrow = 2, ncol = 2, c(exp(theta[1]),
                                        complex(real = theta[3], imaginary = theta[4])*
                                          sqrt(exp(theta[1])*exp(theta[2])),
                                        complex(real = theta[3], imaginary = -theta[4])*sqrt(exp(theta[1])*exp(theta[2])),
                                        exp(theta[2])))
  if (tail(eigen(Sigma, only.values = T)[['values']], 1) <= 0.0001 |
      (exp(theta[7])/exp(theta[8])> 5) | (exp(theta[8])/exp(theta[7]) > 5)) {
    return(10^6)
  }
  cov_mat <- int_covariance_matrix(dist_matrix,
                                    nu1 = exp(theta[5]), 
                                    nu2 = exp(theta[6]),
                                    c11 = exp(theta[1]), 
                                    c22 = exp(theta[2]), 
                                    c12 = complex(real = theta[3], imaginary = theta[4])*
                                      sqrt(exp(theta[1])*exp(theta[2])), 
                                    a1 = exp(theta[7]), a2 = exp(theta[8]),
                                    nugget1 = exp(theta[9]), 
                                    nugget2 = exp(theta[10]))
  c_chol <- base::chol(cov_mat)
  v2 <- .Internal(backsolve(r = c_chol, x = response, 
                            k = ncol(c_chol), upper.tri = T, transpose = T))
  quad_form <- sum(v2^2)
  det_val <-  2* sum(log(diag(c_chol)))
  print(  quad_form+ det_val)
  quad_form+ det_val
}

joint_optim <- optim(par = c(log(var(response1)),log(var(response2)), 
                             cor(response1, response2), 0,
                             log(.57), log(.51), log(.00001),
                             log(.00001),
                             log(.001), log(.001)), fn = likelihood_int,
                     response = response, method = 'L-BFGS-B', 
                     hessian = T, 
                     lower = c(log(var(response1)/1500), log(var(response2)/1500), -1, -1,
                               log(.01),log(.01),
                               log(.000005), log(.000005), 
                               log(.00000002),log(.00000002)),
                     upper = c(log(2), log(2), 1, 1, log(1.5), log(1.5),log(.02),log(.02),
                               log(2), log(2)),
                     control = list(parscale = c(1, 1, 1/10,1/10, 1, 1, 1, 1, 1,1)),
                     dist_matrix = dist_matrix, dist_upper_tri = dist_upper_tri)
#save(joint_optim, file = 'argo_data_analysis_int_test_results.RData')


test25 <- sapply(test_seq, function(z) {
  integrate(spec_dens, h = z, nu1 = nu1_val, nu2 = .25, re_z = 0, im_z = 1,
            a1 = 1, a2 = 1,
            lower = 0, upper = 100)$value*
    norm_constant(d = 1, x = plot_seq[x], nu_1 = nu1_val, nu_2 = .25, a_1 = 1, a_2 = 1)
})
























cov_mat <- real_covariance_matrix(dist_matrix, dist_upper_tri,
                                       nu1 = nu1, 
                                       nu2 = nu2,
                                       c11 = s11, 
                                       c22 = s22, 
                                       c12 = s12, 
                                       a1 = a1, a2 = a2, nugget1 = nugget1,
                                       nugget2 = nugget2)
library(fields)
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



likelihood_re_only <- function(theta, response, dist_matrix, dist_upper_tri) {
  print(exp(theta))
  Sigma <- matrix(nrow = 2, ncol = 2, 
                  complex(real = c(exp(theta[1]), theta[3], theta[3], exp(theta[2])),
                          imaginary = c(0, theta[4], -theta[4], 0)))
  if (tail(eigen(Sigma, only.values = T)[['values']], 1) <= 0.0001) {
    return(10^6)
  }
  cov_mat <- imaginary_covariance_matrix(dist_matrix, dist_upper_tri,
                                         nu1 = exp(theta[5]), 
                                         nu2 = exp(theta[5]),
                                         c11 = exp(theta[1]), 
                                         c22 = exp(theta[2]), 
                                         c12 = complex(real = theta[3], 
                                                       imaginary = theta[4]), 
                                         a1 = exp(theta[6]), a2 = exp(theta[6]),
                                         nugget1 = exp(theta[7]), 
                                         nugget2 = exp(theta[8]))
  c_chol <- base::chol(cov_mat)
  v2 <- .Internal(backsolve(r = c_chol, x = response, 
                            k = ncol(c_chol), upper.tri = T, transpose = T))
  quad_form <- sum(v2^2)
  det_val <-  2* sum(log(diag(c_chol)))
  print(  quad_form+ det_val)
  quad_form+ det_val
}