
r <- seq(-3, 3, length.out = 50)
test_r1 <- function(r,d,h, nu1, nu2) {
 # print(r)
  diff_r <- r < 0
  r_scaled <- ifelse(diff_r, -r, r)
  if (d %% 2 == 0) {
    sign_r <- rep(1, length(r))
  } else {
    sign_r <- r_scaled/r
  }
  # sign_r * (Re(besselJ(r_scaled*h, nu = (d-2)/2) * (1 + complex(imaginary = r))^(-nu1-d/2) *
  #   (1 - complex(imaginary = r))^(-nu2-d/2) * abs(r)))
  # sign_r * (Re(besselJ(r_scaled*h, nu = (d-2)/2) * (1 + complex(imaginary = r))^(-nu1-d/2) * 
  #                (1 - complex(imaginary = r))^(-nu2-d/2) * r))
  v1 <- sign_r * (Re(besselJ(r_scaled*h, nu = (d-2)/2) * (1 + complex(imaginary = r))^(-nu1-d/2) *
                       (1 - complex(imaginary = r))^(-nu2-d/2) * abs(r)))
  # v2 <- sign_r * (Re(besselJ(r_scaled*h, nu = (d-2)/2) * (1 + complex(imaginary = r))^(-nu1-d/2) *
  #                      (1 - complex(imaginary = r))^(-nu2-d/2) * r))
  # plot(v1, v2)
  # all.equal(v1, v2)
  # sum(v1)
  # sum(v2)
  # v2
}
test_r2 <- function(r,d,h, nu1, nu2) {
  diff_r <- r < 0
  r_scaled <- ifelse(diff_r, -r, r)
  if (d %% 2 == 0) {
    sign_r <- rep(1, length(r))
  } else {
    sign_r <- r_scaled/r
  }
  sign_r * Im(besselJ(r_scaled*h, nu = (d-2)/2) * (1 + complex(imaginary = r))^(-nu1-d/2) * 
    (1 - complex(imaginary = r))^(-nu2-d/2) * abs(r))
  # sign_r * Im(besselJ(r_scaled*h, nu = (d-2)/2) * (1 + complex(imaginary = r))^(-nu1-d/2) * 
  #               (1 - complex(imaginary = r))^(-nu2-d/2) * r)
}
test_easy_r1 <- function(r,d,h, nu1, nu2) {
  (Re(exp(complex(imaginary = -r*h)) * (1 + complex(imaginary = r))^(-nu1-d/2) * 
                 (1 - complex(imaginary = r))^(-nu2-d/2)))
}
test_easy_r2 <- function(r,d,h, nu1, nu2) {
  (Im(exp(complex(imaginary = -r*h)) * (1 + complex(imaginary = r))^(-nu1-d/2) * 
        (1 - complex(imaginary = r))^(-nu2-d/2)))
}
#results_re_no_abs <- results_re

test_seq <- seq(-5, 5, by = .1)
grid <-expand.grid(test_seq, test_seq)
results_re <- c()
results_im <- c()
d <- 2
nu1 <- 2.3
nu2 <- .7

for (i in 1:nrow(grid)) {
  h <- sqrt(grid[i,1]^2 + grid[i,2]^2)
  #h <- h * sign(grid[i,1])
  # results_re[i] <- integrate(test_r1, lower = -5, upper = 5, d = d, h = h, 
  #                             nu1 = nu1, nu2=nu2)$value
  # results_im[i] <- integrate(test_r2, lower = -5, upper = 5, d = d, h = h,
  #                           nu1 =nu1, nu2=nu2)$value
  results_re[i] <- integrate(test_r1, lower = -5, upper = 5, d = d, h = h, 
                             nu1 = nu1, nu2=nu2)$value
  results_im[i] <- integrate(test_r2, lower = -5, upper = 5, d = d, h = h, 
                              nu1 = nu1, nu2=nu2)$value
  if (i %% 100 ==0 ) {
    print(i)
  }
}
test_seq_new <- seq(-5, 5, by = .01)
matern_nu_p <- c()
for (i in 1:length(test_seq_new)) {
    results_re[i] <- integrate(test_r1, lower = -20, upper = 20, d = d,
                               h = abs(test_seq_new[i]),
                               nu1 = nu1, nu2=nu2)$value
  matern_nu_p[i] <- integrate(test_r1, lower = -20, upper = 20, d = d,
                              h = abs(test_seq_new[i]),
                              nu1 = nu1/2 + nu2/2, nu2= nu1/2 + nu2/2)$value
  if (i %% 100 ==0 ) {
    print(i)
  }
}
plot(test_seq_new, results_re, type = 'l', ylim = c(0, max(c(results_re, matern_nu_p))))
lines(test_seq_new, matern_nu_p, col = 2)
nu_fixed <- .7
nu_vary <- c(.7, .8, 1, 1.5, 2, 2.5, 3, 3.5)
results_re <- matrix(nrow = length(test_seq_new), ncol = length(nu_vary))
matern_nu_p <- c()
for (i in 1:length(test_seq_new)) {
  for (j in 1:length(nu_vary)) {
    results_re[i,j] <- integrate(test_r1, lower = -20, upper = 20, d = d,
                               h = abs(test_seq_new[i]),
                               nu1 = nu_fixed, nu2=nu_vary[j])$value
  }
  if (i %% 100 ==0 ) {
    print(i)
  }
}
png('d_2_ccov.png', height = 900, width = 1200,pointsize = 12,res = 144)
plot(test_seq_new[test_seq_new >=0], results_re[test_seq_new >=0,1]/ results_re[test_seq_new == 0,1],
     type = 'l', ylim = c(0, 2),
     xlab = 'Lag', ylab = 'Cross-Covariance',
     main = 'nu_1 = .7')
lines(test_seq_new[test_seq_new >=0], results_re[test_seq_new >=0,2]/ results_re[test_seq_new == 0,2], lty = 2)
lines(test_seq_new[test_seq_new >=0], results_re[test_seq_new >=0,3]/ results_re[test_seq_new == 0,3], lty = 3)
lines(test_seq_new[test_seq_new >=0], results_re[test_seq_new >=0,4]/ results_re[test_seq_new == 0,4], lty = 4)
lines(test_seq_new[test_seq_new >=0], results_re[test_seq_new >=0,5]/ results_re[test_seq_new == 0,5], lty = 5)
lines(test_seq_new[test_seq_new >=0], results_re[test_seq_new >=0,6]/ results_re[test_seq_new == 0,6], lty = 6)
lines(test_seq_new[test_seq_new >=0], results_re[test_seq_new >=0,7]/ results_re[test_seq_new == 0,7], col = 2)
#lines(test_seq_new[test_seq_new >=0], results_re[test_seq_new >=0,8]/ results_re[test_seq_new == 0,8], col = 3)
legend('topright', paste0('nu_2 = ', nu_vary[1:7]),
       lty = c(1:6,1), col = c(rep(1, 6),2))
dev.off()

test_seq_new <- seq(-1, 1, by = .001)
nu_vary <- c(2.1, 2.2, 2.3)
results_re <- matrix(nrow = length(test_seq_new), ncol = length(nu_vary))
matern_nu_p <- c()
for (i in 1:length(test_seq_new)) {
  for (j in 1:length(nu_vary)) {
    results_re[i,j] <- integrate(test_r1, lower = -20, upper = 20, d = d,
                                 h = abs(test_seq_new[i]),
                                 nu1 = nu_fixed, nu2=nu_vary[j])$value
  }
  if (i %% 100 ==0 ) {
    print(i)
  }
}
plot(test_seq_new, results_re[,1]/ results_re[test_seq_new == 0,1],
     type = 'l', ylim = c(.8, 1.2))
lines(test_seq_new, results_re[,2]/ results_re[test_seq_new == 0,2], col = 2)
lines(test_seq_new, results_re[,3]/ results_re[test_seq_new == 0,3], col = 3)
lines(test_seq_new, results_re[,4]/ results_re[test_seq_new == 0,4], col = 4)
lines(test_seq_new, results_re[,5]/ results_re[test_seq_new == 0,5], col = 5)
lines(test_seq_new, results_re[,6]/ results_re[test_seq_new == 0,6], col = 6)
lines(test_seq_new, results_re[,7]/ results_re[test_seq_new == 0,7], col = 7)

test_seq_new <- seq(-20, 20, by = .02)
nu_vary <- c(3.5, 3.8, 4.0, 4.8, 10.1)
results_re <- matrix(nrow = length(test_seq_new), ncol = length(nu_vary))
matern_nu_p <- c()
for (i in 1:length(test_seq_new)) {
  for (j in 1:length(nu_vary)) {
    results_re[i,j] <- integrate(test_r1, lower = -10, upper = 10, d = d,
                                 h = abs(test_seq_new[i]),
                                 nu1 = nu_fixed, nu2=nu_vary[j])$value
  }
  if (i %% 100 ==0 ) {
    print(i)
  }
}
plot(test_seq_new, results_re[,1]/ results_re[test_seq_new == 0,1],
     type = 'l')
plot(test_seq_new, results_re[,2]/ results_re[test_seq_new == 0,2],
     type = 'l')
plot(test_seq_new, results_re[,3]/ results_re[test_seq_new == 0,3],
     type = 'l')
plot(test_seq_new, results_re[,4]/ results_re[test_seq_new == 0,4],
     type = 'l')
plot(test_seq_new, results_re[,5]/ results_re[test_seq_new == 0,5],
     type = 'l')

library(ggplot2)
summary(results_im)
ggplot(data = cbind(grid, results_re), aes(x = Var1, y = Var2, fill = results_re))+
  geom_raster()+scale_fill_viridis_b()
ggplot(data = cbind(grid, results_im), aes(x = Var1, y = Var2, fill = results_im))+
  geom_raster()+scale_fill_viridis_b()
test(.3, d = 2, h = 2, nu1 = 1, nu2 = .2)
test(.3, d = 2, h = 2, nu1 = 1, nu2 = 1)
test(.3, d = 1, h = 2, nu1 = .6, nu2 = 1.2)
besselJ(r*h, nu = (d-2)/2) * (1 + complex(imaginary = r))^(-nu1-d/2) * 
  (1 - complex(imaginary = r))^(-nu2-d/2)*r


