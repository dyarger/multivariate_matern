# numerically approximate integral
struve_version <- function(h,nu, d) {
  if(h == 0) {
    return(0)
  }
  sign(h)*(abs(h))^nu * pi * 2^(-nu-d/2) / (gamma(nu + d/2) * cos(pi * nu))*
    (besselI(x = abs(h), nu = nu) - RandomFieldsUtils::struveL(x = abs(h), nu = -nu))
}
library(fields)
a <- c(1,0)
a_perp <- c(0,1)
d <- length(a)
nu = .75
x_seq <- seq(-5, 5, by = .1)
x_grid <- as.matrix(expand.grid(x_seq, x_seq))
x_results_new <- sapply(1:nrow(x_grid), function(y) {
  #x_grid[y,] <- c(1,0)
  a_norm <- sqrt(sum(a_perp^2))
  h_norm <- sqrt(sum(x_grid[y,]^2))
  gamma <- sum(a_perp * x_grid[y,])/a_norm/h_norm
  if (is.nan(gamma)) {gamma = 1}
  s1 <- struve_version(h_norm, nu = nu, d = 2) 
  sign2 <- sign(sum(a_perp * x_grid[y,]))
  sign1 <- sign(sum(a * x_grid[y,]))
  if (sign2==0) {sign2=1}
  if (sign1==0) {sign1=1}
  tv <- abs(gamma) * h_norm
  s2 <- struve_version(tv, nu = nu, d = 2)
 #-s1 * sign1
  sign1 * s2
  tv
  #s2
   #s2
  # abs(gamma) * h_norm
  # abs(gamma)^(-d/2 +1)
  #abs(gamma)
})
png('images/testnu0_75.png')
image.plot(matrix(x_results_new, nrow = length(x_seq)), main = 'nu = 0.75')
dev.off()

image.plot(t(matrix(ifelse(h_vals < 0, p1_vals-test_x_new,
                           -(p1_vals-test_x_new)) * 2, nrow = length(cov_seq))))
image.plot(t(matrix(ifelse(h_vals < 0, p1_vals-test_x_new,
                           -(p1_vals-test_x_new)), nrow = length(cov_seq))))


image.plot(t(matrix(ifelse(h_vals < 0, p1_vals,
                           -(p1_vals)) * 2, nrow = length(cov_seq))))

image.plot(t(matrix(ifelse(h_vals < 0, -test_x_new,
                           +(test_x_new)) * 2, nrow = length(cov_seq))))


image.plot(2*matrix(x_results, nrow = length(x_seq)), main = 'nu = 0.75',
           zlim = c(-1,1))
image.plot(matrix(test_save, nrow = length(x_seq)), main = 'nu = 0.75',zlim = c(-1,1))
image.plot(2*matrix(x_results, nrow = length(x_seq)) - 
             matrix(test_save, nrow = length(x_seq)), main = 'nu = 0.75',
           zlim = c(-1,1))

te <- sapply(x_seq, function(x) {
  gamma = 2
  struve_version(abs(gamma) * x, nu = nu, d = 2)
})
plot(x_seq, te)

x_results2 <- sapply(1:nrow(x_grid), function(y) {
  #y = y+1
  #x_grid[1,] <- c(0,1)
  a_norm <- sqrt(sum(a_perp^2))
  h_norm <- sqrt(sum(x_grid[y,]^2))
  gamma <- sum(a_perp * x_grid[y,])/a_norm/h_norm
  if (is.nan(gamma)) {gamma = 1}
  #struve_version(h_norm, nu = nu, d = 2) #+
  h_norm <- sign(sum(a_perp * x_grid[y,]))
  struve_version(abs(gamma) * h_norm, nu = nu, d = 2)
})

image.plot(matrix(x_results2, nrow = length(x_seq)), main = 'nu = 0.75')
image.plot(matrix(x_results+x_results2 , nrow = length(x_seq)), main = 'nu = 0.75')




