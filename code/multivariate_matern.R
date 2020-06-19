# 
# 
# 
# f <- function(x,h,nu) {
#   return(sin(h*x)/(1 + x^2)^(nu + 1/2))
# }
# 
# f_cos <- function(x,h,nu) {
#   return(cos(h*x)/(1 + x^2)^(nu + 1/2))
# }
# 
# f(1, 1, 1/2)
# 
# h <- 4
# nu <- 3/2
# interval = .05
# x_vec <- seq(0,10, by =interval)
# vec <- c()
# cos_vec <- c()
# for (j in 1:length(x_vec)) {
#   vec[j] <- f(x_vec[j],h, nu)
#   cos_vec[j] <- f_cos(x_vec[j],h, nu)
# } 
# plot(x_vec, vec, ylim = c(-1, 1), type = 'l')
# lines(x_vec, cos_vec, col = 2)
# # intapprox <- sum(vec)/length(x_vec)
# # intapprox_cos <- sum(cos_vec*interval) * (gamma(nu + 1/2) * 2^nu)/(pi^(1/2) * h^nu)
# integrate(f, lower = 0, upper = 100, h = h, nu = nu)$value
# integrate(f_cos, lower = 0, upper = 100, h = h, nu = nu)$value# * (gamma(nu + 1/2) * 2^nu)/(pi^(1/2) * h^nu)
# besselK(h, nu = nu)
# besselI(h, nu = nu)
# f_cos(x = .5, h = 1, nu = nu)
# 
# vals <- expand.grid(h = seq(0, 10, by = .2), nu =seq(.1, 5, by = .2))
# vals_mat <- matrix(nrow = nrow(vals), ncol = 2)
# for (j in 1:nrow(vals)) {
#   h <- vals[j,1]
#   nu <- vals[j,2]
#   vals_mat[j,] <- c(integrate(f, lower = 0, upper = 90, h = h, nu = nu)$value,
#                     integrate(f_cos, lower = 0, upper = 100, h = h, nu = nu)$value
#                   #  besselI(h, nu = nu)
#                     )
#   
# }
# 
# library(ggplot2)
# vals_mat <- as.data.frame(vals_mat)
# colnames(vals_mat) <- c('sin', 'cos')
# vals_mat_long <- tidyr::gather(cbind(vals,vals_mat), type, value, -c(1,2))
# library(viridis)
# 
# ggplot(vals_mat_long, aes(x = h, y = nu, fill = value))+
#   geom_raster()+
#   scale_fill_viridis(trans = 'log')+
#   facet_wrap(~type)
# 
# 
# ggplot(vals_mat_long[vals_mat_long$type == 'sin',], aes(x = h, y = value))+
#   geom_line()+
#   facet_wrap(~nu)
#   #scale_fill_viridis(trans = 'log')+
#   #facet_wrap(~type)
# 
# ggplot(vals_mat_long[vals_mat_long$type == 'cos',], aes(x = h, y = value))+
#   geom_line()+
#   facet_wrap(~nu)


# build cross covariance model given matrix A
nu <- .01
A <- matrix(nrow = 2, ncol = 2,
            complex(real = c(1, .2, .2, 1), 
                    imaginary = c(.2, .1, .2, .1)))
Re(A); Im(A)
i <- 1
j <- 2
z_ij <- as.complex(A[i,] %*% Conj(A[j,]))
Re(z_ij);Im(z_ij)

# or just by defining z_ij
z_ij <- complex(real = .02, imaginary =  .0001)
matern_cov <- function(s,t, nu) {
  (gamma(nu + 1/2) * 2^nu)/(pi^(1/2) * abs(t-s)^nu) * besselK(abs(t-s), nu = nu) 
}
library(RandomFieldsUtils)
struve_version <- function(s,t,nu) {
  abs(t-s)^nu * pi^(3/2)/(cos(pi * nu) * gamma(nu + 1/2)) * 2^(-nu-1) * 
    (besselI(abs(t-s), nu = nu) - RandomFieldsUtils::struveL(abs(t-s), nu = -nu))
}
cross_cov <- function(s,t, nu, z_ij) {
  Re(z_ij) * matern_cov(s,t,nu) -2 * 
    Im(z_ij) * struve_version(s,t,nu)
} # integral is only approximated
test_seq <- seq(-5, 5, by = .1)
grid <- expand.grid(test_seq, test_seq)
cov_val <- sapply(1:nrow(grid), function(x) {
  cross_cov(grid$Var1[x], grid$Var2[x], nu = nu, z_ij = z_ij)
})
library(fields)
image.plot(matrix(cov_val, nrow = length(test_seq)), x = test_seq, y = test_seq,
           xlab = 's', ylab = 't',
           zlim = c(min(cov_val[cov_val > -Inf], na.rm = T),
                    max(cov_val[cov_val < Inf], na.rm = T)))




# sin_version <- function(x,s,t,nu) {
#   return(sin( (t-s)*x)/(1 + x^2)^(nu + 1/2))
# }
# s_test <- .1
# t_test <- .5
# struve_version(s_test, t_test, nu)
# integrate(sin_version, lower = 0, upper = 100, s = s_test, t = t_test, nu = nu)$value
# these confirm that the wolfram-alpha version numerically seems right
