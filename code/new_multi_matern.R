
source('/Users/drewyarger/multi_matern_source.R')
A <- complex(real = rnorm(1), 
             imaginary = rnorm(1))

plot_cc <- function(plot_seq =  seq(-5, 5, by = .05), cons, nu) {
  params <- par()$mfrow
  cov_val_lag <- sapply(1:length(plot_seq), function(x) cross_cov(0, plot_seq[x], nu = nu, z_ij = cons))
  plot(plot_seq, cov_val_lag, type = 'l', xlab = 'lag', main = 'Cross Covariance')
}
# plot_cc(cons = complex(real = 1, imaginary = .2), nu = .6)
# plot_cc(cons = complex(real = .5, imaginary = .2), nu = .6)
# plot_cc(cons = complex(real = .5, imaginary = .5), nu = .6)
# plot_cc(cons = complex(real = .5, imaginary = -.5), nu = .6)
# plot_cc(cons = complex(real = -.5, imaginary = 1), nu = .6)
# plot_cc(cons = complex(real = 0, imaginary = 1), nu = .6)


plot_cc(cons = complex(real = 0, imaginary = 1), nu = 1)
plot_cc(cons = complex(real = .5, imaginary = .5), nu = 1)
plot_cc(cons = complex(real = -.2, imaginary = -.5 ), nu = 1)

plot_cc(cons = complex(real = 0, imaginary = 1), nu = .2)
plot_cc(cons = complex(real = .5, imaginary = .5), nu = .2)
plot_cc(cons = complex(real = -.2, imaginary = -.5 ), nu = .2)
