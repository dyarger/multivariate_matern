
source('code/multi_matern_source.R')
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


plot_function <- function(s,t,nu) {
  if(s-t == 0) {
    return(0)
  }
 sign(t-s) * (abs(t-s))^nu*
    (besselI(abs(t-s), nu = nu) - struve(abs(t-s), -nu))
}

plot_seq <- seq(-8, 8, by = .05)
nu25 <- sapply(1:length(plot_seq), function(x) plot_function(0, plot_seq[x], nu = .25))
nu75 <- sapply(1:length(plot_seq), function(x) plot_function(0, plot_seq[x], nu = .75))
nu125 <- sapply(1:length(plot_seq), function(x) plot_function(0, plot_seq[x], nu = 1.25))
nu01 <- sapply(1:length(plot_seq), function(x) plot_function(0, plot_seq[x], nu = .1))
nu175 <- sapply(1:length(plot_seq), function(x) plot_function(0, plot_seq[x], nu = 1.75))
nu225 <- sapply(1:length(plot_seq), function(x) plot_function(0, plot_seq[x], nu = 2.25))
nu275 <- sapply(1:length(plot_seq), function(x) plot_function(0, plot_seq[x], nu = 2.75))

png('example_fun.png', width = 480*3, height = 480*3)
par(mar=c(5.1, 4.1+5, 4.1, 2.1))
plot(plot_seq, nu01, type = 'l', ylim = c(-1.2, 1.2),# ylim = c(-.8, .8), 
     xlab = 'Lag', ylab = 'Function Value',
     lwd = 4, cex.lab =4, cex.axis = 2)
lines(plot_seq, nu25, type = 'l', col = 2, lwd = 4)
lines(plot_seq, nu75, type = 'l', col = 3, lwd=4)
lines(plot_seq, nu125, type = 'l', col = 4, lwd =4)
lines(plot_seq, nu175, type = 'l', lty = 2, lwd =4)
lines(plot_seq, nu225, type = 'l', lty = 3, lwd =4)
lines(plot_seq, nu275, type = 'l', lty = 4, lwd =4)
legend('bottomleft', c('nu = .01', 'nu = .25', 'nu = .75', 'nu = 1.25',
                        'nu = 1.75', 'nu = 2.25', 'nu = 2.75'),
       lty = c(rep(1, 4),2,3,4), col = c(1:4,1,1,1), cex = 3, lwd = rep(4,7))
dev.off()

plot_matern <- function(s,t,nu) {
   (abs(t-s))^nu*besselK(abs(t-s), nu)
}

plot_seq <- seq(-8, 8, by = .05)
nu25 <- sapply(1:length(plot_seq), function(x) plot_function(0, plot_seq[x], nu = 1.25))
maternnu25 <- sapply(1:length(plot_seq), function(x) plot_matern(0, plot_seq[x], nu = 1.25))

plot(plot_seq, nu25, type = 'l', ylim = c(-1.2, 1.2),# ylim = c(-.8, .8), 
     xlab = 'Lag', ylab = 'Function Value',
     lwd = 4, cex.lab =4, cex.axis = 2)
lines(plot_seq, maternnu25, type = 'l', col = 2, lwd = 4)



par(mar=c(5.1, 4.1 + .5, 4.1, 2.1))
fields::image.plot(pred_vals_section, ylim = c(1, 0), xlab = 'Latitude', ylab = 'Pressure (decibar)', axes=F,
                   legend.cex = 1.8, cex.lab = 1.8, legend.lab = '°C/day', legend.line = 4,
                   legend.mar = 7, axis.args=list(cex.axis=1.6))
axis(2, at = seq(0, 2000, by = 250)/2000, seq(0, 2000, by = 250), cex.axis = 1.6)
axis(1, at=(seq(-50.5, 50.5, by = 10) + 50.5)/103.5, labels=seq(-50.5, 50.5, by = 10),
     cex.axis = 1.6)

plot_cc(cons = complex(real = 0, imaginary = 1), nu = 1)
plot_cc(cons = complex(real = 0, imaginary = 1), nu = .500000001)
plot_cc(cons = complex(real = 0, imaginary = 1), nu = .48)
plot_cc(cons = complex(real = .5, imaginary = .5), nu = 1)
plot_cc(cons = complex(real = -.2, imaginary = -.5 ), nu = 1)

plot_cc(cons = complex(real = 0, imaginary = 1), nu = 5)
plot_cc(cons = complex(real = .5, imaginary = .5), nu = .2)
plot_cc(cons = complex(real = -.2, imaginary = -.5 ), nu = .2)

