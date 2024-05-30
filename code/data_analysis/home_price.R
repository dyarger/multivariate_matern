library(dplyr)
test <- read.csv('~/Downloads/weekly_housing_market_data_most_recent.tsv000', sep = '\t')
head(test)
sf <- test[test$region_name == 'San Diego, CA metro area',]
sf_week <- sf[sf$duration == '1 weeks',]

save(sf_week, file = 'data/housing_SD_1w.RData')

lm1 <- stats::loess(sf_week, formula = median_sale_price~
                      as.numeric(as.Date(period_begin) ))
resid1 <- lm1$residuals
lm2 <- stats::loess(sf_week, formula = inventory~
                      as.numeric(as.Date(period_begin) ))
resid2 <- lm2$residuals
ccf_im <- ccf(resid1, resid2, lag.max = 30 )

source('code/simulation/estimation_source_1d.R')
source('code/multi_matern_source.R')

response <- c(resid1 / sd(resid1), resid2 / sd(resid2))
n <- length(resid1)
n

dist <- matrix(0, n, n)
for (i in 1:n) {
  for (j in 1:n) {
    dist[i,j] <- i - j
  }
}
grid_info <- create_grid_info_1d(2^16, 1800)

test_optim_re <- optim(
  par = c(log(c(1.5, 1.5, .05, .05, .2, .2, var(response[1:n]), var(response[-c(1:n)]))),
          cor(resid1, resid2)),
  fn = ll_fun_re, dist_one = dist, response = response, 
  lower = c(log(c(.001, .001, .00005, .00005,  0.00005, 0.00005, var(response[1:n])/4,  var(response[-c(1:n)])/4)), -1),
  upper = c(log(c(5, 5, 1, 1, 4, 4, var(response[1:n]) * 4,  var(response[-c(1:n)])*4)), 1),
  method = 'L-BFGS-B',
  grid_info = grid_info,
  control = list(parscale = c(rep(1, 8), .1))
)
-test_optim_re$value

test_optim_im <- optim(
  par = c(test_optim_re$par, .01),
  fn = ll_fun, dist_one = dist, response = response, 
  lower = c(log(c(.001, .001, .00005, .00005,  0.00005, 0.00005, 
                  var(response[1:n])/4,  var(response[-c(1:n)])/4)), -1, -1),
  upper = c(log(c(5, 5, 1, 1, 4, 4, var(response[1:n]) * 4, 
                  var(response[-c(1:n)])*4)), 1, 1),
  method = 'L-BFGS-B',
  grid_info = grid_info,
  control = list(parscale = c(rep(1, 8), .1, .1))
)
-test_optim_im$value


test_optim_re_opp <- optim(
  par = c(test_optim_re$par),
  fn = ll_fun_re, dist_one = -dist, response = response, 
  lower = c(log(c(.001, .001, .00005, .00005,  0.00005, 0.00005, var(response[1:n])/4,  var(response[-c(1:n)])/4)), -1),
  upper = c(log(c(5, 5, 1, 1, 4, 4, var(response[1:n]) * 4,  var(response[-c(1:n)])*4)), 1),
  method = 'L-BFGS-B',
  grid_info = grid_info,
  control = list(parscale = c(rep(1, 8), .1))
)
-test_optim_re_opp$value

test_optim_im_opp <- optim(
  par = c(test_optim_re_opp$par, .01),
  fn = ll_fun, dist_one = -dist, response = response, 
  lower = c(log(c(.001, .001, .00005, .00005,  0.00005, 0.00005, 
                  var(response[1:n])/4,  var(response[-c(1:n)])/4)), -1, -1),
  upper = c(log(c(5, 5, 1, 1, 4, 4, var(response[1:n]) * 4, 
                  var(response[-c(1:n)])*4)), 1, 1),
  method = 'L-BFGS-B',
  grid_info = grid_info,
  control = list(parscale = c(rep(1, 8), .1, .1))
)
-test_optim_im_opp$value

save(test_optim_re_opp, test_optim_re, test_optim_im, test_optim_im_opp,
     file = 'results/housing_params.RData')

load('results/housing_params.RData')

if (test_optim_re_opp$value < test_optim_re$value) {
  test_optim_re_old <- test_optim_re
  test_optim_re <- test_optim_re_opp
}

if (test_optim_im_opp$value < test_optim_im$value) {
  test_optim_im_old <- test_optim_im
  test_optim_im <- test_optim_im_opp
  switch <- T
} else{
  switch <- F
}


cc1 <- fft_1d(grid_info = grid_info, nu1 = exp(test_optim_im$par[1]), nu2 = exp(test_optim_im$par[1]), 
              a1 = exp(test_optim_im$par[3]), a2 = exp(test_optim_im$par[3]), 
              Sigma_re = exp(test_optim_im$par[7]) / 
                (exp(test_optim_im$par[7]) + exp(test_optim_im$par[5])),
              Sigma_im = 0) 
plot(cc1, type = 'l', xlim = c(-20, 20))

cc1_re <- fft_1d(grid_info = grid_info, nu1 = exp(test_optim_re$par[1]), nu2 = exp(test_optim_re$par[1]), 
              a1 = exp(test_optim_re$par[3]), a2 = exp(test_optim_re$par[3]), 
              Sigma_re = exp(test_optim_re$par[7]) / 
                (exp(test_optim_re$par[7]) + exp(test_optim_re$par[5])),
              Sigma_im = 0) 

cc2 <- fft_1d(grid_info = grid_info, nu1 = exp(test_optim_im$par[2]), nu2 = exp(test_optim_im$par[2]), 
              a1 = exp(test_optim_im$par[4]), a2 = exp(test_optim_im$par[4]), 
              Sigma_re = exp(test_optim_im$par[8]) / 
                (exp(test_optim_im$par[8]) + exp(test_optim_im$par[6])), 
              Sigma_im = 0)
cc2_re <- fft_1d(grid_info = grid_info, nu1 = exp(test_optim_re$par[2]), nu2 = exp(test_optim_re$par[2]), 
              a1 = exp(test_optim_re$par[4]), a2 = exp(test_optim_re$par[4]), 
              Sigma_re = exp(test_optim_re$par[8]) / 
                (exp(test_optim_re$par[8]) + exp(test_optim_re$par[6])), 
              Sigma_im = 0)


cc12 <- fft_1d(grid_info = grid_info, nu1 = exp(test_optim_im$par[1]), 
               nu2 = exp(test_optim_im$par[2]), 
               a1 = exp(test_optim_im$par[3]), a2 = exp(test_optim_im$par[4]), 
               Sigma_re = sqrt(exp(test_optim_im$par[7])*exp(test_optim_im$par[8])) * test_optim_im$par[9] /
                 (sqrt(exp(test_optim_im$par[7]) + exp(test_optim_im$par[5])) *
                    sqrt(exp(test_optim_im$par[8]) + exp(test_optim_im$par[6]))),
               Sigma_im = sqrt(exp(test_optim_im$par[7])*exp(test_optim_im$par[8])) * test_optim_im$par[10] /
                 (sqrt(exp(test_optim_im$par[7]) + exp(test_optim_im$par[5])) *
                    sqrt(exp(test_optim_im$par[8]) + exp(test_optim_im$par[6])))) 
cc12_re <- fft_1d(grid_info = grid_info, nu1 = exp(test_optim_re$par[1]), 
               nu2 = exp(test_optim_re$par[2]), 
               a1 = exp(test_optim_re$par[3]), a2 = exp(test_optim_re$par[4]), 
               Sigma_re = sqrt(exp(test_optim_re$par[7])*exp(test_optim_re$par[8])) * test_optim_re$par[9] /
                 (sqrt(exp(test_optim_re$par[7]) + exp(test_optim_re$par[5])) *
                    sqrt(exp(test_optim_re$par[8]) + exp(test_optim_re$par[6]))),
               Sigma_im = 0) 

ccf_im <- ccf(resid1, resid2, lag.max = 30 )
library(ggplot2)
theme_set(theme_bw() + theme(text = element_text(size = 16)))

if (switch) {
  est_ccov_im <- rev(cc12[,2])
  est_ccov_re <- rev(cc12_re[,2])
} else {
  est_ccov_im <- cc12[,2]
  est_ccov_re <- cc12_re[,2]
}
line_cex <- .65

ggplot(data = data.frame(date = as.Date(sf_week$period_begin), 
                               acf = c(response[1:n], response[-c(1:n)],
                                       sf_week$median_sale_price, sf_week$inventory), 
                               type = factor(rep(c('Residuals', 'Residuals', 'Median Sale Price ($)',
                                            'Inventory (Units)'), each = n), 
                                            levels = c(c('Median Sale Price ($)',
                                                         'Inventory (Units)', 
                                                         'Residuals'))),
                               variable = factor(rep(c('Median\nSale\nPrice', 'Inventory'), each = n), 
                                                 levels = c('Median\nSale\nPrice', 'Inventory')))) +
  geom_line(aes(x = date, y = acf, color = variable, linetype = variable,
                group = variable), linewidth = line_cex) +
  facet_wrap(~type, ncol = 1, scales  = 'free_y') + 
  #scale_x_continuous(limits =c(-30, 30)) +
  labs(x = 'Date',
       y = '', 
       color = 'Variable', 
       linetype = 'Variable') +
  theme(legend.position = 'right', text = element_text(size = 10),
        strip.text = element_text(size = 8))+
  guides(color = guide_legend(nrow = 2), linetype = guide_legend(nrow = 2))
ggsave('images/redfin_data.png', height = 3, width = 4)

ggplot(data = data.frame(date = as.Date(sf_week$period_begin), 
                         acf = c(response[1:n], response[-c(1:n)],
                                 sf_week$median_sale_price, sf_week$inventory), 
                         type = rep(c('Residuals', 'Residuals', 'Median Sale Price ($)',
                                      'Inventory (Units)'), each = n),
                         variable = rep(c('Median Sale Price', 'Inventory'), each = n))) +
  geom_line(aes(x = date, y = acf, color = variable, linetype = variable,
                group = variable), linewidth = line_cex) +
  facet_wrap(~type, ncol = 3, scales  = 'free_y') + 
  #scale_x_continuous(limits =c(-30, 30)) +
  labs(x = 'Date',
       y = 'Data', 
       color = 'Variable', 
       linetype = 'Variable') +
  theme(legend.position = 'bottom', text = element_text(size = 12),
        strip.text = element_text(size = 10))+
  guides(color = guide_legend(nrow = 1), linetype = guide_legend(nrow = 1))
ggsave('images/redfin_data_jsm.png', height = 3.5, width = 7.5)

var(lm1$residuals)  * exp(test_optim_re$par)[7]
var(lm1$residuals)  * exp(test_optim_im$par)[7]
var(lm2$residuals) * exp(test_optim_re$par)[8]
var(lm2$residuals) * exp(test_optim_im$par)[8]
test_optim_re$par
test_optim_im$par
exp(test_optim_re$par) 
exp(test_optim_im$par)

ggplot(data = rbind(data.frame(lag = ccf_im$lag, acf = ccf_im$acf, type = 'Empirical'),
                    data.frame(lag = cc12[,1], acf = est_ccov_im, type = 'Est Complex'),
                    data.frame(lag = cc12[,1], acf = est_ccov_re, type = 'Est Real')) %>%
         mutate(type = factor(type, levels = c('Empirical', 'Est Real', 'Est Complex')))) +
  geom_line(aes(x = lag, y = acf, color = type, linetype = type), linewidth = line_cex) +
  scale_x_continuous(limits = c(-30, 30)) +
  labs(x = 'Lag h (weeks)',
       y = expression(atop('Cross-correlation','Cor['*Y[1]*'(t+h), '*Y[2]*'(t)]')), 
       color = 'Model', 
       linetype = 'Model') +
  theme(legend.position = 'bottom', text = element_text(size = 12)) +
  guides(color = guide_legend(nrow = 1), linetype = guide_legend(nrow = 1))
ggsave('images/redfin_cc.png', height = 3, width = 4.5)

ggplot(data = rbind(data.frame(lag = ccf_im$lag, acf = ccf_im$acf, type = 'Empirical'),
                    data.frame(lag = cc12[,1], acf = est_ccov_im, type = 'Est Complex'),
                    data.frame(lag = cc12[,1], acf = est_ccov_re, type = 'Est Real')) %>%
         mutate(type = factor(type, levels = c('Empirical', 'Est Real', 'Est Complex')))) +
  geom_line(aes(x = lag, y = acf, color = type, linetype = type), linewidth = line_cex) +
  scale_x_continuous(limits = c(-30, 30)) +
  labs(x = 'Lag (weeks)',
       y = 'Cross-correlation', 
       color = 'Model', 
       linetype = 'Model') +
  theme(legend.position = 'bottom') +
  guides(color = guide_legend(nrow = 2), linetype = guide_legend(nrow = 2))
ggsave('images/redfin_cc_jsm.png', height = 3.5, width = 7.5)


acf1 <- acf(resid1, lag.max = 30)
acf2 <- acf(resid2, lag.max = 30)

pchisq(-2 * (-test_optim_re$value + test_optim_im$value), 
       df = 1, lower.tail = FALSE)
ggplot(data = rbind(data.frame(lag = acf1$lag, acf = acf1$acf, type = 'Empirical'),
                    data.frame(lag = cc1[,1], acf = cc1[,2], type = 'Est Complex'),
                    data.frame(lag = cc1[,1], acf = cc1_re[,2], type = 'Est Real')) %>%
         mutate(type = factor(type, levels = c('Empirical', 'Est Real', 'Est Complex')))) +
  geom_line(aes(x = lag, y = acf, color = type, linetype = type), linewidth = line_cex) +
  scale_x_continuous(limits = c(0, 30)) +
  labs(x = 'Lag (weeks)',
       y = 'Correlation', 
       color = 'Model', 
       linetype = 'Model') +
  theme(legend.position = 'bottom')+
  guides(color = guide_legend(nrow = 2), linetype = guide_legend(nrow = 2))
ggsave('images/redfin_c1.png', height = 3, width = 4)

ggplot(data = rbind(data.frame(lag = acf2$lag, acf = acf2$acf, type = 'Empirical'),
                    data.frame(lag = cc2[,1], acf = cc2[,2], type = 'Est Complex'),
                    data.frame(lag = cc2[,1], acf = cc2_re[,2], type = 'Est Real')) %>%
         mutate(type = factor(type, levels = c('Empirical', 'Est Real', 'Est Complex')))) +
  geom_line(aes(x = lag, y = acf, color = type, linetype = type), linewidth = line_cex) +
  scale_x_continuous(limits = c(0, 30)) +
  labs(x = 'Lag (weeks)',
       y = 'Correlation', 
       color = 'Model', 
       linetype = 'Model') +
  theme(legend.position = 'bottom')+
  guides(color = guide_legend(nrow = 2), linetype = guide_legend(nrow = 2))
ggsave('images/redfin_c2.png', height = 3, width = 4)

test_optim_re_alt <- optim(
  par = c(log(c(1.5, 1.5, .05, .05, .2, .2, var(response[1:n]), var(response[-c(1:n)]))),
          cor(resid1, resid2)),
  fn = ll_fun_alt_re, dist_one = dist, response = response, 
  lower = c(log(c(.001, .001, .00005, .00005,  0.00005, 0.00005, var(response[1:n])/4,  var(response[-c(1:n)])/4)), -1),
  upper = c(log(c(5, 5, 1, 1, 4, 4, var(response[1:n]) * 4,  var(response[-c(1:n)])*4)), 1),
  method = 'L-BFGS-B',
  grid_info = grid_info,
  control = list(parscale = c(rep(1, 8), .1))
)
-test_optim_re_alt$value

test_optim_im_alt <- optim(
  par = c(test_optim_re_alt$par, .01),
  fn = ll_fun_alt, dist_one = dist, response = response, 
  lower = c(log(c(.001, .001, .00005, .00005,  0.00005, 0.00005, 
                  var(response[1:n])/4,  var(response[-c(1:n)])/4)), -1, -1),
  upper = c(log(c(5, 5, 1, 1, 4, 4, var(response[1:n]) * 4, 
                  var(response[-c(1:n)])*4)), 1, 1),
  method = 'L-BFGS-B',
  grid_info = grid_info,
  control = list(parscale = c(rep(1, 8), .1, .1))
)
-test_optim_im_alt$value
