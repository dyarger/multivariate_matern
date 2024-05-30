
files <- list.files('results/sim_results_2d/', pattern = '300')
result_df <- list()
for (i in 1:length(files)) {
  load(paste0('results/sim_results_2d/', files[i]))
  
  vals <- strsplit(files[i], '\\_')[[1]]
  
  result_df[[i]] <- data.frame(model = c('A', 'B', 'C')[as.integer(vals[2])],
                               Re_true = vals[3],
                               Im_true = vals[4], 
                               theta_star = vals[5],
                               sim = strsplit(vals[15], '\\.')[[1]][1], 
                               A_Co_ll = -test_optim_im$value,
                               A_Re_ll = -test_optim_re$value,
                               B_Co_ll = -test_optim_im_psi_fixed$value,
                               B_Re_ll = -test_optim_re_psi_fixed$value,
                               A_Co_Sigma12_Re = test_optim_im$par[9] * sqrt(exp(test_optim_im$par[7])*exp(test_optim_im$par[8])),
                               A_Co_Sigma12_Im = test_optim_im$par[10] * sqrt(exp(test_optim_im$par[7])*exp(test_optim_im$par[8])),
                               A_Re_Sigma12_Re = test_optim_re$par[9] * sqrt(exp(test_optim_re$par[7])*exp(test_optim_re$par[8])),
                               B_Co_Sigma12_Re = test_optim_im_psi_fixed$par[9] * sqrt(exp(test_optim_im_psi_fixed$par[7])*exp(test_optim_im_psi_fixed$par[8])),
                               B_Co_Sigma12_Im = test_optim_im_psi_fixed$par[10] * sqrt(exp(test_optim_im_psi_fixed$par[7])*exp(test_optim_im_psi_fixed$par[8])),
                               B_Re_Sigma12_Re = test_optim_re_psi_fixed$par[9] * sqrt(exp(test_optim_re_psi_fixed$par[7])*exp(test_optim_re_psi_fixed$par[8])))
}
library(dplyr)
df <- dplyr::bind_rows(result_df)
table(df$theta_star)

df %>%
  group_by(model, Re_true, Im_true, theta_star) %>%
  mutate(p_val_A = pchisq(-2 * (A_Re_ll - A_Co_ll), df = 1, lower.tail = FALSE),
         p_val_B = pchisq(-2 * (B_Re_ll - B_Co_ll), df = 1, lower.tail = FALSE)) %>%
  summarize(reject_A = mean(p_val_A < 0.05),
            reject_B = mean(p_val_B < 0.05)) %>%
  arrange(model, theta_star, as.numeric(Im_true), as.numeric(Re_true))

df %>%
  group_by(model, Re_true, Im_true, theta_star) %>%
  summarize_at(vars(contains('Sigma')), function(x) round(mean((x)), 2)) %>%
  arrange(model, theta_star, as.numeric(Im_true), as.numeric(Re_true)) %>%
  filter(model == 'A', theta_star == '0') %>%
  select(starts_with('A'))

df %>%
  group_by(model, Re_true, Im_true, theta_star) %>%
  summarize_at(vars(contains('Sigma')), function(x) round(mean(abs(x)), 2)) %>%
  arrange(model, theta_star, as.numeric(Im_true), as.numeric(Re_true)) %>%
  filter(model == 'A')

df %>%
  group_by(model, Re_true, Im_true, theta_star) %>%
  summarize_at(vars(contains('Sigma')), function(x) round(mean(x), 2)) %>%
  arrange(model, theta_star, as.numeric(Im_true), as.numeric(Re_true)) %>%
  filter(model == 'A', theta_star == 0)

df %>%
  group_by(model, Re_true, Im_true, theta_star) %>%
  summarize_at(vars(contains('Sigma')), function(x) round(sd(x), 2)) %>%
  arrange(model, theta_star, as.numeric(Im_true), as.numeric(Re_true)) %>%
  filter(model == 'A', theta_star == 0)

