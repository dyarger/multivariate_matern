files <- list.files('results/sim_results_1d/', pattern = '300')
result_df <- list()
for (i in 1:length(files)) {
  load(paste0('results/sim_results_1d/', files[i]))
  
  vals <- strsplit(files[i], '\\_')[[1]]
  response <- pred[,1]
  
  result_df[[i]] <- data.frame(model = c('A', 'B', 'C')[as.integer(vals[2])],
                               Re_true = vals[3],
                               Im_true = vals[4], 
                               sim = strsplit(vals[15], '\\.')[[1]][1], 
                               A_Co_ll = -test_optim$value,
                               A_Re_ll = -test_optim_real$value,
                               B_Co_ll = -test_optim_alt$value,
                               B_Re_ll = -test_optim_real_alt$value,
                               C_Co_ll = -test_optim_mm$value,
                               C_Re_ll = -test_optim_mm_re$value,
                               A_Co_Sigma12_Re = test_optim$par[9] * sqrt(exp(test_optim$par[7])*exp(test_optim$par[8])),
                               A_Co_Sigma12_Im = test_optim$par[10] * sqrt(exp(test_optim$par[7])*exp(test_optim$par[8])),
                               A_Re_Sigma12_Re = test_optim_real$par[9] * sqrt(exp(test_optim_real$par[7])*exp(test_optim_real$par[8])),
                               B_Co_Sigma12_Re = test_optim_alt$par[9] * sqrt(exp(test_optim_alt$par[7])*exp(test_optim_alt$par[8])),
                               B_Co_Sigma12_Im = test_optim_alt$par[10] * sqrt(exp(test_optim_alt$par[7])*exp(test_optim_alt$par[8])),
                               B_Re_Sigma12_Re = test_optim_real_alt$par[9] * sqrt(exp(test_optim_real_alt$par[7])*exp(test_optim_real_alt$par[8])),
                               C_Co_Sigma12_Re = test_optim_mm$par[11] * sqrt(exp(test_optim_mm$par[10])*exp(test_optim_mm$par[9])),
                               C_Co_Sigma12_Im = test_optim_mm$par[12] * sqrt(exp(test_optim_mm$par[10])*exp(test_optim_mm$par[9])),
                               C_Re_Sigma12_Re = test_optim_mm_re$par[11] * sqrt(exp(test_optim_mm_re$par[9])*exp(test_optim_mm_re$par[10])),
                               A_Co_resp1_pred12 = sqrt(mean((pred[1:100,2] - response[1:100])^2)),
                               A_Co_resp1_pred1 = sqrt(mean((pred[1:100,3] - response[1:100])^2)),
                               A_Co_resp1_pred2 = sqrt(mean((pred[1:100,4] - response[1:100])^2)),
                               A_Re_resp1_pred12 = sqrt(mean((pred_re[1:100,2] - response[1:100])^2)),
                               A_Re_resp1_pred1 = sqrt(mean((pred_re[1:100,3] - response[1:100])^2)),
                               A_Re_resp1_pred2 = sqrt(mean((pred_re[1:100,4] - response[1:100])^2)),
                               B_Co_resp1_pred12 = sqrt(mean((pred_alt[1:100,2] - response[1:100])^2)),
                               B_Co_resp1_pred1 = sqrt(mean((pred_alt[1:100,3] - response[1:100])^2)),
                               B_Co_resp1_pred2 = sqrt(mean((pred_alt[1:100,4] - response[1:100])^2)),
                               B_Re_resp1_pred12 = sqrt(mean((pred_re_alt[1:100,2] - response[1:100])^2)),
                               B_Re_resp1_pred1 = sqrt(mean((pred_re_alt[1:100,3] - response[1:100])^2)),
                               B_Re_resp1_pred2 = sqrt(mean((pred_re_alt[1:100,4] - response[1:100])^2)),
                               C_Co_resp1_pred12 = sqrt(mean((pred_mm[1:100,2] - response[1:100])^2)),
                               C_Co_resp1_pred1 = sqrt(mean((pred_mm[1:100,3] - response[1:100])^2)),
                               C_Co_resp1_pred2 = sqrt(mean((pred_mm[1:100,4] - response[1:100])^2)),
                               C_Re_resp1_pred12 = sqrt(mean((pred_mm_re[1:100,2] - response[1:100])^2)),
                               C_Re_resp1_pred1 = sqrt(mean((pred_mm_re[1:100,3] - response[1:100])^2)),
                               C_Re_resp1_pred2 = sqrt(mean((pred_mm_re[1:100,4] - response[1:100])^2)),
                               A_Co_resp2_pred12 = sqrt(mean((pred[101:200,2] - response[101:200])^2)),
                               A_Co_resp2_pred1 = sqrt(mean((pred[101:200,3] - response[101:200])^2)),
                               A_Co_resp2_pred2 = sqrt(mean((pred[101:200,4] - response[101:200])^2)),
                               A_Re_resp2_pred12 = sqrt(mean((pred_re[101:200,2] - response[101:200])^2)),
                               A_Re_resp2_pred1 = sqrt(mean((pred_re[101:200,3] - response[101:200])^2)),
                               A_Re_resp2_pred2 = sqrt(mean((pred_re[101:200,4] - response[101:200])^2)),
                               B_Co_resp2_pred12 = sqrt(mean((pred_alt[101:200,2] - response[101:200])^2)),
                               B_Co_resp2_pred1 = sqrt(mean((pred_alt[101:200,3] - response[101:200])^2)),
                               B_Co_resp2_pred2 = sqrt(mean((pred_alt[101:200,4] - response[101:200])^2)),
                               B_Re_resp2_pred12 = sqrt(mean((pred_re_alt[101:200,2] - response[101:200])^2)),
                               B_Re_resp2_pred1 = sqrt(mean((pred_re_alt[101:200,3] - response[101:200])^2)),
                               B_Re_resp2_pred2 = sqrt(mean((pred_re_alt[101:200,4] - response[101:200])^2)),
                               C_Co_resp2_pred12 = sqrt(mean((pred_mm[101:200,2] - response[101:200])^2)),
                               C_Co_resp2_pred1 = sqrt(mean((pred_mm[101:200,3] - response[101:200])^2)),
                               C_Co_resp2_pred2 = sqrt(mean((pred_mm[101:200,4] - response[101:200])^2)),
                               C_Re_resp2_pred12 = sqrt(mean((pred_mm_re[101:200,2] - response[101:200])^2)),
                               C_Re_resp2_pred1 = sqrt(mean((pred_mm_re[101:200,3] - response[101:200])^2)),
                               C_Re_resp2_pred2 = sqrt(mean((pred_mm_re[101:200,4] - response[101:200])^2)))
}
library(dplyr)
df <- dplyr::bind_rows(result_df)



head(df)
dim(df)

df %>%
  group_by(model, Re_true, Im_true) %>%
  mutate(p_val_A = pchisq(-2 * (A_Re_ll - A_Co_ll), df = 1, lower.tail = FALSE),
         p_val_B = pchisq(-2 * (B_Re_ll - B_Co_ll), df = 1, lower.tail = FALSE),
         p_val_C = pchisq(-2 * (C_Re_ll - C_Co_ll), df = 1, lower.tail = FALSE)) %>%
  summarize(reject_A = mean(p_val_A < 0.05),
            reject_B = mean(p_val_B < 0.05),
            reject_C = mean(p_val_C < 0.05)) %>%
  arrange(model, as.numeric(Im_true), as.numeric(Re_true))

library(ggplot2)
df %>%
  group_by(model, Re_true, Im_true) %>%
  mutate(p_val_A = pchisq(-2 * (A_Re_ll - A_Co_ll), df = 1, lower.tail = FALSE),
         p_val_B = pchisq(-2 * (B_Re_ll - B_Co_ll), df = 1, lower.tail = FALSE),
         p_val_C = pchisq(-2 * (C_Re_ll - C_Co_ll), df = 1, lower.tail = FALSE), 
         reject_A = p_val_C < 0.05, 
         reject_B = p_val_C < 0.05, 
         reject_C = p_val_C < 0.05) %>%
  filter(model == 'A') %>%
  ggplot(data = . , aes(x = p_val_A, group = paste(as.numeric(Im_true), as.numeric(Re_true)),
                        color = factor(paste(as.numeric(Re_true), as.numeric(Im_true)),
                                       levels = c('0.4 0', '0 0.4', '0.4 0.4')),
                        linetype = factor(paste(as.numeric(Re_true), as.numeric(Im_true)),
                                          levels = c('0.4 0', '0 0.4', '0.4 0.4')))) +
  stat_ecdf(geom = 'step') +
  geom_abline(slope = 1, intercept = 0) +
  labs(color = expression('True '*sigma[12]), 
       linetype = expression('True '*sigma[12]), 
       x = 'Value', 
       y = 'Empirical CDF\nof p-values') + 
  scale_color_discrete(labels = expression(0.4, 0.0 + 0.4*i,0.4 + 0.4*i)) +
  scale_linetype_discrete(labels = expression(0.4, 0.0 + 0.4*i,0.4 + 0.4*i)) +
  theme_bw() +theme(text = element_text(size = 16))
ggsave(filename = 'images/ecdf_lrt.png', height = 3, width = 6)

df_null <- df %>%
  group_by(model, Re_true, Im_true) %>%
  mutate(p_val_A = pchisq(-2 * (A_Re_ll - A_Co_ll), df = 1, lower.tail = FALSE),
         p_val_B = pchisq(-2 * (B_Re_ll - B_Co_ll), df = 1, lower.tail = FALSE),
         p_val_C = pchisq(-2 * (C_Re_ll - C_Co_ll), df = 1, lower.tail = FALSE), 
         reject_A = p_val_C < 0.05, 
         reject_B = p_val_C < 0.05, 
         reject_C = p_val_C < 0.05) %>%
  filter(model == 'A', Re_true == 0.4, Im_true == 0) 
ks.test(df_null$p_val_A, y = 'punif')



df %>%
  group_by(model, Re_true, Im_true) %>%
  summarize_at(vars(contains('Sigma')), function(x) round(mean(x), 2)) %>%
  arrange(model, as.numeric(Im_true), as.numeric(Re_true)) %>%
  filter(model == 'A')

df %>%
  group_by(model, Re_true, Im_true) %>%
  summarize_at(vars(contains('Sigma')), function(x) round(sd(x), 2)) %>%
  arrange(model, as.numeric(Im_true), as.numeric(Re_true)) %>%
  filter(model == 'A')


df %>%
  group_by(model, Re_true, Im_true) %>%
  summarize_at(vars(contains('resp')), function(x) round(mean(x), 3)) %>%
  arrange(model, as.numeric(Im_true), as.numeric(Re_true)) %>%
  ungroup() %>%
  filter(model == 'A') %>%
  tidyr::pivot_longer(cols = contains('resp')) %>%
  tidyr::separate(name, into = c('model_est', 'type', 'resp', 'pred')) %>%
  dplyr::filter(model_est == 'A') %>%
  dplyr::mutate(model_summary = paste(Re_true, Im_true)) %>%
  dplyr::select(-Re_true, -Im_true, -model, -model_est) %>%
  tidyr::pivot_wider(names_from = model_summary, values_from = value) %>%
  dplyr::mutate(pred = factor(pred, levels = c('pred12', 'pred1', 'pred2')),
                type = factor(type, levels = c('Co', 'Re'))) %>% 
  arrange(resp, pred, type) %>%
  as.data.frame() 


df_pred <- df %>%
  group_by(model, Re_true, Im_true) %>%
  summarize_at(vars(contains('resp')), function(x) round(mean(x), 3)) %>%
  arrange(model, as.numeric(Im_true), as.numeric(Re_true)) %>%
  ungroup() %>%
  filter(model == 'C') %>%
  tidyr::pivot_longer(cols = contains('resp')) %>%
  tidyr::separate(name, into = c('model_est', 'type', 'resp', 'pred')) %>%
  dplyr::filter(model_est == 'C') %>%
  dplyr::mutate(model_summary = paste(Re_true, Im_true)) %>%
  dplyr::select(-Re_true, -Im_true, -model, -model_est) %>%
  tidyr::pivot_wider(names_from = model_summary, values_from = value) %>%
  dplyr::mutate(pred = factor(pred, levels = c('pred12', 'pred1', 'pred2')),
                type = factor(type, levels = c('Co', 'Re'))) %>% 
  arrange(resp, pred, type) %>%
  as.data.frame() 

df_pred

library(ggplot2)
library(tidyr)
df %>%
  select(-ends_with('ll'), -contains('Sigma')) %>%
  pivot_longer(contains('resp')) %>%
  tidyr::separate(name, into = c('model_est', 'type', 'resp', 'pred')) %>%
  dplyr::mutate(model_summary = paste(Re_true, Im_true)) %>%
  filter(Re_true == 0, Im_true == 0.4) %>%
  ggplot() +
  geom_boxplot(aes(x = c(type), group = type, y = value)) +
  facet_grid(resp~pred)
