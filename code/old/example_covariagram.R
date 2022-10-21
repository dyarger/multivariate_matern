
# data motivation
library(tidyverse)
mi_data <- read.delim('https://www.ncei.noaa.gov/pub/data/uscrn/products/hourly02/2021/CRNH0203-2021-KS_Oakley_19_SSW.txt',
                 sep = '', header = F)
# il_data <- read.delim('https://www.ncei.noaa.gov/pub/data/uscrn/products/hourly02/2021/CRNH0203-2021-IL_Shabbona_5_NNE.txt',
#                       sep = '', header = F)
# il_data <- read.delim('https://www.ncei.noaa.gov/pub/data/uscrn/products/hourly02/2021/CRNH0203-2021-VA_Charlottesville_2_SSE.txt',
#                       sep = '', header = F)
il_data <- read.delim('https://www.ncei.noaa.gov/pub/data/uscrn/products/hourly02/2021/CRNH0203-2021-KS_Oakley_19_SSW.txt',
                      sep = '', header = F)
colnames(mi_data) <- colnames(il_data) <- 
  c('WBANNO', 'UTC_DATE', 'UTC_TIME', 'LST_DATE', 'LST_TIME', 'CRX_VN', 'LONGITUDE', 'LATITUDE', 'T_CALC', 'T_HR_AVG', 'T_MAX', 'T_MIN', 'P_CALC', 'SOLARAD', 'SOLARAD_FLAG', 'SOLARAD_MAX', 'SOLARAD_MAX_FLAG', 
    'SOLARAD_MIN', 'SOLARAD_MIN_FLAG', 'SUR_TEMP_TYPE', 'SUR_TEMP', 'SUR_TEMP_FLAG', 'SUR_TEMP_MAX', 'SUR_TEMP_MAX_FLAG', 'SUR_TEMP_MIN', 'SUR_TEMP_MIN_FLAG', 'RH_HR_AVG', 'RH_HR_AVG_FLAG', 'SOIL_MOISTURE_5', 
    'SOIL_MOISTURE_10', 'SOIL_MOISTURE_20', 'SOIL_MOISTURE_50', 'SOIL_MOISTURE_100', 
    'SOIL_TEMP_5', 'SOIL_TEMP_10', 'SOIL_TEMP_20', 'SOIL_TEMP_50', 'SOIL_TEMP_100')
head(mi_data)


mi_data <- mi_data %>%
  mutate(year = substr(UTC_DATE, 1, 4),
         month = substr(UTC_DATE, 5, 6),
         day = substr(UTC_DATE, 7, 8),
         date = as.Date(paste0(year, '-', month, '-', day), 
                        format = '%Y-%m-%d'),
         date_julian = as.numeric(date),
         date_time = date_julian + as.numeric(UTC_TIME)/2400)%>%
  group_by(date) %>%
  summarize(T_HR_AVG = max(T_HR_AVG, na.rm = T))
il_data <- il_data %>%
  mutate(year = substr(UTC_DATE, 1, 4),
         month = substr(UTC_DATE, 5, 6),
         day = substr(UTC_DATE, 7, 8),
         date = as.Date(paste0(year, '-', month, '-', day), 
                        format = '%Y-%m-%d'),
         date_julian = as.numeric(date),
         date_time = date_julian + as.numeric(UTC_TIME)/2400) %>%
  group_by(date) %>%
  summarize(T_HR_AVG = max(SOIL_MOISTURE_10, na.rm = T))
data_comb <- left_join(dplyr::select(mi_data, date, T_HR_AVG) %>%
                         rename(MI = T_HR_AVG), 
                       dplyr::select(il_data, date, T_HR_AVG) %>%
                         rename(IL = T_HR_AVG)) %>%
  filter(IL > -100, MI >  - 100)

ggplot(data = data_comb, aes(x = date, y = MI)) + 
  geom_line()
ggplot(data = data_comb, aes(x = date, y = IL)) + 
  geom_line()

X_val <- cbind(1, sin((1:366) / 366* 2*pi), cos((1:366) / 366* 2*pi))
mean_fun_MI <- solve(crossprod(X_val), crossprod(X_val, data_comb$MI))
mean_fun_IL <- solve(crossprod(X_val), crossprod(X_val, data_comb$IL))
grid <- expand.grid('V1' = 1:nrow(data_comb),
                    'V2' = 1:nrow(data_comb))
MI_val <- data.frame('V' = 1:nrow(data_comb), 
                     'MI' = data_comb$MI - X_val %*% mean_fun_MI)
IL_val <- data.frame('V' = 1:nrow(data_comb), 
                     'IL' = data_comb$IL - X_val %*% mean_fun_IL)

grid_MI <-grid %>%
  left_join(rename(MI_val, V1 = V), by = 'V1') %>%
  left_join(rename(MI_val, V2 = V), by = 'V2') %>%
  mutate(prod = MI.x * MI.y,
         lag_val = (V1 - V2)) %>%
  group_by(lag_val) %>%
  summarize(cov = mean(prod))

ggplot(data = grid_MI, aes(x = lag_val, y = cov)) + 
  geom_point()

grid_IL <-grid %>%
  left_join(rename(IL_val, V1 = V), by = 'V1') %>%
  left_join(rename(IL_val, V2 = V), by = 'V2') %>%
  mutate(prod = IL.x * IL.y,
         lag_val = (V1 - V2)) %>%
  group_by(lag_val) %>%
  summarize(cov = mean(prod))
ggplot(data = grid_MI, aes(x = lag_val, y = cov)) + 
  geom_point()

grid_IL_MI <-grid %>%
  left_join(rename(IL_val, V1 = V), by = 'V1') %>%
  left_join(rename(MI_val, V2 = V), by = 'V2') %>%
  mutate(prod = IL * MI,
         lag_val = (V1 - V2)) %>%
  group_by(lag_val) %>%
  summarize(cov = mean(prod))
ggplot(data = grid_IL_MI, aes(x = lag_val, y = cov)) + 
  geom_point() + 
  scale_x_continuous(limits = c(-200, 200))

df_sum <- data.frame(grid_IL_MI, 'VA_cov' = grid_IL[['cov']],
                     'KS_cov' = grid_MI[['cov']]) %>%
  rename(VA_KS_cov = cov) %>%
  pivot_longer(cols = ends_with('cov'))
label_df <- data.frame('name'  =  c('VA_cov', 'KS_cov', 'VA_KS_cov'),
                       'label' = factor(c('Charlottesville, VA covariogram',
                                   'Oakley, KS covariogram', 
                                   'VA and KS cross-covariogram'),
                                   c('Charlottesville, VA covariogram',
                                     'VA and KS cross-covariogram',
                                     'Oakley, KS covariogram')))



ggplot(data = df_sum %>% filter(lag_val > -21, lag_val < 21) %>%
         left_join(label_df), aes(x = lag_val, y = value)) + 
  geom_point() + 
  geom_line() + 
  facet_wrap(~label, ncol = 2, scales = 'free_y') + 
  labs(x = 'Lag (number of days)',
       y = 'Covariances and cross-covariances')



  