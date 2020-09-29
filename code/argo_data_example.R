

# argo data example

load('data/region_residuals.RData')

head(df_region_resid)

library(dplyr)
df <- df_region_resid %>%
  filter(!is.na(oxy_resids) & !is.na(pressure) & temp_QC != 8 &
           oxy_QC != 8 & psal_QC!=8 & pressure_QC != 8,
         pressure > 280 & pressure < 320)
dim(df)
library(ggplot2)
ggplot(data = df, aes(x = longitude, y = latitude, color = oxy_resids))+
  geom_point() + 
  scale_color_gradient2()

ggplot(data = df, aes(x = longitude, y = latitude, color = temp_resids))+
  geom_point() + 
  scale_color_gradient2()

ggplot(data = df, aes(x = longitude, y = latitude, color = psal_resids))+
  geom_point() + 
  scale_color_gradient2()

cor(df[, c('oxy_resids', 'temp_resids', 'psal_resids')])
library(fields)
distances <- rdist.earth(df[,c('longitude', 'latitude')], miles = F)
distances_lat <- rdist.earth(df[,c('longitude', 'latitude')],
                             cbind(mean(df$longitude), df[,'latitude']), miles = F)
distances_lon <- rdist.earth(df[,c('longitude', 'latitude')],
                             cbind(df[,'longitude'], mean(df$latitude)), miles = F)
lat_sign <- (rep(1, nrow(df)) %*% t(df$latitude)) < (df$latitude %*%t(rep(1, nrow(df))))
lon_sign <- (rep(1, nrow(df)) %*% t(df$longitude)) < (df$longitude %*%t(rep(1, nrow(df))))
distances_lat <- ifelse(as.vector(lat_sign), as.vector(distances_lat),-as.vector(distances_lat))
distances_lon <- ifelse(as.vector(lon_sign), as.vector(distances_lon),-as.vector(distances_lon))


TO <- df$oxy_resids %*% t(df$psal_resids)

# dist_test <- data.frame(dist = as.vector(distances),
#                         dist_lon = as.vector(distances_lon),
#                         dist_lat = as.vector(distances_lat),
#                         TO = as.vector(TO)) %>%
#   mutate(dist_group = round(dist/25) * 25) %>%
#   group_by(dist_group) %>%
#   summarise(TO = mean(TO))

dist_test <- data.frame(dist = as.vector(distances),
                        dist_lon = as.vector(distances_lon),
                        dist_lat = as.vector(distances_lat),
                        TO = as.vector(TO)) %>%
  mutate(dist_group = round(dist_lon/50) * 50) %>%
  group_by(dist_group) %>%
  summarise(TO = mean(TO))


ggplot(data = dist_test, aes(x = dist_group, y = TO))+
  geom_point()+
  coord_cartesian(xlim = c(-900, 900))

ggplot(data = dist_test, aes(x = dist_lon, y = TO))+
  geom_point()

