library(devtools)

devtools::install_bitbucket("davidbolin/ngme",ref="default")

library(R.matlab)
library(ggplot2)

weather <- R.matlab::readMat('bolin_code/article_code/Application/TempPress/weather_data.mat')

weather <- as.data.frame(weather)
colnames(weather) <- c('long', 'lat', 'pres', 'temp')

ggplot(data = weather, aes(x = long, y = lat, color = temp))+
  geom_point()
ggplot(data = weather, aes(x = long, y = lat, color = pres))+
  geom_point()

dist <- fields::rdist.earth(x1 = weather[, c('long', 'lat')],
                            x2 = weather[, c('long', 'lat')])
dim(dist)
cps <- weather$pres %*% t(weather$pres)
dist_df <- data.frame(dist = as.vector(dist), 
                      cp  = as.vector(cps))
ggplot(data=dist_df,aes(dist, cp))+geom_point() + geom_smooth()+
  coord_cartesian(ylim = c(-35000, 35000))
cps <- weather$temp %*% t(weather$temp)
dist_df <- data.frame(dist = as.vector(dist), 
                      cp  = as.vector(cps))
ggplot(data=dist_df,aes(dist, cp))+geom_point() + geom_smooth()+
  coord_cartesian(ylim = c(-35000, 35000))

dist_long <- fields::rdist.earth(x1 = weather[, c('long', 'lat')],
                                 x2 = cbind(weather[, c('long')], mean(weather$lat)))
dist_long2 <- matrix(nrow = nrow(weather),
                    ncol = nrow(weather), weather[, c('lat')]) > 
  matrix(nrow = nrow(weather),
         ncol = nrow(weather), weather[, c('lat')], byrow = T)
dist_long_final <- dist_long
dist_long_final[dist_long2] <- -dist_long_final[dist_long2]

dist_lat <- fields::rdist.earth(x1 = weather[, c('long', 'lat')],
                                 x2 = cbind(mean(weather$long), weather[, c('lat')]))
dist_lat2 <- matrix(nrow = nrow(weather),
                     ncol = nrow(weather), weather[, c('long')]) > 
  matrix(nrow = nrow(weather),
         ncol = nrow(weather), weather[, c('long')], byrow = T)
dist_lat_final <- dist_lat
dist_lat_final[dist_lat2] <- -dist_lat_final[dist_lat2]



cps <- weather$pres %*% t(weather$temp)
dist_df <- data.frame(dist = as.vector(dist_lat_final), 
                      cp  = as.vector(cps))
dist_df
ggplot(data=dist_df,aes(dist, cp))+geom_point() + geom_smooth()+
  coord_cartesian(xlim = c(-100, 100))


