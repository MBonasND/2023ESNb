#####################################################
#####################################################
### Bonas, Castruccio - 2021ESNDiscussion Figures ###
#####################################################
#####################################################

#clear enviroment and load libraries
rm(list = ls())
library(tidyverse)
library(ggpubr)
library(tigris)
library(ggsn)
library(spcov)
library(lubridate)
library(geoR)
library(fields)
library(zoo)
library(reshape2)
library(scales)
library(verification)
library(forecast)
library(pracma)
library(TTR)
library(fda)


#Declare borders
elnino_border = expand.grid(c(-120, -170), c(5, -5))
swap = c(1, 2, 4, 3)
elnino_border =  elnino_border[swap, ]


#load functions
source('functions.R')

area_grid = function(lat, lon)
{
  temp = meshgrid(lon, lat)
  R = earth_radius(temp$Y)
  
  dlat = deg2rad(apply(temp$Y, 2, gradient))
  dlon = deg2rad(apply(temp$X, 1, gradient))
  
  dy = dlat * R
  dx = t(dlon) * R * cos(deg2rad(temp$Y))
  
  area = dy*dx
  return(area)
}

earth_radius = function(lat)
{
  a = 6378137
  b = 6356752.3142
  e2 = 1 - (b^2/a^2)
  
  lat = deg2rad(lat)
  lat_gc = atan((1-e2)*tan(lat))
  
  r = ((a*(1-e2)^0.5) / (1-(e2*cos(lat_gc)^2))^0.5)
  return(r)
}


################
### Figure 1 ###
################

load('merra2_enso_grid.RData')
load('MERRA2_Enso_Anomaly.RData')

#Plots for Data Section
#Extract only Nino-3.4 Region
full_nino_anomaly = cbind(merra2_enso_grid, enso.anomaly)

inside = locations.inside(merra2_enso_grid, elnino_border, as.is = FALSE)
colnames(inside) = c('lon', 'lat')
colnames(full_nino_anomaly) = c('lon', 'lat')

nino34 = data.frame(inside) %>% left_join(data.frame(full_nino_anomaly), by = c('lon', 'lat'))

nino34_index = apply(nino34[,-c(1,2)], 2, mean)


#plot jan 1989 --- coolest time = 97
cols <- c("navy", "blue", "cyan", "green", "yellow", "red", "red4")
ln_gg = ggplot() +
  geom_point(mapping = aes(x = merra2_enso_grid[,1],
                           y = merra2_enso_grid[,2],
                           color = enso.anomaly[,97]),
             shape = 15,
             size = 2) +
  scale_color_gradientn(colours = cols, 
                        limits = c(-5.5,5.5)) +
  geom_point(mapping = aes(x = merra2_enso_grid[!complete.cases(enso.anomaly),1],
                           y = merra2_enso_grid[!complete.cases(enso.anomaly),2]),
             shape = 15,
             size = 2,
             color = 'black') +
  geom_rect(mapping = aes(xmin = -170,
                          xmax = -120,
                          ymin = -5, 
                          ymax = 5),
            color = 'black',
            alpha = 0,
            lwd = 1) +
  theme(plot.title = element_text(hjust = 0.5)) +
  labs(title = 'January 1989', x = '', y = '', color = 'Temperature Anomaly (°C)') + 
  theme(plot.title = element_text(hjust = 0.5, size = 16), 
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 12),
        legend.title = element_text(size = 14), 
        legend.text = element_text(size = 11),
        legend.position = 'bottom') + 
  scale_y_continuous(breaks = seq(-30, 30, 10), expand = c(0,0)) +
  scale_x_continuous(breaks = seq(-180, -60, 20), expand = c(0,0)) +
  theme(legend.key.width = unit(1.5,"cm")) +
  guides(color = guide_colorbar(title.vjust = 1))

#plot nov 2015 --- warmest time = 419
cols <- c("navy", "blue", "cyan", "green", "yellow", "red", "red4")
en_gg = ggplot() +
  geom_point(mapping = aes(x = merra2_enso_grid[,1],
                           y = merra2_enso_grid[,2],
                           color = enso.anomaly[,419]),
             shape = 15,
             size = 2) +
  scale_color_gradientn(colours = cols, 
                        limits = c(-5.5,5.5)) +
  geom_point(mapping = aes(x = merra2_enso_grid[!complete.cases(enso.anomaly),1],
                           y = merra2_enso_grid[!complete.cases(enso.anomaly),2]),
             shape = 15,
             size = 2,
             color = 'black') +
  geom_rect(mapping = aes(xmin = -170,
                          xmax = -120,
                          ymin = -5, 
                          ymax = 5),
            color = 'black',
            alpha = 0,
            lwd= 1) +
  labs(title = 'November 2015', x = '', y = '', color = 'Temperature Anomaly (°C)') + 
  theme(plot.title = element_text(hjust = 0.5, size = 16), 
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 12),
        legend.title = element_text(size = 14), 
        legend.text = element_text(size = 11),
        legend.position = 'bottom') + 
  scale_y_continuous(breaks = seq(-30, 30, 10), expand = c(0,0)) +
  scale_x_continuous(breaks = seq(-180, -60, 20), expand = c(0,0)) +
  theme(legend.key.width = unit(1.5,"cm")) +
  guides(color = guide_colorbar(title.vjust = 1))

#Nino Index Plot
dates = seq(as.Date("1981/1/16"), by = "month", length.out = 480)
nino.gg = ggplot() +
  geom_line(mapping = aes(x = as.yearmon(dates), 
                          y = nino34_index), 
            lwd = 1, 
            lty = 1, 
            color = 'black') +
  geom_hline(yintercept = c(0.5,-0.5),
             color = 'gray50', 
             lwd = 1,
             lty = 2) +
  geom_vline(xintercept = as.yearmon(dates[c(97,419)]) ,
             color = 'indianred1',
             lwd = 1,
             lty =1, 
             alpha = 0.4) +
  geom_point(mapping = aes(x = as.yearmon(dates[c(97,419)])),
             y = nino34_index[c(97,419)],
             color = 'indianred1', 
             shape = 19, 
             size = 2) +
  theme_bw() +
  labs(title = 'Niño3.4 Index', x = '', y = 'Temp. Anomaly (°C)') + 
  theme(plot.title = element_text(hjust = 0.5, size = 16), 
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 12)) + 
  scale_y_continuous(breaks = seq(-2, 3, 1),
                     limits = c(-2,3.25), 
                     sec.axis = sec_axis( trans=~., name='Temp. Anomaly (°C)')) +
  scale_x_continuous(breaks = seq(1980, 2021, 5),
                     limits = c(as.yearmon(dates[1]), as.yearmon(dates[480])))


enln.gg = ggarrange(en_gg, ln_gg,
                    nrow = 1, ncol = 2,
                    labels = c('A', 'B'),
                    common.legend = TRUE,
                    legend = 'bottom')


panel.gg = ggarrange(enln.gg, nino.gg,
                     nrow = 2, ncol = 1,
                     heights = c(1.75,1), 
                     labels = c('', 'C'))




################
### Figure 2 ###
################

load('Full_CMIP6_ENSO_EOF.RData')
load('Full_CMIP6_ENSO_EOF_Betas.RData')
load('cmip6_enso_grid.RData')
load('C6_ENSO_NH_R134.RData')
load('QDist_ENSO_CMIP6.RData')
load('1850_CMIP6_ENSO_Anomaly.RData')
load('1950_CMIP6_ENSO_Anomaly.RData')
load('CMIP6_ENSO_Anomaly.RData')
full_cmip6_anomaly = cbind(cmip6_anomaly_1850, cmip6_anomaly_1950, cmip6_anomaly)

rawData = full_cmipenso_eof
index = 31
tau = 1
trainLen = 1739 + (index-1)*36 #total = 2819
forward = 36


#Generate training/testing/valid sets
sets = cttv(rawData = rawData,
            tau = tau,
            trainLen = trainLen,
            testLen = forward,
            valid.flag = F)


upper = mean.pred + quant_dist[,2,]
lower = mean.pred - quant_dist[,1,]


#load('C6_ENSO_Windows31.RData')
#sd.pred = sapply(1:20, function(x) apply(ensemb.pred[x,,], 2, sd))
#upper = mean.pred + qnorm(0.975)*sd.pred
#lower = mean.pred - qnorm(0.975)*sd.pred

full_grid_forcs = full_cmipenso_eof_betas %*% t(mean.pred)
full_grid = full_cmipenso_eof_betas %*% t(sets$yTest)

full_grid_forcs_upper = full_cmipenso_eof_betas %*% t(upper)
full_grid_forcs_lower = full_cmipenso_eof_betas %*% t(lower)

#fix upper and lower bounds
grid_lower = full_grid_forcs_lower
grid_upper = full_grid_forcs_upper
filt = full_grid_forcs_upper < full_grid_forcs_lower

grid_lower[filt] = full_grid_forcs_upper[filt]
grid_upper[filt] = full_grid_forcs_lower[filt]


filt2 = grid_upper < full_grid_forcs
grid_upper[filt2] = full_grid_forcs[filt2]

filt3 = grid_lower > full_grid_forcs
grid_lower[filt3] = full_grid_forcs[filt3]

borderfix = cbind(rep(360, 7483), rep(0, 7483))
nino3.4_forcs = cbind(cmip6_enso_grid[complete.cases(full_cmip6_anomaly),]-borderfix, full_grid_forcs)
nino3.4_actual = cbind(cmip6_enso_grid[complete.cases(full_cmip6_anomaly),]-borderfix, full_grid)
nino3.4_upper = cbind(cmip6_enso_grid[complete.cases(full_cmip6_anomaly),]-borderfix, grid_upper)
nino3.4_lower = cbind(cmip6_enso_grid[complete.cases(full_cmip6_anomaly),]-borderfix, grid_lower)

inside = locations.inside(cmip6_enso_grid[complete.cases(cmip6_anomaly),]-borderfix, elnino_border, as.is = FALSE)
colnames(inside) = c('lon', 'lat')
colnames(nino3.4_forcs[,1:2]) = c('lon', 'lat')
colnames(nino3.4_actual[,1:2]) = c('lon', 'lat')
colnames(nino3.4_upper[,1:2]) = c('lon', 'lat')
colnames(nino3.4_lower[,1:2]) = c('lon', 'lat')

nino34_f = data.frame(inside) %>% left_join(data.frame(nino3.4_forcs), by = c('lon', 'lat'))
nino34_a = data.frame(inside) %>% left_join(data.frame(nino3.4_actual), by = c('lon', 'lat'))
nino34_u = data.frame(inside) %>% left_join(data.frame(nino3.4_upper), by = c('lon', 'lat'))
nino34_l = data.frame(inside) %>% left_join(data.frame(nino3.4_lower), by = c('lon', 'lat'))

mean((nino34_a[,-c(1,2)] <= nino34_u[,-c(1,2)] & nino34_a[,-c(1,2)] >= nino34_l[,-c(1,2)]))



forc.avg = apply(nino34_f[,-c(1,2)], 2, mean)
actual.avg = apply(nino34_a[,-c(1,2)], 2, mean)
upper.avg = apply(nino34_u[,-c(1,2)], 2, mean)
lower.avg = apply(nino34_l[,-c(1,2)], 2, mean)



dates_cmip6 = seq(as.Date("1850/1/16"), by = "month", length.out = 3012)
forc.dates = as.POSIXct(dates_cmip6[2821:2856])

#Panel A
eof_num = 1
upper = mean.pred[,eof_num] + (quant_dist[,2,eof_num])
lower = mean.pred[,eof_num] - (quant_dist[,1,eof_num])



eof_gg = ggplot() +
  geom_line(mapping = aes(x = forc.dates,
                          y = sets$yTest[,eof_num],
                          color = 'black'),
            lty = 1,
            lwd = 1.5) +
  geom_line(mapping = aes(x = forc.dates, 
                          y = upper,
                          color = 'red'),
            lty = 2,
            lwd = 1.5) +
  geom_line(mapping = aes(x = forc.dates,
                          y = mean.pred[,eof_num],
                          color = 'red1'), 
            lty = 1,
            lwd = 1.5) +
  geom_line(mapping = aes(x = forc.dates,
                          y = lower),
            color = 'red',
            lty = 2,
            lwd = 1.5) + 
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5, size = 16),
        axis.text = element_text(size = 12), 
        axis.title = element_text(size = 14)) +
  labs(y = '', x = 'Time', title = 'CMIP6-MIROC EOF1', color = '') +
  scale_color_manual(values = c('black', 'red', 'red1'), 
                     labels = c('True Value', 'Mean Forecasts', '95% Prediction Interval'),
                     guide = guide_legend(override.aes = list(lty = c(1,1,2)))) +
  theme(legend.position = c(0.25,0.15), legend.text = element_text(size = 14)) +
  theme(legend.key.width = unit(2.5,"cm")) +
  scale_y_continuous(breaks = seq(-1.5,1.5,0.5), limits = c(-1.7, 1.6)) +
  scale_x_datetime(date_labels = '%Y')



#eof_gg

#Panel B
nino_gg = ggplot() +
  geom_line(mapping = aes(x = forc.dates,
                          y = actual.avg,
                          color = 'black'),
            lty = 1,
            lwd = 1.5) +
  geom_line(mapping = aes(x = forc.dates, 
                          y = upper.avg,
                          color = 'red'),
            lty = 2,
            lwd = 1.5) +
  geom_line(mapping = aes(x = forc.dates,
                          y = forc.avg,
                          color = 'red1'), 
            lty = 1,
            lwd = 1.5) +
  geom_line(mapping = aes(x = forc.dates,
                          y = lower.avg),
            color = 'red',
            lty = 2,
            lwd = 1.5) + 
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5, size = 16),
        axis.text = element_text(size = 12), 
        axis.title = element_text(size = 14),
        axis.title.y = element_text(angle = 90, vjust = 0.5)) +
  labs(y = 'Temperature Anomaly (°C)', x = 'Time', title = 'CMIP6-MIROC Niño3.4 Index', color = '') +
  scale_color_manual(values = c('black', 'red','red1' ), 
                     labels = c('True Value', 'Mean Forecasts', '95% Prediction Interval'),
                     guide = guide_legend(override.aes = list(lty = c(1,1,2)))) +
  theme(legend.position = c(0.25,0.2), legend.text = element_text(size = 12)) +
  theme(legend.key.width = unit(2.5,"cm")) +
  scale_y_continuous(breaks = seq(-4,3,1), limits = c(-4, 3), position = 'right') +
  scale_x_datetime(date_labels = '%Y')


#nino_gg

ggarrange(eof_gg, nino_gg, nrow = 1, ncol = 2,
          common.legend = TRUE, legend = 'bottom',
          labels = c('A', 'B'))




################
### Figure 3 ###
################


#run all code above in Figure 2 to obtain this plot

time = 9
grid_lower[grid_lower[,time] < -5,time] = -5


cols = c("navy", "blue", "cyan", "green", "yellow", "red", "red4")
forc_gg = ggplot() +
  geom_point(mapping = aes(x = cmip6_enso_grid[complete.cases(cmip6_anomaly),1]-360,
                           y = cmip6_enso_grid[complete.cases(cmip6_anomaly),2],
                           color = full_grid_forcs[,time]),
             shape = 15,
             size = 3) +
  geom_point(mapping = aes(x = cmip6_enso_grid[!complete.cases(cmip6_anomaly),1]-360,
                           y = cmip6_enso_grid[!complete.cases(cmip6_anomaly),2]),
             color = 'black',
             shape = 15,
             size = 3) +
  geom_rect(mapping = aes(xmin = -170,
                          xmax = -120,
                          ymin = -5, 
                          ymax = 5),
            color = 'black',
            alpha = 0,
            lwd = 1) +
  labs(x = '', y = '', color = 'Temperature Anomaly (°C)', title = 'Forecasts') +
  scale_color_gradientn(colours = cols, 
                        limits = c(-5,5)) +
  theme(plot.title = element_text(hjust = 0.5, size = 16),
        axis.text = element_text(size = 12)) + 
  theme(legend.text = element_text(size = 12), legend.title = element_text(size = 14)) +
  theme(legend.key.width = unit(1.5,"cm")) + 
  scale_y_continuous(breaks = seq(-30,30,10), expand = c(0,0)) +
  scale_x_continuous(breaks = seq(-180, -60, 20), expand = c(0,0))

#Actual
actual_gg = ggplot() +
  geom_point(mapping = aes(x = cmip6_enso_grid[complete.cases(cmip6_anomaly),1]-360,
                           y = cmip6_enso_grid[complete.cases(cmip6_anomaly),2],
                           color = full_grid[,time]),
             shape = 15,
             size = 3) +
  geom_point(mapping = aes(x = cmip6_enso_grid[!complete.cases(cmip6_anomaly),1]-360,
                           y = cmip6_enso_grid[!complete.cases(cmip6_anomaly),2]),
             color = 'black',
             shape = 15,
             size = 3) +
  geom_rect(mapping = aes(xmin = -170,
                          xmax = -120,
                          ymin = -5, 
                          ymax = 5),
            color = 'black',
            alpha = 0,
            lwd = 1) +
  labs(x = '', y = '', color = 'Anomaly °C', title = 'Actual') +
  scale_color_gradientn(colours = cols, 
                        limits = c(-5,5)) +
  theme(plot.title = element_text(hjust = 0.5, size = 16),
        axis.text = element_text(size = 12)) + 
  scale_y_continuous(breaks = seq(-30,30,10), expand = c(0,0)) +
  scale_x_continuous(breaks = seq(-180, -60, 20), expand = c(0,0))


#Forecasts Upper
upper_gg = ggplot() +
  geom_point(mapping = aes(x = cmip6_enso_grid[complete.cases(cmip6_anomaly),1]-360,
                           y = cmip6_enso_grid[complete.cases(cmip6_anomaly),2],
                           color = grid_upper[,time]),
             shape = 15,
             size = 3) +
  geom_point(mapping = aes(x = cmip6_enso_grid[!complete.cases(cmip6_anomaly),1]-360,
                           y = cmip6_enso_grid[!complete.cases(cmip6_anomaly),2]),
             color = 'black',
             shape = 15,
             size = 3) +
  geom_rect(mapping = aes(xmin = -170,
                          xmax = -120,
                          ymin = -5, 
                          ymax = 5),
            color = 'black',
            alpha = 0,
            lwd = 1) +
  labs(x = '', y = '', color = 'Anomaly °C', title = 'Upper 95% PI') +
  scale_color_gradientn(colours = cols, 
                        limits = c(-5,5)) +
  theme(plot.title = element_text(hjust = 0.5, size = 16),
        axis.text = element_text(size = 12)) + 
  scale_y_continuous(breaks = seq(-30,30,10), expand = c(0,0)) +
  scale_x_continuous(breaks = seq(-180, -60, 20), expand = c(0,0))

#Forecasts Lower
lower_gg = ggplot() +
  geom_point(mapping = aes(x = cmip6_enso_grid[complete.cases(cmip6_anomaly),1]-360,
                           y = cmip6_enso_grid[complete.cases(cmip6_anomaly),2],
                           color = grid_lower[,time]),
             shape = 15,
             size = 3) +
  geom_point(mapping = aes(x = cmip6_enso_grid[!complete.cases(cmip6_anomaly),1]-360,
                           y = cmip6_enso_grid[!complete.cases(cmip6_anomaly),2]),
             color = 'black',
             shape = 15,
             size = 3) +
  geom_rect(mapping = aes(xmin = -170,
                          xmax = -120,
                          ymin = -5, 
                          ymax = 5),
            color = 'black',
            alpha = 0,
            lwd = 1) +
  labs(x = '', y = '', color = 'Anomaly °C', title = 'Lower 95% PI') +
  scale_color_gradientn(colours = cols, 
                        limits = c(-5,5)) +
  theme(plot.title = element_text(hjust = 0.5, size = 16),
        axis.text = element_text(size = 12)) + 
  scale_y_continuous(breaks = seq(-30,30,10), expand = c(0,0)) +
  scale_x_continuous(breaks = seq(-180, -60, 20), expand = c(0,0))


ggarrange(forc_gg, actual_gg, upper_gg, lower_gg, nrow=2, ncol = 2, common.legend = TRUE, labels = c('A', 'B','C','D'), legend = 'bottom')







mean(full_grid >= grid_lower & full_grid <= grid_upper)

contained = (full_grid[,time] >= grid_lower[,time] & full_grid[,time] <= grid_upper[,time])

contained_gg = ggplot() +
  geom_point(mapping = aes(x = cmip6_enso_grid[complete.cases(cmip6_anomaly),1]-360,
                           y = cmip6_enso_grid[complete.cases(cmip6_anomaly),2],
                           color = as.factor(contained)),
             shape = 15,
             size = 3.1) +
  geom_point(mapping = aes(x = cmip6_enso_grid[!complete.cases(cmip6_anomaly),1]-360,
                           y = cmip6_enso_grid[!complete.cases(cmip6_anomaly),2]),
             color = 'black',
             shape = 15,
             size = 3) +
  geom_rect(mapping = aes(xmin = -170,
                          xmax = -120,
                          ymin = -5, 
                          ymax = 5),
            color = 'black',
            alpha = 0,
            lwd = 1) +
  labs(x = '', y = '', color = 'Anomaly °C', title = 'Within 95% PI?') +
  scale_color_manual(values = c('gray90', 'darkgreen'))+
  theme(plot.title = element_text(hjust = 0.5, size = 16),
        axis.text = element_text(size = 12)) + 
  scale_y_continuous(breaks = seq(-30,30,10), expand = c(0,0)) +
  scale_x_continuous(breaks = seq(-180, -60, 20), expand = c(0,0))

contained_gg


################
### Figure 4 ###
################

load('merra2_enso_grid.RData')
load('MERRA2_Enso_Anomaly.RData')
load('M2_ESN_TL_Future.RData')
load('MERRA2_ENSO_EOF_Betas.RData')
load('MERRA2_ENSO_EOF.RData')
load('merra2_enso_grid.RData')
load('Full_TL_QDist.RData')

tl.pred = sapply(1:20, function(x) colMeans(ensemb.pred[x,,]))

full_grid_tl = m2enso_eof_betas %*% t(tl.pred)

nino3.4_tl = cbind(merra2_enso_grid[complete.cases(enso.anomaly),], full_grid_tl)
inside = locations.inside(merra2_enso_grid[complete.cases(cmip6_anomaly),], elnino_border, as.is = FALSE)
colnames(inside) = c('Var1', 'Var2')

nino34_tl = data.frame(inside) %>% left_join(data.frame(nino3.4_tl), by = c('Var1', 'Var2'))
tl.avg = apply(nino34_tl[,-c(1,2)], 2, mean)



trainLen = 478 #adjust
validLen = 36 #adjust
testLen = 1 #adjust
forward = 1
tau  = 1
rawData = m2enso_eof
#Generate training/testing/valid sets
sets = cttv(rawData = rawData,
            tau = tau,
            trainLen = trainLen,
            testLen = forward,
            valid.flag = F)


full_grid = m2enso_eof_betas %*% t(sets$yTrain)
full_grid_forcs_upper = m2enso_eof_betas %*% t(tl.pred + full_tl_quant_dist[,2,])
full_grid_forcs_lower = m2enso_eof_betas %*% t(tl.pred - full_tl_quant_dist[,1,])

#fix upper and lower bounds
grid_lower = full_grid_forcs_lower
grid_upper = full_grid_forcs_upper
filt = full_grid_forcs_upper < full_grid_forcs_lower

grid_lower[filt] = full_grid_forcs_upper[filt]
grid_upper[filt] = full_grid_forcs_lower[filt]


filt2 = grid_upper < full_grid_tl
grid_upper[filt2] = full_grid_tl[filt2]

filt3 = grid_lower > full_grid_tl
grid_lower[filt3] = full_grid_tl[filt3]

nino3.4_actual = cbind(merra2_enso_grid[complete.cases(enso.anomaly),], full_grid)
nino3.4_upper = cbind(merra2_enso_grid[complete.cases(enso.anomaly),], grid_upper)
nino3.4_lower = cbind(merra2_enso_grid[complete.cases(enso.anomaly),], grid_lower)

inside = locations.inside(merra2_enso_grid[complete.cases(enso.anomaly),], elnino_border, as.is = FALSE)
colnames(inside) = c('Var1', 'Var2')

nino34_a = data.frame(inside) %>% left_join(data.frame(nino3.4_actual), by = c('Var1', 'Var2'))
nino34_u = data.frame(inside) %>% left_join(data.frame(nino3.4_upper), by = c('Var1', 'Var2'))
nino34_l = data.frame(inside) %>% left_join(data.frame(nino3.4_lower), by = c('Var1', 'Var2'))


actual.avg = apply(nino34_a[,-c(1,2)], 2, mean)
upper.avg = apply(nino34_u[,-c(1,2)], 2, mean)
lower.avg = apply(nino34_l[,-c(1,2)], 2, mean)


actual.avg = apply(nino34_a[,-c(1,2)], 2, mean)
dates = seq(as.Date("1981/1/16"), by = "month", length.out = 480)
old.dates = tail(dates, 120)
old.nino = tail(actual.avg, 120)


new.dates = tail(seq(as.Date("1981/1/16"), by = "month", length.out = (480+36)), 36)
new.nino = tl.avg


#Transfer Learning Plot
ggplot() +
  geom_line(mapping = aes(x = as.yearmon(old.dates), y = old.nino),
            color = 'black', 
            lwd = 2) +
  geom_line(mapping = aes(x = as.yearmon(new.dates), y = new.nino),
            color = 'red',
            lwd = 2) + 
  geom_line(mapping = aes(x = as.yearmon(new.dates), y = upper.avg),
            color = 'red',
            lwd = 2, lty = 2) +
  geom_line(mapping = aes(x = as.yearmon(new.dates), y = lower.avg),
            color = 'red',
            lwd = 2, lty = 2) +
  geom_hline(yintercept = c(-0.5,0.5),
             color = 'gray60', 
             lty = 2,
             lwd = 1.25) +
  labs(x = '', y = 'Temperature Anomaly (°C)', title = 'Future Niño3.4 Index Forecasts') +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5, size = 16), 
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 12)) + 
  scale_y_continuous(breaks = seq(-2.5, 2.5, 0.5),
                     limits = c(-2.6,2.6)) +
  scale_x_continuous(breaks = seq(2011, 2024, 2),
                     limits = c(as.yearmon(old.dates[1]), as.yearmon(new.dates[36])))



#################
### Figure S1 ###
#################

load('merra2_pdo_grid.RData')
load('MERRA2_PDO_Anomaly.RData')
file = paste0('MERRA2_PDO_EOF.RData')
load(file)
#Equivalent of Warm/Cold/Index ENSO Figure

borderfix = cbind(rep(360,19291), rep(0,19291))
m2_pdo_grid = m2_pdo_grid - borderfix

#plot jun 1997
cols <- c("navy", "blue", "cyan", "green", "yellow", "red", "red4")
pospdo_gg = ggplot() +
  geom_point(mapping = aes(x = m2_pdo_grid[,1],
                           y = m2_pdo_grid[,2],
                           color = pdo.anomaly[,198]),
             shape = 15,
             size = 2) +
  scale_color_gradientn(colours = cols, 
                        limits = c(-4,4)) +
  geom_point(mapping = aes(x = m2_pdo_grid[!complete.cases(pdo.anomaly),1],
                           y = m2_pdo_grid[!complete.cases(pdo.anomaly),2]),
             shape = 15,
             size = 2,
             color = 'black') +
  theme(plot.title = element_text(hjust = 0.5)) +
  labs(title = 'June 1997', x = '', y = '', color = 'Temperature Anomaly (°C)') + 
  theme(plot.title = element_text(hjust = 0.5, size = 16), 
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 12),
        legend.title = element_text(size = 14), 
        legend.text = element_text(size = 11),
        legend.position = 'bottom') + 
  scale_y_continuous(breaks = seq(20, 70, 10), expand = c(0,0)) +
  scale_x_continuous(breaks = seq(-240, -120, 20), expand = c(0,0)) +
  theme(legend.key.width = unit(1.5,"cm")) +
  guides(color = guide_colorbar(title.vjust = 1))

#plot nov 2011
cols <- c("navy", "blue", "cyan", "green", "yellow", "red", "red4")
negpdo_gg = ggplot() +
  geom_point(mapping = aes(x = m2_pdo_grid[,1],
                           y = m2_pdo_grid[,2],
                           color = pdo.anomaly[,371]),
             shape = 15,
             size = 2) +
  scale_color_gradientn(colours = cols, 
                        limits = c(-4,4)) +
  geom_point(mapping = aes(x = m2_pdo_grid[!complete.cases(pdo.anomaly),1],
                           y = m2_pdo_grid[!complete.cases(pdo.anomaly),2]),
             shape = 15,
             size = 2,
             color = 'black') +
  labs(title = 'November 2011', x = '', y = '', color = 'Anomaly (°C)') + 
  theme(plot.title = element_text(hjust = 0.5, size = 16), 
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 12),
        legend.title = element_text(size = 14), 
        legend.text = element_text(size = 11),
        legend.position = 'bottom') + 
  scale_y_continuous(breaks = seq(20, 70, 10), expand = c(0,0)) +
  scale_x_continuous(breaks = seq(-240, -120, 20), expand = c(0,0)) +
  theme(legend.key.width = unit(1.5,"cm")) +
  guides(color = guide_colorbar(title.vjust = 1))


#PDO Index Plot
dates = seq(as.Date("1981/1/16"), by = "month", length.out = 480)
pdoindex.gg = ggplot() +
  geom_line(mapping = aes(x = as.yearmon(dates), 
                          y = -m2pdo_eof[,1]), 
            lwd = 1, 
            lty = 1, 
            color = 'black') +
  geom_hline(yintercept = c(0),
             color = 'gray50', 
             lwd = 1,
             lty = 2) +
  geom_vline(xintercept = as.yearmon(dates[c(198,371)]) ,
             color = 'indianred1',
             lwd = 1,
             lty =1, 
             alpha = 0.4) +
  geom_point(mapping = aes(x = as.yearmon(dates[c(198,371)])),
             y = -m2pdo_eof[c(198,371),1],
             color = 'indianred1', 
             shape = 19, 
             size = 2) +
  theme_bw() +
  labs(title = 'Pacific Decadal Oscillation Index', x = '', y = 'Index') + 
  theme(plot.title = element_text(hjust = 0.5, size = 16), 
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 12)) + 
  scale_y_continuous(breaks = seq(-1, 1, 0.5),
                     limits = c(-1,1)) +
  scale_x_continuous(breaks = seq(1980, 2021, 5),
                     limits = c(as.yearmon(dates[1]), as.yearmon(dates[480])))

posneg.gg = ggarrange(pospdo_gg, negpdo_gg,
                      nrow = 1, ncol = 2,
                      labels = c('A', 'B'),
                      common.legend = TRUE,
                      legend = 'bottom')


panel.gg = ggarrange(posneg.gg, pdoindex.gg,
                     nrow = 2, ncol = 1,
                     heights = c(1.75,1), 
                     labels = c('', 'C'))




#################
### Figure S2 ###
#################


load('merra2_amo_grid.RData')
load('MERRA2_AMO_Anomaly.RData')
#Equivalent of Warm/Cold/Index ENSO Figure



#plot jun 1997
cols <- c("navy", "blue", "cyan", "green", "yellow", "red", "red4")
amo.anomaly[amo.anomaly[,141] > 3,141] = 3
amo.anomaly[amo.anomaly[,141] < -3,141] = -3
posamo_gg = ggplot() +
  geom_point(mapping = aes(x = merra2_amo_grid[,1],
                           y = merra2_amo_grid[,2],
                           color = amo.anomaly[,141]),
             shape = 15,
             size = 2) +
  scale_color_gradientn(colours = cols, 
                        limits = c(-3,3)) +
  geom_point(mapping = aes(x = merra2_amo_grid[!complete.cases(amo.anomaly),1],
                           y = merra2_amo_grid[!complete.cases(amo.anomaly),2]),
             shape = 15,
             size = 2,
             color = 'black') +
  theme(plot.title = element_text(hjust = 0.5)) +
  labs(title = 'September 1992', x = '', y = '', color = 'Temperature Anomaly (°C)') + 
  theme(plot.title = element_text(hjust = 0.5, size = 16), 
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 12),
        legend.title = element_text(size = 14), 
        legend.text = element_text(size = 11),
        legend.position = 'bottom') + 
  scale_y_continuous(breaks = seq(0, 65, 10), expand = c(0,0), limits = c(0,65)) +
  scale_x_continuous(breaks = seq(-80, 0, 10), expand = c(0,0)) +
  theme(legend.key.width = unit(1.5,"cm")) +
  guides(color = guide_colorbar(title.vjust = 1))

#plot nov 2011
cols <- c("navy", "blue", "cyan", "green", "yellow", "red", "red4")
amo.anomaly[amo.anomaly[,356] > 3,356] = 3
amo.anomaly[amo.anomaly[,356] < -3,356] = -3
negamo_gg = ggplot() +
  geom_point(mapping = aes(x = merra2_amo_grid[,1],
                           y = merra2_amo_grid[,2],
                           color = amo.anomaly[,356]),
             shape = 15,
             size = 2) +
  scale_color_gradientn(colours = cols, 
                        limits = c(-3,3)) +
  geom_point(mapping = aes(x = merra2_amo_grid[!complete.cases(amo.anomaly),1],
                           y = merra2_amo_grid[!complete.cases(amo.anomaly),2]),
             shape = 15,
             size = 2,
             color = 'black') +
  labs(title = 'August 2010', x = '', y = '', color = 'Temperature Anomaly (°C)') + 
  theme(plot.title = element_text(hjust = 0.5, size = 16), 
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 12),
        legend.title = element_text(size = 14), 
        legend.text = element_text(size = 11),
        legend.position = 'bottom') + 
  scale_y_continuous(breaks = seq(0, 65, 10), expand = c(0,0), limits = c(0,65)) +
  scale_x_continuous(breaks = seq(-80, 0, 10), expand = c(0,0)) +
  theme(legend.key.width = unit(1.5,"cm")) +
  guides(color = guide_colorbar(title.vjust = 1))


posneg.gg = ggarrange(posamo_gg, negamo_gg,
                      nrow = 1, ncol = 2,
                      labels = c('A', 'B'),
                      common.legend = TRUE,
                      legend = 'bottom')



#Calculate Area Averaged Anomalies - Derive AMO Index
#lat = lat[lat<=65] #https://climatedataguide.ucar.edu/climate-data/atlantic-multi-decadal-oscillation-amo
amo.anomaly = amo.anomaly[merra2_amo_grid[,2]<=65,]
lat = unique(merra2_amo_grid[merra2_amo_grid[,2]<=65,2])
lon = unique(merra2_amo_grid[,1])
point_area = area_grid(lat, lon)
tot_area = sum(as.numeric(point_area)[!is.na(amo.anomaly[,1])])

weighted_anomaly = rep(0,480)
for(i in 1:480)
{
  weighted_anomaly[i] = sum(as.numeric(point_area)*amo.anomaly[,i], na.rm = TRUE)/tot_area
}

#plot(weighted_anomaly, type = 'l')


#Plot weighted mean 10-year average
#plot(runMean(weighted_anomaly, 120), type = 'l')

#Lowpass Filter
lowpass.loess <- loess(y ~ x, data = data.frame(x = 1:480, y = weighted_anomaly), span = 0.3) ## control span to define the amount of smoothing




#Nino Index Plot
ten_lp_filter = runMean(predict(lowpass.loess, 1:480), 120)[120:480]
dates = seq(as.Date("1981/1/16"), by = "month", length.out = 480)
pdoindex.gg = ggplot() +
  geom_line(mapping = aes(x = as.yearmon(dates), 
                          y = weighted_anomaly), 
            lwd = 1, 
            lty = 1, 
            color = 'black') +
  geom_hline(yintercept = c(0),
             color = 'gray50', 
             lwd = 1,
             lty = 2) +
  geom_vline(xintercept = as.yearmon(dates[c(141,356)]) ,
             color = 'indianred1',
             lwd = 1,
             lty =1, 
             alpha = 0.4) +
  geom_point(mapping = aes(x = as.yearmon(dates[c(141,356)])),
             y = weighted_anomaly[c(141,356)],
             color = 'indianred1', 
             shape = 19, 
             size = 2) +
  geom_line(mapping = aes(x = as.yearmon(dates[120:480]),
                          y = ten_lp_filter),
            color = 'navy',
            lwd = 1, 
            lty = 1) +
  theme_bw() +
  labs(title = 'Atlantic Multi-Decadal Oscialltion Area Averaged SST Anomaly', x = '', y = 'Temp. Anomaly (°C)') + 
  theme(plot.title = element_text(hjust = 0.5, size = 16), 
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 12)) + 
  scale_y_continuous(breaks = seq(-0.5, 0.5, 0.25),
                     limits = c(-0.6,0.6)) +
  scale_x_continuous(breaks = seq(1980, 2021, 5),
                     limits = c(as.yearmon(dates[1]), as.yearmon(dates[480])))




panel.gg = ggarrange(posneg.gg, pdoindex.gg,
                     nrow = 2, ncol = 1,
                     heights = c(1.75,1), 
                     labels = c('', 'C'))



#################
### Figure S3 ###
#################


#Comparison of M2 and C6 Nino3.4 Indices
load('merra2_enso_grid.RData')
load('MERRA2_Enso_Anomaly.RData')
full_nino_anomaly = cbind(merra2_enso_grid, enso.anomaly)

inside = locations.inside(merra2_enso_grid, elnino_border, as.is = FALSE)
colnames(inside) = c('lon', 'lat')
colnames(full_nino_anomaly) = c('lon', 'lat')

nino34 = data.frame(inside) %>% left_join(data.frame(full_nino_anomaly), by = c('lon', 'lat'))
nino34_index = apply(nino34[,-c(1,2)], 2, mean)




load('cmip6_enso_grid.RData')
load('1950_CMIP6_ENSO_Anomaly.RData')
load('CMIP6_ENSO_Anomaly.RData')
borderfix = cbind(rep(360, 8140), rep(0, 8140))
full_cmip6_anomaly = cbind(cmip6_enso_grid-borderfix,
                           cmip6_anomaly_1950,
                           cmip6_anomaly)


inside = locations.inside(cmip6_enso_grid-borderfix,
                          elnino_border,
                          as.is = FALSE)
colnames(inside) = c('lon', 'lat')


nino34.c = data.frame(inside) %>% left_join(data.frame(full_cmip6_anomaly), by = c('lon', 'lat'))
nino34_index_c = apply(nino34.c[,-c(1,2)], 2, mean)


dates = seq(as.Date("1981/1/16"), by = "month", length.out = 480)
dates_cmip6 = seq(as.Date("1950/1/16"), by = "month", length.out = 1812) #jan81 = 373, Dec20 = 852 



nino.gg = ggplot() +
  geom_line(mapping = aes(x = as.yearmon(dates), 
                          y = nino34_index), 
            lwd = 1, 
            lty = 1, 
            color = 'black') +
  geom_hline(yintercept = c(0.5,-0.5),
             color = 'gray50', 
             lwd = 1,
             lty = 2) +
  theme_bw() +
  labs(title = 'MERRA2 Niño3.4 Index', x = '', y = 'Temp. Anomaly (°C)') + 
  theme(plot.title = element_text(hjust = 0.5, size = 16), 
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 12)) + 
  scale_y_continuous(breaks = seq(-4, 4, 2),
                     limits = c(-4.25,4.25)) +
  scale_x_continuous(breaks = seq(1980, 2021, 5),
                     limits = c(as.yearmon(dates[1]), as.yearmon(dates[480])))

filt = 373:852
nino.gg.c = ggplot() +
  geom_line(mapping = aes(x = as.yearmon(dates_cmip6[filt]), 
                          y = nino34_index_c), 
            lwd = 1, 
            lty = 1, 
            color = 'black') +
  geom_hline(yintercept = c(0.5,-0.5),
             color = 'gray50', 
             lwd = 1,
             lty = 2) +
  theme_bw() +
  labs(title = 'CMIP6-MIROC Niño3.4 Index', x = '', y = 'Temp. Anomaly (°C)') + 
  theme(plot.title = element_text(hjust = 0.5, size = 16), 
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 12)) + 
  scale_y_continuous(breaks = seq(-4, 4, 2),
                     limits = c(-4.25,4.25)) +
  scale_x_continuous(breaks = seq(1980, 2021, 5),
                     limits = c(as.yearmon(dates[1]), as.yearmon(dates[480])))



#Functional boxplot of N34 Index for all 30 CMIPs
inside = locations.inside(merra2_enso_grid[complete.cases(enso.anomaly),],
                          elnino_border,
                          as.is = FALSE)
colnames(inside) = c('Var1', 'Var2')

n34_reals = matrix(NaN, nrow = 30, ncol = 480)
filt = (373:852)+1200
dates_cmip6 = seq(as.Date("1850/1/16"), by = "month", length.out = 3012)
for(real in 30:1)
{
  
  load(file = paste0('Full_TL_CMIP6_ENSO_Anomaly_R', real, '.RData'))
  
  full_cmip6_anomaly = cbind(merra2_enso_grid[complete.cases(enso.anomaly),],
                             full_tl_cmipenso_anomaly[,filt])
  
  nino34.c = data.frame(inside) %>% left_join(data.frame(full_cmip6_anomaly), by = c('Var1', 'Var2'))
  nino34_index_c = apply(nino34.c[,-c(1,2)], 2, mean)
  
  
  n34_reals[real,] = nino34_index_c
  
  print(real)
}



#Functional Boxplots
fit = t(n34_reals)
x = NULL
method = "MBD"
depth = NULL
plot = FALSE
prob = 0.5
color = 6
outliercol = 2
barcol = 4
fullout = FALSE
factor = 1.5
xlim = c(1, nrow(fit))
ylim = c(-5,5)
depth = fbplot(t(n34_reals), plot = FALSE)$depth



tp = dim(fit)[1]
n = dim(fit)[2]
if (length(x) == 0) {
  x = 1:tp
}
dp_s = sort(depth, decreasing = TRUE)
index = order(depth, decreasing = TRUE)
med = depth == max(depth)
medavg = matrix(fit[, med], ncol = sum(med), nrow = tp)
y = apply(medavg, 1, mean)
if (plot) {
  plot(x, y, lty = 1, lwd = 2, col = 1, type = "l", 
       xlim, ylim)
}
for (pp in 1:length(prob)) {
  m = ceiling(n * prob[pp])
  center = fit[, index[1:m]]
  out = fit[, index[(m + 1):n]]
  inf = apply(center, 1, min)
  sup = apply(center, 1, max)
  if (prob[pp] == 0.5) {
    dist = factor * (sup - inf)
    upper = sup + dist
    lower = inf - dist
    outly = (fit <= lower) + (fit >= upper)
    outcol = colSums(outly)
    remove = (outcol > 0)
    colum = 1:n
    outpoint = colum[remove == 1]
    out = fit[, remove]
    woout = fit
    good = woout[, (remove == 0), drop = FALSE]
    maxcurve = apply(good, 1, max)
    mincurve = apply(good, 1, min)
    if (sum(outly) > 0) {
      if (plot) {
        matlines(x, out, lty = 2, col = outliercol, 
                 type = "l", ...)
      }
    }
    barval = (x[1] + x[tp])/2
    bar = which(sort(c(x, barval)) == barval)[1]
    if (plot) {
      lines(c(x[bar], x[bar]), c(maxcurve[bar], sup[bar]), 
            col = barcol, lwd = 2)
      lines(c(x[bar], x[bar]), c(mincurve[bar], inf[bar]), 
            col = barcol, lwd = 2)
    }
  }
  xx = c(x, x[order(x, decreasing = TRUE)])
  supinv = sup[order(x, decreasing = TRUE)]
  yy = c(inf, supinv)
  if (plot) {
    if (prob[pp] == 0.5) {
      polygon(xx, yy, col = color[pp], border = barcol, 
              lwd = 2)
    }
    else {
      polygon(xx, yy, col = color[pp], border = NA)
    }
  }
}
if (plot) {
  lines(x, fit[, index[1]], lty = 1, lwd = 2, col = 1, 
        type = "l")
  lines(x, maxcurve, col = barcol, lwd = 2)
  lines(x, mincurve, col = barcol, lwd = 2)
  if (fullout) {
    if (sum(outly) > 0) {
      if (plot) {
        matlines(x, out, lty = 2, col = outliercol, 
                 type = "l", ...)
      }
    }
  }
}

medcurve = y

c6_fb = ggplot()+
  geom_polygon(mapping = aes(x = as.yearmon(dates_cmip6[filt][xx]),
                             y = yy),
               fill = 'gray60')+
  geom_line(mapping = aes(x = as.yearmon(dates_cmip6[filt]),
                          y = medcurve),
            color = 'black',
            lwd = 1,
            lty = 1) +
  geom_line(mapping = aes(x = as.yearmon(dates_cmip6[filt]),
                          y = maxcurve),
            color = 'blue',
            lwd = 1,
            lty = 1) +
  geom_line(mapping = aes(x = as.yearmon(dates_cmip6[filt]),
                          y = mincurve),
            color = 'blue',
            lwd = 1,
            lty = 1) +
  theme_bw() +
  labs(title = 'Functional Boxplot of CMIP6-MIROC Niño3.4 Indices', x = '', y = 'Temp. Anomaly (°C)') + 
  theme(plot.title = element_text(hjust = 0.5, size = 16), 
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 12)) +
  scale_x_continuous(breaks = seq(1980, 2021, 5),
                     limits = c(as.yearmon(dates_cmip6[filt][1]), as.yearmon(dates_cmip6[filt][480])))



#c6_fb


panel.gg = ggarrange(nino.gg, nino.gg.c, c6_fb,
                     nrow = 3, ncol = 1,
                     labels = c('A', 'B', 'C'))





#################
### Figure S4 ###
#################



#QQ-Plot of Residuals
load(paste0('C6_ENSO_NH_R134.RData'))
load('QDist_ENSO_CMIP6.RData')

file = paste0('Full_CMIP6_ENSO_EOF.RData')
load(file)
rawData = full_cmipenso_eof
index = 31
tau = 1
trainLen = 1739 + (index-1)*36 #total = 2819
forward = 36

#Generate training/testing/valid sets
sets = cttv(rawData = rawData,
            tau = tau,
            trainLen = trainLen,
            testLen = forward,
            valid.flag = F)


# file = paste0('C6_ENSO_Windows', index, '.RData')
# load(file)
# sd.pred = sapply(1:20, function(x) apply(ensemb.pred[x,,], 2, sd))


new.resids = (sets$yTest - mean.pred)

new.resids = data.frame('nr' = as.vector(new.resids))


new.density = ggplot(new.resids, aes(nr)) +
  geom_density(color = 'black', lwd = 1.05) +
  labs(x = bquote('Residuals'), y = 'Density', title = 'Density of CMIP6-MIROC ENSO Residuals') +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5, size = 18),
        axis.text = element_text(size = 14),
        axis.title = element_text(size = 16))

new.qq = ggplot(new.resids, aes(sample = nr)) +
  stat_qq(size = 1.25) +
  stat_qq_line(lwd = 1.05, color = 'red') +
  labs(x = 'Theoretical', y = bquote('Sample'), title = 'QQ-Plot of CMIP6-MIROC ENSO Residuals') +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5, size = 18),
        axis.text = element_text(size = 14),
        axis.title = element_text(size = 16)) +
  scale_x_continuous(breaks = seq(-3,3,1)) +
  scale_y_continuous(position = 'right')


qq = ggarrange(new.density, new.qq, 
               nrow = 1, ncol = 2, 
               labels = c('A', 'B'))

#################
### Figure S5 ###
#################


load('C6_Optim_Params.RData')

nu = optim_params[,5]
pi.w = optim_params[,6]
pi.win = optim_params[,7]
alpha = optim_params[,8]



nu.gg = ggplot() + 
  geom_density(mapping = aes(x = nu),
               lwd = 2) +
  expand_limits(x = c(0,1)) + 
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5, size = 16),
        axis.text = element_text(size = 12), 
        axis.title = element_text(size = 14)) +
  labs(y = 'Density', x = expression(nu), title = '')


piw.gg = ggplot() + 
  geom_density(mapping = aes(x = pi.w),
               lwd = 2) +
  expand_limits(x = c(0,1)) + 
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5, size = 16),
        axis.text = element_text(size = 12), 
        axis.title = element_text(size = 14)) +
  labs(y = 'Density', x = expression(pi[W]), title = '')


piwin.gg = ggplot() + 
  geom_density(mapping = aes(x = pi.win),
               lwd = 2) +
  expand_limits(x = c(0,1)) + 
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5, size = 16),
        axis.text = element_text(size = 12), 
        axis.title = element_text(size = 14)) +
  labs(y = 'Density', x = expression(pi[W^"in"]), title = '') +
  scale_y_continuous(position = 'right')

alpha.gg = ggplot() + 
  geom_density(mapping = aes(x = alpha),
               lwd = 2) +
  expand_limits(x = c(0,1)) + 
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5, size = 16),
        axis.text = element_text(size = 12), 
        axis.title = element_text(size = 14)) +
  labs(y = 'Density', x = expression(alpha), title = '') +
  scale_y_continuous(breaks = seq(0,3,0.5), position = 'right')



ggarrange(nu.gg, alpha.gg, piw.gg, piwin.gg, 
          nrow = 2, ncol = 2, 
          labels = c('A', 'B', 'C', 'D'))





#################
### Figure S6 ###
#################

index = 31

load('Full_CMIP6_PDO_EOF.RData')
load('QDist_PDO_CMIP6.RData')
file = paste0('Full_C6_PDO_Windows', index, '.RData')
load(file)
file = paste0('C6_PDO_Alpha_CV', 0.13, '.RData')
load(file)

#Coverage performance
tau = 1
trainLen = 1895 + (index-1)*36 #total = 839
testLen = 1
forward = 36

#Generate training/testing/valid sets
rawData = full_cmippdo_eof
sets = cttv(rawData, tau, trainLen, testLen = forward)


eof_num = 1
#mean.pred = sapply(1:1, function(x) apply(ensemb.pred[x,,], 2, mean))
upper = mean.pred + (quant_dist[,2])
lower = mean.pred - (quant_dist[,1])




dates_cmip6 = seq(as.Date("1850/1/16"), by = "month", length.out = 3012)
forc.dates = tail(as.POSIXct(dates_cmip6), 36)

eof_gg = ggplot() +
  geom_line(mapping = aes(x = forc.dates,
                          y = sets$yTest,
                          color = 'black'),
            lty = 1,
            lwd = 1.5) +
  geom_line(mapping = aes(x = forc.dates, 
                          y = upper,
                          color = 'red'),
            lty = 2,
            lwd = 1.5) +
  geom_line(mapping = aes(x = forc.dates,
                          y = mean.pred,
                          color = 'red1'), 
            lty = 1,
            lwd = 1.5) +
  geom_line(mapping = aes(x = forc.dates,
                          y = lower),
            color = 'red1',
            lty = 2,
            lwd = 1.5) + 
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5, size = 16),
        axis.text = element_text(size = 12), 
        axis.title = element_text(size = 14)) +
  labs(y = '', x = 'Time', title = 'CMIP6-MIROC PDO Index', color = '') +
  scale_color_manual(values = c('black',  'red', 'red1' ), 
                     labels = c('True Value', 'Mean Forecasts', '95% Prediction Interval'),
                     guide = guide_legend(override.aes = list(lty = c(1,1,2)))) +
  theme(legend.position = c(0.45,0.9), legend.text = element_text(size = 14)) +
  theme(legend.key.width = unit(2.5,"cm")) +
  scale_y_continuous(breaks = seq(-1,1,0.5), limits = c(-1, 1)) +
  scale_x_datetime(date_labels = '%Y')


eof_gg



















