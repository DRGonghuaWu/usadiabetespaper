setwd("E:/1工作/项目/美国糖尿病数据/数据")

library(tidyverse)
library(gnm)
library(CBCgrps)
library(qgcomp)
library(qgcompint)
library(forestploter)
library(splines)
library(gWQS)
colors1 = c("#fee090","#abd9e9","#74add1","#4575b4","#fdae61","#f46d43","#d73027")

source("used functions.R")
Anal_dat = openxlsx::read.xlsx("AnalysisData0916.xlsx")
Incidence_dat = openxlsx::read.xlsx("Incidencedata08_17.xlsx")
Incidence_dat$spatialid = paste0(Incidence_dat$State,Incidence_dat$County)
Incidence_dat$YEAR = as.character(Incidence_dat$YEAR)
Final_dat = Anal_dat %>% left_join(.,Incidence_dat[,c("spatialid","YEAR","Rate")], by = c("spatialid", "YEAR")) %>%
  mutate(Newcase_20 = ifelse(is.na(Newcase_20),round(0.001*Rate*Population_20,0),Newcase_20))
exposures <- c("mean_EC","mean_OC","mean_NO3","mean_SO4","mean_NH4")

#### ========== Descriptive analysis =========== ####
mean(Final_dat$incidence, na.rm = T);sd(Final_dat$incidence, na.rm = T)
quantile(Final_dat$incidence, na.rm = T);quantile(Final_dat$incidence, na.rm = T)[4]-quantile(Final_dat$incidence, na.rm = T)[2]
mean(Final_dat$mortality, na.rm = T);sd(Final_dat$mortality, na.rm = T)
quantile(Final_dat$mortality, na.rm = T);quantile(Final_dat$mortality, na.rm = T)[4]-quantile(Final_dat$mortality, na.rm = T)[2]

mean(Final_dat$PM_25, na.rm = T);sd(Final_dat$PM_25, na.rm = T)
quantile(Final_dat$PM_25, na.rm = T);quantile(Final_dat$PM_25, na.rm = T)[4]-quantile(Final_dat$PM_25, na.rm = T)[2]
mean(Final_dat$mean_EC, na.rm = T);sd(Final_dat$mean_EC, na.rm = T)
quantile(Final_dat$mean_EC, na.rm = T);quantile(Final_dat$mean_EC, na.rm = T)[4]-quantile(Final_dat$mean_EC, na.rm = T)[2]
mean(Final_dat$mean_OC, na.rm = T);sd(Final_dat$mean_OC, na.rm = T)
quantile(Final_dat$mean_OC, na.rm = T);quantile(Final_dat$mean_OC, na.rm = T)[4]-quantile(Final_dat$mean_OC, na.rm = T)[2]
mean(Final_dat$mean_SO4, na.rm = T);sd(Final_dat$mean_SO4, na.rm = T)
quantile(Final_dat$mean_SO4, na.rm = T);quantile(Final_dat$mean_SO4, na.rm = T)[4]-quantile(Final_dat$mean_SO4, na.rm = T)[2]
mean(Final_dat$mean_NO3, na.rm = T);sd(Final_dat$mean_NO3, na.rm = T)
quantile(Final_dat$mean_NO3, na.rm = T);quantile(Final_dat$mean_NO3, na.rm = T)[4]-quantile(Final_dat$mean_NO3, na.rm = T)[2]
mean(Final_dat$mean_NH4, na.rm = T);sd(Final_dat$mean_NH4, na.rm = T)
quantile(Final_dat$mean_NH4, na.rm = T);quantile(Final_dat$mean_NH4, na.rm = T)[4]-quantile(Final_dat$mean_NH4, na.rm = T)[2]

mean(Final_dat$Summer_temp, na.rm = T);sd(Final_dat$Summer_temp, na.rm = T)
quantile(Final_dat$Summer_temp, na.rm = T);quantile(Final_dat$Summer_temp, na.rm = T)[4]-quantile(Final_dat$Summer_temp, na.rm = T)[2]
mean(Final_dat$Winter_temp, na.rm = T);sd(Final_dat$Winter_temp, na.rm = T)
quantile(Final_dat$Winter_temp, na.rm = T);quantile(Final_dat$Winter_temp, na.rm = T)[4]-quantile(Final_dat$Winter_temp, na.rm = T)[2]
mean(Final_dat$Summer_sd, na.rm = T);sd(Final_dat$Summer_sd, na.rm = T)
quantile(Final_dat$Summer_sd, na.rm = T);quantile(Final_dat$Summer_sd, na.rm = T)[4]-quantile(Final_dat$Summer_sd, na.rm = T)[2]
mean(Final_dat$Winter_sd, na.rm = T);sd(Final_dat$Winter_sd, na.rm = T)
quantile(Final_dat$Winter_sd, na.rm = T);quantile(Final_dat$Winter_sd, na.rm = T)[4]-quantile(Final_dat$Winter_sd, na.rm = T)[2]

mean(Final_dat$Population_20, na.rm = T);sd(Final_dat$Population_20, na.rm = T)
quantile(Final_dat$Population_20, na.rm = T);quantile(Final_dat$Population_20, na.rm = T)[4]-quantile(Final_dat$Population_20, na.rm = T)[2]
mean(Final_dat$GDP_PER, na.rm = T);sd(Final_dat$GDP_PER, na.rm = T)
quantile(Final_dat$GDP_PER, na.rm = T);quantile(Final_dat$GDP_PER, na.rm = T)[4]-quantile(Final_dat$GDP_PER, na.rm = T)[2]
mean(Final_dat$GDP_PER, na.rm = T);sd(Final_dat$GDP_PER, na.rm = T)
quantile(Final_dat$GDP_PER, na.rm = T);quantile(Final_dat$GDP_PER, na.rm = T)[4]-quantile(Final_dat$GDP_PER, na.rm = T)[2]
mean(Final_dat$AGEING, na.rm = T);sd(Final_dat$AGEING, na.rm = T) # because only few counties were not aging, so it was included in the first group
quantile(Final_dat$AGEING, na.rm = T);quantile(Final_dat$AGEING, na.rm = T)[4]-quantile(Final_dat$AGEING, na.rm = T)[2]
mean(Final_dat$obesity_20, na.rm = T);sd(Final_dat$obesity_20, na.rm = T) 
quantile(Final_dat$obesity_20, na.rm = T);quantile(Final_dat$obesity_20, na.rm = T)[4]-quantile(Final_dat$obesity_20, na.rm = T)[2]
mean(Final_dat$Inactivity_20, na.rm = T);sd(Final_dat$Inactivity_20, na.rm = T) 
quantile(Final_dat$Inactivity_20, na.rm = T);quantile(Final_dat$Inactivity_20, na.rm = T)[4]-quantile(Final_dat$Inactivity_20, na.rm = T)[2]
mean(Final_dat$poverty_18, na.rm = T);sd(Final_dat$poverty_18, na.rm = T) 
quantile(Final_dat$poverty_18, na.rm = T);quantile(Final_dat$poverty_18, na.rm = T)[4]-quantile(Final_dat$poverty_18, na.rm = T)[2]
mean(Final_dat$INSURE, na.rm = T);sd(Final_dat$INSURE, na.rm = T) 
quantile(Final_dat$INSURE, na.rm = T);quantile(Final_dat$INSURE, na.rm = T)[4]-quantile(Final_dat$INSURE, na.rm = T)[2]

## plot
Fig1_dat <- Final_dat %>% group_by(YEAR) %>%
  summarise(Mort_mean = mean(mortality,na.rm = T),
            Mort_25 = quantile(mortality,na.rm = T)[2],
            Mort_75 = quantile(mortality,na.rm = T)[4],
            Inc_mean = mean(incidence,na.rm = T),
            Inc_25 = quantile(incidence,na.rm = T)[2],
            Inc_75 = quantile(incidence,na.rm = T)[4],
            EC_mean = mean(mean_EC,na.rm = T),
            EC_25 = quantile(mean_EC,na.rm = T)[2],
            EC_75 = quantile(mean_EC,na.rm = T)[4],
            OC_mean = mean(mean_OC,na.rm = T),
            OC_25 = quantile(mean_OC,na.rm = T)[2],
            OC_75 = quantile(mean_OC,na.rm = T)[4],
            SO4_mean = mean(mean_SO4,na.rm = T),
            SO4_25 = quantile(mean_SO4,na.rm = T)[2],
            SO4_75 = quantile(mean_SO4,na.rm = T)[4],
            NO3_mean = mean(mean_NO3,na.rm = T),
            NO3_25 = quantile(mean_NO3,na.rm = T)[2],
            NO3_75 = quantile(mean_NO3,na.rm = T)[4],
            NH4_mean = mean(mean_NH4,na.rm = T),
            NH4_25 = quantile(mean_NH4,na.rm = T)[2],
            NH4_75 = quantile(mean_NH4,na.rm = T)[4]) %>%
  mutate(YEAR = as.numeric(YEAR))
p1.1 = ggplot(data = Fig1_dat, aes(x = YEAR)) +
  geom_line(aes(y = Inc_mean))+
  geom_point(aes(y = Inc_mean),shape = 10) +
  geom_errorbar(aes(ymin = Inc_25, ymax = Inc_75), width = 0.2) +
  geom_hline(yintercept = 9.14, color = "#E64B35CC", lty = 2) +
  scale_x_continuous(name = "", breaks = c(2008,2010,2012,2014,2016,2018)) +
  scale_y_continuous(name = "Diabetes incidence (/1,000)") +
  theme_bw()+
  theme(axis.title = element_text(size = 9))
p1.2 = ggplot(data = Fig1_dat, aes(x = YEAR)) +
  geom_line(aes(y = Mort_mean))+
  geom_point(aes(y = Mort_mean),shape = 10) +
  geom_errorbar(aes(ymin = Mort_25, ymax = Mort_75), width = 0.2) +
  geom_hline(yintercept = 142.38, color = "#E64B35CC", lty = 2) +
  scale_x_continuous(name = "", breaks = c(2008,2010,2012,2014,2016,2018)) +
  scale_y_continuous(name = "Diabetes mortality (/100,000)") +
  theme_bw()+
  theme(axis.title = element_text(size = 9))
p1.3 = ggplot(data = Fig1_dat, aes(x = YEAR)) +
  geom_line(aes(y = EC_mean), color = "#ef1828")+
  geom_line(aes(y = OC_mean), color = "#f88421")+
  geom_line(aes(y = SO4_mean), color = "#ffbc14")+
  geom_line(aes(y = NO3_mean), color = "#00bdcd")+
  geom_line(aes(y = NH4_mean), color = "#006b7b")+
  geom_point(aes(y = EC_mean),shape = 10, color = "#ef1828") +
  geom_point(aes(y = OC_mean),shape = 10, color = "#f88421") +
  geom_point(aes(y = SO4_mean),shape = 10, color = "#ffbc14") +
  geom_point(aes(y = NO3_mean),shape = 10, color = "#00bdcd") +
  geom_point(aes(y = NH4_mean),shape = 10, color = "#006b7b") +
  geom_rect(aes(xmin = 2014.5,xmax=2015.1,ymin=2.8,ymax=2.8),color="#ef1828")+
  annotate("text", x = 2015.4, y = 2.8, label = "Elemental carbon",size = 3,hjust=0)+
  geom_rect(aes(xmin = 2014.5,xmax=2015.1,ymin=2.6,ymax=2.6),color="#f88421")+
  annotate("text", x = 2015.4, y = 2.6, label = "Organic carbon",size = 3,hjust=0)+
  geom_rect(aes(xmin = 2014.5,xmax=2015.1,ymin=2.4,ymax=2.4),color="#ffbc14")+
  annotate("text", x = 2015.4, y = 2.4, label = "Sulfate",size = 3,hjust=0)+
  geom_rect(aes(xmin = 2014.5,xmax=2015.1,ymin=2.2,ymax=2.2),color="#00bdcd")+
  annotate("text", x = 2015.4, y = 2.2, label = "Nitrate",size = 3,hjust=0)+
  geom_rect(aes(xmin = 2014.5,xmax=2015.1,ymin=2.0,ymax=2.0),color="#006b7b")+
  annotate("text", x = 2015.4, y = 2.0, label = "Ammonium",size = 3,hjust=0)+
  scale_x_continuous(name = "", breaks = c(2008,2010,2012,2014,2016,2018)) +
  scale_y_continuous(name = expression(paste("Concentrations of component (μg/ ",m^3,")")), breaks = c(0,0.5,1,1.5,2,2.5,3)) +
  theme_bw()+
  theme(axis.title = element_text(size = 9))
tiff("Figure108_17.tiff", res = 900, width = 6,height = 8, units = "in")
ggpubr::ggarrange(p1.1,p1.2,p1.3, ncol = 1, nrow = 3, align = "v", labels = c("A","B","C"), font.label = list(size = 12),
                  label.x = 0.05, label.y = 1.05) +
  theme(plot.margin = margin(t = 0.5, r = 0, b = 0, l = 0, unit = "cm"))
dev.off()

cor_matrix <- Final_dat[,exposures] %>% setNames(c("EC","OC","SO4","NO3","NH4")) %>% cor()
tiff("FigureS108_17.tiff", res = 900, width = 6,height = 6, units = "in")
corrplot::corrplot(cor_matrix, cl.pos = "b",tl.cex = 0.8,tl.col = "black",tl.srt = 45)
dev.off()

## maps
library(sf)
library(rnaturalearth)
library(rnaturalearthdata)
library(ggspatial)
USA <- read_sf("gadm36_USA_2.shp")
MapDat <- read.csv("MapDat.csv")
USA_dat = USA %>% left_join(MapDat, by = "spatialid")
ps2 = ggplot(data = USA_dat) + 
  geom_sf(aes(fill = M_mort)) +
  annotation_scale(location = "bl") +
  annotation_north_arrow(location = "bl", which_north = "true", height = unit(1,"cm"),width = unit(1,"cm"),
                         pad_x = unit(0.2, "in"), pad_y = unit(0.3, "in"),
                         style = north_arrow_fancy_orienteering) +
  coord_sf(xlim = c(-126, -64), ylim = c(23, 50), expand = FALSE) +
  scale_fill_gradient2(name = "Diabetes mortality\n(/100,000)", low = "#fff5f0",mid = "#fb6a4a",high = "#99000d", na.value = "white",midpoint = 250)+
  theme_bw()+
  theme(legend.position = c(0.8,0.42),
        legend.justification = c(0,1),
        legend.key.size = unit(10,"pt"),
        legend.key.width = unit(20,"pt"),
        legend.key.height = unit(10,"pt"),
        legend.title = element_text(size = 10),
        panel.grid = element_blank())
tiff("FigureS2_mortality_08_17.tiff",res = 900, width = 7,height = 4.5,units = "in")
ps2
dev.off()

ps3 = ggplot(data = USA_dat) + 
  geom_sf(aes(fill = M_PM)) +
  annotation_scale(location = "bl") +
  annotation_north_arrow(location = "bl", which_north = "true", height = unit(1,"cm"),width = unit(1,"cm"),
                         pad_x = unit(0.2, "in"), pad_y = unit(0.3, "in"),
                         style = north_arrow_fancy_orienteering) +
  coord_sf(xlim = c(-126, -64), ylim = c(23, 50), expand = FALSE) +
  scale_fill_gradient2(name = expression(paste(PM[2.5]," (μg/",m^3,")")), low = "#fff5f0",mid = "#fb6a4a",high = "#99000d", na.value = "white",midpoint = 8)+
  theme_bw()+
  theme(legend.position = c(0.8,0.42),
        legend.justification = c(0,1),
        legend.key.size = unit(10,"pt"),
        legend.key.width = unit(20,"pt"),
        legend.key.height = unit(10,"pt"),
        legend.title = element_text(size = 10),
        panel.grid = element_blank())
tiff("FigureS3_PM2.5_08_17.tiff",res = 900, width = 7,height = 4.5,units = "in")
ps3
dev.off()

ps4 = ggplot(data = USA_dat) + 
  geom_sf(aes(fill = M_EC)) +
  annotation_scale(location = "bl") +
  annotation_north_arrow(location = "bl", which_north = "true", height = unit(1,"cm"),width = unit(1,"cm"),
                         pad_x = unit(0.2, "in"), pad_y = unit(0.3, "in"),
                         style = north_arrow_fancy_orienteering) +
  coord_sf(xlim = c(-126, -64), ylim = c(23, 50), expand = FALSE) +
  scale_fill_gradient2(name = expression(paste("EC proportion(%)")), low = "#ffffcc",mid = "#fd8d3c",high = "#800026", na.value = "white",midpoint = 10)+
  theme_bw()+
  theme(legend.position = c(0.78,0.42),
        legend.justification = c(0,1),
        legend.key.size = unit(10,"pt"),
        legend.key.width = unit(20,"pt"),
        legend.key.height = unit(10,"pt"),
        legend.title = element_text(size = 10),
        panel.grid = element_blank())
tiff("FigureS4_EC_08_17.tiff",res = 900, width = 7,height = 4.5,units = "in")
ps4
dev.off()

ps5 = ggplot(data = USA_dat) + 
  geom_sf(aes(fill = M_OC)) +
  annotation_scale(location = "bl") +
  annotation_north_arrow(location = "bl", which_north = "true", height = unit(1,"cm"),width = unit(1,"cm"),
                         pad_x = unit(0.2, "in"), pad_y = unit(0.3, "in"),
                         style = north_arrow_fancy_orienteering) +
  coord_sf(xlim = c(-126, -64), ylim = c(23, 50), expand = FALSE) +
  scale_fill_gradient2(name = expression(paste("OC proportion(%)")), low = "#ffffcc",mid = "#fd8d3c",high = "#800026", na.value = "white",midpoint = 30)+
  theme_bw()+
  theme(legend.position = c(0.78,0.42),
        legend.justification = c(0,1),
        legend.key.size = unit(10,"pt"),
        legend.key.width = unit(20,"pt"),
        legend.key.height = unit(10,"pt"),
        legend.title = element_text(size = 10),
        panel.grid = element_blank())
tiff("FigureS5_OC_08_17.tiff",res = 900, width = 7,height = 4.5,units = "in")
ps5
dev.off()

ps6 = ggplot(data = USA_dat) + 
  geom_sf(aes(fill = M_NO3)) +
  annotation_scale(location = "bl") +
  annotation_north_arrow(location = "bl", which_north = "true", height = unit(1,"cm"),width = unit(1,"cm"),
                         pad_x = unit(0.2, "in"), pad_y = unit(0.3, "in"),
                         style = north_arrow_fancy_orienteering) +
  coord_sf(xlim = c(-126, -64), ylim = c(23, 50), expand = FALSE) +
  scale_fill_gradient2(name = expression(paste(NO[3]^{"-"}~proportion,"\n(%)")), low = "#ffffcc",mid = "#fd8d3c",high = "#800026", na.value = "white",midpoint = 20)+
  theme_bw()+
  theme(legend.position = c(0.78,0.42),
        legend.justification = c(0,1),
        legend.key.size = unit(10,"pt"),
        legend.key.width = unit(20,"pt"),
        legend.key.height = unit(10,"pt"),
        legend.title = element_text(size = 10),
        panel.grid = element_blank())
tiff("FigureS6_NO3_08_17.tiff",res = 900, width = 7,height = 4.5,units = "in")
ps6
dev.off()

ps7 = ggplot(data = USA_dat) + 
  geom_sf(aes(fill = M_SO4)) +
  annotation_scale(location = "bl") +
  annotation_north_arrow(location = "bl", which_north = "true", height = unit(1,"cm"),width = unit(1,"cm"),
                         pad_x = unit(0.2, "in"), pad_y = unit(0.3, "in"),
                         style = north_arrow_fancy_orienteering) +
  coord_sf(xlim = c(-126, -64), ylim = c(23, 50), expand = FALSE) +
  scale_fill_gradient2(name = expression(paste(SO[4]^{"2-"}~proportion,"\n(%)")), 
                       low = "#ffffcc",mid = "#fd8d3c",high = "#800026", na.value = "white",midpoint = 25)+
  theme_bw()+
  theme(legend.position = c(0.78,0.42),
        legend.justification = c(0,1),
        legend.key.size = unit(10,"pt"),
        legend.key.width = unit(20,"pt"),
        legend.key.height = unit(10,"pt"),
        legend.title = element_text(size = 10),
        panel.grid = element_blank())
tiff("FigureS7_SO4_08_17.tiff",res = 900, width = 7,height = 4.5,units = "in")
ps7
dev.off()

ps8 = ggplot(data = USA_dat) + 
  geom_sf(aes(fill = M_NH4)) +
  annotation_scale(location = "bl") +
  annotation_north_arrow(location = "bl", which_north = "true", height = unit(1,"cm"),width = unit(1,"cm"),
                         pad_x = unit(0.2, "in"), pad_y = unit(0.3, "in"),
                         style = north_arrow_fancy_orienteering) +
  coord_sf(xlim = c(-126, -64), ylim = c(23, 50), expand = FALSE) +
  scale_fill_gradient2(name = expression(paste(NH[4]^{"+"}~proportion,"\n(%)")), 
                       low = "#ffffcc",mid = "#fd8d3c",high = "#800026", na.value = "white",midpoint = 10)+
  theme_bw()+
  theme(legend.position = c(0.78,0.42),
        legend.justification = c(0,1),
        legend.key.size = unit(10,"pt"),
        legend.key.width = unit(20,"pt"),
        legend.key.height = unit(10,"pt"),
        legend.title = element_text(size = 10),
        panel.grid = element_blank())
tiff("FigureS8_NH4_08_17.tiff",res = 900, width = 7,height = 4.5,units = "in")
ps8
dev.off()

ps9 = ggplot(data = USA_dat) + 
  geom_sf(aes(fill = M_ST)) +
  annotation_scale(location = "bl") +
  annotation_north_arrow(location = "bl", which_north = "true", height = unit(1,"cm"),width = unit(1,"cm"),
                         pad_x = unit(0.2, "in"), pad_y = unit(0.3, "in"),
                         style = north_arrow_fancy_orienteering) +
  coord_sf(xlim = c(-126, -64), ylim = c(23, 50), expand = FALSE) +
  scale_fill_gradient2(name = "Summer mean \ntemperature(℃)",
                       low = "#006b7b",mid = "#ffffcc",high = "#800026", na.value = "white",midpoint = 20)+
  theme_bw()+
  theme(legend.position = c(0.78,0.42),
        legend.justification = c(0,1),
        legend.key.size = unit(10,"pt"),
        legend.key.width = unit(20,"pt"),
        legend.key.height = unit(10,"pt"),
        legend.title = element_text(size = 10),
        panel.grid = element_blank())
tiff("FigureS9_ST_08_17.tiff",res = 900, width = 7,height = 4.5,units = "in")
ps9
dev.off()

ps10 = ggplot(data = USA_dat) + 
  geom_sf(aes(fill = M_WT)) +
  annotation_scale(location = "bl") +
  annotation_north_arrow(location = "bl", which_north = "true", height = unit(1,"cm"),width = unit(1,"cm"),
                         pad_x = unit(0.2, "in"), pad_y = unit(0.3, "in"),
                         style = north_arrow_fancy_orienteering) +
  coord_sf(xlim = c(-126, -64), ylim = c(23, 50), expand = FALSE) +
  scale_fill_gradient2(name = "Winter mean \ntemperature(℃)",
                       low = "#00bdcd",mid = "#ffffcc",high = "#800026", na.value = "white",midpoint = 0)+
  theme_bw()+
  theme(legend.position = c(0.78,0.42),
        legend.justification = c(0,1),
        legend.key.size = unit(10,"pt"),
        legend.key.width = unit(20,"pt"),
        legend.key.height = unit(10,"pt"),
        legend.title = element_text(size = 10),
        panel.grid = element_blank())
tiff("FigureS10_WT_08_17.tiff",res = 900, width = 7,height = 4.5,units = "in")
ps10
dev.off()
#### ========== Modelling =========== ####
# Single component DID
Inc_PM <- gnm(Newcase_20 ~ mean_PM_iqr + spatialid + YEAR + Winter_temp +
                Summer_temp + Winter_sd + Summer_sd + GDP_PER, eliminate = as.factor(spatialid),
              family = quasipoisson, data=Final_dat, offset = log(Final_dat$Population_20))
Inc_EC <- gnm(Newcase_20 ~ mean_EC_iqr + spatialid + YEAR + Winter_temp +
                Summer_temp + Winter_sd + Summer_sd + GDP_PER, eliminate = as.factor(spatialid),
              family = quasipoisson, data=Final_dat, offset = log(Final_dat$Population_20))
Inc_OC <- gnm(Newcase_20 ~ mean_OC_iqr + spatialid + YEAR + Winter_temp +
                Summer_temp + Winter_sd + Summer_sd + GDP_PER, eliminate = as.factor(spatialid),
              family = quasipoisson, data=Final_dat, offset = log(Final_dat$Population_20))
Inc_NO3 <- gnm(Newcase_20 ~ mean_NO3_iqr + spatialid + YEAR + Winter_temp +
                 Summer_temp + Winter_sd + Summer_sd + GDP_PER, eliminate = as.factor(spatialid),
               family = quasipoisson, data=Final_dat, offset = log(Final_dat$Population_20))
Inc_SO4 <- gnm(Newcase_20 ~ mean_SO4_iqr + spatialid + YEAR + Winter_temp +
                 Summer_temp + Winter_sd + Summer_sd + GDP_PER, eliminate = as.factor(spatialid),
               family = quasipoisson, data=Final_dat, offset = log(Final_dat$Population_20))
Inc_NH4 <- gnm(Newcase_20 ~ mean_NH4_iqr + spatialid + YEAR + Winter_temp +
                 Summer_temp + Winter_sd + Summer_sd + GDP_PER, eliminate = as.factor(spatialid),
               family = quasipoisson, data=Final_dat, offset = log(Final_dat$Population_20))
round(Calc_IR(Inc_PM, component = "mean_PM_iqr"),3) # 1.054 0.585 1.525 0.000
round(Calc_IR(Inc_EC, component = "mean_EC_iqr"),3) # 0.209 -0.511  0.933  0.571
round(Calc_IR(Inc_OC, component = "mean_OC_iqr"),3) # -1.835 -2.678 -0.985  0.000
round(Calc_IR(Inc_SO4, component = "mean_SO4_iqr"),3) # 1.968 1.194 2.748 0.000
round(Calc_IR(Inc_NO3, component = "mean_NO3_iqr"),3) # -0.420 -1.112  0.278  0.238
round(Calc_IR(Inc_NH4, component = "mean_NH4_iqr"),3) # 0.741 0.026 1.462 0.042



Mort_PM <- gnm(Deaths_20 ~ mean_PM_iqr + spatialid + YEAR + Winter_temp +
                 Summer_temp + Winter_sd + Summer_sd + GDP_PER, eliminate = as.factor(spatialid),
               family = quasipoisson, data=Final_dat, offset = log(Final_dat$Population_20))
Mort_EC <- gnm(Deaths_20 ~ mean_EC_iqr + spatialid + YEAR + Winter_temp +
                 Summer_temp + Winter_sd + Summer_sd + GDP_PER, eliminate = as.factor(spatialid),
               family = quasipoisson, data=Final_dat, offset = log(Final_dat$Population_20))
Mort_OC <- gnm(Deaths_20 ~ mean_OC_iqr + spatialid + YEAR + Winter_temp +
                 Summer_temp + Winter_sd + Summer_sd + GDP_PER, eliminate = as.factor(spatialid),
               family = quasipoisson, data=Final_dat, offset = log(Final_dat$Population_20))
Mort_NO3 <- gnm(Deaths_20 ~ mean_NO3_iqr + spatialid + YEAR + Winter_temp +
                  Summer_temp + Winter_sd + Summer_sd + GDP_PER, eliminate = as.factor(spatialid),
                family = quasipoisson, data=Final_dat, offset = log(Final_dat$Population_20))
Mort_SO4 <- gnm(Deaths_20 ~ mean_SO4_iqr + spatialid + YEAR + Winter_temp +
                  Summer_temp + Winter_sd + Summer_sd + GDP_PER, eliminate = as.factor(spatialid),
                family = quasipoisson, data=Final_dat, offset = log(Final_dat$Population_20))
Mort_NH4 <- gnm(Deaths_20 ~ mean_NH4_iqr + spatialid + YEAR + Winter_temp +
                  Summer_temp + Winter_sd + Summer_sd + GDP_PER, eliminate = as.factor(spatialid),
                family = quasipoisson, data=Final_dat, offset = log(Final_dat$Population_20))
round(Calc_IR(Mort_PM, component = "mean_PM_iqr"),3) # 3.414 (2.750,4.083) P<0.001
round(Calc_IR(Mort_EC, component = "mean_EC_iqr"),3) # 11.205 (10.017,12.406) P<0.001
round(Calc_IR(Mort_OC, component = "mean_OC_iqr"),3) # 0.800 (-0.475,2.092) P=0.220
round(Calc_IR(Mort_SO4, component = "mean_SO4_iqr"),3)  # 9.686 (8.503,10.882) P<0.001
round(Calc_IR(Mort_NO3, component = "mean_NO3_iqr"),3)  # 3.862 (2.825,4.909) P<0.001
round(Calc_IR(Mort_NH4, component = "mean_NH4_iqr"),3)  # 10.185 (9.054,11.328) P<0.001


# Multiple components DID
Final_dat_inc1 <- Final_dat %>% filter(!is.na(Newcase_20) & !is.na(Winter_temp) & !is.na(GDP_PER))
Inc_qgcomp1 <- qgcomp.noboot(Newcase_20 ~ spatialid + YEAR + Winter_temp + Summer_temp + Winter_sd + 
                               Summer_sd + mean_EC + mean_OC + mean_SO4 + mean_NO3 + mean_NH4 + GDP_PER + offset(log(Population_20)),
                             data = Final_dat_inc1, expnms = exposures, family = poisson, q = 4)
Inc_qgcomp2 <- qgcomp.boot.quasi(Newcase_20 ~ spatialid + YEAR + Winter_temp + Summer_temp + Winter_sd + 
                                   Summer_sd + mean_EC + mean_OC + mean_SO4 + mean_NO3 + mean_NH4 + GDP_PER + offset(log(Population_20)),
                                 data = Final_dat_inc1, expnms = exposures, family = poisson, q = 4, 
                                 B = 100, seed = 123,parallel = TRUE,parplan = T)

(exp(summary(Inc_qgcomp2)$coefficients["psi1", "Estimate"])-1)*100
(exp(summary(Inc_qgcomp2)$coefficients["psi1", "Lower CI"])-1)*100
(exp(summary(Inc_qgcomp2)$coefficients["psi1", "Upper CI"])-1)*100 # -2.18%(-3.69%,-0.64%) P=0.006

Final_dat_mort1 <- Final_dat %>% filter(!is.na(Deaths_20) & !is.na(Winter_temp) & !is.na(GDP_PER))
Mort_qgcomp1 <- qgcomp.noboot(Deaths_20 ~ spatialid + YEAR + Winter_temp + Summer_temp + Winter_sd + 
                                Summer_sd + mean_EC + mean_OC + mean_SO4 + mean_NO3 + mean_NH4 + GDP_PER + offset(log(Population_20)),
                              data = Final_dat_mort1, expnms = exposures, family = poisson, q = 4)
Mort_qgcomp2 <- qgcomp.boot.quasi(Deaths_20 ~ spatialid + YEAR + Winter_temp + Summer_temp + Winter_sd + 
                                    Summer_sd + mean_EC + mean_OC + mean_SO4 + mean_NO3 + mean_NH4 + GDP_PER + offset(log(Population_20)),
                                  data = Final_dat_mort1, expnms = exposures, family = poisson, q = 4, B = 100, 
                                  seed = 123,parallel = TRUE,parplan = T)

Mort_qgcomp1
(exp(summary(Mort_qgcomp2)$coefficients["psi1", "Estimate"])-1)*100
(exp(summary(Mort_qgcomp2)$coefficients["psi1", "Lower CI"])-1)*100
(exp(summary(Mort_qgcomp2)$coefficients["psi1", "Upper CI"])-1)*100  # 3.58%(1.84%, 5.36%)   08-17年结果

# the coefficients and bootstrap confidence intervals of psi are not different between poisson and quasipoisson distribution,
# so that the following weights were estimated under the assumption on poisson distribution. However, the confidence intervals
# for components in qgcomp were also estimated under the assumption of quasi-poisson distribution
# Mort_qgcomp3 <- qgcomp.boot.quasi(Deaths_20 ~ STATE + YEAR + Winter_temp + Summer_temp + Winter_sd + 
#                                    Summer_sd + mean_EC + mean_OC + mean_SO4 + mean_NO3 + mean_NH4,
#                                  data = Final_dat, expnms = exposures, family = quasipoisson, q = 4, B = 200, seed = 123)

Mort_res <- as.data.frame(matrix(NA, nrow = 5, ncol = 6)) %>% setNames(c("Component","Weight","IR","IR_L","IR_U","P"))
Mort_res$Component <- c("Elemental carbon","Organic carbon","Sulfate","Nitrate","Ammonium")
Mort_res$Weight <- c(-0.893,-0.107,0.519,0.194,0.287)
Mort_res$id <- 1:5
tiff("weight_plot08_17.tiff", res = 900, width = 8, height = 3, units = "in")
ggplot(data = Mort_res, aes(x = id)) +
  geom_rect(aes(xmin = 5-0.3, xmax = 5+0.3,ymin=0,ymax=0.519),fill = "#ef1828",color = "black") +
  geom_rect(aes(xmin = 4-0.3, xmax = 4+0.3,ymin=0,ymax=0.287),fill = "#ef1828",color = "black") +
  geom_rect(aes(xmin = 3-0.3, xmax = 3+0.3,ymin=0,ymax=0.194),fill = "#ef1828",color = "black") +
  geom_rect(aes(xmin = 2-0.3, xmax = 2+0.3,ymin=-0.107,ymax=0),fill = "#00bdcd",color = "black") +
  geom_rect(aes(xmin = 1-0.3, xmax = 1+0.3,ymin=-0.893,ymax=0),fill = "#00bdcd",color = "black") +
  geom_hline(yintercept = 0, lty = 2, color = "#911310") +
  annotate("text", x = 5, y = -0.9, label = "Positive coefficients: 0.053", hjust = 0, size = 4) +
  annotate("text", x = 4.6, y = -0.9, label = "Negative coefficients: -0.018", hjust = 0, size = 4) +
  scale_x_continuous(name = "",breaks = c(1,2,3,4,5),labels = c("Elemental carbon","Organic carbon",
                                                                "Nitrate","Ammonium","Sulfate"))+
  scale_y_continuous(name = "Weight", breaks = c(-0.9,-0.5,0,0.3,0.6),limits = c(-0.95,0.65), expand = c(0,0))+
  coord_flip()+
  theme_bw()+
  theme(axis.text.x = element_text(size = 12),axis.text.y = element_text(size = 12),
        axis.title.x = element_text(size = 13))
dev.off()

## nonlinear effect
{
  library(dlnm)
  #pm2.5
  PM25_ns <- onebasis(Final_dat_mort1$PM_25, fun="ns", df=3)
  PMfit <- glm(Deaths_20 ~ PM25_ns + spatialid + YEAR + Winter_temp +
                 Summer_temp + Winter_sd + Summer_sd + GDP_PER,
               family = quasipoisson, data=Final_dat_mort1, offset = log(Final_dat_mort1$Population_20))
  pred.PM <- crosspred(PM25_ns,PMfit,by=0.2,cumul = F,cen = 6.6)
  
  fitPM <- pred.PM[["matRRfit"]]%>% data.frame()%>%
    mutate(.,PM25=seq(3.2,17,by=0.2))
  lowPM <- pred.PM[["matRRlow"]]%>% data.frame()%>%
    mutate(.,PM25=seq(3.2,17,by=0.2))
  highPM <- pred.PM[["matRRhigh"]]%>% data.frame()%>%
    mutate(.,PM25=seq(3.2,17,by=0.2))
  histPM <- data.frame(fitPM,lowPM[,1],highPM[,1]) %>% setNames(c("RR","PM25","RRLow","RRHigh"))
  # EC
  EC_ns <- onebasis(Final_dat_mort1$mean_EC, fun="ns", df=3)
  ECfit <- glm(Deaths_20 ~ EC_ns + spatialid + YEAR + Winter_temp +
                 Summer_temp + Winter_sd + Summer_sd + GDP_PER,
               family = quasipoisson, data=Final_dat_mort1, offset = log(Final_dat_mort1$Population_20))
  pred.EC <- crosspred(EC_ns,ECfit,by=0.02,cumul = F,cen = 0.24)
  
  fitEC <- pred.EC[["matRRfit"]]%>% data.frame()%>%
    mutate(.,EC25=seq(0.06,1.64,by=0.02))
  lowEC <- pred.EC[["matRRlow"]]%>% data.frame()%>%
    mutate(.,EC25=seq(0.06,1.64,by=0.02))
  highEC <- pred.EC[["matRRhigh"]]%>% data.frame()%>%
    mutate(.,EC25=seq(0.06,1.64,by=0.02))
  histEC <- data.frame(fitEC,lowEC[,1],highEC[,1]) %>% setNames(c("RR","EC25","RRLow","RRHigh"))
  # OC
  OC_ns <- onebasis(Final_dat_mort1$mean_OC, fun="ns", df=3)
  OCfit <- glm(Deaths_20 ~ OC_ns + spatialid + YEAR + Winter_temp +
                 Summer_temp + Winter_sd + Summer_sd + GDP_PER,
               family = quasipoisson, data=Final_dat_mort1, offset = log(Final_dat_mort1$Population_20))
  pred.OC <- crosspred(OC_ns,OCfit,by=0.04,cumul = F,cen = 2.05)
  
  fitOC <- pred.OC[["matRRfit"]]%>% data.frame()%>%
    mutate(.,OC25=seq(0.45,3.21,by=0.04))
  lowOC <- pred.OC[["matRRlow"]]%>% data.frame()%>%
    mutate(.,OC25=seq(0.45,3.21,by=0.04))
  highOC <- pred.OC[["matRRhigh"]]%>% data.frame()%>%
    mutate(.,OC25=seq(0.45,3.21,by=0.04))
  histOC <- data.frame(fitOC,lowOC[,1],highOC[,1]) %>% setNames(c("RR","OC25","RRLow","RRHigh"))
  # SO4
  SO4_ns <- onebasis(Final_dat_mort1$mean_SO4, fun="ns", df=3)
  SO4fit <- glm(Deaths_20 ~ SO4_ns + spatialid + YEAR + Winter_temp +
                  Summer_temp + Winter_sd + Summer_sd + GDP_PER,
                family = quasipoisson, data=Final_dat_mort1, offset = log(Final_dat_mort1$Population_20))
  pred.SO4 <- crosspred(SO4_ns,SO4fit,by=0.12,cumul = F,cen = 1.16)
  
  fitSO4 <- pred.SO4[["matRRfit"]]%>% data.frame()%>%
    mutate(.,SO425=seq(0.2,4.88,by=0.12))
  lowSO4 <- pred.SO4[["matRRlow"]]%>% data.frame()%>%
    mutate(.,SO425=seq(0.2,4.88,by=0.12))
  highSO4 <- pred.SO4[["matRRhigh"]]%>% data.frame()%>%
    mutate(.,SO425=seq(0.2,4.88,by=0.12))
  histSO4 <- data.frame(fitSO4,lowSO4[,1],highSO4[,1]) %>% setNames(c("RR","SO425","RRLow","RRHigh"))
  # NO3
  NO3_ns <- onebasis(Final_dat_mort1$mean_NO3, fun="ns", df=3)
  NO3fit <- glm(Deaths_20 ~ NO3_ns + spatialid + YEAR + Winter_temp +
                  Summer_temp + Winter_sd + Summer_sd + GDP_PER,
                family = quasipoisson, data=Final_dat_mort1, offset = log(Final_dat_mort1$Population_20))
  pred.NO3 <- crosspred(NO3_ns,NO3fit,by=0.1,cumul = F,cen = 0.1)
  
  fitNO3 <- pred.NO3[["matRRfit"]]%>% data.frame()%>%
    mutate(.,NO325=seq(0.1,2.8,by=0.1))
  lowNO3 <- pred.NO3[["matRRlow"]]%>% data.frame()%>%
    mutate(.,NO325=seq(0.1,2.8,by=0.1))
  highNO3 <- pred.NO3[["matRRhigh"]]%>% data.frame()%>%
    mutate(.,NO325=seq(0.1,2.8,by=0.1))
  histNO3 <- data.frame(fitNO3,lowNO3[,1],highNO3[,1]) %>% setNames(c("RR","NO325","RRLow","RRHigh"))
  # NH4
  NH4_ns <- onebasis(Final_dat_mort1$mean_NH4, fun="ns", df=3)
  NH4fit <- glm(Deaths_20 ~ NH4_ns + spatialid + YEAR + Winter_temp +
                  Summer_temp + Winter_sd + Summer_sd + GDP_PER,
                family = quasipoisson, data=Final_dat_mort1, offset = log(Final_dat_mort1$Population_20))
  pred.NH4 <- crosspred(NH4_ns,NH4fit,by=0.08,cumul = F,cen = 0.1)
  
  fitNH4 <- pred.NH4[["matRRfit"]]%>% data.frame()%>%
    mutate(.,NH425=seq(0.1,2.1,by=0.08))
  lowNH4 <- pred.NH4[["matRRlow"]]%>% data.frame()%>%
    mutate(.,NH425=seq(0.1,2.1,by=0.08))
  highNH4 <- pred.NH4[["matRRhigh"]]%>% data.frame()%>%
    mutate(.,NH425=seq(0.1,2.1,by=0.08))
  histNH4 <- data.frame(fitNH4,lowNH4[,1],highNH4[,1]) %>% setNames(c("RR","NH425","RRLow","RRHigh"))
  
  
  PM_plot <- ggplot(histPM, aes(x=PM25)) + 
    geom_line(aes(y=RR), color="black")+ 
    geom_ribbon(aes(ymin=RRLow,ymax=RRHigh),alpha=0.2)+
    geom_hline(yintercept = 1, color = "red", lty = 2) + 
    scale_x_continuous(name = expression(paste(PM[2.5]," (",mu,"g/",m^3,")")), breaks = c(3,6,9,12,15)) +
    theme_bw()
  EC_plot <- ggplot(histEC, aes(x=EC25)) + 
    geom_line(aes(y=RR), color="black")+ 
    geom_ribbon(aes(ymin=RRLow,ymax=RRHigh),alpha=0.2)+
    geom_hline(yintercept = 1, color = "red", lty = 2) + 
    scale_x_continuous(name = expression(paste("Elemental carbon (",mu,"g/",m^3,")")),breaks = c(0.08,0.4,0.72,1.04,1.36)) +
    theme_bw()
  OC_plot <- ggplot(histOC, aes(x=OC25)) + 
    geom_line(aes(y=RR), color="black")+ 
    geom_ribbon(aes(ymin=RRLow,ymax=RRHigh),alpha=0.2)+
    geom_hline(yintercept = 1, color = "red", lty = 2) + 
    scale_x_continuous(name = expression(paste("Organic carbon (",mu,"g/",m^3,")")),breaks = c(0.6,1.2,1.8,2.4,3)) +
    theme_bw()
  SO4_plot <- ggplot(histSO4, aes(x=SO425)) + 
    geom_line(aes(y=RR), color="black")+ 
    geom_ribbon(aes(ymin=RRLow,ymax=RRHigh),alpha=0.2)+
    geom_hline(yintercept = 1, color = "red", lty = 2) + 
    scale_x_continuous(name = expression(paste("Sulfate (",mu,"g/",m^3,")")),breaks = c(0.3,1.2,2.1,3.0,3.9,4.8)) +
    theme_bw()
  NO3_plot <- ggplot(histNO3, aes(x=NO325)) + 
    geom_line(aes(y=RR), color="black")+ 
    geom_ribbon(aes(ymin=RRLow,ymax=RRHigh),alpha=0.2)+
    geom_hline(yintercept = 1, color = "red", lty = 2) + 
    scale_x_continuous(name = expression(paste("Nitrate (",mu,"g/",m^3,")")),breaks = c(0.1,0.6,1.1,1.6,2.1,2.6)) +
    theme_bw()
  NH4_plot <- ggplot(histNH4, aes(x=NH425)) + 
    geom_line(aes(y=RR), color="black")+ 
    geom_ribbon(aes(ymin=RRLow,ymax=RRHigh),alpha=0.2)+
    geom_hline(yintercept = 1, color = "red", lty = 2) + 
    scale_x_continuous(name = expression(paste("Ammonium (",mu,"g/",m^3,")")),breaks = c(0.1,0.5,0.9,1.3,1.7,2.1)) +
    theme_bw()
  
  png("single_nonlinear.png",res = 900, width = 6, height = 6, units = "in")
  ggpubr::ggarrange(PM_plot,EC_plot,OC_plot,SO4_plot,NO3_plot,NH4_plot, ncol = 2, nrow = 3)
  dev.off()
  
  Mort_qgcomp_ns1 <- qgcomp.boot.quasi(Deaths_20 ~ spatialid + YEAR + Winter_temp + Summer_temp + Winter_sd + 
                                         Summer_sd + ns(mean_EC) + ns(mean_OC) + ns(mean_SO4) + ns(mean_NO3) + ns(mean_NH4) + GDP_PER + offset(log(Population_20)),
                                       data = Final_dat_mort1, expnms = exposures, family = poisson, q = 4, 
                                       B = 100, seed = 123, degree = 3,parallel = TRUE,parplan = T)
  Mort_qgcomp_nsplot1 <- plot(Mort_qgcomp_ns1, suppressprint = T, pointwiseref = 4)
  
  png("mixture_nonlinear08_17.png",res = 900, width = 6, height = 4, units = "in")
  Mort_qgcomp_nsplot1
  dev.off()
}



## modification effect
# temperature modification
Mort_ST_qgcomp1 <- qgcomp.modi.noboot(Deaths_20 ~ spatialid + YEAR + Winter_temp + Summer_temp + Winter_sd + GDP_PER + 
                                        Summer_sd + mean_EC + mean_OC + mean_SO4 + mean_NO3 + mean_NH4 + offset(log(Population_20)),
                                      data = Final_dat_mort1, 
                                      expnms = exposures, emmvar = "Summer_temp", family = poisson(), q = 4)
Mort_ST_qgcomp2 <- qgcomp.modi.boot(Deaths_20 ~ spatialid + YEAR + Winter_temp + Summer_temp + Winter_sd + GDP_PER + 
                                      Summer_sd + mean_EC + mean_OC + mean_SO4 + mean_NO3 + mean_NH4 + offset(log(Population_20)),
                                    data = Final_dat_mort1, 
                                    expnms = exposures, emmvar = "Summer_temp", family = poisson(), q = 4,
                                    B=100, seed = 123,parallel = T, parplan = T)
getstrateffects(Mort_ST_qgcomp1, emmval=23.92) # 0.039374   0.003406 0.032698  0.04605   11.56 < 2.2e-16
getstrateffects(Mort_ST_qgcomp1, emmval=23.92-1) # 0.0444507  0.0034613 0.037667 0.051235  12.842 < 2.2e-16
getstrateffects(Mort_ST_qgcomp1, emmval=23.92-2) # 0.0495274  0.0036005 0.042471 0.056584  13.756 < 2.2e-16
getstrateffects(Mort_ST_qgcomp1, emmval=23.92-3) # 0.0546041  0.0038144 0.047128  0.06208  14.315 < 2.2e-16
getstrateffects(Mort_ST_qgcomp1, emmval=23.92+1) # 0.0342973  0.0034386 0.027558 0.041037  9.9741 < 2.2e-16
getstrateffects(Mort_ST_qgcomp1, emmval=23.92+2) # 0.0292206  0.0035567  0.02225 0.036192  8.2157 2.22e-16
getstrateffects(Mort_ST_qgcomp1, emmval=23.92+3) # 0.0241439  0.0037522  0.01679 0.031498  6.4346 1.238e-10

Mort_WT_qgcomp1 <- qgcomp.modi.noboot(Deaths_20 ~ spatialid + YEAR + Winter_temp + Summer_temp + Winter_sd + GDP_PER + 
                                        Summer_sd + mean_EC + mean_OC + mean_SO4 + mean_NO3 + mean_NH4 + offset(log(Population_20)),
                                      data = Final_dat_mort1, 
                                      expnms = exposures, emmvar = "Winter_temp", family = poisson(), q = 4)
Mort_WT_qgcomp2 <- qgcomp.modi.boot(Deaths_20 ~ spatialid + YEAR + Winter_temp + Summer_temp + Winter_sd + GDP_PER + 
                                      Summer_sd + mean_EC + mean_OC + mean_SO4 + mean_NO3 + mean_NH4 + offset(log(Population_20)),
                                    data = Final_dat_mort1, 
                                    expnms = exposures, emmvar = "Winter_temp", family = poisson(), q = 4,
                                    B=100, seed = 123,parallel = T, parplan = T)
getstrateffects(Mort_WT_qgcomp1, emmval=1.68) # 0.0387318 0.0035135 0.031845 0.045618 11.024 < 2.2e-16
getstrateffects(Mort_WT_qgcomp1, emmval=1.68-1) # 0.04000913 0.0036154 0.033005 0.047177 11.089 < 2.2e-16
getstrateffects(Mort_WT_qgcomp1, emmval=1.68-2) # 0.0414509  0.0037432 0.034114 0.048788  11.073 < 2.2e-16
getstrateffects(Mort_WT_qgcomp1, emmval=1.68-3) # 0.0428105  0.0038946 0.035177 0.050444  10.992 < 2.2e-16
getstrateffects(Mort_WT_qgcomp1, emmval=1.68+1) # 0.037372 0.003440 0.03063 0.044144 10.864 < 2.2e-16
getstrateffects(Mort_WT_qgcomp1, emmval=1.68+2) # 0.0360126  0.0033967 0.029355  0.04267  10.602 < 2.2e-16
getstrateffects(Mort_WT_qgcomp1, emmval=1.68+3) # 0.0346530  0.0033846 0.028019 0.041287  10.238 < 2.2e-16
# seasonal temperature plot
{
  Tempmodi_dat = as.data.frame(matrix(NA,nrow = 28, ncol = 7)) %>% setNames(c("Season","model","temp","Tempid","ir","ir_l","ir_u"))
  Tempmodi_dat$Season = rep(c("Summer","Winter"), each = 14)
  Tempmodi_dat$model = c(rep("TWFE",times=7),rep("QGC-TWFE",times=7),rep("TWFE",times=7),rep("QGC-TWFE",times=7))
  Tempmodi_dat$temp = rep(c(1,2,3,4,5,6,7),times=4)
  Tempmodi_dat$temp = factor(Tempmodi_dat$temp, levels = c(1,2,3,4,5,6,7),labels = c("T","T-1℃","T-2℃","T-3℃","T+1℃","T+2℃","T+3℃"))
  Tempmodi_dat$Tempid = c(0.7,0.8,0.9,1.0,1.1,1.2,1.3,0.7,0.8,0.9,1.0,1.1,1.2,1.3,
                          1.7,1.8,1.9,2.0,2.1,2.2,2.3,1.7,1.8,1.9,2.0,2.1,2.2,2.3)
  Tempmodi_dat$ir = c(3.070294,
                      3.325386,3.581110,3.837467,
                      2.815832,2.561997,2.308790,
                      (exp(0.039374)-1)*100, # st
                      (exp(0.0444507)-1)*100,(exp(0.0495274)-1)*100,(exp(0.0546041)-1)*100, #st-
                      (exp(0.0342973)-1)*100,(exp(0.0292206)-1)*100, (exp(0.0241439)-1)*100,#st+
                      3.604455,
                      3.882834,4.161960,4.441836,
                      3.326823,3.049934,2.773788,
                      (exp(0.0387318)-1)*100, # wt
                      (exp(0.04000913)-1)*100,(exp(0.0414509)-1)*100,(exp(0.0428105)-1)*100, # wt-
                      (exp(0.037372)-1)*100,(exp(0.0360126)-1)*100,(exp(0.0346530)-1)*100) # wt+
  Tempmodi_dat$ir_l = c(2.38143,
                        2.660761,2.910318,3.130471,
                        2.075743,1.749207,1.407393,
                        (exp(0.032698)-1)*100, # st
                        (exp(0.037667)-1)*100,(exp(0.042471)-1)*100,(exp(0.047128)-1)*100, # st-
                        (exp(0.027558)-1)*100,(exp(0.02225)-1)*100,(exp(0.01679)-1)*100, # st+
                        2.938499,
                        3.204802,3.463365 ,3.714784,
                        2.664132,2.381710,2.091586,
                        (exp(0.031845)-1)*100, # wt
                        (exp(0.033005)-1)*100,(exp(0.034114)-1)*100,(exp(0.035177)-1)*100, # wt-
                        (exp(0.03063)-1)*100,(exp(0.029355)-1)*100,(exp(0.028019)-1)*100) # wt+
  Tempmodi_dat$ir_u = c(3.763793,
                        3.994315,4.256274,4.549309,
                        3.561286,3.381281,3.218199,
                        (exp(0.04605)-1)*100, # st
                        (exp(0.051235)-1)*100,(exp(0.056584)-1)*100,(exp(0.06208)-1)*100, # st-
                        (exp(0.041037)-1)*100,(exp(0.036192)-1)*100,(exp(0.031498)-1)*100, # st+
                        4.274719,
                        4.56532,4.865272,5.173986,
                        3.993791,3.722519,3.460547,
                        (exp(0.045618)-1)*100, # wt
                        (exp(0.047177)-1)*100,(exp(0.048788)-1)*100,(exp(0.050444)-1)*100, # wt-
                        (exp(0.044144)-1)*100,(exp(0.04267)-1)*100,(exp(0.041287)-1)*100) # wt+
}

compareEffect = function(beta,var){
  beta1 = beta[1]
  beta2 = beta[2]
  var1 = var[1]
  var2 = var[2]
  pooled_beta = ((beta1/var1)+(beta2/var2))/(1/var1+1/var2)
  CochranQ = ((beta1-pooled_beta)^2)/var1 + ((beta2-pooled_beta)^2)/var2
  P = pchisq(CochranQ,df = 1, lower.tail = F)
  return(P)
}
{
  mst_beta1 = log((3.070294/100)+1)
  mst_beta2 = log((3.325386/100)+1);mst_beta3 = log((3.581110/100)+1);mst_beta4 = log((3.837467/100)+1)
  mst_beta5 = log((2.815832/100)+1);mst_beta6 = log((2.561997/100)+1);mst_beta7 = log((2.308790/100)+1)
  mst_beta1_var = ((log((3.763793/100)+1)-log((2.38143/100)+1))/(2*1.96))^2
  mst_beta2_var = ((log((3.994315/100)+1)-log((2.660761/100)+1))/(2*1.96))^2
  mst_beta3_var = ((log((4.256274/100)+1)-log((2.910318/100)+1))/(2*1.96))^2
  mst_beta4_var = ((log((4.549309/100)+1)-log((3.130471/100)+1))/(2*1.96))^2
  mst_beta5_var = ((log((3.561286/100)+1)-log((2.075743/100)+1))/(2*1.96))^2
  mst_beta6_var = ((log((3.381281/100)+1)-log((1.749207/100)+1))/(2*1.96))^2
  mst_beta7_var = ((log((3.218199/100)+1)-log((1.407393/100)+1))/(2*1.96))^2
  mwt_beta1 = log((3.604455/100)+1)
  mwt_beta2 = log((3.882834/100)+1);mwt_beta3 = log((4.161960/100)+1);mwt_beta4 = log((4.441836/100)+1)
  mwt_beta5 = log((3.326823/100)+1);mwt_beta6 = log((3.049934/100)+1);mwt_beta7 = log((2.773788/100)+1)
  mwt_beta1_var = ((log((4.274719/100)+1)-log((2.938499/100)+1))/(2*1.96))^2
  mwt_beta2_var = ((log((4.56532/100)+1)-log((3.204802/100)+1))/(2*1.96))^2
  mwt_beta3_var = ((log((4.865272/100)+1)-log((3.463365/100)+1))/(2*1.96))^2
  mwt_beta4_var = ((log((5.173986/100)+1)-log((3.714784/100)+1))/(2*1.96))^2
  mwt_beta5_var = ((log((3.993791/100)+1)-log((2.664132/100)+1))/(2*1.96))^2
  mwt_beta6_var = ((log((3.722519/100)+1)-log((2.381710/100)+1))/(2*1.96))^2
  mwt_beta7_var = ((log((3.460547/100)+1)-log((2.091586/100)+1))/(2*1.96))^2
  compareEffect(beta = c(mst_beta1,mst_beta2), var = c(mst_beta1_var,mst_beta2_var)) 
  compareEffect(beta = c(mst_beta1,mst_beta3), var = c(mst_beta1_var,mst_beta3_var)) 
  compareEffect(beta = c(mst_beta1,mst_beta4), var = c(mst_beta1_var,mst_beta4_var)) 
  compareEffect(beta = c(mst_beta1,mst_beta5), var = c(mst_beta1_var,mst_beta5_var)) 
  compareEffect(beta = c(mst_beta1,mst_beta6), var = c(mst_beta1_var,mst_beta6_var)) 
  compareEffect(beta = c(mst_beta1,mst_beta7), var = c(mst_beta1_var,mst_beta7_var)) 
  compareEffect(beta = c(mwt_beta1,mwt_beta2), var = c(mwt_beta1_var,mwt_beta2_var)) 
  compareEffect(beta = c(mwt_beta1,mwt_beta3), var = c(mwt_beta1_var,mwt_beta3_var)) 
  compareEffect(beta = c(mwt_beta1,mwt_beta4), var = c(mwt_beta1_var,mwt_beta4_var)) 
  compareEffect(beta = c(mwt_beta1,mwt_beta5), var = c(mwt_beta1_var,mwt_beta5_var)) 
  compareEffect(beta = c(mwt_beta1,mwt_beta6), var = c(mwt_beta1_var,mwt_beta6_var)) 
  compareEffect(beta = c(mwt_beta1,mwt_beta7), var = c(mwt_beta1_var,mwt_beta7_var))
}
compareEffect(beta = c(0.039374,0.0444507), var = c(0.003406^2,0.0034613^2)) 
compareEffect(beta = c(0.039374,0.0495274), var = c(0.003406^2,0.0036005^2))
compareEffect(beta = c(0.039374,0.0546041), var = c(0.003406^2,0.0038144^2))
compareEffect(beta = c(0.039374,0.0342973), var = c(0.003406^2,0.0034386^2))
compareEffect(beta = c(0.039374,0.0292206), var = c(0.003406^2,0.0035567^2))
compareEffect(beta = c(0.039374,0.0241439), var = c(0.003406^2,0.0037522^2))

compareEffect(beta = c(0.0387318,0.04000913), var = c(0.0035135^2,0.0036154^2)) 
compareEffect(beta = c(0.0387318,0.0414509), var = c(0.0035135^2,0.0037432^2))
compareEffect(beta = c(0.0387318,0.0428105), var = c(0.0035135^2,0.0038946^2))
compareEffect(beta = c(0.0387318,0.037372), var = c(0.0035135^2,0.003440^2))
compareEffect(beta = c(0.0387318,0.0360126), var = c(0.0035135^2,0.0033967^2))
compareEffect(beta = c(0.0387318,0.0346530), var = c(0.0035135^2,0.0033846^2))

temp.1 = ggplot(Tempmodi_dat[Tempmodi_dat$model == "TWFE",],aes(x=Tempid)) +
  geom_point(aes(y=ir, color = temp)) +
  geom_errorbar(aes(ymin = ir_l, ymax = ir_u, color = temp), width = 0.03, linewidth = 0.8) +
  geom_hline(yintercept = 0, color = "black", lty = 2) +
  scale_x_continuous(name = "", breaks = c(1,2), labels = c("","")) +
  scale_y_continuous(name = "IR%") +
  scale_color_manual(values = colors1) +
  theme_bw()+
  theme(axis.title = element_text(size = 9),
        legend.title = element_blank(),
        legend.background = element_blank())
temp.2 = ggplot(Tempmodi_dat[Tempmodi_dat$model == "QGC-TWFE",],aes(x=Tempid)) +
  geom_point(aes(y=ir, color = temp)) +
  geom_errorbar(aes(ymin = ir_l, ymax = ir_u, color = temp), width = 0.03, linewidth = 0.8) +
  geom_hline(yintercept = 0, color = "black", lty = 2) +
  annotate("text",x=0.9, y = (exp(0.042471)-1)*100 -0.4, label = "*", size = 6)+
  annotate("text",x=1, y = (exp(0.047128)-1)*100 -0.4, label = "*", size = 6)+
  annotate("text",x=1.2, y = (exp(0.02225)-1)*100 -0.4, label = "*", size = 6)+
  annotate("text",x=1.3, y = (exp(0.01679)-1)*100 -0.4, label = "*", size = 6)+
  scale_x_continuous(name = "", breaks = c(1,2), labels = c("Summer","Winter")) +
  scale_y_continuous(name = "IR%") +
  scale_color_manual(values = colors1) +
  theme_bw()+
  theme(axis.title = element_text(size = 9),
        legend.title = element_blank(),
        legend.background = element_blank())
tiff("Temp_modi08_17.tiff",res = 900, width = 4, height = 5, units = "in")
ggpubr::ggarrange(temp.1,temp.2, nrow = 2, align = "v", labels = c("TWFE model","QGC-TWFE model"),
                  font.label = list(size = 10), common.legend = T, legend = "right",
                  label.x = 0.1, label.y = 1.05, hjust = 0) +
  theme(plot.margin = margin(t=0.5,r=0,b=0,l=0, unit = "cm"))
dev.off()
# Poverty level
Mort_poverty_qgcomp1 <- qgcomp.modi.noboot(Deaths_20 ~ spatialid + YEAR + Winter_temp + Summer_temp + Winter_sd + GDP_PER + 
                                             Summer_sd + mean_EC + mean_OC + mean_SO4 + mean_NO3 + mean_NH4 + offset(log(Population_20)),
                                           data = Final_dat_mort1, expnms = exposures, emmvar = "Povertycat", family = poisson(), q = 4)

Mort_poverty_qgcomp2 <- qgcomp.modi.boot(Deaths_20 ~ spatialid + YEAR + Winter_temp + Summer_temp + Winter_sd + GDP_PER + 
                                           Summer_sd + mean_EC + mean_OC + mean_SO4 + mean_NO3 + mean_NH4 + offset(log(Population_20)),
                                         data = Final_dat_mort1, expnms = exposures, emmvar = "Povertycat", family = poisson(), q = 4,
                                         B=100, seed = 123, parallel = TRUE, parplan = T)
##08-17年
#                           Estimate  Std. Error  Lower CI  Upper CI  Z value  Pr(>|z|)
#psi1                      0.0336167  0.0119585  0.010179  0.057055   2.8111   0.004937
#PovertycatGroup2:mixture  0.0028304  0.0151081 -0.026781  0.032442   0.1873   0.851394
# Obesity rate
Mort_obese_qgcomp1 <- qgcomp.modi.noboot(Deaths_20 ~ spatialid + YEAR + Winter_temp + Summer_temp + Winter_sd + GDP_PER + 
                                           Summer_sd + mean_EC + mean_OC + mean_SO4 + mean_NO3 + mean_NH4 + offset(log(Population_20)),
                                         data = Final_dat_mort1[!is.na(Final_dat_mort1$Obesitycat),], expnms = exposures, emmvar = "Obesitycat", family = poisson(), q = 4)
Mort_obese_qgcomp2 <- qgcomp.modi.boot(Deaths_20 ~ spatialid + YEAR + Winter_temp + Summer_temp + Winter_sd + GDP_PER + 
                                         Summer_sd + mean_EC + mean_OC + mean_SO4 + mean_NO3 + mean_NH4 + offset(log(Population_20)),
                                       data = Final_dat_mort1[!is.na(Final_dat_mort1$Obesitycat),], expnms = exposures, emmvar = "Obesitycat", family = poisson(), q = 4,
                                       B=100, seed = 123,parallel = TRUE, parplan = T)
##08-17年
#                           Estimate Std. Error  Lower CI  Upper CI  Z value  Pr(>|z|)
#(Intercept)               4.535167   0.021883  4.492278 4.5780558 207.2502 < 2.2e-16
#psi1                     -0.020623   0.014039 -0.048138 0.0068918  -1.4690    0.1418
#ObesitycatGroup2:mixture  0.094094   0.015226  0.064251 0.1239364   6.1797 6.422e-10
# physical inactivity
Mort_inact_qgcomp1 <- qgcomp.modi.noboot(Deaths_20 ~ spatialid + YEAR + Winter_temp + Summer_temp + Winter_sd + GDP_PER + 
                                           Summer_sd + mean_EC + mean_OC + mean_SO4 + mean_NO3 + mean_NH4 + offset(log(Population_20)),
                                         data = Final_dat_mort1[!is.na(Final_dat_mort1$Inactivitycat),], expnms = exposures, emmvar = "Inactivitycat", family = poisson(), q = 4)
Mort_inact_qgcomp2 <- qgcomp.modi.boot(Deaths_20 ~ spatialid + YEAR + Winter_temp + Summer_temp + Winter_sd + GDP_PER + 
                                         Summer_sd + mean_EC + mean_OC + mean_SO4 + mean_NO3 + mean_NH4 + offset(log(Population_20)),
                                       data = Final_dat_mort1[!is.na(Final_dat_mort1$Inactivitycat),], expnms = exposures, emmvar = "Inactivitycat", family = poisson(), q = 4,
                                       B=100, seed = 123,parallel = TRUE, parplan = T)
##08-17年
#                             Estimate    Std. Error  
#psi1                         0.03569786  0.01243682  
#InactivitycatGroup2:mixture  0.00057577  0.01467023 
# aging degree
Mort_aging_qgcomp1 <- qgcomp.modi.noboot(Deaths_20 ~ spatialid + YEAR + Winter_temp + Summer_temp + Winter_sd + GDP_PER + 
                                           Summer_sd + mean_EC + mean_OC + mean_SO4 + mean_NO3 + mean_NH4 + offset(log(Population_20)),
                                         data = Final_dat_mort1, expnms = exposures, emmvar = "Agingcat", family = poisson(), q = 4)
Mort_aging_qgcomp2 <- qgcomp.modi.boot(Deaths_20 ~ spatialid + YEAR + Winter_temp + Summer_temp + Winter_sd + GDP_PER + 
                                         Summer_sd + mean_EC + mean_OC + mean_SO4 + mean_NO3 + mean_NH4 + offset(log(Population_20)),
                                       data = Final_dat_mort1, expnms = exposures, emmvar = "Agingcat", family = poisson(), q = 4,
                                       B=100, seed = 123,parallel = TRUE, parplan = T)
##08-17年
#                        Estimate   Std. Error  
#psi1                    0.2879063  0.0098316  
#AgingcatAging2:mixture -0.4276947  0.0145413
# health insurance coverage
Mort_insure_qgcomp1 <- qgcomp.modi.noboot(Deaths_20 ~ spatialid + YEAR + Winter_temp + Summer_temp + Winter_sd + GDP_PER + 
                                            Summer_sd + mean_EC + mean_OC + mean_SO4 + mean_NO3 + mean_NH4 + offset(log(Population_20)),
                                          data = Final_dat_mort1[!is.na(Final_dat_mort1$Insurecat),], expnms = exposures, emmvar = "Insurecat", family = poisson(), q = 4)
Mort_insure_qgcomp2 <- qgcomp.modi.boot(Deaths_20 ~ spatialid + YEAR + Winter_temp + Summer_temp + Winter_sd + GDP_PER + 
                                          Summer_sd + mean_EC + mean_OC + mean_SO4 + mean_NO3 + mean_NH4 + offset(log(Population_20)),
                                        data = Final_dat_mort1[!is.na(Final_dat_mort1$Insurecat),], expnms = exposures, emmvar = "Insurecat", family = poisson(), q = 4,
                                        B=100, seed = 123,parallel = TRUE, parplan = T)
##08-17年
#                         Estimate   Std. Error  
#psi1                     0.0062715  0.0129745  
#InsurecatGroup2:mixture  0.0534909  0.0143125

Modif_weight <- as.data.frame(matrix(NA, nrow = 5, ncol = 11)) %>% 
  setNames(c("Component","Lower poverty","Higher poverty","Lower obesity","Higher obesity","Lower Inactivity","Higher Inactivity","Lower insurance","Higher insurance","Aging","Aged"))
Modif_weight$Component <- c("Elemental carbon","Organic carbon","Sulfate","Nitrate","Ammonium")
Modif_weight$`Lower poverty` <- c(-0.712,-0.288,0.602,0.179,0.219)
Modif_weight$`Higher poverty` <- c(-1,0.0267,0.5703,0.0804,0.3227)
Modif_weight$`Lower obesity` <- c(-0.827,-0.173,0.554,0.154,0.293)
Modif_weight$`Higher obesity` <- c(-0.9391,-0.0609,0.6270,0.0962,0.2768)
Modif_weight$`Lower Inactivity` <- c(-0.772,-0.228,0.588,0.193,0.218)
Modif_weight$`Higher Inactivity` <- c(-1,0.014,0.560,0.082,0.344)
Modif_weight$`Lower insurance` <- c(-0.842,-0.158,0.211,0.470,0.319)
Modif_weight$`Higher insurance` <- c(-0.5860,-0.0887,-0.3254,0.829,0.171)
Modif_weight$Aging <- c(-1,0.0383,0.4342,0.2587,0.2688)
Modif_weight$Aged <- c(-0.835,-0.165,0.309,0.354,0.337)
Modif_weight <- Modif_weight %>% gather(., key = "group", value = "weight",-Component) %>%
  mutate(group = case_when(group == "Aged" ~ 10, group == "Aging" ~ 9, 
                           group == "Higher insurance" ~ 8, group == "Lower insurance" ~ 7,
                           group == "Higher Inactivity" ~ 6, group == "Lower Inactivity" ~ 5, 
                           group == "Higher obesity" ~ 4, group == "Lower obesity" ~ 3, 
                           group == "Higher poverty" ~ 2, group == "Lower poverty" ~ 1),
         component = case_when(Component == "Elemental carbon" ~ 5,
                               Component == "Organic carbon" ~ 4,
                               Component == "Sulfate" ~ 1,
                               Component == "Nitrate" ~ 2,
                               Component == "Ammonium" ~ 3))

Inter_res1 <- as.data.frame(matrix(NA, nrow = 5, ncol = 4)) %>%
  setNames(c("group","IR","IR_L","IR_U"))
Inter_res1$group <- c(1,2,3,4,5)
Inter_res1$IR <- c(0.169,-0.063,0.338,-0.089,-0.273)
Inter_res1$IR_L <- c(-0.095,-0.251,0.15,-0.344,-0.521)
Inter_res1$IR_U <- c(0.434,0.124,0.527,0.168,-0.023)

Inter_res2 <- as.data.frame(matrix(NA, nrow = 5, ncol = 4)) %>%
  setNames(c("group","IR","IR_L","IR_U"))
Inter_res2$group <- c(1,2,3,4,5)
Inter_res2$IR <- c(exp(0.0034468)-1,
                   exp(0.00005966)-1,
                   exp(0.0044611)-1,
                   exp(0.0061610)-1,
                   exp(-0.0017182)-1)
Inter_res2[1,3:4] <- c(exp(0.0034468-1.96*0.0151081)-1,exp(0.0034468+1.96*0.0151081)-1)
Inter_res2[2,3:4] <- c(exp(0.00005966-1.96*0.015226)-1,exp(0.00005966+1.96*0.015226)-1)
Inter_res2[3,3:4] <- c(exp(0.0044611-1.96*0.01467023)-1,exp(0.0044611+1.96*0.01467023)-1)
Inter_res2[4,3:4] <- c(exp(0.0061610-1.96*0.0143125)-1,exp(0.0061610+1.96*0.0143125)-1)
Inter_res2[5,3:4] <- c(exp(-0.0017182-1.96*0.0145413)-1,exp(-0.0017182+1.96*0.0145413)-1)
Inter_res2$IR <- 100*Inter_res2$IR
Inter_res2$IR_L <- 100*Inter_res2$IR_L
Inter_res2$IR_U <- 100*Inter_res2$IR_U


weightplot <- ggplot(Modif_weight, aes(x = component, y = group, fill = weight)) +
  geom_raster() +
  geom_hline(yintercept = c(0.5,2.5,4.5,6.5,8.5,10.5),linewidth = 1.8, color = "white") +
  geom_hline(yintercept = c(1.5,3.5,5.5,7.5,9.5), linewidth = 0.8, color = "white") +
  geom_vline(xintercept = c(1.5,2.5,3.5,4.5), linewidth = 0.8, color = "white") +
  
  scale_fill_gradient2(low = "#00b0eb",high = "#e20612",mid = "grey97") +
  scale_x_continuous(name = "", breaks = c(1,2,3,4,5), labels = c("SO4","NO3","NH4","OC","EC"),
                     expand = c(0,0)) +
  scale_y_continuous(name = "", limits = c(0.5,10.5), breaks = c(1.5,3.5,5.5,7.5,9.5), 
                     labels = c("Poverty","Obesity","Inactivity","Insurance","Aging"),
                     sec.axis = sec_axis(~.,
                                         breaks = c(1,2,3,4,5,6,7,8,9,10),
                                         labels = c(expression(paste(Lower^ref,"",sep = "")),"Higher",
                                                    expression(paste(Lower^ref,"",sep = "")),"Higher",
                                                    expression(paste(Lower^ref,"",sep = "")),"Higher",
                                                    expression(paste(Lower^ref,"",sep = "")),"Higher",
                                                    expression(paste(Aging^ref,"",sep = "")),"Aged")),
                     expand = c(0,0)) +
  coord_flip()+
  theme_bw() +
  theme(panel.grid.major=element_blank(), #去掉网格线中的竖线
        panel.grid.minor=element_blank(), #去掉网格线中的横线,
        axis.text.y.left = element_text(angle = 90, hjust = 0.5, vjust = 4, margin = unit(c(0,0,0,-0.5),"cm")),
        panel.border = element_blank(),
        axis.text.y = element_text(size = 10),
        axis.text.x.bottom = element_text(size = 10))

interactplot1 <- ggplot(Inter_res1, aes(x = group)) +
  geom_point(aes(y = IR),size = 2)+
  geom_errorbar(aes(ymin = IR_L, ymax = IR_U), width = 0.1, linewidth = 0.8) +
  geom_hline(yintercept = 0, lty = 2, color = "red") +
  annotate("text", x = 3, y = 0.1, label = "※", size = 5)+
  annotate("text", x = 5, y = -0.571, label = "※", size = 5)+
  annotate("text", x = 0.6, y = -0.55, label = "※: significant difference",size=4.5,hjust=0)+
  scale_x_continuous(name = "", limits = c(0.5,5.5),breaks = c(1,2,3,4,5), labels = c("Poverty","Obesity","Inactivity","Insurance","Aging"),expand = c(0,0)) +
  scale_y_continuous(name = "", limits = c(-0.6,0.6), breaks = c(-0.6,-0.3,0,0.3,0.6), labels = c("-0.6%","-0.3%","0.0%","0.3%","0.6%"))+
  theme_bw()+
  theme(axis.text.y = element_text(size = 10),
        axis.text.x.bottom = element_text(size = 10))
interactplot2 <- ggplot(Inter_res2, aes(x = group)) +
  geom_point(aes(y = IR),size = 2)+
  geom_errorbar(aes(ymin = IR_L, ymax = IR_U), width = 0.1, linewidth = 0.8) +
  geom_hline(yintercept = 0, lty = 2, color = "red") +
  scale_x_continuous(name = "", limits = c(0.5,5.5),breaks = c(1,2,3,4,5), labels = c("Poverty","Obesity","Inactivity","Insurance","Aging"),expand = c(0,0)) +
  scale_y_continuous(name = "",limits = c(-3,4.5), breaks = c(-3,-1.5,0,1.5,3,4.5), labels = c("-3.0%","-1.5%","0.0%","1.5%","3.0%","4.5%"))+
  theme_bw()+
  theme(axis.text.y = element_text(size = 10),
        axis.text.x.bottom = element_text(size = 10))

sociomodi.1 = ggpubr::ggarrange(interactplot1,interactplot2, nrow = 2, labels = c("A","B"),
                                label.x = 0.15,label.y = 1.08, font.label = list(size = 12), hjust = 0)+
  theme(plot.margin = margin(t=0.5,r=0.5,b=0,l=0,unit = "cm"))

tiff("Modification_effect08_17.tiff", res = 900, width = 11, height = 5.5,units = "in")
ggpubr::ggarrange(sociomodi.1,weightplot,widths = c(1.15,1.85), common.legend = F, labels = "C",
                  label.x = 0.94, label.y = 1.,font.label = list(size = 12))+
  theme(plot.margin = margin(t=0.5,r=0,b=0,l=0,unit = "cm"))
dev.off()




#### sensitivity analysis ####
sen_qgcomp_q5 <- qgcomp.boot.quasi(Deaths_20 ~ spatialid + YEAR + Winter_temp + Summer_temp + Winter_sd + 
                                     Summer_sd + mean_EC + mean_OC + mean_SO4 + mean_NO3 + mean_NH4 + GDP_PER + offset(log(Population_20)),
                                   data = Final_dat_mort1, expnms = exposures, family = poisson, q = 5, B = 100, seed = 123,parallel = TRUE, parplan = T)
(exp(summary(sen_qgcomp_q5)$coefficients["psi1", "Estimate"])-1)*100
(exp(summary(sen_qgcomp_q5)$coefficients["psi1", "Lower CI"])-1)*100
(exp(summary(sen_qgcomp_q5)$coefficients["psi1", "Upper CI"])-1)*100
sen_qgcomp_q5$pval

sen_qgcomp_q10 <- qgcomp.boot.quasi(Deaths_20 ~ spatialid + YEAR + Winter_temp + Summer_temp + Winter_sd + 
                                      Summer_sd + mean_EC + mean_OC + mean_SO4 + mean_NO3 + mean_NH4 + GDP_PER + offset(log(Population_20)),
                                    data = Final_dat_mort1, expnms = exposures, family = poisson, q = 10, B = 100, seed = 123,parallel = TRUE, parplan = T)
(exp(summary(sen_qgcomp_q10)$coefficients["psi1", "Estimate"])-1)*100
(exp(summary(sen_qgcomp_q10)$coefficients["psi1", "Lower CI"])-1)*100
(exp(summary(sen_qgcomp_q10)$coefficients["psi1", "Upper CI"])-1)*100
sen_qgcomp_q10$pval


sen_gwqs <- gwqs(Deaths_20 ~ wqs + spatialid + YEAR + Winter_temp +
                   Summer_temp + Winter_sd + Summer_sd + GDP_PER + offset(log(Population_20)),
                 data = Final_dat_mort1, mix_name = exposures, b1_pos = T, b1_constr = T, 
                 q = 4, b = 100, validation = 0.3, family = quasipoisson, seed = 123)
round(Calc_IR(sen_gwqs, component = "wqs"),3)

exposures2 <- c("mean_NO3","mean_SO4","mean_NH4")
sen_gwqs_pos <- gwqs(Deaths_20 ~ wqs + spatialid + YEAR + Winter_temp +
                       Summer_temp + Winter_sd + Summer_sd + GDP_PER + offset(log(Population_20)),
                     data = Final_dat_mort1, mix_name = exposures2, b1_pos = T, b1_constr = T, 
                     q = 4, b = 100, validation = 0.3, family = quasipoisson, seed = 123)
round(Calc_IR(sen_gwqs_pos, component = "wqs"),3)

P_dif_tb <- as.data.frame(matrix(NA, nrow = 5, ncol = 3))
colnames(P_dif_tb) <- c("beta","se","group")
P_dif_tb$group <- c("Primary","q5","q10","wqs5","wqs3")
P_dif_tb[1,1] <- Mort_qgcomp2$psi
P_dif_tb[1,2] <- (Mort_qgcomp2$ci[2]-Mort_qgcomp2$ci[1])/(2*1.96)
P_dif_tb[2,1] <- sen_qgcomp_q5$psi
P_dif_tb[2,2] <- (sen_qgcomp_q5$ci[2]-sen_qgcomp_q5$ci[1])/(2*1.96)
P_dif_tb[3,1] <- sen_qgcomp_q10$psi
P_dif_tb[3,2] <- (sen_qgcomp_q10$ci[2]-sen_qgcomp_q10$ci[1])/(2*1.96)
P_dif_tb[4,1] <- sen_gwqs$fit$coefficients['wqs']
P_dif_tb[4,2] <- summary(sen_gwqs)$coefficients['wqs','Std. Error']
P_dif_tb[5,1] <- sen_gwqs_pos$fit$coefficients['wqs']
P_dif_tb[5,2] <- summary(sen_gwqs_pos)$coefficients['wqs','Std. Error']
P_dif_tb$group <- factor(P_dif_tb$group, levels = c("Primary","q5","q10","wqs5","wqs3"))

set.seed(1)
dif_tb <- mvmeta::mvmeta(beta ~ group, S = se^2, data = P_dif_tb, method = "ml")
summary(dif_tb)

# E value
# https://www.evalue-calculator.com/

# assumption checking
# calculate relative change rate
Ec <- Final_dat_mort1 %>% group_by(spatialid) %>% summarize(SO4_Ec = mean(mean_SO4),
                                                            NO3_Ec = mean(mean_NO3),
                                                            NH4_Ec = mean(mean_NH4),
                                                            OC_Ec = mean(mean_OC),
                                                            EC_Ec = mean(mean_EC),
                                                            Wintertemp_Ec = mean(Winter_temp),
                                                            Summertemp_Ec = mean(Summer_temp),
                                                            WSD_Ec = mean(Winter_sd),
                                                            SSD_Ec = mean(Summer_sd),
                                                            GDP_Ec = mean(GDP_PER))
res <- Final_dat_mort1 %>% left_join(.,Ec, by = "spatialid") %>%
  mutate(SO4_RC = 100*(mean_SO4 - SO4_Ec)/SO4_Ec,
         NO3_RC = 100*(mean_NO3 - NO3_Ec)/NO3_Ec,
         NH4_RC = 100*(mean_NH4 - NH4_Ec)/NH4_Ec,
         OC_RC = 100*(mean_OC - OC_Ec)/OC_Ec,
         EC_RC = 100*(mean_EC - EC_Ec)/EC_Ec,
         Wintertemp_RC = 100*(Winter_temp - Wintertemp_Ec)/Wintertemp_Ec,
         Summertemp_RC = 100*(Summer_temp - Summertemp_Ec)/Summertemp_Ec,
         WSD_RC = 100*(Winter_sd - WSD_Ec)/WSD_Ec,
         SSD_RC = 100*(Summer_sd - SSD_Ec)/SSD_Ec,
         GDP_RC = 100*(GDP_PER - GDP_Ec)/GDP_Ec) %>%
  select(c(STATE,YEAR, SO4_RC,NO3_RC,NH4_RC,OC_RC,EC_RC
           ,Wintertemp_RC,Summertemp_RC,WSD_RC,SSD_RC,GDP_RC))

plot_dat <- res %>% group_by(YEAR) %>% summarise(mSO4_RC = mean(SO4_RC),
                                                 mNO3_RC = mean(NO3_RC),
                                                 mNH4_RC = mean(NH4_RC),
                                                 mOC_RC = mean(OC_RC),
                                                 mEC_RC = mean(EC_RC),
                                                 mWintertemp_RC = mean(Wintertemp_RC),
                                                 mSummertemp_RC = mean(Summertemp_RC),
                                                 mWSD_RC = mean(WSD_RC),
                                                 mSSD_RC = mean(SSD_RC),
                                                 mGDP_RC = mean(GDP_RC))
RES <- res %>% left_join(., plot_dat, by = "YEAR") %>%
  mutate(SO4_did = SO4_RC - mSO4_RC,
         NO3_did = NO3_RC - mNO3_RC,
         NH4_did = NH4_RC - mNH4_RC,
         OC_did = OC_RC - mOC_RC,
         EC_did = mEC_RC - EC_RC,
         Wintertemp_did = Wintertemp_RC - mWintertemp_RC,
         Summertemp_did = Summertemp_RC - mSummertemp_RC,
         WSD_did = WSD_RC - mWSD_RC,
         SSD_did = SSD_RC - mSSD_RC,
         GDP_did = GDP_RC - mGDP_RC)

cor.test(RES$SO4_did,RES$Wintertemp_did, method = "pearson")
cor.test(RES$SO4_did,RES$Summertemp_did, method = "pearson")
cor.test(RES$SO4_did,RES$WSD_did, method = "pearson")
cor.test(RES$SO4_did,RES$SSD_did, method = "pearson")
cor.test(RES$SO4_did,RES$GDP_did, method = "pearson")

cor.test(RES$NO3_did,RES$Wintertemp_did, method = "pearson")
cor.test(RES$NO3_did,RES$Summertemp_did, method = "pearson")
cor.test(RES$NO3_did,RES$WSD_did, method = "pearson")
cor.test(RES$NO3_did,RES$SSD_did, method = "pearson")
cor.test(RES$NO3_did,RES$GDP_did, method = "pearson")

cor.test(RES$NH4_did,RES$Wintertemp_did, method = "pearson")
cor.test(RES$NH4_did,RES$Summertemp_did, method = "pearson")
cor.test(RES$NH4_did,RES$WSD_did, method = "pearson")
cor.test(RES$NH4_did,RES$SSD_did, method = "pearson")
cor.test(RES$NH4_did,RES$GDP_did, method = "pearson")

cor.test(RES$OC_did,RES$Wintertemp_did, method = "pearson")
cor.test(RES$OC_did,RES$Summertemp_did, method = "pearson")
cor.test(RES$OC_did,RES$WSD_did, method = "pearson")
cor.test(RES$OC_did,RES$SSD_did, method = "pearson")
cor.test(RES$OC_did,RES$GDP_did, method = "pearson")

cor.test(RES$EC_did,RES$Wintertemp_did, method = "pearson")
cor.test(RES$EC_did,RES$Summertemp_did, method = "pearson")
cor.test(RES$EC_did,RES$WSD_did, method = "pearson")
cor.test(RES$EC_did,RES$SSD_did, method = "pearson")
cor.test(RES$EC_did,RES$GDP_did, method = "pearson")

