## Create a simple map of Sierra Leone and Guinea with pie charts showing rodent
## capture composition.
##_2 adds in new dataset from EFC
##_3 only has Mn LASV prevalence with pie size scaled by trap rate

##require(ggthemes)
##require(ggplot2)
library(scatterpie)
library(ggmap)
require(rgdal)
require(sf)

storage.fold <- '../../Storage'
## Data used to plot focal WA countries
foc.shp <- st_read(dsn = paste(storage.fold, '/Shapefiles/West_Africa', sep = ''),
                   layer = 'foc', quiet = TRUE)
foc.shp.ogr <- readOGR(dsn = paste(storage.fold, '/Shapefiles/West_Africa', sep = ''),
                       layer = 'foc', verbose = FALSE)

## PREEMPT dataset
pre.dat <- read.csv(file = '../Trap_Data_PRE_V2_1/Data/Clean/PRE_Trap_Data.csv') 
##pre.dat <- subset(pre.dat, Sp=='Mn')
pre.dat$Source <- 'PRE'
pre.dat <- data.table::setDT(pre.dat)
agg.pre.dat <- pre.dat[,.(Longitude = mean(Longitude), Latitude = mean(Latitude),
                          Tot.Rr = sum(Sp=='Rr', na.rm = TRUE),
                          Tot.Mm = sum(Sp=='Mm', na.rm = TRUE),
                          Tot.Ct = sum(Sp=='Ct', na.rm = TRUE),
                          Tot.Me = sum(Sp=='Me', na.rm = TRUE),
                          Tot.Pd = sum(Sp=='Pd', na.rm = TRUE),
                          Tot.Co = sum(Sp=='Co', na.rm = TRUE),
                          Tot.Mn = sum(Sp=='Mn', na.rm = TRUE),
                          Tot.Ur = sum(Sp=='Ur', na.rm = TRUE),
                          Tot.Ls = sum(Sp=='Ls', na.rm = TRUE),
                          Tot.Ms = sum(Sp=='Ms', na.rm = TRUE),
                          Tot.Hs = sum(Sp=='Hs', na.rm = TRUE),
                          Tot.Pr = sum(Sp=='Pr', na.rm = TRUE),
                          Tot.Gs = sum(Sp=='Gs', na.rm = TRUE),
                          Tot.Cb = sum(Sp=='Cb', na.rm = TRUE),
                          Tot.Cg = sum(Sp=='Cg', na.rm = TRUE),
                          Tot.Gk = sum(Sp=='Gk', na.rm = TRUE),
                          Tot.Lf = sum(Sp=='Lf', na.rm = TRUE),
                          Uni.Sp = length(unique(Sp[!is.na(Sp)])),
                          TotTraps = .N,
                          LASV.Pos = sum(Lassa==1, na.rm = TRUE),
                          LASV.Test = sum(Rodent==1)),
                       by = list(Village)]
agg.pre.dat

## Fill in prevalence and trap rate
agg.pre.dat$LASV.Prop <- with(agg.pre.dat, LASV.Pos / LASV.Test)

sp.list <- as.character(unique(pre.dat$Sp[!is.na(pre.dat$Sp)]))
agg.pre.dat <- as.data.frame(agg.pre.dat)
agg.pre.dat[,'TotCaptures'] = rowSums(agg.pre.dat[,paste0('Tot.', sp.list)])
shannon.diversity = 0
for(sp in sp.list){
    agg.pre.dat[,paste0('TS.', sp)] <- agg.pre.dat[, paste0('Tot.', sp)] / agg.pre.dat[,'TotTraps'] 
    agg.pre.dat[,paste0('p.', sp)] <- agg.pre.dat[, paste0('Tot.', sp)] / agg.pre.dat$TotCaptures
    p.sp = agg.pre.dat[,paste0('p.', sp)]
    shan = p.sp*log(p.sp)
    shan[is.na(shan)] = 0
    shannon.diversity = shannon.diversity -  shan
}
agg.pre.dat$Shannon <- shannon.diversity

FIX THIS TO JUST BE MN LASV PREVALENCE
## Graph of Shannon diversity vs LASV prevalence in rodents
p = ggplot(data = agg.pre.dat) +
    geom_point(shape = 21, aes(x = Shannon, y = LASV.Prop, size = TotCaptures),
               fill = 'darkorange', alpha = 0.5) + 
    scale_size(range = c(1,10), breaks = c(0,50,100,200,300)) +
    xlab('Shannon diversity index') + ylab('LASV Prevalence') +
    ggtitle('PREEMPT Dataset')
p$labels$size = 'Total captures'
p
ggsave(filename = 'Figures/PRE_Shannon_vs_LASV.png')


agg.pre.dat$TS.MR <- with(agg.pre.dat, Tot.Mn / (Tot.Mn + Tot.Rr))
agg.pre.dat$TS.nMR <- with(agg.pre.dat, Tot.Rr / (Tot.Mn + Tot.Rr))





lat <- 9.84 - 1.5
lon <- -11.03333 - 0.75
height = 2
width = 2
xmin = lon - width
xmax = lon + width
ymin = lat - height
ymax = lat + height

sl_borders <- c(bottom = ymin,
                 top = ymax,
                 left = xmin,
                 right = xmax)
map <- get_stamenmap(sl_borders, zoom = 8, maptype = "terrain")

## Map of diversity 
g1 = ggmap(map) +
    xlab('Longitude') + ylab('Latitude') +
    geom_point(shape = 21, mapping = aes(x = Longitude, y = Latitude, size = Uni.Sp),
               data = agg.pre.dat,
               fill = 'darkorange', alpha = 0.5) + 
    scale_size(range = c(1,10),breaks=c(0, 4,6,8,10))
g1

## Map of Shannon diversity 
g1 = ggmap(map) +
    xlab('Longitude') + ylab('Latitude') +
    geom_point(shape = 21, mapping = aes(x = Longitude, y = Latitude, size = Shannon),
               data = agg.pre.dat,
               fill = 'darkorange', alpha = 0.5) + 
    scale_size(range = c(2,10),breaks=c(0,0.5,1,1.5,2))
g1



b <- 0.1
m <- 0.001
agg.pre.dat$focal.rodents <- agg.pre.dat$Tot.Mn + agg.pre.dat$Tot.Rr
agg.pre.dat$focal.rodents1 <- log10(agg.pre.dat$focal.rodents + 1)
plot.agg.pre.dat <- subset(agg.pre.dat, focal.rodents > 0)
plot.agg.pre.dat$radius = (plot.agg.pre.dat$focal.rodents1) / 10
plot.agg.pre.dat$frac.Mn <- with(plot.agg.pre.dat, Tot.Mn / focal.rodents)
plot.agg.pre.dat$frac.Rr <- with(plot.agg.pre.dat, Tot.Rr / focal.rodents)
g1 = ggmap(map) +
    xlab('Longitude') + ylab('Latitude')
g1
g2 = g1 +
    geom_scatterpie(aes(x=Longitude, y=Latitude, group=Village, r = radius),
                    data=plot.agg.pre.dat, pie_scale = 0.1,
                    cols=c('frac.Mn', 'frac.Rr'),
                    color='black', size = 0.15, alpha = 0.5) + coord_equal() +
        scale_fill_discrete(name = 'Proportion',
                        labels = c('Mn',
                                   'Rr'))
g2 +  geom_scatterpie_legend(radius = log10(c(1,20,200) + 1)/10, x=-13.25, y=6.7,
                             labeller = function(x){10^(10*x) - 1}) +
    annotate("text", x=-13.125, y=7.2, label= "Mn + Rr Captures")
ggsave(filename = 'Figures/Map_Mn_vs_Rr.png')


### --- 



mod.1 <- glm(TS.Mn~Latitude, plot.agg.pre.dat, family = quasibinomial,
             weights = TotTraps)
summary(mod.1)
## Mn increases with latitude**

mod.1 <- glm(TS.Rr~Latitude, plot.agg.pre.dat, family = quasibinomial,
             weights = TotTraps)
summary(mod.1)
## Rr decreases with latitude 

mod.2 <- glm(TS.Mn~TS.Rr, plot.agg.pre.dat, family = quasibinomial,
             weights = TotTraps)
summary(mod.2)

