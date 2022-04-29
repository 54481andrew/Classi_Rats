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

## Load local packages
library(raster)
library(plotmo)
library(dismo)

## Load packages 
library(ggplot2)
library(reshape)
library(plotrix)
library(parallel)
library(plyr)
library(geosphere)
library(gbm)
library(lubridate)
library(zoo)
library(MASS)

res <- 5 ## Choose resolution of predictors (1~ 01 degrees, 5~05 degrees
omit.preds <- c('Pcv','Ncv', 'TS.House', 'TS.InTown')

## - Specify paths to EFC and PREEMPT datasets
efc.data.path = '../Trap_Data_EFC_V2/Data/Clean/Aggregated_EFC_Capture_Data.csv'
preempt.data.path = '../Trap_Data_PRE_V2/Data/Clean/Aggregated_PREEMPT_Capture_Data.csv'

## - Specify folder that contains predictor stacks
raster.data.path <- "../../Storage/Raster_Data"
storage.fold <- '../../Storage'

## Data used to plot focal WA countries
foc.shp <- st_read(dsn = paste(storage.fold, '/Shapefiles/West_Africa', sep = ''),
                   layer = 'foc', quiet = TRUE)
foc.shp.ogr <- readOGR(dsn = paste(storage.fold, '/Shapefiles/West_Africa', sep = ''),
                       layer = 'foc', verbose = FALSE)

## Create rodent.data, the concatenated EFC / PRE dataset
source("Tools/Prep_Data.r")

## Add on MODIS predictors to dataset
## TODO: ADD IN CORRECT LC TYPE (YEAR) FROM MODIS 
source("Tools/Prep_Features.r")





rodent.data <- data.table::setDT(rodent.data)
agg.data = rodent.data[,.(Tot = sum(Tot), NumPosVirus = sum(NumPosVirus),
                          NumTestVirus = sum(NumTestVirus),
               TotTraps = sum(TotTraps)),
            by = list(Site, Longitude, Latitude,
                      Sp, Country, Source)]

rodent.data <- as.data.frame(rodent.data)
agg.data <- as.data.frame(agg.data)
agg.data1 <- merge(agg.data, unique(rodent.data[,c('Site', hab.pred.names)]), 
                   all.x = TRUE, all.y = FALSE)

temp.names <- names(agg.data1)
focal.cols <- c('Tot', 'NumTestVirus', 'NumPosVirus', 'Sp')
agg.data2 = reshape(agg.data1, timevar = 'Sp', idvar = temp.names[!(temp.names %in% focal.cols)],
            direction = 'wide')

## Fill in TS for each species and compute shannon diversity index
sp.list <- as.character(unique(agg.data1$Sp))
agg.data2[,'TotCaptures'] = rowSums(agg.data2[,paste0('Tot.', sp.list)], na.rm = TRUE)


## -- compute shannon diversity index
ts.names <- c()
p.names <- c()
shannon.diversity = 0
for(sp in sp.list){
    tot.name <- paste0('Tot.', sp)
    ts.name <- paste0('TS.', sp)
    ts.names <- c(ts.names, ts.name)
    vals <- agg.data2[,tot.name]
    agg.data2[is.na(vals), tot.name] <- 0
    agg.data2[, ts.name] <- agg.data2[,tot.name] / agg.data2$TotTraps

    p.name = paste0('p.', sp)
    p.names = c(p.names, p.name)
    agg.data2[,p.name] <- agg.data2[, tot.name] / agg.data2$TotCaptures
    p.sp = agg.data2[,paste0('p.', sp)]
    shan = p.sp*log(p.sp)
    shan[is.na(shan)] = 0
    shannon.diversity = shannon.diversity -  shan

}

agg.data2$Shannon <- shannon.diversity
agg.data2$NumPosVirus = rowSums(agg.data2[,paste0('NumPosVirus.', sp.list)], na.rm = TRUE)
agg.data2$NumTestVirus = rowSums(agg.data2[,paste0('NumTestVirus.', sp.list)], na.rm = TRUE)
agg.data2$LASV.Prop = with(agg.data2, NumPosVirus / NumTestVirus)
agg.data2$LASV.Prop.Mna = with(agg.data2, NumPosVirus.Mna / NumTestVirus.Mna)


## Add predictor indicating Rr presence 
agg.data2$RrPres <- 1*(agg.data2$TS.Rra > 0)
agg.data2$MnPres <- 1*(agg.data2$TS.Mna > 0)


## --- LASV presence in Mn vs Shannon diversity

## Graph of Shannon diversity vs LASV prevalence in rodents
p = ggplot(data = agg.data2) +
    geom_point(shape = 21, aes(x = Shannon, y = LASV.Prop.Mna, size = TotCaptures),
               fill = 'darkorange', alpha = 0.5) + 
    scale_size(range = c(1,10), breaks = c(0,50,100,200,300)) +
    xlab('Shannon diversity index') + ylab('LASV Prevalence in Mn') +
    ggtitle('PREEMPT and EFC Datasets')
p$labels$size = 'Total captures'
p
ggsave(filename = 'Figures/Both/Shannon_vs_LASV.png')

mod = glm(LASV.Prop.Mna~Shannon, family = quasibinomial, data = agg.data2,
    weights = TotCaptures)
summary(mod)
## Coefficients:
##             Estimate Std. Error t value Pr(>|t|)   
## (Intercept)  -2.8812     0.7815  -3.687  0.00169 **
## Shannon       0.6181     0.6078   1.017  0.32262   


mod = glm(LASV.Prop.Mna~Shannon, family = quasibinomial,
          data = subset(agg.data2, LASV.Prop.Mna > 0),
    weights = TotCaptures)
summary(mod)
## Coefficients:
##             Estimate Std. Error t value Pr(>|t|)    
## (Intercept)  -4.2666     0.5654  -7.546 0.000132 ***
## Shannon       1.9098     0.4442   4.300 0.003569 ** 



## --- Effect of Rr on shannon diversity

## Graph of Shannon diversity vs Rr Trap success
p = ggplot(data = agg.data2) +
    geom_point(shape = 21, aes(x = TS.Rra, y = Shannon, size = TotTraps),
               fill = 'darkorange', alpha = 0.5) + 
    scale_size(range = c(1,10), breaks = c(0, 500, 1000, 2000)) +
    xlab('TS: Rr') + ylab('Shannon diversity index') +
    ggtitle('PREEMPT and EFC Datasets')
p$labels$size = 'Total traps'
p
ggsave(filename = 'Figures/Both/Shannon_vs_Rr.png')

mod <- glm(Shannon~TS.Rra, data = agg.data2, family = gaussian)
summary(mod)
## Coefficients:
##             Estimate Std. Error t value Pr(>|t|)    
## (Intercept)  1.29324    0.07727  16.736 1.94e-15 ***
## TS.Rra      -9.21758    3.92911  -2.346   0.0269 *  

mod <- glm(Shannon~RrPres, data = agg.data2, family = gaussian)
summary(mod)

## 

## --- Effect of Mn on shannon diversity

## Graph of Shannon diversity vs Mn Trap success
p = ggplot(data = agg.data2) +
    geom_point(shape = 21, aes(x = TS.Mna, y = Shannon, size = TotTraps),
               fill = 'darkorange', alpha = 0.5) + 
    scale_size(range = c(1,10), breaks = c(0, 500, 1000, 2000)) +
    xlab('TS: Mn') + ylab('Shannon diversity index') +
    ggtitle('PREEMPT and EFC Datasets')
p$labels$size = 'Total traps'
p
ggsave(filename = 'Figures/Both/Shannon_vs_Mn.png')

mod <- glm(Shannon~TS.Mna, data = agg.data2, family = gaussian)
summary(mod)
## Coefficients:
##             Estimate Std. Error t value Pr(>|t|)    
## (Intercept)  1.22994    0.09233  13.322    4e-13 ***
## TS.Mna      -1.11357    3.60850  -0.309     0.76    


mod <- glm(Shannon~MnPres, data = agg.data2, family = gaussian)
summary(mod)

## Coefficients:
##             Estimate Std. Error t value Pr(>|t|)    
## (Intercept)   0.7317     0.1071   6.830    3e-07 ***
## MnPres        0.6429     0.1237   5.198    2e-05 ***

## --- Rr on LASV

## 
p = ggplot(data = agg.data2[agg.data2$NumTestVirus.Mna > 0,]) +
    geom_point(shape = 21, aes(x = TS.Rra, y = LASV.Prop.Mna, size = NumTestVirus.Mna),
               fill = 'darkorange', alpha = 0.5) + 
    scale_size(range = c(1,10), breaks = c(0,50,100,200,300)) +
    xlab('TS: Rr') + ylab('LASV Prevalence in Mn') +
    ggtitle('PREEMPT and EFC Datasets')
p$labels$size = 'Mn tested'
p
ggsave(filename = 'Figures/Both/Rr_vs_MnLASV.png')

mod = glm(LASV.Prop.Mna~ TS.Rra, weights = NumTestVirus.Mna, data = agg.data2,
          family = quasibinomial)
summary(mod)

mod = glm(LASV.Prop.Mna~ RrPres, weights = NumTestVirus.Mna, data = agg.data2,
          family = quasibinomial)
summary(mod)

## --- Mn vs Rr

## Graph of Shannon diversity vs LASV prevalence in rodents
p = ggplot(data = agg.data2) +
    geom_point(shape = 21, aes(x = TS.Rra, y = TS.Mna, size = TotTraps),
               fill = 'darkorange', alpha = 0.5) + 
    scale_size(range = c(1,10), breaks = c(0,500,1000,2000,3000)) +
    xlab('TS: Rr') + ylab('TS: Mn') +
    xlim(0, 0.0375) + 
    ggtitle('PREEMPT and EFC Datasets')
p$labels$size = 'Total traps'
p
ggsave(filename = 'Figures/Both/TS_Mn_vs_Rr_Truncated.png')


p = ggplot(data = agg.data2) +
    geom_point(shape = 21, aes(x = TS.Rra, y = TS.Mna, size = TotTraps),
               fill = 'darkorange', alpha = 0.5) + 
    scale_size(range = c(1,10), breaks = c(0,500,1000,2000,3000)) +
    xlab('TS: Rr') + ylab('TS: Mn') +
    ggtitle('PREEMPT and EFC Datasets')
p$labels$size = 'Total traps'
p
ggsave(filename = 'Figures/Both/TS_Mn_vs_Rr.png')

b <- 0.1
m <- 0.001
agg.data2$focal.rodents <- agg.data2$Tot.Mna + agg.data2$Tot.Rra
agg.data2$focal.rodents1 <- log10(agg.data2$focal.rodents + 1)
plot.agg.data2 <- subset(agg.data2, focal.rodents > 0)
plot.agg.data2$radius = (plot.agg.data2$focal.rodents1) / 10
plot.agg.data2$frac.Mna <- with(plot.agg.data2, Tot.Mna / focal.rodents)
plot.agg.data2$frac.Rra <- with(plot.agg.data2, Tot.Rra / focal.rodents)
g1 = ggmap(map) +
    xlab('Longitude') + ylab('Latitude')
g1
g2 = g1 +
    geom_scatterpie(aes(x=Longitude, y=Latitude, group=Site, r = radius),
                    data=plot.agg.data2, pie_scale = 0.1,
                    cols=c('frac.Mna', 'frac.Rra'),
                    color='black', size = 0.15, alpha = 0.5) + coord_equal() +
        scale_fill_discrete(name = 'Proportion',
                        labels = c('Mn',
                                   'Rr'))
g2 +  geom_scatterpie_legend(radius = log10(c(1,20,200) + 1)/10, x=-13.25, y=6.7,
                             labeller = function(x){10^(10*x) - 1}) +
    annotate("text", x=-13.125, y=7.2, label= "Mn + Rr Captures") +
    ggtitle('Mn : Rr in EFC and PREEMPT sites')
ggsave(filename = 'Figures/Both/Map_Mn_vs_Rr.png')



mod = glm(TS.Mna~TS.Rra, family = quasibinomial,
          data = agg.data2,
    weights = TotTraps)
summary(mod)
## Coefficients:
##             Estimate Std. Error t value Pr(>|t|)
## (Intercept)   -3.164      4.104  -0.771    0.448
## TS.Rra      -481.239   3826.619  -0.126    0.901

mod = glm(TS.Mna~RrPres, family = quasibinomial,
          data = agg.data2,
    weights = TotTraps)
summary(mod)
## Coefficients:
##             Estimate Std. Error t value Pr(>|t|)    
## (Intercept)  -3.0554     0.1200 -25.456  < 2e-16 ***
## RrPres       -2.4773     0.4365  -5.676 5.71e-06 ***
## ---
## Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

## (Dispersion parameter for quasibinomial family taken to be 11.92822)

mod = glm(TS.Mna~TS.Rra + TS.Mer, family = quasibinomial,
          data = agg.data2,
    weights = TotTraps)
summary(mod)


## --- Try another approach using SLN reshuffling method

## This will fit the function k1*exp(k2*x) to the data, and assess whether the
## resulting likelihood is better than the background null models

model.fun <- function(x, parms){parms[1]*exp(x*parms[2])}
ll.fun <- function(parms){
    ##parms <- c(0.05, -50)
    preds <- model.fun(agg.data2$TS.Rra, parms) + 1e-5
    ll.tot <- -sum(dbinom(x = agg.data2$Tot.Mna, size = agg.data2$TotTraps, prob = preds,
                         log = TRUE))
    return(ll.tot)
}

out = optim(par = c(0.005, -5000), f = ll.fun)

plot(TS.Mna~TS.Rra, agg.data2)
lines(xseq, model.fun(xseq, out$par))
ll.star <- out$value

fit.dat <- data.frame(x = xseq, y = model.fun(xseq, out$par))

p = ggplot(data = agg.data2) +
    geom_point(shape = 21, aes(x = TS.Rra, y = TS.Mna, size = TotTraps),
               fill = 'darkorange', alpha = 0.5) + 
    scale_size(range = c(1,10), breaks = c(0,500,1000,2000,3000)) +
    xlab('TS: Rr') + ylab('TS: Mn') +
    ggtitle('PREEMPT and EFC Datasets') +
    geom_line(aes(x = x, y = y), data = fit.dat)
p$labels$size = 'Total traps'
p
ggsave(filename = 'Figures/Both/TS_Mn_vs_Rr_Fit.png')


## Get background likelihood by randomly pairing TS.Mna and TS.Rra
store <- c()
for(ii in 1:10000){
    res <- sample(1:nrow(agg.data2))
    new.data <- cbind(agg.data2[res,c('Tot.Mna', 'TS.Mna','TotTraps')], TS.Rra = agg.data2[,'Tot.Rra'])
    model.fun <- function(x, parms){parms[1]*exp(x*parms[2])}
    ll.fun <- function(parms){
        ##parms <- c(0.05, -50)
        preds <- model.fun(new.data$TS.Rra, parms) + 1e-5
        ll.tot <- -sum(dbinom(x = new.data$Tot.Mna, size = new.data$TotTraps, prob = preds,
                              log = TRUE))
        return(ll.tot)
    }
    out.res = optim(par = c(0.005, -5000), f = ll.fun)
    store <- rbind(store, data.frame(ll = out.res$value))
}
names(store) = 'Null_LL'

ggplot(data = store) + geom_histogram(aes(x = Null_LL))


a = quantile(-store$Null_LL, probs = seq(0,1,by = 0.01))
f <- approxfun(x = unlist(a, use.names = FALSE), y = seq(0,1,by = 0.01))


plot(TS.Mna~TS.Rra, new.data)


## --- 

## Map of Shannon diversity 
lat <- 9.84 
lon <- -11.03333 
height = 4
width = 4
xmin = lon - width
xmax = lon + width
ymin = lat - height
ymax = lat + height

sl_borders <- c(bottom = ymin,
                 top = ymax,
                 left = xmin,
                 right = xmax)
map <- get_stamenmap(sl_borders, zoom = 8, maptype = "terrain")

g1 = ggmap(map)
g1

g1 = ggmap(map) +
    xlab('Longitude') + ylab('Latitude') +
    geom_point(shape = 21, mapping = aes(x = Longitude, y = Latitude, size = Shannon),
               data = agg.data2,
               fill = 'darkorange', alpha = 0.5) + 
    scale_size(range = c(2,10),breaks=c(0,0.5,1,1.5,2)) +
    ggtitle('Shannon diversity in EFC and PREEMPT sites')
g1
ggsave(filename = 'Figures/Both/Map_shannon_diversity.png')



require(factoextra)

rownames(agg.data2) = paste(agg.data2$Site)
pr.out = princomp(x = agg.data2[,p.names], cor = TRUE)

png(filename = 'Figures/Both/Community_scree.png')
plot(cumsum(pr.out$sdev^2) / sum(pr.out$sdev^2))
dev.off()



fviz_pca_biplot(pr.out, repel = TRUE, habillage = agg.data2$Country)
ggsave('Figures/Both/Community_biplot12.png', width = 14, height = 7)

fviz_pca_biplot(pr.out, axes = c(3,4), repel = TRUE, habillage = agg.data2$Country)
ggsave('Figures/Both/Community_biplot34.png', width = 7, height = 7)

fviz_pca_biplot(pr.out, axes = c(5,6), repel = TRUE, habillage = agg.data2$Country)
ggsave('Figures/Both/Community_biplot56.png', width = 7, height = 7)


f.out <- fviz_nbclust(x = agg.data2[,p.names], FUNcluster = kmeans)
k.out = kmeans(x = agg.data2[,p.names], centers = 3, nstart = 25)


fviz_cluster(k.out, data = agg.data2[,p.names])
## 1: Mn
## 2: Rr
## 3. Me


## --- Create map of community
k.out = kmeans(x = agg.data2[,p.names], centers = 3, nstart = 25)
agg.out = cbind(agg.data2, cluster =  k.out$cluster)
g1 = ggmap(map) +
    xlab('Longitude') + ylab('Latitude')
g1
g2 = g1 + geom_point(data = agg.out, aes(x = Longitude, y = Latitude,
                                         fill = as.factor(cluster)),
                     size = 5, pch = 21, color = 'black')
g2


## Total rodent abundance
ggplot(data = agg.data2, aes(x = TotCaptures, y = reorder(Site, TotCaptures), fill = Country)) +
    geom_bar(stat = 'identity') + ylab('') + xlab('Rodents captured')
ggsave(filename = 'Figures/Both/Rodent_Captures.png')

ggplot(data = agg.data2, aes(x = TotCaptures / TotTraps, y = reorder(Site, TotCaptures / TotTraps),
                             fill = Country)) +
    geom_bar(stat = 'identity') + ylab('') + xlab('Rodents captured per trap')
ggsave(filename = 'Figures/Both/Rodent_Captures_Per_Trap.png')

## Rodent species
a = colSums(agg.data2[,ts.names]*agg.data2$TotTraps)
sp.out = data.frame(Species = sp.list, Captured = unlist(a, use.names = FALSE))
ggplot(data = sp.out, aes(y = reorder(Species, Captured), x = Captured)) +
    geom_bar(stat = 'identity')
ggsave(filename = 'Figures/Both/Rodent_Total_Captures.png')


## --- Incorporate LC predictors


pc.vars <- c('Elev', 'Tmu', 'Pmu', 'Pc', 'Pm', 'Pop', 'Evergreen_Broadleaf_Forest',
             'Woody_Savanna', 'Savannas', 'Grasslands')
pc.out <- prcomp(x = agg.data2[,c(pc.vars)], scale.=TRUE)
pc.sd <- pc.out$sdev
cu.pc.var <- cumsum(pc.sd^2) / sum(pc.sd^2)
plot(cu.pc.var)

agg.data3 = cbind(agg.data2, pc.out$x)
out = glm(Shannon ~  PC1 + PC2 + PC3 + PC4 + PC5 + PC6, family = gaussian,
    data = agg.data3)
summary(out)

step.out = stepAIC(out)
summary(step.out)

plot(Shannon ~ PC4, agg.data3)

## --- 

pc.vars = hab.pred.names
pc.out <- prcomp(x = agg.data2[,c(pc.vars)], scale.=TRUE)
pc.sd <- pc.out$sdev
cu.pc.var <- cumsum(pc.sd^2) / sum(pc.sd^2)
plot(cu.pc.var)

agg.data3 = cbind(agg.data2, pc.out$x)
pc.names = paste('Shannon ~ ', paste(paste0('PC', 1:10), collapse = ' + '))
out = glm(as.formula(pc.names), family = gaussian,
    data = agg.data3)
summary(out)

step.out = stepAIC(out)
summary(step.out)

plot(Shannon ~ PC8, agg.data3)
lines(step.out$fitted)

fviz_pca_biplot(pc.out, axes = c(1,2), repel = TRUE, habillage = agg.data2$Country)

