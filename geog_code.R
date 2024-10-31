#Load libraries
library(ggmap)
library(mapdata)
library(viridis)
library(ggplot2)
library(grid)
library(rgdal)
library(ggalt)
library(ggthemes)
library(sp)
library(rgeos)
library(raster)
library(devtools)
library(data.table)
library(RColorBrewer)
library(rgdal)
library(rgeos)

devtools::install_github("eliocamp/ggalt@new-coord-proj",force=TRUE)

#####################

# import conservation evidence freshwater data (fwmarinemammal_data.csv)
cedat <- fread(choose.files())

########################
#rename dataframe and only extract columns that we need
mapdata = unique(cedat[,list(pageid,lat,long,country,syn)])
mapdata$lat = as.numeric(mapdata$lat)
mapdata$long = as.numeric(mapdata$long)
mapdata$country = as.character(mapdata$country)

#####################
#check there are no errors in lat and longitude
#these are reviews
unique(mapdata[(lat==0 | is.na(lat)==TRUE),pageid])
unique(mapdata[(long==0 | is.na(long)==TRUE),pageid])

###
#this code below is a bit complex but what it essentially does it to take the conservation evidence data
#determines how many studies fit within each grid cell (4x4 degrees in this case) and creates a dataframe
#containing this data. It first works out which grid cells each study is in, then aggregates to find the number of studies in each grid cell
data.df=mapdata
data.df$LONGITUDE = data.df$long
data.df$LATITUDE = data.df$lat
#can set cellsize to be bigger or small cells, at the moment it is 4x4 degrees
ji <- function(xy, origin=c(0,0), cellsize=c(4,4)) {
  t(apply(xy, 1, function(z) cellsize/2+origin+cellsize*(floor((z - origin)/cellsize))))
}
JI <- ji(cbind(data.df$LONGITUDE, data.df$LATITUDE))
data.df$X <- JI[, 1]
data.df$Y <- JI[, 2]
data.df$Cell <- paste(data.df$X, data.df$Y)
counts <- by(data.df, data.df$Cell, function(d) c(d$X[1], d$Y[1], nrow(d)))
counts.m <- matrix(unlist(counts), nrow=3)
rownames(counts.m) <- c("X", "Y", "Count")
count.max <- max(counts.m["Count",])
counts.m2=data.frame(cbind(t(counts.m)),stringsAsFactors = FALSE)
counts.m2$X = as.numeric(counts.m2$X)
counts.m2$Y = as.numeric(counts.m2$Y)
counts.m2$Count = as.numeric(counts.m2$Count)
counts.m3 <- counts.m2[which(is.na(counts.m2$X)==FALSE&is.na(counts.m2$Y)==FALSE),]

#we now aggregate the data to find number of studies in each unique grid cell 
counts.m4 = aggregate(Count ~ X + Y, data=counts.m3,sum)
max(counts.m4$Count) #work out the highest number in each grid cell so we know what the 
                     #range should be in the number of studies on our maps below.


#linear model
#counts studies multiple times if in different ecoregions
studiespereco <- data.table(Marine_Ecoregion=unique(cedat$Marine_Ecoregion),studies=sapply(unique(cedat$Marine_Ecoregion), function(x){length(unique(cedat[Marine_Ecoregion==x,pageid]))}))
sum(studiespereco$studies)
length(unique(cedat$pageid))

#import dataset on number of species in each ecoregion obtained from overlaying MEOW dataset and IUCN range maps and merge with studies data
map.df_meow <- fread(choose.files())
ecostudiesmams <- merge(studiespereco,map.df_meow,by.x="Marine_Ecoregion",by.y="ECOREGION",all.y=TRUE)
ecostudiesmams <- unique(ecostudiesmams[,list(Marine_Ecoregion,studies,num.species)])
sum(ecostudiesmams$studies,na.rm=TRUE)
ecostudiesmams[is.na(num.species),num.species:=0]
ecostudiesmams[is.na(studies),studies:=0]

### Check assumptions of the model
#install.packages("DHARMa")
library(DHARMa)

### geographic model if poisson - results show geog model cannot use poisson distribution
Geogmodel1 <- glm(studies~num.species, family="poisson",data=ecostudiesmams)
testDispersion(Geogmodel1)
simulationOutput <- simulateResiduals(fittedModel = Geogmodel1, plot = F)
residuals(simulationOutput)
plot(simulationOutput)
testZeroInflation(simulationOutput)

### ZERO-INFLATED POISSON REGRESSION suggests that the excess zeros are generated from a separate process from the count values - not sure if this is best option for this data

### Tests for the quasipoisson distribution for geog model - symmetrical
Geogmodelqp <- glm(studies~num.species, family="quasipoisson",data=ecostudiesmams)
dev_residuals <- residuals(Geogmodelqp, type = "deviance")
# Plot deviance residuals against predicted values
plot(fitted(Geogmodelqp), dev_residuals, main = "Deviance Residuals vs. Fitted Values")
summary(glm(studies~num.species, family="quasipoisson",data=ecostudiesmams))
