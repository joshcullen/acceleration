# Download and extract environmental covariates for extraction

library(elevatr)
library(tidyverse)
library(sf)
library(raster)
library(lubridate)
library(sp)
library(rgee)
library(raster)
library(rasterVis)
library(googledrive)




###################
### Import data ###
###################

dat<- read.csv("Giant Armadillo state estimates.csv", as.is = T)
dat<-  dat %>% 
  rename(x = easting, y = northing) %>% 
  mutate(across(c('z.map','z.post.thresh','z.post.max'), factor,
                levels = c("Slow-Turn","Slow-Unif","Exploratory","Transit","Unclassified"))
  )
dat$month<- month.abb[month(dat$date)]
dat$month<- factor(dat$month, levels = month.abb[c(5:12,1)])
dat$season<- ifelse(dat$month %in% c(month.abb[3:5]), "Fall",
                    ifelse(dat$month %in% c(month.abb[6:8]), "Winter",
                           ifelse(dat$month %in% c(month.abb[9:11]), "Spring", "Summer")))
dat$season<- factor(dat$season, levels = c("Fall","Winter","Spring","Summer"))





#############################
### Topographic variables ###
#############################

# Create dummy raster for extent

rast<- raster(ext=extent(c(min(dat$x) - 3000, max(dat$x) + 3000, min(dat$y) - 3000,
                             max(dat$y + 3000))),
                crs = "+init=epsg:32721",
                res = 30)
values(rast)<- 0

#download elevation data
dem<- get_elev_raster(rast, z = 12) #18 m res
dem2<- crop(dem, rast)
plot(dem2)
points(dat$x, dat$y)


# Visualize different features using elev (i.e., slope, aspect, hillshade)
slope = terrain(dem2, opt='slope')
aspect = terrain(dem2, opt='aspect')
hill = hillShade(slope, aspect, 40, 270)

#plot elev
plot(dem2, main="DEM", col=terrain.colors(25,alpha=0.7))
#plot slope
plot(slope, main="Slope", col=topo.colors(25,alpha=0.7))
#plot aspect
plot(aspect, main="Aspect", col=rainbow(25,alpha=0.7))
#plot DEM w/ hillshade
plot(hill,
     col=grey(1:100/100),  # create a color ramp of grey colors for hillshade
     legend=FALSE,         # no legend, we don't care about the grey of the hillshade
     main="DEM",
     axes=FALSE)           # makes for a cleaner plot, if the coordinates aren't necessary
plot(dem2, 
     axes=FALSE,
     col=terrain.colors(12, alpha=0.35), add=TRUE)
points(dat$x, dat$y, cex = 0.25)


## Export DEM layer
writeRaster(dem2, filename='giantarm_dem.tif', format="GTiff", overwrite=TRUE)




####################################################
### Tasseled Cap Indices via Google Earth Engine ###
####################################################

### Part of code taken from Oeser et al 2020
# https://code.earthengine.google.com/?accept_repo=users/julianoeser/Oeser_et_al_2019_Habitat_metrics

ee_Initialize(email = "joshcullen10@gmail.com", drive = TRUE)

# Date range for giant armadillos
range(dat$date)  #2019-05-19 to 2020-01-24


# Define bounds for Region of Interest (ROI) based on data extent
bounds<- st_as_sf(data.frame(rasterToPoints(dem2)), coords = c("x","y"), 
                  crs = "+init=epsg:32721") %>%  #convert to sf object
  st_bbox() %>%  #extract bounding box
  st_as_sfc() %>%  #convert into polygon
  sf_as_ee()  #convert into GEE format 



## Define functions to apply on imagery

# This functions masks clouds and cloud shadows (bits 3 and 5)
maskL8sr<- function(image) {
  cloudShadowBitMask<- bitwShiftL(1, 3)
  cloudsBitMask<- bitwShiftL(1, 5)
  qa<- image$select('pixel_qa')  # Get the pixel QA band
  # Both flags should be set to zero, indicating clear conditions
  mask<- qa$bitwiseAnd(cloudShadowBitMask)$eq(0)$
    And(qa$bitwiseAnd(cloudsBitMask)$eq(0))
  return(image$updateMask(mask))
}


# This function calculates Tasseled Cap indices
add_TC = function(image) {
  
  img = image$select(c('B', 'G', 'R', 'NIR', 'SWIR1', 'SWIR2'))
  
  # coefficients for Landsat surface reflectance (Crist 1985)
  # A TM tasseled cap equivalent transformation for reflectance factor data.
  # Remote Sensing of Environment, 17(3), 301-306.
  
  brightness_c= ee$Image(list(c(0.2043, 0.4158, 0.5524, 0.5741, 0.3124, 0.2303)))
  greenness_c= ee$Image(list(c(-0.1603, 0.2819, -0.4934, 0.7940, -0.0002, -0.1446)))
  wetness_c= ee$Image(list(c(0.0315,  0.2021,  0.3102,  0.1594, -0.6806, -0.6109)))
  
  brightness = img$multiply(brightness_c)
  greenness = img$multiply(greenness_c)
  wetness = img$multiply(wetness_c)
  
  brightness = brightness$reduce(ee$call("Reducer.sum"))$toFloat()  #needs to be float, not double
  greenness = greenness$reduce(ee$call("Reducer.sum"))$toFloat()
  wetness = wetness$reduce(ee$call("Reducer.sum"))$toFloat()
  
  return(image$addBands(ee$Image(brightness)$
                        addBands(greenness)$
                        addBands(wetness)$
                        rename(c('brightness','greenness','wetness'))))
}



## Access and map functions to Landsat 8 imagery (images collected every 16 days)

l8bands = c('B2', 'B3', 'B4', 'B5', 'B6', 'B7', 'B1', 'B10', 'B11', 'pixel_qa')
l8band_names = c('B', 'G', 'R', 'NIR', 'SWIR1', 'SWIR2', 'UB', 'T', 'T2','pixel_qa')

l8<- ee$ImageCollection('LANDSAT/LC08/C01/T1_SR')$
  select(l8bands, l8band_names)$
  filterBounds(bounds)$
  filterDate('2019-05-01', '2020-01-31')$
  sort("system:time_start", TRUE)$  # Sort the collection in chronological order
  map(function(x) x$reproject("EPSG:32721"))

print(l8$size()$getInfo())  #check number of images
ee_print(l8)
print(l8$first()$getInfo())

## Calculate TC indices and assign to separate objects
tc = l8$map(maskL8sr)$
  map(add_TC)$
  select(c('brightness', 'greenness', 'wetness'))
ee_print(tc)
print(tc$first()$getInfo())


brightness = tc$select('brightness')
greenness = tc$select('greenness')
wetness = tc$select('wetness')




### Export data locally

setwd("~/Documents/Snail Kite Project/Data/R Scripts/acceleration/TC_Brightness")
tc_brightness<- ee_imagecollection_to_local(
  ic = brightness,
  scale = 30,
  region = bounds,
  via = 'drive'
)

setwd("~/Documents/Snail Kite Project/Data/R Scripts/acceleration/TC_Greenness")
tc_greenness<- ee_imagecollection_to_local(
  ic = greenness,
  scale = 30,
  region = bounds,
  via = 'drive'
)

setwd("~/Documents/Snail Kite Project/Data/R Scripts/acceleration/TC_Wetness")
tc_wetness<- ee_imagecollection_to_local(
  ic = wetness,
  scale = 30,
  region = bounds,
  via = 'drive')
)






#############################################################
### Mosaic and Average Dynamic Covariates by Month/Season ###
#############################################################


#### Brightness ####

setwd("~/Documents/Snail Kite Project/Data/R Scripts/acceleration/TC_Brightness")

#load files as raster brick
bright<- list.files(getwd(), pattern = "*.tif$")
bright.stack<- stack(bright)
bright.stack2<- flip(bright.stack, direction = 'y')
names(bright.stack2)<- names(bright.stack)
bright.brick<- bright.stack2

#change values x == 0 to NA (these are masked pixels)
raster::values(bright.brick)[raster::values(bright.brick) == 0] <- NA

#rename images by date
names(bright.brick)<- gsub(pattern = "LC08_22.07._", replacement = "", names(bright.brick))


#visualize original rasters
breaks<- seq(0, 15000, by=10)
cols<- colorRampPalette(RColorBrewer::brewer.pal(9, 'Reds'))(length(breaks)-1)
rasterVis::levelplot(bright.brick[[1:6]], at=breaks, col.regions=cols, main="Brightness")



### Mosaic RasterBrick over space and time

mosaic.bright<- list()
tmp1<- gsub(pattern = "LC08_22.07._", replacement = "", names(bright.stack2))
ind.m<- format(as.Date(tmp1, format = "%Y%m%d"), format = "%m") %>% 
  as.numeric()
ord.months<- c(5:12,1)


for (i in 1:length(unique(ind.m))) {
  ind<- which(ind.m %in% ord.months[i])  #using unique(ind.m) instead of ord.months results in wrong order
  bright.sub<- bright.stack2[[ind]]
  
  # time1<- gsub(pattern = "LC08_22.07._", replacement = "", names(bright.stack))[i]
  ind.pr<- stringr::str_sub(names(bright.stack2)[ind], 6, 11)
  
  if (length(unique(ind.pr)) > 1) {
    tile.list<- vector("list", length = length(unique(ind.pr)))
    
    for (j in 1:length(unique(ind.pr))) {
      ind.tiles<- grep(pattern = unique(ind.pr)[j], names(bright.stack2)[ind])
      tile.list[[j]]<- mean(bright.sub[[ind.tiles]], na.rm=TRUE)
    }
    
    tile.list$fun <- mean
    tile.list$na.rm <- TRUE
    
    tile.mosaic<- do.call(raster::mosaic, tile.list)
    
  } else {
    tile.mosaic<- mean(bright.sub, na.rm=TRUE)
  }
  
  mosaic.bright[[i]]<- tile.mosaic
}


coordinates(dat) <- ~x+y


mosaic.bright<- raster::brick(mosaic.bright)
names(mosaic.bright)<- month.abb[c(5:12,1)]
rasterVis::levelplot(mosaic.bright, at=breaks, col.regions=cols, main="Bright") +
  layer(sp.points(dat, cex=0.1, pch = 16, col = alpha("black", 0.5)))



### If wanting to aggregate into seasons:

# May is Fall
# Jun/Jul/Aug are Winter
# Sep/Oct/Nov are Spring
# Dec/Jan are Summer

season.ind<- c("Fall", rep("Winter",3), rep("Spring",3), rep("Summer",2))
bright.season<- stackApply(mosaic.bright, season.ind, fun = mean)
names(bright.season)<- c("Fall","Winter","Spring","Summer")
rasterVis::levelplot(bright.season, at=breaks, col.regions=cols, main="Bright") +
  layer(sp.points(dat, cex=0.1, pch = 16, col = alpha("black", 0.5)))




### Try plotting in ggplot2 so points can be plotted by month and season

mosaic.bright.df<- as.data.frame(mosaic.bright, xy=TRUE)
mosaic.bright.df<- pivot_longer(mosaic.bright.df, cols = -c(x,y), names_to = "month",
                              values_to = "bright")
mosaic.bright.df$month<- factor(mosaic.bright.df$month, levels = month.abb[c(5:12,1)])
mosaic.bright.df<- filter(mosaic.bright.df, bright > 0)  #only keep values from 0-1

dat<- as.data.frame(dat)


ggplot() +
  geom_raster(data = mosaic.bright.df, aes(x, y, fill = bright)) +
  geom_point(data = dat, aes(x, y, color = id), size = 0.5, alpha = 0.75) + 
  theme_bw() +
  scale_fill_gradientn(colors = cols, na.value = "transparent", limits = c(0,10000)) +
  scale_color_brewer(palette = "Dark2") +
  coord_fixed() +
  facet_wrap(~ month) +
  theme(strip.text = element_text(size = 12, face = "bold"))




bright.season.df<- as.data.frame(bright.season, xy=TRUE)
bright.season.df<- pivot_longer(bright.season.df, cols = -c(x,y), names_to = "season",
                              values_to = "bright")
bright.season.df$season<- factor(bright.season.df$season,
                               levels = c("Fall","Winter","Spring","Summer"))



ggplot() +
  geom_raster(data = bright.season.df, aes(x, y, fill = bright)) +
  geom_point(data = dat, aes(x, y, color = id), size = 0.5, alpha = 0.75) + 
  theme_bw() +
  scale_fill_gradientn(colors = cols, na.value = "transparent", limits = c(0,10000)) +
  scale_color_brewer(palette = "Dark2") +
  coord_fixed() +
  facet_wrap(~ season) +
  theme(strip.text = element_text(size = 12, face = "bold"))






#### Greenness ####

setwd("~/Documents/Snail Kite Project/Data/R Scripts/acceleration/TC_Greenness")

#load files as raster brick
green<- list.files(getwd(), pattern = "*.tif$")
green.stack<- stack(green)
green.stack2<- flip(green.stack, direction = 'y')
names(green.stack2)<- names(green.stack)
green.brick<- green.stack2

#change values x == 0 to NA (these are masked pixels)
raster::values(green.brick)[raster::values(green.brick) == 0] <- NA

#rename images by date
names(green.brick)<- gsub(pattern = "LC08_22.07._", replacement = "", names(green.brick))


#visualize original rasters
breaks<- seq(-2000, 5000, by=1)
cols<- colorRampPalette(RColorBrewer::brewer.pal(9, 'Greens'))(length(breaks)-1)
rasterVis::levelplot(green.brick[[1:6]], at=breaks, col.regions=cols, main="Greenness")



### Mosaic RasterBrick over space and time

mosaic.green<- list()
tmp1<- gsub(pattern = "LC08_22.07._", replacement = "", names(green.stack2))
ind.m<- format(as.Date(tmp1, format = "%Y%m%d"), format = "%m") %>% 
  as.numeric()
ord.months<- c(5:12,1)


for (i in 1:length(unique(ind.m))) {
  ind<- which(ind.m %in% ord.months[i])  #using unique(ind.m) instead of ord.months results in wrong order
  green.sub<- green.stack2[[ind]]
  
  # time1<- gsub(pattern = "LC08_22.07._", replacement = "", names(green.stack))[i]
  ind.pr<- stringr::str_sub(names(green.stack2)[ind], 6, 11)
  
  if (length(unique(ind.pr)) > 1) {
    tile.list<- vector("list", length = length(unique(ind.pr)))
    
    for (j in 1:length(unique(ind.pr))) {
      ind.tiles<- grep(pattern = unique(ind.pr)[j], names(green.stack2)[ind])
      tile.list[[j]]<- mean(green.sub[[ind.tiles]], na.rm=TRUE)
    }
    
    tile.list$fun <- mean
    tile.list$na.rm <- TRUE
    
    tile.mosaic<- do.call(raster::mosaic, tile.list)
    
  } else {
    tile.mosaic<- mean(green.sub, na.rm=TRUE)
  }
  
  mosaic.green[[i]]<- tile.mosaic
}


coordinates(dat) <- ~x+y


mosaic.green<- raster::brick(mosaic.green)
names(mosaic.green)<- month.abb[c(5:12,1)]
rasterVis::levelplot(mosaic.green, at=breaks, col.regions=cols, main="Green") +
  layer(sp.points(dat, cex=0.1, pch = 16, col = alpha("black", 0.5)))



### If wanting to aggregate into seasons:

# May is Fall
# Jun/Jul/Aug are Winter
# Sep/Oct/Nov are Spring
# Dec/Jan are Summer

season.ind<- c("Fall", rep("Winter",3), rep("Spring",3), rep("Summer",2))
green.season<- stackApply(mosaic.green, season.ind, fun = mean)
names(green.season)<- c("Fall","Winter","Spring","Summer")
rasterVis::levelplot(green.season, at=breaks, col.regions=cols, main="Green") +
  layer(sp.points(dat, cex=0.1, pch = 16, col = alpha("black", 0.5)))




### Try plotting in ggplot2 so points can be plotted by month and season

mosaic.green.df<- as.data.frame(mosaic.green, xy=TRUE)
mosaic.green.df<- pivot_longer(mosaic.green.df, cols = -c(x,y), names_to = "month",
                                values_to = "green")
mosaic.green.df$month<- factor(mosaic.green.df$month, levels = month.abb[c(5:12,1)])

dat<- as.data.frame(dat)


ggplot() +
  geom_raster(data = mosaic.green.df, aes(x, y, fill = green)) +
  geom_point(data = dat, aes(x, y, color = id), size = 0.5, alpha = 0.75) + 
  theme_bw() +
  scale_fill_gradientn(colors = cols, na.value = "transparent", limits = c(-2000,5000)) +
  scale_color_brewer(palette = "Dark2") +
  coord_fixed() +
  facet_wrap(~ month) +
  theme(strip.text = element_text(size = 12, face = "bold"))




green.season.df<- as.data.frame(green.season, xy=TRUE)
green.season.df<- pivot_longer(green.season.df, cols = -c(x,y), names_to = "season",
                                values_to = "green")
green.season.df$season<- factor(green.season.df$season,
                                 levels = c("Fall","Winter","Spring","Summer"))



ggplot() +
  geom_raster(data = green.season.df, aes(x, y, fill = green)) +
  geom_point(data = dat, aes(x, y, color = id), size = 0.5, alpha = 0.75) + 
  theme_bw() +
  scale_fill_gradientn(colors = cols, na.value = "transparent", limits = c(-2000,5000)) +
  scale_color_brewer(palette = "Dark2") +
  coord_fixed() +
  facet_wrap(~ season) +
  theme(strip.text = element_text(size = 12, face = "bold"))








#### Wetness ####

setwd("~/Documents/Snail Kite Project/Data/R Scripts/acceleration/TC_Wetness")

#load files as raster brick
wet<- list.files(getwd(), pattern = "*.tif$")
wet.stack<- stack(wet)
wet.stack2<- flip(wet.stack, direction = 'y')
names(wet.stack2)<- names(wet.stack)
wet.brick<- wet.stack2

#change values x == 0 to NA (these are masked pixels)
raster::values(wet.brick)[raster::values(wet.brick) == 0] <- NA

#rename images by date
names(wet.brick)<- gsub(pattern = "LC08_22.07._", replacement = "", names(wet.brick))


#visualize original rasters
breaks<- seq(-7000, 200, by=10)
cols<- colorRampPalette(RColorBrewer::brewer.pal(9, 'Blues'))(length(breaks)-1)
rasterVis::levelplot(wet.brick[[1:6]], at=breaks, col.regions=cols, main="Wetness")




### Mosaic RasterBrick over space and time

mosaic.wet<- list()
tmp1<- gsub(pattern = "LC08_22.07._", replacement = "", names(wet.stack2))
ind.m<- format(as.Date(tmp1, format = "%Y%m%d"), format = "%m") %>% 
  as.numeric()
ord.months<- c(5:12,1)


for (i in 1:length(unique(ind.m))) {
  ind<- which(ind.m %in% ord.months[i])  #using unique(ind.m) instead of ord.months results in wrong order
  wet.sub<- wet.stack2[[ind]]
  
  # time1<- gsub(pattern = "LC08_22.07._", replacement = "", names(wet.stack))[i]
  ind.pr<- stringr::str_sub(names(wet.stack2)[ind], 6, 11)
  
  if (length(unique(ind.pr)) > 1) {
    tile.list<- vector("list", length = length(unique(ind.pr)))
    
    for (j in 1:length(unique(ind.pr))) {
      ind.tiles<- grep(pattern = unique(ind.pr)[j], names(wet.stack2)[ind])
      tile.list[[j]]<- mean(wet.sub[[ind.tiles]], na.rm=TRUE)
    }
    
    tile.list$fun <- mean
    tile.list$na.rm <- TRUE
    
    tile.mosaic<- do.call(raster::mosaic, tile.list)
    
  } else {
    tile.mosaic<- mean(wet.sub, na.rm=TRUE)
  }
  
  mosaic.wet[[i]]<- tile.mosaic
}


coordinates(dat) <- ~x+y


mosaic.wet<- raster::brick(mosaic.wet)
names(mosaic.wet)<- month.abb[c(5:12,1)]
rasterVis::levelplot(mosaic.wet, at=breaks, col.regions=cols, main="Wet") +
  layer(sp.points(dat, cex=0.1, pch = 16, col = alpha("black", 0.5)))



### If wanting to aggregate into seasons:

# May is Fall
# Jun/Jul/Aug are Winter
# Sep/Oct/Nov are Spring
# Dec/Jan are Summer

season.ind<- c("Fall", rep("Winter",3), rep("Spring",3), rep("Summer",2))
wet.season<- stackApply(mosaic.wet, season.ind, fun = mean)
names(wet.season)<- c("Fall","Winter","Spring","Summer")
rasterVis::levelplot(wet.season, at=breaks, col.regions=cols, main="Wet") +
  layer(sp.points(dat, cex=0.1, pch = 16, col = alpha("black", 0.5)))




### Try plotting in ggplot2 so points can be plotted by month and season

mosaic.wet.df<- as.data.frame(mosaic.wet, xy=TRUE)
mosaic.wet.df<- pivot_longer(mosaic.wet.df, cols = -c(x,y), names_to = "month",
                               values_to = "wet")
mosaic.wet.df$month<- factor(mosaic.wet.df$month, levels = month.abb[c(5:12,1)])

dat<- as.data.frame(dat)


ggplot() +
  geom_raster(data = mosaic.wet.df, aes(x, y, fill = wet)) +
  geom_point(data = dat, aes(x, y, color = id), size = 0.5, alpha = 0.75) + 
  theme_bw() +
  scale_fill_gradientn(colors = cols, na.value = "transparent", limits = c(-5000,250)) +
  scale_color_brewer(palette = "Dark2") +
  coord_fixed() +
  facet_wrap(~ month) +
  theme(strip.text = element_text(size = 12, face = "bold"))




wet.season.df<- as.data.frame(wet.season, xy=TRUE)
wet.season.df<- pivot_longer(wet.season.df, cols = -c(x,y), names_to = "season",
                               values_to = "wet")
wet.season.df$season<- factor(wet.season.df$season,
                                levels = c("Fall","Winter","Spring","Summer"))



ggplot() +
  geom_raster(data = wet.season.df, aes(x, y, fill = wet)) +
  geom_point(data = dat, aes(x, y, color = id), size = 0.5, alpha = 0.75) + 
  theme_bw() +
  scale_fill_gradientn(colors = cols, na.value = "transparent", limits = c(-5000,250)) +
  scale_color_brewer(palette = "Dark2") +
  coord_fixed() +
  facet_wrap(~ season) +
  theme(strip.text = element_text(size = 12, face = "bold"))




### Export aggregated Tasseled Cap Transformation data

setwd("~/Documents/Snail Kite Project/Data/R Scripts/acceleration")


## Monthly

#if saving individual layers
# writeRaster(mosaic.bright, filename=names(mosaic.bright), bylayer=TRUE, format="GTiff")

#if saving RasterBrick as raster format
writeRaster(mosaic.bright, filename = 'GiantArm_tcbright_monthly.grd', format="raster",
            overwrite=TRUE)
writeRaster(mosaic.green, filename = 'GiantArm_tcgreen_monthly.grd', format="raster",
            overwrite=TRUE)
writeRaster(mosaic.wet, filename = 'GiantArm_tcwet_monthly.grd', format="raster",
            overwrite=TRUE)

#if saving RasterBrick as geoTIFF
# writeRaster(mosaic.bright, filename='GiantArm_tcbright_monthly.tif', format="GTiff",
#             overwrite=TRUE, options=c("INTERLEAVE=BAND","COMPRESS=LZW"))



## Season

#if saving individual layers
# writeRaster(bright.season, filename=names(bright.season), bylayer=TRUE, format="GTiff")

#if saving RasterBrick as raster format
writeRaster(bright.season, filename = 'GiantArm_tcbright_season.grd', format="raster",
            overwrite=TRUE)
writeRaster(green.season, filename = 'GiantArm_tcgreen_season.grd', format="raster",
            overwrite=TRUE)
writeRaster(wet.season, filename = 'GiantArm_tcwet_season.grd', format="raster",
            overwrite=TRUE)

#if saving RasterBrick as geoTIFF
# writeRaster(bright.season, filename='GiantArm_tcbright_season.tif', format="GTiff",
#             overwrite=TRUE, options=c("INTERLEAVE=BAND","COMPRESS=LZW"))