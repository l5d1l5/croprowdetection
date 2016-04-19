library(rgeos)
library(rgdal)
library(raster)
library(bfastSpatial)
library(sp)
library(igraph)
library(RcppCNPy)
library(spatstat)
library(maptools)

## SET WORKING DIRECTORY AND SETTINGS ##
setwd('D:/Sugarcane_Project/201601_Sugar_Bacolod_sugarcanfields_zone_1/orthomosaics/')
# Create output dir
dir.create('output')

## LOAD SUGARCANE FIELDS ##
ras <- stack('Sugar_Bacolod_sugarcanefields_onefield_ortho.tif')

## CACLULATE VEGETATION INDICES ##
# NDVI calculator
ViCalc <- function(x, y) {
  ndvi <- (y - x) / (y + x)
  return(ndvi)
}

# Calculate NDVI
ndvi <- overlay(x = ras$Sugar_Bacolod_sugarcanefields_onefield_ortho.1, 
                y = ras$Sugar_Bacolod_sugarcanefields_onefield_ortho.4, fun = ViCalc)
names(ndvi) <- "NDVI"

# Inspect NDVI histoggram to determine treshold value for soil pixels
par(mfrow=c(1,2)) 
plot(ndvi, main = "Map")
hist(ndvi, main = "NDVI frequency")

tresholdmax = 0.10
areatreshold = 0.02

#--------------------------------------------------------------------------------------#
## IMAGE SEGMENTATION ##

# Veg extract function
VegExtract <- function(inputlayer) {
  veg <- inputlayer >= tresholdmax
  veg[veg == 0] <- 0
  veg <- areaSieve(veg, thresh = areatreshold, directions = 4)
  return(veg)
}

# Extract the vegetation from NDVI layer
veg <- VegExtract(ndvi)

# Plot extracted vegetation pixels
par(mfrow=c(1,1)) 
objectsLegend1 <- c("Veg Pixels")
plot(veg, col = 'Black', legend=FALSE, main = 'Plot of extracted veg pixels') 
legend("bottomright", legend=objectsLegend1, cex=0.9, fill='Green', ncol=1)


# Write to file
writeRaster(veg, 'output/vegNA.tif', overwrite=TRUE, prj=TRUE, format = 'Gtiff',
            options=c("COMPRESS=NONE", "TFW=YES, PROFILE=BASELINE"))

system(paste('-t_srs +proj=utm +zone=51 +datum=WGS84 gdalwarp output/vegNA.tif output/vegNAout.tif'))



## Manual re-saving in Qgis still needed to save tif as a rendered image in order for it to load into python line detection module ##

# Run python script
system(paste("python D:/Sugarcane_Project/Scripts/HoughLineDetection.py"))

## LOAD PYTHON RESULTS ##
# Load lines tiff
detectedlines <- raster('croprowsP.tif')
detectedlines_bin <- detectedlines > 0.1 # binary detected lines

writeRaster(detectedlines_bin, 'output/lines_bin.tif', overwrite = TRUE) # write lines to file

# Load point coordinates
linesNptxt <- read.table('output/coordinates.txt', sep = '') # load numpy array
linesNptxt <- round(linesNptxt, digits = 2)
colnames(linesNptxt) <- c('p1x','p1y','p2x','p2y')
linesNptxt <- transform(linesNptxt, id = as.numeric(interaction(p1x, p2y, drop=TRUE)))

# Convert to spatialpoints
epsg <- crs(ras)
linesPoints1 <- SpatialPointsDataFrame(linesNptxt[,1:2],
                                                linesNptxt,    #the R object to convert
                                                proj4string = epsg)   # assign a CRS
linesPoints2 <- SpatialPointsDataFrame(linesNptxt[,3:4],
                                       linesNptxt,    #the R object to convert
                                       proj4string = epsg)   # assign a CRS 

points <- rbind(linesPoints1, linesPoints2)

# write a shapefile
writeOGR(points, getwd(),
         "output/all_points", driver="ESRI Shapefile")

#--------------------------------------------------------------------------------------#
## EXTRACT LINES ##
begin.coord <- data.frame(lon=c(-linesNptxt[,2]), lat=c(linesNptxt[,1]))
end.coord <- data.frame(lon=c(-linesNptxt[,4]), lat=c(linesNptxt[,3]))

## raw list to store Lines objects
rawlist <- vector("list", nrow(begin.coord))

for (i in seq_along(rawlist)) {
  rawlist[[i]] <- Lines(list(Line(rbind(begin.coord[i, ], end.coord[i,]))), as.character(i))
}

splines <- SpatialLines(rawlist, proj4string = epsg)
splinesdf <- SpatialLinesDataFrame(splines, linesNptxt, match.ID = TRUE)

# write a shapefile
writeOGR(splinesdf, getwd(),
         "output/all_lines", driver="ESRI Shapefile")

