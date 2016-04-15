library(rgeos)
library(rgdal)
library(raster)
library(bfastSpatial)
library(sp)
library(igraph)


setwd('D:/Sugarcane_Project/201601_Sugar_Bacolod_sugarcanfields_zone_1/orthomosaics/')
# Create output dir
dir.create('output')

# Read rasterstack
ras <- stack('Sugar_Bacolod_sugarcanefields_onefield_ortho.tif')


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
writeRaster(veg, 'output/vegNA.tif', format='GTiff', overwrite=TRUE, progress='window', options=c("COMPRESS=NONE", "TFW=YES"))

## Manual re-saving in Qgis still needed to save tif as a rendered image in order for it to load into python line detection module ##

# Run python script
system(paste("python D:/Sugarcane_Project/Scripts/HoughLineDetection.py"))

# Load results from Python script
detectedlines <- raster('croprowsP.tif')
detectedlines_bin <- detectedlines > 0.1 # binary detected lines

writeRaster(detectedlines_bin, 'output/lines_bin.tif', overwrite = TRUE)

