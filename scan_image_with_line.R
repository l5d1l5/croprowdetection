# Linescan over image

# Load libraries
packages <- c('rgeos', 'rgdal', 'raster', 'bfastSpatial', 'sp', 'igraph', 'spatstat', 'maptools')
lapply(packages, require, character.only = TRUE)

# setwd
setwd('D:/Sugarcane_Project/201601_Sugar_Bacolod_sugarcanfields_zone_1/orthomosaics/')
# Create output dir
dir.create('output')

## LOAD SUGARCANE FIELD ##
ras <- stack('Sugar_Bacolod_sugarcanefields_onefield_ortho_square.tif')

## CACLULATE VEGETATION INDEX ##
# NDVI calculator
ViCalc <- function(x, y) {
  ndvi <- (y - x) / (y + x)
  return(ndvi)
}

# Calculate NDVI and write to raster
ndvi <- overlay(x = ras$Sugar_Bacolod_sugarcanefields_onefield_ortho_square.1, 
                y = ras$Sugar_Bacolod_sugarcanefields_onefield_ortho_square.4, fun = ViCalc)
names(ndvi) <- "NDVI"

writeRaster(ndvi, 'output/ndvi.tif', overwrite=TRUE, prj=TRUE, format = 'GTiff',
            options=c("COMPRESS=NONE", "TFW=YES, PROFILE=BASELINE"))

## Manual re-saving in Qgis still needed to save tif as a rendered image in order for it to load into 
## python line detection module

# Run Python script
system(paste("python C:/Users/darellvdv/Documents/Python_Scripts/generate_houghlines.py"))

# Load Python output
ndvi <- raster('output/ndvi.tif')
linesNptxt <- read.table('output/coordinates.txt', sep = '') # load numpy array
linesNptxt <- round(linesNptxt, digits = 2)
colnames(linesNptxt) <- c('p1x','p1y','p2x','p2y')
linesNptxt <- transform(linesNptxt, id = as.numeric(interaction(p1x, p2y, drop=TRUE)))

# Angle calculator from hough transform numpy array
anglecalc <- function(numpy_array){
  
  # Calculate angle
  dx <- numpy_array[,3] - numpy_array[,1]
  dy <- numpy_array[,4] - numpy_array[,2]
  
  rads <- atan2(dy,dx)
  
  # Select rads based on quantile (80% of data) and avarage
  q1 <- quantile(rads, 0.10)
  q2 <- quantile(rads, 0.90)
  
  radssub <- as.data.frame(subset(rads, rads >= q1 & rads <= q2))
  radavg <- mean(radssub[,1])
  
  # Convert rad
  radconv <- pi - (-radavg)
  radconv <- pi - radconv
  
  rad2deg <- function(rad) {(rad * 180) / (pi)}
  angle <- rad2deg(radconv)
  
  return(angle)
}

# Set Parameters
angle <- anglecalc(linesNptxt)

coord1ymax <- extent(ndvi)@ymax
coord1xmin <- extent(ndvi)@xmin
coord1 <- cbind(coord1ymax, coord1xmin)

coord2ymin <- extent(ndvi)@ymin
coord2xmin <- extent(ndvi)@xmin
coord2 <- cbind(coord2ymin, coord2xmin)

xmin <- extent(ndvi)@xmin
xmax <- extent(ndvi)@xmax

ymin <- extent(ndvi)@ymin
ymax <- extent(ndvi)@ymax

thresh <- 0.5

for(i in seq((angle - 5), (angle + 5), 0.2)){
  for(j in seq(xmin, xmax)){
    
  }
}
    






for(i in seq((angle -5),(angle +5),0.2){
  for j in seq(xmin, xmax){
    # shift lane
    extract(ndvi,line,fun=sum)
    if: (> thresh) { break
      }
  }
}
      line <- elide(line,center, rotate= i)
      ((j, (0.5 + ymax + ymin)))
    )