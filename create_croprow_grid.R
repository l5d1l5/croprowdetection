library(rgeos)
library(rgdal)
library(raster)


###############
## LOAD DATA ##
###############
setwd('D:/Sugarcane_Project/201601_Sugar_Bacolod_sugarcanfields_zone_1/orthomosaics/')

# Load binary raster
vegbin <- raster('output/vegNA.tif')

# Load point coordinates from Python
linesNptxt <- read.table('output/coordinates.txt', sep = '') # load numpy array
linesNptxt <- round(linesNptxt, digits = 2)
colnames(linesNptxt) <- c('p1x','p1y','p2x','p2y')
linesNptxt <- transform(linesNptxt, id = as.numeric(interaction(p1x, p2y, drop=TRUE)))

#####################
## CALCULATE ANGLE ##
#####################

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
  
  return(radconv)
}

# Set parameters
line_length <- 250;
angle <- anglecalc(linesNptxt)
epsg <- crs(vegbin)
line_amount <- 60

#####################
## INCREASE EXTENT ##
#####################

extentincrease <- function(percentage, raster) {
  xminnew <- xmin(vegbin) - ((xmin(vegbin) / 1000) * percentage)
  xmaxnew <- xmax(vegbin) + ((xmax(vegbin) / 1000) * percentage)
  
  yminnew <- ymin(vegbin) - ((ymin(vegbin) / 1000) * percentage)
  ymaxnew <- ymax(vegbin) + ((ymax(vegbin) / 1000) * percentage)
  
  ex <- extent(xminnew, xmaxnew, ymin(raster), ymaxnew)
  area <- setExtent(raster, ex, keepres=TRUE, snap=TRUE)
  
  return(area)
}

newex <- extentincrease(0.20, vegbin)

#################
## CREATE GRID ##
#################

# Create grid function
gridcreate <- function(rasterextent, epsg, line_amount, line_length, angle, spacing) {
  
  #initialize list to to 2 empty matrices of 2 by 3
  clist <-  matrix(NA, nrow=line_amount, ncol=2)
  names(clist) <- c('x','y')
  
  # Rad 2 degreee
  rad2deg <- function(rad) {(rad * 180) / (pi)}
  angledeg <- rad2deg(angle)
  angledegsup <- 180 + angledeg
  
  if(angledegsup <= 210 & angledegsup >= 150){
    
    # Create y points
    for(i in 1:line_amount){
      clist[[i,2]]<- (ymin(rasterextent) + (i*spacing))
    }
    
    # Create x values
    for (i in 1:line_amount){
      clist[[i,1]] <- xmin(rasterextent) 
    }} else {
      for(i in 1:line_amount){
        clist[[i,2]]<- ymin(rasterextent)
      }
      for (i in 1:line_amount){
        clist[[i,1]] <- (xmin(rasterextent) + (i*spacing))
      }
    }
  
  # Generate line
  x1 = clist[,1]
  y1 = clist[,2]
  x2 = x1 + line_length * cos(angle)
  y2 = y1 + line_length * sin(angle)
  
  begin.coord <- cbind(x1, y1)
  end.coord <- cbind(x2, y2)
  
  points <- rbind(begin.coord, end.coord)
  colnames(points) <- c('x','y')
  points <- transform(points, id = as.numeric(interaction(x, y, drop=TRUE)))
  
  ## raw list to store Lines objects
  rawlist <- vector("list", nrow(begin.coord))
  
  for (i in seq_along(rawlist)) {
    rawlist[[i]] <- Lines(list(Line(rbind(begin.coord[i, ], end.coord[i,]))), as.character(i))
  }
  
  splines <- SpatialLines(rawlist, proj4string = epsg)
  data <- points[1:line_amount,]
  splinesdf <- SpatialLinesDataFrame(splines, data, match.ID = FALSE)
  
  return(splinesdf)
}

grid <- gridcreate(newex, epsg, 145, 225, angle, 1.8)


# write a shapefile
writeOGR(grid, getwd(),
         "temp/60_linesnewex", driver="ESRI Shapefile")

###################
## CREATE BUFFER ##
###################

system(paste('ogr2ogr.exe "[temporary file]" D:/Sugarcane_Project/201601_Sugar_Bacolod_sugarcanfields_zone_1/orthomosaics/temp/60_lines.shp 60_lines -dialect sqlite -sql "SELECT ST_Buffer( geometry , 0.4 ),* FROM 60_lines'))

buf <- readOGR('temp/buffer.shp', 'buffer')
ndvi <- raster('ndvi.tif')

rowsndvi <- mask(ndvi, buf)
# Write to file
writeRaster(rowsndvi, 'output/rowsndvi.tif', overwrite=TRUE, prj=TRUE, format = 'GTiff',
            options=c("COMPRESS=NONE", "TFW=YES"))

# Inspect row histoggram to determine treshold value for soil pixels
par(mfrow=c(1,2)) 
plot(rowsndvi, main = "Map")
hist(rowsndvi, main = "NDVI frequency")

gaps <- rowsndvi <= -0.12
gaps[gaps == 0] <- NA

# Write to file
writeRaster(gaps, 'output/gaps.tif', overwrite=TRUE, prj=TRUE, format = 'GTiff',
            options=c("COMPRESS=NONE", "TFW=YES"))
