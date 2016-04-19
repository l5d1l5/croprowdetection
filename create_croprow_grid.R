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
  radavg <- mean(rads)
  radconv <- pi - (-radavg)
  radconv <- pi - radconv
  
  return(radconv)
}

# Set parameters
line_length <- 250;
angle <- anglecalc(linesNptxt)
epsg <- crs(vegbin)
line_amount <- 60

#################
## CREATE GRID ##
#################

# Create grid function
gridcreate <- function(rasterextent, epsg, line_amount, line_length, angle, spacing) {
  
  #initialize list to to 2 empty matrices of 2 by 3
  clist <-  matrix(NA, nrow=line_amount, ncol=2)
  names(clist) <- c('x','y')
  
  # Create x points
  for(i in 1:line_amount){
    clist[[i,2]]<-ymin(rasterextent)
  }
  
  # Create y values
  for (i in 1:line_amount){
    clist[[i,1]] <- (xmin(rasterextent) + (i*spacing))
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

grid <- gridcreate(vegbin, epsg, 45, 225, angle, 1.8)


# write a shapefile
writeOGR(grid, getwd(),
         "temp/60_lines", driver="ESRI Shapefile")

###################
## CREATE BUFFER ##
###################

system(paste('ogr2ogr.exe "[temporary file]" D:/Sugarcane_Project/201601_Sugar_Bacolod_sugarcanfields_zone_1/orthomosaics/temp/60_lines.shp 60_lines -dialect sqlite -sql "SELECT ST_Buffer( geometry , 0.4 ),* FROM 60_lines'))


