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
  
  list <- list('Rads' = radconv, 'Degree' = angle)
  
  return(list)
}

## Set parameters ##
line_length <- 120;
angle <- as.numeric(anglecalc(linesNptxt)[1])
epsg <- crs(ndvi)
line_amount <- 250
spacing <- 0.045

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
  
  # Generate begin coords based on degree
  if(angledegsup <= 210 & angledegsup >= 150){
    
    # Create y points
    for(i in 1:line_amount){
      clist[[i,2]]<- (ymin(ndvi) + (i*spacing))
    }
    # Create x values
    for (i in 1:line_amount){
      clist[[i,1]] <- xmin(ndvi) 
    }} else {
      
      # Create x values
      for(i in 1:line_amount){
        clist[[i,2]]<- ymin(ndvi)
      }
      # Create y values
      for (i in 1:line_amount){
        clist[[i,1]] <- (xmin(ndvi) + (i*spacing))
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
  
  # set rotatations for lines
  #rotatelist <- list(-1,-0.5,0,0.5,1)
  
  # Rotate spatial lines and unlist
  #gridElide <- lapply(rotatelist, function(x) elide(splines, rotate = x, center=apply(bbox(splines), 1, mean)))
  
  #ll0 <- lapply(gridElide , function(x) `@`(x , "lines"))
  #ll1 <- lapply( unlist( ll0 ) , function(y) `@`(y,"Lines"))
  #Sl <- SpatialLines( list( Lines( unlist( ll1 ) , ID = 1 ))) # ERROR: makes 1 feature of gridElide
  #proj4string(Sl) <- epsg
  
  
  # Extract NDVI values
  #gridVal <- lapply(gridElide, function(x) extract(ndvi, x, method = 'simple', fun = 'sum', na.rm = TRUE, df = TRUE, sp = FALSE))
  
  data <- points[1:line_amount,]
  splinesdf <- SpatialLinesDataFrame(splines, data, match.ID = FALSE)
  
  splinesdfVal <- extract(ndvi, splinesdf, method = 'simple', fun = 'sum', na.rm = TRUE, df = TRUE, sp = TRUE)
  
  return(splinesdfVal)
}

system.time(grid <- gridcreate(ndvi, epsg, line_amount, line_length, angle, spacing))

# write a shapefile
writeOGR(grid, getwd(),
         "temp/scanlines_splinesdfVal", driver="ESRI Shapefile")


# Subset grid based on high NDVI values (top 10% of data)
q <- quantile(grid$ndvi, 0.90)
grid_croprows <- grid[grid$ndvi >= q,]

# write a shapefile
writeOGR(grid_croprows, getwd(),
         "temp/scanlines_gridcroprows", driver="ESRI Shapefile")

# Apply elide on crop rows in order to test other angles
# set rotatations for lines
rotatelist <- list(-0.2, -0.1, 0.1, 0.2)

# Rotate spatial lines and unlist
gridElide <- lapply(rotatelist, function(x) elide(grid_croprows, rotate = x, center=apply(bbox(grid_croprows), 1, mean)))
# Extract ndvi from rotated lines
gridValdf <- lapply(gridElide, function(x) extract(ndvi, x, method = 'simple', fun = 'sum', na.rm = TRUE, df = TRUE, sp = FALSE))
gridValsp <- lapply(gridElide, function(x) extract(ndvi, x, method = 'simple', fun = 'sum', na.rm = TRUE, df = TRUE, sp = TRUE))

# Combine ndvi of angles to ndvi of straight line and check highest value
anglesndvi <- as.data.frame(gridValdf)
keep <- c('ndvi', 'ndvi.1', 'ndvi.2', 'ndvi.3')
anglesndvi <- anglesndvi[keep]
straightndvi <- as.data.frame(grid_croprows$ndvi)
allanglesndvi <- cbind(straightndvi, anglesndvi)
colnames(allanglesndvi) <- c('0', '-0.2','-0.1','0.1','0.2')

# Compare rotated angles with grid_croprows
maxndvi <- colnames(allanglesndvi)[apply(allanglesndvi,1,which.max)]
allanglesndvi['MaxNDVIangle'] <- maxndvi

croprowone <- elide(grid_croprows[1])


# write a shapefile
writeOGR(grid, getwd(),
         "temp/scanlines_test_5000_0002", driver="ESRI Shapefile")

# Calculate the number of cores
#no_cores <- detectCores() - 1
# Initiate cluster
#cl <- makeCluster(no_cores)
#linevalMC <- parLapply(cl, ndvi, fun = extract(ndvi, grid, method = 'simple', fun = 'sum', na.rm = TRUE, df = TRUE))

cropmintwo <- elide(grid_croprows, rotate = -0.2, center=apply(bbox(grid_croprows), 1, mean))
proj4string(cropmintwo) <- epsg
cropminone <- elide(grid_croprows, rotate = -0.1, center=apply(bbox(grid_croprows), 1, mean))
proj4string(cropminone) <- epsg

# write a shapefile
writeOGR(cropminone, getwd(),
         "temp/cropminone", driver="ESRI Shapefile")


#----------------------------------------------------------#
rotatelist <- list(-1,-0.5,0,0.5,1)

gridElide <- lapply(rotatelist, function(x) elide(splines, rotate = x, center=apply(bbox(splines), 1, mean)))

#  Get the Lines objects which contain multiple 'lines'
ll0 <- lapply(gridElide , function(x) `@`(x , "lines") )

#  Extract the individual 'lines'
ll1 <- lapply( unlist( ll0 ) , function(y) `@`(y,"Lines") )

#  Combine them into a single SpatialLines object
Sl <- SpatialLines( list( Lines( unlist( ll1 ) , ID = 1)))
proj4string(Sl) <- epsg

sldf <- SpatialLinesDataFrame(Sl, data = gridValsp, match.id = FALSE)

gridVal <- lapply(gridElide, function(x) extract(ndvi, x, method = 'simple', fun = 'sum', na.rm = TRUE, df = TRUE, sp = TRUE))

# write a shapefile
writeOGR(gridElide, getwd(),
         "temp/scanlines_gridElide", driver="ESRI Shapefile")

#-------------------------------------------------------------#













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