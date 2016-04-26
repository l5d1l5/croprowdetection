setwd('D:/Sugarcane_Project/201601_Sugar_Bacolod_sugarcanfields_zone_1/orthomosaics/')
ras <- raster('output/shapefield_sample_ndvi_ren.tif')

# Load centroidpoint coordinates
centroids <- read.table('output/centroidscord.txt', sep = '') # load numpy array
centroids <- round(centroids, digits = 4)
colnames(centroids) <- c('x','y')
# Create random ID
centroids <- transform(centroids, id = as.numeric(interaction(x, y, drop=TRUE)))

y = -centroids[,1]
x = centroids[,2]

xy <- cbind(x,y)

# Convert to spatialpoints
#epsg <- crs(ras)
centroid_points <- SpatialPointsDataFrame(xy,
                                       centroids,    #the R object to convert
                                       proj4string = epsg)   # assign a CRS

# write a shapefile
writeOGR(centroid_points, getwd(),
         "output/centroid_points", driver="ESRI Shapefile")

#-------------------------------------------------------------------#
# Load line point coordinates to calculate avarage angle
linesNptxt <- read.table('output/coordinates.txt', sep = '') # load numpy array
linesNptxt <- round(linesNptxt, digits = 4)
colnames(linesNptxt) <- c('p1x','p1y','p2x','p2y')


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

angle <- anglecalc(linesNptxt)
line_length <- 550

#-------------------------#
# create lines based on point centroids

# Generate line
x1 = centroids[,2]
y1 = centroids[,1]
x2 = x1 + line_length * cos(-angle)
y2 = y1 + line_length * sin(-angle)

begin.coord <- cbind(x1, -y1)
end.coord <- cbind(x2, -y2)

points <- rbind(begin.coord, end.coord)
colnames(points) <- c('x','y')
points <- transform(points, id = as.numeric(interaction(x, y, drop=TRUE)))

## raw list to store Lines objects
rawlist <- vector("list", nrow(begin.coord))

for (i in seq_along(rawlist)) {
  rawlist[[i]] <- Lines(list(Line(rbind(begin.coord[i, ], end.coord[i,]))), as.character(i))
}

splines <- SpatialLines(rawlist, proj4string = epsg)
data <- points[1:nrow(centroids),]
splinesdf <- SpatialLinesDataFrame(splines, data, match.ID = FALSE)


# write a shapefile
writeOGR(splinesdf, getwd(),
         "temp/linesfromcetrnoids", driver="ESRI Shapefile")



