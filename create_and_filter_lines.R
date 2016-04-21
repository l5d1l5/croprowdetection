## LOAD PYTHON RESULTS ##
# Load lines tiff
detectedlines <- raster('croprowsP.tif')
detectedlines_bin <- detectedlines > 0.1 # binary detected lines

writeRaster(detectedlines_bin, 'output/lines_bin.tif', overwrite = TRUE) # write lines to file

#--------------------------------------------------------------------------------------#
# Load point coordinates
linesNptxt <- read.table('output/coordinates.txt', sep = '') # load numpy array
linesNptxt <- round(linesNptxt, digits = 4)
colnames(linesNptxt) <- c('p1x','p1y','p2x','p2y')

#--------------------------------------------------------------------------------------#
# Add angles
dx <- linesNptxt[,3] - linesNptxt[,1]
dy <- linesNptxt[,4] - linesNptxt[,2]

rads <- atan2(dy,dx)

linesNptxt <- cbind(linesNptxt, rads)

# Create random ID
linesNptxt <- transform(linesNptxt, id = as.numeric(interaction(p1x, p2y, drop=TRUE)))

# Filter out extremes on both ends based on percentiles
q1 <- quantile(rads, 0.10)
q2 <- quantile(rads, 0.90)

linesNptxt <- subset(linesNptxt, rads >= q1 & rads <= q2)

#--------------------------------------------------------------------------------------#

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
begin.coord <- data.frame(lon=c(linesNptxt[,1]), lat=c(-linesNptxt[,2]))
end.coord <- data.frame(lon=c(linesNptxt[,3]), lat=c(-linesNptxt[,4]))

## raw list to store Lines objects
rawlist <- vector("list", nrow(begin.coord))

for (i in seq_along(rawlist)) {
  rawlist[[i]] <- Lines(list(Line(rbind(begin.coord[i, ], end.coord[i,]))), as.character(i))
}

splines <- SpatialLines(rawlist, proj4string = epsg)
splinesdf <- SpatialLinesDataFrame(splines, linesNptxt, match.ID = FALSE)

# write a shapefile
writeOGR(splinesdf, getwd(),
         "output/all_lines", driver="ESRI Shapefile")

#--------------------------------------------------------------------------------------#
# Merge lines

linesdf <- as.data.frame(splinesdf)

linesmerg <- gLineMerge(splinesdf, byid=FALSE)
