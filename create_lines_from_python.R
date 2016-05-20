# Load libraries
packages <- c('rgeos', 'rgdal', 'raster', 'bfastSpatial', 'sp', 'igraph', 'spatstat', 'maptools')
lapply(packages, require, character.only = TRUE)


## LOAD PYTHON RESULTS ##
ras <- raster('output/vegNA.tif')
epsg <- crs(ras)

# setwd
setwd('D:/Sugarcane_Project/201601_Sugar_Bacolod_sugarcanfields_zone_1/orthomosaics/')

linecoords <- read.table('rotatedlines.txt', sep = '')

## EXTRACT LINES ##
begin.coord <- data.frame(lon=c(linecoords[,1]), lat=c(linecoords[,2]))
end.coord <- data.frame(lon=c(linecoords[,3]), lat=c(linecoords[,4]))

## raw list to store Lines objects
rawlist <- vector("list", nrow(begin.coord))

for (i in seq_along(rawlist)) {
  rawlist[[i]] <- Lines(list(Line(rbind(begin.coord[i,], end.coord[i,]))), as.character(i))
}

splines <- SpatialLines(rawlist, proj4string = epsg)
splinesdf <- SpatialLinesDataFrame(splines, linecoords, match.ID = FALSE)


# write a shapefile
writeOGR(splinesdf, getwd(),
         "output/f3test_all", driver="ESRI Shapefile")