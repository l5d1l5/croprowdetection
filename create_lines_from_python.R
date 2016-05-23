options(warn=-1)

# Load libraries
packages <- c('rgdal', 'raster', 'sp')
lapply(packages, require, character.only = TRUE)

# setwd
#setwd('D:/Sugarcane_Project/201601_Sugar_Bacolod_sugarcanfields_zone_1/orthomosaics/')

# Fetch command line arguments
myArgs <- commandArgs(trailingOnly = TRUE)

## LOAD PYTHON RESULTS ##
ras <- raster(paste(myArgs,'/ndvi_f6.tif', sep = ''))
epsg <- crs(ras)

linecoords <- read.table(paste(myArgs,'/rotatedlines.txt', sep = ''))

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
writeOGR(splinesdf, (paste(myArgs, sep = '')), '/croprows', driver="ESRI Shapefile")

cat('Shapefile created')