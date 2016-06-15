options(warn=-1)

# Load libraries
packages <- c('rgdal', 'raster', 'sp', 'rgeos')
suppressMessages(lapply(packages, require, character.only = TRUE))

# setwd
#setwd('D:/Sugarcane_Project/201601_Sugar_Bacolod_sugarcanfields_zone_1/orthomosaics/')

# Fetch command line arguments
myArgs <- commandArgs(trailingOnly = TRUE)

print('loading Python results into R')

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


# Load temp data
#splinesdf <- readOGR("D:/Sugarcane_project/Ecuador_output", "croprows")
#ras <- raster("D:/Sugarcane_project/Ecuador_output/ndvi_f6.tif")

# Clip spatial lines to raster extent
splinesclipped <- crop(splinesdf, ras)

print("Spatial lines created!.. now generating buffer")

# Buffer spatial lines
buffer <- gBuffer(splinesdf, byid = TRUE, width = 0.5, quadsegs=5, capStyle="ROUND",
           joinStyle="ROUND")


writeOGR(buffer, (paste(myArgs, sep = '')), "/buffer", driver="ESRI Shapefile")
# write a shapefile
writeOGR(splinesclipped, (paste(myArgs, sep = '')), '/croprows', driver="ESRI Shapefile")

print("Calling GDAL from system to perfrom raster cutting")

# System call to GDAL to clip raster to buffer
system(paste("gdalwarp -q -cutline D:/Sugarcane_project/Ecuador_output/buffer.shp -crop_to_cutline -of GTiff D:/Sugarcane_project/Ecuador_output/ndvi_f6.tif D:/Sugarcane_Project/Ecuador_output/ndvi_f6_buffer5.tif"))
croprowsclip <- raster("D:/Sugarcane_Project/Ecuador_output/ndvi_f6_buffer5.tif")

print("Raster sucessfully clipped!")
print("Detecting crop gaps..")

# Set threshold to detect crop gaps inside rows
cropgaps_all <- croprowsclip <= 0
# Extend cropgaps and clump pixels
#cropgapsex <- extend(cropgaps_all, c(1,1))
cropgapsclump <- clump(cropgaps_all, directions = 8) 
# Build a frequancy table
fr <- as.data.frame(freq(cropgapsclump))
# Determine regions with less than 15 pixels clumped
excluderegions <- fr$value[which(fr$count <= 100)]

# make a new raster to be sieved
cropgaps <- cropgapsclump
# assign NA to all clumps whose IDs are found in excludeID
cropgaps[cropgapsclump %in% excluderegions] <- NA
cropgaps[cropgaps >= 1] <- 1

print("Crop gaps detected! Now wring gaps to raster file")

writeRaster(cropgaps, "D:/Sugarcane_Project/Ecuador_output/cropgaps.tif", format = "GTiff", overwrite = TRUE)

print('Gaps created')