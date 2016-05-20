# -*- coding: utf-8 -*-
"""
Created on Mon May 09 10:11:46 2016

@author: darellvdv
"""
import numpy as np
import os
import math
from osgeo import gdal
from osgeo.gdalconst import GA_ReadOnly, GDT_Float32
from matplotlib import pyplot as plt
import cv2

from skimage.segmentation import clear_border
from skimage.measure import label
from skimage.measure import regionprops
from skimage.color import label2rgb
from skimage import img_as_ubyte

import scipy.misc
#import scipy.ndimage
#from scipy.stats import itemfreq

#------------------------------------------------------------------------------#
## STEP 1 LOAD SUGARCANE RASTER
#------------------------------------------------------------------------------#

# Set working directory
path = os.chdir('D:/Sugarcane_Project/201601_Sugar_Bacolod_sugarcanfields_zone_1/orthomosaics')

# Open Ortho and print metadata
filename = 'Sugar_Bacolod_sugarcanefields_onefield_ortho.tif'
dataSource = gdal.Open(filename, GA_ReadOnly)
if not dataSource:
    print "THE FOLLOWING RASTER FAILED TO LOAD:\n", filename
else:
    print "\nInformation about " + filename 
    print "Driver: ", dataSource.GetDriver().ShortName,"/", \
    dataSource.GetDriver().LongName
    print "Size is ",dataSource.RasterXSize,"x",dataSource.RasterYSize, \
    'x',dataSource.RasterCount
    
    print '\nProjection is: ', dataSource.GetProjection()
    
    print "\nInformation about the location of the image and the pixel size:"
    geotransform = dataSource.GetGeoTransform()
    if not geotransform is None:
        print 'Origin = (',geotransform[0], ',',geotransform[3],')'
        print 'Pixel Size = (',geotransform[1], ',',geotransform[5],')'
    
#------------------------------------------------------------------------------#
## STEP 2 CALCULATE NDVI
#------------------------------------------------------------------------------#
    
# Read data into an array
band1Arr = dataSource.GetRasterBand(1).ReadAsArray(0,0,dataSource.RasterXSize, dataSource.RasterYSize)
band4Arr = dataSource.GetRasterBand(4).ReadAsArray(0,0,dataSource.RasterXSize, dataSource.RasterYSize)
print type(band1Arr)

# set the data type
band1Arr=band1Arr.astype(np.float32)
band4Arr=band4Arr.astype(np.float32)

# Derive the NDVI
mask = np.greater(band1Arr+band4Arr,0)

# set np.errstate to avoid warning of invalid values (i.e. NaN values) in the divide 
with np.errstate(invalid='ignore'):
    ndvi = np.choose(mask,(float('NaN'),(band4Arr-band1Arr)/(band4Arr+band1Arr)))
print "NDVI min and max values", ndvi.min(), ndvi.max()
# Check the real minimum value
print ndvi[ndvi>-99].min() # gives error since turning -99 to NaN!

# Write the result to disk
driver = gdal.GetDriverByName('GTiff')
outDataSet=driver.Create('output/ndvi_f5.tif', dataSource.RasterXSize, dataSource.RasterYSize, 1, GDT_Float32)
outBand = outDataSet.GetRasterBand(1)
outBand.WriteArray(ndvi,0,0)
outBand.SetNoDataValue(float('NaN'))

# set the projection and extent information of the dataset
outDataSet.SetProjection(dataSource.GetProjection())
outDataSet.SetGeoTransform(dataSource.GetGeoTransform())

# Flush to save
outBand.FlushCache()
outDataSet.FlushCache()

#------------------------------------------------------------------------------#
## STEP 3 MANUAL RESAVING IN QGIS 
#------------------------------------------------------------------------------#

# Renders image, however contrast not good yet
maskedndvi = np.ma.masked_array(ndvi, np.isnan(ndvi)) # mask NaN value
#scipy.misc.imsave('output/ndvi_ren_test2.png', maskedndvi)
#image = scipy.misc.toimage(maskedndvi, low=0, high=255, cmin=(np.min(maskedndvi)), cmax=(np.max(maskedndvi))).save('output/ndvi_ren_f5.png')
image = scipy.misc.toimage(maskedndvi, low=0, high=255, cmin=(np.min(maskedndvi)), cmax=(np.max(maskedndvi)))

image = img_as_ubyte(image)
# Manual re-saving in qgis is still needed to get a correct rendered image -> look for solution

#------------------------------------------------------------------------------#
## STEP 4 INTERACTIVE CROPPING
#------------------------------------------------------------------------------#

# import the necessary packages
#import argparse

# initialize the list of reference points and boolean indicating
# whether cropping is being performed or not
refPt = []
cropping = False
 
def click_and_crop(event, x, y, flags, param):
    # grab references to the global variables
    global refPt, cropping
 
    # if the left mouse button was clicked, record the starting
    # (x, y) coordinates and indicate that cropping is being
    # performed
    if event == cv2.EVENT_LBUTTONDOWN:
        refPt = [(x, y)]
        cropping = True
 
    # check to see if the left mouse button was released
    elif event == cv2.EVENT_LBUTTONUP:
         # record the ending (x, y) coordinates and indicate that
         # the cropping operation is finished
         refPt.append((x, y))
         cropping = False
 
         # draw a rectangle around the region of interest
         cv2.rectangle(image, refPt[0], refPt[1], (0, 255, 0), 2)
         cv2.imshow("image", image)
  
# construct the argument parser and parse the arguments
#ap = argparse.ArgumentParser()
#ap.add_argument("-i", "--image", required=True, help="Path to the image")
#args = vars(ap.parse_args())

# load the image, clone it, and setup the mouse callback function
#image = cv2.imread('output/ndvi_ren_test2.png')
clone = image.copy()
cv2.namedWindow("image", cv2.WINDOW_NORMAL)
cv2.setMouseCallback("image", click_and_crop)
 
# keep looping until the 'q' key is pressed
while True:

    # display the image and wait for a keypress
    cv2.imshow("image", image)
    key = cv2.waitKey(1) & 0xFF
 
    # if the 'r' key is pressed, reset the cropping region
    if key == ord("r"):
        image = clone.copy()
 
    # if the 'c' key is pressed, break from the loop
    elif key == ord("c"):
        break
 
# if there are two reference points, then crop the region of interest
# from teh image and display it
if len(refPt) == 2:
    roi = clone[refPt[0][1]:refPt[1][1], refPt[0][0]:refPt[1][0]]
    cv2.imshow("ROI", roi)
    cv2.waitKey(0)

# close all open windows
cv2.destroyAllWindows()

cv2.imwrite('output/roi.tiff', roi)
#------------------------------------------------------------------------------#
## STEP 4 APPLY ADAPTIVE THRESHOLDING AND CLEANING &
## APPLY HOUGH LINES TRANSFORMATION
#------------------------------------------------------------------------------# 

# READ RASTER (needs to be an square NDVI image of just the sugarcane field)
imagepath = 'D:/Sugarcane_Project/201601_Sugar_Bacolod_sugarcanfields_zone_1/orthomosaics/output/roi.tiff'
image = cv2.imread(imagepath, 0)
try:
    if not image:
        print "NO IMAGE LOADED IN:\n", imagepath
except ValueError:
    pass

image = img_as_ubyte(roi)

def adaptivethresh_and_hough(image):
    '''
    image = square NDVI image of field
    
    Applies adaptive thresholding and cleaning on NDVI image to prepare for hough line transformation
    ''' 

    # apply adaptive opencv threshold
    th3 = cv2.adaptiveThreshold(image,255,cv2.ADAPTIVE_THRESH_MEAN_C,\
                cv2.THRESH_BINARY,95,2)

    # Plot NDVI and adaptive thresholding            
    plt.subplot(1,2,1),plt.imshow(image,cmap='nipy_spectral')
    plt.title('NDVI Image'), plt.xticks([]), plt.yticks([])
    plt.subplot(1,2,2),plt.imshow(th3,cmap = 'gray')
    plt.title('Adaptive thresholding'), plt.xticks([]), plt.yticks([])

    # Remove noise bij opening
    kernel = np.ones((3,3),np.uint8)
    opening = cv2.morphologyEx(th3, cv2.MORPH_OPEN, kernel)

    # close image using scikit
    #bw = closing(th3 > opening, square(14))

    # remove artifacts connected to image border
    cleared = opening.copy()
    clear_border(cleared)

    # Plot before and after cleaning
    plt.subplot(1,2,1),plt.imshow(th3,cmap='gray')
    plt.title('Before cleaning'), plt.xticks([]), plt.yticks([])
    plt.subplot(1,2,2),plt.imshow(cleared,cmap = 'gray')
    plt.title('After cleaning'), plt.xticks([]), plt.yticks([])

    # label image regions
    crop_regions = label(cleared)
    image_label_overlay = label2rgb(crop_regions, image=image)

    # Remove noise by area
    region_size = np.bincount(crop_regions.ravel())
    region_mask = region_size > 200
    region_mask[0] = 0
    regions_cleaned = region_mask[crop_regions]

    # Re-label cleaned areas
    crop_regions_relab = label(regions_cleaned)
    regions = regionprops(crop_regions_relab)
    region_size_new = np.bincount(crop_regions_relab.ravel())

    # Plot before and after labeling and cleaning on area
    plt.subplot(1,2,1),plt.imshow(cleared,cmap='gray')
    plt.title('Before labeling'), plt.xticks([]), plt.yticks([])
    plt.subplot(1,2,2),plt.imshow(crop_regions_relab,cmap = 'nipy_spectral')
    plt.title('After labeling'), plt.xticks([]), plt.yticks([])
    
    # Convert image to CV image
    cv_image = img_as_ubyte(regions_cleaned)

    # Apply hough line transformation
    lines = cv2.HoughLinesP(cv_image,rho=1,theta=np.pi/180,threshold=200,lines=np.array([]),
                        minLineLength=100,maxLineGap=5) # TO DO: MAKE SURE ONLY 180 RANGE IS RETURNED and minlinelength automatic adjust
    for line in lines:
        x1,y1,x2,y2 = line[0]
        cv2.line(cv_image,(x1,y1),(x2,y2),(50,255,10),2)
     
    # Extract only the coordinates from NP array
    coordinates = lines[0:,0,]
    np.save('coordinates.npy', coordinates)
    np.savetxt('coordinates.txt', coordinates)
    
    return coordinates

coordinates = adaptivethresh_and_hough(image)

#row_centers = ndimage.measurements.center_of_mass(regions_cleaned,crop_regions_relab,itemfreq(crop_regions_relab)[1:,0])
#np.savetxt('centroidscord.txt', row_centers)

#------------------------------------------------------------------------------#
## STEP 5 USE HOUGH LINES TRANSFORMATION TO ESTIMATE CROP ROW ANGLE
#------------------------------------------------------------------------------# 


def anglecalc(hough_lines):
    '''
    hough_lines = output of OpenCV hough lines transformation
    
    Calculates angle based on middle percentile of hough lines
    '''
    
    # Calculate angle
    dx = hough_lines[0:,2] - hough_lines[0:,0]
    dy = hough_lines[0:,3] - hough_lines[0:,1]
    
    rads = np.arctan2(dy,dx)
    
    # Select rads based on quantile (80% of data) and avarage
    q1, q2 = np.percentile(rads, 10), np.percentile(rads, 90)
    if np.mean(rads) > 0: # CAN GIVE ERRORS BECAUSE RADS IS IN MINUS AS OF NOW
        radssubq1, radssubq2 = np.mean((rads[rads[0:,] >= q1])), np.mean((rads[rads[0:,] <= q2]))
    else:
        radssubq1, radssubq2 = np.mean((rads[rads[0:,] <= q1])), np.mean((rads[rads[0:,] >= q2]))
    radm = (radssubq1 + radssubq2) / 2
    
    # Convert rads 
    radconv = np.pi - (-radm)
    radconv = np.pi - radconv
    
    angle = (radconv * 180) / (np.pi)
    
    radangle = (radconv, angle)
    
    return radangle

angles = anglecalc(coordinates)
#------------------------------------------------------------------------------#
## STEP 6 CREATE A GRID BASED ON ANGLE
#------------------------------------------------------------------------------# 

#imgun = cv2.imread('ndvi_ren.tif',cv2.IMREAD_UNCHANGED)
#img = cv2.imread('ndvi_ren.tif',cv2.IMREAD_GRAYSCALE)

#cv2.imshow('ndvi',img)
#cv2.waitKey(0)
#cv2.destroyAllWindows()

#pixels = img[57, 23:87]

line_amount = 6000 # must be able to divided by 40!
line_length = 5000
angle = anglecalc(coordinates)[0]
spacing = 1 # change this to automaticcaly detect cell size

def gridcreate(line_amount, line_length, extent, angle, spacing):
    '''
    line_amount = amount of lines to be created
    line_length = total length from starting point for each line
    angle = angle as determined by hough transform
    spacing = spacing between lines (pixel size of raster)
    
    Creates grid without spatial reference
    '''
    
    # Initialize array of 2 by 3
    w, h = 2, line_amount
    begin_coord = [[0 for x in range(w)] for y in range(h)]
    
    begin_coord_array = np.zeros((h, w))
    end_coord_array = np.zeros((h, w))
    
    # check if angle is positive or negative and adjust xmin
    if anglecalc(coordinates)[1] < 0:
        xmin = 0
        ymax = 0
        
        for i in range(line_amount):
            begin_coord_array[i][0] = (xmin + (i * spacing))
            
        for i in range(line_amount):
            begin_coord_array[i][1] = ymax
        
    else:
        xmin = extent.shape[1]
        ymax = extent.shape[0]
                
        for i in range(line_amount):
            begin_coord_array[i][0] = (xmin - (i * spacing))
            
        for i in range(line_amount):
            begin_coord_array[i][1] = ymax
    
    # Calulcate end points based on length and angle 
    end_coord = [[0 for x in range(w)] for y in range(h)]
    
    for i in range(line_amount):
        end_coord_array[i][0] = begin_coord_array[i][0] + line_length * math.cos(angle)
        
    for i in range(line_amount):
        end_coord_array[i][1] = begin_coord_array[i][1] - line_length * math.sin(angle)
        
    return begin_coord_array, end_coord_array
    
grid = gridcreate(line_amount, line_length, ndvi, angle, spacing)
grid = np.array(grid)

# Plot the grid
plt.imshow(ndvi)
plt.plot([grid[0,:,0], grid[1,:,0]], [grid[0,:,1], grid[1,:,1]], 'ro-')

plt.show()

gridstack = np.hstack((begin_coord_array, end_coord_array))
        
#------------------------------------------------------------------------------#
## STEP 6 EXTRACT VALUES FROM LINES
#------------------------------------------------------------------------------# 
def extractvalues(line_grid, raster):
    '''
    line_grid = grid of lines as created with the gridcreate tool
    raster = raster image of which values need to be extracted
    
    extracts values based on coordinates as created with the gridcreate tool
    Samples at each pixel using nearest neighbour and line length
    '''
    
    # Extract all coordinates for later plotting
    x0, y0 = line_grid[0,:,0], line_grid[0,:,1]
    x1, y1 = line_grid[1,:,0], line_grid[1,:,1]
    
    length = (int(np.hypot(x1[1]-x0[1], y1[1]-y0[1])))
    
    # Check if line length is correct
    if length != line_length:
        print 'Error: line length not in same range'
        return
        
    else:
        w, h = line_length, line_amount
        xvalues_array = np.zeros((h, w))
        for i in range(line_amount):
            xvalues_array[i] = np.linspace(line_grid[0][i][0], line_grid[1][i][0], line_length)
    
        yvalues_array = np.zeros((h, w))
        for i in range(line_amount):
            yvalues_array[i] = np.linspace(line_grid[0][i][1], line_grid[1][i][1], line_length)
            
        # Transpose image and add 0 values to increase extent according to grid
        imaget = np.transpose(raster)
        # determine extension of image canvas based on grid
        if anglecalc(coordinates)[1] < 0:
            plus = (((xvalues_array[(line_amount-1),(line_length-1)]) + 2) - ndvi.shape[0])
        else:
            plus = (((xvalues_array[(0),(0)]) + 2) + ndvi.shape[0])
        newcol0 = np.zeros((plus, ndvi.shape[0]))
        imaget2 = np.append(imaget, newcol0, axis = 0)
        newcol1 = np.zeros((imaget2.shape[0], plus))
        imaget2 = np.append(imaget2, newcol1, axis = 1)
        
        # Calculate values underneath lines
        try:
            zi = imaget2[xvalues_array.astype(np.int), yvalues_array.astype(np.int)]
        except IndexError:
            zi = 'null'
            
        zitrans = zi.T
    
    return zitrans

 # sets -99 values to nan
#rasteri = ndvi.astype('int')

# Calculate the values underneath the lines
rasterf = ndvi.astype('float')
values = extractvalues(grid, rasterf) # Non integer warning since changing ndvi calc no data to float('NaN'), therefore convert to float


# Get pixel extent of most outer line:
#xmax = xvalues_array[(line_amount-1),(line_length-1)]
#ymin = yvalues_array[0,0]

# create column column
#new_colx = np.zeros(((int(xmax) + 2 - ndvi.shape[0]), ndvi.shape[1]))


# Extract line from coordinates
x0, y0 = grid[0,:,0], grid[0,:,1]
x1, y1 = grid[1,:,0], grid[1,:,1]

#-- Plot grid and value lines
fig, axes = plt.subplots(nrows=2)
axes[0].imshow(raster)
axes[0].plot([x0, x1], [y0, y1], 'ro-')
axes[0].axis('image')

axes[1].plot(values)

plt.show()

#------------------------------------------------------------------------------#
## STEP 7 GET LINES WITH HIGHEST VALUE
#------------------------------------------------------------------------------# 

# Calculate lines with highest values
maskedvalues = np.ma.masked_array(values, np.isnan(values)) # mask NaN values
valuestotal = np.zeros((line_amount, 1))
for i in range(line_amount):
    valuestotal[i] = np.sum(maskedvalues[:,i])

#valuestotal_max = valuestotal > np.percentile(valuestotal, 75)
lineam = (line_amount / 30) # interval of 40 lines (differs per crop spacing! -> look for automatic detection)
valuestotal_interval = ((valuestotal.reshape(lineam, 30).T) > np.percentile((valuestotal.reshape(lineam, 30).T), 
                         99, axis=0)).flatten('F')
                         


# Select lines from top 25 percentile for rotating
best_lines = [val for is_good, val in zip(valuestotal_interval, (valuestotal.tolist())) if is_good]

# Create array with sorted begin and end coordinates
x0, y0 = grid[0,:,0], grid[0,:,1]
begxy = np.vstack((x0, y0)).T
x1, y1 = grid[1,:,0], grid[1,:,1]
endxy = np.vstack((x1, y1)).T
allresults = np.array(np.hstack((begxy, endxy)), dtype=float)

# Get the best lines with coordinates
bestlines_coord = allresults[valuestotal_interval == True]
bestlines_coord = bestlines_coord.astype(np.float)

# Plot the grid with best lines
plt.imshow(ndvi)
plt.plot([bestlines_coord[:,0], bestlines_coord[:,2]], [bestlines_coord[:,1], bestlines_coord[:,3]], '-ro')

plt.show()

#------------------------------------------------------------------------------#
## STEP 8 ROTATE LINES AROUND CENTROID
#------------------------------------------------------------------------------# 

def rotate_lines(bestlines_coord, degreelist):
    '''
    line_grid_max = grid of max value lines as created with the gridcreate tool and subsetted
    degrees = list of degrees to rotate lines with
    
    Rotate lines the given angle about their centers. 
    '''
    # Convert angle from degrees to radians and create list
    #degreelist = degreelist
    theta = []
    for degree in degreelist:
        theta.append(math.radians(degree))       
    cosang = []
    sinang = []
    for theta in theta:
        cosang.append(math.cos(theta))
        sinang.append(math.sin(theta))
    
    cosin = np.vstack((cosang, sinang)).T
    
    # Find logical center (avg x and avg y) of entire line 
    cxcy = np.zeros((bestlines_coord.shape[0], 2))
    for i in range(len(bestlines_coord)):
        cxcy[i,0] = (bestlines_coord[i,0] + bestlines_coord[i,2]) / 2
        cxcy[i,1] = (bestlines_coord[i,1] + bestlines_coord[i,3]) / 2
    
    # Rotate each around whole polyline's center point
    tx1, ty1 = bestlines_coord[0,0] - cxcy[0,0], bestlines_coord[0,1] - cxcy[0,1]   
    tx2, ty2 = bestlines_coord[0,2] - cxcy[0,0], bestlines_coord[0,3] - cxcy[0,1] 

    r1x = [] 
    r1y = [] 
    r2x = [] 
    r2y = []
    
    p1x = []
    p1y = [] 
    p2x = [] 
    p2y = []

    # Create 5 times rotation for each point
    for j in range(len(degreelist)):
        r1x.append(tx1 * cosin[j,0] + ty1 * cosin[j,1])
        r1y.append(-tx1 * cosin[j,1] + ty1 * cosin[j,0])
        r2x.append(tx2 * cosin[j,0] + ty2 * cosin[j,1])
        r2y.append(-tx2 * cosin[j,1] + ty2 * cosin[j,0])
    
    # add each rotated values to cxcy list
    for i in xrange(len(bestlines_coord)):
        for j in range(len(degreelist)):
            p1x.append(r1x[j] + cxcy[i,0])   
            p1y.append(r1y[j] + cxcy[i,1])
            p2x.append(r2x[j] + cxcy[i,0])
            p2y.append(r2y[j] + cxcy[i,1])  
        
        rotatedlines = np.vstack((p1x, p1y, p2x, p2y)).T
    
    return(rotatedlines)

degreelist = [359.8, 359.9, 0, 0.1, 0.2]
rotatedlines = rotate_lines(bestlines_coord, degreelist)

# Plot the rotated lines
plt.imshow(ndvi)
plt.plot([rotatedlines[:,0], rotatedlines[:,2]], [rotatedlines[:,1], rotatedlines[:,3]], '-ro')

plt.show()

#-------------------------------------------------------------------#
# Get values from rotated lines
raster = ndvi

w, h = line_length, len(rotatedlines)
xvalues_array = np.zeros((h, w))
for i in range(len(rotatedlines)):
    xvalues_array[i] = np.linspace(rotatedlines[i,0], rotatedlines[i,2], line_length)
    
yvalues_array = np.zeros((h, w))
for i in range(len(rotatedlines)):
    yvalues_array[i] = np.linspace(rotatedlines[i,1], rotatedlines[i,3], line_length)
            
# Transpose image and add 0 values to increase extent according to grid
imaget = np.transpose(raster)
# determine extension of image canvas based on grid
if anglecalc(coordinates)[1] < 0:
    plus = (((xvalues_array[(line_amount-1),(line_length-1)]) + 2) - ndvi.shape[0])
else:
    plus = (((xvalues_array[(0),(0)]) + 2) + ndvi.shape[0])
newcol0 = np.zeros((plus, ndvi.shape[0]))
imaget2 = np.append(imaget, newcol0, axis = 0)
newcol1 = np.zeros((imaget2.shape[0], plus))
imaget2 = np.append(imaget2, newcol1, axis = 1)
        
# Calculate values underneath lines
try:
    zi = imaget2[xvalues_array.astype(np.int), yvalues_array.astype(np.int)]
except IndexError:
    zi = 'null'
            
valuesrotated = zi.T

# Get line with highest rotated values

# Calculate lines with highest values
maskedvalues = np.ma.masked_array(valuesrotated, np.isnan(valuesrotated)) # mask NaN values
valuestotalro = np.zeros((len(rotatedlines), 1))
for i in range(len(rotatedlines)):
    valuestotalro[i] = np.sum(maskedvalues[:,i])


lineam = (len(rotatedlines) / 5) # interval of 40 lines (differs per crop spacing! -> look for automatic detection)
valuestotal_intervalrotated = ((valuestotalro.reshape(lineam, 5).T) > np.percentile((valuestotalro.reshape(lineam, 5).T), 
                         99, axis=0)).flatten('F')
                                                
# Select lines from top 25 percentile for rotating
resultlines = rotatedlines[valuestotal_intervalrotated == True]

# Plot the grid with best lines
plt.imshow(ndvi)
plt.plot([resultlines[:,0], resultlines[:,2]], [resultlines[:,1], resultlines[:,3]], '-ro')

plt.show()



#------------------------------------------------------------------------------#
## STEP 9 TRANSLATE IMAGE COORDINATES TO REAL WORLD COORDINATES
#------------------------------------------------------------------------------# 

# Flip Y axis
projected_lines = np.copy(gridstack)
projected_lines[:,1] = abs(projected_lines[:,1] - ndvi.shape[0])
projected_lines[:,3] = abs(projected_lines[:,3] - ndvi.shape[0])

# Get origin & cellsize of input image
proj_im = gdal.Open(filename)
geoinformation = proj_im.GetGeoTransform()
xres = abs(geoinformation[1])
yres = abs(geoinformation[5])
xorigin= abs(geoinformation[0])
yorigin= abs(geoinformation[3]) - (ndvi.shape[0] * yres)

# Correct coordinates
projected_lines[:,0] = projected_lines[:,0]*xres+xorigin
projected_lines[:,1] = projected_lines[:,1]*yres+yorigin

projected_lines[:,2] = projected_lines[:,2]*xres+xorigin
projected_lines[:,3] = projected_lines[:,3]*yres+yorigin

# Save projected points
#name='projectedlines.csv'
#np.savetxt(name, projected_lines ,delimiter=',',header='x0,y0,x1,y1', comments='')

np.savetxt('rotatedlines.txt', projected_lines)

# Create spatial lines in R
import subprocess

# Define command and arguments
command = 'Rscript'
path2script = 'C:\Users\darellvdv\Documents\croprowdetection\create_lines_from_python.R'

# Build subprocess command
cmd = [command, path2script]

# check_output will run the command and store to result
x = subprocess.check_output(cmd, universal_newlines=True)

