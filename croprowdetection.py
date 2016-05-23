# -*- coding: utf-8 -*-
"""
Created on Mon May 23 10:31:57 2016

Crop row detection script based on sugarcane fields

@author: darellvdv
"""
### CROP ROW DETECTION ###

# Load modules
import numpy as np
import os
import sys
import errno
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
import argparse

import warnings
warnings.filterwarnings("ignore")

# construct the argument parser and parse the arguments
ap = argparse.ArgumentParser()
ap.add_argument("-i", "--image", required=True, help="Path to the image")
ap.add_argument("-c", "--config", required=True, help="Path to the config file")
ap.add_argument("-o", "--output", required=True, help="Path to the output folder")
ap.add_argument("-p", "--plot", required=True, help="Plotting yes or no")
args = vars(ap.parse_args())

# Load settings file
config = args["config"]
import config

# Load parameters
line_amount = config.line_amount # must be able to divided by 40!
line_length = config.line_length
spacing = config.spacing # change this to automaticcaly detect cell size
degreelist = config.degreelist #[359.8, 359.9, 0, 0.1, 0.2] # degrees to rotate best lines

# Get filename and create output dir
filename = args["image"]
output = args["output"]

def make_sure_path_exists(path):
    try:
        os.makedirs(path)
    except OSError as exception:
        if exception.errno != errno.EEXIST:
            raise
make_sure_path_exists(output)
output = os.chdir(output)

# Load functions
def click_and_crop(event, x, y, flags, param):
    '''
    Function to click and crop a square of the ndvi raster
    '''
    
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
         cv2.rectangle(image, refPt[0], refPt[1], (255, 0, 0), 2)
         cv2.imshow("image", image)

        
def adaptivethresh_and_hough(image):
    '''
    image = square NDVI image of field
    
    Applies adaptive thresholding and cleaning on NDVI image to prepare for hough line transformation
    ''' 

    # apply adaptive opencv threshold
    th3 = cv2.adaptiveThreshold(image,255,cv2.ADAPTIVE_THRESH_MEAN_C,\
                cv2.THRESH_BINARY,95,2)

    # Remove noise bij opening
    kernel = np.ones((3,3),np.uint8)
    opening = cv2.morphologyEx(th3, cv2.MORPH_OPEN, kernel)

    # close image using scikit
    #bw = closing(th3 > opening, square(14))

    # remove artifacts connected to image border
    cleared = opening.copy()
    clear_border(cleared)

    # label image regions
    crop_regions = label(cleared)

    # Remove noise by area
    region_size = np.bincount(crop_regions.ravel())
    region_mask = region_size > 200
    region_mask[0] = 0
    regions_cleaned = region_mask[crop_regions]

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
    cv2.imwrite('detected_lines.png', cv_image)
    
    return coordinates
    
    
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
    q1, q2 = np.percentile(rads, 20), np.percentile(rads, 80)
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
    if length != (config.line_length):
        print 'Warning: line length not in same range'
        pass
        
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
            plus = (((xvalues_array[(line_amount-1),(line_length-1)]) + 2) - int(ndvi.shape[0]))
        else:
            plus = (((xvalues_array[(0),(0)]) + 2) + int(ndvi.shape[0]))
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
#------------------------------------------------------------------------------#
## STEP 1 LOAD CROP ORTHOPHOTO
#------------------------------------------------------------------------------#

print 'STEP 1: LOADING IMAGE...'

# Open Ortho and print metadata
filename = args["image"]
dataSource = gdal.Open(filename, GA_ReadOnly)
if not dataSource:
    print "THE FOLLOWING RASTER FAILED TO LOAD:\n", filename
    sys.exit()

else:
    print "IMAGE SUCESSFULLY LOADED!"
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
           
print 'STEP 2: CALCULATING NDVI...'

# Read data into an array
band1Arr = dataSource.GetRasterBand(1).ReadAsArray(0,0,dataSource.RasterXSize, dataSource.RasterYSize)
band4Arr = dataSource.GetRasterBand(4).ReadAsArray(0,0,dataSource.RasterXSize, dataSource.RasterYSize)

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
#print ndvi[ndvi>-99].min() # gives error since turning -99 to NaN!

# Write the result to disk
driver = gdal.GetDriverByName('GTiff')
outDataSet=driver.Create('ndvi_f6.tif', dataSource.RasterXSize, dataSource.RasterYSize, 1, GDT_Float32)
outBand = outDataSet.GetRasterBand(1)
outBand.WriteArray(ndvi,0,0)
outBand.SetNoDataValue(float('NaN'))

# set the projection and extent information of the dataset
outDataSet.SetProjection(dataSource.GetProjection())
outDataSet.SetGeoTransform(dataSource.GetGeoTransform())

# Flush to save
outBand.FlushCache()
outDataSet.FlushCache()

print "SUCESSFULLY CALCULATED NDVI AND SAVED TO DISC!"

#----------------------------------------------------------------------------#
## STEP 3: CONVERT TO GRAYSCALE
#----------------------------------------------------------------------------#

print "STEP 3: NOW CONVERTING NDVI TO GRAYSCALE"

maskedndvi = np.ma.masked_array(ndvi, np.isnan(ndvi)) # mask NaN value
image = scipy.misc.toimage(maskedndvi, cmin=(np.min(maskedndvi)), cmax=(np.max(maskedndvi)))
image = img_as_ubyte(image)

#------------------------------------------------------------------------------#
## STEP 4 INTERACTIVE CROPPING
#------------------------------------------------------------------------------#

print "STEP 4 INTERACTIVE CROPPING..."
print "Now crop image; draw rectangle, c = crop, r = reset, q = quit"

# initialize the list of reference points and boolean indicating
# whether cropping is being performed or not
refPt = []
cropping = False

# load the image, clone it, and setup the mouse callback function
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

print "SUCCESFULLY CROPPED IMAGE FOR HOUG LINE DETECTION!"

#------------------------------------------------------------------------------#
## STEP 5 APPLY ADAPTIVE THRESHOLDING AND CLEANING &
## APPLY HOUGH LINES TRANSFORMATION
#------------------------------------------------------------------------------# 

print "STEP 5 DETECTING HOUGH LINES..."

image = img_as_ubyte(roi)
coordinates = adaptivethresh_and_hough(image)

print 'HOUGH LINES DETECTED!'

#------------------------------------------------------------------------------#
## STEP 6 USE HOUGH LINES TRANSFORMATION TO ESTIMATE CROP ROW ANGLE
#------------------------------------------------------------------------------# 
print "STEP 6 CALCULATING ANGLE..."

angles = anglecalc(coordinates)

print 'ANGLE CALCULATED! The angles are (rads, degrees):'
print angles
#------------------------------------------------------------------------------#
## STEP 7 CREATE A GRID BASED ON ANGLE
#------------------------------------------------------------------------------# 
print "STEP 7 CREATING GRID..."

grid = gridcreate(line_amount, line_length, ndvi, (anglecalc(coordinates)[0]), spacing)
grid = np.array(grid)

# Plot the grid
#if args["plot"] == 'yes':
#    print 'GRID CREATED! ..now showing plot'
#    plt.imshow(ndvi)
#    plt.plot([grid[0,:,0], grid[1,:,0]], [grid[0,:,1], grid[1,:,1]], 'ro-')
    
#    plt.show()
#else:
#    print 'GRID CREATED!'
#    pass

#------------------------------------------------------------------------------#
## STEP 8 EXTRACT VALUES FROM LINES
#------------------------------------------------------------------------------# 
print "STEP 8 EXTRACTING VALUES FROM GRID..."

# Calculate the values underneath the lines
values = extractvalues(grid, ndvi) 

print 'VALUES EXTRACTED!'

#------------------------------------------------------------------------------#
## STEP 9 GET LINES WITH HIGHEST VALUE
#------------------------------------------------------------------------------# 
print "STEP 9 DETERMINING LINES WITH HIGHEST VALUES..."

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
if args["plot"] == 'yes':
    print 'BEST LINES DETERMINED! ..now showing plot'
    plt.imshow(ndvi)
    plt.plot([bestlines_coord[:,0], bestlines_coord[:,2]], [bestlines_coord[:,1], bestlines_coord[:,3]], '-ro')

    plt.show()

else:
    print 'BEST LINES DETERMINED!'
    pass

#------------------------------------------------------------------------------#
## STEP 10 ROTATE LINES AROUND CENTROID
#------------------------------------------------------------------------------# 
print "STEP 10 ROTATING LINES WITH HIGHEST VALUES AND EXTRACT NEW VALUES..."

rotatedlines = rotate_lines(bestlines_coord, degreelist)

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
    plus = (((xvalues_array[(line_amount-1),(line_length-1)]) + 2) - int(ndvi.shape[0]))
else:
    plus = (((xvalues_array[(0),(0)]) + 2) + int(ndvi.shape[0]))
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
if args["plot"] == 'yes':
    print 'BEST ROTATED DETERMINED! ..now showing plot'
    plt.imshow(ndvi)
    plt.plot([resultlines[:,0], resultlines[:,2]], [resultlines[:,1], resultlines[:,3]], '-ro')

    plt.show()

else:
    print 'BEST ROTATED DETERMINED!'
    pass

#------------------------------------------------------------------------------#
## STEP 11 TRANSLATE IMAGE COORDINATES TO REAL WORLD COORDINATES
#------------------------------------------------------------------------------# 
print "STEP 10 TRANSLATING IMAGE COORDINATES TO REAL WORLD COORDINATES AND OUTPUTTING SHAPEFILE USING R..."

#begincoord1 = grid[0:1,:,0].T
#begincoord2 = grid[0:1,:,1].T
#begincoord = np.hstack((begincoord1, begincoord2))
#endcoord1 = grid[0,:,0:1]
#endcoord2 = grid[1,:,0:1]
#endcoord = np.hstack((endcoord1, endcoord2))

#gridstack = np.hstack((begincoord, endcoord))

# Flip Y axis
projected_lines = np.copy(resultlines)
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
x = subprocess.call(cmd, universal_newlines=True)

print "ALL DONE!"