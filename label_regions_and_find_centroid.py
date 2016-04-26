# -*- coding: utf-8 -*-
"""
Created on Fri Apr 22 16:30:13 2016

@author: darellvdv
"""

import os
import gdal
import math
import cv2
from matplotlib import pyplot as plt
import matplotlib.patches as mpatches
#from cv2 import *
import numpy as np
from math import atan2, pi

from skimage import data
from skimage.filters import threshold_otsu, threshold_adaptive
from skimage.filters import threshold_otsu
from skimage.segmentation import clear_border
from skimage.measure import label
from skimage.morphology import closing, square
from skimage.measure import regionprops
from skimage.color import label2rgb
from skimage.draw import ellipse
from skimage.measure import label, regionprops
from skimage.transform import rotate

from scipy import ndimage
from scipy.stats import itemfreq

# Set working directory
path = os.chdir('D:/Sugarcane_Project/201601_Sugar_Bacolod_sugarcanfields_zone_1/orthomosaics/output')
print os.getcwd(); # Prints the working directory

# READ RASTER (needs to be an NDVI image)
image = cv2.imread('ndvi_ren.tif',0)

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

cv2.imwrite('crop_row.tif', cleared)

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

#centroids = props.centroid(regions)

image_label_overlay = label2rgb(crop_regions_relab, image=image)


row_centers = ndimage.measurements.center_of_mass(regions_cleaned,crop_regions_relab,itemfreq(crop_regions_relab)[1:,0])

cv2.imwrite('crop_regions_relab.tif', crop_regions_relab)


np.savetxt('centroidscord.txt', row_centers)

fig, ax = plt.subplots()
ax.imshow(image, cmap=plt.cm.gray)

for props in regions:
    y0, x0 = props.centroid
    orientation = props.orientation
    x1 = x0 + math.cos(orientation) * 0.5 * props.major_axis_length
    y1 = y0 - math.sin(orientation) * 0.5 * props.major_axis_length
    x2 = x0 - math.sin(orientation) * 0.5 * props.minor_axis_length
    y2 = y0 - math.cos(orientation) * 0.5 * props.minor_axis_length

    ax.plot((x0, x1), (y0, y1), '-r', linewidth=2.5)
    ax.plot((x0, x2), (y0, y2), '-r', linewidth=2.5)
    ax.plot(x0, y0, '.g', markersize=15)

    minr, minc, maxr, maxc = props.bbox
    bx = (minc, maxc, maxc, minc, minc)
    by = (minr, minr, maxr, maxr, minr)
    ax.plot(bx, by, '-b', linewidth=2.5)

ax.axis((0, 600, 600, 0))
plt.show()


#-----------------------------------------
# Calulcate lines from skikit regions cleaned

from skimage import img_as_ubyte

cv_image = img_as_ubyte(regions_cleaned)

lines = cv2.HoughLinesP(cv_image,rho=0.9,theta=np.pi/180,threshold=500,lines=np.array([]),
                        minLineLength=100,maxLineGap=25)
for line in lines:
   x1,y1,x2,y2 = line[0]
   cv2.line(cv_image,(x1,y1),(x2,y2),(50,255,10),2)
   
cv2.imwrite('croprows_regionscleaned.tif',cv_image)
# Extract only the coordinates from NP array
coordinates = lines[0:,0,]
np.save('coordinates.npy', coordinates)
np.savetxt('coordinates.txt', coordinates)

cv_image = img_as_ubyte(regions_cleaned)

lines = cv2.HoughLinesP(cv_image,rho=0.9,theta=np.pi/34,threshold=500,lines=np.array([]),
                        minLineLength=100,maxLineGap=10)
for line in lines:
   x1,y1,x2,y2 = line[0]
   cv2.line(cv_image,(x1,y1),(x2,y2),(50,255,10),2)
   
cv2.imwrite('croprows_regionscleaned2.tif',cv_image)

import sys
sys.exit()