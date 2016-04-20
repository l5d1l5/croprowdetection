# -*- coding: utf-8 -*-
"""
Created on Thu Apr 14 13:05:18 2016

@author: darellvdv
"""

import os
import cv2
import numpy as np
from math import atan2, pi

# Set working directory
path = os.chdir('D:/Sugarcane_Project/201601_Sugar_Bacolod_sugarcanfields_zone_1/orthomosaics/output')
print os.getcwd(); # Prints the working directory

# HOUGH LINES P
img = cv2.imread('vegNArenderd.tif')
gray = cv2.cvtColor(img,cv2.COLOR_BGR2GRAY)
lines = cv2.HoughLinesP(gray,rho=0.9,theta=np.pi/180,threshold=500,lines=np.array([]),
                        minLineLength=800,maxLineGap=25)
for line in lines:
   x1,y1,x2,y2 = line[0]
   cv2.line(img,(x1,y1),(x2,y2),(0,255,10),2)
   
cv2.imwrite('croprowsP2.tif',img)
# Extract only the coordinates from NP array
coordinates = lines[0:,0,]
np.save('coordinates.npy', coordinates)
np.savetxt('coordinates.txt', coordinates)

# Calculate angle
#linepoints = lines[0:1653]

#for i in linepoints:
   # rads = atan2(i[0,1] - i[0,3], i[0,0] - i[0,2])
   # rads %= 2*pi
   # degr = 180 * rads / pi

#degrconv = 180 - degr # supplementary angle

# Read geoTIFFS
#from osgeo.gdal_array import *