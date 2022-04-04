from ij import IJ, ImagePlus
from ij.gui import PointRoi, PolygonRoi, Roi
from ij.plugin.frame import RoiManager
import math
import csv
from ij.io import DirectoryChooser

# Get current ImagePlus
image = IJ.getImage()
RM = RoiManager()
rm = RM.getRoiManager()

dc_in = DirectoryChooser("Choose directory with locations.csv file")
path_in = dc_in.getDirectory()

imageTitle=image.getTitle()
fname = 'locations'+imageTitle[11:14]+'.csv'

locations=csv.reader(open(path_in+fname))
next(locations)

counters=[]
zslices=[]
# Create ROI
roi = PointRoi()
for row in locations:
    roi.addPoint(image,int(row[2]),int(row[1]))
    zslices.append(int(row[3]))
    if int(row[0])<100:
        counters.append(int(row[0]))
    else:
        counters.append(0)
    
roi.setCounters([256*((zslice+1)*4-2)+counter for zslice,counter in zip(zslices,counters)])
    
rm.addRoi(roi)
    
rm.runCommand("Associate", "true")	 
rm.runCommand("Show All with labels")
