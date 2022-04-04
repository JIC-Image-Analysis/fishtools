from ij import IJ
from ij.gui import PointRoi, PolygonRoi, Roi
from ij.plugin.frame import RoiManager
import math
import csv
from ij.io import DirectoryChooser

# Get current ImagePlus
image = IJ.getImage()
RM = RoiManager()
rm = RM.getRoiManager()

dc_in = DirectoryChooser("Choose directory with annotation file")
path_in = dc_in.getDirectory()

imageTitle=image.getTitle()
fname = 'coordinates'+imageTitle[11:14]+'.csv'

locations=csv.reader(open(path_in+fname))
next(locations)
# Create ROI
roi = PointRoi()
for row in locations:
    roi.setCounter(int(row[0]))
    roi.addPoint(image,int(row[1]),int(row[2]))

rm.add(image, roi, 0)
rm.runCommand("Associate", "true")	 
rm.runCommand("Show All with labels")
