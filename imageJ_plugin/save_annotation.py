from ij import IJ
from ij.gui import PolygonRoi, Roi
import math
import csv
from ij.io import DirectoryChooser

# Get current ImagePlus
image = IJ.getImage()

dc = DirectoryChooser("Choose directory to save")
path = dc.getDirectory()

imageTitle=image.getTitle()
fname = 'coordinates'+imageTitle[11:14]+'.csv'

# Get current ROI
roi = image.getRoi()
if roi is not None:
  # Get ROI points
  polygon = roi.getPolygon()
  n_points = polygon.npoints
  x = polygon.xpoints
  y = polygon.ypoints

  f = open(path+fname, 'wb')
  try:
    writer = csv.writer(f)
    writer.writerow(('counter','x', 'y'))
    for i in range(n_points):
      writer.writerow((str(roi.getCounter(i)),str(x[i]), str(y[i])))
  finally:
    f.close()
