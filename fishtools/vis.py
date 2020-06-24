import numpy as np

from fishtools.utils import force_to_2d_rgb


def maxproj_composite_from_fishimage(fishimage, probe_channel=0):
     probe_maxproj = np.max(fishimage.probes[probe_channel], axis=2)
     nuclear_maxproj = np.max(fishimage.nuclei, axis=2)

     blank = np.zeros_like(probe_maxproj)

     return np.dstack([probe_maxproj, blank, nuclear_maxproj])



   