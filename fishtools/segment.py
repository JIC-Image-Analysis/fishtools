import numpy as np

import skimage.filters
import skimage.exposure


def cell_mask_from_fishimage(fishimage, params, probe_channel=0):

    ks = params.ks
    bs = params.bs
    sigma = params.sigma

    minproj = np.min(fishimage.probes[probe_channel], axis=2)
    eq = skimage.exposure.equalize_adapthist(minproj, kernel_size=(ks, ks))
    eq_nuclear_proj_smoothed = skimage.filters.gaussian(eq, sigma=sigma)
    thresh_image = skimage.filters.threshold_local(eq_nuclear_proj_smoothed, block_size=bs)
    result = (eq_nuclear_proj_smoothed > thresh_image)

    return result