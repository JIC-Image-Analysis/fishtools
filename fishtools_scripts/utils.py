import numpy as np
from dtoolbioimage import scale_to_uint8

import logging


logger = logging.getLogger("fishtools")


def clipped_image_difference_uint8(im1, im2):
    return np.clip(im1.astype(np.int16) - im2, 0, 255).astype(np.uint8)


def extract_nuclei(annotated_im):

    ar1 = annotated_im[:,:,1]
    ar2 = np.maximum(annotated_im[:,:,0], annotated_im[:,:,2])

    result = clipped_image_difference_uint8(ar1, ar2)
    
    return scale_to_uint8(result > 30)


def force_to_2d_rgb(im):
    if len(im.shape) == 2:
        return scale_to_uint8(np.dstack(3 * [im]))
    
    rdim, cdim, zdim = im.shape
    if zdim == 3:
        return scale_to_uint8(im)

    raise ValueError("Can't handle that image type")


def crop_to_non_empty(im):
    """
    Sample pixel at 30, 30 to determine background colour, then return an
    image cropped to a bounding box based on colour other than that background.
    """

    notmask = tuple(im[30,30,:])
    rr, cc, zz = np.where(im != notmask)
    r, h = min(rr), max(rr) - min(rr)
    c, w = min(cc), max(cc) - min(cc)

    logging.debug(f"crop to {r},{c}, {w},{h}")

    return im[r:r+h,c:c+w]


def select_near_colour(im, colour, tolerance=1):
    """
    Given an RGB image, return a mask where values are True if the pixel value
    is within tolerance of the given colour.
    """

    within_tolerance = np.abs(im - colour) < tolerance

    return within_tolerance.all(axis=2)