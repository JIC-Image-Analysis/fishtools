import numpy as np
from dtoolbioimage import scale_to_uint8


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