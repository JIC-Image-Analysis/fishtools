import numpy as np
import skimage.segmentation
from PIL import ImageDraw, ImageFont
import PIL.Image as pilImage

from dtoolbioimage import Image as dbiImage, scale_to_uint8

from fishtools.utils import force_to_2d_rgb


def maxproj_composite_from_fishimage(fishimage, probe_channel=0):
     probe_maxproj = np.max(fishimage.probes[probe_channel], axis=2)
     nuclear_maxproj = np.max(fishimage.nuclei, axis=2)

     blank = np.zeros_like(probe_maxproj)

     return np.dstack([probe_maxproj, blank, nuclear_maxproj])


def merge_images(im1, im2, p=0.5):
    return p * force_to_2d_rgb(im1) + (1-p) * force_to_2d_rgb(im2)


def merge_composite_from_fishimage(fishimage, cell_mask):
     composite = maxproj_composite_from_fishimage(fishimage)
     merge_composite = merge_images(cell_mask, composite, p=0.1).view(dbiImage)

     return merge_composite


def visualise_counts(maxproj, scaled_cell_regions, centroids):
    maxproj = scale_to_uint8(maxproj)
    
    centroids_by_cell = {
        idx: {tuple(p) for p in scaled_cell_regions.rprops[idx].coords} & set(centroids)
        for idx in scaled_cell_regions.labels
    }

    counts_by_cell = {
        idx: len(centroids)
        for idx, centroids in centroids_by_cell.items()
    }

    boundary_image = scale_to_uint8(
         skimage.segmentation.mark_boundaries(maxproj, scaled_cell_regions)
     )
    boundary_pilimage = pilImage.fromarray(boundary_image)
    draw = ImageDraw.ImageDraw(boundary_pilimage)
    font = ImageFont.truetype("Microsoft Sans Serif.ttf", size=16)

    for idx, count in counts_by_cell.items():
        r, c = map(int, scaled_cell_regions.rprops[idx].centroid)
        draw.text((c-20, r-12), f"[{idx}] {count}", font=font)

    return boundary_pilimage