import json

import numpy as np

import skimage.draw
import skimage.filters
import skimage.exposure
import skimage.segmentation
import scipy.ndimage

from dtoolbioimage.segment import Segmentation
from dtoolbioimage import Image as dbiImage


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


def nuc_cell_mask_from_fishimage(fishimage, params):

    ks = params.ks
    bs = params.bs
    sigma = params.sigma

    minproj = np.min(fishimage.nuclei, axis=2)
    eq = skimage.exposure.equalize_adapthist(minproj, kernel_size=(ks, ks))
    eq_nuclear_proj_smoothed = skimage.filters.gaussian(eq, sigma=sigma)
    thresh_image = skimage.filters.threshold_local(eq_nuclear_proj_smoothed, block_size=bs)
    result = (eq_nuclear_proj_smoothed > thresh_image)

    return result


def label_image_from_coords(label_coords, dim):
    label_img = np.zeros(dim, dtype=np.uint16)

    current_label = 1

    for label, points in label_coords.items():
        for p in points:
            r, c = p
            rr, cc = skimage.draw.disk((r, c), 12)
            try:
                label_img[rr, cc] = current_label
                current_label += 1
            except IndexError:
                pass
    
    return label_img


def label_coords_from_points_fpath(fpath):
    with open(fpath) as fh:
        label_coords = json.load(fh)

    return label_coords


def label_image_from_points_fpath(fpath, dim):
    label_coords = label_coords_from_points_fpath(fpath)
    return label_image_from_coords(label_coords, dim)


def filter_segmentation_by_region_list(segmentation, region_ids):
    rids_not_in_files = segmentation.labels - set(region_ids)
    trimmed_segmentation = segmentation.copy()

    for rid in rids_not_in_files:
        trimmed_segmentation[np.where(trimmed_segmentation == rid)] = 0

    return Segmentation.from_array(trimmed_segmentation)


def segmentation_from_nuclear_channel_and_markers(fishimage, label_img, params):

    nucmask = nuc_cell_mask_from_fishimage(fishimage, params)
    # print(nucmask.shape, label_img.shape)
    assert nucmask.shape == label_img.shape

    n_segmentation = skimage.segmentation.watershed(
        -scipy.ndimage.distance_transform_edt(nucmask),
        markers=label_img,
        mask=nucmask
    ).view(Segmentation)

    return n_segmentation


def segmentation_from_cellmask_and_label_image(cell_mask, label_img):
    noholes = skimage.morphology.remove_small_holes(cell_mask, area_threshold=150)

    segmentation = skimage.segmentation.watershed(
        -scipy.ndimage.distance_transform_edt(noholes),
        markers=label_img,
        mask=noholes
    )

    return Segmentation.from_array(segmentation)


def filter_segmentation_by_label_coords(segmentation, label_coords):

    valid_labels = {
        segmentation[tuple(p)]
        for p in label_coords["1"]
    }
    
    return filter_segmentation_by_region_list(segmentation, valid_labels)


def scale_segmentation(cell_regions, maxproj):
    scaled_cell_regions = Segmentation.from_array(
        skimage.transform.resize(
            cell_regions,
            maxproj.shape,
            anti_aliasing=False,
            order=0,
            preserve_range=True
        ).astype(int)
    )

    return scaled_cell_regions


def get_filtered_segmentation(dataitem, params):
    nuc_label_image = segmentation_from_nuclear_channel_and_markers(
        dataitem.fishimage,
        skimage.measure.label(dataitem.scaled_markers),
        params
    )
    nuc_label_image.pretty_color_image.view(dbiImage).save("nuc_label_img.png")

    segmentation = segmentation_from_cellmask_and_label_image(
        dataitem.cell_mask(params),
        nuc_label_image
    )

    scaled_good_mask = scale_segmentation(dataitem.good_mask, dataitem.maxproj)
    labelled_points = skimage.measure.label(scaled_good_mask)
    rprops = skimage.measure.regionprops(labelled_points)
    region_centroids = [r.centroid for r in rprops]
    icentroids = [(int(r), int(c)) for r, c in region_centroids]
    good_regions = [segmentation[r, c] for r, c in icentroids]

    filtered_segmentation = filter_segmentation_by_region_list(
        segmentation,
        good_regions
    )

    return filtered_segmentation
