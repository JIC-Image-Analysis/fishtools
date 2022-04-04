import os
import logging
import pathlib

from dataclasses import dataclass
from typing import List
from types import SimpleNamespace

import parse
import numpy as np
import skimage.measure

import pandas as pd

from dtoolbioimage import (
    Image as dbiImage,
    Image3D,
    ImageDataSet,
    scale_to_uint8
)

from utils import extract_nuclei, crop_to_non_empty, select_near_colour
from segment import scale_segmentation, cell_mask_from_fishimage
from probes import find_probe_locations_3d


logger = logging.getLogger("fishtools")


@dataclass
class FISHImage(object):

    probes: List[Image3D]
    nuclei: Image3D

    @classmethod
    def from_ids_im_sn(cls, ids, image_name, series_name, nuclear_channel_first):
        channels = list(ids.planes_index[image_name][series_name][0].keys())
        n_channels = len(channels)
        n_probe_channels = n_channels - 1

        if nuclear_channel_first:
            nuclei = ids.get_stack(image_name, series_name, 0, 0)
            probe_channels = range(1, 1+n_probe_channels)
        else:
            nuclei = ids.get_stack(image_name, series_name, 0, n_channels-1)
            probe_channels = range(0, n_probe_channels)

        probes = []
        for n in probe_channels:
            probes.append(ids.get_stack(image_name, series_name, 0, n))


        return cls(probes, nuclei)

    @classmethod
    def from_stack_fpath(cls, fpath):
        return cls([Image3D.from_file(fpath)], None)


@dataclass
class ManualAnnotatedImage(object):
    raw_image: dbiImage

    @classmethod
    def from_fpath(cls, fpath):
        im = dbiImage.from_file(fpath)

        return cls(im)

    @property
    def marked_cell_mask(self):
        return extract_nuclei(self.raw_image)


def get_specs(config):
    diriter = pathlib.Path(config.annotation_dirpath).iterdir()
    fnameiter = (fpath.name for fpath in diriter)

    try:
        annotation_type = config.annotation_type
    except KeyError:
        annotation_type = 'csv'

    if annotation_type == 'csv':
        logger.debug(f"Matching with {config.annotation_csv_template}")
        all_specs = [
            parse.parse(config.annotation_csv_template, fname)
            for fname in fnameiter
        ]
    else:
        logger.debug(f"Matching with {config.annotation_template}")
        all_specs = [
            parse.parse(config.annotation_template, fname)
            for fname in fnameiter
        ]

    valid_specs = [
        spec.named
        for spec in all_specs
        if spec is not None
    ]

    return valid_specs


def get_slice(config, spec):
    fname = config.slice_template.format(**spec)
    fpath = os.path.join(config.annotation_dirpath, fname)
    im = dbiImage.from_file(fpath)

    regionim = (im[:, :, 3] == 255)
    r = skimage.measure.regionprops(skimage.measure.label(regionim))[0]
    rmin, cmin, rmax, cmax = r.bbox
    return np.s_[rmin:rmax, cmin:cmax]

def sliceImage_from_template_and_spec(template, config, spec, sl):
    fname = template.format(**spec)
    fpath = os.path.join(config.annotation_dirpath, fname)
    im = dbiImage.from_file(fpath)
    return im[sl]

def dotImage_from_template_and_spec(template, config, spec, sl):
    fname = template.format(**spec)
    fpath = os.path.join(config.annotation_dirpath, fname)
    im = dbiImage.from_file(fpath)
    maxflatten = np.max(im[:, :, :3], axis=2)
    return maxflatten[sl]

def mask_from_template_and_spec(template, config, spec, sl):
    fname = template.format(**spec)
    fpath = os.path.join(config.annotation_dirpath, fname)
    im = dbiImage.from_file(fpath)
    maxflatten = np.max(im[:, :, :3], axis=2)
    return (maxflatten > 0)[sl]



class MultiAnnotationDataItem(object):

    def __init__(self, config, spec, use_deconv=False):
        self.config = config
        self.ids = ImageDataSet(self.config.ids_uri)

        try:
            nuclear_channel_first = self.config.nuclear_channel_first
        except KeyError:
            nuclear_channel_first = True
        
        image_name = self.config.image_name_template.format(**spec)
        series_name = self.config.series_name_template.format(**spec)
        self.fishimage = FISHImage.from_ids_im_sn(
            self.ids, image_name, series_name, nuclear_channel_first
        )

        if use_deconv:
            fname = self.config.deconv_fname_template.format(**spec)
            fpath = os.path.join(self.config.deconv_dirpath, fname)
            self.deconv_stack = Image3D.from_file(fpath)

            self.deconv_stack = np.clip(self.deconv_stack, 0, 10000)

            if self.deconv_stack.shape != self.fishimage.nuclei.shape:
                logger.warning("Deconv stack doesn't match shape, trimming")
                rdim, cdim, zdim = self.fishimage.nuclei.shape
                self.deconv_stack = self.deconv_stack[:rdim,:cdim,:zdim]
                
            self.probe_stack = self.deconv_stack
        else:
            self.probe_stack = self.fishimage.probes[0]
            

        try:
            self.annotation_type = self.config.annotation_type
        except KeyError:
            self.annotation_type = 'csv'

        if self.annotation_type == 'csv':
            fname = config.annotation_csv_template.format(**spec)
            fpath = os.path.join(config.annotation_dirpath, fname)
            df = pd.read_csv(fpath)
            self.df_good = df[df['counter']==self.config.good_value]
            self.df_bad = df[df['counter']==self.config.bad_value]
            try:
                self.df_nuc = df[df['counter']==self.config.nuc_value]
            except KeyError:
                self.df_nuc = None
            self.sliceImage = self.maxproj
        else:
            sl = get_slice(config, spec)
            
            self.good_mask_attr = mask_from_template_and_spec(
                config.good_template,
                config,
                spec,
                sl
            )
            self.bad_mask_attr = mask_from_template_and_spec(
                config.bad_template,
                config,
                spec,
                sl
            )
            self.nuc_mask_attr = mask_from_template_and_spec(
                config.nuc_template,
                config,
                spec,
                sl
            )
            self.good_dotImage = dotImage_from_template_and_spec(
                config.good_template,
                config,
                spec,
                sl
            )
            self.bad_dotImage = dotImage_from_template_and_spec(
                config.bad_template,
                config,
                spec,
                sl
            )
            self.nuc_dotImage = dotImage_from_template_and_spec(
                config.nuc_template,
                config,
                spec,
                sl
            )
            self.sliceImage = sliceImage_from_template_and_spec(
                config.slice_template,
                config,
                spec,
                sl
            )


    @property
    def maxproj(self):
        return np.max(self.probe_stack, axis=2).view(dbiImage)

    @property
    def good_mask(self):
        if self.annotation_type == 'csv':
            df_good_matrix = self.maxproj*0
            df_dim = df_good_matrix.shape
            for index, row in self.df_good.iterrows():
                for top in range(-5,6):
                    if (row["y"]+top)>=0 and (row["y"]+top)<df_dim[0]:
                        for right in range(-5,6):
                            if (row["x"]+right)>=0 and (row["x"]+right)<df_dim[1]:
                                df_good_matrix[row["y"]+top,row["x"]+right] = 1
            return df_good_matrix
        else:
            return self.good_mask_attr

    @property
    def bad_mask(self):
        if self.annotation_type == 'csv':
            df_bad_matrix = self.maxproj*0
            df_dim = df_bad_matrix.shape
            for index, row in self.df_bad.iterrows():
                for top in range(-3,4):
                    if (row["y"]+top)>=0 and (row["y"]+top)<df_dim[0]:
                        for right in range(-3,4):
                            if (row["x"]+right)>=0 and (row["x"]+right)<df_dim[1]:
                                df_bad_matrix[row["y"]+top,row["x"]+right] = 1
            return df_bad_matrix
        else:
            return self.bad_mask_attr

    @property
    def nuc_mask(self):
        if self.annotation_type == 'csv':
            df_nuc_matrix = self.maxproj*0
            df_dim = df_nuc_matrix.shape
            if self.df_nuc is not None:
                for index, row in self.df_nuc.iterrows():
                    for top in range(-6,7):
                        if (row["y"]+top)>=0 and (row["y"]+top)<df_dim[0]:
                            for right in range(-6,7):
                                if (row["x"]+right)>=0 and (row["x"]+right)<df_dim[1]:
                                    df_nuc_matrix[row["y"]+top,row["x"]+right] = 1
            return df_nuc_matrix
        else:
            return self.nuc_mask_attr

    @property
    def nuc_centroids(self):
        if self.annotation_type == 'csv':
            icentroids = [(row["y"], row["x"]) for index, row in self.df_nuc.iterrows()]
        else:
            scaled_nuc_mask = scale_segmentation(self.nuc_mask, self.maxproj) # scale nuc dots image to original image size
            labelled_points = skimage.measure.label(scaled_nuc_mask) # Give each (non-touching) dot a different integer label
            rprops = skimage.measure.regionprops(labelled_points) # Calculate the geometric properties (including centroids) of each of these labelled regions
            region_centroids = [r.centroid for r in rprops] # Take the centroids from the properties for each region
            icentroids = [(int(r), int(c)) for r, c in region_centroids] # Round the centroids to integer coordinates
        return icentroids
    

    @property
    def icentroids(self):
        if self.annotation_type == 'csv':
            icentroids = [(row["y"], row["x"]) for index, row in self.df_good.iterrows()]
        else:
            scaled_good_mask = scale_segmentation(self.good_mask, self.maxproj) # scale good dots image to original image size
            labelled_points = skimage.measure.label(scaled_good_mask) # Give each (non-touching) dot a different integer label
            rprops = skimage.measure.regionprops(labelled_points) # Calculate the geometric properties (including centroids) of each of these labelled regions
            region_centroids = [r.centroid for r in rprops] # Take the centroids from the properties for each region
            icentroids = [(int(r), int(c)) for r, c in region_centroids] # Round the centroids to integer coordinates
        return icentroids


    @property
    def all_mask(self):
        #if self.annotation_type == 'csv':
        #    df_all = self.df_good.append(self.df_bad)
        #    df_all_matrix = self.maxproj*0
        #    for index, row in df_all.iterrows():
        #        df_all_matrix[row["y"],row["x"]] = 1
        #    return df_all_matrix
        #else:
        return self.bad_mask ^ self.good_mask

    @property
    def scaled_markers(self):
        if self.annotation_type == 'csv':
            scaled_markers = self.all_mask
        else:
            scaled_markers = scale_segmentation(self.all_mask, self.maxproj)
        return scaled_markers

    @property
    def scaled_good_dotImage(self):
        if self.annotation_type == 'csv':
            return self.good_mask*255
        else:
            scaled_markers = scale_segmentation(self.good_dotImage, self.maxproj)
            return scaled_markers

    @property
    def scaled_bad_dotImage(self):
        if self.annotation_type == 'csv':
            return self.bad_mask*255
        else:
            scaled_markers = scale_segmentation(self.bad_dotImage, self.maxproj)
            return scaled_markers

    @property
    def scaled_nuc_dotImage(self):
        if self.annotation_type == 'csv':
            return self.nuc_mask*255
        else:
            scaled_markers = scale_segmentation(self.nuc_dotImage, self.maxproj)
            return scaled_markers

    @property
    def scaled_sliceImage(self):
        if self.annotation_type == 'csv':
            allDotMask=((self.bad_mask | self.good_mask | self.nuc_mask)<1)
            intensity_scaled_sliceImage = self.sliceImage/np.max(self.sliceImage)*65535
            scaled_dim1 = (intensity_scaled_sliceImage) * allDotMask
            info = np.iinfo(self.sliceImage.dtype)
            return np.dstack((scaled_dim1,scaled_dim1,scaled_dim1))
        else:
            allDotMask=~(self.bad_mask | self.good_mask | self.nuc_mask)
            scaledDotMask=scale_segmentation(allDotMask, self.maxproj)
            scaled_dim1 = scale_segmentation(self.sliceImage[:,:,0], self.maxproj)*scaledDotMask
            scaled_dim2 = scale_segmentation(self.sliceImage[:,:,1], self.maxproj)*scaledDotMask
            scaled_dim3 = scale_segmentation(self.sliceImage[:,:,2], self.maxproj)*scaledDotMask
            return np.dstack((scaled_dim1,scaled_dim2,scaled_dim3))

    def cell_mask(self, params):
        cell_mask = cell_mask_from_fishimage(
            self.fishimage, params
        ).view(dbiImage)

        return cell_mask

    def probe_locs_2d(self, thresh=100, ball_size=1):
        rdim, cdim, zdim = self.probe_stack.shape
        probe_locs_3d = find_probe_locations_3d(self.probe_stack[1:(rdim-1),1:(cdim-1),1:(zdim-1)], thresh, ball_size)
        probe_locs_2d = [(r, c) for r, c, z in probe_locs_3d]

        return probe_locs_2d

    def probe_locs_3d(self, thresh=100, ball_size=1):
        rdim, cdim, zdim = self.probe_stack.shape
        probe_locs_3d = find_probe_locations_3d(self.probe_stack[1:(rdim-1),1:(cdim-1),1:(zdim-1)], thresh, ball_size)

        return probe_locs_3d


def load_multiannotation_di(config, spec, use_deconv=False):

    di = MultiAnnotationDataItem(config, spec, use_deconv)

    return di
