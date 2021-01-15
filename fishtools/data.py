import os
import logging
import pathlib

from dataclasses import dataclass
from typing import List
from types import SimpleNamespace

import parse
import numpy as np
import skimage.measure

from dtoolbioimage import (
    Image as dbiImage,
    Image3D,
    ImageDataSet,
    scale_to_uint8
)

from fishtools.utils import extract_nuclei, crop_to_non_empty, select_near_colour
from fishtools.segment import scale_segmentation, cell_mask_from_fishimage
from fishtools.probes import find_probe_locations_3d


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


@dataclass
class DataItem(object):
    fishimage: FISHImage
    deconv_stack: Image3D
    annotation: dbiImage
    config: SimpleNamespace

    @property
    def good_mask(self):
        good_mask = select_near_colour(self.cropped_im, self.config.good_col)
        return good_mask
    
    @property
    def bad_mask(self):
        bad_mask = select_near_colour(self.cropped_im, self.config.bad_col)
        return bad_mask
    
    @property
    def all_mask(self):
        return self.bad_mask ^ self.good_mask
    
    def __post_init__(self):
        small_crop = self.annotation[10:-10,10:-10]
        self.cropped_im = crop_to_non_empty(small_crop)
        
    @property
    def maxproj(self):
        return np.max(self.deconv_stack, axis=2).view(dbiImage)

    @property
    def scaled_markers(self):
        scaled_markers = scale_segmentation(self.all_mask, self.maxproj)
        return scaled_markers
    
    def cell_mask(self, params):
        cell_mask = cell_mask_from_fishimage(self.fishimage, params).view(dbiImage)
        return cell_mask
    
    def probe_locs_2d(self, thresh=100):
        probe_locs_3d = find_probe_locations_3d(self.deconv_stack, thresh)
        probe_locs_2d = [(r, c) for r, c, z in probe_locs_3d]
        
        return probe_locs_2d

    
class DataLoader(object):
    
    def __init__(self, config):
        self.config = SimpleNamespace(**config)
        self.ids = ImageDataSet(self.config.ids_uri)
        
    def load_by_specifier(self, **kwargs):
        fname = self.config.deconv_fname_template.format(**kwargs)
        fpath = os.path.join(self.config.deconv_dirpath, fname)
        deconv_stack = Image3D.from_file(fpath)
        
        nuclear_channel_first = True
        image_name = self.config.image_name_template.format(**kwargs)
        series_name = self.config.series_name_template.format(**kwargs)
        fishimage = FISHImage.from_ids_im_sn(self.ids, image_name, series_name, nuclear_channel_first)
        
        annotation_fname = self.config.annotation_template.format(**kwargs)
        annotation_fpath = os.path.join(self.config.annotation_dirpath, annotation_fname)
        annotation = dbiImage.from_file(annotation_fpath)
        
        return DataItem(fishimage, deconv_stack, annotation, self.config)
    

def get_specs(config):
    diriter = pathlib.Path(config.annotation_dirpath).iterdir()
    fnameiter = (fpath.name for fpath in diriter)

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
    # FIXME - put this in config
    # template = "{expid}-{expid}.png"
    fname = config.slice_template.format(**spec)
    fpath = os.path.join(config.annotation_dirpath, fname)
    im = dbiImage.from_file(fpath)

    regionim = (im[:, :, 3] == 255)
    r = skimage.measure.regionprops(skimage.measure.label(regionim))[0]
    rmin, cmin, rmax, cmax = r.bbox
    return np.s_[rmin:rmax, cmin:cmax]


def mask_from_template_and_spec(template, config, spec, sl):
    fname = template.format(**spec)
    fpath = os.path.join(config.annotation_dirpath, fname)
    im = dbiImage.from_file(fpath)
    maxflatten = np.max(im[:, :, :3], axis=2)
    return (maxflatten > 0)[sl]


class MultiAnnotationDataLoader(DataLoader):

    def load_by_specifier(self, **kwargs):
        pass


class MultiAnnotationDataItem(object):

    def __init__(self, config, spec):
        self.config = config
        self.ids = ImageDataSet(self.config.ids_uri)

        nuclear_channel_first = True
        image_name = self.config.image_name_template.format(**spec)
        series_name = self.config.series_name_template.format(**spec)
        self.fishimage = FISHImage.from_ids_im_sn(
            self.ids, image_name, series_name, nuclear_channel_first
        )

        fname = self.config.deconv_fname_template.format(**spec)
        fpath = os.path.join(self.config.deconv_dirpath, fname)
        self.deconv_stack = Image3D.from_file(fpath)

        self.deconv_stack = np.clip(self.deconv_stack, 0, 10000)

        if self.deconv_stack.shape != self.fishimage.nuclei.shape:
            logger.warning("Deconv stack doesn't match shape, trimming")
            rdim, cdim, zdim = self.fishimage.nuclei.shape
            self.deconv_stack = self.deconv_stack[:rdim,:cdim,:zdim]

        sl = get_slice(config, spec)
        self.good_mask = mask_from_template_and_spec(
            config.good_template,
            config,
            spec,
            sl
        )
        self.bad_mask = mask_from_template_and_spec(
            config.bad_template,
            config,
            spec,
            sl
        )
        self.nuc_mask = mask_from_template_and_spec(
            config.nuc_template,
            config,
            spec,
            sl
        )

    @property
    def all_mask(self):
        # return self.nuc_mask
        return self.bad_mask ^ self.good_mask

    @property
    def maxproj(self):
        return np.max(self.deconv_stack, axis=2).view(dbiImage)

    @property
    def scaled_markers(self):
        # print(f"Scaling to {self.maxproj.shape}")
        # print(f"Diag {self.deconv_stack.shape}, {self.fishimage.nuclei.shape}")
        scaled_markers = scale_segmentation(self.all_mask, self.maxproj)
        return scaled_markers

    def cell_mask(self, params):
        cell_mask = cell_mask_from_fishimage(
            self.fishimage, params
        ).view(dbiImage)

        return cell_mask

    def probe_locs_2d(self, thresh=100):
        probe_locs_3d = find_probe_locations_3d(self.deconv_stack, thresh)
        probe_locs_2d = [(r, c) for r, c, z in probe_locs_3d]

        return probe_locs_2d


def load_multiannotation_di(config, spec):

    di = MultiAnnotationDataItem(config, spec)

    return di