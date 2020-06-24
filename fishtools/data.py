from dataclasses import dataclass
from typing import List

from dtoolbioimage import Image, Image3D

from fishtools.utils import extract_nuclei


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


@dataclass
class ManualAnnotatedImage(object):
    raw_image: Image

    @classmethod
    def from_fpath(cls, fpath):
        im = Image.from_file(fpath)

        return cls(im)

    @property
    def marked_cell_mask(self):
        return extract_nuclei(self.raw_image)