from types import SimpleNamespace

import click

import skimage.measure

from fishtools.data import DataLoader, DataItem
from fishtools.segment import (
    segmentation_from_nuclear_channel_and_markers,
    segmentation_from_cellmask_and_label_image
)
from fishtools.vis import visualise_counts

@click.command()
def main():
    config = {
        "deconv_dirpath": "/Users/mhartley/Dropbox/fca alleles project/smFISH Venus probes/smFISH Venus deconvolved/fca-3",
        "deconv_fname_template": "fca3-FLCVenus-VenusRNA{n}_fca3-FLCVenus-VenusRNA{n}_0_qmle_ch01.tif",
        "ids_uri": "/Users/mhartley/data_repo/fish_test_ids/",
        "image_name_template": "fca3-FLCVenus-VenusRNA{n}",
        "series_name_template": "fca3-FLCVenus-VenusRNA{n}.czi #1",
        "annotation_dirpath": "local-data/fca3",
        "annotation_template": "fca3_{n}.png",
        "bad_col": (10, 14, 96),
        "good_col": (103, 20, 0)
    }

    dl = DataLoader(config)
    n = 14
    dataitem = dl.load_by_specifier(n=n)

    params = SimpleNamespace(ks=11, bs=191, sigma=2)

    nuc_label_image = segmentation_from_nuclear_channel_and_markers(
        dataitem.fishimage,
        skimage.measure.label(dataitem.scaled_markers),
        params
    )

    segmentation = segmentation_from_cellmask_and_label_image(
        dataitem.cell_mask(params),
        nuc_label_image
    )

    # import skimage.measure

    # scaled_good_mask = scale_segmentation(dataitem.good_mask, dataitem.maxproj)
    # labelled_points = skimage.measure.label(scaled_good_mask)
    # rprops = skimage.measure.regionprops(labelled_points)
    # region_centroids = [r.centroid for r in rprops]
    # icentroids = [(int(r), int(c)) for r, c in region_centroids]

    vis = visualise_counts(dataitem.maxproj, segmentation, dataitem.probe_locs_2d(75))

    vis.save(f"vis{n}.png")



if __name__ == "__main__":
    main()