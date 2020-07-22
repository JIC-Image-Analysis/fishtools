import logging
import pathlib
from types import SimpleNamespace

import click
import parse
import skimage.measure
import pandas as pd

from fishtools.config import Config
from fishtools.data import DataLoader
from fishtools.segment import segmentation_from_nuclear_channel_and_markers, segmentation_from_cellmask_and_label_image, scale_segmentation, filter_segmentation_by_region_list
from fishtools.vis import visualise_counts
from fishtools.probes import get_counts_by_cell


logger = logging.getLogger("fishtools")


def get_filtered_segmentation(dataitem, params):
    nuc_label_image = segmentation_from_nuclear_channel_and_markers(
        dataitem.fishimage,
        skimage.measure.label(dataitem.scaled_markers),
        params
    )

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


def get_specs(config):
    diriter = pathlib.Path(config.annotation_dirpath).iterdir()
    fnameiter = (fpath.name for fpath in diriter)

    logger.debug(f"Matching with {config.annotation_template}")
    
    specs = [
        parse.parse(config.annotation_template, fname).named
        for fname in fnameiter
    ]

    return specs


def process_dataitem(dataitem, spec, params, config):

    probe_locs = dataitem.probe_locs_2d(params.probethresh)
    filtered_segmentation = get_filtered_segmentation(dataitem, params)
    vis = visualise_counts(
        dataitem.maxproj,
        filtered_segmentation,
        probe_locs

    )

    output_dirpath = pathlib.Path(config.output_dirpath)
    output_fname = "vis{n}.png".format(**spec)
    vis.save(output_dirpath/"images"/output_fname)

    areas_by_cell = {
        l: int(filtered_segmentation.rprops[l].area)
        for l in filtered_segmentation.labels
    }
    counts_by_cell = get_counts_by_cell(filtered_segmentation, probe_locs)

    measurements = [
        {
            "label": l,
            "pixelarea": areas_by_cell[l],
            "probecount": counts_by_cell[l]
        }
        for l in areas_by_cell
    ]

    df = pd.DataFrame(measurements)
    csv_output_fname = "results{n}.csv".format(**spec)
    df.to_csv(output_dirpath/"csv"/csv_output_fname, index=False)


@click.command()
@click.argument('config_fpath')
def main(config_fpath):

    logging.basicConfig(level=logging.INFO)
    config = Config(config_fpath)
    params = SimpleNamespace(**config.params)

    dl = DataLoader(config.raw_config)

    specs = get_specs(config)

    # for spec in specs:

    # spec = {'n': 3}

    for spec in specs:
        logger.info("Processing n={n}".format(**spec))

        try:
            dataitem = dl.load_by_specifier(**spec)
            process_dataitem(dataitem, spec, params, config)
        except FileNotFoundError as err:
            logger.warning(f"Couldn't load: {err}")




if __name__ == "__main__":
    main()
