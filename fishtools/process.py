import logging
import pathlib
from types import SimpleNamespace

import click
import parse
import dtoolcore
import skimage.measure
import pandas as pd

from dtoolbioimage import Image as dbiImage

from fishtools.config import Config
from fishtools.data import DataLoader, get_specs
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


def process_dataitem(dataitem, spec, params, config, output_ds):

    probe_locs = dataitem.probe_locs_2d(params.probethresh)
    filtered_segmentation = get_filtered_segmentation(dataitem, params)
    vis = visualise_counts(
        dataitem.maxproj,
        filtered_segmentation,
        probe_locs
    )

    # FIXME
    output_fname = "vis{expid}.png".format(**spec)
    image_abspath = output_ds.prepare_staging_abspath_promise(
        f"images/{output_fname}")
    vis.save(image_abspath)

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
    # FIXME
    csv_output_fname = "results{expid}.csv".format(**spec)
    csv_abspath = output_ds.prepare_staging_abspath_promise(
        f"csv/{csv_output_fname}")
    df.to_csv(csv_abspath, index=False)

    return df


@click.command()
@click.argument('config_fpath')
def process_from_config(config_fpath):

    logging.basicConfig(level=logging.INFO)
    config = Config(config_fpath)
    params = SimpleNamespace(**config.params)

    dl = DataLoader(config.raw_config)

    all_specs = get_specs(config)

    specs = all_specs

    from fishtools.data import load_multiannotation_di
    from dtoolbioimage import Image as dbiImage

    readme_str = config.as_readme_format()

    dfs = []
    with dtoolcore.DataSetCreator(
        config.output_name,
        config.output_base_uri
    ) as output_ds:
        for spec in specs:
            logger.info("Processing n={expid}".format(**spec))
            try:
                # FIXME - naming!
                dataitem = load_multiannotation_di(config, spec)
                df = process_dataitem(
                    dataitem, spec, params, config, output_ds)
                df['expid'] = spec['expid']
                dfs.append(df)
            except FileNotFoundError as err:
                logger.warning(f"Couldn't load: {err}")

        summary_output_abspath = output_ds.prepare_staging_abspath_promise(
            f"summary.csv")
        pd.concat(dfs).to_csv(summary_output_abspath, index=False)

        output_ds.put_readme(readme_str)


if __name__ == "__main__":
    main()
