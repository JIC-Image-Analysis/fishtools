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
    image_abspath = output_ds.prepare_staging_abspath_promise(f"images/{output_fname}")
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
    csv_abspath = output_ds.prepare_staging_abspath_promise(f"csv/{csv_output_fname}")
    df.to_csv(csv_abspath, index=False)

    return df


def diagnostics(dataitem, spec, config, params):

    import os
    import numpy as np
    from fishtools.segment import nuc_cell_mask_from_fishimage
    import skimage.measure
    from dtoolbioimage import scale_to_uint8

    ncm = nuc_cell_mask_from_fishimage(dataitem.fishimage, params)
    ncm.view(dbiImage).save("ncm.png")

    template = "{expid}-{expid}.png"
    # template = "{expid}-good.png"
    fname = template.format(**spec)
    fpath = os.path.join(config.annotation_dirpath, fname)
    im = dbiImage.from_file(fpath)

    # template = "{expid}-{expid}.png"
    template = "{expid}-good.png"
    fname = template.format(**spec)
    fpath = os.path.join(config.annotation_dirpath, fname)
    imgood = dbiImage.from_file(fpath)

    im.save("floop.png")
    print(im[0, 0, :])
    print(im[480, 360, :])
    regionim = (im[:,:,3] == 255)
    r = skimage.measure.regionprops(skimage.measure.label(regionim))[0]
    rmin, cmin, rmax, cmax = r.bbox
    sliceme = np.s_[rmin:rmax,cmin:cmax]

    imgood[sliceme].save("SLIGAES.png")


    dataitem.scaled_markers.view(dbiImage).save("scaled_markers.png")

    rdim, cdim = dataitem.maxproj.shape
    canvas = np.dstack(3 * [scale_to_uint8(dataitem.maxproj)])
    canvas[np.where(dataitem.scaled_markers)] = 0, 255, 0
    canvas.view(dbiImage).save("canvas.png")



@click.command()
@click.argument('config_fpath')
def main(config_fpath):

    logging.basicConfig(level=logging.INFO)
    config = Config(config_fpath)
    params = SimpleNamespace(**config.params)

    dl = DataLoader(config.raw_config)

    all_specs = get_specs(config)

    import random
    # specs = random.sample(all_specs, 10)
    specs = all_specs

    # print(specs)
    # import sys; sys.exit(0)

    # print(specs)
    # from dtoolbioimage import ImageDataSet
    # ids = ImageDataSet(config.ids_uri)
    # print(ids.all_possible_stack_tuples())
    # specs = get_specs(config)

    # spec = specs[0]
    # df = dl.load_by_specifier(**spec)


    from fishtools.data import load_wubbly
    from dtoolbioimage import Image as dbiImage

    # TODO - some visual info dumping
    # di = load_wubbly(config, specs[0])
    # diagnostics(di, specs[0], config, params)
    # di.maxproj.view(dbiImage).save("max.png")

    # filt = get_filtered_segmentation(di, params)

    # filt.pretty_color_image.view(dbiImage).save("seg.png")

    # process_dataitem(di, spec, params, config, None)


    readme_str = config.as_readme_format()

    dfs = []
    with dtoolcore.DataSetCreator(
        config.output_name,
        config.output_base_uri
    ) as output_ds:
        for spec in specs:
            # FIXME
            logger.info("Processing n={expid}".format(**spec))
            try:
                # dataitem = dl.load_by_specifier(**spec)
                # FIXME - naming!
                dataitem = load_wubbly(config, spec)
                df = process_dataitem(dataitem, spec, params, config, output_ds)
                df['expid'] = spec['expid']
                dfs.append(df)
            except FileNotFoundError as err:
                logger.warning(f"Couldn't load: {err}")

        summary_output_abspath = output_ds.prepare_staging_abspath_promise(f"summary.csv")
        pd.concat(dfs).to_csv(summary_output_abspath, index=False)

        output_ds.put_readme(readme_str)


if __name__ == "__main__":
    main()
