import logging
import pathlib
from types import SimpleNamespace

import click
import parse
import dtoolcore
import skimage.measure
import pandas as pd

from dtoolbioimage import Image as dbiImage

from config import Config
from data import get_specs
from segment import segmentation_from_cellmask_and_label_image, scale_segmentation, filter_segmentation_by_region_list
from vis import visualise_counts
from probes import get_counts_by_cell


logger = logging.getLogger("fishtools")


def get_filtered_segmentation(dataitem, params):

    segmentation = segmentation_from_cellmask_and_label_image(
        dataitem.cell_mask(params),
        skimage.measure.label(dataitem.scaled_markers)
    )

    good_regions = [segmentation[r, c] for r, c in dataitem.icentroids] # Find the segmentation values (==labels) at the pixel centroids of each of the good dots

    nuc_regions = [segmentation[r, c] for r, c in dataitem.nuc_centroids] # Find the segmentation values (==labels) at the pixel centroids of each of the nuc dots

    good_nuc_regions = []
    for i in good_regions:
        if sum([j==i for j in nuc_regions])==1:
            good_nuc_regions.append(i)

    filtered_segmentation = filter_segmentation_by_region_list(
        segmentation,
        good_nuc_regions
    )

    return filtered_segmentation


def visualise_marker_positions(dataitem):

    r = dataitem.bad_mask
    g = dataitem.good_mask
    b = dataitem.nuc_mask

    import numpy as np
    from dtoolbioimage import scale_to_uint8
    merged = np.dstack([r, g, b])
    scaled = scale_to_uint8(scale_segmentation(merged, dataitem.maxproj))
    maxproj_rgb = scale_to_uint8(np.dstack(3 * [dataitem.maxproj]))
    scaled.view(dbiImage).save("scaled.png")
    v = 0.5 * scaled + 0.5 * maxproj_rgb
    v.view(dbiImage).save("v.png")

    visd = np.array(vis).view(dbiImage)
    (0.5 * visd + 0.5 * scaled).view(dbiImage).save("svis.png")


def process_dataitem(dataitem, spec, params, config, output_ds):

    probe_locs = dataitem.probe_locs_2d(params.probethresh, params.ball_size)
    filtered_segmentation = get_filtered_segmentation(dataitem, params)
    
    centroids_by_cell = {
        idx: {tuple(p) for p in filtered_segmentation.rprops[idx].coords} & set(probe_locs)
        for idx in filtered_segmentation.labels
    }
    
    vis = visualise_counts(
        dataitem,
        filtered_segmentation,
        centroids_by_cell
    )

    # FIXME
    output_fname = "vis{expid}.png".format(**spec)
    image_abspath = output_ds.prepare_staging_abspath_promise(
        f"images/{output_fname}")
    vis.save(image_abspath)


    probe_locs_3d = dataitem.probe_locs_3d(params.probethresh, params.ball_size)
    
    centroids_3d_by_cell={}
    for idx, centroids in centroids_by_cell.items():
        p3d=set()
        for p in probe_locs_3d:
            if tuple(p[0:2]) in centroids:
                p3d.add(tuple(p))
        centroids_3d_by_cell[idx] = p3d
        
    

    import numpy as np
    locations = {'label':[], 'xpoint':[], 'ypoint':[], 'zpoint':[]}
    for idx, centroids in centroids_3d_by_cell.items():
        for p in centroids:
            locations['label'].append(idx)
            locations['xpoint'].append(p[0]+int(params.ball_size/2))
            locations['ypoint'].append(p[1]+int(params.ball_size/2))
            locations['zpoint'].append(p[2]+1+int(params.ball_size/2))

    
    df = pd.DataFrame(locations)
    # FIXME
    csv_locations_output_fname = "locations{expid}.csv".format(**spec)
    csv_locations_abspath = output_ds.prepare_staging_abspath_promise(
        f"csv_locations/{csv_locations_output_fname}")
    df.to_csv(csv_locations_abspath, index=False)

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

    all_specs = get_specs(config)

    specs = all_specs

    from data import load_multiannotation_di
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
                use_deconv = config.deconv_fname_template
                use_deconv = config.deconv_dirpath
                use_deconv = config.use_deconv
            except KeyError:
                use_deconv = False
                logger.info("Not using deconvolution. Option not set or missing optional arguments: use_deconv, deconv_dirpath, deconv_fname_template")
                
            try:
                dataitem = load_multiannotation_di(config, spec, use_deconv)
                df = process_dataitem(
                    dataitem, spec, params, config, output_ds)
                df['expid'] = spec['expid']
                dfs.append(df)
            except (FileNotFoundError, IndexError) as err:
                logger.warning(f"Couldn't load: {err}")

        summary_output_abspath = output_ds.prepare_staging_abspath_promise(
            f"summary.csv")
        pd.concat(dfs).to_csv(summary_output_abspath, index=False)

        output_ds.put_readme(readme_str)


if __name__ == "__main__":
    main()
