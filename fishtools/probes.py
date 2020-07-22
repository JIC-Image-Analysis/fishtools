import skimage.measure
import skimage.morphology

from dtoolbioimage import scale_to_uint8


def find_probe_locations_3d(stack, thresh=100):
    selem = skimage.morphology.ball(1)
    th_stack = skimage.morphology.white_tophat(stack, selem)
    scaled_stack = scale_to_uint8(th_stack)
    labelled = skimage.measure.label(scaled_stack > thresh)
    rprops = skimage.measure.regionprops(labelled)

    centroids_int = [
        list(map(int, r.centroid))
        for r in rprops
    ]

    return centroids_int


def get_counts_by_cell(cell_regions, probe_centroids):
    centroids_by_cell = {
        idx: {tuple(p) for p in cell_regions.rprops[idx].coords} & set(probe_centroids)
        for idx in cell_regions.labels
    }

    counts_by_cell = {
        idx: len(centroids)
        for idx, centroids in centroids_by_cell.items()
    }

    return counts_by_cell
