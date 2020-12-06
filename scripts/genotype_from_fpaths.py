import os
import pathlib

import click
import parse

from fishtools.config import Config


def is_image(filename, image_exts=['.czi']):

    _, ext = os.path.splitext(filename)

    return ext in image_exts


@click.command()
@click.argument('config_fpath')
def main(config_fpath):

    config = Config(config_fpath)

    dirpath = pathlib.Path(config.images_root_dirpath)

    dirpaths_fns = []
    for dirpath, dirnames, filenames in os.walk(dirpath):
        for fn in filenames:
            if is_image(fn):
                dirpaths_fns.append((dirpath, fn))

    expid_to_genotype = {}
    image_name_template = "Experiment-{expid:d}.czi"
    for dirpath, fn in dirpaths_fns:
        result = parse.parse(image_name_template, fn)
        expid = result.named['expid']
        expid_to_genotype[expid] = os.path.basename(dirpath)

    
    for expid, genotype in expid_to_genotype.items():
        print(f"{expid}\t{genotype}")




if __name__ == "__main__":
    main()
