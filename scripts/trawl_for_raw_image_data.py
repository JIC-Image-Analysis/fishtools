import os
import shutil
import logging

import click
import dtoolcore


logger = logging.getLogger(__file__)


def is_image(filename, image_exts=['.czi']):
    
    _, ext = os.path.splitext(filename)

    return ext in image_exts


@click.command()
@click.argument('root_dirpath')
@click.argument('output_base_uri')
@click.argument('output_name')
def main(root_dirpath, output_base_uri, output_name):
    
    logging.basicConfig(level=logging.INFO)


    dirpaths_fns = []
    for dirpath, dirnames, filenames in os.walk(root_dirpath):
        for fn in filenames:
            if is_image(fn):
                dirpaths_fns.append((dirpath, fn))


    with dtoolcore.DataSetCreator(output_name, output_base_uri) as output_ds:
        for dirpath, filename in dirpaths_fns:
            basedir = os.path.basename(dirpath)
            relpath = f"{basedir}/{filename}"
            src_abspath = os.path.join(dirpath, filename)
            dst_abspath = output_ds.prepare_staging_abspath_promise(relpath)

            shutil.copy(src_abspath, dst_abspath)





        

if __name__ == "__main__":
    main()