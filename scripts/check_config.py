import logging
import pathlib

import click
import parse

from fishtools.data import DataLoader
from fishtools.config import Config


logger = logging.getLogger("fishtools")


def get_specs(config):
    diriter = pathlib.Path(config.annotation_dirpath).iterdir()
    fnameiter = (fpath.name for fpath in diriter)

    logger.debug(f"Matching with {config.annotation_template}")

    specs = [
        parse.parse(config.annotation_template, fname).named
        for fname in fnameiter
    ]

    return specs


@click.command()
@click.argument('config_fpath')
def main(config_fpath):

    logging.basicConfig(level=logging.INFO)

    config = Config(config_fpath)

    spec = get_specs(config)

    dl = DataLoader(config.raw_config)

    print(dl.ids.all_possible_stack_tuples())


if __name__ == "__main__":
    main()