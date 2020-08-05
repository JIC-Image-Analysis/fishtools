import logging
import pathlib

import click
import parse

from fishtools.data import DataLoader, get_specs
from fishtools.config import Config


logger = logging.getLogger("fishtools")


@click.command()
@click.argument('config_fpath')
def main(config_fpath):

    logging.basicConfig(level=logging.INFO)

    config = Config(config_fpath)

    specs = get_specs(config)

    print(specs)
    dl = DataLoader(config.raw_config)

    print(dl.ids.all_possible_stack_tuples())


if __name__ == "__main__":
    main()