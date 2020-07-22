import pathlib

import click
import parse
import pandas as pd


def load_and_add_metadata(fpath, fname_template):
    fname_metadata = parse.parse(fname_template, fpath.name).named
    df = pd.read_csv(fpath)
    for k, v in fname_metadata.items():
        df[k] = v

    return df


@click.command()
@click.argument("csv_dirpath")
@click.argument("output_fpath")
def main(csv_dirpath, output_fpath):
    fname_template = "results{n}.csv"

    dirpath = pathlib.Path(csv_dirpath)
    df = pd.concat(
        [load_and_add_metadata(fpath, fname_template)
         for fpath in dirpath.iterdir()]
    )

    renamer = {
        "n": "image_number"
    }
    colname_mapper = {cname: cname for cname in df.columns}
    colname_mapper.update(renamer)

    df = df.rename(columns=colname_mapper)

    df.to_csv(output_fpath, index=False)


if __name__ == "__main__":
    main()
