import io

import ruamel.yaml


class Config(object):

    def __init__(self, config_fpath):
        yaml = ruamel.yaml.YAML()

        with open(config_fpath) as fh:
            self.raw_config = yaml.load(fh)

    def __getattr__(self, name):
        return self.raw_config[name]

    def as_readme_format(self):
        yaml = ruamel.yaml.YAML()

        readme_dict = {
            "config": self.raw_config
        }

        with io.StringIO() as f:
            yaml.dump(readme_dict, f)
            readme_content = f.getvalue()

        return readme_content
