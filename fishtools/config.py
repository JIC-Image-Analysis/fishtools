import ruamel.yaml


class Config(object):

    def __init__(self, config_fpath):

        yaml = ruamel.yaml.YAML()

        with open(config_fpath) as fh:
            self.raw_config = yaml.load(fh)

    def __getattr__(self, name):
        return self.raw_config[name]