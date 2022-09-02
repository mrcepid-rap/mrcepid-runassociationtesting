import dxpy
import yaml
from typing import List, TypedDict
from abc import ABC, abstractmethod


class BaseOptions(TypedDict):
    pass


# An abstract class that allows for identical loading of 'modules' as I develop them
class ModuleLoader(ABC):

    def __init__(self, output_prefix: str, yaml_file: dict):
        self.outputs = []
        self.yaml_file = yaml_file
        self.output_prefix = output_prefix
        self.inputs = self._load_yaml()

    def get_outputs(self) -> List[str]:
        return self.outputs

    @abstractmethod
    def start_module(self) -> None:
        pass

    def _load_yaml(self) -> BaseOptions:
        dxpy.download_dxfile(dxpy.DXFile(self.yaml_file).get_id(), 'options.yaml')
        options = yaml.load(open('options.yaml'), Loader=yaml.SafeLoader)
        return options
