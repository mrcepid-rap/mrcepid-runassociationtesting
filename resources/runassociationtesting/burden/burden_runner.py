from abc import ABC, abstractmethod
from typing import List


class BurdenRunner(ABC):

    def __init__(self, output_prefix: str, yaml_file: str):
        self.yaml_file = yaml_file

    def get_outputs(self) -> List[str]:
        return self.outputs

    @abstractmethod
    def run_tool(self):
        pass

