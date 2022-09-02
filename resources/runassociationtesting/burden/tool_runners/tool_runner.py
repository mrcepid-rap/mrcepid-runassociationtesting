from abc import ABC, abstractmethod

class ToolRunner(ABC):

    def __init__(self):
        pass

    @abstractmethod
    def run_tool(self):
        pass
