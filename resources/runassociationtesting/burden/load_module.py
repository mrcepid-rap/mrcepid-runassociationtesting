from ..module_loader import *
import burden_ingester


class LoadModule(ModuleLoader):

    def start_module(self) -> None:


        ingested_data = burden_ingester.IngestData(self.inputs, )


