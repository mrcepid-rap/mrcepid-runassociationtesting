from ..ingest_data import *


class BurdenIngest(IngestData):

    def __init__(self, parsed_options: dict, required_options: set):
        super().__init__(parsed_options, required_options)

        # Put additional options required by this specific package here
