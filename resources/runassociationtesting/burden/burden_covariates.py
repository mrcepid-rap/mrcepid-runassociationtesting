from ..covariate_processor import *


class BurdenCovariates(CovariateProcessor):

    def __init__(self, ingested_data: IngestData, parsed_options: dict):
        super().__init__(ingested_data, parsed_options)

        # Put additional covariate processing specific to this module here
