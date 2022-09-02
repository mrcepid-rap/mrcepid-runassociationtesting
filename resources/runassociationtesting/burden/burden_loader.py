from ..module_loader import *
from ..association_pack import *
import burden_ingester
import burden_covariates

from tool_runners.glm_runner import *
from tool_runners.bolt_runner import *
from tool_runners.saige_runner import *
from tool_runners.regenie_runner import *
from tool_runners.staar_runner import *


class LoadModule(ModuleLoader):

    def start_module(self) -> None:

        ingested_data = burden_ingester.IngestData(self.inputs, self.required_options)
        processed_covariates = burden_covariates.CovariateProcessor(ingested_data, self.inputs)

    # Just defines possible tools useable by this module
    def check_tools(self, input_tool, association_pack: AssociationPack):

        module_tools = {'bolt': BOLTRunner(association_pack),
                        'saige': SAIGERunner(association_pack),
                        'regenie': REGENIERunner(association_pack),
                        'staar': STAARRunner(association_pack),
                        'glm': GLMRunner(association_pack)}
        if input_tool in module_tools:
            return module_tools
        else:
            raise dxpy.AppError(f'Tool – {input_tool} – not support. Please try a different input tool!')

