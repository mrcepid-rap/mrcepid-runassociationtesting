from abc import ABC, abstractmethod
from typing import List

import pandas as pd

from runassociationtesting.burden.burden_ingester import BurdenAssociationPack


class ToolRunner(ABC):

    def __init__(self, association_pack: BurdenAssociationPack, output_prefix: str):
        self._association_pack = association_pack
        self._output_prefix = output_prefix
        self._outputs = []

    def get_outputs(self) -> List[str]:
        return self._outputs

    # These two methods help the different tools in defining the correct field names to include in outputs
    @staticmethod
    def _define_field_names_from_pandas(field_one: pd.Series) -> List[str]:

        # Test what columns we have in the 'SNP' field so we can name them...
        field_one = field_one['SNP'].split("-")
        field_names = ['ENST']
        if len(field_one) == 2:  # This is the bare minimum, always name first column ENST, and second column 'var1'
            field_names.append('var1')
        elif len(field_one) == 3:  # This could be the standard naming format... check that column [2] is MAF/AC
            if 'MAF' in field_one[2] or 'AC' in field_one[2]:
                field_names.extend(['MASK', 'MAF'])
            else:  # This means we didn't hit on MAF in column [2] and a different naming convention is used...
                field_names.extend(['var1', 'var2'])
        else:
            for i in range(2, len(field_one) + 1):
                field_names.append('var%i' % i)

        return field_names

    @staticmethod
    def _define_field_names_from_tarball_prefix(tarball_prefix: str, variant_table: pd.DataFrame) -> pd.DataFrame:
        tarball_prefix_split = tarball_prefix.split("-")
        if len(tarball_prefix_split) == 2:  # This could be the standard naming format. Check that column [1] is MAF/AC
            if 'MAF' in tarball_prefix_split[1] or 'AC' in tarball_prefix_split[1]:
                field_names = ['MASK', 'MAF']
                variant_table[field_names] = tarball_prefix_split
            else:  # This means we didn't hit on MAF/AC in column [2] and a different naming convention is used...
                field_names = ['var1', 'var2']
                variant_table[field_names] = tarball_prefix_split
        else:
            for i in range(1, len(tarball_prefix_split) + 1):
                field_name = 'var%i' % i
                variant_table[field_name] = tarball_prefix_split[i - 1]

        return variant_table

    @abstractmethod
    def run_tool(self) -> None:
        pass
