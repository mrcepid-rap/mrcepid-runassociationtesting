#!/usr/bin/env python3

import json
import pandas as pd

from pathlib import Path

with Path('test_data/expected_results.tsv').open('r') as expected_tsv,\
        Path('test_data/expected_results.json').open('w') as out_json:

    expected_table = pd.read_csv(expected_tsv, header=[0, 1], sep='\t')
    final_array = []

    for row in expected_table.iterrows():
        data = row[1]
        current_dict = {'test_name': data.test_name.test_name, 'parameters': {}, 'outputs': {}}
        for parameter in data.parameters.items():
            if pd.isna(parameter[1]):
                current_dict['parameters'][parameter[0]] = None
            elif type(parameter[1]) == float:
                current_dict['parameters'][parameter[0]] = f'{parameter[1]:0.0f}'
            else:
                current_dict['parameters'][parameter[0]] = parameter[1]
        for output in data.outputs.items():
            if pd.isna(output[1]):
                current_dict['outputs'][output[0]] = None
            else:
                current_dict['outputs'][output[0]] = output[1]
        final_array.append(current_dict)

    json.dump(final_array, out_json, indent=1)
