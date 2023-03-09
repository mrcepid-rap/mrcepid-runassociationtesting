#!/usr/bin/env python3

import pandas as pd

test = pd.read_csv('test_data/expected_results.tsv', header=[0, 1], sep='\t')

final_array = []

print('[')
for row in test.iterrows():
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
    print(f'{current_dict},')

print(']')
