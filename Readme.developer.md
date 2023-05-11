# MRCEPID-RunAssociationTesting Developer Readme

## Testing

To run tests for this applet:

1. `cd` into the current module `test/` directory.

2. Check to make sure the `test_data/` folder has been built properly. You can find a list of required files for testing
in the header of the `runassotiationtesting_test.py` script. If these files are not available, please run two scripts:
   * the R script `generate_test_data.R` will generate all required data files
   * The python script can regenerate the input parameters required for the `test_load_general_options()` test so long
   as the file `expected_results.tsv` is still in the `test_data/` directory.
   * The `transcripts.tsv.gz` file cannot be regenerated here, but a description of how to generate this file is 
   available in the [QC_Workflow](https://github.com/mrcepid-rap/QC_workflow) repository in this project.

3. Acquire the test_launch.py script available on github: 

```
git clone https://github.com/mrcepid-rap/mrcepid-testing
cd mrcepid-testing/
```

4. Run the testing script

```{commandline}
# --script is the pytest compatible script
# --files are the test data required for testing
# --root_dir is the path to the root directory containing the source code for mrcepid-runassociationtesting (could be `..` based on instructions above)
# --modules are modules required for the current test. A branch (e.g., v1.1.0) of a given module can be requested using syntax like: general_utilities:v1.1.0 
./test_launch.py --script /path/to/mrcepid-runassociationtesting/test/runassociationtesting_test.py \ 
   --files /path/to/mrcepid-runassociationtesting/test/test_data/ \
   --root_dir /path/to/mrcepid-runassociationtesting/ \
   --modules general_utilities
   --add_opts mode:burden input_args:
```

Please see the [mrcepid-testing repository](https://github.com/mrcepid-rap/mrcepid-testing) for more details on running 
tests.

## Implementing your own modules

T.B.D.
