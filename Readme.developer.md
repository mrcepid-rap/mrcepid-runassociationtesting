# MRCEPID-RunAssociationTesting Developer Readme

## Testing

### Running Tests

Tests of this module have been implemented using `pytest`. To run tests follow these instructions:

1. `cd` into the `test/` directory. All tests assume you are executing code from there. 

2. Check to make sure the `test_data/` folder has been built properly. You can find a list of required files for testing
in the header of the `runassotiationtesting_test.py` script. If these files are not available, please run two scripts:
   * the R script `generate_test_data.R` will generate all required data files
   * The python script can regenerate the input parameters required for the `test_load_general_options()` test so long
   as the file `expected_results.tsv` is still in the `test_data/` directory.
   * The `transcripts.tsv.gz` file cannot be regenerated here, but a description of how to generate this file is 
   available in the [QC_Workflow](https://github.com/mrcepid-rap/QC_workflow) repository in this project.

3. Run the `test_launch.py` script like:

```{commandline}
# --script is the pytest compatible script
# --files are the test data required for testing
# --modules are modules required for the current test. A branch (e.g., v1.1.0) of a given module can be requested using syntax like: general_utilities:v1.1.0 
./test_launch.py --script runassociationtesting_test.py --files test_data/ --modules general_utilities
```
       
   This script will:

   1. Collate the resources and files the directory above (`..`; e.g., from mrcepid-runassociationtesting) into a 'test' 
   applet on the DNANexus platform. The applet will automatically be placed into your current project with a name 
   like `runassociationtesting_test_<TIMESTAMP>` where timestamp is a timestamp of the current test start time in the format `YYYYMMDDmmss`.
   2. Run the test applet with the script provided to `--script`.
   3. Retrieve the `pytest` log indicating test(s) pass/fail.
   4. Tear down the temporary files/folders that were created in your project. 

Test results will automatically be downloaded at `test_data/pytest.<TIMESTAMP>.log`. For more information on how 
`test_launch.py`, please see the documentation included in the script.

### Test Implementation

As this applet runs via the DNANexus environment, tests cannot be run on a local machine. Thus, we made several 
modifications to mrcepid-runassociationtesting to enable testing functionality across all modules:

1. Added a `test()` function to `mrcepid-runassociationtesting.py` that allows for testing behaviour
2. Added additional command-line inputs to the applet (`-itesting_script`) and (`-itesting_directory`) that are used by
the `test_launch.py` script detailed above to put `mrcepid-runassociationtesting.py` into testing mode.
3. Created a testing repository with a 'blank' module: [mrcepid-test_loader](https://github.com/mrcepid-rap/mrcepid-test_loader).
This module has no additional functionality beyond that implemented in `mrcepid-runassociationtesting.py`. This module 
was required as mrcepid-runassociationtesting has no way to actually trigger loading of options / data without a module
that implements its interfaces (see [Implementing your own modules](#implementing-your-own-modules) for more information).

## Implementing your own modules

T.B.D.
