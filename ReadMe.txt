1. Put the contents to a new directory.
The folder "codes" contains the scripts.
The folder "data" contains data that are stored in RData format.
The folder "output" will save the outputs files. Subfolders "primary" and "sensAnaly" hold the outputs for the primary and the sensitivity analysis, respectively. 
Outputs for "fit_CxInc", "proj_CxInc", and "CEA" are saved in separate folders, with the date of running the script as suffix of the folder. These folders will be created automatically when running the main script "0_main.R".


2. Please change the working folder "folder_main" in "codes/0_main.R" before running the scripts.
The scripts require some external R packages. Please install the required packages if needed.

3. Parameter inference could take a long time, depending on the performance of the computer.
The folder "output/primary/260430_fit_CxInc" contains the parameter sets that are previously inferred, which is saved in the file "260430_nonCeCx_PSApara.xlsx".
To use this file of parameters for projecting cancer incidence and running CEA, one can run the script "0_main_use_fitted_params.R" directly.
This script skips the codes for estimating parameters for calibrating to non-cervical cancer incidence.
