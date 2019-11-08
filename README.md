# IBMAbsLemFer
C++ and R codes for running the microbial individual-based model and plotting Figures 1 to 4.

To test a new value of argument in one of the codes:
1- change value in “description file”
2- check path in "script" (if you creates a new folder for the new test, make sure name of folder in path is the right one)

To run one test on the cluster:
1- change directory to "result" folder
2- execute the 3 lines of code in "to-compile-on-cluster" (without bash)
3- change directory back to main folder
4- run code (e.g. "condor_submit description_file")

To run the code to get "communs" file:
1- adjust min(IDcluster), max(IDcluster), max(IDjob) in "scriptcreationScores”
2- “bash scriptcreationScores” (directory = main folder) (will run "creationScores.cpp")
The file “scoresCommuns.txt” will be created in the main folder.

The code "non-spatial-PDMP" was ran on computer, while codes "spatial-PDMP-pairwise-contest" and "spatial-PDMP-one-strain" were ran on cluster to test several values of "Diff" and "phi*" respectively. To run the 2 latter on computer, delete "description file" and "to-compile-on-cluster", and create a file "script" like "scripttest" in "non-spatial-PDMP".

The outputs used for figures can be found in the "plots-..." folder ("result" folders are eempty).
