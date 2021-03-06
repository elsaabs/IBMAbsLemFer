# IBMAbsLemFer
C++ and R codes to run the 3 PDMP models and plot Figures 1 to 4.

"non-spatial-PDMP" is the code of the non spatial PDMP (used for Figure 1).
"spatial-PDMP-pairwise-contest" is the code of the spatial PDMP with 2 coexisting/competing microbial strains (one large resident and one rare mutant) (used for Figures 2 and 3).
"spatial-PDMP-one-strain" is the spatial PDMP with 1 microbial strain to get average decomposition rate and soil C stock values per phi value (used for Figure 4).

All models should contain "totaux" and "communs" folders in "result" folder.
Each simulation generates one output file going into "totaux" folder, which contains total number of individuals (M) and molecules (C,D,Z) at each writing step (can be larger than time step) (example of use in figures 1 and 2).
Each simulation also generates one output file going into "communs" folder, which contains outputs calculated and defined in "main.cpp". A "scoresCommuns" output file can be generated after the end of simulation (see below for how) - to combine info of several simulations that differ only by the random generator (but same parameter values) (example of use in figures 3 and 4).
Spatial models additionally generate for each simulation one output file going into "cases", which contains the number of individuals (M) and molecules (C,D,Z) in each microsite of the whole lattice to plot grids (example of use in figure 2).

To test a new value of argument in one of the codes:
1- change value in “description file”.
2- check path in "script" (if you creates a new folder for the new test, make sure name of folder in path is the right one).

To run one test on the cluster:
1- change directory to "result" folder.
2- execute the 3 lines of code in "to-compile-on-cluster" (without bash).
3- change directory back to main folder.
4- run code (e.g. "condor_submit description_file").

To generate a "scoresCommuns" file:
1- adjust min(IDcluster), max(IDcluster), max(IDjob) in "scriptcreationScores”.
2- “bash scriptcreationScores” (directory = main folder) (will run "creationScores.cpp").
The file “scoresCommuns.txt” will be created in the main folder.

The code "non-spatial-PDMP" was ran on computer, while codes "spatial-PDMP-pairwise-contest" and "spatial-PDMP-one-strain" were ran on cluster to test several values of "Diff" and "phi*" respectively. To run the 2 latter on computer, delete "description file" and "to-compile-on-cluster", and create a file "script" like "scripttest" in "non-spatial-PDMP".

The outputs used for figures can be found in the "plots-..." folder ("result" folders are empty).

Note: check path in R files. The written ones correspond to codes folders saved on Desktop.
