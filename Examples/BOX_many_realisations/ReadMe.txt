Example run for BOX mode with 5 realisations.
Code computes the statistics for every scale and realisation,
at the end combines the outputs for all realisations according to Real_template in parameter file.

To recreate this example:
- remove everything in Results/ directory
- copy here the compiled executable (Avcorr.exe) located in ../../
- run with bash script (firstly ensuring if it has permissions (chmod +700 ./run)):
    ./run N
    where N - number of threads

The code runs for 15 scales, spaced between 2 and 80 Mpc/h.
Bash script runs the code and plots Xi_2 and S_3, using Plot_Xi.py and Plot_S.py python scripts.

Data file is located in Data/, results are in Results/ directory
The Results/*merrout* files contain statistics results, while Results/*nnR_XXX_*.fndhst contain histogram counts for XXX-th scale
For more information, check main readme file.
