Example run for angular mode, one catalog.

To recreate this example:
- remove everything in Results/ directory
- copy here the compiled executable (Avcorr.exe) located in ../../
- run with bash script (firstly ensuring if it has permissions (chmod +700 ./run)):
    ./run N
    where N - number of threads

The code runs for 25 scales, spaced between 0.1 and 5 deg, on 60 deg x 60 deg catalog.
Bash script runs the code and plots w_2 and S_3, using Plot_w.py and Plot_S.py python scripts.
By w_J we call angular version of xi_J, S_J name is not changed (notation)

Data file is located in Data/, results are in Results/ directory
The Results/*merrout* files contain statistics results, while Results/*nnR_XXX_*.fndhst contain histogram counts for XXX-th scale
For more information, check main readme file.
