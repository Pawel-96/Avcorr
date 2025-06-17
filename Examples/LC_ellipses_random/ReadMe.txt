Example run for LC_ellipses mode, one catalog.

To recreate this example:
- remove everything in Results/ directory
- copy here the compiled executable (Avcorr.exe) located in ../../
- run with bash script (firstly ensuring if it has permissions (chmod +700 ./run)):
    ./run N
    where N - number of threads

The code runs for 15x15 [r_parallel, r_perpendicular] ellipsoid axes, spaced between 2 and 70 Mpc/h.
Bash script runs the code and plots Xi_2 and S_3, using Plot_heatmap_Xi.py, Plot_heatmap_Xi_extended.py and Plot_heatmap_S.py python scripts.
For interpretation, check: https://arxiv.org/abs/2506.01762

Data file and random file is located in Data/ and Randoms/, respectively.
These are created with LCto3D.py script converting 2D ra,dec,z lightcone from ../Angular/Data/ into 3D cartesian
comoving coordinates. Results are in Results/ directory
The Results/*merrout* files contain statistics results, while Results/*nnR_XXX_*.fndhst contain histogram counts for XXX-th scale
For more information, check main readme file.
