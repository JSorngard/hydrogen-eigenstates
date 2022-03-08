# hydrogen-eigenstates
Program calculating the energies of the eigenstates of hydrogen, as well as writing the specified wavefunction to a .txt file.

It does this using the B-spline method to solve a matrix eigenvalue problem. The quantum numbers n and l are hardcoded in the program, and found at the start.

In order to run this program gfortran, lapack and blas must be present. Compile the program by running "./compile". First compile must be made executable by running "chmod u+x compile". After compilation a new executable file will have been made called "eigensplines.out".
