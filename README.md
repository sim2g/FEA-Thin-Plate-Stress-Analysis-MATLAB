# FEA-Thin-Plate-Stress-Analysis-MATLAB

## Description
A finite element analysis (FEA) program to solve for displacement and strain of a thin plate in plane stress, using both 3-noded and 6-noded triangular elements. The results of this program are validated using ANSYS Mechanical APDL.


<p align="center">
  <img src="./Project Resources/MATLAB_Deformation.png" alt="graph" width="600"/>
</p>

## Dependencies
To use this script, ensure you have MATLAB installed on your system. No additional libraries are required as the script uses standard MATLAB libraries.

## Usage
To use these functions, run the file ```Thin_Plate_Solver.m``` and the program will return the Strain and Nodal Displacement in each dimension, as well as graphs of the deformed plate for both the 3-noded and 6-noded case.

The geometry and loading conditions of the thin plate in plane stress can be modified in ```Thin_Plate_Solver.m``` for your use case.
