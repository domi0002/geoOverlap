# geoOverlap
 1D/2D Overset with OpenFOAM(Conservative) - tested with OpenFOAM-v1812
 
 
External Dependencies(included)
-------------------------------
Octree search: https://github.com/jbehley/octree  (No compilation required)

sparseMatrix: https://github.com/uestla/Sparse-Matrix (OpenFOAM based compilation required inside folder sparseMatrix)

itsol:https://github.com/huiscliu/itsol (standard compilation required inside folder itsol-master)
 
Case Directory
---------------
Overlapping 1D/2D regions case can be found in run/

Meshing
--------
For 1D, meshes are with empty patches on the y and the z directions. Navigate to run/meshes/ and there are two regions region1 and region2. Execute blockMesh inside each of these regions. Once done with meshing, run the script ./createLinks.sh(Just once!). This establishes a link between the meshes inside meshes/ and the root folder inside constant/ 

for 2D, several mesh refinements can be generated by issuing the command
./switchMesh.sh <refinementLevel>
 
where <refinementLeve> is either tiny, coarse, medium or fine.

Running
--------

1D:
onedOverset (basic run)
onedOverset -cons (with a correction)

2D:
twoDOverset  <options>

Various options are as follows:
-neumannrhs:  Solve using Neumann BCs on left and top boundaries (You also need to change the 0/T files appropriately). The default option in 0/ folder is the neumann BCs.
-inv:         Use inverse distance interpolation instead of polynomial interplation.
-cons      :  Activate conservative correction


Outputs
-------
1. For 1D: Solution on both meshes written in solution.0 and solution.1
2. For 1D: The matrix and the RHS also written in matrix.txt
3. L2 norm of the error and conservation error also reported for each run

