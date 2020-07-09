# geoOverlap
 1D Overset with OpenFOAM

Meshing
--------
Meshes are 1D with empty patches on the y and the z directions. Navigate to meshes/ and there are two regions region1 and region2. Execute blockMesh inside each of these regions.

Running
--------
onedOverset (basic run)

onedOverset -cons (with a correction)

Outputs
-------
1. Solution on both meshes written in solution.0 and solution.1
2. The matrix and the RHS also written in matrix.txt
3. L2 norm of the error and conservation error also reported for each run

