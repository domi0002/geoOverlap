# geoOverlap
 1D Overset with OpenFOAM - tested with OpenFOAM-v1812


Case Directory
---------------
Overlapping 1D region case can be found in run/

Meshing
--------
Meshes are 1D with empty patches on the y and the z directions. Navigate to run/meshes/ and there are two regions region1 and region2. Execute blockMesh inside each of these regions.
Once done with meshing, run the script ./createLinks.sh(Just once!). This establishes a link between the meshes inside meshes/ and the root folder inside constant/ 

Running
--------
onedOverset (basic run)

onedOverset -cons (with a correction)

Outputs
-------
1. Solution on both meshes written in solution.0 and solution.1
2. The matrix and the RHS also written in matrix.txt
3. L2 norm of the error and conservation error also reported for each run

