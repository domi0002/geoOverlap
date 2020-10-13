#!/bin/bash

if [ $# -eq 0 ]; then
	echo " Usage : ./switchMesh  <meshType>, where meshType is one of "
	echo "         tiny, coarse, medium, fine "
	exit
fi

cd meshes
	cd region1
	blockMesh -dict system/blockMeshDict.$1
	cd ../
	cd region2
	blockMesh -dict system/blockMeshDict.$1
	transformPoints -origin '(0.45 0.45 0)' -rotate-angle '((0 0 1) 25)'
        transformPoints -translate '(-0.2 -0.2 0)'	
        transformPoints -translate '( 0 0.5 0)'	
cd ..
cd ..
rm -rf 0
cp -r 0.org.NBC 0

cp -r meshes/region1/constant/polyMesh constant/region1/.
cp -r meshes/region2/constant/polyMesh constant/region2/.

