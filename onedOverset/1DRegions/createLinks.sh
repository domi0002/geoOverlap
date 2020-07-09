#!/bin/bash

cd constant/region1
ln -s -f ../../meshes/region1/constant/polyMesh .

cd ../../constant/region2
ln -s -f ../../meshes/region2/constant/polyMesh .

