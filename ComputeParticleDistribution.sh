#! /bin/bash
set -e
# call Scilab to generate the coordinates of particles
echo "exec('./ComputeParticleDistribution.sce', -1)" | ~/Programs/Unzip_Installed/Scilab/bin/scilab-cli
