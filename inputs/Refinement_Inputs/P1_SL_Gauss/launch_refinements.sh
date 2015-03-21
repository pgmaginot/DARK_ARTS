#!/bin/csh

bsub -n 1 openmpi-mpirun ~/Research/dark_arts_local/build/dark_arts Refinement_8_Cell.xml
bsub -n 1 openmpi-mpirun ~/Research/dark_arts_local/build/dark_arts Refinement_16_Cell.xml
bsub -n 1 openmpi-mpirun ~/Research/dark_arts_local/build/dark_arts Refinement_32_Cell.xml
bsub -n 1 openmpi-mpirun ~/Research/dark_arts_local/build/dark_arts Refinement_64_Cell.xml
bsub -n 1 openmpi-mpirun ~/Research/dark_arts_local/build/dark_arts Refinement_128_Cell.xml
bsub -n 1 openmpi-mpirun ~/Research/dark_arts_local/build/dark_arts Refinement_256_Cell.xml
bsub -n 1 openmpi-mpirun ~/Research/dark_arts_local/build/dark_arts Refinement_512_Cell.xml
