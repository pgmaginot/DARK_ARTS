# run the initial files

# bsub -n 1 openmpi-mpirun dark_arts ../inputs/SL_Lobatto_P1_4_Cells.xml 

# sleep "1" 

# bsub -n 1 openmpi-mpirun dark_arts ../inputs/SL_Gauss_P1_4_Cells.xml 

# sleep "1" 

# bsub -n 1 openmpi-mpirun dark_arts ../inputs/SL_Lobatto_P2_4_Cells.xml 

# sleep "1" 

# bsub -n 1 openmpi-mpirun dark_arts ../inputs/SL_Gauss_P2_4_Cells.xml 

# sleep "1" 

# bsub -n 1 openmpi-mpirun dark_arts ../inputs/SL_Lobatto_P3_4_Cells.xml 

# sleep "1" 

# bsub -n 1 openmpi-mpirun dark_arts ../inputs/SL_Gauss_P3_4_Cells.xml 

# sleep "1" 

# bsub -n 1 openmpi-mpirun dark_arts ../inputs/SL_Lobatto_P4_4_Cells.xml 

# sleep "1" 

# bsub -n 1 openmpi-mpirun dark_arts ../inputs/SL_Gauss_P4_4_Cells.xml 

# sleep "30" 

# P1 Refinements

bsub -n 1 openmpi-mpirun dark_arts ../inputs/Refinement_Inputs/P1_SL_Lobatto/Refinement_8_Cell.xml 

sleep "1" 

bsub -n 1 openmpi-mpirun dark_arts ../inputs/Refinement_Inputs/P1_SL_Lobatto/Refinement_16_Cell.xml 

sleep "1" 

bsub -n 1 openmpi-mpirun dark_arts ../inputs/Refinement_Inputs/P1_SL_Lobatto/Refinement_32_Cell.xml 

sleep "1" 

bsub -n 1 openmpi-mpirun dark_arts ../inputs/Refinement_Inputs/P1_SL_Lobatto/Refinement_64_Cell.xml 

sleep "1" 

bsub -n 1 openmpi-mpirun dark_arts ../inputs/Refinement_Inputs/P1_SL_Lobatto/Refinement_128_Cell.xml 

sleep "1" 

bsub -n 1 openmpi-mpirun dark_arts ../inputs/Refinement_Inputs/P1_SL_Lobatto/Refinement_256_Cell.xml 

sleep "1" 

bsub -n 1 openmpi-mpirun dark_arts ../inputs/Refinement_Inputs/P1_SL_Lobatto/Refinement_512_Cell.xml 

sleep "1"  

# Gauss 

bsub -n 1 openmpi-mpirun dark_arts ../inputs/Refinement_Inputs/P1_SL_Gauss/Refinement_8_Cell.xml 

sleep "1" 

bsub -n 1 openmpi-mpirun dark_arts ../inputs/Refinement_Inputs/P1_SL_Gauss/Refinement_16_Cell.xml 

sleep "1" 

bsub -n 1 openmpi-mpirun dark_arts ../inputs/Refinement_Inputs/P1_SL_Gauss/Refinement_32_Cell.xml 

sleep "1" 

bsub -n 1 openmpi-mpirun dark_arts ../inputs/Refinement_Inputs/P1_SL_Gauss/Refinement_64_Cell.xml 

sleep "1" 

bsub -n 1 openmpi-mpirun dark_arts ../inputs/Refinement_Inputs/P1_SL_Gauss/Refinement_128_Cell.xml 

sleep "1" 

bsub -n 1 openmpi-mpirun dark_arts ../inputs/Refinement_Inputs/P1_SL_Gauss/Refinement_256_Cell.xml 

sleep "1" 

bsub -n 1 openmpi-mpirun dark_arts ../inputs/Refinement_Inputs/P1_SL_Gauss/Refinement_512_Cell.xml 

sleep "1" 

# P2 Refinements

bsub -n 1 openmpi-mpirun dark_arts ../inputs/Refinement_Inputs/P2_SL_Lobatto/Refinement_8_Cell.xml 

sleep "1" 

bsub -n 1 openmpi-mpirun dark_arts ../inputs/Refinement_Inputs/P2_SL_Lobatto/Refinement_16_Cell.xml 

sleep "1" 

bsub -n 1 openmpi-mpirun dark_arts ../inputs/Refinement_Inputs/P2_SL_Lobatto/Refinement_32_Cell.xml 

sleep "1" 

bsub -n 1 openmpi-mpirun dark_arts ../inputs/Refinement_Inputs/P2_SL_Lobatto/Refinement_64_Cell.xml 

sleep "1" 

bsub -n 1 openmpi-mpirun dark_arts ../inputs/Refinement_Inputs/P2_SL_Lobatto/Refinement_128_Cell.xml 

sleep "1" 

bsub -n 1 openmpi-mpirun dark_arts ../inputs/Refinement_Inputs/P2_SL_Lobatto/Refinement_256_Cell.xml 

sleep "1" 

bsub -n 1 openmpi-mpirun dark_arts ../inputs/Refinement_Inputs/P2_SL_Lobatto/Refinement_512_Cell.xml 

sleep "1" 

# Gauss 

bsub -n 1 openmpi-mpirun dark_arts ../inputs/Refinement_Inputs/P2_SL_Gauss/Refinement_8_Cell.xml 

sleep "1" 

bsub -n 1 openmpi-mpirun dark_arts ../inputs/Refinement_Inputs/P2_SL_Gauss/Refinement_16_Cell.xml 

sleep "1" 

bsub -n 1 openmpi-mpirun dark_arts ../inputs/Refinement_Inputs/P2_SL_Gauss/Refinement_32_Cell.xml 

sleep "1" 

bsub -n 1 openmpi-mpirun dark_arts ../inputs/Refinement_Inputs/P2_SL_Gauss/Refinement_64_Cell.xml 

sleep "1" 

bsub -n 1 openmpi-mpirun dark_arts ../inputs/Refinement_Inputs/P2_SL_Gauss/Refinement_128_Cell.xml 

sleep "1" 

bsub -n 1 openmpi-mpirun dark_arts ../inputs/Refinement_Inputs/P2_SL_Gauss/Refinement_256_Cell.xml 

sleep "1" 

bsub -n 1 openmpi-mpirun dark_arts ../inputs/Refinement_Inputs/P2_SL_Gauss/Refinement_512_Cell.xml 

sleep "1" 

# P3 Refinements

bsub -n 1 openmpi-mpirun dark_arts ../inputs/Refinement_Inputs/P3_SL_Lobatto/Refinement_8_Cell.xml 

sleep "1" 

bsub -n 1 openmpi-mpirun dark_arts ../inputs/Refinement_Inputs/P3_SL_Lobatto/Refinement_16_Cell.xml 

sleep "1" 

bsub -n 1 openmpi-mpirun dark_arts ../inputs/Refinement_Inputs/P3_SL_Lobatto/Refinement_32_Cell.xml 

sleep "1" 

bsub -n 1 openmpi-mpirun dark_arts ../inputs/Refinement_Inputs/P3_SL_Lobatto/Refinement_64_Cell.xml 

sleep "1" 

bsub -n 1 openmpi-mpirun dark_arts ../inputs/Refinement_Inputs/P3_SL_Lobatto/Refinement_128_Cell.xml 

sleep "1"  

bsub -n 1 openmpi-mpirun dark_arts ../inputs/Refinement_Inputs/P3_SL_Lobatto/Refinement_256_Cell.xml 

sleep "1" 

bsub -n 1 openmpi-mpirun dark_arts ../inputs/Refinement_Inputs/P3_SL_Lobatto/Refinement_512_Cell.xml 

sleep "1" 

# Gauss 

bsub -n 1 openmpi-mpirun dark_arts ../inputs/Refinement_Inputs/P3_SL_Gauss/Refinement_8_Cell.xml 

sleep "1" 

bsub -n 1 openmpi-mpirun dark_arts ../inputs/Refinement_Inputs/P3_SL_Gauss/Refinement_16_Cell.xml 

sleep "1" 

bsub -n 1 openmpi-mpirun dark_arts ../inputs/Refinement_Inputs/P3_SL_Gauss/Refinement_32_Cell.xml 

sleep "1" 

bsub -n 1 openmpi-mpirun dark_arts ../inputs/Refinement_Inputs/P3_SL_Gauss/Refinement_64_Cell.xml 

sleep "1" 

bsub -n 1 openmpi-mpirun dark_arts ../inputs/Refinement_Inputs/P3_SL_Gauss/Refinement_128_Cell.xml 

sleep "1" 

bsub -n 1 openmpi-mpirun dark_arts ../inputs/Refinement_Inputs/P3_SL_Gauss/Refinement_256_Cell.xml 

sleep "1" 

bsub -n 1 openmpi-mpirun dark_arts ../inputs/Refinement_Inputs/P3_SL_Gauss/Refinement_512_Cell.xml 

sleep "1" 

# P4 Refinements

bsub -n 1 openmpi-mpirun dark_arts ../inputs/Refinement_Inputs/P4_SL_Lobatto/Refinement_8_Cell.xml 

sleep "1"  

bsub -n 1 openmpi-mpirun dark_arts ../inputs/Refinement_Inputs/P4_SL_Lobatto/Refinement_16_Cell.xml 

sleep "1" 

bsub -n 1 openmpi-mpirun dark_arts ../inputs/Refinement_Inputs/P4_SL_Lobatto/Refinement_32_Cell.xml 

sleep "1" 

bsub -n 1 openmpi-mpirun dark_arts ../inputs/Refinement_Inputs/P4_SL_Lobatto/Refinement_64_Cell.xml 

sleep "1" 

bsub -n 1 openmpi-mpirun dark_arts ../inputs/Refinement_Inputs/P4_SL_Lobatto/Refinement_128_Cell.xml 

sleep "1" 

bsub -n 1 openmpi-mpirun dark_arts ../inputs/Refinement_Inputs/P4_SL_Lobatto/Refinement_256_Cell.xml 

sleep "1" 

bsub -n 1 openmpi-mpirun dark_arts ../inputs/Refinement_Inputs/P4_SL_Lobatto/Refinement_512_Cell.xml 

sleep "1" 

# Gauss 

bsub -n 1 openmpi-mpirun dark_arts ../inputs/Refinement_Inputs/P4_SL_Gauss/Refinement_8_Cell.xml 

sleep "1" 

bsub -n 1 openmpi-mpirun dark_arts ../inputs/Refinement_Inputs/P4_SL_Gauss/Refinement_16_Cell.xml 

sleep "1" 

bsub -n 1 openmpi-mpirun dark_arts ../inputs/Refinement_Inputs/P4_SL_Gauss/Refinement_32_Cell.xml 

sleep "1"  

bsub -n 1 openmpi-mpirun dark_arts ../inputs/Refinement_Inputs/P4_SL_Gauss/Refinement_64_Cell.xml 

sleep "1" 

bsub -n 1 openmpi-mpirun dark_arts ../inputs/Refinement_Inputs/P4_SL_Gauss/Refinement_128_Cell.xml 

sleep "1" 

bsub -n 1 openmpi-mpirun dark_arts ../inputs/Refinement_Inputs/P4_SL_Gauss/Refinement_256_Cell.xml 

sleep "1" 

bsub -n 1 openmpi-mpirun dark_arts ../inputs/Refinement_Inputs/P4_SL_Gauss/Refinement_512_Cell.xml 

sleep "1" 



