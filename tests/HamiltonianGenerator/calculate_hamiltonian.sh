#!/bin/bash

./test_hamiltonian_generator.x -2

n_terms=$(./test_hamiltonian_generator.x -1)
echo "Number of hamiltonian terms:" $n_terms

seq 0 $(($n_terms - 1)) | parallel -j12 ./run_one.sh

#for ((i=0;i<$n_terms;i++))
#do
#   ./test_hamiltonian_generator.x $i > "./hamiltonian/term_${i}.dat"
#   echo $i
#done
