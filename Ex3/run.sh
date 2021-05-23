#!/bin/bash
Nx=("100")
dt=("0.01" "0.1" "1" "10")
num_iter=("100" "1000" "10000")
print_mod="10"
task=("1" "2")
for a in "${Nx[@]}"
do
	for b in "${dt[@]}"
	do
		for c in "${num_iter[@]}"
		do
			for d in "${task[@]}"
			do
				echo "python3 task1_hickel.py $a $b $c $print_mod $task Nx_${a}_dt_${b}_numIter_${c}_task_${d}.txt"
				python3 task1_hickel.py $a $b $c $print_mod $task Nx_${a}_dt_${b}_numIter_${c}_task_${d}.txt
			done
		done
	done
done

