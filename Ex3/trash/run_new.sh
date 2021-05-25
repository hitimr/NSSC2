#!/bin/bash
Nx=("100")
dt=("0.01" "0.1" "1" "20" "49" "51")
#dt=("49" "51")
num_iter=("100" "500" "1000" "5000" "10000" "20000")
print_mod="1"
task=("1" "2")
for a in "${Nx[@]}"
do
	for b in "${dt[@]}"
	do
		for c in "${num_iter[@]}"
		do
			for d in "${task[@]}"
			do
				if [ $d -eq 1 ]
				then
				echo "python3 task1_1.py $a $b $c $print_mod Nx_${a}_dt_${b}_numIter_${c}_task_${d}.txt"
				python3 task1_1.py $a $b $c $print_mod Nx_${a}_dt_${b}_numIter_${c}_task_${d}.txt
				else
				echo "python3 task1_2.py $a $b $c $print_mod Nx_${a}_dt_${b}_numIter_${c}_task_${d}.txt"
                                python3 task1_2.py $a $b $c $print_mod  Nx_${a}_dt_${b}_numIter_${c}_task_${d}.txt
				fi


			done
		done
	done
done

