#!/bin/bash
Nx=("100")
dt=("0.1" "1" "20")
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

echo "python3 task1_1.py 100 49 5000 1 Nx_100_dt_49_numIter_5000_task_1.txt"
python3 task1_1.py 100 49 5000 1 Nx_100_dt_49_numIter_5000_task_1.txt
echo "python3 task1_2.py 100 49 5000 1 Nx_100_dt_49_numIter_5000_task_2.txt"
python3 task1_2.py 100 49 5000 1 Nx_100_dt_49_numIter_5000_task_1.txt

