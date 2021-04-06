mpirun -np 2 valgrind --leak-check=full \
         --show-leak-kinds=all \
         --track-origins=yes \
         --verbose \
         --log-file=../out/valgrind-out.txt \
         ../out/build/jacobiMPI 100 10