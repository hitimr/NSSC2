valgrind --leak-check=full \
         --show-leak-kinds=all \
         --track-origins=yes \
         --log-file=../out/valgrind-out.txt \
         ../out/build/jacobiMPI 100 10