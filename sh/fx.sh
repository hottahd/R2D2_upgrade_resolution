##PJM -L "rscgrp=fx-extra"
##PJM -L "rscgrp=fx-middle"
#PJM -L "rscgrp=fx-small"
#PJM -L "node=8"
#PJM --mpi "proc=32"
#PJM -L "elapse=01:00:00"
#PJM -j
#PJM -o "log.txt"
#PJM --step
#------ Program execution -------#
#export PARALLEL=12
export OMP_NUM_THREADS=12
#rm -rf prof
#rm -rf data
#mpiexec ./a.out
#fapp -C -d prof -I mpi -I hwm mpiexec ./a.out

mpiexec ./a.out
