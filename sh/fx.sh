##PJM -L "rscgrp=fx-extra"
#PJM -L "rscgrp=fx-large"
##PJM -L "rscgrp=fx-small"
#PJM -L "node=128"
#PJM --mpi "proc=512"
#PJM -L "elapse=06:00:00"
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
