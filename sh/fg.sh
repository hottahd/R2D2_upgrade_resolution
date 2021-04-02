#!/bin/bash
#PJM -N "upgrade"
#PJM -L "node=1x32x16"
#PJM -L "rscunit=rscunit_ft01"
##PJM -L "rscgrp=small"
#PJM -L "rscgrp=large"
#PJM -L "elapse=5:00:00"
#PJM -o out.txt
#PJM -e err.txt
#PJM --mpi "max-proc-per-node=4"
#PJM --mpi "proc=2048"
#PJM --step

export PARALLEL=12
export OMP_NUM_THREADS=${PARALLEL}

# execute job
mpiexec -n 2048 --of log.txt ./a.out

