#!/bin/bash

# mpich2
mpiexec -machinefile $1 -np `expr $4 \* $5 + 1` ../dist/Release/GNU-Linux-x86/par_conv $2 $3 $4 $5
