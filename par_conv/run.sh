#!/bin/bash

# mpich2
mpiexec -machinefile $1 -np `expr $2 \* $3 + 1` ./dist/Debug/GNU-Linux-x86/par_conv $2 $3

