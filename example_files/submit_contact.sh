#!/bin/bash -l
#
#$ -l h_rt=6:00:00
#$ -j y
#$ -N Contacts
#$ -pe omp 16
#$ -V

module load gcc/11.2.0

OMP_NUM_THREADS=16
export OMP_NUM_THREADS

time ./contact



