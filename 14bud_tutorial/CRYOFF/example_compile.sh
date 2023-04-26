#!/bin/bash

#module load intel/18.0.2 impi/18.0.2 mkl/18.0.2

mpiifort -fpp cry3.0.0.f90 -DMPI -O3 -mavx2 -axAVX,core-AVX2,core-AVX512 -lmkl_lapack95_lp64 -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -liomp5 -lpthread -lmpi -o cry.300.x

