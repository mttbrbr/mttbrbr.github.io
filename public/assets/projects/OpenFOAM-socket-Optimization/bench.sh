#!/bin/sh
# Build every variant and run all timing experiments. Output: CSV on stdout. ~8 minutes.
set -e
mkdir -p bin
b() { name=$1; cc=$2; shift 2; $cc -std=c11 "$@" -DFLAGS="\"$cc $*\"" ldu_cg.c -o bin/$name -lm; }
b gcc_O0            gcc   -O0
b gcc_O2_generic    gcc   -O2 -march=x86-64 -mtune=generic
b gcc_O3            gcc   -O3
b gcc_O3_native     gcc   -O3 -march=native
b gcc_O3_native512  gcc   -O3 -march=native -mprefer-vector-width=512
b gcc_fastmath      gcc   -O3 -march=native -ffast-math
b clang_O3          clang -O3
b clang_O3_native   clang -O3 -march=native
b clang_fastmath    clang -O3 -march=native -ffast-math
b dp                gcc   -O3 -march=native -DPREC_DP
b sp                gcc   -O3 -march=native -DPREC_SP
b spdp              gcc   -O3 -march=native -DPREC_SPDP
echo "kind,precision,flags,n,cells,MB,ns_per_cell_iter,checksum,procs"
for rep in 1 2 3; do
  # 1. compiler flags, solver kernel, in cache (24^3) and in RAM (96^3)
  for n in 24 96; do
    for v in gcc_O0 gcc_O2_generic gcc_O3 gcc_O3_native gcc_O3_native512 gcc_fastmath clang_O3 clang_O3_native clang_fastmath; do
      ./bin/$v $n time | sed 's/$/,1/'
    done
  done
  # 2. precision vs mesh size, solver and field kernels, one process
  for n in 16 24 32 48 64 96 128 160; do
    for v in dp sp spdp; do ./bin/$v $n time | sed 's/$/,1/'; ./bin/$v $n fields | sed 's/$/,1/'; done
  done
  # 3. precision with 4 processes at once (memory bandwidth shared), 96^3 each
  for v in dp sp spdp; do
    for mode in time fields; do
      for k in 0 1 2 3; do taskset -c $k ./bin/$v 96 $mode | sed 's/$/,4/' & done; wait
    done
  done
done
