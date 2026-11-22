CC     = gcc
MPICC  = mpicc
CFLAGS = -O3 -march=native
NP    ?= 4

all: pingpong halo_big halo_small decompose

pingpong: pingpong.c ; $(MPICC) $(CFLAGS) $< -o $@
halo_big: halo.c decompose.c ; $(MPICC) $(CFLAGS) $< -o $@ -lm
halo_small: halo.c decompose.c ; $(MPICC) $(CFLAGS) -DNX=800 -DNY=400 -DITERS=1000 $< -o $@ -lm
decompose: decompose.c ; $(CC) $(CFLAGS) $< -o $@ -lm

run: all
	mpirun -np 2 --cpu-set 0,1 --bind-to hwthread ./pingpong
	mpirun -np $(NP) --bind-to core ./halo_big
	mpirun -np $(NP) --bind-to core ./halo_small

clean: ; rm -f pingpong halo_big halo_small decompose
.PHONY: all run clean
