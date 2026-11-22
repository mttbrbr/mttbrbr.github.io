CC      ?= gcc
CFLAGS  = -O3 -march=native -std=c11 -Wall
LDLIBS  = -lm
SIZES   = 64 128 256 512 724 1024 1448 2048 2896

all: bin_soa_ij bin_soa_ji bin_aos_ij bin_aos_ji bin_soa_ij_pad

bin_soa_ij: locality.c ; $(CC) $(CFLAGS) $< -o $@ $(LDLIBS)
bin_soa_ji: locality.c ; $(CC) $(CFLAGS) -DORDER_JI $< -o $@ $(LDLIBS)
bin_aos_ij: locality.c ; $(CC) $(CFLAGS) -DLAYOUT_AOS $< -o $@ $(LDLIBS)
bin_aos_ji: locality.c ; $(CC) $(CFLAGS) -DLAYOUT_AOS -DORDER_JI $< -o $@ $(LDLIBS)
bin_soa_ij_pad: locality.c ; $(CC) $(CFLAGS) -DPAD $< -o $@ $(LDLIBS)

# Sweep grid size for every variant; output is CSV
bench: all
	@echo "layout,order,N,working_set_MB,steps,seconds,ns_per_cell_update,checksum_u,checksum_p"
	@for n in $(SIZES); do for b in bin_soa_ij bin_aos_ij bin_soa_ji bin_aos_ji bin_soa_ij_pad; do ./$$b $$n; done; done

clean: ; rm -f bin_*
.PHONY: all bench clean
