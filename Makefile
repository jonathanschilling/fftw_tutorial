
CC=gcc

LIBS=-lfftw3 -lm

.PHONY: all clean test

TARGETS=test_util \
        test_1d_c2c \
        test_1d_c2r \
        test_1d_r2c \
        test_1d_redft00 \
        test_1d_redft10 \
        test_1d_redft01 \
        test_1d_redft11 \
        test_1d_rodft00 \
        test_1d_rodft10 \
        test_1d_rodft01 \
        test_1d_rodft11 \
        test_2d_c2c \
        test_2d_c2r \
        test_2d_r2c \
        test_2d_r2r_e10_o10 \
        test_2d_r2r_true2d \
        app_magn_axis \
        app_magn_axis_stellsym \
        app_flux_surface

all: $(TARGETS)

clean:
	rm -f test_*
	rm -f app_*
	rm -f axis_R*.dat
	rm -f lcfs_*.dat

test: $(addsuffix _test,$(TARGETS))

%_test: %
	./$*

$(TARGETS): % : src/%.c src/util.h
	$(CC) src/$@.c -o $@ $(LIBS)
