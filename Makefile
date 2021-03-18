
CC=gcc

LIBS=-lfftw3 -lm

.PHONY: all clean

all: test_1d_c2c \
     test_1d_c2r

clean:
	rm test_1d_c2c

test_1d_c2c: src/test_1d_c2c.c src/util.h
	${CC} src/test_1d_c2c.c -o test_1d_c2c ${LIBS}

test_1d_c2r: src/test_1d_c2r.c src/util.h
	${CC} src/test_1d_c2r.c -o test_1d_c2r ${LIBS}

