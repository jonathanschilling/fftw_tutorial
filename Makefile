
CC=gcc

LIBS=-lfftw3 -lm

.PHONY: all clean

all: test_1d_c2c \
     test_1d_c2r \
     test_1d_r2c \
     test_1d_redft00 \
     test_1d_redft10

clean:
	rm test_1d_c2c

test_1d_c2c: src/test_1d_c2c.c src/util.h
	${CC} src/test_1d_c2c.c -o test_1d_c2c ${LIBS}

test_1d_c2r: src/test_1d_c2r.c src/util.h
	${CC} src/test_1d_c2r.c -o test_1d_c2r ${LIBS}

test_1d_r2c: src/test_1d_r2c.c src/util.h
	${CC} src/test_1d_r2c.c -o test_1d_r2c ${LIBS}

test_1d_redft00: src/test_1d_redft00.c src/util.h
	${CC} src/test_1d_redft00.c -o test_1d_redft00 ${LIBS}

test_1d_redft10: src/test_1d_redft10.c src/util.h
	${CC} src/test_1d_redft10.c -o test_1d_redft10 ${LIBS}


