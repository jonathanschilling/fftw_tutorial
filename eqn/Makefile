PNGLATEX=./pnglatex

.PHONY: all clean

all: complex_dft.png \
     c2c_product.png \
     euler.png phi.png \
     c2c_product_real.png \
     hermitian.png \
     hermitian_sum.png \
     redft00.png \
     redft10.png \
     redft01.png \
     redft11.png \
     rodft00.png \
     rodft10.png \
     rodft01.png \
     rodft11.png \
     sin_cos_add.png \
     true_2d_idct.png

clean:
	rm *.png

%.png: %.tex
	$(PNGLATEX) $*.tex
