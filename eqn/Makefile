PNGLATEX=./pnglatex

.PHONY: all clean

all: complex_dft.png c2c_product.png euler.png phi.png

clean:
	rm *.png

%.png: %.tex
	$(PNGLATEX) $*.tex
