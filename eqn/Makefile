PNGLATEX=./pnglatex

.PHONY: all clean

all: complex_dft.png

clean:
	rm *.png

%.png: %.tex
	$(PNGLATEX) $*.tex
