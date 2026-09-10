# Build the manual. The option reference is generated from the source table in
# src/ADZE_pfile.cpp so it cannot document options the program does not have.

PDFLATEX ?= pdflatex
MANUAL   = ADZE_Manual

all: $(MANUAL).pdf

options.tex: ../src/ADZE_pfile.cpp gen_options.py
	./gen_options.py ../src/ADZE_pfile.cpp options.tex

$(MANUAL).pdf: $(MANUAL).tex options.tex $(wildcard figures/*.png)
	$(PDFLATEX) -interaction=nonstopmode $(MANUAL).tex >/dev/null
	$(PDFLATEX) -interaction=nonstopmode $(MANUAL).tex >/dev/null   # table of contents

install: $(MANUAL).pdf
	cp $(MANUAL).pdf ../$(MANUAL).pdf

clean:
	rm -f $(MANUAL).aux $(MANUAL).log $(MANUAL).out $(MANUAL).toc options.tex

.PHONY: all install clean
