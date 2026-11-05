# Build the manual. The option reference is generated from the source table in
# src/ADZE_pfile.cpp so it cannot document options the program does not have.

PDFLATEX ?= pdflatex
MANUAL   = ADZE_Manual

all: $(MANUAL).pdf

options.tex: ../src/ADZE_pfile.cpp gen_options.py
	./gen_options.py ../src/ADZE_pfile.cpp options.tex

# Every other .tex file here is a chapter the manual \input{}s (options.tex is
# generated, so it has its own rule). Listing them by wildcard means a new
# chapter is a dependency without editing this file.
CHAPTERS = $(filter-out $(MANUAL).tex options.tex,$(wildcard *.tex))

# LaTeX resolves cross-references and the table of contents from the previous
# pass, so a fixed number of passes can leave them stale: after a clean build,
# two were not always enough. Run until the log no longer asks for a rerun, at
# most five times. Then fail if any reference is still undefined, so that a
# \ref to a removed section stops the build instead of printing "??".
$(MANUAL).pdf: $(MANUAL).tex $(CHAPTERS) options.tex $(wildcard figures/*.png)
	@n=0; while :; do \
	  n=$$((n+1)); \
	  $(PDFLATEX) -interaction=nonstopmode -halt-on-error $(MANUAL).tex >/dev/null \
	    || { tail -20 $(MANUAL).log; exit 1; }; \
	  grep -q 'Rerun to get' $(MANUAL).log || break; \
	  if [ $$n -ge 5 ]; then echo "manual: references still changing after $$n passes"; exit 1; fi; \
	done; \
	echo "manual: $$n passes"
	@if grep -E "Reference .* undefined|Citation .* undefined" $(MANUAL).log; then \
	  echo "manual: undefined references (see above)"; rm -f $(MANUAL).pdf; exit 1; fi

install: $(MANUAL).pdf
	cp $(MANUAL).pdf ../$(MANUAL).pdf

clean:
	rm -f $(MANUAL).aux $(MANUAL).log $(MANUAL).out $(MANUAL).toc options.tex

.PHONY: all install clean
