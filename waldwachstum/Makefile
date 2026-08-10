DOC = waldwachstum
SCRIPT_PLOTS = pic/plots.jl
SCRIPT_STAMM = pic/stammanalyse.jl
BIB_FILE     = literature.bib

.PHONY: all clean scripts

all: $(DOC).pdf

$(SCRIPT_PLOTS):
	cd pic && julia plots.jl

$(SCRIPT_STAMM):
	cd pic && julia stammanalyse.jl

$(DOC).pdf: $(DOC).tex $(SCRIPT_PLOTS) $(SCRIPT_STAMM) $(BIB_FILE)
	pdflatex -shell-escape $(DOC).tex
	biber $(DOC)
	-texindy -L german-duden $(DOC).idx
	pdflatex -shell-escape $(DOC).tex

clean:
	rm -f *.aux *.bbl *.bcf *.blg *.idx *.ilg *.ind *.log *.out *.run.xml *.toc *.synctex.gz
