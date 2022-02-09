PAPER := mischOderReinbestand
PICDIR := pic
XFIGDIR := xfig
RPLOTDIR := r

FIG := $(wildcard $(XFIGDIR)/*.fig)
FIG_AS_PDF := $(subst $(XFIGDIR)/, $(PICDIR)/, $(FIG:.fig=.pdf))

RPLOTS := $(wildcard $(RPLOTDIR)/*.plot.r)
RPLOTS_AS_PDF := $(subst $(RPLOTDIR)/, $(PICDIR)/, $(RPLOTS:.plot.r=.pdf))

PICTURES := $(FIG_AS_PDF) $(RPLOTS_AS_PDF)

.PHONY: all purge clean

all: $(PAPER).pdf

$(PAPER).pdf: $(PICTURES) $(wildcard *.tex) $(PAPER).bbl
	pdflatex $(PAPER)
	if grep "Citation.*undefined" $(PAPER).log; then bibtex $(PAPER); pdflatex $(PAPER); fi
	if grep "Writing index file" $(PAPER).log; then makeindex $(PAPER); pdflatex $(PAPER) ; fi
	LATEX_COUNT=0
	while grep "Rerun (to get|LaTeX)|undefined references|Linenumber reference failed" $(PAPER).log  \
	   && [ $LATEX_COUNT -lt 5 ] ; do      \
	  LATEX_COUNT=`expr $LATEX_COUNT + 1` ; \
	  pdflatex $(PAPER) ;                    \
	done

$(PAPER).bbl: literature.bib
	if [ ! -f "$(PAPER).aux" ] ; then pdflatex $(PAPER) ; fi
	bibtex $(PAPER)

#$(PAPER).ind: $(PAPER).idx
#	makeindex $(PAPER)

$(PICTURES): | $(PICDIR)

$(PICDIR):
	mkdir -p $(PICDIR)

$(PICDIR)/%.pdf : $(XFIGDIR)/%.fig
	fig2dev -L pdf $< $@

$(PICDIR)/%.pdf : $(RPLOTDIR)/%.plot.r
	Rscript --vanilla $< PDFNAME=$@

purge:
	$(RM) $(PAPER).aux $(PAPER).bbl $(PAPER).blg $(PAPER).log $(PAPER).out $(PAPER).toc $(PAPER).idx $(PAPER).ilg $(PAPER).ind
	$(RM) $(PICTURES)
#	if [ -d "$(PICDIR)" ] && [ `ls -A "$(PICDIR)" | wc -l` -eq 0 ]; then rmdir $(PICDIR); fi
	if [ -d "$(PICDIR)" ] && [ ! -f "$(PICDIR)/"* ]; then rmdir $(PICDIR); fi

clean: purge
	$(RM) $(PAPER).pdf
	$(RM) $(PAPER).tar.xz

tar:
	tar --exclude="$(PAPER).tar.xz" -cJf $(PAPER).tar.xz .
