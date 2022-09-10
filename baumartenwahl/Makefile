PAPER := baumartenwahl
PICDIR := pic
XFIGDIR := xfig
RPLOTDIR := r
TABSDIR := tabs
MAPSDIR := maps
WEBPICDIR := webPic

FIG := $(wildcard $(XFIGDIR)/*.fig)
FIG_AS_PDF := $(subst $(XFIGDIR)/, $(PICDIR)/, $(FIG:.fig=.pdf))

RPLOTS := $(wildcard $(RPLOTDIR)/*.plot.r)
RPLOTS_AS_PDF := $(subst $(RPLOTDIR)/, $(PICDIR)/, $(RPLOTS:.plot.r=.pdf))

MAPS := $(wildcard $(MAPSDIR)/*Get.sh)
MAPS_AS_PDF := $(subst $(MAPSDIR)/, $(PICDIR)/, $(MAPS:Get.sh=.pdf))

WEBPIC := $(wildcard $(WEBPICDIR)/*Get.sh)
WEBPIC_AS_JPG := $(subst $(WEBPICDIR)/, $(PICDIR)/, $(WEBPIC:Get.sh=.jpg))

PICTURES := $(FIG_AS_PDF) $(RPLOTS_AS_PDF) $(MAPS_AS_PDF) $(WEBPIC_AS_JPG)

TABLES := $(wildcard $(TABSDIR)/*.csv)

.PHONY: all purge clean rmDownloads

all: $(PAPER).pdf

$(PAPER).pdf: $(PICTURES) $(TABLES) $(wildcard *.tex) $(PAPER).bbl
	pdflatex $(PAPER)
	LATEX_COUNT=0; \
	while grep "Warning: Citation" $(PAPER).log \
	  && [ $$LATEX_COUNT -lt 5 ] ; do      \
	  LATEX_COUNT=`expr $$LATEX_COUNT + 1` ; \
	  bibtex $(PAPER); pdflatex $(PAPER) ;                    \
	done
#	if grep "Citation.*undefined" $(PAPER).log; then bibtex $(PAPER); pdflatex $(PAPER); bibtex $(PAPER); fi
	if grep "Writing index file" $(PAPER).log; then makeindex $(PAPER); pdflatex $(PAPER) ; fi
	LATEX_COUNT=0; \
	while grep "Rerun \(to get\|LaTeX\)\|undefined references\|Linenumber reference failed" $(PAPER).log  \
	   && [ $$LATEX_COUNT -lt 5 ] ; do      \
	  LATEX_COUNT=`expr $$LATEX_COUNT + 1` ; \
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

$(PICDIR)/%.pdf : $(MAPSDIR)/%Get.sh
	( cd $(MAPSDIR) ; bash $(<:$(MAPSDIR)/%=%) ../$@ )

$(PICDIR)/%.jpg : $(WEBPICDIR)/%Get.sh
	( cd $(WEBPICDIR) ; bash $(<:$(WEBPICDIR)/%=%) ../$@ )

purge:
	$(RM) $(PAPER).aux $(PAPER).bbl $(PAPER).blg $(PAPER).log $(PAPER).out $(PAPER).toc $(PAPER).idx $(PAPER).ilg $(PAPER).ind
	$(RM) $(PICTURES)
#	if [ -d "$(PICDIR)" ] && [ `ls -A "$(PICDIR)" | wc -l` -eq 0 ]; then rmdir $(PICDIR); fi
	if [ -d "$(PICDIR)" ] && [ ! -f "$(PICDIR)/"* ]; then rmdir $(PICDIR); fi

clean: purge
	$(RM) $(PAPER).pdf
	$(RM) $(PAPER).tar.xz

rmDownloads:
	$(RM) $(MAPSDIR)/waelderDesOstalpenraumesORIG.jpg
	$(RM) $(MAPSDIR)/gadm40_AUT_1.kmz
	$(RM) $(WEBPICDIR)/fichtenHerkuenfteORIG.jpg
	$(RM) $(WEBPICDIR)/tannerVerbreitungskarteORIG.png

tar: $(PAPER).tar.xz

$(PAPER).tar.xz: $(shell find * -type f ! -name $(PAPER).tar.xz)
	touch $(PAPER).tar.xz
	tar --exclude="$(PAPER).tar.xz" -cJf $(PAPER).tar.xz .
