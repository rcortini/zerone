all: document.pdf ../tutorial.pdf

document.pdf: document.tex document.bib
	pdflatex document.tex
	bibtex document.aux
	pdflatex document.tex
	pdflatex document.tex

document.bib: pubmed2bibtex.xsl pmid.xml
	xsltproc pubmed2bibtex.xsl pmid.xml > $@

pmid.xml: getref.py pmid.txt
	./getref.py pmid.txt > $@

../tutorial.pdf: tutorial.tex
	pdflatex tutorial.tex
	mv tutorial.pdf ..
	#rm tutorial.aux tutorial.log tutorial.out
