
name=metro0

default: $(name).pdf

dvipdf: $(name).tex $(name).bib
	latex $(name).tex
	bibtex $(name).aux
	latex $(name).tex
	latex $(name).tex
	latex $(name).tex
	dvipdf $(name).dvi

$(name).pdf: $(name).tex $(name).bib
	pdflatex $(name).tex
	bibtex $(name).aux
	pdflatex $(name).tex
	pdflatex $(name).tex
	pdflatex $(name).tex

clean:
	rm -vf $(name).aux
	rm -vf $(name).log
	rm -vf $(name).dvi
	rm -vf $(name).pdf
	rm -vf $(name).out
	rm -vf $(name).bbl
	rm -vf $(name).blg
