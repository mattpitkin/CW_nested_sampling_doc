default: pdf

pdf: *.tex bibliography/biblio.bib
	pdflatex '\def\usepdffigs{true}\input{main.tex}'
	bibtex main
	pdflatex '\def\usepdffigs{true}\input{main.tex}'
	pdflatex '\def\usepdffigs{true}\input{main.tex}'

# compile with .png versions of figures 
png: *.tex bibliography/biblio.bib
	pdflatex main.tex
	bibtex main
	pdflatex main.tex
	pdflatex main.tex

# create a postscript verions
ps: *.tex bibliography/biblio.bib
	latex main.tex
	bibtex main
	latex main.tex
	latex main.tex
	dvips main

# use the postscript version to create a compressed pdf version
compressed: ps
	ps2pdf main.ps

clean:
	@echo "Cleaning directory of backups and logs"
	rm -f *~ *.log *.aux *.dvi *.out *.backup *.bbl *.blg *.ps *Notes.bib

