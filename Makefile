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

ps: *.tex bibliography/biblio.bib layout.md
	latex main.tex
	bibtex main
	latex main.tex
	latex main.tex
	dvips main
	ps2pdf main.ps

compressed: *.tex bibliography/biblio.bib
	pdflatex '\def\usepdffigs{true}\input{main.tex}'
	bibtex main
	pdflatex '\def\usepdffigs{true}\input{main.tex}'
	pdflatex '\def\usepdffigs{true}\input{main.tex}'
	pdf2ps main.pdf main.ps
	ps2pdf -dPDFSETTINGS=/screen -dDownsampleColorImages=true -dColorImageResolution=200 -dColorImageDownsampleType=/Bicubic main.ps maincomp.pdf

clean:
	@echo "Cleaning directory of backups and logs"
	rm -f *~ *.log *.aux *.dvi *.out *.backup *.bbl *.blg *.ps *Notes.bib

