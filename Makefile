# make file to automate all the annoying commands

default: pdf

pdf: *.tex bibliography/biblio.bib layout.md
	local_build.py --latex pdflatex --filename main.tex --build-dir .

ps: *.tex bibliography/biblio.bib layout.md
	local_build.py --latex latex --filename main.tex --build-dir .
	dvips main

clean:
	@echo "Cleaning directory of backups and logs"
	rm -f *~ *.log *.aux *.dvi *.out *.backup *.bbl *.blg *.ps *Notes.bib

