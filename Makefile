# make file to automate all the annoying commands

default: pdf

pdf: nested_sampling_doc.tex nested_sampling_doc.bib
	pdflatex nested_sampling_doc
	bibtex nested_sampling_doc
	pdflatex nested_sampling_doc
	pdflatex nested_sampling_doc

ps: nested_sampling_doc.tex nested_sampling_doc.bib
	latex nested_sampling_doc
	bibtex nested_sampling_doc
	latex nested_sampling_doc
	latex nested_sampling_doc
	dvips nested_sampling_doc
	ps2pdf nested_sampling_doc.ps

clean:
	@echo "Cleaning directory of backups and logs"
	rm -f *~ *.log *.aux *.dvi *.out *.backup *.bbl *.blg *.ps *Notes.bib

