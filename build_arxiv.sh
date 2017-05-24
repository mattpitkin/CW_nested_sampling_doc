figs=`find figures -name '*.pdf'`
caps=`find figures -name caption.tex`

tar -cvzf arXiv.tar.gz abstract.tex acknowledgements.tex appendix0*.tex emulateapj.cls main.tex main.bbl posttitle.tex section0*.tex $figs $caps
