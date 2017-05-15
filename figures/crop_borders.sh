#!/usr/bin/bash

# crop whitespace borders around all figures
pdffiles=`find . -name "*.pdf"`

for f in $pdffiles; do y=${f%.*}; z=${y#*/}; echo ${z}; pdfcrop $f ${z}_trimmed.pdf; pdftops -eps ${z}_trimmed.pdf; mv -f ${z}_trimmed.pdf $f; mv -f ${z}_trimmed.eps ${z}.eps; done;
