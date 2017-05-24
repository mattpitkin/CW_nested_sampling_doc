for f in $pdffiles; do echo $f; convert -density 150 -flatten -trim $f ${f}_new.png; convert -density 150 ${f}_new.png ${f}_new.pdf; rm -f ${f}_new.png; mv -f ${f}_new.pdf ${f}; done
