target : html/manual.html manual.txt
#html/default.css

html/manual.html : manual.xml pairingfns.xml elementfns.xml basics.xml
	#db2html -o html -d "custom.dsl#html" manual.xml
	xmlto -o html -m custom.xsl html manual.xml
	xmlto -o html -m custom-nochunks.xsl html-nochunks manual.xml

#html/default.css : default.css
	#cp default.css html/

manual.txt : html/manual.html
	html2text -nobs -style pretty html/manual.html > manual.txt

clean :
	-rm -rf html/*.html manual.txt
