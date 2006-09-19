SRCFILES := *.xml

target : manual.html manual.txt chunked/index.html pretty/index.html

chunked/index.html : $(SRCFILES) custom.xsl
	xmlto -o chunked -m custom.xsl html manual.xml

pretty/index.html : $(SRCFILES) custom-pretty.xsl makeover
	xmlto -o pretty -m custom-pretty.xsl html manual.xml
	ls pretty/*.html | xargs -n 1 tidy -utf8 -m -i
	./makeover

manual.html : $(SRCFILES) custom-nochunks.xsl
	xmlto -m custom-nochunks.xsl html-nochunks manual.xml

manual.txt : manual.html
	html2text -nobs -style pretty manual.html > manual.txt

clean :
	-rm -rf pretty/*.html chunked/*.html manual.html manual.txt
