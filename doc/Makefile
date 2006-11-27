# Makefile for generating manual in various formats from DocBook sources
.PHONY: target genxml

target : pretty/index.html

GENFILES:=gen/*.xml
SRCFILES:=*.xml $(GENFILES)

genxml $(GENFILES) : ../include/*.h extract
	-rm $(GENFILES)
	-mkdir gen
	cat ../include/*.h | ./extract
	for a in gen/*.*.xml; do b=$${a%.*.xml}.xml; cat $$a $$b > tmp; mv tmp $$b ; rm $$a; done

chunked/index.html : $(SRCFILES) custom.xsl
	xmlto -o chunked -m custom.xsl html manual.xml

pretty/index.html : $(SRCFILES) custom-pretty.xsl makeover index.html
	xmlto -o pretty -m custom-pretty.xsl html manual.xml
	-ls pretty/*.html | xargs -n 1 tidy -utf8 -m -i -q
	./makeover

manual.html : $(SRCFILES) custom-nochunks.xsl
	xmlto -m custom-nochunks.xsl html-nochunks manual.xml
	-tidy -utf8 -imq manual.html

manual.txt : manual.html
	html2text -nobs -style pretty manual.html > manual.txt

manual.pdf : $(SRCFILES)
	docbook2pdf manual.xml

pdf : $(SRCFILES)
	docbook2pdf manual.xml

clean :
	-rm -rf pretty/*.html chunked manual.html manual.txt gen forsite.tar.bz2 manual.pdf manual.tex

forsite: pretty/* manual.html manual.txt chunked/* manual.pdf
	tar cvfj forsite.tar.bz2 $^
