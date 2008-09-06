.PHONY: target clean sync gendoc

target: manual.html

TXTFILES:=preface.txt quickstart.txt basics.txt sigex.txt \
  pairingfns.txt elementfns.txt paramfns.txt randomfns.txt \
  auxfns.txt alloc.txt bundle.txt internal.txt security.txt \
  contributors.txt

GENFILES:=gen/*.txt

gendoc $(GENFILES) : ../include/*.h extract
	-rm $(GENFILES)
	-mkdir gen
	cat ../include/*.h | ./extract
	for a in gen/*.*.txt; do b=$${a%.*.txt}.txt; cat $$a $$b > tmp; mv tmp $$b ; rm $$a; done

manual.xml: $(TXTFILES) $(GENFILES)
	( for FILE in $(TXTFILES) ; do cat $$FILE ; echo ; done ) | asciidoc -f macros.ad -d book -b docbook - > $@

chunked/index.html : manual.xml custom.xsl
	xmlto -o chunked -m custom.xsl html manual.xml

manual: manual.xml custom-pretty.xsl pretty.css
	xmlto -m custom-pretty.xsl -o manual html manual.xml
	sed -i 's/xmlns:fo[^ ]*//g' manual/*.html
	ls manual/*.html | xargs -n 1 tidy -utf8 -m -i -q
	./makeover
	cp pretty.css manual/

manual.html : manual.xml custom-nochunks.xsl
	xmlto -m custom-nochunks.xsl html-nochunks manual.xml
	-tidy -utf8 -imq $@

manual.txt : manual.html
	html2text -nobs -style pretty manual.html > manual.txt

manual.pdf: manual.xml
	docbook2pdf manual.xml

clean:
	-rm -rf manual.xml manual.html manual

sync: manual/* manual.html manual.txt chunked/* manual.pdf
	#rsync -r $^ blynn@tl1.stanford.edu:pbc/doc/
