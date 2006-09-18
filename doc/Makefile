target : html/manual.html manual.txt
#html/default.css

html/manual.html : manual.xml
	db2html -o html -d "custom.dsl#html" -u manual.xml
	cp default.css html/

#html/default.css : default.css
	#cp default.css html/

manual.txt : html/manual.html
	html2text -nobs -style pretty html/manual.html > manual.txt

clean :
	-rm -rf html manual manual.junk manual.txt
