<!DOCTYPE style-sheet PUBLIC "-//James Clark//DTD DSSSL Style Sheet//EN" [
<!ENTITY html-ss
  PUBLIC "-//Norman Walsh//DOCUMENT DocBook HTML Stylesheet//EN" CDATA dsssl>
<!ENTITY print-ss
  PUBLIC "-//Norman Walsh//DOCUMENT DocBook Print Stylesheet//EN" CDATA dsssl>
]>

<style-sheet>
<style-specification id="print" use="print-stylesheet">
<style-specification-body> 

;; don't remember where I found the original version of this file (Ben)

;; customize the print stylesheet

</style-specification-body>
</style-specification>
<style-specification id="html" use="html-stylesheet">
<style-specification-body> 

;; customize the html stylesheet

(define %stylesheet% "default.css")
(define (chunk-section-depth) 0)
(define %stylesheet-type% "text/css")

</style-specification-body>
</style-specification>
<external-specification id="print-stylesheet" document="print-ss">
<external-specification id="html-stylesheet" document="html-ss">
</style-sheet>
