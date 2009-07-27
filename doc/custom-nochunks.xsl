<?xml version='1.0'?>
<xsl:stylesheet xmlns:xsl="http://www.w3.org/1999/XSL/Transform"
                xmlns:fo="http://www.w3.org/1999/XSL/Format"
		                version="1.0">
<xsl:param name="html.stylesheet" select="'default.css'"/>
<xsl:param name="generate.toc" select="'book toc'"/>
<xsl:output method="html" encoding="UTF-8" indent="no"
doctype-public="-//W3C//DTD HTML 4.01 Transitional//EN"
/>
<xsl:template name="user.footer.navigation">
<script type="text/javascript">
var gaJsHost = (("https:" == document.location.protocol) ? "https://ssl." : "http://www.");
document.write(unescape("%3Cscript src='" + gaJsHost + "google-analytics.com/ga.js' type='text/javascript'%3E%3C/script%3E"));
</script>
<script type="text/javascript">
try{
var pageTracker = _gat._getTracker("UA-1901330-5");
pageTracker._trackPageview();
} catch(err) {}
</script>
</xsl:template>
</xsl:stylesheet>
