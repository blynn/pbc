<?xml version='1.0'?>
<xsl:stylesheet xmlns:xsl="http://www.w3.org/1999/XSL/Transform"
                xmlns:fo="http://www.w3.org/1999/XSL/Format"
		                version="1.0">
<!--
To chunk by chapter only:
<xsl:param name="chunk.section.depth" select="0"></xsl:param>
-->
<xsl:param name="chunker.output.encoding" select="'UTF-8'"></xsl:param>
<xsl:param name="chunker.output.doctype-public" select="'-//W3C//DTD HTML 4.01 Transitional//EN'"></xsl:param>

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
