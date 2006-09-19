// From my own website(!)
//TODO: only do this for links in the table of contents menu

function find_selflink() {
    var a = document.links;
    var i = 0;
    while (i < a.length) {
	if (a[i].href == document.URL) {
	    a[i].className = "currentlink";
	    /*
	    var s_new = document.createElement("span");
	    s_new.className = "currentlink";
	    s_new.appendChild(a[i]);
	    a[i].parentNode.replaceChild(s_new, a[i]);
	    */
	}
	i++;
    }
}

find_selflink();
