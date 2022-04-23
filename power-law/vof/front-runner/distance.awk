BEGIN {
    FS=" |:"
}
{
    if ($1 == "#") {
        # get column indices of the relevant fields (X,Y)
	for (i = 1; i <= NF; i++) {
	    if ($i == "T")
		iT = $(i-1);
	    else if ($i == "X")
		iX = $(i-1);
	    else if ($i == "Y")
		iY = $(i-1);
	}
    }
    else {
	T = $iT
	if (T > 0.011 && T < 0.989) {
            # the cell contains an interface
	    x = $iX
	    y = $iY
	    print $1,y;
	    fflush (stdout);
	}
    }
} 
