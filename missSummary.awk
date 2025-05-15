#!/bin/awk -f

BEGIN {
	print "RECORD\tPOSITION\tPROPORTION";
}
{
	if (substr($0, 0, 1) == "#") {
		numRecords = 1;
		next;
	} else {
		total = 0;
		for (i = 9; i <= NF; i++) {
			if (substr($i, 1, 1) == "." && substr($i, 3, 1) == ".") {
				total += 1;
			}
		}
		print numRecords"\t"$2"\t"(total / (NF - 9));
		numRecords += 1; 
	}
}
