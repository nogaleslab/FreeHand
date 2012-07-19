#!/bin/csh

set l = $1

foreach file ($l*)
	
	tail -n +4 "$file"

end

