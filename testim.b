#!/bin/csh -f

set pix = $1
set tot = $2
setenv IMAGIC_BATCH 1
echo "! "
echo "! "
echo "! ====================== "
echo "! IMAGIC ACCUMULATE FILE "
echo "! ====================== "
echo "! "
echo "! "
echo "! IMAGIC program: testim -----------------------------------------------"
echo "! "
/opt/qb3/imagic-070813/stand/testim.e <<EOF
test,1,$tot
$pix,$pix
REAL
BLOBS
EOF
