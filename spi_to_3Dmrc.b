#!/bin/csh -f

set input=$1

setenv IMAGIC_BATCH 1
echo "! "
echo "! "
echo "! ====================== "
echo "! IMAGIC ACCUMULATE FILE "
echo "! ====================== "
echo "! "
echo "! "
echo "! IMAGIC program: em2em ------------------------------------------------"
echo "! "
/opt/qb3/imagic-070813/stand/em2em.e <<EOF
SPI
SINGLE_FILE
IMAGIC
2D
${input:r}
spi
${input:r}
1,1
EOF
echo "! "
echo "! IMAGIC program: headers ----------------------------------------------"
echo "! "
/opt/qb3/imagic-070813/stand/headers.e <<EOF
${input:r}
WRITE
INDEX
61
IN
0
EOF
echo "! "
echo "! IMAGIC program: em2em ------------------------------------------------"
echo "! "
/opt/qb3/imagic-070813/stand/em2em.e <<EOF
IM
MRC
3
${input:r}
${input:r}.mrc
NO
EOF
