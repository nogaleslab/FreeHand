#!/bin/csh -f

set input = $1

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
IM
2D
${input:r}
spi
${input:r}
1,1
EOF
