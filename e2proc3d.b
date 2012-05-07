#!/bin/csh -f

set file = $1

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
IMAGIC
SPI
SINGLE_FILE
3
${file:r}
${file:r}.spi
LINUX
NO
EOF
