#!/bin/csh

set f = $1

echo ${f}_plots_CC_v101_merge.mrc
totsumstack.exe << eot
${f}_plots_CC_v101_merge.mrc
${f}_averageplot_CC_v101.mrc
eot
