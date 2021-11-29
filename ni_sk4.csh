#!/bin/csh -f

# Usage :  
#  prompt> linac_sk4.sh [output file] [posx] [posy] [posz]
#

## assume skdetsim exists in the same directory
set DIR=`echo $0:h` 
if ("$DIR" == "$0") set DIR="."

set CARD="ni_sk4_example_zbs.card"

if ($#argv>1) then
    set FNAME_OUT=$1
else
    echo 'Usage : linac_sk4.sh [output file] [run number]'
    exit(0)
endif

$DIR/skdetsim $CARD $FNAME_OUT $2 $3 $4
