#!/bin/csh -f

# Usage :  
#  prompt> skdetsim.sh [card name] [output file] [input file] [# of sukipped event]
#

#
# skdetsim_high seems to have much memory size, so need to unlimit memory use in csh.
#
#limit memoryuse unlimited
#limit vmemoryuse unlimited

## assume skdetsim exists in the same directory
set DIR=`echo $0:h` 
if ("$DIR" == "$0") set DIR="."

## SUBRUN number is input from enviromment variable.
setenv MC_SUBNUM 0

echo $#argv

if ($#argv>0) then
    set CARD=$1
else
    set CARD=supersim.card
endif

if ($#argv>1) then
    set FNAME_OUT=$2
else
    set FNAME_OUT=dummy.output
endif

if ($#argv>2) then
    set FNAME_IN=$3
else
    set FNAME_IN=dummy.input
endif

if ($#argv>3) then
$DIR/skdetsim_high $CARD $FNAME_OUT $FNAME_IN $4
else
$DIR/skdetsim_high $CARD $FNAME_OUT $FNAME_IN
endif

