#!/bin/csh -f

# Usage :  
#  prompt> skdetsim.sh [card name] [output file] [input file] [# of sukipped event]
#

## assume skdetsim exists in the same directory

set DIR=`echo $0:h` 
if ("$DIR" == "$0") set DIR="."

## SUBRUN number is input from enviromment variable.
## This should be editted by user.

setenv MC_SUBNUM 0

## Random seed file for NEUT.
## This should be controlled by user.

#setenv RANFILE /home/atmpd/skdetsim/skdetsim-v13p90_mar16/random.tbl.000
setenv RANFILE ./random.tbl.000

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






