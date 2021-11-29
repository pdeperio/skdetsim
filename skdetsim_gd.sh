#!/bin/bash
#
#   sample shell script to run skdetsim_gd (SK-Gd skdetsim with Geant4 ext.)
#

# Usage :  
#  prompt> skdetsim_gd.sh [card name] [output file] [input file] [# of sukipped event]
#

## assume skdetsim exists in the same directory

#DIR=.
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"


source /usr/local/sklib_gcc4.8.5/geant4.10.03.p01/bin/geant4.sh
export G4NEUTRONHPDATA="/usr/local/sklib_gcc4.8.5/geant4.10.03.p01/share/Geant4-10.3.1/data/G4NDL4.5"

# set SKDETSIMDIR to skdetsim directory for "gdcap_gamma.dat" and "Gd155(7).dat"
export SKDETSIMDIR=${DIR}



## SUBRUN number is input from enviromment variable.
export MC_SUBNUM=${MC_SUBNUM-0}


if [ $# -gt 0 ] ; then
    CARD=$1
else
    CARD=sk6_example.card
fi

if [ $# -gt 1 ] ; then
    FNAME_OUT=$2
else
    FNAME_OUT=dummy.output
fi

if [ $# -gt 2 ] ; then
    FNAME_IN=$3
else
    FNAME_IN=dummy.input
fi

if [ $# -gt 3 ] ; then
   $DIR/skdetsim $CARD $FNAME_OUT $FNAME_IN $4
else
   $DIR/skdetsim $CARD $FNAME_OUT $FNAME_IN
fi










