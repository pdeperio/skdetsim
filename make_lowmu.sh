#!/bin/csh -f

# Usage :  
#  prompt> apdetsim.sh [card name] [output file] [input file] [# of sukipped event]
#
 
setenv RFLIST  /tmp/rflist.$$.`hostname`

## SUBRUN number is input from enviromment variable.

ln -s /cern/95a/lib/xsneut95.dat xsneut.dat
ln -s /home/atmpd/const/chetc.dat .

setenv MC_SUBNUM $1

setenv NUMBER $1

cat <<! >! $RFLIST
20{{"lowmumc$1.dat",LOCAL,,APP,,,"recl=5670 status=new"}}
!

./solaris_sparc/make_lowmu supersim.card $NUMBER

/bin/rm $RFLIST
/bin/rm xsneut.dat
/bin/rm chetc.dat



