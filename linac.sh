#!/bin/csh -f

# Usage :  
#     imake_boot (after edit Imakefile and modify like Imakefile_linac)
#  prompt> linac.sh [card name] [output file]  [run_num] [input file] [# of sukipped event]
#
 
setenv RFLIST  /tmp/rflist.$$.`hostname`

## SUBRUN number is input from enviromment variable.

ln -s /cern/95a/lib/xsneut95.dat xsneut.dat
ln -s /home/atmpd/const/chetc.dat .

setenv MC_SUBNUM 0

# echo $#argv
if ($#argv>3) then
cat <<! >! $RFLIST
10{{"$4",LOCAL,,RED,,,"recl=5670  form=formatted",}}
20{{"$2",LOCAL,,RPL,,,"recl=5670  status=new"}}
!

else
cat <<! >! $RFLIST
20{{"$2",LOCAL,,RPL,,,"recl=5670 status=new"}}
!

endif

if ($#argv>4) then
./solaris_sparc/skdetsim $1 $5 $3
else
./solaris_sparc/skdetsim $1 0 $3
endif

/bin/rm $RFLIST
/bin/rm xsneut.dat
/bin/rm chetc.dat



