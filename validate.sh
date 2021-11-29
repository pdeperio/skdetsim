#!/bin/bash
#
#------------------------------------------------------------------
# Instructions:
# 1) Check out a fresh copy of skdetsim from svn and compile it.
# 2) Enter fresh skdetsim directory and run this script with the 
#    argument 'svn' to indicate that this is a fresh version of
#    the skdetsim code (the first parameter is the SK geometry):
#    % cd ~/skdetsim
#    % validate.sh sk[]/all svn
# 3) Copy the resultant subdirectory '~/skdetsim/validation/' and all
#    of its contents to the directory that contains your working 
#    copy of skdetsim:
#    % cp -r ~/skdetsim/validation ~/my_skdetsim/.
# 4) Run this script again, this time without the 'svn' argument:
#    % cd ~/my_skdetsim
#    % validate.sh sk[]/all
#------------------------------------------------------------------
# 
# This script runs skdetsim in "validate" mode.  This mode produces a
# text dump of the PMT info.  If one runs this from a fresh svn
# checkout, pass it the parameter "svn".
#
# One can then make changes to the code, and run it again without the 
# "svn" parameter.  This will produce a dump of the current PMT info, 
# and compare to the version generated earlier from the repository-based
# code.
#
# The goal is to let the coder know when the changes he has made are
# changing the simulation.
#
# Note that due to differences in libraries and compilers, the svn-based
# tube dump cannot be kept in the repository, or it will always show
# changes even from a freshly checked out skdetsim.  The user must
# generate his own reference file locally for this to work.
#
# Two parameters are allowed.
#
#   1) "skver" = "sk1", "sk2", "sk3", "sk4", or "all"
#       Indicates which geometry to use. This parameter is required.
#       If "all" is given, the script will run over all geometries
#       in order.
#
#   2) "svn"
#      Indicates that this is a reference run based on fresh svn code.  
#      If omitted, the script will assume you want to compare your
#      current executable with a previously generated reference file.
#
# All output files are deposited into the 'validation' subdirectory.
#====================================================================
# -M. Litos   Jan.21, 2009
# based on original validation code by A. Habig
#====================================================================

if [ $# -lt 1 -o $# -gt 2 ] || [ "$1" != "sk1" -a "$1" != "sk2" -a "$1" != "sk3" -a "$1" != "sk4" -a "$1" != "all" ] ||  [ $# -eq 2 -a "$2" != "svn" ]; then
    echo $"Usage:" >&2
    echo $"       validate.sh skver [svn]" >&2
    echo $"where 'skver' = 'sk1', 'sk2', 'sk3', 'sk4', or 'all' " >&2
    echo $"and 'svn' indicates that this is a fresh copy of skdetsim" >&2
    exit 1
fi

export RFLIST=validate.rflist

# make subdir, to keep the clutter down
if [ ! -d validation ]; then
    mkdir validation
    for (( i=1; i<=4; i++ ))
      do
      ln -s ../validate_le_sk$i.card validation/validate_le_sk$i.card
      ln -s ../validate_ap_sk$i.card validation/validate_ap_sk$i.card
    done
fi

cd validation

# see the end for the main function.  It's functionized so we can do things 
# a couple times, for atmpd and le

# function which does all the work
# arguments: $1=skver $2=[svn] $3=ap/le
doit() {

rm -f validation.dat validation.zbs  # clobber existing scratch 
                                     # validation output files

if [ "$2" == svn ]; then
    if [ "$1" == all ]; then
	for (( i=1; i<=4; i++ ))
	  do
	  echo $"Creating initial $3 sk$i validation data from clean svn build" >&2
	  if [ -f sk$i.$SVNTUBEDAT ]; then
	      echo $"sk$i.$SVNTUBEDAT already exists, replacing it" >&2
	      rm -f sk$i.$SVNTUBEDAT sk$i.$SVNZBS
	  fi
# do initial run to create reference tube file from clean svn build
	  echo $"Running GEANT to make $3 sk$i tube data file, output log is in sk$i.$SVNLOG" >&2
	  rm -f sk$i.$SVNLOG
	  ../skdetsim validate_$3_sk$i.card validation.zbs > sk$i.$SVNLOG 2>&1
	  mv validation.dat sk$i.$SVNTUBEDAT
	  mv validation.zbs sk$i.$SVNZBS
	done
    else
	echo $"Creating initial $3 $1 validation data from clean svn build" >&2
	if [ -f $1.$SVNTUBEDAT ]; then
	    echo $"$1.$SVNTUBEDAT already exists, replacing it" >&2
	    rm -f $1.$SVNTUBEDAT $1.$SVNZBS
	fi
# do initial run to create reference tube file from clean svn build
	echo $"Running GEANT to make $3 $1 tube data file, output log is in $1.$SVNLOG" >&2
	rm -f $1.$SVNLOG
	../skdetsim validate_$3_$1.card validation.zbs > $1.$SVNLOG 2>&1
	mv validation.dat $1.$SVNTUBEDAT
	mv validation.zbs $1.$SVNZBS
    fi
else    
    if [ "$1" == all ]; then
	for (( i=1; i<=4; i++ ))
	  do
# create tube file to compare to pre-existing reference tube file
# First, test for existance of reference tube file
	  if [ -f sk$i.$SVNTUBEDAT ]; then
# then run the skdetsim to create tube file
	      if [ -f sk$i.$TUBEDAT ]; then
		  echo $"sk$i.$TUBEDAT already exists, replacing it" >&2
		  rm -f sk$i.$TUBEDAT sk$i.$ZBS
	      fi
	      echo $"Running GEANT to make $3 sk$i tube data file, output log is in sk$i.$LOG" >&2
	      rm -f sk$i.$LOG
	      ../skdetsim validate_$3_sk$i.card validation.zbs > sk$i.$LOG 2>&1
	      mv validation.dat sk$i.$TUBEDAT
	      mv validation.zbs sk$i.$ZBS	      
# then check to see if things differ
	      rm -f validation_$3_sk$i_diff.txt
	      diff sk$i.$SVNTUBEDAT sk$i.$TUBEDAT > validation_$3_sk$i_diff.txt
	      if [ ! "$?" == 0 ]; then
		  echo 
		  echo $"Uh-oh, the $3 sk$i files are different, see validation_$3_sk$i_diff.txt for details" >&2
		  return 1
	      else
		  echo $"No $3 sk$i differences found!" >&2
		  rm -f validation_$3_sk$i_diff.txt
	      fi	      
	  else
	      echo $"no sk$i.$SVNTUBEDAT reference file exists" >&2
	      echo $"re-run this script with validate.sh svn" >&2
	      echo $"to create one.  " >&2
	      return 1
	  fi
	done
    else
# create tube file to compare to pre-existing reference tube file
# First, test for existance of reference tube file
	if [ -f $1.$SVNTUBEDAT ]; then
# then run the skdetsim to create tube file
	    if [ -f $1.$TUBEDAT ]; then
		echo $"$1.$TUBEDAT already exists, replacing it" >&2
		rm -f $1.$TUBEDAT $1.$ZBS
	    fi
	    echo $"Running GEANT to make $3 $1 tube data file, output log is in $1.$LOG" >&2
	    rm -f $1.$LOG
	    ../skdetsim validate_$3_$1.card validation.zbs > $1.$LOG 2>&1
	    mv validation.dat $1.$TUBEDAT
	    mv validation.zbs $1.$ZBS	    
# then check to see if things differ
	    rm -f validation_$3_$1_diff.txt
	    diff $1.$SVNTUBEDAT $1.$TUBEDAT > validation_$3_$1_diff.txt
	    if [ ! "$?" == 0 ]; then
		echo 
		echo $"Uh-oh, the $3 $1 files are different, see validation_$3_$1_diff.txt for details" >&2
		return 1
	    else
		echo $"No $3 $1 differences found!" >&2
		rm -f validation_$3_$1_diff.txt
	    fi	    
	else
	    echo $"no $1.$SVNTUBEDAT reference file exists" >&2
	    echo $"re-run this script with validate.sh svn" >&2
	    echo $"to create one.  " >&2
	    return 1
	fi
    fi
fi
return 0
}

# This is at the end so that the function above gets defined so it 
# can be called

# do the LE run
SVNTUBEDAT=validation_le_svn.dat
TUBEDAT=validation_le_current.dat
SVNLOG=validate_le_svn.log
LOG=validate_le_current.log
SVNZBS=validation_le_svn.zbs
ZBS=validation_le_current.zbs
doit "$1" "$2" le

# do the ATMPD run
SVNTUBEDAT=validation_ap_svn.dat
TUBEDAT=validation_ap_current.dat
SVNLOG=validate_ap_svn.log
LOG=validate_ap_current.log
SVNZBS=validation_ap_svn.zbs
ZBS=validation_ap_current.zbs
RANFILE=${NEUT_ROOT}/src/neutmspl/random.tbl # Required if using NEUT pion scattering model with fixed seed
doit "$1" "$2" ap
