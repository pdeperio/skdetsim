#!/bin/csh -f
#
#   csh script for compiling SK-Gd skdetsim with G4 extension (skdetsim_gd)
#   please use this script to compile.
#
#   Usage:
#     %  ./make_gd.csh clean
#     %  ./make_gd.csh 
#
#   Note:
#     Currently this works only on Kamioka sukap machine.
#     Geant4.10.03.p01
#


#
#  Check if you are using gcc4/gfortran
#
setenv LANG C
setenv LC_ALL C
gcc -v |& grep 4.8 > /dev/null
if ( $? != 0 ) then;
   echo " Please use gcc4.8 version "
   exit 1;
endif

#
#  enable G4 extension
#
setenv WithG4 1

#
# simulate neutron 
#
setenv NEUTRON 1


#
# if you want to enlarge memory size for atmpd use,
# use one of these flags depending on your use
#
#setenv MIDDLE 1
#setenv HIGH   1


set cwd = `pwd`
setenv G4WORKDIR "${cwd}/geant4_work/"
setenv G4TMPDIR "$G4WORKDIR/tmp/"



set g4top = "/usr/local/sklib_gcc4.8.5/geant4.10.03.p01"
set g4ver = `cd ${g4top}/share ;  ls -d Geant4-*`

setenv G4LEDATA "${g4top}/share/${g4ver}/data/G4EMLOW6.50"
setenv G4LEVELGAMMADATA "${g4top}/share/${g4ver}/data/PhotonEvaporation4.3.2"
setenv G4LIB "${g4top}/lib64/${g4ver}/"
setenv G4INCLUDE "${g4top}/include/Geant4"
setenv G4INSTALL "${g4top}/share/${g4ver}/geant4make"
setenv G4LIB_BUILD_SHARED "1"
setenv G4LIB_USE_ZLIB "1"
setenv G4NEUTRONHPDATA "${g4top}/share/${g4ver}/data/G4NDL4.5"
setenv G4NEUTRONXSDATA "${g4top}/share/${g4ver}/data/G4NEUTRONXS1.4"
setenv G4PIIDATA "${g4top}/share/${g4ver}/data/G4PII1.3"
setenv G4RADIOACTIVEDATA "${g4top}/share/${g4ver}/data/RadioactiveDecay5.1.1"
setenv G4REALSURFACEDATA "${g4top}/share/${g4ver}/data/RealSurface1.0"
setenv G4SAIDXSDATA "${g4top}/share/${g4ver}/G4SAIDDATA1.1"
setenv G4SYSTEM "Linux-g++"
setenv G4UI_USE_TCSH "1"
setenv G4VIS_USE_OPENGLX "1"

setenv CLHEP_BASE_DIR    /home/skofl/sklib_gcc4.8.5/clhep-2.3.4.3
setenv CLHEP_INCLUDE_DIR $CLHEP_BASE_DIR/include
setenv CLHEP_LIB_DIR     $CLHEP_BASE_DIR/lib

if ( ! ${?LD_LIBRARY_PATH} ) then
  setenv LD_LIBRARY_PATH ${CLHEP_LIB_DIR}
else
  setenv LD_LIBRARY_PATH ${CLHEP_LIB_DIR}:${LD_LIBRARY_PATH}
endif

setenv LIBRARY_PATH "$G4TMPDIR/$G4SYSTEM"


make $1

