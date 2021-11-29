#
# In case of building locally, set SKOFL_ROOT variable 
#      setenv SKOFL_ROOT   ... 
#  or directly set in this makefile 
# SKOFL_ROOT = /skofl
#

#
# if you want run sk1 MC in lowe tune,
# please set the following parameter;
#
# SK1LE := true
#

ifdef WithG4
G4OPTIMISE := true

name := gtneut
G4TARGET := $(name)
G4EXLIB := true
G4_GTNEUT := $(G4TMPDIR)/$(G4SYSTEM)/gtneut/libgtneut.so
G4_LIBDIR := $(G4TMPDIR)/$(G4SYSTEM)

CCFLAGS   := -O3 -march=native
endif

ifndef SKOFL_ROOT
  SKOFL_ROOT = ../..
endif

FC=g77

.PHONY: all
ifdef WithG4
all: libg4 bin
endif

# ROOT support
CPPFLAGS += -I$(shell root-config --incdir) 
EXTRALIBS = $(shell root-config --glibs)
ROOTLIBDIR = $(shell root-config --libdir)

ifdef WithG4
CPPFLAGS += -I./include/ANNRI_Gd #-I./include/myReaction
# define for G4 compilation
include ./architecture.gmk
include ./binmake.gmk
endif

include $(SKOFL_ROOT)/config.gmk
# check if UseGFORTRAN flag is ON
#ifndef UseGFORTRAN
#   TMP := $(shell echo " ERROR: make sure you are using gcc4/gfortran version SKOFL library" 1>&2 )
#   TMP := $(shell echo "        SKOFL_ROOT=$(SKOFL_ROOT) " 1>&2 )
#   THIS_MUST_CRASH
#endif

LOCAL_INC	= -g -O0 -I$(CERN_ROOT)/include/geant321 

ifdef NEUT_ROOT
LOCAL_INC += -I$(NEUT_ROOT)/include
endif

LOCAL_LIBS	= $(INITOBJ) $(GOBJ) $(DOBJ) -L. -l $(LIBNAME) -lwtlib_5.1 $(CERN_ROOT)/lib/libgeant321.a\
							-L$(SKOFL_ROOT)/lib -lConnectionTableReader

LIBALIAS = skdetsim
APVER_MAJOR = 0
APVER_MINOR = 0
LEVER       = 14
VER = $(LEVER).$(APVER_MAJOR)$(APVER_MINOR)
EXEVER = $(LEVER).$(APVER_MAJOR)$(APVER_MINOR)

LIBNAME = $(LIBALIAS)_$(VER)

FORTRANDEFINES += -DCERNLIB_TYPE -DDETAIL_2NDARY -DREDUCED_PHOTON -DTUNED_1PE -DSKDETSIM_VERSION=$(VER) -DAPDETSIM_VERSION=$(APVER_MAJOR).$(APVER_MINOR) # -DFADC
#FORTRANDEFINES += -DDSSWTRIGSTOREINFO
#FORTRANDEFINES += -DUPMU
ifdef WithG4
FORTRANDEFINES += -DWithG4
endif

#
#  Objects
#


INITOBJ = sginit.o

SGO =	sgffin.o sggeom.o sgcard.o sgabsc.o sgpmt.o  sgpmta.o sgpmta_old.o \
        sgbst.o  sggvcd.o sgabsc.o sgpich.o sgpiel.o sgpiie.o \
        sgpint.o sgpist.o sgtyvk.o sgtyvk_odsegm.o sgtyvk_util.o \
	sgrefpmt.o skdonuts.o skdonuts_sk1.o sggvus.o sggvus_escale.o\
        sggtxt.o sggeom_sk1.o sggeom_sk2.o sggeom_sk3.o cntlcard.o \
	skidentify.o sktrginfo.o skmkruninf.o sgcnt_photon.o sgcnt_aphoton.o \
	skdcy_muminus.o muminus.o skdcy_muplus.o sgpintneut.o getneutxs.o \
        sgnuceff.o rotatevect.o efpncrsh_pnsep.o sggeom_dome.o

DETO =	chersg.o pmtqesg.o rfpmsg.o wtrsg.o  wlngsg.o wsrtsg.o tyvksg.o \
	drefsg.o effnsg.o refsg.o

DSO =	 ds.o dsnois.o dsinip.o dsres.o dsres_multi.o dstrig.o dstrig_multi.o dsprsv.o \
	dsbdel.o dsanti.o dseff_qbee.o dssatu_sk4.o \
        dscard.o dshigh.o dssubchk.o sktsig_le.o sktsig_le_sk3.o dsnois_multi.o dscnt.o \
	dsswtrg.o dsdigi_sk4.o dsdigi_anti_sk4.o dstrgnum.o dshigh_sk4.o dsanti_sk4.o rootSinglePe.o dscptqztoraw.o \
	skqtfitscat.o

SEG =   seg.o seg_openf.o seg_closef.o 363tools.o seg_packprt.o

OBJ =	vecbst.o rfbssg.o dsthr.o dsthr_sk3.o dsthr_sk4.o dsthr_sk5.o dsreso_qb.o \
	rndark.o skfakehd.o aprdcd.o \
	skgamma.o sgmies.o wmiesg.o sgcardbnk.o  skmkmcver.o \
	tunedat.o apdssbread.o apfill2nd.o apvcmpmk.o wsymscsg.o wasyscsg.o \
	skmkapmcstop.o wtmisg.o gauss.o wabssg.o wraysg.o save_muon_path.o \
	nemknebk.o do_validate.o sgpope_read.o sgbsdist.o popefrysg.o \
	skphonuc.o skphonucxs.o disttodet.o \
	wtrgad.o wtgdeff.o
#        fadcwrinit.o fadcrecord.o fadcwrite.o fadcdark.o 

ifdef LASER
GOBJ =  gukine.o gutrak.o gutrev.o gustep.o guphad.o guhadr.o guout.o \
        grun.o gtckov.o ggckov.o gtmuon.o gskpho.o gltrac.o \
        gpart.o gcalor_unix.o skrn1pe.o skrn1pe_sk2.o skrn1pe_sk3.o skrn1pe_sk5.o \
	tr2fac.o trabs.o trabs_sk4.o gtnino.o lsrgen.o rngausyoo.o 
else
ifdef WTSCTT
GOBJ =  gukine.o gutrak.o gutrev.o gustep.o guphad.o guhadr.o guout.o \
        grun.o gtckov.o ggckov.o gtmuon.o gskpho.o gltrac.o \
        gpart.o gcalor_unix.o skrn1pe.o skrn1pe_sk2.o skrn1pe_sk3.o skrn1pe_sk5.o \
	tr2fac.o trabs.o trabs_sk4.o gtnino.o lsrgen.o rngausyoo.o
else
GOBJ =  gukine.o gutrak.o gutrev.o gustep.o guphad.o guhadr.o guout.o \
        grun.o gtckov.o ggckov.o gtmuon.o gtgama.o gthadr.o gskpho.o gltrac.o \
        gpart.o gcalor_unix.o skrn1pe.o skrn1pe_sk2.o skrn1pe_sk3.o skrn1pe_sk5.o \
	tr2fac.o trabs.o trabs_sk4.o gtneut.o klksregen.o 
endif
endif



ifdef SK1LE
DOBJ =
else
DOBJ =  gdrela.o gdrele_h2o.o gdrelm_h2o.o gdrelx_h2o.o gdrelp_h2o.o
endif

ACR =	sgacredge.o sgacrefl.o sgacrefr.o sgacrylic.o sggapedge.o \
	sgrefacr.o acrsg.o acrnsg.o cheracsg.o sgacrnor.o \
        sggeom.o sgpmtnum_sk1.o sgpmtnum_sk2.o sgpmtnum_sk3.o\
	sgpmt.o sgbst.o sktest.o sgholebst.o 

ifdef UPMU
FORTRANDEFINES += -DUPMU
endif

ifdef NEUTRON
FORTRANDEFINES += -DNEUTRON
FORTRANDEFINES += -DDSSWTRIGSTOREINFO
endif

ifdef NICKEL
FORTRANDEFINES += -DNICKEL
endif

ifdef ESCALE 
FORTRANDEFINES += -DESCALE
endif

ifdef LINAC
FORTRANDEFINES += -DLINAC
OBJLIN = array_locate_lin.o getreflect_lin.o lingeom1.o linsum2.o \
	 refstn.o rfpip.o rfpip1_5.o rfpip1_7.o rfwin.o \
	 rffra.o rfpip1.o rfpip1_6.o rfpip2.o
else
OBJLIN =
endif

ifdef DTG
FORTRANDEFINES += -DDTG
OBJDTG = dtgsum.o
else
OBJDTG =
endif

ifdef LASER
FORTRANDEFINES += -DLASER
endif

ifdef WTSCTT
FORTRANDEFINES += -DWTSCTT
endif

ifdef NEUT_ROOT
FORTRANDEFINES += -DSITE_HAS_NEUT

  #LDLIBS += -L$(NEUT_ROOT)/lib/Linux_pc -lneutcore -lnuceff -lneutcore -lnuccorrspl -lpartnuck -lskmcsvc
  #LDLIBS += $(NEUT_ROOT)/src/neutsmpl/rmarin_dum.o

  NEUTCOREVER = 5.4.0
  NUCEFFVER   = 5.4.0
  NUCCORVER   = 5.4.0
  PARTNUCKVER = 5.4.0
  SKMCSVCVER  = 5.4.0
  ZBSFNSVER   = 5.4.0
  SPECFUNCVER = 5.4.0
  RADCORRVER  = 5.4.0
  TAUOLAVER   = 5.4.0
  N1P1HVER    = 5.4.0
  HT2P2HVER   = 5.4.0

  NEUTINCDIR  = ${NEUT_ROOT}/include
  NEUTLIBDIR  = ${NEUT_ROOT}/lib
  NEUTCOREDIR = ${NEUT_ROOT}/src/neutcore
  NEUTCLASSDIR= ${NEUT_ROOT}/src/neutclass


  NEUTDISOBJS = ${NEUTCOREDIR}/structm.o ${NEUTCOREDIR}/pdfset.o \
		${NEUTCOREDIR}/grv94di.o  ${NEUTCOREDIR}/grv98_lo.o

  NEUTCLSOBJS = ${NEUTCLASSDIR}/neutctrl.so    ${NEUTCLASSDIR}/neutpart.so \
		${NEUTCLASSDIR}/neutvect.so    ${NEUTCLASSDIR}/neutvtx.so \
		${NEUTCLASSDIR}/neutfsipart.so ${NEUTCLASSDIR}/neutfsivert.so \
		${NEUTCLASSDIR}/neutnucfsistep.so \
		${NEUTCLASSDIR}/neutnucfsivert.so \
		${NEUTCLASSDIR}/neutfill.o \
		${NEUTCLASSDIR}/neutrootfortwrapper.o \
		${NEUTCLASSDIR}/NeutRootHandlers.o

  NEUTINC  = -I$(NEUT_ROOT)/include
  NEUTLIBS = ${NEUTDISOBJS} ${NEUTCLSOBJS} \
		-L${NEUTLIBDIR} -ltauola_${TAUOLAVER} \
		-lneutcore_${NEUTCOREVER} -lspecfunc_${SPECFUNCVER} \
		-lnuccorspl_${NUCCORVER} \
		-lnuceff_${NUCEFFVER} \
		-lneutcore_${NEUTCOREVER} \
		-lpartnuck_${PARTNUCKVER} -lskmcsvc_${SKMCSVCVER} \
		-lN1p1h_${N1P1HVER}       -lHT2p2h_${HT2P2HVER}   \
		-lradcorr_${RADCORRVER} -lzbsfns_${ZBSFNSVER}     \
		`cernlib jetset74 pdflib804 photos202`

LDLIBS += -L$(NEUTLIBS)

#LDLIBS += $(NEUT_ROOT)/src/neutsmpl/rmarin_dum.o
LDLIBS += $(CLIB)
endif

LOBJS= $(SGO) $(DETO) $(DSO) $(SEG) $(OBJ) $(ACR) $(OBJLIN) $(OBJDTG) $(GOBJ) $(DOBJ) 
OBJS = $(INITOBJ) $(LOBJS)

SUBOBJ  = apdssub.o

EXEC_FILES = skdetsim make_lowmu

.PHONY:  clean setup includes install.includes depend lib install.lib exec install.exec

emptyrule : skdetsim_svn_version skidentify.o lib exec

# These switches select the Zebra store size for the executable
# If you do
#    setenv HIGH 1  
# before compiling (need a clean compile!)
#    then larger zebra store will be selected with ifdefs in the source
#   and an executable named skdetsim_high will be created.
#    setenv MIDDLE 1  
#  will use a medium size store.
# The default is the smallest store size.

ifdef WithG4
SKDETSIM_RUN_PATH = $(A_LIBDIR):$(SKOFL_LIBDIR):$(ROOTLIBDIR):/usr/local/sklib_gcc4.8.2/xerces-c-3.1.1/lib:$(G4TMPDIR):$(G4LIB)/..:$(CLHEP_LIB_DIR)
else
SKDETSIM_RUN_PATH = $(A_LIBDIR):$(SKOFL_LIBDIR):$(ROOTLIBDIR):/usr/local/sklib_gcc4.8.2/xerces-c-3.1.1/lib
endif

make_lowmu : make_lowmu.o skdetsim_svn_version $(OBJS) 
	LD_RUN_PATH=$(SKDETSIM_RUN_PATH) $(CXX) $(CXXFLAGS) $(CPPFLAGS) -o make_lowmu make_lowmu.o $(LDFLAGS) $(G4_GTNEUT) $(LDLIBS) $(LDLIBS2)

skdetsim : skdetsim.o skdetsim_svn_version $(OBJS) 
	LD_RUN_PATH=$(SKDETSIM_RUN_PATH) $(CXX) $(CXXFLAGS) $(CPPFLAGS) -o skdetsim skdetsim.o $(LDFLAGS) $(G4_GTNEUT) $(LDLIBS) $(LDLIBS2)

tmp:
	@echo G4TMPDIR  $(G4TMPDIR)

ifdef HIGH
 FORTRANDEFINES += -DHIGH
 EXEC_FILES += skdetsim_high
else
 ifdef MIDDLE
  FORTRANDEFINES += -DMIDDLE
  EXEC_FILES += skdetsim_middle
 endif
endif 

#
#  Targets 
#


skdetsim_svn_version :
	./skdetsim_svn_version.sh

clean::
	$(RM) *.o *~ *.a *.so core* fort.* dummy.* $(OBJS) $(EXEC_FILES) skdetsim_svn_version.h
	$(RM) -f skdetsim skdetsim_middle skdetsim_high
ifdef WithG4
	$(RM) -rf $(G4_LIBDIR)
 endif 

setup::

includes::
	./skdetsim_svn_version.sh

install.includes::

depend::

lib:: lib$(LIBNAME).a

install.lib:: $(LIBDIR)lib$(LIBNAME).a

exec:: $(EXEC_FILES)

install.exec:: 

#
#  Rules for building library 
#

.PHONY:  lib$(LIBNAME).a $(LIBDIR)lib$(LIBNAME).a

lib$(LIBNAME).a : skdetsim_svn_version $(OBJS) 
	$(RM) $@
	$(AR) $@ $(LOBJS) 
	$(RANLIB) $@

$(LIBDIR)lib$(LIBNAME).a : lib$(LIBNAME).a
	$(RM) $@
	$(INSTALL_LIB) $< $@

skdetsim_high: skdetsim_svn_version $(OBJS)
	rm skdetsim.o
	mv skdetsim skdetsim_high


skdetsim_middle: skdetsim_svn_version $(OBJS)
	rm skdetsim.o
	mv skdetsim skdetsim_middle
