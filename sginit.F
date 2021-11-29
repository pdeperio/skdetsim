************************************************************************
*     -----------------
      SUBROUTINE SGINIT 
*     -----------------
*
*     (Purpose)
*       Initialize ZBS system, reading control parameters from CARDs,
*       define detector geometries, and so on.    
*     (Input)
*       none
*
*     (Output)
*       none
*
*     (Creation Date and Author)
*       1992.04.24 ; First version by J.Kanzaki
*       1997.02.26 ; KZINIT -> APKZINIT by Y.Hayato
*                    Also add COMMON KZBS entry
*       Merged verstion cretated from
*        apdetsim/sggeom.F 1.4 and ledetsim_4.0/sggeom.F
*
************************************************************************

      implicit none

#include "skmc.h"
#include "skparm.h"
#include "gupht.h"
#include "guhigh.h"
#include "guanti.h"
#include "dsdata.h"
#include "dsmulti.h"

#include "geant321/gcphys.inc"
#include "geant321/gctrak.inc"

#ifdef FADC
#include "fadcwrbnk.h"
#endif

      integer    NPDIV,NPSTOREX
      PARAMETER (NPDIV=20)

#ifdef HIGH
      PARAMETER (NPSTOREX=100000000)
#else
#ifdef MIDDLE
      PARAMETER (NPSTOREX=10000000)
#else
      PARAMETER (NPSTOREX=1000000)
#endif
#endif

      integer       IXSTOR,IXDIV(NPDIV),IFENCE,LZBS(NPSTOREX)
      COMMON /KZBS/ IXSTOR,IXDIV       ,IFENCE,LZBS

      real bratio
      integer mode

      DIMENSION BRATIO(6),MODE(6)
      DATA BRATIO/ 0., 0., 0., 0., 0., 0./
      DATA MODE/ 0, 0, 0, 0, 0, 0/

*     common initialisation

      CALL VZERO(NOPH,   MAXRNG)
      CALL VZERO(TGESK,  MAXPM*MAXRNG)
      CALL VZERO(NQGESK, MAXPM*MAXRNG)
      CALL VZERO(TGESKA, MAXPMA*MAXRNG)
      CALL VZERO(NQGESKA,MAXPMA*MAXRNG)

      call vzero(pmt_id,MXGATE)
      call vzero(nphoto_gate,MXGATE)
      call vzero(timing_gate,MXGATE)
      call vzero(width_photons,MXGATE)
      call vzero(list_gate,MAXPM*50)
      call vzero(num_gate,MAXPM)
      num_gate_all = 0

      call clrmuver

*     Initialization of Zebra Bank System

      CALL KZOPTN('NO-ZEBRA-INIT')
C     CALL KZINIT
C     CALL APKZINIT(NPSTOREX)
      CALL KZLIMIT(NPSTOREX)

*     Initialisation of GEANT3 variables

      CALL GINIT

*     Initialization of DPHYS1 ( = IABAN ) for 
*     electron tracking
*
      DPHYS1 = 0

*     Print job and environment information
      CALL SKIDENTIFY(6)
      
*     Input data cards and initialize detector simulation parameter

      CALL CNTLCARD
      CALL SGCARD
      CALL SGFFIN
      CALL DSINIP
      CALL DSCARD
      CALL APRDCD

*     Fill missing channel information into COMMON
      Call skmisch

C-T.O Initialize for Flash ADC
#ifdef FADC
      if (IFADCI) then
*   FADC initialization
         Call fadcwrinit
         Call fadccontabrd
      endif
#endif

*     Initialization of ZEBRA core divisions and link areas

      CALL GZINIT

*     Creation of the 'particle' data structure JPART

      CALL GPART
      CALL GSDK(9,BRATIO,MODE)

*     Creation of the 'material' data structure JMATE

      CALL GMATE

*     User code for the description of the geometrical setup and
*     the sensitive detectors and for the creation of data structures
*     JVOLUM, JTMED, JROTM and JSETS

      CALL SGGEOM

*     Preparation of cross-section and energy loss tables for 
*     all used materials

      CALL GPHYSI

      RETURN
      END
