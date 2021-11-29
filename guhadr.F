CDECK  ID>, GUHADR.
*CMZ :  3.15/01 26/06/91  11.45.52  by  Federico Carminati
*-- Author :
      SUBROUTINE GUHADR
C.
C.    ******************************************************************
C.    *                                                                *
C.    *       User routine to generate one hadronic interaction        *
C.    *                                                                *
C.    *    ==>Called by : GTHADR,GTNEUT                                *
C.    *                                                                *
C.    ******************************************************************
C.
      IMPLICIT NONE

#include "geant321/gckine.inc"
#include "geant321/gcphys.inc"
C-ADD Y.H./96/11/16 (Add GCTRAK)<-TO get momentum
#include "geant321/gctrak.inc"
#include "apcard.h"

C local definition

      real rdummy

C.
C.    ------------------------------------------------------------------
C.
C          GHEISHA only if IHADR<3 (default)
C          GHEISHA and HADRIN/NUCRIN if IHADR=3
C
c     CALL GHEISH

C      IF(IPART.EQ.8.OR.IPART.EQ.9)THEN

C-IF PION and MOMENTUM < RAPITHR(MeV) THEN USE ORIGINAL
C in case IAPION==1.
C
C If IAPIMOD==1, use cross sections calculated by NEUT
C If IAPIMOD==2, use cross sections calculated by NEUT
C                and NEUT-NUCEFF for final state

      IF(     (IPART.EQ.8.OR.IPART.EQ.9)
     $   .and.(VECT(7).lt.RAPITHR/1000.)
     $   .and.(IAPION.eq.1))THEN

         if (IAPIMOD.eq.1) then
            CALL SGPINTNEUT
#ifdef SITE_HAS_NEUT
         else if (IAPIMOD.eq.2) then
            CALL SGNUCEFF
#endif
         else
            CALL SGPINT
         end if

C         ZINTHA=-LOG(RNDM())
         CALL GRNDM(RDUMMY,1)
         ZINTHA=-LOG(RDUMMY)

      ELSE
         IF (IAPHADR.eq.IFLUKA) THEN
            CALL FLUFIN
            GO TO 100
         ENDIF
         IF (IAPHADR.eq.IGHEIS) THEN
            CALL GHEISH
            GO TO 100
         ENDIF
         IF (IAPHADR.eq.ICALOR) THEN
            CALL GCALOR
         ENDIF


      ENDIF
 100  CONTINUE
C
      END
