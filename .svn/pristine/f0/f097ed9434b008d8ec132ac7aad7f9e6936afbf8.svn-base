CDECK  ID>, GUPHAD.
*CMZ :  3.14/14 23/05/90  18.26.11  by  Rene Brun
*-- Author :
      SUBROUTINE GUPHAD
C.
C.    ******************************************************************
C.    *                                                                *
C.    *       User routine to compute Hadron. inter. probabilities     *
C.    *                                                                *
C.    *    ==>Called by : GTHADR,GTNEUT                                *
C.    *                                                                *
C.    ******************************************************************
C.
      IMPLICIT NONE

#include "geant321/gckine.inc"
#include "geant321/gcphys.inc"
#include "geant321/gctrak.inc"
#include "apcard.h"

C local definitions

      real p, piwint
      integer ip
C.
C.    ------------------------------------------------------------------
C.
C             GPHADR for TATINA and CASCADE84
C             GPGHEI for GHEISHA
C
c     CALL GPGHEI

C---      IF(IPART.EQ.8.OR.IPART.EQ.9)THEN
C-IF PION and MOMENTUM < RAPITHR(MeV) THEN USE ORIGINAL
C in case IAPION==1.

      IF(     (IPART.EQ.8.OR.IPART.EQ.9)
     $     .and.(VECT(7).lt.RAPITHR/1000.)
     $     .and.(IAPION.eq.1))THEN
         IP = IPART
         P  = VECT(7)*1000.
         CALL SGPIST(P,IP,PIWINT)
         STEPHA = 1./PIWINT
         SHADR  = STEPHA*ZINTHA
      ELSE
         IF (IAPHADR.eq.IFLUKA) THEN
            CALL FLDIST
            GO TO 100
         ENDIF
         IF (IAPHADR.eq.IGHEIS) THEN
            CALL GPGHEI
            GO TO 100
         ENDIF
         IF (IAPHADR.eq.ICALOR) THEN
            CALL CALSIG
         ENDIF
      ENDIF
 100  CONTINUE
      END
