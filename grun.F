CDECK  ID>, GRUN.   
*CMZ :  3.15/01 03/02/92  19.29.07  by  Federico Carminati
*-- Author :
      SUBROUTINE GRUN
*.
*.    ******************************************************************
*.    *                                                                *
*.    *       Steering routine to process all the events               *
*.    *                                                                *
*.    *    ==>Called by : <USER>, main program                         *
*.    *       Author    R.Brun  *********                              *
*.    *                                                                *
*.    ******************************************************************
*.
      IMPLICIT NONE
#include "gcbank.h"
#include "geant321/gcflag.inc"
#include "geant321/gcunit.inc"
#include "geant321/gctime.inc"
#include "skdetsim_parm.h"
C
C local defs for implicit none
C Hey!  Is timnow being used before it's initialized?  If you know about 
C this, please check it
      real timnow
      integer ifirst
C
      SAVE IFIRST
      DATA IFIRST/0/
*.
*.    ------------------------------------------------------------------
*.
*             Keep starting time
*
      IF(IFIRST.EQ.0)THEN
         IFIRST=1
         CALL TIMEL(TIMINT)
      ENDIF

#ifdef WithG4
*     initialize geant4
c      print *,"grun (MHARADA), GDWTR, USEG4", GD_WTR, USEG4 
      IF(USEG4==1) CALL g4_initialize(G4_INITI_RAND, GD_WTR, GD_MODEL, GD_CONC)
*
#endif

   10 IF(IEVENT.LT.NEVENT) THEN
         IEVENT=IEVENT+1
*
*             Initialises event partition
*
         CALL GTRIGI
*
*             Process one event (trigger)
*
         CALL GTRIG
*
*             Clear event partition
*
         CALL GTRIGC
*
         IF(IEORUN.EQ.0) THEN
*
*             Check time left
*
            IF(ITIME.LE.0)GO TO 10
            IF(MOD(IEVENT,ITIME).NE.0)GO TO 10
            CALL TIMEL(TIMNOW)
C no time limit
            TIMNOW = 1.E10 
            IF(TIMNOW.GT.TIMEND)GO TO 10
            WRITE(CHMAIL,10000)TIMEND
            CALL GMAIL(0,2)
            IEORUN = 1
         ENDIF
      ENDIF
*

#ifdef WithG4
*     finalize geant4
      IF (USEG4==1) CALL g4_finalize
#endif

10000 FORMAT(5X,'***** THE JOB STOPS NOW BECAUSE THE TIME LEFT IS LESS',
     +' THAN     ',F8.3,' SECONDS *****')
      END
