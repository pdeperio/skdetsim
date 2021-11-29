      SUBROUTINE KLKSREGEN
C.
C.    ******************************************************************
C.    *                                                                *
C.    *   Generates a daughter KS from the parent KL with same         *
C.    *   kinematics as the parent (only particle id is changed)       *
C.    *                                                                *
C.    *   ==>Called by : GTNEUT                                        *
C.    *      Authors     C.REGIS                                       *
C.    *                                                                *
C.    ******************************************************************
C.
#include "geant321/gcbank.inc"
#include "geant321/gccuts.inc"
#include "geant321/gckine.inc"
#include "geant321/gcking.inc"
#include "geant321/gconsp.inc"
#include "geant321/gctrak.inc"
#include "geant321/gcphys.inc"
#include "geant321/gcstak.inc"
#include "geant321/gctmed.inc"
#include "geant321/gctrak.inc"
#if defined(CERNLIB_USRJMP)
#include "geant321/gcjump.inc"
#endif
#if defined(CERNLIB_DEBUG)
#include "geant321/gcunit.inc"
#endif

      include "klksregen.inc"
      
      NGKINE=1
      DO K=1,NGKINE
         GKIN(1,K)=VECT(4)*VECT(7)
         GKIN(2,K)=VECT(5)*VECT(7)
         GKIN(3,K)=VECT(6)*VECT(7)
         GKIN(4,K)=GETOT
         GKIN(5,K)=16
         TOFD(K)=0.
         IFLGK(K)=0
         GPOS(1,K) = VECT(1)
         GPOS(2,K) = VECT(2)
         GPOS(3,K) = VECT(3)
         print*,'cregis klksregen.F make daughter',(gkin(J,K),J=1,5),':'
*     +        ,(gpos(J,K),J=1,3)
      ENDDO

      CALL UCTOH('KLKS',KCASE,4,4)

      END
