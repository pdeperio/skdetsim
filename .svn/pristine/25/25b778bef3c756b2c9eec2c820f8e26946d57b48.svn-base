***********************************************************************
*     -----------------
      SUBROUTINE SGPINTNEUT
*     -----------------
*
*     (Purpose)
*       GET THE KIND OF INTERACTION(INTKND) FOR THE GIVEN MOMENTUM(P)
*
*     (Input)
*       none
*
*     (Output)
*       none
*
*     (Creation Date and Author)
*       2010.07.09 ; First version   by  P.de Perio
*
************************************************************************
#include "geant321/gccuts.inc"
#include "geant321/gckine.inc"
#include "geant321/gcking.inc"
#include "geant321/gctrak.inc"
#include "skcnst.h"

      REAL tot,absb,inel,scx,dcx,elas
      REAL PABS,PINEL,PSCX,PDCX,PELAS,PROBSUM

      DIMENSION DIR(3)

      IF(GEKIN.LT.CUTHAD)THEN
         ISTOP = 1
         NMEC = NMEC +1
         LMEC(NMEC) = 5
         
         IF(IPART.EQ.8) then
            CALL GDECAY

C        Placeholder proton to track captured pi-
         ELSE

            NGKINE=1
            
            DO I=1,3
               GKIN(I,NGKINE)=0
               GPOS(I,NGKINE)=VECT(I)
            end do
            GKIN(5,NGKINE) = 14
            TOFD(NGKINE) = 0.
            IFLGK(NGKINE) = -1
                        
         endif

         GOTO 999
      ENDIF

      P=VECT(7)*1000.
      T=GEKIN*1000.
      IF(T.GT.380.) then
         write (*,*) 'SGPINTNEUT: Error NEUT cross sections not valid for 
     &tpi (',T,') > 380, setting TPI=380 MeV'
         stop
      end if

C-----Get interaction cross sections
      Call getneutxs(P,tot,absb,inel,scx,dcx,elas)

C      WRITE(*,*) "SGPINTNEUT: ", T, tot, absb, inel, scx, dcx, elas

C-----Calculate THE RATIO OF EACH INTERACTION

      PABS=absb/tot
      PINEL=inel/tot
      PSCX=scx/tot
      PDCX=dcx/tot
      PELAS=elas/tot

      PROBSUM=PABS+PINEL+PSCX+PDCX+PELAS
      IF (PROBSUM.ne.1.) then
         PABS=PABS/PROBSUM
         PINEL=PINEL/PROBSUM
         PSCX=PSCX/PROBSUM
         PDCX=PDCX/PROBSUM
         PELAS=PELAS/PROBSUM
      ENDIF

      PINEL=PABS+PINEL
      PSCX=PINEL+PSCX
      PDCX=PSCX+PDCX
      PELAS=PDCX+PELAS

C      write(*,*) PABS,PINEL,PSCX,PDCX,PELAS

C-----GENERATE RANDOM NUMBER
C      R=RNDM()
      CALL GRNDM(R,1)
      IF(R.GE.0..AND.R.LT.PABS)   GOTO 100
      IF(R.GE.PABS.AND.R.LT.PINEL) GOTO 300
      IF(R.GE.PINEL.AND.R.LT.PSCX) GOTO 400
      IF(R.GE.PSCX.AND.R.LT.PDCX) GOTO 500
      IF(R.GE.PDCX) GOTO 200
      RETURN
C     
C---ABSORPTION
C     
 100  ISTOP = 1
      NMEC = NMEC + 1
      LMEC(NMEC) = 12

C     Generate placeholder proton to store in secondaries stack
      NGKINE=1

      DO J=1,3
         GKIN(J,NGKINE)=0
         GPOS(J,NGKINE)=VECT(J)
      end do
      GKIN(5,NGKINE) = 14
      TOFD(NGKINE) = 0.
      IFLGK(NGKINE) = -1

      KCASE = NAMEC(12)

      GOTO 999
C     
C---ELASTIC SCATTERING
C     
 200  ISTOP = 1

      CALL SGPIEL(P,COST)

      NGKINE=1

      KCASE = NAMEC(13)
      NMEC = NMEC + 1
      LMEC(NMEC) = 13

      GOTO 600
C     
C---INELASTICSCATTERING (PARTICLE KIND DOES NOT CHANGE)
C     
 300  ISTOP = 1
      CALL SGPIIE(P,COST)

      NGKINE=1

      KCASE = NAMEC(12)
      NMEC = NMEC + 1
      LMEC(NMEC) = 12

      GOTO 600
C     
C---INELASTICSCATTERING (CHARGE EXCHANGE)
C     
 400  ISTOP=1

      CALL SGPICH(P,COST)
      NGKINE=1
      CONE = ACOS(COST)*180./PI
      CALL RNCONE(VECT(4),CONE, GKIN(1,NGKINE))
      DO 410 J=1,3
         GKIN(J,NGKINE)=GKIN(J,NGKINE)*P/1000.
         GPOS(J,NGKINE)=VECT(J)
 410  CONTINUE
      GKIN(5,NGKINE) = 7.
      TOFD(NGKINE) = 0.

      KCASE = NAMEC(12)
      NMEC = NMEC + 1
      LMEC(NMEC) = 12

      GOTO 999
C     
C---INELASTICSCATTERING (DOUBLE CHARGE EXCHANGE)
C     
 500  ISTOP=1

      CALL SGPICH(P,COST)
      NGKINE=1
      CONE = ACOS(COST)*180./PI
      CALL RNCONE(VECT(4),CONE, GKIN(1,NGKINE))
      DO 510 J=1,3
         GKIN(J,NGKINE)=GKIN(J,NGKINE)*P/1000.
         GPOS(J,NGKINE)=VECT(J)
 510  CONTINUE
      IF (IPART.EQ.8.) THEN
         GKIN(5,NGKINE) = 9.
      ELSE IF (IPART.EQ.9.) THEN
         GKIN(5,NGKINE) = 8.
      END IF
      TOFD(NGKINE) = 0.

      KCASE = NAMEC(12)
      NMEC = NMEC + 1
      LMEC(NMEC) = 12

      GOTO 999

C     
C   Correct pion for energy lost and scattered angle
C     
 600  CONTINUE
      CONE = ACOS(COST)*180./PI
      CALL RNCONE(VECT(4),CONE,DIR)
      DO 610 I=1,3
         VECT(I+3)=DIR(I)
 610  CONTINUE
      VECT(7)=P/1000.
      GETOT = SQRT(VECT(7)**2+AMASS**2)
      GEKIN = GETOT-AMASS
      CALL GEKBIN

      DO 620 J=1,3
         GKIN(J,NGKINE)=DIR(J)*P/1000.
         GPOS(J,NGKINE)=VECT(J)
 620  CONTINUE

      GKIN(5,NGKINE) = IPART
      TOFD(NGKINE) = 0.


 999  RETURN
      END
