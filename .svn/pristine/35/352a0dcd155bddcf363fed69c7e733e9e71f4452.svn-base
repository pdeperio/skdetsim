***********************************************************************
*     -----------------
      SUBROUTINE SGPINT
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
*       1980.12.12 ; First version   by  T.Mashimo
*       1994.11.28 ; Modified for SK by  Koshio,Y
*
************************************************************************
#include "geant321/gccuts.inc"
#include "geant321/gckine.inc"
#include "geant321/gcking.inc"
#include "geant321/gctrak.inc"
#include "skcnst.h"


      DIMENSION ADATA(4,16),BDATA(4),SLOPE(4)
      DIMENSION DIR(3)

      DATA ADATA/170.,320.,280.,39.
     &          ,270.,370.,340.,49.
     &          ,340.,410.,340.,60.
     &          ,330.,410.,290.,66.
     &          ,225.,385.,285.,72.
     &          ,130.,340.,300.,78.
     &          , 80.,300.,300.,83.
     &          , 50.,290.,290.,88.
     &          , 35.,280.,280.,90.
     &          , 20.,275.,275.,95.
     &          , 15.,270.,270.,95.
     &          , 10.,260.,260.,98.
     &          ,  5.,260.,260.,98.
     &          ,  0.,260.,260.,98.
     &          ,  0.,260.,260.,98.
     &          ,  0.,260.,260.,98./



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
      IF(T.GT.750.) T=750.

      NOINT=IFIX((T-50.)/50.)+1
C     NOINT IS THE NUMBER OF INTERVAL TO WHICH THE MOMENTUM BELONGS.

      DIFF=T-50.*FLOAT(NOINT)
C     DIFF IS THE DIFFERENCE BETWEEN THE MOMENTUM P AND THE MOMENTUM
C     AT THE LEFT END OF THE INTERVAL(NO.NOINT).

      DO 10 I=1,4
         SLOPE(I)=(ADATA(I,NOINT+1)-ADATA(I,NOINT))/50.
         BDATA(I)=ADATA(I,NOINT)+SLOPE(I)*DIFF
 10   CONTINUE
C     BDATA IS THE CROSS SECTION OF THE INTERACTION AT THE MOMENTUM P

      TOTAL=BDATA(1)+BDATA(2)+BDATA(3)+BDATA(4)

C-----DIVIDE THE (0,1) TO THE FOUR PARTS ACCORDING TO THE RATIO OF THE
C     EACH INTERACTION

      DEL1=BDATA(1)/TOTAL
      DEL2=DEL1+BDATA(2)/TOTAL
      DEL3=DEL2+BDATA(3)/TOTAL

C-----GENERATE RANDOM NUMBER
C      R=RNDM()
      CALL GRNDM(R,1)
      IF(R.GE.0..AND.R.LT.DEL1)   GOTO 100
      IF(R.GE.DEL1.AND.R.LT.DEL2) GOTO 200
      IF(R.GE.DEL2.AND.R.LT.DEL3) GOTO 300
      IF(R.GE.DEL3)               GOTO 400
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

      GOTO 500
C     
C---INELASTICSCATTERING (PARTICLE KIND DOES NOT CHANGE)
C     
 300  ISTOP = 1
      CALL SGPIIE(P,COST)

      NGKINE=1

      KCASE = NAMEC(12)
      NMEC = NMEC + 1
      LMEC(NMEC) = 12

      GOTO 500
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
c     2011-05-23: tmw23 assigning the position to the initial pion step position
         GPOS(J,NGKINE)=VECT(J)
 410  CONTINUE
      GKIN(5,NGKINE) = 7.
      TOFD(NGKINE) = 0.

      KCASE = NAMEC(12)
      NMEC = NMEC + 1
      LMEC(NMEC) = 12

      GOTO 999
C     
C   Correct pion for energy lost and scattered angle
C     
 500  CONTINUE
      CONE = ACOS(COST)*180./PI
      CALL RNCONE(VECT(4),CONE,DIR)
      DO 510 I=1,3
         VECT(I+3)=DIR(I)
 510  CONTINUE
      VECT(7)=P/1000.
      GETOT = SQRT(VECT(7)**2+AMASS**2)
      GEKIN = GETOT-AMASS
      CALL GEKBIN

      DO 520 J=1,3
         GKIN(J,NGKINE)=DIR(J)*P/1000.
         GPOS(J,NGKINE)=VECT(J)
 520  CONTINUE

      GKIN(5,NGKINE) = IPART
      TOFD(NGKINE) = 0.

 999  RETURN
      END
