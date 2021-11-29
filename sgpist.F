***********************************************************************
*     ------------------------------
      SUBROUTINE SGPIST(P,IP,PIWINT)
*     ------------------------------
*
*     (Purpose)
*       Calculate interaction probability for PI in water
*
*     (Input)
*       P      ; Momemtum of the particle (MeV)
*       IP     ; Particle identification number
*                8 --- Pion +
*                9 --- Pion -
*
*     (Output)
*       PIWINT ; Interaction probability ( / cm )
*
*     (Creation Date and Author)
*       1980.09.30 ; First version   by  K.Arisaka
*       1994.11.28 ; Modified for SK by  Koshio,Y
*       2010.07.09 ; Can use total cross section tuned to data and NEUT 
*                    calculation below 100 MeV/c (by P.de Perio)
*
************************************************************************
C
C     ROH ; H-NUCLEON DENCITY IN WATER (/ CM**2)
C     ROO ; O-NUCLEUS DENCITY IN WATER (/ CM**2)

#include "apcard.h"

      DATA ROH,ROO/6.689E22,3.344E22/

      REAL SIGEL,SIGINEL,SIGCX,SIGTOT
      INTEGER PICHARGE,NUCCHARGE,IERROR

      SIGEL=0
      SIGINEL=0
      SIGCX=0
      
      PIWINT=0.

      if (IAPIMOD.eq.1) then
         PIWINT=(PIPINT(P,IP)*ROH + PIOINTNEUT(P,IP)*ROO)*1.E-24
C         PIWINT=PIOINTNEUT(P,IP)

      else if (IAPIMOD.eq.2) then

         NUCCHARGE=1
         if (IP.eq.8) then
            PICHARGE=1
         else if (IP.eq.9) then
            PICHARGE=-1
         end if
         
         call EFPNCRSH_PNSEP(P,PICHARGE,NUCCHARGE,SIGEL,SIGINEL,SIGCX,IERROR)

         SIGTOT=(SIGEL+SIGINEL+SIGCX)/1000

         PIWINT=(SIGTOT*ROH + PIOINTNEUT(P,IP)*ROO)*1.E-24

      else
         PIWINT=(PIPINT(P,IP)*ROH + PIOINT(P,IP)*ROO)*1.E-24
C         PIWINT=PIOINT(P,IP)

      end if

      RETURN
      END
C
C
      FUNCTION PIPINT(P,IP)
C
C     J. KANZAKI     80-09-30
C
C---CALCULATE TOTAL CROSS SECTION FOR PI-P INTERACTION
C
C     IP= 8 ; PI+
C     IP= 9 ; PI-
C
C      UNIT ; BARN
C
      DIMENSION SIG(39,4,2)
      DATA SIG/  6.2, 26.1, 68.5,168.0,191.1,123.4, 73.8, 46.7, 35.6,
     &          27.8, 21.7, 20.3, 17.4, 15.7, 15.5, 18.0, 20.7, 23.6,
     &          25.1, 26.1, 26.9, 28.0, 29.5, 31.7, 34.4, 37.8, 40.0,
     &          40.9, 40.2, 38.6, 36.5, 34.8, 33.1, 31.8, 30.9, 30.1,
     &          29.6, 29.2, 29.0,
     &          29.0, 29.4, 30.1, 30.8, 30.1, 29.2, 28.6, 28.3, 28.2,
     &          28.0, 27.8, 28*0.0,
     &          27.8, 25.7, 25.1, 24.8, 35*0.0,
     &          24.8, 23.5, 23.3, 23.1, 23.1, 23.3, 33*0.0,
     &          27.9, 55.4, 65.9, 45.0, 32.0, 26.4, 27.1, 28.5, 30.4,
     &          36.7, 43.0, 46.5, 36.7, 38.1, 41.4, 52.7, 61.2, 57.9,
     &          47.2, 39.9, 37.1, 18*0.0,
     &          37.1, 36.0, 34.4, 34.7, 35.8, 35.7, 34.7, 33.3, 32.8,
     &          31.9, 31.6, 31.2, 31.4, 30.9, 30.6, 24*0.0,
     &          30.6, 28.1, 27.5, 26.5, 26.3, 26.0, 25.1, 25.4, 25.2,
     &          30*0.0,
     &          25.2, 24.8, 24.4, 24.3, 24.2, 34*0.0/

      DIMENSION PDIV(4)
      DATA PDIV/ 50., 200., 2000., 10000./

      DIMENSION PLIM(5,2)
      DATA PLIM/ 100., 2000., 4000., 10000., 60000.,
     &           200., 1200., 4000., 20000., 60000./

      DIMENSION SIGASY(2)
      DATA SIGASY/ 23.3, 24.2/

      PIPINT=0.0

      IC = IP - 7
     
      IF(P.LT.PLIM(1,IC)) RETURN
     
      DO 10 I=2,5
         IF(P.GT.PLIM(I,IC)) GO TO 10
         N=I-1
         GO TO 20
   10 CONTINUE
      PIPINT=SIGASY(IC)/1000.
      RETURN
     
   20 PDIV0=PDIV(N)
      PLOW=PLIM(N,IC)
C     IP -> IPX not to modify value in IP   K.Okumura  
C      IP=INT((P-PLOW)/PDIV0)+1
C      PIPINT=SIG(IP,N,IC)+(SIG(IP+1,N,IC)-SIG(IP,N,IC))
C     &     *(P-PLOW-PDIV0*FLOAT(IP-1))/PDIV0
      IPX=INT((P-PLOW)/PDIV0)+1
      PIPINT=SIG(IPX,N,IC)+(SIG(IPX+1,N,IC)-SIG(IPX,N,IC))
     &     *(P-PLOW-PDIV0*FLOAT(IPX-1))/PDIV0
      PIPINT=PIPINT/1000.0
     
      RETURN
      END
C
C
      FUNCTION PIOINT(P,IP)
C
C     J. KANZAKI     80-09-30
C
C---CALCULATE TOTAL CROSS SECTION FORT PI-O INTERACTION
C
C     IP= 8 ; PI+
C     IP= 9 ; PI-
C
C     UNIT   ; BARN
C     REF.  A.S.CARROLL ET AL., PHYS. REV. C14(1976)635
C
      DATA E0,E1,C,GAM2,RR/ 1186.7, 2.36, 0.31, 129.6, 3.55/
      DATA HP/ 197.3/, PI/ 3.1415927/
      DATA OXM,PIM,PRM/ 14900., 139.57, 939./

      CMM=OXM**2*P**2/(PIM**2+OXM**2+2.*OXM*SQRT(P**2+PIM**2))
      E=SQRT(PIM**2+CMM)+SQRT(OXM**2+CMM)-OXM+PRM
      CMM=SQRT(CMM)
      GAM1=C*CMM
      GAM=GAM1+GAM2
      ALM=RR*CMM/HP+1
      ALAM=HP/CMM/10.

      SIG=ATAN(2.*(E0+E1*ALM**2-E)/GAM)-ATAN(2.*(E0-E)/GAM)
      PIOINT=SIG*2.*PI*ALAM**2*GAM1/E1

      RETURN
      END


      FUNCTION PIOINTNEUT(P,IP)
C
C     P. de Perio     2010-07-09
C
C---CALCULATE TOTAL CROSS SECTION FORT PI-O INTERACTION
C
C     IP= 8 ; PI+
C     IP= 9 ; PI- 
C
C     UNIT   ; BARN
C     REF.  A.S.CARROLL ET AL., PHYS. REV. C14(1976)635
C           Parameterization re-fitted to pi+ data including
C           Saunders et al., Phys. Rev. C53(1996)1745.
C        
C           Use NEUT calculation below 100 MeV/c where no data
C           exists
C
C      These are the newly fitted values using abovementioned data:
C     
C      DATA E0,E1,C,GAM2,RR/ 1178.79, 8.8193, 0.31, 129.59, 21.4709/
C      DATA HP/ 197.3/, PI/ 3.1415927/
C      DATA OXM,PIM,PRM/ 14900., 139.57, 939./
C
C      The above function is evaluated externally and input via
C      table below.

      real dum,tot

      call getneutxs(P,tot,dum,dum,dum,dum,dum)
      PIOINTNEUT=tot/1000.

      RETURN
      END
