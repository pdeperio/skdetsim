*
* $Id: gtneut.F,v 1.1.1.1 1995/10/24 10:21:44 cernlib Exp $
*
* $Log: gtneut.F,v $
* Revision 1.1.1.1  1995/10/24 10:21:44  cernlib
* Geant
*
*
#include "geant321/pilot.h"
*CMZ :  3.21/02 29/03/94  15.41.24  by  S.Giani
*-- Author :
      SUBROUTINE GTNEUT
C.
C.    ******************************************************************
C.    *                                                                *
C.    *   Neutral hadron type track. Computes step size and propagates *
C.    *    particle through step.                                      *
C.    *                                                                *
C.    *   ==>Called by : GTRACK                                        *
C.    *      Authors     R.Brun, F.Bruyant  *********                  *
C.    *                                                                *
C.    ******************************************************************
C.
#include "geant321/gcbank.inc"
#include "geant321/gccuts.inc"
#include "geant321/gckine.inc"
#include "geant321/gcking.inc"
#include "geant321/gconsp.inc"
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
#if !defined(CERNLIB_SINGLE)
      PARAMETER (EPSMAC=1.E-6)
#endif
#if defined(CERNLIB_SINGLE)
      PARAMETER (EPSMAC=1.E-11)
#endif

#ifdef WithG4
#include "gdneutron.h"
#include "skdetsim_parm.h"
#endif

c cregis include klksregeneration
      include "klksregen.inc"
C.
C.    ------------------------------------------------------------------
*
* *** Particle below energy threshold ? Special treatment
*

#ifdef WithG4
C     ------------------------------------------------------------------
C       parameters added by M.Harada(original is by Y.Ashida, R.Akutsu,
C       Ka Ming 
C     ------------------------------------------------------------------

      LOGICAL insideVolume, needGd, needG4Water
      INTEGER ITER, testNumSec, ITER2, 
     +        gdCaptureCount, gdElaCount, gdInElaCount,  HCapCount
      REAL Displacement, radius, dz
      INTEGER numSecondary, neutronCount, gammaCount, processNumber
      REAL*8 aRndm, secVec(9,100), captureTime, sumGamma, PP
     +       totXsec, wtrXsec, 
     +       HcapXsec, HelaXsec, HinelaXsec, 
     +       OcapXsec, OelaXsec, OinelaXsec, 
     +       GdcapXsec, GdelaXsec, GdinelaXsec
      DATA neutronCount/0/, gdCaptureCount/0/, radius/20000./
     +     dz/20000./, gdElaCount/0/, gdInElaCount/0/, HCapCount/0/

C     ------------------------------------------------------------------


C     ------------------------------------------------------------------
C       Activate Gd-simulation if satisfy several conditions. 
C        1. inside the inner water tank volume <-- no need now
C        2. neutron energy < 20 MeV 
C        3. GD_WTR == 1 as defined in the input card file 
C        4. current material is water 
C        added by M.Harada(original is by Y.Ashida, R.Akutsu, Ka Ming
C     ------------------------------------------------------------------

      insideVolume = .TRUE.
      needG4Water = .FALSE.

C      print *,"check (MHARADA), GDWTR, USEG4, needGd", GD_WTR, USEG4, needGd 
      IF ((VECT(1)*VECT(1)+VECT(2)*VECT(2)>radius*radius) .OR. 
     $    abs(VECT(3))>dz) insideVolume=.FALSE.
C     if insideVolume is used, Gd-capture only happens in ID. Want OD events as well?
      needGd = .FALSE.
      IF (GEKIN<0.02 .AND. USEG4==1 .AND. NUMED==2) needG4Water = .TRUE.
      IF (GEKIN<0.02 .AND. USEG4==1 .AND. GD_WTR==1
     &   .AND. NUMED==2) 
     $   needGd = .TRUE.

C     ------------------------------------------------------------------
#endif

C     ------------------------------------------------------------------
C       If particle below energy threshold, go to special treatment. 
C       Added by M.Harada(original is by Y.Ashida, R.Akutsu, Ka Ming
C     ------------------------------------------------------------------

      IF (GEKIN.LE.CUTNEU) THEN
         GEKIN  = 0.
         GETOT  = AMASS
         VECT(7)= 0.
         ISTOP  = 2
         NMEC = NMEC + 1
         LMEC(NMEC) = 30
         IF (IHADR.EQ.0.AND.AMASS.GT.0.) THEN
            IF (IDCAY.NE.0) THEN
               GAMMA  = GETOT/AMASS
               TOFG = TOFG +GAMMA*SUMLIF/CLIGHT
               SUMLIF = 0.
               IF (TOFG.GE.TOFMAX) THEN
                  ISTOP = 4
                  NMEC = NMEC + 1
                  LMEC(NMEC) = 22
                  NGKINE = 0
               ELSE
                  NMEC = NMEC + 1
                  LMEC(NMEC) = 5
                  ISTOP =1
                  CALL GDECAY
               ENDIF
            ENDIF
            GO TO 999
         ENDIF
         IPROC = 12
         GO TO 40
      ENDIF

*
* *** Compute step size
*
      IPROC  = 103
      STEP   = STEMAX
*
*  **   Step limitation due to hadron interaction ?
*
      IF (IHADR.GT.0) THEN

#ifdef WithG4
C       Get random number following uniform distribution for
C       determining interaction type ('aRndm').
        CALL g4_getrndm(aRndm)
C        write(*,*) 'aRndm = ',aRndm

C       Initialize parameters for each cross section. 
C       (n-H/O/Gd, and g3Sig for the others)
        HcapXsec = 0
        HelaXsec = 0
        HinelaXsec = 0
        OcapXsec = 0
        OelaXsec = 0
        OinelaXsec = 0
        GdcapXsec = 0 
        GdelaXsec = 0
        GdinelaXsec = 0
        g3Sig = 0

C
C  ***  If the particle is neutron (IPART==13), load GEANT4 cross sections.
C       If you want default GEANT3 simulation, comment out all the parts
C       regarding GEANT4 and read GUPHAD instead. 
C
C       Why having this condition? gtneut only handles neutron anyway
        IF (IPART.EQ.13) THEN 


C         GEANT3 cross section if necessary  
          wtrXsec = 0.
          if (needG4Water) then
C         Load neutron-water cross section, and calculate interaction length. 
            CALL g4_gethxsec(GEKIN,HcapXsec,HelaXsec,HinelaXsec) 
            CALL g4_getoxsec(GEKIN,OcapXsec,OelaXsec,OinelaXsec) 
            wtrXsec = HcapXsec+HelaXsec+HinelaXsec+OcapXsec+OelaXsec
     $          +OinelaXsec
            if (wtrXsec.gt.0) SHADR = ZINTHA/wtrXsec

c          write(*,*) 'E_n [GeV] = ',GEKIN
c          write(*,*) 'sigma_G3water [/cm] = ',g3Sig
c          write(*,*) 'sigma_G4water [/cm] = ',wtrXsec

C         One possible problem here is the "impurity" in salt, e.g. sulfur atoms. Okay to completely neglect?
C         Get cross sections for neutron-Gd interactions and modify step size.
            if (needGd) CALL g4_getgdxsec(GEKIN,GdcapXsec,GdelaXsec,
     $                                     GdinelaXsec)
            if (wtrXsec.gt.0 .AND. needGd) SHADR = SHADR*wtrXsec
     $           /(wtrXsec+GdcapXsec+GdelaXsec+GdinelaXsec)
          else 
#if !defined(CERNLIB_USRJMP)
            CALL GUPHAD
#endif
#if defined(CERNLIB_USRJMP)
            CALL JUMPT0(JUPHAD)
#endif
          endif
C
C  ***  If the particle is not neutron, or no GEANT4, load GEANT3 cross sections.
C
        ELSE
#endif  
C#endif  // WithG4

#if !defined(CERNLIB_USRJMP)
         CALL GUPHAD

#endif
#if defined(CERNLIB_USRJMP)
         CALL JUMPT0(JUPHAD)
#endif

#ifdef WithG4
        ENDIF
#endif

        IF (SHADR.LT.STEP) THEN
            IF (SHADR.LE.0.) SHADR = PREC
            STEP  = SHADR
            IPROC = 12
        ENDIF
      ENDIF
*
*  **   Step limitation due to decay ?
*
      IF ((IDCAY.NE.0).AND.(AMASS.GT.0.)) THEN
         SDCAY = SUMLIF*VECT(7)/AMASS
         IF (SDCAY.LT.STEP) THEN
            STEP  = SDCAY
            IPROC = 5
         ENDIF
      ENDIF
*

*
*  ** Step limitation due to KL->KS regeneration (C.REGIS)
*
      IF (IPART.EQ.10 .AND. KLKSENABLE) THEN
         KLKSREGL = KLKSLIF*VECT(7)/AMASS
         IF (KLKSREGL.LT.STEP) THEN
            STEP  = KLKSREGL
            IPROC = 42
         ENDIF
      ENDIF      


c      if (ipart.eq.10 .or. ipart.eq.16) print*,
c     +     'cregis gtneut', shadr, sdcay
c     +     ,klksregl,ipart

      IF (STEP.LT.0.) STEP = 0.
*
*  **   Step limitation due to geometry ?
*
      IF (STEP.GE.SAFETY) THEN
         CALL GTNEXT
         IF (IGNEXT.NE.0) THEN
            STEP = SNEXT + PREC
            IPROC = 0
            INWVOL = 2
            NMEC = NMEC + 1
            LMEC(NMEC)  = 1
         ENDIF
*
*        Update SAFETY in stack companions, if any
         IF (IQ(JSTAK+3).NE.0) THEN
            DO 10 IST = IQ(JSTAK+3),IQ(JSTAK+1)
               JST = JSTAK +3 +(IST-1)*NWSTAK
               Q(JST+11) = SAFETY
   10       CONTINUE
            IQ(JSTAK+3) = 0
         ENDIF
*
      ELSE
         IQ(JSTAK+3) = 0
      ENDIF
*
* *** Linear transport
*
      IF (INWVOL.EQ.2) THEN
         DO 20 I = 1,3
            VECTMP  = VECT(I) +STEP*VECT(I+3)
            IF(VECTMP.EQ.VECT(I)) THEN
*
* *** Correct for machine precision
*
               IF(VECT(I+3).NE.0.) THEN
                  VECTMP = VECT(I)+ABS(VECT(I))*SIGN(1.,VECT(I+3))*
     +            EPSMAC
                  IF(NMEC.GT.0) THEN
                     IF(LMEC(NMEC).EQ.104) NMEC=NMEC-1
                  ENDIF
                  NMEC=NMEC+1
                  LMEC(NMEC)=104
#if defined(CERNLIB_DEBUG)
                  WRITE(CHMAIL, 10000)
                  CALL GMAIL(0,0)
                  WRITE(CHMAIL, 10100) GEKIN, NUMED, STEP, SNEXT
                  CALL GMAIL(0,0)
10000 FORMAT(' Boundary correction in GTNEUT: ',
     +       '    GEKIN      NUMED       STEP      SNEXT')
10100 FORMAT(31X,E10.3,1X,I10,1X,E10.3,1X,E10.3,1X)
#endif
               ENDIF
            ENDIF
            VECT(I) = VECTMP
   20    CONTINUE
      ELSE
         DO 30 I = 1,3
            VECT(I)  = VECT(I) +STEP*VECT(I+3)
   30    CONTINUE
      ENDIF
*
      SLENG = SLENG +STEP
*
* *** Update time of flight
*
      KLKSLIF = KLKSLIF -STEP*AMASS/VECT(7)
      SUMLIF = SUMLIF -STEP*AMASS/VECT(7)
      TOFG = TOFG +STEP*GETOT/(VECT(7)*CLIGHT)
      IF (TOFG.GE.TOFMAX) THEN
         ISTOP = 4
         NMEC  = NMEC +1
         LMEC(NMEC) = 22
         GO TO 999
      ENDIF
*
* *** Update interaction probabilities
*
      IF (IHADR.GT.0) ZINTHA = ZINTHA*(1.-STEP/SHADR)
*
* *** apply the selected process if any
*
   40 IF (IPROC.EQ.0) GO TO 999
      NMEC = NMEC +1
      LMEC(NMEC) = IPROC
*
*  **   Hadron interaction ?
*

#ifdef WithG4
C     ------------------------------------------------------------------
C       Hadron interaction (IPROC==12) 
C        - Determine first Water or GdWater based on needGd. 
C        - If processNumber > 0, CALL g4_run().  
C        - processNumber 
C           = 0  : not-neutron interaction with water (by GEANT3)
C           = 10 : neutron-H/O    capture             (by GEANT4)
C           = 11 : neutron-H/O    elastic             (by GEANT4)
C           = 12 : neutron-H/O    inelastic           (by GEANT4)
C           = 20 : neutron-H/O/Gd capture             (by GEANT4)
C           = 21 : neutron-H/O/Gd elastic             (by GEANT4)
C           = 22 : neutron-H/O/Gd inelastic           (by GEANT4)
C       Added by M.Harada(original is by Y.Ashida, R.Akutsu, Ka Ming
C     ------------------------------------------------------------------
#endif

      IF (IPROC.EQ.12) THEN

#ifdef WithG4
        processNumber = 0

C ***   1. Interact with pure water
C
        if ((wtrXsec.gt.0).AND.(needGd.eqv..false.)) then 

C         total-Xsec = not-n + n-H/O-{cap,ela,inela}
C         g3Sig = 0 for all cases?
          totXsec = g3Sig+wtrXsec 
C         not-neutron interaction with water
C         What is the not-neutron interaction? Legacy from Ashida's code?
          if (g3Sig/totXsec.gt.aRndm) then 
            processNumber = 0
C         neutron-H/O capture
          else if ((g3Sig+HcapXsec+OcapXsec)/totXsec.gt.aRndm) then 
            processNumber = 10
C         neutron-H/O elastic
          else if ((g3Sig+HcapXsec+OcapXsec+HelaXsec+OelaXsec)/
     $       totXsec.gt.aRndm) then 
            processNumber = 11
C         neutron-H/O inelastic
          else
            processNumber = 12
          endif

C ***   2. Interact with Gd water
C
        else if ((wtrXsec.gt.0) .AND. (needGd.eqv..true.)) then 

C         total-Xsec = not-n + n-H/O/Gd-{cap,ela,inela}
C          print *, wtrXsec, GdcapXsec, GdelaXsec, GdinelaXsec
          totXsec = g3Sig+wtrXsec+GdcapXsec+GdelaXsec+GdinelaXsec 
C         not-neutron interaction with water
          if (g3Sig/totXsec.gt.aRndm) then 
            processNumber = 0
C         neutron-H/O/Gd capture
          else if ((g3Sig+HcapXsec+OcapXsec+GdcapXsec)/totXsec.gt.aRndm) 
     $       then 
            processNumber = 20
C            print *, GdcapXsec, HcapXsec, OcapXsec
C         neutron-H/O/Gd elastic
          else if ((g3Sig+HcapXsec+OcapXsec+GdcapXsec+HelaXsec+OelaXsec+
     $       GdelaXsec)/totXsec.gt.aRndm) then 
            processNumber = 21
C         neutron-H/O/Gd inelastic
          else
            processNumber = 22
          endif

        endif


C ***   Get secondary products. 
C
C        write(*,*) 'random number = ',aRndm
C        write(*,*) 'processNumber (G3) = ',processNumber
        if (needG4Water) then
          if (processNumber.gt.0) then
            DO ITER=1,100
              DO ITER2=1,9
                secVec(ITER2,ITER) = 0
              ENDDO      
            ENDDO

            CALL g4_run(VECT(1),VECT(2),VECT(3),
     $                  VECT(4)*VECT(7),VECT(5)*VECT(7),VECT(6)*VECT(7),
     $                  secVec,processNumber)
            GO TO 590
          endif
        endif
#endif

#if !defined(CERNLIB_USRJMP)
c         if (ipart.ne.13) print*, 'cregis gtneut.f called guhadr: ipart'
c     +        , IPART
         CALL GUHADR
#endif
#if defined(CERNLIB_USRJMP)
         CALL JUMPT0(JUHADR)
#endif
*          Check time cut-off for decays at rest
         IF (LMEC(NMEC).EQ.5) THEN
            TOFG   = TOFG +SUMLIF/CLIGHT
            SUMLIF = 0.
            IF (TOFG.GE.TOFMAX) THEN
               ISTOP = 4
               LMEC(NMEC) = 22
               NGKINE = 0
            ENDIF
         ENDIF
*
*  **   Decay ?
*
      ELSE IF (IPROC.EQ.5) THEN
         ISTOP = 1
c         print*,'cregis gtneut.f called gdecay: ipart', IPART
         CALL GDECAY

*     Apr 2010 C.Regis 
*     Kl->Ks regeneration 
*     stop the tracking of KL and generate a daughter KS with same
c     kinematics as parent KL
      ELSE IF (IPROC.EQ.42) THEN
         ISTOP = 1
c         print*,'cregis gtneut.f called regeneration: ipart', IPART
         CALL KLKSREGEN
         
      ENDIF

#ifdef WithG4
      GO TO 600

C     ------------------------------------------------------------------


C     ------------------------------------------------------------------
C       Return secondary products to GEANT3 GKIN bank. 
C       Added by M.Harada(original is by Y.Ashida, R.Akutsu, Ka Ming
C     ------------------------------------------------------------------

C *** for the neutron-water interaction loaded from GEANT4
C      - process only limited particles 
C         + neutron   : secVec(1,ITER)==2112
C         + proton    : secVec(1,ITER)==2212
C         + electron  : secVec(1,ITER)==11
C         + positron  : secVec(1,ITER)==-11
C         + gamma-ray : secVec(1,ITER)==22
C         + deuteron  : secVec(1,ITER)==10020
C         + triton    : secVec(1,ITER)==10030
C         + 3He       : secVec(1,ITER)==20030
C         + alpha     : secVec(1,ITER)==20040
C

  590 continue
C  Change the way to accept seconadry for elastic process
C  Neutron is put back to original track instead of secondary
C  Occasionally reaches the maximum number of steps (MAXNST=10000), in this case put to seconadry
      numSecondary = 0
      ISTOP = 1
C  ZINTHA is number of interaction lengths to the next hadronic interaction
C  Properly initialize for next interaction
      ZINTHA=GARNDM(6)
      SLHADR=SLENG
      STEPHA=BIG
      DO ITER=1,100

        IF (secVec(1,ITER)==2112  .OR. secVec(1,ITER)==2212  .OR. !Neutron and proton 
     $      secVec(1,ITER)==11    .OR. secVec(1,ITER)==-11   .OR. !electron and positron 
     $      secVec(1,ITER)==22    .OR. secVec(1,ITER)==10020 .OR.  
     $      secVec(1,ITER)==10030 .OR. secVec(1,ITER)==20030 .OR. 
     $      secVec(1,ITER)==20040) THEN
C         Put the first neutron from elastic process to current track 
          IF (ITER==1.AND.secVec(1,ITER)==2112.AND. 
     $        NSTEP.LT.MAXNST .AND. 
     $        (processNumber==21.OR.processNumber==11)) THEN
               PP = SQRT(secVec(6,ITER)*secVec(6,ITER)+
     $                    secVec(7,ITER)*secVec(7,ITER)+
     $                    secVec(8,ITER)*secVec(8,ITER))
               VECT(4) = secVec(6,ITER)/PP
               VECT(5) = secVec(7,ITER)/PP
               VECT(6) = secVec(8,ITER)/PP
               VECT(7) = PP
               GEKIN = secVec(9,ITER)-
     $                 SQRT(secVec(9,ITER)*secVec(9,ITER)-PP*PP)
               GETOT = secVec(9,ITER)
               ISTOP = 0
C               write(*,*) 'ITER = ',ITER
C               write(*,*) 'Process = ',processNumber
C               write(*,*) 'Particle = NEUTRON' 
C               write(*,*) 'PP [GeV] = ',PP
C               write(*,*) 'GEKIN [GeV] = ',GEKIN
C               write(*,*) 'GETOT [GeV] = ',GETOT
C               write(*,*) 'TOFG = ',TOFG
          ELSE

               numSecondary = numSecondary + 1

               if (secVec(1,ITER)==2112) then 
                  GKIN(5,numSecondary) = 13
               else if (secVec(1,ITER)==2212) then 
                  GKIN(5,numSecondary) = 14
               else if (secVec(1,ITER)==11) then 
                  GKIN(5,numSecondary) = 3
               else if (secVec(1,ITER)==-11) then 
                  GKIN(5,numSecondary) = -3
               else if (secVec(1,ITER)==22) then 
                  GKIN(5,numSecondary) = 1 
               else if (secVec(1,ITER)==10020) then 
                  GKIN(5,numSecondary) = 45 
               else if (secVec(1,ITER)==10030) then 
                  GKIN(5,numSecondary) = 46
               else if (secVec(1,ITER)==20030) then 
                  GKIN(5,numSecondary) = 49 
               else if (secVec(1,ITER)==20040) then 
                  GKIN(5,numSecondary) = 47
               endif

               GKIN(1,numSecondary) = secVec(6,ITER)
               GKIN(2,numSecondary) = secVec(7,ITER)
               GKIN(3,numSecondary) = secVec(8,ITER)
               GKIN(4,numSecondary) = secVec(9,ITER)
               TOFD(numSecondary) = 0
               GPOS(1,numSecondary) = VECT(1)
               GPOS(2,numSecondary) = VECT(2)
               GPOS(3,numSecondary) = VECT(3)
c               write(*,*) 'ITER = ',ITER
c               write(*,*) 'Process = ',processNumber
c               write(*,*) 'Particle = ',secVec(1,ITER)
c               write(*,*) 'Energy [GeV] = ',secVec(9,ITER)
c          write(*,*) 'Momentum-X [GeV] = ',secVec(6,ITER)
c          write(*,*) 'Momentum-Y [GeV] = ',secVec(7,ITER)
c          write(*,*) 'Momentum-Z [GeV] = ',secVec(8,ITER)
          ENDIF

        ENDIF

      ENDDO
C     LMEC is the geant3 interaction code
C     Distinguish the capture, elastic and inleastic interactions
      NGKINE = numSecondary
      NMEC = NMEC + 1
      IF (processNumber==10 .OR. processNumber==20) THEN 
         LMEC(NMEC) = 18
         KCASE = NAMEC(18)
      ELSEIF (processNumber==11 .OR. processNumber==21) THEN 
         LMEC(NMEC) = 19
         KCASE = NAMEC(19)
      ELSEIF (processNumber==12 .OR. processNumber==22) THEN 
         LMEC(NMEC) = 20
         KCASE = NAMEC(20)
      ENDIF

C      write(*,*) '  NGKINE = ',NGKINE
C      write(*,*) ''
C      write(*,*) 'LMEC   = ',LMEC(NMEC)

C *** for the interactions other than neutron-water 
C     (or all the interactions in G3 mode) 
C
  600 continue

C      NMEC = NMEC + 1
C      LMEC(NMEC) = 12
C      IF (NGKINE > 0) THEN
C        DO ITER=1,NGKINE
C          if (GKIN(5,ITER).eq.45 .AND. ISTOP.eq.1 .AND. NUMED.eq.2) then
C            LMEC(NMEC) = 18
C            HCapCount = HCapCount + 1
C            write(*,*) 'HcaptureCount = ',HCapCount
C          endif
C        ENDDO
C      ENDIF

C     ------------------------------------------------------------------
#endif

*                                                             END GTNEUT
  999 END
