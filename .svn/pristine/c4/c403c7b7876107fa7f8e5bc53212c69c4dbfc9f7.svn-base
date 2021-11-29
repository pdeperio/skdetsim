*
* $Id: gtnino.F,v 1.1.1.1 1995/10/24 10:21:44 cernlib Exp $
*
* $Log: gtnino.F,v $
* Revision 1.1.1.1  1995/10/24 10:21:44  cernlib
* Geant
*
*
#ifndef LASER
#include "geant321/pilot.h"
#endif
*CMZ :  3.21/04 01/11/94  16.25.59  by  S.Ravndal
*-- Author :
      SUBROUTINE GTNINO
C.
C.    ******************************************************************
C.    *                                                                *
C.    *       Tracking routine for "GEANTINO" type tracks              *
C.    *                                                                *
C.    *    ==>Called by : GTRACK                                       *
C.    *       Author    R.Brun, F.Bruyant                              *
C.    *                                                                *
C.    ******************************************************************
C.

#ifndef LASER

#include "geant321/gconsp.inc"
#include "geant321/gctmed.inc"
#include "geant321/gctrak.inc"
#if defined(CERNLIB_DEBUG)
#include "geant321/gcunit.inc"
#endif
#if !defined(CERNLIB_SINGLE)
      PARAMETER (EPSMAC=1.E-6)
#endif
#if defined(CERNLIB_SINGLE)
      PARAMETER (EPSMAC=1.E-11)
#endif
C.
C.    ------------------------------------------------------------------
C.
      STEP = STEMAX
      IPROC = 103
      CALL GTNEXT
      IF (IGNEXT.NE.0) THEN
         STEP = SNEXT + PREC
         IPROC=0
         INWVOL= 2
         NMEC =  1
         LMEC(1)=1
*
* *** Linear transport
*
         DO 10 I = 1,3
            VECTMP  = VECT(I) +STEP*VECT(I+3)
            IF(VECTMP.EQ.VECT(I)) THEN
*
* *** Correct for machine precision
*
               IF(VECT(I+3).NE.0.) THEN
                  VECTMP = VECT(I)+ABS(VECT(I))*SIGN(1.,VECT(I+3))*
     +            EPSMAC
*                 IF(NMEC.GT.0) THEN
*                    IF(LMEC(NMEC).EQ.32) NMEC=NMEC-1
*                 ENDIF
*                 NMEC=NMEC+1
*                 LMEC(NMEC)=32
#if defined(CERNLIB_DEBUG)
                  WRITE(CHMAIL, 10000)
                  CALL GMAIL(0,0)
                  WRITE(CHMAIL, 10100) GEKIN, NUMED, STEP, SNEXT
                  CALL GMAIL(0,0)
10000 FORMAT(' Boundary correction in GTNINO: ',
     +       '    GEKIN      NUMED       STEP      SNEXT')
10100 FORMAT(31X,E10.3,1X,I10,1X,E10.3,1X,E10.3,1X)
#endif
               ENDIF
            ENDIF
            VECT(I) = VECTMP
   10    CONTINUE
      ELSE
         DO 20 I = 1,3
            VECT(I)  = VECT(I) +STEP*VECT(I+3)
   20    CONTINUE
      ENDIF
      SLENG = SLENG +STEP
*
* *** Update time of flight
*
      TOFG = TOFG +STEP/CLIGHT
*
      IF(IPROC.NE.0) THEN
         NMEC=NMEC+1
         LMEC(NMEC)=IPROC
      ENDIF

#else

      DOUBLE PRECISION PI,TWOPI,PIBY2,DEGRAD,RADDEG,CLIGHT,BIG,EMASS
      DOUBLE PRECISION EMMU,PMASS,AVO
*
      PARAMETER (PI=3.14159265358979324D0)
      PARAMETER (TWOPI=6.28318530717958648D0)
      PARAMETER (PIBY2=1.57079632679489662D0)
      PARAMETER (DEGRAD=0.0174532925199432958D0)
      PARAMETER (RADDEG=57.2957795130823209D0)
      PARAMETER (CLIGHT=29979245800.D0)
      PARAMETER (BIG=10000000000.D0)
      PARAMETER (EMASS=0.0005109990615D0)
      PARAMETER (EMMU=0.105658387D0)
      PARAMETER (PMASS=0.9382723128D0)
      PARAMETER (AVO=0.60221367D0)
*
      COMMON/GCTMED/NUMED,NATMED(5),ISVOL,IFIELD,FIELDM,TMAXFD,STEMAX
     +      ,DEEMAX,EPSIL,STMIN,CFIELD,PREC,IUPD,ISTPAR,NUMOLD
      COMMON/GCTLIT/THRIND,PMIN,DP,DNDL,JMIN,ITCKOV,IMCKOV,NPCKOV
C
      PARAMETER (MAXMEC=30)
      COMMON/GCTRAK/VECT(7),GETOT,GEKIN,VOUT(7),NMEC,LMEC(MAXMEC)
     + ,NAMEC(MAXMEC),NSTEP ,MAXNST,DESTEP,DESTEL,SAFETY,SLENG
     + ,STEP  ,SNEXT ,SFIELD,TOFG  ,GEKRAT,UPWGHT,IGNEXT,INWVOL
     + ,ISTOP ,IGAUTO,IEKBIN, ILOSL, IMULL,INGOTO,NLDOWN,NLEVIN
     + ,NLVSAV,ISTORY
      PARAMETER (MAXME1=30)
      COMMON/GCTPOL/POLAR(3), NAMEC1(MAXME1)
C
      PARAMETER (EPSMAC=1.E-6)
C.
C#include "geant321/gcflag.inc"
      COMMON/GCFLAG/IDEBUG,IDEMIN,IDEMAX,ITEST,IDRUN,IDEVT,IEORUN
     +        ,IEOTRI,IEVENT,ISWIT(10),IFINIT(20),NEVENT,NRNDM(2)
      COMMON/GCFLAX/BATCH, NOLOG
      LOGICAL BATCH, NOLOG
C.    ------------------------------------------------------------------
C.
 
C      write(*,*) '&&&&&&&&& GTNINO &&&&&&&&&&&&&&&&', ievent
      ISTOP = 2
      CALL LSRGEN
C      RETURN

#endif

      END
