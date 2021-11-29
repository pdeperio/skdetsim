*
* $Id: gdrele_h2o.F,v 1.3 2003-01-31 20:41:05 habig Exp $
*
* $Log: not supported by cvs2svn $
* Revision 1.2  2003/01/31 20:30:59  habig
* Fixed gcbank.h.  It had been redefined in the main program to make for a
* larger bank.  Changed it to be like the standard geant include with a
* larger bank, instead of abusing an equivalence.
*
* Updated all files which include gcbank.h to use geant321 headers and use
* implicit none.
*
* Revision 1.1  2003/01/10 23:01:11  kameda
*
* Copy from apdetsim module. No modification.
*
* Revision 1.2  2001/09/17 13:28:48  itow
*
* dedx parameters in a,m,x0,x1 are changed to PLB 26 version same as in EGS.
* Also bug in drele_he0.F is fixed
*
* Revision 1.1  2001/08/27 12:23:36  itow
*
* New dedx caclulation routines made by Kim are installed by Itow.
*
* Revision 1.1.1.1  1995/10/24 10:21:24  cernlib
* Geant
*
*
*CMZ :  3.21/02 29/03/94  15.41.21  by  S.Giani
*-- Author :
      SUBROUTINE GDRELE_H2O(EEL,CHARGE,JMA,DEDX)
C.
C.    ******************************************************************
C.    *                                                                *
C.    *       Compute ion losses for electron/positron                 *
C.    *                                                                *
C.    *    ==>Called by : GDRELA                                       *
C.    *       Author    G.Patrick *********                            *
C.    *                                                                *
C.    ******************************************************************
C. 2001.JUN.14 created from gdrele.F by H.Kim to get better dedx in water

      IMPLICIT NONE

CHK   NEW DE/DX PARAMETERS FOR ONLY WATER

#include "gcbank.h"
#include "geant321/gconsp.inc"
#include "geant321/gccuts.inc"

C local defines for implicit none
      integer jma, jprob
      real gam,eel,gam2,t,dedx,tcme,bet2,charge,y,d,d2,d3,d4,f,potl
      real fac,c,x0,x1,aa,x,del,cons

C Parameters of dedx for liquid water taken from PRB 26 6067(1082)
C Partialy taken from BNL-33571
C
      REAL WATLQ(11)
      data WATLQ /0.55509, 75.0, 1.00, 21.47, 2.203,3.502,0.2400
     &     ,2.5,0.2065,3.007, 0.097 /
CYI
C
      DATA CONS/0.153536E-3/
C.
C.    ------------------------------------------------------------------
C.
      JPROB=LQ(JMA-4)
C
      GAM=EEL/EMASS + 1.
      GAM2=GAM*GAM
      T=GAM-1.
      DEDX=0.
      IF(T.LE.0.)GO TO 99
      TCME=DCUTE/EMASS
      BET2=1.-1./GAM2
C     ------------------------------
      IF(CHARGE.GT.0.) THEN
         Y=1./(1.+GAM)
         D=TCME
         IF(T.LT.TCME) D=T
         D2=D*D/2.
         D3=2.*D2*D/3.
         D4=D2*D2
         F=LOG(T*D)-BET2*(T+2.*D-Y*(3.*D2
     *    +Y*(D-D3+Y*(D2-T*D3+D4))))/T
C
      ELSE
        D=TCME
         IF(T.LT.2.*TCME) D=0.5*T
         F=-1.-BET2+LOG((T-D)*D)+T/(T-D)
     *    +(0.5*D*D+(1.+2.*T)*LOG(1.-D/T))/GAM2
      ENDIF

CHK   NEW DENSITY EFFECT
C      POTL=Q(JPROB+11)
      POTL=LOG(WATLQ(2)*1.E-9/REAL(EMASS))
C      FAC=Q(JPROB+12)
      FAC=WATLQ(1)
C      C=Q(JPROB+13)
      C=WATLQ(6)
C      X0=Q(JPROB+14)
      X0=WATLQ(7)
C      X1=Q(JPROB+15)
      X1=WATLQ(8)
C      AA=Q(JPROB+16)
      AA=WATLQ(9)
C
      X=LOG(GAM2-1.)/4.6052
      DEL=0.
      IF(X.GT.X0)THEN
         DEL=4.6052*X-C
         IF(X.LE.X1)DEL=DEL+AA*(X1-X)**WATLQ(10)
      ENDIF
CHK
C
      DEDX=CONS*FAC*(LOG(2.*T+4.)-2.*POTL+F-DEL)/BET2
      IF(DEDX.LT.0.)DEDX=0.
C
  99  RETURN
      END

