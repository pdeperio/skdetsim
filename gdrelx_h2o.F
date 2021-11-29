CDECK  ID>, GDRELX.
*CMZ :  3.21/03 06/10/94  18.22.43  by  S.Giani
*-- Author :
      SUBROUTINE GDRELX_H2O(A,Z,DENS,T,HMASS,DEDX)
C.
C.    ******************************************************************
C.    *                                                                *
C.    *  Calculates the mean 1/DENS*dE/dx of a particle with kinetic   *
C.    *  energy T in an element of atomic number Z, atomic weight A    *
C.    *  and density DENS ( the density is just used for the           *
C.    *  calculation of the density effect in the case of high T).     *
C.    *  The routine reproduces the experimental and/or tabulated      *
C.    *  energy losses rather well down to T -> 0.                     *
C.    *  Simple parametrization is used for  T .le. T2L=2 MeV (see     *
C.    *  H.H.Andersen,J.F.Ziegler:Hydrogen stopping powers and         *
C.    *  ranges in all element,Pergamon Press,1977.).                  *
C.    *  For T .gt. T2L=2 MeV the corrected Bethe-Bloch stopping       *
C.    *  power / restricted energy loss formula is used.               *
C.    *                                                                *
C.    *                                                                *
C.    *    ==>Called by : GDRELA                                       *
C.    *       Author    L.Urban    *********                           *
C.    *                                                                *
C.    ******************************************************************
C. 2001.JUN.14 created from gdrelx.F by H.Kim to get better dedx
      IMPLICIT NONE
#include "geant321/gconsp.inc"
#include "geant321/gccuts.inc"
#include "geant321/gcunit.inc"
#include "gcbank.h"
#include "geant321/gcjloc.inc"
C local defs for implicit none
      integer nlmat,i,j,iz
      real z,fac,a,t1lim,hmass,t2lim,t,p,e,bet2,eta,eta2,b2g2,tmax,poti
      real cc,xa,xm,x0,x1,aa,x,delta,potsq
      real f1,f2,f3,f4,f5,ce,tupp,dedx,tmpwght,dens,sbb,corbb
      real emprot
      real amukev,d,t1l,t2l,b,c,cecof,eps

CHK      PARAMETER (AMUKEV=931494.32,D=0.000153537,T1L=0.00001,T2L=0.002)
      PARAMETER (AMUKEV=931494.32,D=0.000153536,T1L=0.00001,T2L=0.002)
      PARAMETER (EMPROT=0.9382723)
      DIMENSION B(6,92),C(6,92),CECOF(6)
*
      DATA ((B(I,J),I=1,6),J=1,19) /
     +     1.262,1.44,242.6,12000.,0.1159,18.8,
     +     1.229,1.397,484.5,5873.,0.05225,41.7,
     +     1.411,1.6,725.6,3013.,0.04578,47.6,
     +     2.248,2.59,966.,153.8,0.03475,62.7,
     +     2.474,2.815,1206.,1060.,0.02855,76.0,
     +     2.631,2.989,1445.,957.2,0.02819,77.3,
     +     2.954,3.35,1683.,1900.,0.02513,86.7,
     +     2.652,3.,1920.,2000.,0.0223,97.7,
     +     2.085,2.352,2157.,2634.,0.01816,120.,
     +     1.951,2.199,2393.,2699.,0.01568,139.,
     +     2.542,2.869,2628.,1854.,0.01472,148.,
     +     3.792,4.293,2862.,1009.,0.01397,156.,
     +     4.154,4.739,2766.,164.5,0.02023,162.,
     +     4.15,4.7,3329.,550.,0.01321,165.,
     +     3.232,3.647,3561.,1560.,0.01267,172.,
     +     3.447,3.891,3792.,1219.,0.01211,180.,
     +     5.047,5.714,4023.,878.6,0.01178,185.,
     +     5.731,6.5,4253.,530.,0.01123,194.,
     +     5.151,5.833,4482.,545.7,0.01129,193./
      DATA ((B(I,J),I=1,6),J=20,38) /
     +     5.521,6.252,4710.,553.3,0.01112,196.,
     +     5.201,5.884,4938.,560.9,0.009995,218.,
     +     4.862,5.496,5165.,568.5,0.009474,230.,
     +     4.48,5.055,5391.,952.3,0.009117,239.,
     +     3.983,4.489,5616.,1336.,0.008413,259.,
     +     3.469,3.907,5725.,1461.,0.008829,270.,
     +     3.519,3.963,6065.,1243.,0.007782,280.,
     +     3.14,3.535,6288.,1372.,0.007361,296.,
     +     3.553,4.004,6205.,555.1,0.008763,310.,
     +     3.696,4.175,4673.,387.8,0.02188,322.,
     +     4.21,4.75,6953.,295.2,0.006809,320.,
     +     5.041,5.697,7173.,202.6,0.006725,324.,
     +     5.554,6.3,6496.,110.,0.009689,330.,
     +     5.323,6.012,7611.,292.5,0.006447,338.,
     +     5.874,6.656,7395.,117.5,0.007684,340.,
     +     5.611,6.335,8046.,365.2,0.006244,349.,
     +     6.411,7.25,8262.,220.,0.006087,358.,
     +     5.694,6.429,8478.,292.9,0.006087,358.,
     +     6.339,7.159,8693.,330.3,0.006003,363./
      DATA ((B(I,J),I=1,6),J=39,57)  /
     +     6.407,7.234,8907.,367.8,0.005889,370.,
     +     6.734,7.603,9120.,405.2,0.005765,378.,
     +     6.902,7.791,9333.,442.7,0.005587,390.,
     +     6.425,7.248,9545.,480.2,0.005367,406.,
     +     6.799,7.671,9756.,517.6,0.005315,410.,
     +     6.108,6.887,9966.,555.1,0.005151,423.,
     +     5.924,6.677,10180.,592.5,0.004919,443.,
     +     5.238,5.9,10380.,630.,0.004758,458.,
     +     5.623,6.354,7160.,337.6,0.01394,466.,
     +     5.814,6.554,10800.,355.5,0.004626,471.,
     +     6.23,7.024,11010.,370.9,0.00454,480.,
     +     6.41,7.227,11210.,386.4,0.004474,487.,
     +     7.5,8.48,8608.,348.,0.009074,494.,
     +     6.979,7.871,11620.,392.4,0.004402,495.,
     +     7.725,8.716,11830.,394.8,0.004376,498.,
     +     8.231,9.289,12030.,397.3,0.004384,497.,
     +     7.287,8.218,12230.,399.7,0.004447,490.,
     +     7.899,8.911,12430.,402.1,0.004511,483.,
     +     8.041,9.071,12630.,404.5,0.00454,480./
      DATA ((B(I,J),I=1,6),J=58,76)  /
     +     7.489,8.444,12830.,406.9,0.00442,493.,
     +     7.291,8.219,13030.,409.3,0.004298,507.,
     +     7.098,8.,13230.,411.8,0.004182,521.,
     +     6.91,7.786,13430.,414.2,0.004058,537.,
     +     6.728,7.58,13620.,416.6,0.003976,548.,
     +     6.551,7.38,13820.,419.,0.003877,562.,
     +     6.739,7.592,14020.,421.4,0.003863,564.,
     +     6.212,6.996,14210.,423.9,0.003725,585.,
     +     5.517,6.21,14400.,426.3,0.003632,600.,
     +     5.219,5.874,14600.,428.7,0.003498,623.,
     +     5.071,5.706,14790.,433.,0.003405,640.,
     +     4.926,5.542,14980.,433.5,0.003342,652.,
     +     4.787,5.386,15170.,435.9,0.003292,662.,
     +     4.893,5.505,15360.,438.4,0.003243,672.,
     +     5.028,5.657,15550.,440.8,0.003195,682.,
     +     4.738,5.329,15740.,443.2,0.003186,684.,
     +     4.574,5.144,15930.,442.4,0.003144,693.,
     +     5.2,5.851,16120.,441.6,0.003122,698.,
     +     5.07,5.704,16300.,440.9,0.003082,707./
      DATA ((B(I,J),I=1,6),J=77,92)  /
     +     4.945,5.563,16490.,440.1,0.002965,735.,
     +     4.476,5.034,16670.,439.3,0.002871,759.,
     +     4.856,5.46,18320.,438.5,0.002542,755.,
     +     4.308,4.843,17040.,487.8,0.002882,756.,
     +     4.723,5.311,17220.,537.,0.002913,748.,
     +     5.319,5.982,17400.,586.3,0.002871,759.,
     +     5.956,6.7,17800.,677.,0.00266,765.,
     +     6.158,6.928,17770.,586.3,0.002812,775.,
     +     6.204,6.979,17950.,586.3,0.002776,785.,
     +     6.181,6.954,18120.,586.3,0.002748,793.,
     +     6.949,7.82,18300.,586.3,0.002737,796.,
     +     7.506,8.448,18480.,586.3,0.002727,799.,
     +     7.649,8.609,18660.,586.3,0.002697,808.,
     +     7.71,8.679,18830.,586.3,0.002641,825.,
     +     7.407,8.336,19010.,586.3,0.002603,837.,
     +     7.29,8.204,19180.,586.3,0.002573,847./
      DATA C/92*0.,92*0.,92*0.,92*0.,92*0.,92*0./
      DATA CECOF/0.42237,0.0304,-0.00038,3.858,-0.1668,0.00158/
      DATA EPS/0.000001/

CYI  New dE/dx parameters for lq. water taken from R.M.Sternheimer et al
CYI  11 parameters mean
CYI  Z/A, I(eV), density, plazma freq., rho, -C, X0, X1, a, m, delta_max
C Parameters of dedx for liquid water taken from PRB 26 6067(1082)
C Partialy taken from BNL-33571
C
      REAL WATLQ(11)
      data WATLQ /0.55509, 75.0, 1.00, 21.47, 2.203,3.502,0.2400
     &     ,2.5,0.2065,3.007, 0.097 /
CYI
CHK
      REAL TMPDEDX,MIDDEDX,TMPAA,TMPZZ,TMPWGH,TMPCE
      INTEGER L
CHK
*     ------------------------------------------------------------------
*      in the case of non-integer Z the low energy parameters
*      and the ionization potential are taken at INT(Z) !
*
      IZ=INT(Z+EPS)
      IF((IZ.LT.1).OR.(IZ.GT.92)) GOTO 10
*
*     Calculate coefficients C(I,J) if it has not been done already
*
      IF(C(1,IZ).LT.EPS) THEN
         FAC=AVO/A
         C(1,IZ)=FAC*AMUKEV**0.5*B(1,IZ)
         C(2,IZ)=FAC*AMUKEV**0.45*B(2,IZ)
         C(3,IZ)=FAC*B(3,IZ)/AMUKEV
         C(4,IZ)=B(4,IZ)/AMUKEV
         C(5,IZ)=AMUKEV*B(5,IZ)
*
*                     POTI=16.E-9*Z**0.9
*
CHK         POTI=B(6,IZ)*1.E-9
         
*
*
C         C(6,IZ)=POTI
         C(6,IZ)=WATLQ(2)*1.E-9
CHK
*
      ENDIF
*     ----------------------------------------------------------------
      T1LIM=HMASS*T1L/EMPROT
      T2LIM=HMASS*T2L/EMPROT
*
*     Calculate dE/dx
*     ---> for T .le. T1LIM (very low energy)
*
CHK
      IF(T.LE.T2LIM) THEN
*
*     ---> for T .gt. T2LIM ( "high " energy , Bethe-Bloch formula)
*
         print *,'ERROR OF T VALUE IN GDRELX_H2O.F'
         stop
      ELSE
CHK
            P=SQRT(T*(T+2.*HMASS))
            E=T+HMASS
            BET2=(P/E)**2
            ETA=P/HMASS
            ETA2=ETA*ETA
*+++ new line follows.....
            B2G2=ETA*ETA
*+++ end of correction
            TMAX=2.*EMASS*T*(T+2.*HMASS)
*+++  correction of the next line
*           TMAX=TMAX/(HMASS**2+EMASS**2+EMASS*(T+HMASS))
            TMAX=TMAX/(HMASS**2+EMASS**2+2.*EMASS*E)
*+++ end of correction
*
*         density correction
*
CHK
            POTI=WATLQ(2)*1.E-9
            CC=WATLQ(6)
            XA=CC/4.6052
            XM=WATLQ(10)
            X0=WATLQ(7)
            X1=WATLQ(8)
            AA=WATLQ(9)
            X=LOG10(ETA)
            DELTA=0.
            IF(X.GT.X0) THEN
               DELTA=4.6052*X-CC
               IF(X.LT.X1) DELTA=DELTA+AA*(X1-X)**XM
            ENDIF
CHK
*
*         calculate shell correction
*
            POTSQ=POTI*POTI
            IF(ETA.GT.0.13) THEN
               F1=1./ETA2
               F2=F1*F1
               F3=F1*F2
               F4=(F1*CECOF(1)+F2*CECOF(2)+F3*CECOF(3))*1.E+12
               F5=(F1*CECOF(4)+F2*CECOF(5)+F3*CECOF(6))*1.E+18
               CE=F4*POTSQ+F5*POTSQ*POTI
            ELSE
               ETA2=0.0169
               F1=1./ETA2
               F2=F1*F1
               F3=F1*F2
               F4=(F1*CECOF(1)+F2*CECOF(2)+F3*CECOF(3))*1.E+12
               F5=(F1*CECOF(4)+F2*CECOF(5)+F3*CECOF(6))*1.E+18
               CE=F4*POTSQ+F5*POTSQ*POTI
               CE=CE*LOG(T/T2LIM)/LOG(0.0079/T2LIM)
            ENDIF
*
            F1=D*WATLQ(1)/BET2
*
*         stopping power or restricted dE/dx ?
*
*+++  correction of the next few lines

            TUPP=DCUTM
            IF(TMAX.LT.DCUTM) TUPP=TMAX
            F2=LOG(2.*EMASS*B2G2/POTI)+LOG(TUPP/POTI)
     +         -BET2*(1.+TUPP/TMAX)
*+++ end of correction
            DEDX=F1*(F2-DELTA-2.*CE/Z)
c            DEDX=F1*(F2-DELTA)
*
*
CHK REPLACEMENT BASED ON OLD METHODS(CONTINUOUS CONDITION AT T2LIM)
C
C            TAU=T2LIM/HMASS
C            SL=C(2,IZ)*TAU**0.45
C            SH=C(3,IZ)*LOG(1.+C(4,IZ)/TAU+C(5,IZ)*TAU)/TAU
C            MIDDEDX=SL*SH/(SL+SH)
         NLMAT = Q(JMA+11)
         NLMAT  = IABS(NLMAT)
         MIDDEDX=0.
         TMPDEDX=0.
         DO 20 L=1,NLMAT
            TMPAA=Q(JMIXT+L)
            TMPZZ=Q(JMIXT+NLMAT+L)
            TMPWGHT=Q(JMIXT+2*NLMAT+L)
         CALL GDRELX(TMPAA,TMPZZ,DENS*TMPWGHT,T2LIM,REAL(EMMU),TMPDEDX)
            IF(TMPDEDX.LT.0.)TMPDEDX=0.
            MIDDEDX=MIDDEDX+TMPWGHT*TMPDEDX
            
   20    CONTINUE

*
            TMAX=2.*EMASS*T2LIM*(T2LIM+2.*HMASS)
*+++  correction of the next line
*           TMAX=TMAX/(HMASS**2+EMASS**2+EMASS*(T2LIM+HMASS))
            TMAX=TMAX/(HMASS**2+EMASS**2+2.*EMASS*(T2LIM+HMASS))
*+++  end of correction
            BET2=T2LIM*(T2LIM+2.*HMASS)/(T2LIM+HMASS)**2
            SBB=2.*(LOG(TMAX/POTI)-BET2)
            SBB=D*WATLQ(1)*SBB/BET2
C            CORBB=(ST/SBB-1.)*T2LIM
            CORBB=(MIDDEDX/SBB-1.)*T2LIM
*
            DEDX=DEDX*(1.+CORBB/T)
CHK
*
      ENDIF

*
      RETURN
   10 WRITE (CHMAIL,10000) Z
      CALL GMAIL(1,1)
10000 FORMAT(5X,'***GDRELP  (Z.LT.1.).OR.(Z.GT.92.) ==> Z=',F10.1)

      END





