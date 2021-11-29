************************************************************************
*     -------------------------------------
      FUNCTION WABSSG(BLAMDA)
*     -------------------------------------
*
*     (Purpose)
*       Get the Atteniation Coeff. of absorption of LE and SK-III ATMPD.
*       To support SK-I & II ATMPD absorption coeff, please call wlngsg.
*       Especially, this calculate the absorption coefficient of short
*       wavelength for SK-I&II LE. Also, pleae call wlngsg or popefrysg
*       for long wavelength part of SK-I&II LE absorption coefficient.
*       Cross point of wlngsg and popefrysg is defined wtrsg.F for SK-II
*
*       Called by wtrsg.F, sgabsc.F
*
*     (Input)
*       BLAMDA : Wave Length of Cherencov Light ( nano Meter )
*
*     (Output)
*       WABSSG  : Attenuation Coeff. of absorption ( 1/Meter )
*
*     (Creation Date and Author)
*       ????.??.?? ; First version by someone
*       2007. 4. 4   : the SK-III part of wlngsg.F and wtrsg.F merged and modified by B. Yang
*          The following log is from wlngsg.F
*           1993. 8.18 ; First version by Koshio
*           2001. 7.01 ; Modified by Shibata  (for ATMPD)
*           2003. 5.17 ; Modified by Yoo      (for LE)
*           2003.11.01 ; Modified by Parker   (for LE)
*           2006.4.18  : Updated for SK-III by K Bays
*
*       2007. 4. 5   : SK-III ATMPD absorption tuned by B. Yang
*       2007. 4.13   : SK-III LOWE absorption tuned by B. Yang
*       2007. 9.28   : SK-III water coefficient retuned due to bug about black sheet by B. Yang
*       2008. 3.14   : SK-III water parameter retuned with Feb, 08 data by B. Yang
*       2009. 5.16   : SK-III water parameter as a function of time by B.Yang
*       2009. 9.04   : SK-IV monthly water parameters installed by T.McLachlan
*                      Absorption function shape was set by Nov'08 data and
*                      modified by multiplicative factor for each month.
*       2012.12      : Changes by T. Akiri
*       2013. 8.15   : Add Gd-loaded effect by P.Weatherly
*
************************************************************************
      IMPLICIT NONE

#include "skdetsim_parm.h"
#include "skhead.h"
#include "./guwtr.h"
#include "rayfac.h"
#include "dsparm.h"


      INTEGER I
      REAL BLAMDA, COEFF,POPEFRYSG
      REAL A(21),B(21),POWER(21)
      REAL PRELAMDA,PREVALUE
      REAL WABSSG, WTGDEFF
      REAL BLAMDA0,VSAPFAC,SAPFAC,LAPFAC,VLAPFAC,SAPPOW
      REAL LEPPFRFAC,BLAMDASPL,BLAMDALNH,POPEFRYL,LEPFHPOW,LEPFHFAC
C      REAL BLAMDASP,LELONGFAC,LEPPFRFAC,LELONGPOW
C      INTEGER*4 PRETUNE, PREGEO
      
      DATA A/.0279,.0375,.0496,.0531,.0573,.0645,
     &     .0713,.0803,.107,.157,.241,.288,.306,.316,.327,
     &     .346,.407,.423,.445,.494,.638/

      DATA B/500.,510.,520.,530.,540.,550.,
     &     560.,570.,580.,590.,600.,610.,620.,630.,640.,
     &     650.,660.,670.,680.,690.,700./

      DATA POWER/14.933,14.933,14.402,3.580,4.072,
     &     6.451,5.563,6.716,16.505,22.429,25.498,10.779,
     &     3.728,2.010,2.173,3.643,10.635,2.564,3.422,7.155,
     &     17.778/

C---
C---BLAMDA0 : cross point of VERY SHORT and SHORT wave range. This is calculated from the following 
CC---VSAPFAC *( (-0.268e-5)*BLAMDA0+0.001569 ) = ((BLAMDA0/500.055)**12.303)*SAPFAC
C---VSAPFAC *( (-0.268e-5)*BLAMDA0+0.001569 ) = ((BLAMDA0/500.0)**SAPPOW)*SAPFAC
C---VSAPFAC : ATMPD ABSORPTION FACTOR under BLAMDA0. Tuend at Apr. 4. 2007
C---SAPFAC  : ATMPD ABSORPTION FACTOR under 500nm.
C---SAPPOW  : POWER OF ATMPD ABSORPTION Coefficient Function under 500nm over BLAMDA0.
C---LAPFAC  : ATMPD ABSORPTION FACTOR under 700nm. NOT tuned yet(SK-II value)
C---VLAPFAC : ATMPD ABSORPTION FACTOR over 700nm. NOT tuned yet(SK-II value)
C---But QE of PMT is simulated in 280nm and 700nm.
C---So in at least TRACK_QE_PHOTON mode, wavelength
C---under 280nm and over 700nm is meaningless.
C---
C        PARAMETER(BLAMDA0 = 395.993)
C        PARAMETER(VSAPFAC = 2.5)
C        PARAMETER(SAPFAC  = 0.0224)

      PARAMETER(BLAMDA0 = 376.783)
      PARAMETER(VSAPFAC = 2.673)
      PARAMETER(SAPFAC  = 0.02232)
      PARAMETER(SAPPOW  = 9.555)
      PARAMETER(LAPFAC  = 0.8)
      PARAMETER(VLAPFAC = 0.8)
C---
C---BLAMDASP  : cross point of popefry and TUNED LONG wavelength fitting curve
C---LELONGFAC : LE ABSORPTION FACTOR under BLAMDASP. Tuend at Apr. 4. 2007
C---LEPPFRFAC : Factor of LE Pope/Fry ABSORPTION over BLAMDASP. This is fixed to 1.0 
C---
C          PARAMETER(BLAMDASP  = 415.0)
C          PARAMETER(LELONGFAC = 2.623)	! tuned at Apr 07, data taken in Feb 07
C          PARAMETER(LELONGFAC = 3.177)	! tuned at Sep 07, data taken in Feb 07

      PARAMETER(LEPPFRFAC = 1.0 )

CC moved into rayfac.h or tr2fac.F, they became common variables except for LELONGPOW at Jun. 2008
CCCCCCC tuned at Mar 08, data taken in Feb 08
C          PARAMETER(BLAMDASP  = 463.918)
C          PARAMETER(LELONGFAC = 0.8013058)
C          PARAMETER(LELONGPOW = 10.6296)

CC New parameter introduced to connect the POPE/FRY's and tuned longer walength absorption coefficient function
      PARAMETER(BLAMDASPL = 450.003) ! cross point of popefry and TUNED LONG wavelength fitting curve at LELONGFAC=1.01279
      PARAMETER(BLAMDALNH = 489.997) ! cross point of popefry and TUNED LONG wavelength fitting curve at LELONGFAC=0.666415
      PARAMETER(POPEFRYL  = 0.009220563) ! =POPEFRYSG(BLAMDASPL)
      PARAMETER(LEPFHPOW  = 24.1328) ! power of POPE/FRY's around 510nm
      PARAMETER(LEPFHFAC  = 0.02250674207) ! =(BLAMDALNH/B(1))**LELONGPOW*A(1)

C---- Get the Attenuation Coefficient in long wave length

C      DATA PRELAMDA,PREVALUE,PRETUNE,PREGEO/0.,0.,0,0/
C      SAVE PRELAMDA,PREVALUE,PRETUNE,PREGEO
      DATA PRELAMDA,PREVALUE/0.,0./
      SAVE PRELAMDA,PREVALUE

C      IF (BLAMDA .EQ. PRELAMDA .AND. WATER_TUNE .EQ. PRETUNE
C     &    .AND. SK_GEOMETRY .EQ. PREGEO ) THEN

      IF (BLAMDA .EQ. PRELAMDA) THEN
         WABSSG = PREVALUE
         RETURN
      END IF

    
*============
C---- AP_TUNE
*============
C---------------------------------------------------------------------------
********* This function don't support SK-I&II ATMPD mode.
********* To support SK-I & II ATMPD absorption coeff, please call wlngsg
C---------------------------------------------------------------------------
      IF (WATER_TUNE .EQ. AP_TUNE) THEN

         IF(BLAMDA.LE.B(1)) THEN
            IF (BLAMDA.LE. BLAMDA0) THEN
               COEFF = ((-0.268e-5)*BLAMDA+0.001569)*VSAPFAC
            ELSE
C            COEFF = ((BLAMDA/500.055)**12.303)*SAPFAC
               COEFF = ((BLAMDA/500.0)**SAPPOW)*SAPFAC
            END IF
         ELSE IF(BLAMDA.GE.B(21)) THEN
            COEFF = (BLAMDA/B(21))**POWER(21)*A(21)*LAPFAC
         ELSE
            DO 100 I=2,21
               IF(BLAMDA.GT.B(I)) GOTO 100
               COEFF = (BLAMDA/B(I))**POWER(I)*A(I)*VLAPFAC
               GOTO 101
 100        CONTINUE
 101        CONTINUE
         ENDIF
         WABSSG = AABL*COEFF

*=======================
*---- LE_TUNE & UNI_TUNE
*=======================

C      ELSE IF (WATER_TUNE.eq.LE_TUNE) then
      ELSE
c         print *,'checking le geometry',sk_geometry
         IF(SK_GEOMETRY.eq.SK_I .OR. SK_GEOMETRY.eq.SK_II ) THEN
C-------------------------------------------------------------------------------------
********* At SK-I&II of LE, this give the absorption coefficient of short wavelength.
C-------------------------------------------------------------------------------------

            WABSSG = absfac/BLAMDA**4

         ELSE IF (SK_GEOMETRY.eq.SK_III) then
C------------------------------------------------------------------------------------
********* At SK-III of LE, this give the absorption coefficient of both short
********* and long wavelength. And here first setting for SK-III was removed here,
********* but those remain in wlngsg.F
C------------------------------------------------------------------------------------

               IF (BLAMDA.GT.BLAMDASP) THEN ! Mar 08

CCCCC This part is from wtrsg.F of SK-I&II 
                  CALL SGPOPE_READ
                  IF (BLAMDA.LT. BLAMDASPL) THEN ! Jun 08
                     COEFF = LEPPFRFAC * POPEFRYL
                  ELSE
                     COEFF = LEPPFRFAC * POPEFRYSG(BLAMDA)
                  ENDIF

CC            COEFF = LEPPFRFAC * POPEFRYSG(BLAMDA)

               ELSEIF(BLAMDA.LE.BLAMDALNH) THEN	! Jun 08
C Mar 08            COEFF = LELONGFAC * (BLAMDA/B(1))**POWER(1)*A(1)
                  COEFF = LELONGFAC * (BLAMDA/B(1))**LELONGPOW*A(1)

               ELSE             ! Jun08
                  COEFF = LELONGFAC * (BLAMDA/BLAMDALNH)**LEPFHPOW * LEPFHFAC
               ENDIF

          ELSE IF(SK_GEOMETRY.eq.SK_IV) then

C------------------------------------------------------------------------------------
********* SK-IV is similar to the SK-III, so now it's temporal.
********* tuned at Mar 08, data taken in Feb 08
C------------------------------------------------------------------------------------

C--       Tuned by Thomas McLachlan 2009 August
C--       See tr2fac.F for tuned variables

               BLAMDASP  = 463.918  

               IF (BLAMDA.GT.BLAMDASP) THEN ! Mar 08
                  CALL SGPOPE_READ
                  COEFF = LEPPFRFAC * POPEFRYSG(BLAMDA)
               ELSE
                  COEFF = absfit * abssfac * A(1) * (BLAMDA/500.0)**abspow
               ENDIF

          ELSE IF(SK_GEOMETRY.eq.SK_V) then
c     Just copy from SK-IV J.Kameda 2019/10/29
             
               BLAMDASP  = 474.75

               IF (BLAMDA.GT.BLAMDASP) THEN ! Mar 08
                  CALL SGPOPE_READ
                  COEFF = LEPPFRFAC * POPEFRYSG(BLAMDA)
               ELSE
                  COEFF = absfit * abssfac * A(1) * (BLAMDA/500.0)**abspow
               ENDIF
          ELSE 
               
           WRITE(0,*) 'ERROR : SK_GEOMETRY is not usual,',SK_GEOMETRY 
           STOP               

          ENDIF

#ifndef WTSCTT
         if (SK_GEOMETRY .ge. SK_IV) then
            WABSSG = absfit*absfac/BLAMDA**4 + COEFF

            if(IWATERP.eq.201208) then
* Tarek Akiri, Dec. 2012
* New function with parameter from 16/04/2009
* Use Pope-Fry data from 480 nm
               if(BLAMDA .LT. 480) then
                  WABSSG = (6.96317e-27)*(BLAMDA**9.01) + 182.9/BLAMDA**2
               else 
                  WABSSG = POPEFRYSG(BLAMDA)
               endif
            endif
         else
            WABSSG = absfac/BLAMDA**4 + COEFF
         endif
#else
         if (SK_GEOMETRY .eq. SK_V) then !just copy from SK-IV
            WABSSG = absfit*absfac/BLAMDA**4 + COEFF*ablfac
         elseif (SK_GEOMETRY .eq. SK_IV) then
            WABSSG = absfit*absfac/BLAMDA**4 + COEFF*ablfac            
         else
            WABSSG = absfac/BLAMDA**4 + COEFF*ablfac
         endif
#endif

C      ELSE
C        WRITE(0,*) 'ERROR : WATER_TUNE is not usual,', WATER_TUNE
C        STOP
      ENDIF


C---- Apply Gd-loaded effect by P.Weatherly

      if(WATER_TUNE.eq.11) then
         WABSSG = WABSSG + WTGDEFF(BLAMDA)
      endif

C---- Get the Attenuation Coefficiency

      PRELAMDA = BLAMDA
      PREVALUE = WABSSG
C      PRETUNE = WATER_TUNE
C      PREGEO = SK_GEOMETRY

      RETURN
      END
