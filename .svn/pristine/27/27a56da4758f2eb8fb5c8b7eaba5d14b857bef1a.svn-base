************************************************************************
*     -------------------------------------
      FUNCTION WLNGSG(BLAMDA)
*     -------------------------------------
*
*     (Purpose)
*       Get the Atteniation Coeff. of long wavelength in Water
*
*     (Input)
*       BLAMDA : Wave Length of Cherencov Light ( nano Meter )
*
*     (Output)
*       WLNGSG  : Attenuation Coeff. ( 1/Meter )
*
*     (Creation Date and Author)
*       1993. 8.18 ; First version by Koshio
*       2001. 7.01 ; Modified by Shibata  (for ATMPD)
*       2003. 5.17 ; Modified by Yoo      (for LE)
*       2003.11.01 ; Modified by Parker   (for LE)
*       2006.4.18  : Updated for SK-III by K Bays
*
************************************************************************
      IMPLICIT NONE

#include "skdetsim_parm.h"
#include "skhead.h"

C--- for LE tune factor by Parker
      real long_absorption_factor

      integer i
      real blamda, wlngsg,coeff
      real A(21),B(21),POWER(21)

      DATA A/.0279,.0375,.0496,.0531,.0573,.0645,
     &       .0713,.0803,.107,.157,.241,.288,.306,.316,.327,
     &       .346,.407,.423,.445,.494,.638/

      DATA B/500.,510.,520.,530.,540.,550.,
     &       560.,570.,580.,590.,600.,610.,620.,630.,640.,
     &       650.,660.,670.,680.,690.,700./

      DATA POWER/14.933,14.933,14.402,3.580,4.072,
     &           6.451,5.563,6.716,16.505,22.429,25.498,10.779,
     &           3.728,2.010,2.173,3.643,10.635,2.564,3.422,7.155,
     &           17.778/


      real LEABS(4)
      real LECOE(5)

C---- Get the Attenuation Coefficient in long wave length

    
*============
C---- AP_TUNE
*============
      IF (WATER_TUNE .EQ. AP_TUNE) THEN
         IF(BLAMDA.LE.350.) THEN
            COEFF = (-0.268e-5)*BLAMDA+0.001569
         ELSE IF(BLAMDA.LE.B(1)) THEN
            IF (BLAMDA.LE.371.) THEN
               COEFF = (-0.268e-5)*BLAMDA+0.001569
            ELSE
               COEFF = ((BLAMDA/500.055)**12.303)*0.0224
            END IF
         ELSE IF(BLAMDA.GE.B(21)) THEN
            COEFF = (BLAMDA/B(21))**POWER(21)*A(21)
            coeff=coeff*(0.8)
         ELSE
            DO 100 I=2,21
               IF(BLAMDA.GT.B(I)) GOTO 100
               COEFF = (BLAMDA/B(I))**POWER(I)*A(I)
               coeff=coeff*(0.8)
               GOTO 101
 100        CONTINUE
 101        CONTINUE
         ENDIF
        

*============
*---- LE_TUNE
*============


      ELSE IF (WATER_TUNE.eq.LE_TUNE) then
c         print *,'checking le geometry',sk_geometry
         IF(SK_GEOMETRY.eq.SK_I) THEN
C----------------SK-I from ledetsim------------------------------
            IF(BLAMDA.LE.350.) THEN
               COEFF = 0.
            ELSE IF(BLAMDA.LE.B(1)) THEN
               COEFF = (BLAMDA/B(1))**POWER(1)*A(1)
            ELSE IF(BLAMDA.GE.B(21)) THEN
               COEFF = (BLAMDA/B(21))**POWER(21)*A(21)
            ELSE
               DO 190 I=2,21
                  IF(BLAMDA.GT.B(I)) GOTO 190
                  COEFF = (BLAMDA/B(I))**POWER(I)*A(I)
                  GOTO 191
 190           CONTINUE
 191           CONTINUE
            ENDIF
            coeff=coeff*0.8     ! should be tuned, 98-9-18 K.I
C----------------SK-I from ledetsim------------------------------

         ELSE IF(SK_GEOMETRY.eq.SK_II) then
*----------------SK-II Tuning------------------------------------
*  Here, we are temporally using the ABS which was taken 
*  at the almost end of SK-1. (run 010322-010327 2001/jun27-jul01)
*
            LEABS(1) = 0.50E8   ! == ABS(337nm) <-- take this as a reference
            LEABS(2) = 0.50E8   ! == ABS(371nm)
            LEABS(3) = 0.86E8   ! == ABS(400nm)
            LEABS(4) = 1.47E8   ! == ABS(420nm)

*---- SK-1 coefficients upto 500nm
            LECOE(1) = (350/B(1))**POWER(1)*A(1)
            LECOE(2) = (371/B(1))**POWER(1)*A(1)
            LECOE(3) = (400/B(1))**POWER(1)*A(1)
            LECOE(4) = (420/B(1))**POWER(1)*A(1)
            LECOE(5) = (500/B(1))**POWER(1)*A(1)
*
*  assume the wavelength dependence of LEABS kept for the moment  
* 
*
*  Note : This version of wlngsg.F (LE_TUNE) basically include only measured 
*  parameter points to avoid any unreasonalble corrections in the very 
*  long wavelength region.
*
************************************************************************

            IF (BLAMDA.LE.350.) THEN
               COEFF = 0. 
            ELSE IF(BLAMDA.LE.B(1)) THEN
               IF (BLAMDA.LE.371.) THEN

                  COEFF = LECOE(1)
     +                 +((LECOE(2)-LECOE(1))+(LEABS(2)-LEABS(1))/BLAMDA**4)
     /                 /(371-350)*(BLAMDA-350)

               ELSE IF (BLAMDA.LE.400.) THEN

                  COEFF = LECOE(2)+(LEABS(2)-LEABS(1))/BLAMDA**4
     +                 +((LECOE(3)-LECOE(2))+(LEABS(3)-LEABS(2))/BLAMDA**4)
     /                 /(400-371)*(BLAMDA-371)

               ELSE IF (BLAMDA.LE.420.) THEN

                  COEFF = LECOE(3)+(LEABS(3)-LEABS(1))/BLAMDA**4
     +                 +((LECOE(4)-LECOE(3))+(LEABS(4)-LEABS(3))/BLAMDA**4)
     /                 /(420-400)*(BLAMDA-400)

               ELSE
 
                  COEFF = LECOE(4) +(LEABS(4)-LEABS(1))/BLAMDA**4
     +                 +((LECOE(5)-LECOE(4))-LEABS(4)/BLAMDA**4)
     /                 /(500-420)*(BLAMDA-420)

               ENDIF
               
            ELSE IF(BLAMDA.GE.B(21)) THEN
               COEFF = (BLAMDA/B(21))**POWER(21)*A(21)
            ELSE
               DO 200 I=2,21
                  IF(BLAMDA.GT.B(I)) GOTO 200
                  COEFF = (BLAMDA/B(I))**POWER(I)*A(I)
                  GOTO 201
 200           CONTINUE
 201           CONTINUE
            ENDIF
            long_absorption_factor = 1.0
            coeff=coeff*long_absorption_factor
                          !End sk-II
*----------------SK-II Tuning------------------------------------
              
      ELSE IF(SK_GEOMETRY.eq.SK_III) then
*----------------SK-III Tuning------------------------------------
*  Temporarily using sk-II data, which itself was supposed to be temporary?

            LEABS(1) = 0.50E8   ! == ABS(337nm) <-- take this as a reference
            LEABS(2) = 0.50E8   ! == ABS(371nm)
            LEABS(3) = 0.86E8   ! == ABS(400nm)
            LEABS(4) = 1.47E8   ! == ABS(420nm)

*---- SK-1 coefficients upto 500nm
            LECOE(1) = (350/B(1))**POWER(1)*A(1)
            LECOE(2) = (371/B(1))**POWER(1)*A(1)
            LECOE(3) = (400/B(1))**POWER(1)*A(1)
            LECOE(4) = (420/B(1))**POWER(1)*A(1)
            LECOE(5) = (500/B(1))**POWER(1)*A(1)
*
*  assume the wavelength dependence of LEABS kept for the moment  
* 
*
*  Note : This version of wlngsg.F (LE_TUNE) basically include only measured 
*  parameter points to avoid any unreasonalble corrections in the very 
*  long wavelength region.
*
************************************************************************

            IF (BLAMDA.LE.350.) THEN
               COEFF = 0. 
            ELSE IF(BLAMDA.LE.B(1)) THEN
               IF (BLAMDA.LE.371.) THEN

                  COEFF = LECOE(1)
     +                 +((LECOE(2)-LECOE(1))+(LEABS(2)-LEABS(1))/BLAMDA**4)
     /                 /(371-350)*(BLAMDA-350)

               ELSE IF (BLAMDA.LE.400.) THEN

                  COEFF = LECOE(2)+(LEABS(2)-LEABS(1))/BLAMDA**4
     +                 +((LECOE(3)-LECOE(2))+(LEABS(3)-LEABS(2))/BLAMDA**4)
     /                 /(400-371)*(BLAMDA-371)

               ELSE IF (BLAMDA.LE.420.) THEN

                  COEFF = LECOE(3)+(LEABS(3)-LEABS(1))/BLAMDA**4
     +                 +((LECOE(4)-LECOE(3))+(LEABS(4)-LEABS(3))/BLAMDA**4)
     /                 /(420-400)*(BLAMDA-400)

               ELSE
 
                  COEFF = LECOE(4) +(LEABS(4)-LEABS(1))/BLAMDA**4
     +                 +((LECOE(5)-LECOE(4))-LEABS(4)/BLAMDA**4)
     /                 /(500-420)*(BLAMDA-420)

               ENDIF
               
            ELSE IF(BLAMDA.GE.B(21)) THEN
               COEFF = (BLAMDA/B(21))**POWER(21)*A(21)
            ELSE
               DO 300 I=2,21
                  IF(BLAMDA.GT.B(I)) GOTO 300
                  COEFF = (BLAMDA/B(I))**POWER(I)*A(I)
                  GOTO 301
 300           CONTINUE
 301           CONTINUE
            ENDIF
            long_absorption_factor = 1.0
            coeff=coeff*long_absorption_factor
                          !End sk-III
*----------------SK-III Tuning------------------------------------
           
      ELSE IF(SK_GEOMETRY.eq.SK_IV) then
*----------------SK-IV Tuning-------------------------------------
*
*   Future tuning for SK-IV, same as SK-III temtatively
*
*  Temporarily using sk-II data, which itself was supposed to be temporary?

            LEABS(1) = 0.50E8   ! == ABS(337nm) <-- take this as a reference
            LEABS(2) = 0.50E8   ! == ABS(371nm)
            LEABS(3) = 0.86E8   ! == ABS(400nm)
            LEABS(4) = 1.47E8   ! == ABS(420nm)

*---- SK-1 coefficients upto 500nm
            LECOE(1) = (350/B(1))**POWER(1)*A(1)
            LECOE(2) = (371/B(1))**POWER(1)*A(1)
            LECOE(3) = (400/B(1))**POWER(1)*A(1)
            LECOE(4) = (420/B(1))**POWER(1)*A(1)
            LECOE(5) = (500/B(1))**POWER(1)*A(1)
*
*  assume the wavelength dependence of LEABS kept for the moment  
* 
*
*  Note : This version of wlngsg.F (LE_TUNE) basically include only measured 
*  parameter points to avoid any unreasonalble corrections in the very 
*  long wavelength region.
*
************************************************************************

            IF (BLAMDA.LE.350.) THEN
               COEFF = 0. 
            ELSE IF(BLAMDA.LE.B(1)) THEN
               IF (BLAMDA.LE.371.) THEN

                  COEFF = LECOE(1)
     +                 +((LECOE(2)-LECOE(1))+(LEABS(2)-LEABS(1))/BLAMDA**4)
     /                 /(371-350)*(BLAMDA-350)

               ELSE IF (BLAMDA.LE.400.) THEN

                  COEFF = LECOE(2)+(LEABS(2)-LEABS(1))/BLAMDA**4
     +                 +((LECOE(3)-LECOE(2))+(LEABS(3)-LEABS(2))/BLAMDA**4)
     /                 /(400-371)*(BLAMDA-371)

               ELSE IF (BLAMDA.LE.420.) THEN

                  COEFF = LECOE(3)+(LEABS(3)-LEABS(1))/BLAMDA**4
     +                 +((LECOE(4)-LECOE(3))+(LEABS(4)-LEABS(3))/BLAMDA**4)
     /                 /(420-400)*(BLAMDA-400)

               ELSE
 
                  COEFF = LECOE(4) +(LEABS(4)-LEABS(1))/BLAMDA**4
     +                 +((LECOE(5)-LECOE(4))-LEABS(4)/BLAMDA**4)
     /                 /(500-420)*(BLAMDA-420)

               ENDIF
               
            ELSE IF(BLAMDA.GE.B(21)) THEN
               COEFF = (BLAMDA/B(21))**POWER(21)*A(21)
            ELSE
               DO 400 I=2,21
                  IF(BLAMDA.GT.B(I)) GOTO 400
                  COEFF = (BLAMDA/B(I))**POWER(I)*A(I)
                  GOTO 401
 400           CONTINUE
 401           CONTINUE
            ENDIF
            long_absorption_factor = 1.0
            coeff=coeff*long_absorption_factor

        ELSE IF(SK_GEOMETRY.eq.SK_V) then
            
C       copy from SK-IV

            IF (BLAMDA.LE.350.) THEN
               COEFF = 0. 
            ELSE IF(BLAMDA.LE.B(1)) THEN
               IF (BLAMDA.LE.371.) THEN

                  COEFF = LECOE(1)
     +                 +((LECOE(2)-LECOE(1))+(LEABS(2)-LEABS(1))/BLAMDA**4)
     /                 /(371-350)*(BLAMDA-350)

               ELSE IF (BLAMDA.LE.400.) THEN

                  COEFF = LECOE(2)+(LEABS(2)-LEABS(1))/BLAMDA**4
     +                 +((LECOE(3)-LECOE(2))+(LEABS(3)-LEABS(2))/BLAMDA**4)
     /                 /(400-371)*(BLAMDA-371)

               ELSE IF (BLAMDA.LE.420.) THEN

                  COEFF = LECOE(3)+(LEABS(3)-LEABS(1))/BLAMDA**4
     +                 +((LECOE(4)-LECOE(3))+(LEABS(4)-LEABS(3))/BLAMDA**4)
     /                 /(420-400)*(BLAMDA-400)

               ELSE
 
                  COEFF = LECOE(4) +(LEABS(4)-LEABS(1))/BLAMDA**4
     +                 +((LECOE(5)-LECOE(4))-LEABS(4)/BLAMDA**4)
     /                 /(500-420)*(BLAMDA-420)

               ENDIF
               
            ELSE IF(BLAMDA.GE.B(21)) THEN
               COEFF = (BLAMDA/B(21))**POWER(21)*A(21)
            ELSE
               DO 500 I=2,21
                  IF(BLAMDA.GT.B(I)) GOTO 500
                  COEFF = (BLAMDA/B(I))**POWER(I)*A(I)
                  GOTO 501
 500           CONTINUE
 501           CONTINUE
            ENDIF
            long_absorption_factor = 1.0
            coeff=coeff*long_absorption_factor
                          !End sk-IV            
         ENDIF                  ! end of sk_geometry flag
         
      ENDIF    ! end of WATER_TUNE flag
      



C---- Get the Attenuation Coefficiency

      WLNGSG = COEFF


      RETURN
      END
