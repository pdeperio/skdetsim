************************************************************************
*     -------------------------------------
      FUNCTION WTRSG(BLAMDA)
*     -------------------------------------
*
*     (Purpose)
*       Get the Atteniation Coefficiency in Water
*
*     (Input)
*       BLAMDA : Wave Length of Cherencov Light ( nano Meter )
*
*     (Output)
*       WTRSG  : Attenuation Coeff. ( 1/Meter )
*
*     (Creation Date and Author)
*       1993. 8.18 ; First version by Koshio
*       2001. 7.01 ; Modified by Shibata 
*       2007. 4.04 ; Modified by Yang
*
************************************************************************

      IMPLICIT NONE
#include "skhead.h"
#include "skdetsim_parm.h"
#include "./guwtr.h"

      REAL WTRSG,POPEFRYSG,WLNGSG
      REAL WABSSG,WRAYSG,WTMISG,WMIESG,WSRTSG
      REAL BLAMDA,PRELAMDA,PREVALUE
      REAL WSYMSCSG,WASYSCSG
      REAL COEFFM,COEFFL,COEFFS,COEFF

      DATA PRELAMDA,PREVALUE/0.,0./
      SAVE PRELAMDA,PREVALUE

      IF (BLAMDA .EQ. PRELAMDA) THEN
        WTRSG=PREVALUE
        RETURN
      END IF

C---- Just copy from SK-V for temporal number
      IF (SK_GEOMETRY.eq.SK_V ) THEN ! just copy from SK-IV

         WTRSG = WSYMSCSG(BLAMDA) + WASYSCSG(BLAMDA) + WABSSG(BLAMDA)

C---- Future tuning for SK-IV, same as SK-III temtatively                  
      else IF (SK_GEOMETRY.eq.SK_IV ) THEN

        WTRSG = WSYMSCSG(BLAMDA) + WASYSCSG(BLAMDA) + WABSSG(BLAMDA)

C---- Separate SK-III and before SK-III. By B. Yang
      else IF (SK_GEOMETRY.eq.SK_III ) THEN

C        WTRSG = WABSSG(BLAMDA) + WSYMSCSG(BLAMDA) + WASYSCSG(BLAMDA)
C---- To avoid unexpected very small difference between SK-I & II calculation form and above calculation form,
C---- addition order changed like following. The rounding off seems to give rise to the difference. 
        WTRSG = WSYMSCSG(BLAMDA) + WASYSCSG(BLAMDA) + WABSSG(BLAMDA)


      ELSE                      ! for SK-I and SK-II

C---- Get the Attenuation Coefficient in long wave length

        IF (WATER_TUNE .EQ. AP_TUNE) THEN
          COEFFL = WLNGSG(BLAMDA) ! absorption only
        ELSE IF (WATER_TUNE .EQ. LE_TUNE) THEN
           IF (SK_GEOMETRY.eq.SK_II ) THEN
              CALL SGPOPE_READ  ! TO GET POPE/FRY ATTENUATION COEFFICIENTS
              IF (BLAMDA.GT.414) THEN
                 COEFFL = POPEFRYSG(BLAMDA)
              ELSE
                 COEFFL = WLNGSG(BLAMDA)
              ENDIF
           ELSE IF (SK_GEOMETRY.eq.SK_I ) THEN
              COEFFL = WLNGSG(BLAMDA) ! absorption only
           ENDIF
        ENDIF

C---- Get the Attenuation Coefficient in short wave length

        IF (WATER_TUNE .EQ. AP_TUNE) THEN
C     ---- correct the Attenuation Coefficient in inverse meter
          COEFFL = COEFFL * (AABL+ASCL)
C     ---- correct the Attenuation Coefficient in inverse meter
C                ( Rejected short wave length absorption )
          COEFFS = WSRTSG(BLAMDA)
          COEFFS = (ARAS+APSS)*COEFFS
C---- Get the Attenuation Coefficient for MIE scattering
          COEFFM = WTMISG(BLAMDA)
C     ---- correct the Attenuation Coefficient for MIE scattering
C                  ( MIE Scatteing = Total Mie Scatteirng )
          COEFFM = (1.0)*COEFFM
        ELSE
          COEFFS = WABSSG(BLAMDA)+WRAYSG(BLAMDA)+WMIESG(BLAMDA)
          COEFFM=0
        END IF

C---- Get the Attenuation Length in meter

        COEFF = COEFFL + COEFFS + COEFFM
        WTRSG = COEFF
      ENDIF

      PRELAMDA=BLAMDA
      PREVALUE=WTRSG

      RETURN
      END
