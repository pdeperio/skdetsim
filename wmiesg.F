************************************************************************
*     -------------------------------------
      FUNCTION WMIESG(BLAMDA)
*     -------------------------------------
*
*     (Purpose)
*       Get the Atteniation Coeff. of short wavelength in Water
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
*                     WMIESG = Asymmmetric Mie Scattering   
*
************************************************************************

      IMPLICIT NONE

      real COEFF, BLAMDA, WMIESG

#include "skdetsim_parm.h"
#include "rayfac.h"
      IF (WATER_TUNE .EQ. LE_TUNE) THEN
         WMIESG = miefac/BLAMDA**4
         RETURN
      END IF

C---- Unit of ALAMDA is Meter, so Convert it into Nano Meter

cc      BLAMDA=ALAMDA*1.0E+9
cc    print*,'blam',BLAMDA

C---- Get the Attenuation Coefficient in short wave length

C     ---- Wave Length =>  515NM

      IF(BLAMDA.GE.515.) THEN
        COEFF=0.

C      ELSE

C     ---- Wave Length =< 472 NM

       ELSE IF(BLAMDA.GE.472.0) THEN
  
         COEFF = 0.6015e23/BLAMDA**10.0

       ELSE 

C     ---- Wave Length =< 515 NM

C        COEFF = 2.332e8/BLAMDA**4
         COEFF = 7.148e14/BLAMDA**7.0          

      ENDIF

C---- Get the Attenuation Coefficiency

      WMIESG = COEFF

      RETURN
      END
