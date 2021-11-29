************************************************************************
*     -------------------------------------
      FUNCTION WASYSCSG(BLAMDA)
*     -------------------------------------
*
*     (Purpose)
*       Get the Atteniation Coefficiency of asymmetric scattering in Water
*       Scattered direction is only forward. Related angle distribution is
*       simulated by sgmies.F. This is called at only SK-III mode.
*
*     (Input)
*       BLAMDA : Wave Length of Cherencov Light ( nano Meter )
*
*     (Output)
*       WASYSCSG  : Attenuation Coeff. of asymmetric scattering ( 1/Meter )
*                     Asymmetric scattering = Asymmetric mie scattering 
*
*     (Creation Date and Author)
*
*       Apr. 5. 2007 ; First version by B. Yang
*       Apr. 8. 2007 ; Meaningless Safe guard commented out by B. Yang
*       Sep. 4. 2009 ; SK-IV monthly water parameters installed by T.McLachlan
*                      Asymmetric function shape was set by Apr'09 data and
*                      modified by multiplicative factor for each month.
*
************************************************************************

      IMPLICIT NONE
#include "skhead.h"
#include "skdetsim_parm.h"
#include "./guwtr.h"
#include "./rayfac.h"
#include "dsparm.h"

      REAL BLAMDA,PRELAMDA,PREVALUE
C      INTEGER*4 PRETUNE
      REAL WASYSCSG

C      DATA PRELAMDA,PREVALUE,PRETUNE/0.,0.,0/
C      SAVE PRELAMDA,PREVALUE,PRETUNE
      DATA PRELAMDA,PREVALUE/0.,0./
      SAVE PRELAMDA,PREVALUE

C      IF (BLAMDA .EQ. PRELAMDA .AND. WATER_TUNE .EQ. PRETUNE) THEN
      IF (BLAMDA .EQ. PRELAMDA) THEN
        WASYSCSG = PREVALUE
        RETURN
      END IF
      
      IF(WATER_TUNE.EQ.AP_TUNE) THEN
        WASYSCSG = AMIS/BLAMDA
      ELSE

         if (SK_GEOMETRY .eq. SK_V) then ! copy from SK-IV

          WASYSCSG = asyfit*(1+ (miefac/BLAMDA**4)*(BLAMDA-miesfac)**2)
            
         else if (SK_GEOMETRY .eq. SK_IV) then

C        SKIV Tuned by month. See tr2fac.F for tuned variables. Thomas McLachlan.

          WASYSCSG = asyfit*(1+ (miefac/BLAMDA**4)*(BLAMDA-miesfac)**2)

          if(IWATERP.eq.201208) then
* Tarek Akiri, Dec. 2012
* New function with parameter from 16/04/2009
             WASYSCSG = 387.013/(BLAMDA**2.5)
          endif
         else

          WASYSCSG = miefac/BLAMDA**4
         endif

      ENDIF

      PRELAMDA = BLAMDA
      PREVALUE = WASYSCSG
C      PRETUNE = WATER_TUNE

      RETURN
      END
