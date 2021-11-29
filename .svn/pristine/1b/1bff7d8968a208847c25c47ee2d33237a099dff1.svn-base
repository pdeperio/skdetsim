************************************************************************
*     -------------------------------------
      FUNCTION WRAYSG(BLAMDA)
*     -------------------------------------
*
*     (Purpose)
*      return Raleigh coefficient     
*
*     (Input)
*       BLAMDA : Wave Length of Cherencov Light ( nano Meter )
*
*     (Output)
*       WRAYSG  : Attenuation Coeff. ( 1/Meter )
*
*     (Creation Date and Author)
*     2006.4.18 updated K Bays 
*
************************************************************************

      IMPLICIT NONE

#include "skdetsim_parm.h"
#include "rayfac.h"
#include "skhead.h"

      real wraysg, blamda
C     parameter(RAYPOW=4.774)  ! from SK-1 data, number not implemented yet

      IF (WATER_TUNE .EQ. AP_TUNE) THEN
         write(6,*) 'WRAYSG is not used for AP_TUNE'
         stop
      ENDIF
      WRAYSG = rayfac/BLAMDA**4
 
cc         WRAYSG = rayfac/BLAMDA**RAYPOW
     

      RETURN
      END
