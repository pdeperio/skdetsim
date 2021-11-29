************************************************************************
*     -------------------------------------
      FUNCTION WTRGAD(BLAMDA)
*     -------------------------------------
*
*     (Purpose)
*       Get the Atteniation Coefficiency in Gd component of EGADS Gd Solution
*
*     (Input)
*       BLAMDA : Wave Length of Cherencov Light ( nano Meter )
*
*     (Output)
*       WTRGAD  : Attenuation Coeff. of Gd component ( 1/Meter )
*
*     (Creation Date and Author)
*       2013. 6.11 ; Coded By T.Yano
*
************************************************************************

      IMPLICIT NONE
#include "skhead.h"
#include "skdetsim_parm.h"
#include "./guwtr.h"

      REAL WTRGAD,DENSITY,BASELINE,BCKGND,GADATTEN, BLAMDA
      DENSITY=0.188617 !   err. 0.000517603
      BASELINE =-0.000492909 !   err. 4.70343e-05
      BCKGND = 2.26306e+12 !   err. 4.70343e-05

      IF (BLAMDA .ge. 337) THEN
          WTRGAD=0
          RETURN
      ELSE IF (BLAMDA .lt. 337) THEN

          GADATTEN =
     $   0.0118147 *exp(-1*(BLAMDA- 244.642 )**2/(2* 0.569769 **2))+
     $   0.02514   *exp(-1*(BLAMDA- 246.342 )**2/(2* 0.365549 **2))+
     $   0.0240045 *exp(-1*(BLAMDA- 252.607 )**2/(2* 0.478219 **2))+
     $   0.0231366 *exp(-1*(BLAMDA- 305.804 )**2/(2* 0.410027 **2))+
     $   0.0430073 *exp(-1*(BLAMDA- 311.462 )**2/(2* 0.447006 **2))+
     $   0.608729  *exp(-1*(BLAMDA- 273.307 )**2/(2* 0.3747   **2))+
     $   0.266735  *exp(-1*(BLAMDA- 274.373 )**2/(2* 0.364533 **2))+
     $   0.367687  *exp(-1*(BLAMDA- 276.099 )**2/(2* 0.436379 **2))+
     $   0.0634436 *exp(-1*(BLAMDA- 279.065 )**2/(2* 0.391507 **2))

          WTRGAD = log(10.)
     $     *(10*(DENSITY*GADATTEN+BASELINE+BCKGND/((BLAMDA)**6)))
          RETURN
      END IF
      RETURN
      END
