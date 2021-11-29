************************************************************************
*     -------------------------------------
      FUNCTION WSYMSCSG(BLAMDA)
*     -------------------------------------
*
*     (Purpose)
*       Get the atteniation coefficiencts of symmetric scattering in water
*       Forward and backward scattered angle distribution is the same.
*       Related angle distribution is 1+cos(theta)^2.
*       (theta:angle between incident and scatterd direction)
*       This is called at only SK-III mode.
*
*
*     (Input)
*       BLAMDA : Wave Length of Cherencov Light ( nano Meter )
*
*     (Output)
*       WSYMSCSG  : Attenuation Coeff. of symmetric scattering ( 1/Meter )
*                    Symmetric scattering = Rayleigh + Symmetric Mie
*
*     (Creation Date and Author)
*
*       Apr. 5. 2007 ; First version by B. Yang
*                      SK-III ATMPD tunning result included
*       Apr. 8. 2007 ; SK-III LOWE tunning result included by B. Yang
*                      Meaningless Safe guard commented out
*       Sep.28. 2007 ; SK-III water coefficient functions retuned
*		       due to bug about black sheet by B. Yang
*       Mar.14. 2008 ; SK-III water parameter retuned with Feb, 08 data by B. Yang
*       Sep. 4. 2009 ; SK-IV monthly water parameters installed by T.McLachlan
*                      Symmetric function now tuned by 2 parameters
*
*       Dec.    2012 ; Change by T.Akiri
*       Jun.18. 2013 ; SK-IV Gd Effect added for WATER_TUNE == 10
*                      Same as UNI_TUNE otherwise, by P.Weatherly
************************************************************************

      IMPLICIT NONE
#include "skhead.h"
#include "skdetsim_parm.h"
#include "./guwtr.h"
#include "./rayfac.h"
#include "dsparm.h"

      REAL BLAMDA,PRELAMDA,PREVALUE
C      INTEGER*4 PRETUNE
      REAL WSYMSCSG, WTGDEFF

C      DATA PRELAMDA,PREVALUE,PRETUNE/0.,0.,0/
C      SAVE PRELAMDA,PREVALUE,PRETUNE
      DATA PRELAMDA,PREVALUE/0.,0./
      SAVE PRELAMDA,PREVALUE

C      IF (BLAMDA .EQ. PRELAMDA .AND. WATER_TUNE .EQ. PRETUNE) THEN
      IF (BLAMDA .EQ. PRELAMDA ) THEN
        WSYMSCSG = PREVALUE
        RETURN
      END IF
      
      IF(WATER_TUNE.EQ.AP_TUNE) THEN
C        WSYMSCSG = ARAS*(1.177E8/BLAMDA**4 + 66.96E16/BLAMDA**8 )	! Apr 07 ( data taken Feb 07 )
        WSYMSCSG = ARAS* (0.7761E8/BLAMDA**4 ) * (1.0 + 1.363E5/BLAMDA**2 )	! Sep 07 ( data taken Feb 07 )
C      ELSE IF(WATER_TUNE.EQ.LE_TUNE) THEN
      ELSE
C        WSYMSCSG = (rayfac/BLAMDA**4) * ( 1.0 + 5.890E9/BLAMDA**4 )	! Apr 07 ( data taken Feb 07 )
C        WSYMSCSG = (rayfac/BLAMDA**4) * ( 1.0 + 1.058E5/BLAMDA**2 )	! Sep 07 ( data taken Feb 07 )


       if (SK_GEOMETRY .eq. SK_IV) then

C       SKIV Tuned by month. See tr2fac.F for tuned variables. Thomas McLachlan.

        WSYMSCSG = (rayfac/BLAMDA**4) * (  1.0 + raysfac/BLAMDA**2 )

        if(IWATERP.eq.201208) then
* Tarek Akiri, Dec. 2012
* New function with parameter from 16/04/2009
* Use old fcn below 330nm because of behaviour of new fcn
          if(BLAMDA .GT. 330) then
             WSYMSCSG = 6.00317e+10/(BLAMDA**5)*EXP(-1.44e29/(BLAMDA**12))
          else
             WSYMSCSG = (rayfac/BLAMDA**4) * (  1.0 + raysfac/BLAMDA**2 )
          endif
        endif
       
       else if (SK_GEOMETRY .eq. SK_V) then  ! just copy from SK-IV

        WSYMSCSG = (rayfac/BLAMDA**4) * (  1.0 + raysfac/BLAMDA**2 )

       else
        WSYMSCSG = (rayfac/BLAMDA**4) * (  1.142925 + 8.266263E4/BLAMDA**2 )	! Mar 08 ( data taken Feb 08 )
       endif

C      ELSE
C        WRITE(0,*) 'ERROR : WATER TUNE is not usual,', WATER_TUNE
C        stop
      ENDIF

C by P.Weatherly for Gd-loaded water
      if(WATER_TUNE.eq.10) then
         WSYMSCSG = WSYMSCSG + WTGDEFF(BLAMDA)
      endif		

      PRELAMDA = BLAMDA
      PREVALUE = WSYMSCSG
C      PRETUNE = WATER_TUNE

      RETURN
      END
