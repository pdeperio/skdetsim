************************************************************************
*     -------------------------------------
      FUNCTION WTGDEFF(BLAMDA)
*     -------------------------------------
*
*     (Purpose)
*       Get the Effect of Gd2(SO4)3 on Atteniation Coefficient
*       from EGADS Data for  337 nm <= BLAMDA <= 700 nm
*		For BLAMDA > 700, Effect set to zero.
*
*	 	 Only gets called in wabssg.F and wsymscsg.F	
*		 Usees WTRGAD for BLAMDA < 337 nm
*
*     (Input)
*       BLAMDA : Wave Length of Cherencov Light ( nano Meter )
*
*     (Output)
*       WTGDEFF  : Attenuation Coeff. Difference of Gd component ( 1/Meter )
*
*     (Creation Date and Author)
*       2013. 6.18 ; Coded By P.Weatherly
*
************************************************************************

      IMPLICIT NONE
#include "skhead.h"
#include "skdetsim_parm.h"
#include "./guwtr.h"

      real gdeff, mylam, blamda, m, y0, wtgdeff
      real x1,x2,y1,y2,yumatch,WTRGAD	
      real*4 dif(7)          !For Linear extrap y 
      Integer i, flag
      real*4 wavelen(7) !This will be the wl points for linear extrap 	             



       yumatch = 0.04801/0.02422 ! 0.0480068/0.0242243 Matches WTRGAD to UDEAL at 337
       flag = -9
           
      DATA dif/.0480068, .008, 0.00103009, 0.00183509, 0.00172512,
     & .000242403, 0.000103571 /
      DATA wavelen/337., 375., 405., 432., 452., 473., 550. /	

      i=0
      m = 0.
      y0 =0. 
      mylam = 0.
      gdeff = 0.
      wtgdeff = 0.
      
      if ( blamda .lt. 337. ) flag = 3
      if ( blamda .gt. 700. ) flag = -1
      if ( blamda.ge.337. .and. blamda.le.550. ) flag = 1
      if ( blamda.gt.550. .and. blamda.le.700. ) flag = 5

      if(BLAMDA .lt. 337.) gdeff = yumatch*WTRGAD(BLAMDA)
     

C  Need to fix this to loop over index on fit points.
      if ( flag .eq. 1 ) then
	
        DO i=1,6
	   if (BLAMDA.ge.wavelen(i) .and. BLAMDA.le.wavelen(i+1)) then
		x1=wavelen(i)
		x2=wavelen(i+1)
		y1=dif(i)
		y2=dif(i+1)
           endif
	enddo
        
        m = (y2-y1)/(x2-x1)
        y0 = y1 -m*x1
        gdeff = y0 + m*BLAMDA
        flag = 3
        
      endif 
      if (flag .eq. 5) then
          gdeff = EXP(-1.*(BLAMDA-550)/50.)
          gdeff = 0.000103571*gdeff
	  flag = 3
c          print *, "2718 Blambda, gdeff", blamda, gdeff

      endif
      if ( gdeff .lt. 0. ) gdeff = 0.  !the difference should be positive, 

      if ( flag .ne. 3 ) gdeff = 0. !only should pass if flag == 3
      wtgdeff = gdeff
      RETURN
      END
