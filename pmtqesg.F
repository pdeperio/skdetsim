C   29/10/92 210291112  MEMBER NAME  PMTQE    *.FORT     M  E2FORT
      REAL FUNCTION PMTQESG(WAVELN)
*
*     Quantum efficiency of a 20-inch PMT (R3600)
*
*     Input : WAVEL ; Wavelength of a photon (nm)
*     Return value  ; Quantum efficiency (fraction of 1)
*
*     Author: M. Mori   29-OCT-1992
*
*     mesured value was changed by M.Shiozawa 15-FEB-1997
*

      implicit none

#include "skhead.h"
#include "skdetsim_parm.h"

      real WAVELN

      integer nmax,i
      PARAMETER (NMAX=20)

      real WL(NMAX),QE_LE_SK1(NMAX),QE(NMAX)
      DATA WL/  280.,  300.,  320.,  340.,  360.,  380.,  400.,  420.,
     +          440.,  460.,  480.,  500.,  520.,  540.,  560.,  580.,
     +          600.,  620.,  640.,  660./
      DATA QE   /  0.0,  1.39, 8.54, 16.9, 20.3, 20.6, 21.1, 20.2,
     + 	        18.8, 16.7, 14.0, 11.6, 8.06,  4.32, 2.646, 1.462,
     +         0.756, 0.508, 0.158, 0.000/
      DATA QE_LE_SK1/ 0.622,  6.80, 15.30, 20.72, 21.80, 20.82, 20.24,
     +     18.82, 17.38, 15.30, 12.70, 11.28, 8.06, 4.32, 2.646, 1.462,
     +     0.756, 0.508, 0.158, 0.000/

      if (REFLECTION_TUNE.eq.LE_TUNE .AND. SK_GEOMETRY.eq.SK_I) THEN
         do i=1,nmax
            QE(i)=QE_LE_SK1(i)
         enddo
      endif

      PMTQESG = 0.
      IF(WAVELN.LT.WL(1) .OR. WAVELN.GT.WL(NMAX)) RETURN

      DO 10 I=1,NMAX-1
        IF(WAVELN.LE.WL(I+1)) THEN
          PMTQESG = 
     &          QE(I) + (QE(I+1)-QE(I))/(WL(I+1)-WL(I))*(WAVELN-WL(I))
          PMTQESG = PMTQESG*0.01
          RETURN
        ENDIF
   10 CONTINUE

      RETURN
      END
