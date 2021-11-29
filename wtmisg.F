************************************************************************
*     -------------------------------------
      FUNCTION WTMISG(BLAMDA)
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
*       2001. 7.01 ; Created by Shibata 
*        ----  Difine Mie Scatteing -----
*               WTMISG = Total Mie Scattering 
*                        Total Mie = Symmmetric Mie + Asymmetric Mie  
*
************************************************************************

      IMPLICIT NONE

      real COEFF, BLAMDA, WTMISG

C---- Unit of ALAMDA is Meter, so Convert it into Nano Meter

cc      BLAMDA=ALAMDA*1.0E+9
cc    print*,'blam',BLAMDA

C---- Get the Attenuation Coefficient in short wave length

C     ---- Wave Length => 515 NM

      IF(BLAMDA.GE.515.) THEN
        COEFF=0.

C     ---- Wave Length =< 515 NM

      ELSE

         COEFF = 0.6015e23/BLAMDA**10       

      ENDIF

C---- Get the Attenuation Coefficiency

C     WMIESG = COEFF
      WTMISG = COEFF
    
      RETURN
      END
