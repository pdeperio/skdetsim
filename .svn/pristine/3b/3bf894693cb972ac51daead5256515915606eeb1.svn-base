************************************************************************
*     -------------------------------------
      FUNCTION WSRTSG(BLAMDA)
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
*
************************************************************************

      IMPLICIT NONE

      real COEFF, BLAMDA, WSRTSG

C---- Unit of ALAMDA is Meter, so Convert it into Nano Meter

cc      BLAMDA=ALAMDA*1.0E+9
cc    print*,'blam',BLAMDA

C---- Get the Attenuation Coefficient in short wave length

C     ---- Wave Length => 515 NM

      IF(BLAMDA.GE.515.) THEN
        COEFF=0.

C     ---- Wave Length =< 515 NM

      ELSE
        COEFF = 2.332e8/BLAMDA**4

      ENDIF

C---- Get the Attenuation Coefficiency

      WSRTSG = COEFF

      RETURN
      END




