************************************************************************
*     -----------------
*     INCLUDE "rayfac.h"
*     -----------------
*     COMMON VARIABLES for LOWE water scattering parameters
*     rayfac : symmetric scattering (Rayleigh + Symmetric Mie)
*     miefac : asymmetric scattering (Asymmetric Mie)
*     absfac : absorption for short wavelength
*     ablfac : absorption for long wavelength (for water scattering analysis)
*     LELONGFAC (LE) : ABSORPTION parameter for longer wavelength
*                      under BLAMDASP
*     BLAMDASP  (LE) : cross point of popefry and
*                      TUNED LONG wavelength fitting curve
*     raysfac : second symmetric scattering factor    | Added for new
*     miesfac : second asymmetric scattering factor   | functions tuned
*     abssfac : second absorption factor (lelongfac)  | by Thomas McLachlan
*     abspow  : absorption function power (lelongpow) | for SKIV SKDETSIM
*     absfit  : multiplicative factor to vary ABS fn. | 
*     asyfit  : multiplicative factor to vary ASY fn. | Aug, 2009
*
************************************************************************
#ifndef WTSCTT
      real*4 rayfac,miefac,absfac,LELONGFAC,LELONGPOW,BLAMDASP,
     &   raysfac,miesfac,abssfac,abspow,absfit,asyfit 
      common /rayfac/rayfac,miefac,absfac,LELONGFAC,BLAMDASP,
     &   raysfac,miesfac,abssfac,abspow,absfit,asyfit 
      PARAMETER(LELONGPOW = 10.6296)
#else
      real*4 rayfac,miefac,absfac,ablfac,LELONGFAC,LELONGPOW,BLAMDASP,
     &   raysfac,miesfac,abssfac,abspow,absfit,asyfit 
      common /rayfac/rayfac,miefac,absfac,ablfac,LELONGFAC,BLAMDASP,
     &   raysfac,miesfac,abssfac,abspow,absfit,asyfit 
      PARAMETER(LELONGPOW = 10.6296)
#endif
