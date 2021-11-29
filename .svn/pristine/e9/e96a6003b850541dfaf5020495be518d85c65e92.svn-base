************************************************************************
*     -----------------
      SUBROUTINE GUTRAK
*     -----------------
*
*     (Purpose)
*
*     (Input)
*       none
*
*     (Output)
*       none
*
*     (Creation Date and Author)
*       1992.12.03 ; First version by Y.Koshio
*
************************************************************************

c cregis include klksregeneration 
      include "klksregen.inc" 
#include "geant321/gconsp.inc"       
 
      DIMENSION RNDM(5) 
       
      CALL GRNDM(RNDM,1) 
      KLKSLIF = -CLIGHT*KLKSLIFE_H20*LOG(RNDM(1)) 
c end regeneration

      CALL GTRACK

*** Added 7/00 Alec Habig
*** If the particle is a penetrating particle - mu or charged pi
*** then store information about where and why it finished tracking
      CALL CHECK_FOR_STOPS

      RETURN
      END
