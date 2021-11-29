***********************************************************************
*     -------------------------
      SUBROUTINE SGPIIE(P,COST)
*     -------------------------
*
*     (Purpose)
*       Calculate PION inelastic scattering
*
*     (Input)
*       P      ; Momemtum of the particle (MeV)
*
*     (Output)
*       P      ; Momemtum of the particle after interaction (MeV)
*       COST   ; 
*
*     (Creation Date and Author)
*       1980.10.29 ; First version   by  K.Arisaka
*       1994.11.28 ; Modified for SK by  Koshio,Y
*
************************************************************************
#include "skcnst.h"

      DIMENSION PIN(3),EV(3),DIR(3)
      DATA EV/1.,0.,0./
      DATA DIR/1.,0.,0./
      DATA AM/139.5669/
      DATA PM/938.2796/

      E=SQRT(P**2+AM**2)
      BC=P/(E+PM)
      GC=1./SQRT(1.-BC**2)
      S=AM**2+PM**2+2.*PM*E
      PC=PM*P/SQRT(S)
C     
C---SET NEW MOMENTUM
C     
 5    CALL RNGAUS(0.,50.,RPC)
      PCR=PC-ABS(RPC)
      IF(PCR.LT.0.) GO TO 5
      PC=PCR
C     
C---SET DIRECTION
C     
      SIG=0.87
      IF(PC.GT.250.) SIG=0.44
      CALL RNSIG(SIG,DTH)

      CONE = DTH*180./PI
      CALL RNCONE(DIR,CONE,PIN)

      DO 10 I=1,3
         PIN(I)=PIN(I)*PC
 10   CONTINUE
C     
C---DO BOOST
C     
      CALL VECBST(PIN,AM,EV,GC)
C     
      PP=PIN(1)**2+PIN(2)**2+PIN(3)**2
      P=SQRT(PP)
      COST=PIN(1)/P
C     
      RETURN
      END
