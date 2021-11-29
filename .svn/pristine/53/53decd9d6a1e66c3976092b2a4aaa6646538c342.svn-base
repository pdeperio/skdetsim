***********************************************************************
*     -------------------------
      SUBROUTINE SGPIEL(P,COST)
*     -------------------------
*
*     (Purpose)
*       Calculate PION elastic scattering
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

      DIMENSION PIN(3),EV(3)
      DATA EV/1.,0.,0./
      DATA AM/139.5669/
      DATA PM/938.2796/

      E=SQRT(P**2+AM**2)
      BC=P/(E+PM)
      GC=1./SQRT(1.-BC**2)
      S=AM**2+PM**2+2.*PM*E
      PC=PM*P/SQRT(S)

 1    CONTINUE

      CALL RNEXP(0.0001054,RND)
      T = - RND
      IF(T.LT.-5.E+4) GO TO 1

      COST=T/2./PC**2+1.
      IF(COST.LT.-1.) GO TO 1
      SINT=SQRT(1.-COST**2)
C     
C---CONSIDER FERMI MORTION & BINDING ENERGY
C     
      CALL RNGAUS(PC,10.,RPC)
      PC=RPC
      IF(PC.LT.1.) PC=1.

      PIN(1)=PC*COST
      PIN(2)=PC*SINT
      PIN(3)=0.

      CALL VECBST(PIN,AM,EV,GC)

      PP=PIN(1)**2+PIN(2)**2+PIN(3)**2
      P=SQRT(PP)
      COST=PIN(1)/P

      RETURN
      END
