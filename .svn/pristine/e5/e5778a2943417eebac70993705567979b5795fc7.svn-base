************************************************************************
*     -----------------------
      FUNCTION EFFNSG(ALAMDA)
*     -----------------------
*
*     (Purpose)
*       Calculate effective index of refraction
*         (light velocity shoud be group velocity)
*             Neff(l)=N(l)-l*N'(l)
*
*     (Input)
*       Wave length (nm)
*
*     (Output)
*       effective index of refraction
*             Vg=C/Neff
*
*     (Creation Date and Author)
*       1995.12.01 ; First version by K.Inoue
*
************************************************************************

      EFFNSG=REFSG(ALAMDA)-ALAMDA*DREFSG(ALAMDA)
      RETURN
      END
