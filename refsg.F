************************************************************************
*     -----------------------
      FUNCTION REFSG(ALAMDA)
*     -----------------------
*
*     (Purpose)
*       Calculate index of refraction
*
*     (Input)
*       Wave length (nm)
*
*     (Output)
*       index of refraction
*
*     (Creation Date and Author)
*       1995.12.01 ; First version by T.Yamaguchi
*            data taken from J.Phys.Chem.Ref.Data
*                    vol.14 No.4 1985
*
************************************************************************

      REAL WL,TEMP,PRES
      
      PARAMETER(TEMP = 13.0,PRES = 1.01325)
      PARAMETER(WLA=0.018085,TEMPB=19.993,PRESB=1.01325)

      PARAMETER(A1=5.743534E-3,A2=1.769238,A3=-2.797222E-2)
      PARAMETER(A4=8.715348E-3,A5=-1.413942E-3)
      PARAMETER(B1=-8.454823E-5,B2=-2.787742E-5,B3=2.608176E-6)
      PARAMETER(B4=-2.050671E-6,B5=1.019989E-6,B6=-2.611919E-6)
      PARAMETER(B7=8.194989E-9,B8=-8.10707E-9,B9=4.877274E-8)
      PARAMETER(C1=8.419632E-6,C2=1.941681E-5,C3=-7.762524E-8)
      PARAMETER(C4=4.371257E-8,C5=7.089664E-9,C6=-2.240384E-8)

      WL=ALAMDA/1000.
      REFSG = SQRT(A1/(WL**2 - WLA)+A2+A3*WL**2+A4*WL**4+A5*WL**6)
     &+(B1+B2*WL**2+B3*WL**4)*(TEMP-TEMPB)
     &+(B4+B5*WL**2+B6*WL**4)*(TEMP-TEMPB)**2
     &+(B7+B8*WL**2+B9*WL**4)*(TEMP-TEMPB)**3
     &+(C1+C2*WL**2+(C3+C4*WL**2)*TEMP)*(PRES-PRESB)
     &+(C5+C6*WL**2)*(PRES-PRESB)**2

      RETURN

      END
