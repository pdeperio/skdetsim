************************************************************************
*     -----------------------
      FUNCTION DREFSG(ALAMDA)
*     -----------------------
*
*     (Purpose)
*       Calculate differential index of refraction
*
*     (Input)
*       Wave length (nm)
*
*     (Output)
*       none
*
*     (Creation Date and Author)
*       1995.12.01 ; First version by K.Inoue
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
	real clamda
	
	clamda = alamda
	if(clamda.le.140.or.clamda.ge.1000.)then
	  write(*,*)'ERROR in drefsg'
	  stop
	end if
      WL=ALAMDA/1000.
      DREFSG= ((-2.*A1*WL/(WL**2-WLA)**2+2.*A3*WL+4.*A4*WL**3
     &          +6.*A5*WL**5)/2./
     &      SQRT(A1/(WL**2 - WLA)+A2+A3*WL**2+A4*WL**4+A5*WL**6)
     &     +(2.*B2*WL+4.*B3*WL**3)*(TEMP-TEMPB)
     &     +(2.*B5*WL+4.*B6*WL**3)*(TEMP-TEMPB)**2
     &     +(2.*B8*WL+4.*B9*WL**3)*(TEMP-TEMPB)**3
     &     +(2.*C2*WL+2.*C4*WL*TEMP)*(PRES-PRESB)
     &     +2.*C6*WL*(PRES-PRESB)**2)/1000.

      RETURN
      END
