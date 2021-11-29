C set absorption coefficient as well as Rayleigh and Mie
C scattering coefficents as a function of water transparency
C
      SUBROUTINE trabs_sk4(wt,absfitfac)

      IMPLICIT NONE

c table to modify absorption coefficient 
      INTEGER    ntab
      PARAMETER (ntab=31)
      REAL factor(ntab)/0.55,0.60,0.65,0.70,0.725,0.75,0.775,
     ,     0.80,0.825,0.85,0.875,0.90,0.925,0.95,0.975,1.00,
     ,     1.025,1.05,1.075,1.10,1.125,1.15,1.175,1.20,
     ,     1.25,1.30,1.35,1.40,1.45,1.50,1.55/
      REAL watert(ntab)/14815.497070,14586.312500,14258.296875,
     ,     13958.295898,13812.751953,13722.947266,13556.947266,
     ,     13420.690430,13250.173828,13159.061523,13034.945312,
     ,     12907.804688,12767.036133,12656.541992,12559.250977,
     ,     12459.154297,12303.238281,12250.059570,12089.152344,
     ,     11971.659180,11901.373047,11800.449219,11705.000977,
     ,     11622.160156,11429.535156,11219.772461,11027.351562,
     ,     10876.535156,10706.962891,10580.900391,10415.204102/


      REAL    wt,absfitfac

      INTEGER i
      REAL    delx,del1,del2

C do linear extrapolation of absorption coefficient if I must
      IF (wt .GT. watert(1)) THEN
         delx=watert(1)-watert(2)
         del1=watert(1)-wt
         del2=wt-watert(2)
         absfitfac=(factor(1)*del2+factor(2)*del1)/delx
         print *,'Water transparency out of range: wt=',wt,'abs=',absfitfac
         RETURN
      ENDIF
      IF (wt .LT. watert(ntab)) THEN
         delx=watert(ntab-1)-watert(ntab)
         del1=watert(ntab-1)-wt
         del2=wt-watert(ntab)
         absfitfac=(factor(ntab-1)*del2+factor(ntab)*del1)/delx
         print *,'Water transparency out of range: wt=',wt,'abs=',absfitfac
         RETURN
      ENDIF

C do linear interpolation of absorption coefficient if I can
      i=1
      DO WHILE((i .LT. ntab) .AND. (watert(i+1) .GT. wt))
         i=i+1
      ENDDO
      delx=watert(i)-watert(i+1)
      del1=watert(i)-wt
      del2=wt-watert(i+1)
      absfitfac=(factor(i)*del2+factor(i+1)*del1)/delx
      print *, '***********************************************'
     /     //'**********'
      print *,'trabs_sk4 : wt=',wt,' absfit factor=',absfitfac
      print *, '***********************************************'
     /     //'**********'
      RETURN
      END
