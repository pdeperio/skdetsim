C set absorption coefficient as well as Rayleigh and Mie
C scattering coefficents as a function of water transparency
C
      SUBROUTINE trabs(wt,abs,ray,mie)
         IMPLICIT NONE

c standard parameters at 83 m water transparency
         REAL       abs83,ray83,mie83
         PARAMETER (abs83=0.61E8,ray83=0.8E8,mie83=1E8)

c table to modify absorption coefficient 
         INTEGER    ntab
C for enelf version 7
C         PARAMETER (ntab=15)
C         REAL factor(ntab)/ 0.20, 0.30,0.40,0.50,0.60,0.70,0.80,
C     ,            0.90,1.00,1.05, 1.20,1.35,1.50,1.80,2.10/
C         REAL watert(ntab)/10412,10099,9827,9567,9317,9053,8815,
C     ,            8566,8300,8178, 7863,7567,7292,6778,6246/
         PARAMETER (ntab=24)
         REAL factor(ntab)/ 0.20, 0.25, 0.30, 0.35,0.40,0.45,
     ,       0.50,0.55,0.60,0.65, 0.70, 0.75, 0.80,0.90,1.00,
     ,       1.10,1.20,1.30,1.40, 1.55, 1.70, 1.85,2.00,2.20/
         REAL watert(ntab)/10430,10280,10160,10080,9918,9796,
     ,       9660,9560,9406,9301, 9145, 9021, 8893,8631,8300,
     ,       8074,7888,7710,7532, 7241, 6997, 6752,6462,6110/

         REAL    wt,abs,ray,mie

         INTEGER i
         REAL    delx,del1,del2

C set scattering coefficients constant
         ray=ray83
         mie=mie83

C do linear extrapolation of absorption coefficient if I must
         IF (wt .GT. watert(1)) THEN
            delx=watert(1)-watert(2)
            del1=watert(1)-wt
            del2=wt-watert(2)
            abs=abs83*(factor(1)*del2+factor(2)*del1)/delx
            print *,'Water transparency out of range: wt=',wt,'abs=',abs
            RETURN
         ENDIF
         IF (wt .LT. watert(ntab)) THEN
            delx=watert(ntab-1)-watert(ntab)
            del1=watert(ntab-1)-wt
            del2=wt-watert(ntab)
            abs=abs83*(factor(ntab-1)*del2+factor(ntab)*del1)/delx
            print *,'Water transparency out of range: wt=',wt,'abs=',abs
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
         abs=abs83*(factor(i)*del2+factor(i+1)*del1)/delx
         RETURN
      END
