************************************************************************
*     -----------------------------
      FUNCTION TYVKSG(WAVELN)
*     -----------------------------
*
*     (Purpose)
*       Get the Reflectance of TYVEK
*
*     (Input)
*       Wave Length (meter) ; WAVELN
*
*     (Output)
*       Reflectance (%)     ; TYVKSG
*
*     (Creation Date and Author)
*       1995.04.26 ; First version by Y.Koshio
*
************************************************************************
C
      PARAMETER (NMAX=11)
      DIMENSION WL(NMAX),RF(NMAX)
      DATA WL/ 300. , 325. , 350. , 375. , 400. , 425. , 450. , 500.,
     +         550. , 575. , 600.  /
      DATA RF/ 0.691, 0.710, 0.839, 0.879, 0.888, 0.895, 0.899, 0.897,
     +         0.894, 0.893, 0.892 /

      BLAMDA = WAVELN * 1.0E+09
      IF(BLAMDA.LT.WL(1))THEN
         TYVKSG=RF(1) + (RF(2)-RF(1))/(WL(2)-WL(1))*(BLAMDA-WL(1))
         IF(TYVKSG.LE.0.)TYVKSG=0.

      ELSE IF(BLAMDA.GT.WL(NMAX)) THEN
         TYVKSG=RF(NMAX)

      ELSE
         DO 10 I=1,NMAX-1
            IF(BLAMDA.LE.WL(I+1)) THEN
               TYVKSG=RF(I)+
     &              (RF(I+1)-RF(I))/(WL(I+1)-WL(I))*(BLAMDA-WL(I))
               GOTO 999
            ENDIF
 10      CONTINUE

      ENDIF

 999  RETURN
      END
