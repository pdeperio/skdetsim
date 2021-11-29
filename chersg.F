************************************************************************
*     ---------------------
      FUNCTION CHERSG(BETA)
*     ---------------------
*
*     (Purpose)
*       Return num of chernkov photon per 1cm
*       This table is made by /skam/pro/example/cherenkov/mktabch.F
*       (instead of GNCKOV)
*
*     (Input)
*       BETA (P/GETOT)
*
*     (Output)
*       num of chernkov photon per 1cm
*
*     (Creation Date and Author)
*       1996. 2.28 ; First version by Y.Koshio
*
************************************************************************

#include "geant321/gckine.inc"

      DIMENSION BETALST(28),CHERNUM(28)

      DATA BETALST /0.73, 0.74, 0.75, 0.76, 0.77, 0.78, 0.79,
     &     0.80, 0.81, 0.82, 0.83, 0.84, 0.85, 0.86, 0.87,
     &     0.88, 0.89, 0.90, 0.91, 0.92, 0.93, 0.94, 0.95,
     &     0.96, 0.97, 0.98, 0.99, 1.00 /

      DATA CHERNUM /0., 1.07477, 13.2025, 35.5941, 57.2964, 78.1693,
     &     98.2547, 117.592, 136.217, 154.164, 171.467, 188.156,
     &     204.260, 219.804, 234.816, 249.319, 263.336, 276.888,
     &     289.996, 302.679, 314.955, 326.841, 338.354, 349.509,
     &     360.321, 370.803, 380.969, 390.832 /

      IF(BETA.LE.BETALST(1)) THEN
         CHERSG = 0.

      ELSE IF(BETA.GT.BETALST(28)) THEN
         WRITE(6,*) 'Warning BETA GT 1. ; in CHERSG '
         CHERSG = 0.

      ELSE
         DO 100 I=2,28
            IF(BETA.GT.BETALST(I)) GOTO 100
            RBETA = (BETA - BETALST(I-1)) * 100.
            CHERSG = CHERNUM(I-1)*(1.-RBETA)+CHERNUM(I)*RBETA
            CHERSG = CHERSG*(CHARGE**2)
            GOTO 101
 100     CONTINUE
      ENDIF

 101  CONTINUE

      RETURN
      END
