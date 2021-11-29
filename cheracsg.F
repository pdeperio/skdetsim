************************************************************************
*     ---------------------
      FUNCTION CHERACSG(BETA)
*     ---------------------
*
*     (Purpose)
*       Return num of chernkov photon per 1cm for Acrylic medium
*       
*
*     (Input)
*       BETA (P/GETOT)
*
*     (Output)
*       num of chernkov photon per 1cm for acrylic
*
*     (Creation Date and Author)
*       2002. 8.30 ; First version by C.Sterner
*
************************************************************************
      implicit none
#include <geant321/gckine.inc>

      integer max
      parameter(max=37)
      real BETALST(max), CHERNUM(max)
      real BETA, cheracsg, RBETA
      integer isave, I
      save isave, BETALST, CHERNUM

      data CHERNUM/0.0,0.4657580,15.31776,36.63940,61.00795,84.35686,
     &     106.7122,128.1296,148.6610,168.3543,187.2547,205.4040,
     &     222.8418,239.6045,255.7266,271.2404,286.1761,300.5620,
     &     314.4248,327.7897,340.6799,353.1180,365.1247,376.7198,
     &     387.9217,398.7483,409.2159,419.3405,429.1366,438.6184,
     &     447.7992,456.6917,465.3077,473.6586,481.7551,489.6076,
     &     497.2257/

c data input
      if(isave.eq.0) then
         isave=1
         BETALST(1)=.64
         do i=2,max
            BETALST(i) = BETALST(i-1)+0.01
c            write(*,*)  BETALST(i), CHERNUM(i)
c            write(*,*) 'ckovacr.dat was read'
         enddo
      endif


      IF(BETA.LE.BETALST(1)) THEN
         CHERACSG = 0.

      ELSE IF(BETA.GT.BETALST(max)) THEN
         WRITE(6,*) 'Warning BETA GT 1. ; in CHERSG '
         CHERACSG = 0.
      ELSE
         DO 100 I=2,max
            IF(BETA.GT.BETALST(I)) GOTO 100
            RBETA = (BETA - BETALST(I-1)) * 100.
            CHERACSG = CHERNUM(I-1)*(1.-RBETA)+CHERNUM(I)*RBETA
            CHERACSG = CHERACSG*(CHARGE**2)
            GOTO 101
 100     CONTINUE
      ENDIF

 101  CONTINUE

      RETURN
      END
      
