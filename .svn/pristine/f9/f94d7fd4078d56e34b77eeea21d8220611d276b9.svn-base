C routine to generate double gaussian random distribution      
C Version 1.0  M. Smy (April 4)
C
      SUBROUTINE rntwogauss(x,sigma,rnd)
         IMPLICIT NONE

         REAL      sigshift,sigsig,as,fac1,fac2
         PARAMETER (sigshift=0.56)
         PARAMETER (as=-0.26)
c     <<<<<<< gauss.F
C         PARAMETER (sigsig=sqrt(1-sigshift*sigshift))
C         PARAMETER (fac1=sqrt(1+as),fac2=sqrt(1-as))
c     =======
c         PARAMETER (sigsig=sqrt(1-sigshift*sigshift))
c         PARAMETER (fac1=sqrt(1+as),fac2=sqrt(1-as))
c     >>>>>>> 1.3

         REAL     x,sigma,rnd

         REAL     peak,sig,x0,rdummy

         sigsig=sqrt(1-sigshift*sigshift)
         fac1=sqrt(1+as)
         fac2=sqrt(1-as)

c         peak=rndm()
         CALL GRNDM(RDUMMY,1)
         peak=rdummy
         IF (peak .GT. 0.5) THEN
            sig=fac1*sigsig*sigma
            x0=x-sigshift*sigma
         ELSE
            sig=fac2*sigsig*sigma
            x0=x+sigshift*sigma
         ENDIF
         CALL rngausk2(x0,sig,rnd)
         RETURN
      END

C routine to generate asymmetric gaussian random distribution
C test version 0.0 Y.HENG (Nov 18, 2009)
C working version 1.0 Y.HENG/Y.OBAYASHI (Dec. 17, 2009)
      SUBROUTINE rntwogauss_asy(x, sigmaE, sigmaL, rnd)
      IMPLICIT NONE
      
      REAL x, sigmaE, sigmaL, rnd
      REAL rdummy, facE, peak

      facE = sigmaE / (sigmaE + sigmaL)
      CALL GRNDM(RDUMMY,1)
      peak = RDUMMY

      IF(peak .lt. facE) THEN
         CALL rngausk2(x, sigmaE, rnd)
         IF(rnd .gt. x) rnd = 2.0*x - rnd
      ELSE
         CALL rngausk2(x, sigmaL, rnd)
         IF(rnd .lt. x) rnd = 2.0*x - rnd
      ENDIF
      
      RETURN
      END

************************************************************************
*     ---------------------------
      subroutine rngausk2(x,sigma,rnd)
*     ---------------------------
*
*     (Purpose)
*       MAKE GAUSSIAN DISTRIBUTED RANDUM NUMBER
*
*     (Input)
*       X : MEAN
*       SIGMA : SIGMA
*
*     (Output)
*       RND : GAUSSIAN
*
*     (Creation Date and Author)
*     2003. 4. 02 : C. Sterner modify to remove storage of second gaussian
*                      random number.  This is done to ensure every call
*                      is indenpendent of past history.  Modification
*                      needed to use random seeds to restart events.
************************************************************************

 10   CALL GRNDM(RDUMMY,1)
      Z1=2.*RDUMMY-1.
      CALL GRNDM(RDUMMY,1)
      Z2=2.*RDUMMY-1.
      
      S=Z1*Z1+Z2*Z2
      IF(S.GT.1.)GOTO 10
      GRAN1=Z1*SQRT(-2.*LOG(S)/S)
      GRAN=GRAN1
c     
      RND=X+SIGMA*GRAN
c
      RETURN
      END


