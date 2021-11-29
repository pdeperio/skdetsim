***********************************************************************
*     ----------------------------------
      SUBROUTINE RNGAUSYOO(XMEAN,SIGMA,RND)
*     ----------------------------------
*     
*     ( purpose )
*       generate gaussian distributed randum number
*       For Water Scattering laser beam width
*
*     ( input ) 
*       XMEAN : MEAN OF DISTRIBUTION
*       SIGMA : SIGMA OF DISTRIBUTION
*
*     ( output ) 
*       RND   : RANDUM ANGLE 
*
*     ( creation date and author )
*       ????.??.?? First version by J.H. Yoo?
*
***********************************************************************


      real a,rdummy,sigma,xmean,rnd,ang
      integer I

      A=0.0

      DO 10 I=1,12
         Call GRNDM(RDUMMY,1)
         A=A+RDUMMY
 10   CONTINUE
C      write(*,*)'A    : ',A      
      RND = (A-6.)*SIGMA+XMEAN
      ang = sqrt(RND**2/(3610.**2))
      ang = sqrt(ang)
      RND = atan(ang)
C      write(*,*)'ang  ;',ang
C      write(*,*)'ANGLE:',RND
      RETURN
      END

