************************************************************************
*     -------------------------------------
      SUBROUTINE RNDARK(PEDARK)
*     -------------------------------------
*
*     (Purpose)
*        return dark noise p.e randomly 
*             pe< 0     : exp(0.6*pe)
*        0   <pe< 5.57  : use table
*        5.57<pe<30     : exp(-2.78*pe)
*        30  <pe        : exp(-0.067*pe)
*     (Input)
*        None
*
*     (Output)
*        PEDARK  ; random number in dark noise p.e. distribution 
*     
*     (Creation Date and Author)
*       1997.01.03 ; First version by T.Hayakawa
*
************************************************************************
      implicit none
c arg
      real    pedark ! output

c const 
      real      P_NEGATIVE
      parameter(P_NEGATIVE=0.00781042)
      real      P_GT5                , P_GT30
      parameter(P_GT5     =0.96      , P_GT30=0.99974772)
c     parameter(P_GT5     =0.95421719, P_GT30=0.99974772)

      real      SLOPE0        , SLOPE1        , SLOPE2
      parameter(SLOPE0=  0.6  , SLOPE1= -0.278, SLOPE2= -0.067)
c      parameter(SLOPE0=  0.511, SLOPE1= -0.278, SLOPE2= -0.067)

c local variable
      real    pe(100),pe1(50),pe2(50)

      EQUIVALENCE (pe( 1),pe1(1))
      EQUIVALENCE (pe(51),pe2(1))
       
      real    x,dummy
      real*8  xx
      integer i
      integer icall/0/
      real    const0,const1,const2

c function
      real    RNDM

      DATA pe1/0.017699, 0.098531, 0.125937, 0.152355, 0.178772, 
     &         0.203579, 0.221802, 0.240025, 0.258248, 0.276470, 
     &         0.294693, 0.313075, 0.331522, 0.349969, 0.368417, 
     &         0.386864, 0.405710, 0.425543, 0.445376, 0.465210, 
     &         0.485043, 0.504981, 0.525241, 0.545502, 0.565762, 
     &         0.586023, 0.606093, 0.625742, 0.645390, 0.665038, 
     &         0.684686, 0.704181, 0.723133, 0.742084, 0.761036, 
     &         0.779988, 0.798939, 0.817475, 0.835985, 0.854496, 
     &         0.873007, 0.891517, 0.909880, 0.928118, 0.946356, 
     &         0.964594, 0.982832, 1.001103, 1.019914, 1.038725/
      DATA pe2/1.057537, 1.076348, 1.095159, 1.114437, 1.133878, 
     &         1.153319, 1.172760, 1.192201, 1.212579, 1.233584, 
     &         1.254589, 1.275595, 1.296600, 1.319222, 1.342156, 
     &         1.365090, 1.388024, 1.412304, 1.438058, 1.463811, 
     &         1.489564, 1.517882, 1.547948, 1.578013, 1.609295, 
     &         1.643891, 1.678486, 1.715516, 1.756550, 1.797584, 
     &         1.846027, 1.894933, 1.952585, 2.013362, 2.083043, 
     &         2.165416, 2.261218, 2.375183, 2.514669, 2.695738, 
     &         2.925270, 3.220021, 3.601475, 4.084537, 4.696189, 
     &         5.473693, 6.532504, 7.998651, 100.    , 100.     /     

c program
      if(icall.eq.0)then
         icall  = 1
         const1 = -SLOPE1/(exp(pe(96)*SLOPE1)-exp(30.*SLOPE1))
         const2 = -SLOPE2/(exp(  30. *SLOPE2)- 0.            )
         const0 =  SLOPE0
      endif

C      X = RNDM(DUMMY)
      CALL GRNDM(X,1)

      if(x.lt.P_NEGATIVE)then
 3       continue
C        X = RNDM(DUMMY)
         CALL GRNDM(X,1)
         if(x.le.0.)goto 3
         pedark = log(x)/SLOPE0
         if(pedark.le.-20.)pedark=-20
      elseif(x.ge.P_GT5.and.x.le.P_GT30)then
 1       continue
C            X = RNDM(DUMMY)
            CALL GRNDM(X,1)
            x = x*SLOPE1/const1+exp(pe(96)*SLOPE1)
         if(x.le.0.)goto 1
         pedark = log(x)/SLOPE1
      elseif(x.gt.P_GT30)then
 2       continue
C            X = RNDM(DUMMY)
            CALL GRNDM(X,1)
            x = x*SLOPE2/const2+exp(30.*SLOPE2)
         if(x.le.0.)goto 2
         pedark = log(x)/SLOPE2
         if(pedark.ge.250.)pedark=250.
      else
         i =  int(x*100.)
         if(i.eq.0)then
            pedark=(x-P_NEGATIVE)*pe(1)
         else
C            x=rndm(dummy)
            CALL GRNDM(X,1)
            pedark= pe(i)+ x*(pe(i+1)-pe(i))
c            pedark= pe(i)+ (x-float(i)*0.01)*(pe(i+1)-pe(i))
         endif
      endif

      RETURN
      END
