******************************************************
*     -----------------
      SUBROUTINE SGGVUS
*     -----------------
*
*     (Purpose)
*        User routine for Generating Vector
*
******************************************************

      implicit none

#ifdef NICKEL

#include "nicpara.h"

      integer i,nvct
      real posx,posy,posz
      REAL POS(3,4),PMOM(3,4),TIM(4),EGAM(4),P(3), rn(4)
      integer IP(4),MODE,NGAM

      CALL GAMNI(MODE,NGAM,EGAM)
c      write(6,*)'EGAM=',(egam(i),i=1,ngam)
c      WRITE(6,*) 'NGAM,EGAM=',Ngam,Egam(1)

 11   continue
      call grndm(rn,4) 
      if(9.*sqrt(rn(1)).lt.2.5.and.17.*rn(3).gt.5.) goto 11

      posx= sorx+9.*sqrt(rn(1))*cos(2.*3.1415927*rn(4))
      posy= sory+9.*sqrt(rn(1))*sin(2.*3.1415927*rn(4))
      posz=  sorz+17.*rn(3)-8.5 

      DO I=1,NGAM /* If you use gamni, please use this line */

         POS(1,I) = posx
         POS(2,I) = posy
         POS(3,I) = posz
         TIM(I) = 0.
         IP(I) = 22
         CALL RNDIR(P)
         PMOM(1,I) = P(1) * EGAM(I)
         PMOM(2,I) = P(2) * EGAM(I)
         PMOM(3,I) = P(3) * EGAM(I)

      ENDDO

      NVCT = NGAM  /* If you use gamni, please use this line */
      CALL VCMPMK(NVCT,POS,PMOM,TIM,IP)

#else
#ifdef DTG

      INTEGER IPART,IDP(2),NVCT
      REAL POS(3,4),P(3),PMOM(3,4),TIM(4),PP

#include "dtginfo.h"

      real m_elec,m_gamm
      integer bob_flag
      real rands(5),testr,rlengthx,rlengthy,rlengthz
      real posx, posy, posz

C      print*,'SGGVUS!!!' ! S.Ito

 101  continue


      print*,'ENTERING DTG SGGVUS!!!'

      read(4,*,END=232) bob_flag,m_elec,m_gamm
      write(*,*) 'bob_flag,m_elec,m_gamm: ',bob_flag,m_elec,m_gamm

c      m_elec = 20000.0
c      m_gamm = 0.0

 110  continue
      
C      print*,'SGGVUS X' ! S.Ito

      call ranmar(rands,5)
      rlengthx = rands(1)*100.0

****  Distribute x over 24 cm sigma
      testr = 0.01663* exp(-1.0*rlengthx**2*0.000868)
      if(rands(4).gt.testr) then
         goto 110
      endif

 120  continue

C      print*,'SGGVUS Y' ! S.Ito

      call ranmar(rands,5)
      rlengthy = rands(1)*100.0

****  Distribute y over 24 cm sigma
      testr = 0.01663* exp(-1.0*rlengthy**2*0.000868)
      if(rands(4).gt.testr) then
         goto 120
      endif

 130  continue

C      print*,'SGGVUS Z' ! S.Ito

      call ranmar(rands,5)
      rlengthz = rands(1)*100.0

****  Distribute z over 40 cm sigma
      testr = 0.009976* exp(-1.0*rlengthz**2*0.000313)
      if(rands(4).gt.testr) then
         goto 130
      endif
      
      call ranmar(rands,5)
      if (rands(1).gt.0.5) then
         posx = dtg_x(dtg_runnum) + (rlengthx/2.0)
      else
         posx = dtg_x(dtg_runnum) - (rlengthx/2.0)
      endif
      if (rands(2).gt.0.5) then
         posy = dtg_y(dtg_runnum) + (rlengthy/2.0)
      else
         posy = dtg_y(dtg_runnum) - (rlengthy/2.0)
      endif
      if (rands(3).gt.0.5) then
         posz = dtg_z(dtg_runnum) + (rlengthz/2.0)
      else
         posz = dtg_z(dtg_runnum) - (rlengthz/2.0)
      endif

c      posx = 0.0
c      posy = 0.0
c      posz = 0.0
      POS(1,1) = posx
      POS(2,1) = posy
      POS(3,1) = posz
      POS(1,2) = posx
      POS(2,2) = posy
      POS(3,2) = posz

      PP = 8.

      TIM(1) = 0.
      TIM(2) = 0.

      CALL RNDIR(P)

c      P(1) = 0.0
c      P(2) = 0.0
c      P(3) = -1.0

      PMOM(1,1) = P(1) * m_elec/1000.0
      PMOM(2,1) = P(2) * m_elec/1000.0
      PMOM(3,1) = P(3) * m_elec/1000.0

      CALL RNDIR(P)
      PMOM(1,2) = P(1) * m_gamm/1000.0
      PMOM(2,2) = P(2) * m_gamm/1000.0
      PMOM(3,2) = P(3) * m_gamm/1000.0
      IDP(1) = 11 ! electron
      IDP(2) = 22 ! gamma
      write(*,*) 'sgguvs dtg    (x,y,z)',pos(1,1),pos(2,1),pos(3,1)
      write(*,*) '  electron (px,py,pz)',pmom(1,1),pmom(2,1),pmom(3,1)
      write(*,*) '  gamma    (px,py,pz)',pmom(1,2),pmom(2,2),pmom(3,2)
      write(*,*) ' '

c      IPART = IDP(1)
c      CALL VCMAKE(POS,P,IPART)

      NVCT = 2 
      CALL APVCMPMK(NVCT,POS,PMOM,TIM,IDP)

      RETURN

 232  Continue
      write (*,*) 'Input E file is rewound'
      rewind 4
      goto 101

#else

      INTEGER IPART
      REAL POS(3),PMOM(3),PP

#ifdef LINAC
#include "linpara.h"
#include "lininfo.h"
      real rn(2),xtrg,ytrg

  100 call grndm(rn,2)
      
      xtrg = 1.2*sqrt(rn(1))*cos(2.*3.1415927*rn(2))
      ytrg = 1.2*sqrt(rn(1))*sin(2.*3.1415927*rn(2))
      if (ytrg.ge.0.5) go to 100

      pos(1)=sorx+xtrg
      pos(2)=sory+ytrg
      pos(3)=sorz+10.0
c      write(*,*) 'sgguvs',pos(1),pos(2),pos(3)

      PP = pinit
      PMOM(1)=0.
      PMOM(2)=0.
      PMOM(3)=-pp
      IPART =11

#else
      CALL RNPOS(POS,150.)
      CALL RNDIR(PMOM)
      PP = 8.
      PMOM(1)=PMOM(1) * PP
      PMOM(2)=PMOM(2) * PP
      PMOM(3)=PMOM(3) * PP
      IPART =22

#endif

      CALL VCMAKE(POS,PMOM,IPART)

#endif
#endif
      RETURN
      END

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      SUBROUTINE GAMNI(MODE,NGAM,EGAM)
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
C     
C     ++ FUNCTION TO RETURN RANDOM GAMMA ENERGY FROM NI(N,GAMMA)NI
C     M.NAKAHATA 07-NOV-86
C     ++ modified by '77 NI(N,GAMMA)NI data
C     T.Maruyama 12-Oct-'95
C--INPUT
C     NOTHING
C--OUTPUT
C     MODE    : DECAY MODE
C     NGAM    : NUMBER OF GAMMA
C     EGAM(I) : ENERGY OF "I"TH GAMMA  (MEV)
C     
      DIMENSION GAME(4,37),BRA(37),RATIO(37),MULTI(37),EGAM(4)
      DATA GAME/9.000, 0., 0., 0.,  8.534, 0.466, 0.  ,0.  ,
     &     8.122, 0.878 , 0., 0. ,  7.698, 1.302 ,0.  ,0.  ,
     &     7.266, 0.856, 0.878, 0., 6.584, 1.950 ,0.466,0.   ,
     &     6.106, 2.894 , 0., 0.,   6.106, 2.555 ,0.340,0.   ,
     &     5.974, 3.026 , 0., 0.,   5.974, 2.686 ,0.340,0.   ,
     &     5.817, 3.182 , 0., 0.,   5.817, 2.842 ,0.340,0.   ,
     &     5.621, 3.039 , 0., 0.,   5.437, 3.563 ,0. , 0.   ,
     &     5.437, 2.686, 0.878,0.,  5.313, 3.347 ,0.340, 0. ,
     &     5.313, 3.221, 0.466,0.,  5.313, 3.687 ,0.878, 0. ,
     &     5.270, 3.266, 0.466,0.,  4.977, 3.555 ,0.466 ,0. ,
     &     4.859, 4.141, 0.   ,0.,  4.859, 3.675 ,0.466 ,0. ,
     &     4.284, 4.716, 0. , 0. ,  4.031, 4.969 ,0.    ,0. ,
     &     3.930, 4.192, 0.878,0.,  7.820, 0.    ,0.    ,0. ,
     &     7.537, 0.283, 0.   ,0.,  6.720, 1.031 ,0.    ,0. ,
     &     6.720, 0.816, 0.284,0.,  6.635, 1.185 ,0.    ,0. ,
     &     5.696, 2.124, 0.   ,0.,  6.838, 0.    ,0.    ,0. ,
     &     6.682, 0.156, 0.   ,0.,  6.320, 0.362 ,0.156 ,0. ,
     &     5.837, 0.846, 0.156,0.,  5.837,0.483,0.362, 0.156,
     &     5.515, 1.168 , 0.156 , 0.   /
      DATA MULTI/1, 2, 2, 2, 3, 3, 2, 3, 2, 3, 2, 3, 2, 2, 3, 3,
     &           3, 3, 3, 3, 2, 3, 2, 2, 3, 1, 2, 2, 3, 2, 2, 1,
     &           2, 3, 3, 4, 3/
      DATA BRA/ 34.26, 16.64, 3.023, 0.832, 0.2015, 1.781, 0.2825,
     &     1.271, 0.3274, 0.2381, 0.5051, 1.900, 0.104, 0.1957,
     &     0.2723, 0.3085, 0.6041, 0.2185, 0.1885, 0.1690, 0.2307,
     &     0.7638, 0.2340, 0.2990, 0.1610, 8.220, 4.800, 0.1485,
     &     0.1054, 0.1455, 0.9465, 10.31, 0.1764, 0.4285, 0.3984,
     &     0.3114, 0.1012/
      DATA ICALL/0/

      DIMENSION rn(1)
    
      IF(ICALL.NE.0)GO TO 90
      ICALL=1
C     -- CALCULATE INTEGRAL RATIO
      TOTBRA=0.
      DO 10 I=1,37
         TOTBRA=TOTBRA+BRA(I)
   10 CONTINUE
      DO 20 K=1,36
         SUBBRA=0.
         DO 21 I=1,K
            SUBBRA=SUBBRA+BRA(I)
   21    CONTINUE
         RATIO(K)=SUBBRA/TOTBRA
   20 CONTINUE
      RATIO(37)=1.
C     
C     -- MAKE VECTOR
C     
   90 CONTINUE

      call grndm(rn,1) 
      RAN=RN(1)
      I=1
  100 IF(RAN.LT.RATIO(I))GO TO 110
      I=I+1
      GO TO 100
C     
  110 MODE=I
      NGAM=MULTI(MODE)
      DO 200 I=1,4
  200 EGAM(I)=GAME(I,MODE)

      RETURN
      END
