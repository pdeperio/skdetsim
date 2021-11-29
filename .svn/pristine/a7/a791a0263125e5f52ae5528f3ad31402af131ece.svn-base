      SUBROUTINE SKGAMMA(ngam,pgam,egam)

      integer ngam,mode,mode2
      real    pgam(3,2),egam(2),dirgam(3)
      real    gam1(7),gam21(4),gam22(4)
      dimension rndm(1)
      data gam1/ 10.07, 9.23, 8.31, 7.30, 6.32, 2.31, 5.30/
      data gam21/ 5.27, 5.30, 5.27, 2.31/
      data gam22/ 4.56, 3.93, 3.30, 1.64/

      do jj=1,2
         do ii=1,3
         pgam(ii,jj)=0.
         end do
         egam(jj)=0.
      end do
      ngam=0
      mode2=0
      mode=0
      call grndm(rndm,1)
C----------------------------------------------------------------------
C  added a correction factor 0.628(=24.5/39) to make room for the
C  ground-state transition for 15N. [ Physics Reports 354 (2001) 371 ]
C    by S.Nakayama (September 1, 2009)
C----------------------------------------------------------------------
      if(rndm(1).lt.0.684*0.628)mode2=4
      if(rndm(1).lt.0.670*0.628)mode2=3
      if(rndm(1).lt.0.647*0.628)mode2=2
      if(rndm(1).lt.0.623*0.628)mode2=1
      if(rndm(1).lt.0.599*0.628)mode=7
      if(rndm(1).lt.0.520*0.628)mode=6
      if(rndm(1).lt.0.518*0.628)mode=5
      if(rndm(1).lt.0.128*0.628)mode=4
      if(rndm(1).lt.0.088*0.628)mode=3
      if(rndm(1).lt.0.059*0.628)mode=2
      if(rndm(1).lt.0.030*0.628)mode=1

      if(mode.gt.0)then
         ngam=1
         egam(1)=gam1(mode)*0.001
         call rndir(dirgam)
         pgam(1,1)=dirgam(1)*egam(1)
         pgam(2,1)=dirgam(2)*egam(1)
         pgam(3,1)=dirgam(3)*egam(1)
      else if(mode2.gt.0)then
         ngam=2
         egam(1)=gam21(mode2)*0.001
         egam(2)=gam22(mode2)*0.001
         do ii=1,2
            call rndir(dirgam)
            pgam(1,ii)=dirgam(1)*egam(ii)
            pgam(2,ii)=dirgam(2)*egam(ii)
            pgam(3,ii)=dirgam(3)*egam(ii)
         end do
      end if
      return
      end
