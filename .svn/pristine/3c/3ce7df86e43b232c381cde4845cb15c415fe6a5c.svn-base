************************************************************************
*     -----------------
      SUBROUTINE SGPMTA
*     -----------------
*
*     (Purpose)
*
*     (Input)
*       none
*
*     (Output)
*       none
*
*     (Creation Date and Author)
*       1992.12.02 ; First version by Y.Koshio
*
************************************************************************
      implicit none
#include "skhead.h"
#include "geant321/gctrak.inc"
#include "geant321/gcvolu.inc"
#include "geant321/gckine.inc"
#include "skparm.h"
#include "geopmt.h"
#include "skmc.h"
#include "gupht.h"
#include "dsparma.h"

#include "geopma.h"
#include "guanti.h"
#include "guxyzpm.h" 
#include "skdetsim_parm.h"
#include "skdetsim_test.h"
#include "geotnk.h"

C Note on NAME, CNAME - NAME(15) from gcvolu.inc should be an array of ints
C If you want to look at it as a character*4, you need to call UHTOC first
C I did not change all the commented-out instances of this.  So, if you
C comment them back in, please take this extra step.
C See http://nucleus.usask.ca/Geant/Exercises/idiosyncracies.html
C midway down the page.


! skdonuts stuff
      real norm
      real inplane_vector(3)
      real dir(3)
      real pmadir(3)
      real costh,r
      real theta

      ! for pythagorean theorem
      real a,b,c
      real cosg, x 
      real gam
      real pmt_efficiency
      real coseffask

      character*4 cname
      character*4 cname_arr(15)
      integer nct,ncb,nopmta,locat,ictr,icabminde,icabde,i,ierr
      real r1,r2,rmindebug ,rdebug,r1new,r2new,acoreff,rdummy
      DIMENSION NCT(16)
      DATA NCT/ 0,1,3,6,9,12,16,20,24,28,32,35,38,41,43,44/
      DIMENSION NCB(16)
      DATA NCB/ 0,1,3,6,10,14,18,22,26,30,34,38,42,45,47,48 /

      INTEGER*4 NATOP1(14,2)
      DATA NATOP1 /
     $             9999,1347,1363,1377,1392,1408,1424,1440,1456,
     $             1472,1487,1503,1517,9999,
     $             9999,1334,1349,1364,1379,1395,1411,1427,1443,
     $             1459,1474,1489,1504,9999/

      INTEGER*4 NATOP2(24,2)
      DATA NATOP2 /
     $             9999,1545,1546,1547,1548,1549,1550,1551,1552,1553,
     $             1555,9999,9999,1556,1566,1558,1559,1560,1561,1562,
     $             1563,1573,1565,9999,
     $             9999,1288,1289,1290,1291,1292,1293,1294,1295,1296,
     $             1297,9999,9999,1298,1299,1300,1301,1302,1303,1304,
     $             1305,1306,1307,9999/

      INTEGER*4 NATOP3(10,2)
      DATA NATOP3 /
     $             9999,9999,1393,1409,1425,1441,1457,1473,9999,9999,
     $             9999,9999,1378,1394,1410,1426,1442,1458,9999,9999/


      INTEGER*4 NATOP4(16,2)
      DATA NATOP4 /
     $             9999,1567,1568,1569,1571,1577,1572,9999,
     $             9999,1574,9999,1575,1576,9999,9999,9999,
     $             9999,1276,1277,1278,1279,1280,1281,9999,
     $             9999,1282,1283,1284,1285,1286,1287,9999/

      INTEGER*4 NABOT1(14,2)
      DATA NABOT1 /1637,1651,1666,1682,1698,1714,1730,1746,1762,1778,
     $             1794,1810,1825,1839,
     $             1624,1638,1653,1669,1685,1701,1717,1733,1749,1765,
     $             1781,1797,1812,1826/

      INTEGER*4 NABOT2(24,2)
      DATA NABOT2 /
     $             1852,1853,1854,1855,1856,1857,1858,1859,1860,1861,
     $             1862,1863,9999,1864,1865,1866,1867,1868,1869,1870,
     $             1871,1872,1873,9999,
     $             9999,1590,1591,1592,1593,1594,1595,1596,1597,1598,
     $             1599,9999,1600,1601,1602,1603,1604,1605,1606,1607,
     $             1608,1609,1610,1611/
      
      INTEGER*4 NABOT3(10,2)
      DATA NABOT3 /1667,1683,1699,1715,1731,1747,1763,1779,1795,1811,
     $             1652,1668,1684,1700,1716,1732,1748,1764,1780,1796/

      INTEGER*4 NABOT4(16,2)
      DATA NABOT4 /
     $             1874,1875,1876,1877,1878,1879,1880,1881,
     $             9999,9999,1882,1883,1884,1885,9999,9999,
     $             9999,9999,1578,1579,1580,1581,9999,9999,
     $             1582,1583,1584,1585,1586,1587,1588,1589/
      

      integer IDATA(2)
      real    RDATA(2)
      EQUIVALENCE (IDATA(1),RDATA(1))
      integer iseg
      
*  Set parameter
 
      if( sk_geometry .lt. 4 ) then 
         call sgpmta_old
         goto 999 
      endif 


 1000 format(' ',15a5)
 1100 format(' ',15i5)

      NOPMTA=99999
      LOCAT=9999

*  If PMT is on the WALL

      CALL UHTOC(NAMES(4),4,CNAME,4)
      IF(CNAME.EQ."WANT".AND.NLEVEL.EQ.6) THEN
         NOPMTA = NUMBER(5) + (NUMBER(6)-1)*17
         
         r1 = sqrt((vect(1)-gxyzpma(1,nopmta))**2
     &             +(vect(2)-gxyzpma(2,nopmta))**2)
         r2 = abs(vect(3)-gxyzpma(3,nopmta))
         LOCAT=1
         
*  If PMT is on the TOP or BOTTOM

         GO TO 5800
      ELSE
         IF(CNAME.EQ."ANTB".AND.(NLEVEL.EQ.6.OR.NLEVEL.EQ.7))THEN
            
C*  for TOP
C         IF(NUMBER(4).EQ.1) NOPMTA = NTWALA + NOPMTA
C*  for BOTTOM
C         IF(NUMBER(4).EQ.2) NOPMTA = NLTOPA + NOPMTA

C*  for TOP
            IF(NUMBER(4).EQ.1) THEN
               CALL UHTOC(NAMES(5),4,CNAME,4)
               IF(CNAME.EQ."ANTC".AND.NLEVEL.EQ.7) THEN
                  NOPMTA = NUMBER(6)+(NUMBER(7)-1)*12+
     $                 NCT(NUMBER(7))+1307
                  IF (NOPMTA.eq.1332) NOPMTA=1333
                  IF (NOPMTA.eq.1502) NOPMTA=1503
                  IF (NOPMTA.eq.1530) NOPMTA=1531
                  LOCAT=2
                  GOTO 5000
               ENDIF
               
               IF(CNAME.EQ."ANT1".AND.NLEVEL.EQ.6)THEN
                  IF((NUMBER(5).EQ.1).OR.(NUMBER(5).EQ.2)) THEN
                     NOPMTA=NATOP1(NUMBER(6),NUMBER(5))
                     LOCAT=3
                  ENDIF
                  GOTO 5000
               ENDIF
               
               IF(CNAME.EQ."ANT2".AND.NLEVEL.EQ.7)THEN
                  IF ((NUMBER(5).eq.1).OR.(NUMBER(5).EQ.2)) THEN
                     NOPMTA=NATOP2(NUMBER(6)+(NUMBER(7)-1)*12,
     $                    NUMBER(5))
                     LOCAT=4
                  ENDIF
                  GOTO 5000
               ENDIF
               
               IF (CNAME.EQ."ANT3".AND.NLEVEL.EQ.6)THEN
                  IF ((NUMBER(5).EQ.1).OR.(NUMBER(5).EQ.2)) THEN
                     NOPMTA=NATOP3(NUMBER(6),NUMBER(5))
                     LOCAT=5
                  ENDIF
                  GO TO 5000
               ENDIF
               
               IF (CNAME.EQ."ANT4".AND.NLEVEL.EQ.7)THEN
                  IF ((NUMBER(5).EQ.1).OR.(NUMBER(5).EQ.2)) THEN
                     NOPMTA=NATOP4(NUMBER(6)+(NUMBER(7)-1)*8,
     $                    NUMBER(5))
                     LOCAT=6
                     GO TO 5000
                  ENDIF
               ENDIF
               write(*,*) "Incorrect Location!"
               GO TO 900
            ENDIF
C
C*  for BOTTOM

            IF(NUMBER(4).EQ.2) THEN
               CALL UHTOC(NAMES(5),4,CNAME,4)
               IF(CNAME.EQ."ANTC".AND.NLEVEL.EQ.7) THEN
                  NOPMTA = NUMBER(6)+(NUMBER(7)-1)*12+
     $                 NCB(NUMBER(7))+1611
                  LOCAT=2
                  GOTO 5000
               ENDIF
               
               IF(CNAME.EQ."ANT1".AND.NLEVEL.EQ.6)THEN
                  IF((NUMBER(5).EQ.1).OR.(NUMBER(5).EQ.2)) THEN
                     NOPMTA=NABOT1(NUMBER(6),NUMBER(5))
                     LOCAT=3
                  ENDIF
                  GOTO 5000
               ENDIF
               
               IF(CNAME.EQ."ANT2".AND.NLEVEL.EQ.7)THEN
                  IF ((NUMBER(5).eq.1).OR.(NUMBER(5).EQ.2)) THEN
                     NOPMTA=NABOT2(NUMBER(6)+(NUMBER(7)-1)*12,
     $                    NUMBER(5))
                     LOCAT=4
                  ENDIF
                  GOTO 5000
               ENDIF
               
               IF (CNAME.EQ."ANT3".AND.NLEVEL.EQ.6)THEN
                  IF ((NUMBER(5).EQ.1).OR.(NUMBER(5).EQ.2)) THEN
                     NOPMTA=NABOT3(NUMBER(6),NUMBER(5))
                     LOCAT=5
                  ENDIF
                  GO TO 5000
               ENDIF

               IF (CNAME.EQ."ANT4".AND.NLEVEL.EQ.7)THEN
                  IF ((NUMBER(5).EQ.1).OR.(NUMBER(5).EQ.2)) THEN
                     NOPMTA=NABOT4(NUMBER(6)+(NUMBER(7)-1)*8,
     $                    NUMBER(5))
                     LOCAT=6
                  ENDIF
                  GO TO 5000
               ENDIF
               write(*,*) "Incorrect Location!"
               GO TO 900
            ENDIF
            write(*,*) "Incorrect Location!"
            GO TO 900
         ENDIF
      ENDIF
         
 900  IF (NOPMTA.ne.9999) THEN
         IF (NLEVEL.eq.4) GO TO 999
         write(*,*) "SGPMTA:??What?? Dumping(and skip)/NOPMTA=",NOPMTA
         DO 910 ICTR=1,5
 910     CONTINUE
         DO 920 ICTR=1,7
            write(*,*) "Number(",ICTR,")=",NUMBER(ICTR)
 920     CONTINUE
      ENDIF
      GO TO 999
            

 5000 CONTINUE            

*  for top & bottom
      IF (nopmta.gt.MAXPMA) GO TO 900
      r1 = abs(vect(1)-gxyzpma(1,nopmta))
      r2 = abs(vect(2)-gxyzpma(2,nopmta))
 5500 CONTINUE


C
C--- search for minimum distance PMT
C

      IF ((IAGEQAL.eq.2).and.((NUMBER(4).eq.1).or.(NUMBER(4).eq.2)))
     $     THEN
C         write(*,*) "TOP/bottom check!"
         RMINDEBUG=99999999.9999
         ICABMINDE=99999
         IF (NUMBER(4).EQ.1) THEN
            DO 10 ICABDE=1276,1577
               RDEBUG= (vect(1)-gxyzpma(1,ICABDE))**2
     $              +(vect(2)-gxyzpma(2,ICABDE))**2
               IF (RDEBUG.le.RMINDEBUG) THEN
                  RMINDEBUG=RDEBUG
                  ICABMINDE=ICABDE
               ENDIF
 10         CONTINUE
         ENDIF
         IF (NUMBER(4).EQ.2) THEN
            DO 110 ICABDE=1578,MAXPMA
               RDEBUG= (vect(1)-gxyzpma(1,ICABDE))**2
     $              +(vect(2)-gxyzpma(2,ICABDE))**2
               IF (RDEBUG.le.RMINDEBUG) THEN
                  RMINDEBUG=RDEBUG
                  ICABMINDE=ICABDE
               ENDIF
 110        CONTINUE
         ENDIF
         
C         IF (NOPMTA.ne.ICABMINDE) THEN
            IF (ICABMINDE.ne.99999) THEN
                r1new = abs(vect(1)-gxyzpma(1,ICABMINDE))
                r2new = abs(vect(2)-gxyzpma(2,ICABMINDE))
                NOPMTA=ICABMINDE
                r1=r1new
                r2=r2new
                GO TO 6000
            ELSE
               write(*,*) "SGPMTA:Illegal Cable/Location!"
               GO TO 900
            ENDIF
C         ENDIF
C---  
         
         
*  Check if Cherenkov is in insensitive regions

      ELSE 
         CALL UHTOC(NAMES(3),4,CNAME,4)
         IF(CNAME.EQ."WMED".AND.NLEVEL.EQ.3) GOTO 999
         CALL UHTOC(NAMES(4),4,CNAME,4)
         IF(CNAME.EQ."ANTB".AND.NLEVEL.EQ.4) GOTO 999
      ENDIF
      
C      WRITE(*,*) "P"
      GO TO 999
      
 5800 CONTINUE
!
!-Search for the nearest cable
! on the wall
!
      IF (IAGEQAL.ne.2)  GO TO 6000
      RMINDEBUG=99999999.9999
      ICABMINDE=99999
      DO 20 ICABDE=1,1275
         RDEBUG = (sqrt((vect(1)-gxyzpma(1,ICABDE))**2
     $             +(vect(2)-gxyzpma(2,ICABDE))**2))**2
     $           +(vect(3)-gxyzpma(3,ICABDE))**2

         IF (RDEBUG.le.RMINDEBUG) THEN
            RMINDEBUG=RDEBUG
            ICABMINDE=ICABDE
         ENDIF
 20   CONTINUE

      IF (NOPMTA.ne.ICABMINDE) THEN
         IF (ICABMINDE.ne.99999) THEN
            r1new = sqrt((vect(1)-gxyzpma(1,nopmta))**2
     &             +(vect(2)-gxyzpma(2,nopmta))**2)
            r2new = abs(vect(3)-gxyzpma(3,nopmta))
            if ((r1new.le.30.).and.(r2new.le.30.)) THEN         
               write(*,"(A10,F8.0,A1,F8.0,A1)") 
     $              "Location:(",vect(1),",",vect(2),")"
               write(*,"(A14,I5,A1,F8.0,A1,F8.0,"
     $                 //"A13,I5,A1,F8.0,A1,F8.0,A2)")
     $              "SGPMTA:NOPMTA=",NOPMTA,
     $              "(",r1,",",r2,",) /NEAREST=",ICABMINDE,
     $              "(",r1new,",",r2new,",)"
               write(*,"(A14,I5,A1,F8.0,A1,F8.0,"
     $                 //"A13,I5,A1,F8.0,A1,F8.0,A2)")
     $              "SGPMTA:NOPMTA=",NOPMTA,
     $              "(",gxyzpma(1,NOPMTA),",",gxyzpma(2,NOPMTA),
     $              ",) /NEAREST=",ICABMINDE,
     $              "(",gxyzpma(1,ICABMINDE),
     $              ",",gxyzpma(2,ICABMINDE),",)"
               NOPMTA=ICABMINDE
               GO TO 6000
            ENDIF
         ELSE
            write(*,*) "SGPMTA:Illegal Cable/Location!"
            GO TO 900
         ENDIF
      ENDIF



 6000 CONTINUE

      IF(NOPMTA.EQ.0)THEN
         write(6,*) ' ' 
         print*,nlevel,(vect(i),i=1,3)
         do i=1,15
            CALL UHTOC(NAMES(i),4,CNAME_ARR(i),4)
         end do
         write(6,1000)(cname_arr(i),i=1,15)
         write(6,1100)(number(i),i=1,15)
         WRITE(6,*) 'Warning in SGPMTA'
c     STOP
         goto 999
      ENDIF
      
      IF(NOPMTA.GT.MAXPMA) GOTO 999
      
*     Check if the photon is in anti PMT?
!   
!  -- wall tubes
!    r1 = distance in x-y plane 
!    r2 = z-distance
!
!  -- top / bottom tubes
!    r1 = x-distance 
!    r2 = y-distance
!
! 


      
c     print*,'Catch'
*     Input Data into Common
      
      if( r1.gt.wlssz .or. r2.gt.wlssz ) then ! nearest cable is on the wall 
         call sgtyvk
         goto 999
!     endif
      else

      ! 8inch tube's radius in cm  is (4.*2.54)
!     IF ((r1**2+r2**2).le.(4.*2.54)**2) THEN

         pmadir(1)  = 0.
         pmadir(1)  = 0.
         pmadir(3)  = 0.
         ! set tube direction
         if( nopmta .le. 1275  ) then ! wall pmt
            ! normal faces outward
            pmadir(1)  = gxyzpma(1, nopmta ) / RINTK + RMED
            pmadir(1)  = gxyzpma(2, nopmta ) / RINTK + RMED

         else if( nopmta .gt. 1275  .and. nopmta .le. 21577 ) then ! top pmt
            pmadir(3)  = 1.
         else if( nopmta .gt. 1577 ) then ! bottom pmt
            pmadir(3)  = -1.
         endif 

!      compute intersection with sphere ( project back )
!      due to lame geometry, photon his is already confined
!      to tank surface in this routine, so we are already in-plane 
         inplane_vector(1) =  vect(1) - gxyzpma(1, nopmta )
         inplane_vector(2) =  vect(2) - gxyzpma(2, nopmta )
         inplane_vector(3) =  vect(3) - gxyzpma(3, nopmta )

         norm = inplane_vector(1)**2 
     &        + inplane_vector(2)**2
     &        + inplane_vector(3)**2
         norm = sqrt( norm )


         ! normalize
         inplane_vector(1) = inplane_vector(1) / ( norm )
         inplane_vector(2) = inplane_vector(2) / ( norm )
         inplane_vector(3) = inplane_vector(3) / ( norm )

         ! get direction of incoming photon
         dir(1) = vect(4)
         dir(2) = vect(5)
         dir(3) = vect(6)
         
         ! dot-product of vector to hit location from tube
         ! ceterpoint with incoming particle direction
         costh = (dir(1)*inplane_vector(1))
     &         + (dir(2)*inplane_vector(2))
     &         + (dir(3)*inplane_vector(3))

         ! r is 4.0 inch radius of pmt sphere converted to 
         theta = acos( costh )
         r = 4.0 * 2.54
         
         ! from law of cosines
         a = 1.
         b = -2*norm*costh
         c = norm**2-r**2
         
         !
         ! ::  x represents distance from  
         !     hit point in plane to pmt sphere
         !     
         call grndm(rdummy,1)
         if( b**2 - 4*a*c  > 0 ) then 
            x = -b + sqrt( b**2 - 4*a*c )
            x = x / (2*a) 
            ! cosine gamma, the angle with the normal 
            ! at the PMT surface
            cosg = (x**2 + r**2 - norm**2)/ (2*x*r) 
            gam = acos( cosg )
           !write(6,*) ' Inside PMT'  , b**2 - 4*a*c, r1, cosg, costh 
            pmt_efficiency = coseffask( cosg ) 
            if( rdummy .lt. pmt_efficiency ) THEN
               ACOREFF=APMTCOR
            else 
               goto 999
            endif
         else
            ! this track doesn't intersect 
            ! the pmt sphere
            ! use the angle with the normal to the plane
            ! --- need a separate function for the WLS?
            norm = pmadir(1)**2 
     &           + pmadir(2)**2
     &           + pmadir(3)**2
            norm = sqrt( norm )

            costh = (dir(1)*pmadir(1))
     &            + (dir(2)*pmadir(2))
     &            + (dir(3)*pmadir(3))
            costh = abs(costh)/norm

           !write(6,*) ' Inside WLS'  , b**2 - 4*a*c, r1, costh, acos(costh) 
            pmt_efficiency = coseffask( costh ) 
            if( rdummy .lt. pmt_efficiency ) THEN
               ACOREFF=AWLSCOR
               tofg = tofg + 3.0e-9 ! delay by a few ns if its on the WLS
            else 
               goto 999
            endif
          endif ! end of intersection with pmt test
      ENDIF ! end of test to see if point landed in the PMT-WLS area
      
      CALL GRNDM(RDUMMY,1)
      IF (RDUMMY.lt.ACOREFF) THEN

       if( sk_geometry.eq.sk_i  .or. sk_geometry.eq.sk_ii .or. 
     +     sk_geometry.eq.sk_iii ) then

         IF(NQGESKA(NOPMTA,ITRA).EQ.0)THEN
            TGESKA(NOPMTA,ITRA)=TOFG
         ELSE
            IF(TOFG.LE.TGESKA(NOPMTA,ITRA)) 
     $           TGESKA(NOPMTA,ITRA)=TOFG
         ENDIF
         NQGESKA(NOPMTA,ITRA)=NQGESKA(NOPMTA,ITRA)+1
         
         IF(OUTPUT_HBOOK .ge. OUT_TEST) then
            TEST_Colla_phot=TEST_Colla_phot+1
            TEST_Hitmapa(NOPMTA)=TEST_Hitmapa(NOPMTA)+1.
         endif
       else if( sk_geometry.ge.sk_iv )then
       
        if( od_photon .eq. count_all_od_photon )then 
         nqgeska(1, 1 ) = nqgeska(1,1)+1
         iseg     = nqgeska(1,1)
         idata(1) = nopmta
         rdata(2) = tofg
         call kzrep2('MCPHOTAN',1, iseg , 2, idata ) 
        else
          call sgcnt_aphoton( nopmta )
        endif

       endif ! end of check on skgeometry

      ENDIF
 999  RETURN
      END
      
************************************************************************
*     ------------------------------
      real function COSEFFASK(cosang)
*     ------------------------------
*
*     (Purpose)
*       EFFICIENCY OF PMT
*
*     (Input)
*       COSANG : DIRECTION ( COS(TH) )
*
*     (Output)
*       coseffask : EFFICIENCY
*      
*    added for anti-pmts by rvw ( from sklib )
*    -- needs updating to correct function, currently using
*    -- function for ID tubes at SK-II 
*
************************************************************************
      implicit none
#include "skhead.h"
#include "skprint.h"
      real x, cosang
      real sk1_a(4), sk2_a(4)

      data sk1_a /0.3542135,0.5104711,0.1160429,0.0118482/
      data sk2_a /0.205349,0.523981,0.389951,-0.131959/

      if ((cosang .ge. 0.0) .and. (cosang .le. 1.0)) then
         x = cosang
      else if (cosang .gt. 1.0) then
         if ( skprint(sklib_io,1) ) then
            write(6,*) 'coseffask: abnormal cosang!!'
            write(6,*) 'cosang = ', cosang
         endif
         x = 1.0
      else
         if ( skprint(sklib_io,1) ) then
            write(6,*) 'coseffask: abnormal cosang!!'
            write(6,*) 'cosang = ', cosang
         endif
         x = 0.0
      end if

      ! we only need one geometry - lets use the one with 
      ! acryllic case, less acceptance closer to od tubes?
         coseffask = sk2_a(1) + sk2_a(2)*x + sk2_a(3)*x**2 +
     &        sk2_a(4)*x**3


      RETURN
      END







