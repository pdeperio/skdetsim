************************************************************************
*     -----------------
      SUBROUTINE SGPMTA_OLD
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

C Note on NAME, CNAME - NAME(15) from gcvolu.inc should be an array of ints
C If you want to look at it as a character*4, you need to call UHTOC first
C I did not change all the commented-out instances of this.  So, if you
C comment them back in, please take this extra step.
C See http://nucleus.usask.ca/Geant/Exercises/idiosyncracies.html
C midway down the page.
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

c      print*,nlevel
c      write(6,1000)(names(i),i=1,15)
c      write(6,1100)(number(i),i=1,15)
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
            
C         IF(NAMES(5).EQ."ANTC".AND.NLEVEL.EQ.7)THEN
C            NOPMTA = NUMBER(6)+(NUMBER(7)-1)*12+32+NC(NUMBER(7))
C         ELSE IF(NAMES(5).EQ."ANT1".AND.NLEVEL.EQ.6)THEN
C            IF(NUMBER(5).EQ.1)THEN
C               IF(NUMBER(6).EQ. 1) NOPMTA = 70
C               IF(NUMBER(6).EQ. 2) NOPMTA = 84
C               IF(NUMBER(6).GE.3.AND.NUMBER(6).LE.10) 
C     &              NOPMTA = 99 + (NUMBER(6)-3)*16
C               IF(NUMBER(6).EQ.11) NOPMTA = 226
C               IF(NUMBER(6).EQ.12) NOPMTA = 240
C            ELSE IF(NUMBER(5).EQ.2)THEN
C               IF(NUMBER(6).EQ. 1) NOPMTA = 57
C               IF(NUMBER(6).EQ. 2) NOPMTA = 71
C               IF(NUMBER(6).GE.3.AND.NUMBER(6).LE.10) 
C     &              NOPMTA = 86 + (NUMBER(6)-3)*16
C               IF(NUMBER(6).EQ.11) NOPMTA = 213
C               IF(NUMBER(6).EQ.12) NOPMTA = 227
C            ENDIF
C
C         ELSE IF(NAMES(5).EQ."ANT2".AND.NLEVEL.EQ.7)THEN
C            IF(NUMBER(5).EQ.1)THEN
C               NOPMTA = 264 + NUMBER(6) + (NUMBER(7)-1)*10
C            ELSE IF(NUMBER(5).EQ.2)THEN
C               NOPMTA =  12 + NUMBER(6) + (NUMBER(7)-1)*10
C            ENDIF
C
C         ELSE IF(NAMES(5).EQ."ANT3".AND.NLEVEL.EQ.6)THEN
C            IF(NUMBER(5).EQ.1)THEN
C               IF(NUMBER(6).GE.1.AND.NUMBER(6).LE.8) 
C     &              NOPMTA = 100 + (NUMBER(6)-1)*16
C            ELSE IF(NUMBER(5).EQ.2)THEN
C               IF(NUMBER(6).GE.1.AND.NUMBER(6).LE.8) 
C     &              NOPMTA =  85 + (NUMBER(6)-1)*16
C            ENDIF
C
C         ELSE IF(NAMES(5).EQ."ANT4".AND.NLEVEL.EQ.7)THEN
C            IF(NUMBER(5).EQ.1)THEN
C               NOPMTA = 284 + NUMBER(6) + (NUMBER(7)-1)*6
C            ELSE IF(NUMBER(5).EQ.2)THEN
C               NOPMTA = NUMBER(6) + (NUMBER(7)-1)*6
C            ENDIF
C            
C         ENDIF
C
C*  for TOP
C         IF(NUMBER(4).EQ.1) NOPMTA = NTWALA + NOPMTA
C
C*  for BOTTOM
C         IF(NUMBER(4).EQ.2) NOPMTA = NLTOPA + NOPMTA
C
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
C         write(*,*) "NLEVEL=",NLEVEL
         IF (NLEVEL.eq.4) GO TO 999
         write(*,*) "SGPMTA:??What?? Dumping(and skip)/NOPMTA=",NOPMTA
         DO 910 ICTR=1,5
c     CALL UHTOC(NAMES(ICTR),4,CNAME_ARR(ICTR),4)
c     write(*,'(A,I,A,A)') "Names(",ICTR,")=",CNAME_ARR(ICTR)
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
C      if (LOCAT.gt.2) THEN
C         if ((r1.le.30.).and.(r2.le.30.)) THECN
C            write(*,*) "ANTI FRINGE REGION HIT! / ",nopmta
C         ELSE
C            write(*,*) "ANTI FRINGE REGION OUT OF BOUNDS! / ",nopmta
C         ENDIF
C      ENDIF
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
C               write(*,*) "r1new=",r1new,"/r2new=",r2new
C               if ((r1new.le.30.).and.(r2new.le.30.)) THEN
C                  IF (ABS(ICABMINDE-NOPMTA).le.16) GO TO 6000
C                  write(*,"(A10,F8,A1,F8,A1)") 
C     $                 "Location:(",vect(1),",",vect(2),")"
C                  write(*,"(A14,I5,A1,F8,A1,F8,A13,I5,A1,F8,A1,F8,A2)")
C     $                 "SGPMTA:NOPMTA=",NOPMTA,
C     $                 "(",r1,",",r2,",) /NEAREST=",ICABMINDE,
C     $                 "(",r1new,",",r2new,",)"
C                  write(*,"(A14,I5,A1,F8,A1,F8,A13,I5,A1,F8,A1,F8,A2)")
C     $                 "SGPMTA:NOPMTA=",NOPMTA,
C     $                 "(",xyzpma(1,NOPMTA),",",xyzpma(2,NOPMTA),
C     $                 ",) /NEAREST=",ICABMINDE,
C     $                 "(",xyzpma(1,ICABMINDE),
C     $                 ",",xyzpma(2,ICABMINDE),",)"
C                  NOPMTA=ICABMINDE
C                  GO TO 900
                NOPMTA=ICABMINDE
                r1=r1new
                r2=r2new
C                if ((r1.le.30.).and.(r2.le.30.)) THEN
C                write(*,*) "CORRECT CABLE=",NOPMTA,"(",r1,",",r2,")"
C                ENDIF

                GO TO 6000
C               ENDIF
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
C
C-Search for the nearest cable
C
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
C               GO TO 900
               GO TO 6000
            ENDIF
         ELSE
            write(*,*) "SGPMTA:Illegal Cable/Location!"
            GO TO 900
         ENDIF
      ENDIF



 6000 CONTINUE

      IF(NOPMTA.EQ.0)THEN
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
      if(r1.gt.wlssz.or.r2.gt.wlssz) then
         call sgtyvk
         goto 999
      endif
      
c     print*,'Catch'
*     Input Data into Common
      
      IF ((r1**2+r2**2).le.(4.*2.54)**2) THEN
         ACOREFF=APMTCOR
      ELSE
         ACOREFF=AWLSCOR
      ENDIF
      
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
      
      
