************************************************************************
*     -------------------------------------------------
      SUBROUTINE SGPMTNUM(NAMES,n5,n6,n7,n8,NOPMT,iflag)
*     -------------------------------------------------
*
*     (Purpose)
**     This routine is used to determine which pmt site a GEANT 
*     volume corresponds to.
*     This routine determines the pmt cable number from the GEANT volume's 
*     position in the geometrical heirarchy, as defined by NAMES 
*     and NUMBER(*).  As the Acrylic cover, Water gap, and inner pmt at a 
*     location are all daughters of the same mother, this routine returns the
*     same cable number for all three volumes at a pmt site.  
*       iflag notes if this volume corresponds to an irregularly 
*     removed pmt, i.e. if this pmt is defined in GEANT but is not in SK-II.
*     
*     (Input)
*       NAMES ; name of geometry (5th array)
*       n5    ; number(5) in gcvol
*       n6    ; number(6) in gcvol
*       n7    ; number(7) in gcvol
*       n8    ; number(8) in gcvol
*
*     (Output)
*       NOPMT ; cable number
*       iflag ; 1:pmt exist, 0:pmt missing
*
*     (Creation Date and Author)
*       2002.10.23 ; First version by Y.Koshio
*       2006.4.18  ; added SK-III case temporarily with an 'or'. needs
*                    more work.  K Bays
*
************************************************************************

      implicit none
#include "skdetsim_parm.h"
#include "skparm.h"
#include "geopmt.h"
#include "skhead.h"


      character*4 names
      integer n5,n6,n7,n8,nopmt,iflag
      character*256 fname

      integer i1_1,i1_2,i1_3,i1_4
      parameter (i1_1=1,i1_2=4,i1_3=6,i1_4=18)
      integer i2top_1,i2bot_1,i2_2,i2_3
      parameter (i2top_1=1,i2bot_1=2,i2_2=4,i2_3=18)
      integer i3_1,i3_2,i3_3
      parameter (i3_1=1,i3_2=4,i3_3=9)
      integer i4_1,i4_2,i4_3
      parameter (i4_1=1,i4_2=12,i4_3=3)
      integer i5_1,i5_2,i5_3
      parameter (i5_1=2,i5_2=3,i5_3=3)
      integer i6_1,i6_2,i6_3
      parameter (i6_1=2,i6_2=4,i6_3=6)
      integer i7_top, i7_bot
      parameter (i7_top=137,i7_bot=125)

      integer i20_1,i20_2,i20_3
      parameter (i20_1=1,i20_2=2,i20_3=18)
      integer i21_1,i21_2
      parameter (i21_1=1,i21_2=18)
      integer i22_1,i22_2
      parameter (i22_1=1,i22_2=18)

      integer isave,i,ityp,i5,i6,i7,i8,ipmt,
     &     itop1(i1_1,i1_2,i1_3,i1_4),itop2(i2top_1,i2_2,i2_3),
     &     itop3(i3_1,i3_2,i3_3),itop4(i4_1,i4_2,i4_3),
     &     itop5(i5_1,i5_2,i5_3),itop6(i6_1,i6_2,i6_3),
     &     itop7(i7_top),itop20(i20_1,i20_2,i20_3),itop21(i21_1,i21_2),
     &     itop22(i22_1,i22_2),
     &     ibot1(i1_1,i1_2,i1_3,i1_4),ibot2(i2bot_1,i2_2,i2_3),
     &     ibot3(i3_1,i3_2,i3_3),ibot4(i4_1,i4_2,i4_3),
     &     ibot5(i5_1,i5_2,i5_3),ibot6(i6_1,i6_2,i6_3),
     &     ibot7(i7_bot)

      save isave, itop1,itop2,itop3,itop4,itop5,itop6,itop7,itop20,
     &     itop21,itop22,ibot1,ibot2,ibot3,ibot4,ibot5,ibot6,ibot7

C     Check detector geometry
      IF(SK_GEOMETRY.eq.SK_I) then !Good Old SK-I
         IF(NAMES.EQ."PMPW") THEN
            NOPMT = N6 + (N7-1)*51
            
         ELSE IF(NAMES.EQ."PMPT") THEN
            NOPMT = ISEQT(N7,N8)
            
         ELSE IF(NAMES.EQ."PMPB") THEN
            NOPMT = ISEQB(N7,N8)
         ENDIF

         IFLAG=1
         IF(NOPMT.EQ.0) THEN
            IFLAG=0
         ENDIF
      ELSE IF((SK_GEOMETRY.eq.SK_II) .OR.
     &        (SK_GEOMETRY.eq.SK_III) .OR.
     &        (SK_GEOMETRY.eq.SK_IV) .OR.
     &        (SK_GEOMETRY.eq.SK_V) ) then                 ! after SK-II
 
         if(NAMES.NE.'PMPW'.and.NAMES.NE.'PMP1'.and.
     &     NAMES.NE.'PMP2'.and.NAMES.NE.'PMP3'.and.
     &     NAMES.NE.'PMP4'.and.NAMES.NE.'PMP5'.and.
     &     NAMES.NE.'PBXT'.and.NAMES.NE.'PBX3'.and.
     &     NAMES.NE.'PBX4'.and.NAMES.NE.'PBX5'.and.
     &     NAMES.NE.'PBX6'.and.NAMES.NE.'PBX7'.and.
     &     NAMES.NE.'PXT8'.and.NAMES.NE.'PBX9'.and.
     &     NAMES.NE.'PXT0'.and.NAMES.NE.'PXTT'.and.
     &     NAMES.NE.'PBXB'.and.NAMES.NE.'PBB3'.and.
     &     NAMES.NE.'PBB4'.and.NAMES.NE.'PBB5'.and.
     &     NAMES.NE.'PBB6'.and.NAMES.NE.'PBB7'.and.
     &     NAMES.NE.'PXB8') then
         print*,'Error in sgpmtnum wrong NAMES',names
         stop
      endif

      if(isave.eq.0) then
         call findconsts(fname,'Conversion from Geant Volume to Cable number',
     &        'skdetsim_sk2_pmtinf.dat',0)
         open(99,file=fname,status='old')
         do i=1,1617
            read(99,*) ipmt,ityp,i5,i6,i7,i8
            if(ityp.eq.1) then
               itop1(i5,i6,i7,i8) = ipmt
            else if(ityp.eq.2) then
               itop2(i5,i6,i7) = ipmt
            else if(ityp.eq.3) then
               itop3(i5,i6,i7) = ipmt
            else if(ityp.eq.4) then
               itop4(i5,i6,i7) = ipmt
            else if(ityp.eq.5) then
               itop5(i5,i6,i7) = ipmt
            else if(ityp.eq.6) then
               itop6(i5,i6,i7) = ipmt
            else if(ityp.eq.7) then
               itop7(i5) = ipmt

            else if(ityp.eq.20) then
               itop20(i5,i6,i7) = ipmt
            else if(ityp.eq.21) then
               itop21(i5,i6) = ipmt
            else if(ityp.eq.22) then
               itop22(i5,i6) = ipmt

            else if(ityp.eq.8) then
               ibot1(i5,i6,i7,i8) = ipmt
            else if(ityp.eq.9) then
               ibot2(i5,i6,i7) = ipmt
            else if(ityp.eq.10) then
               ibot3(i5,i6,i7) = ipmt
            else if(ityp.eq.11) then
               ibot4(i5,i6,i7) = ipmt
            else if(ityp.eq.12) then
               ibot5(i5,i6,i7) = ipmt
            else if(ityp.eq.13) then
               ibot6(i5,i6,i7) = ipmt
            else if(ityp.eq.14) then
               ibot7(i5) = ipmt
            endif
         enddo
         close(99)

         isave=1
      endif

** for barrel
      if(NAMES.eq.'PMPW') then
         NOPMT = 204 + (n7-1)*1836 + (n8-1)/3*51
     &        + ((n6-1)*6+1) + mod((n8-1),3)*2
         if(mod(((n8-1)/3),2).eq.0) then
            NOPMT = NOPMT + 2
         else
            NOPMT = NOPMT + 1
         endif

      else if(NAMES.eq.'PMP1') then
         NOPMT = 5712 + (n7-1001)/3*51
     &        + ((n6-1)*6+1) + mod((n7-1001),3)*2
         if(mod(((n7-1001)/3),2).eq.0) then
            NOPMT = NOPMT + 2
         else
            NOPMT = NOPMT + 1
         endif

      else if(NAMES.eq.'PMP2') then
         NOPMT = 5814 + (n7-2001)/3*51
     &        + ((n6-1)*6+1) + mod((n7-2001),3)*2
         if(mod(((n7-2001)/3),2).eq.0) then
            NOPMT = NOPMT + 2
         else
            NOPMT = NOPMT + 1
         endif
         if(NOPMT.ge.7651) then
            NOPMT = NOPMT - 7650
         endif

      else if(NAMES.eq.'PMP3') then
** for edge correction on bottom part
         if(n6.eq.1) then
            NOPMT = 9580
         else if(n6.eq.2) then
            NOPMT = 9733
         else if(n6.eq.3) then
            NOPMT = 9994
         else if(n6.eq.4) then
            NOPMT = 10225
         else if(n6.eq.5) then
            NOPMT = 10812
         else
            NOPMT = (n6-3001)*51 + 1
         endif

      else if(NAMES.eq.'PMP4') then
         NOPMT = (n6-4000)*51 - 1

      else if(NAMES.eq.'PMP5') then
** for edge correction on top part
         if(n6.eq.1) then
            NOPMT = 7652
         else if(n6.eq.2) then
            NOPMT = 7710
         else if(n6.eq.3) then
            NOPMT = 7737
         else if(n6.eq.4) then
            NOPMT = 7767
         else if(n6.eq.5) then
            NOPMT = 7869
         else if(n6.eq.6) then
            NOPMT = 8026
         else if(n6.eq.7) then
            NOPMT = 8069
         else if(n6.eq.8) then
            NOPMT = 8524
         else if(n6.eq.9) then
            NOPMT = 8525
         else if(n6.eq.10) then
            NOPMT = 8848
         else if(n6.eq.11) then
            NOPMT = 9023
         else if(n6.eq.12) then
            NOPMT = 9180
         else if(n6.eq.13) then
            NOPMT = 9216
         else if(n6.eq.14) then
            NOPMT = 9250
         else if(n6.eq.15) then
            NOPMT = 9396
         else if(n6.eq.16) then
            NOPMT = 9397
         else
            NOPMT = (n6-5000)*51
         endif

** for top
      else if(NAMES.eq.'PBXT') then
         if(n5.gt.i1_1.or.n6.gt.i1_2.or.n7.gt.i1_3.or.n8.gt.i1_4)then
            print*,'Error in sgpmtnum ',n5,n6,n7,n8
            stop
         endif
         NOPMT = itop1(n5,n6,n7,n8)

      else if(NAMES.eq.'PBX3') then
         if(n5.gt.i2top_1.or.n6.gt.i2_2.or.n7.gt.i2_3) then
            print*,'Error in sgpmtnum ',n5,n6,n7
            stop
         endif
         NOPMT = itop2(n5,n6,n7)

      else if(NAMES.eq.'PBX4') then
         if(n5.gt.i3_1.or.n6.gt.i3_2.or.n7.gt.i3_3) then
            print*,'Error in sgpmtnum ',n5,n6,n7
            stop
         endif
         NOPMT = itop3(n5,n6,n7)

      else if(NAMES.eq.'PBX5') then
         if(n5.gt.i4_1.or.n6.gt.i4_2.or.n7.gt.i4_3) then
            print*,'Error in sgpmtnum ',n5,n6,n7
            stop
         endif
         NOPMT = itop4(n5,n6,n7)

      else if(NAMES.eq.'PBX6') then
         if(n5.gt.i5_1.or.n6.gt.i5_2.or.n7.gt.i5_3) then
            print*,'Error in sgpmtnum ',n5,n6,n7
            stop
         endif
         NOPMT = itop5(n5,n6,n7)

      else if(NAMES.eq.'PBX7') then
         if(n5.gt.i6_1.or.n6.gt.i6_2.or.n7.gt.i6_3) then
            print*,'Error in sgpmtnum ',n5,n6,n7
            stop
         endif
         NOPMT = itop6(n5,n6,n7)

      else if(NAMES.eq.'PXT8') then
         if(n5.gt.i7_top) then
            print*,'Error in sgpmtnum ',n5
            stop
         endif
         NOPMT = itop7(n5)

      else if(NAMES.eq.'PBX9') then
         if(n5.gt.i20_1.or.n6.gt.i20_2.or.n7.gt.i20_3) then
            print*,'Error in sgpmtnum ',n5,n6,n7
            stop
         endif
         NOPMT = itop20(n5,n6,n7)
      else if(NAMES.eq.'PXT0') then
         if(n5.gt.i21_1.or.n6.gt.i21_2) then
            print*,'Error in sgpmtnum ',n5,n6,n7
            stop
         endif
         NOPMT = itop21(n5,n6)

      else if(NAMES.eq.'PXTT') then
         if(n5.gt.i22_1.or.n6.gt.i22_2) then
            print*,'Error in sgpmtnum ',n5,n6,n7
            stop
         endif
         NOPMT = itop22(n5,n6)

** for bottom
      else if(NAMES.eq.'PBXB') then
         if(n5.gt.i1_1.or.n6.gt.i1_2.or.n7.gt.i1_3.or.n8.gt.i1_4)then
            print*,'Error in sgpmtnum ',n5,n6,n7,n8
            stop
         endif
         NOPMT = ibot1(n5,n6,n7,n8)

      else if(NAMES.eq.'PBB3') then
         if(n5.gt.i2bot_1.or.n6.gt.i2_2.or.n7.gt.i2_3) then
            print*,'Error in sgpmtnum ',n5,n6,n7
            stop
         endif
         NOPMT = ibot2(n5,n6,n7)

      else if(NAMES.eq.'PBB4') then
         if(n5.gt.i3_1.or.n6.gt.i3_2.or.n7.gt.i3_3) then
            print*,'Error in sgpmtnum ',n5,n6,n7
            stop
         endif
         NOPMT = ibot3(n5,n6,n7)

      else if(NAMES.eq.'PBB5') then
         if(n5.gt.i4_1.or.n6.gt.i4_2.or.n7.gt.i4_3) then
            print*,'Error in sgpmtnum ',n5,n6,n7
            stop
         endif
         NOPMT = ibot4(n5,n6,n7)

      else if(NAMES.eq.'PBB6') then
         if(n5.gt.i5_1.or.n6.gt.i5_2.or.n7.gt.i5_3) then
            print*,'Error in sgpmtnum ',n5,n6,n7
            stop
         endif
         NOPMT = ibot5(n5,n6,n7)

      else if(NAMES.eq.'PBB7') then
         if(n5.gt.i6_1.or.n6.gt.i6_2.or.n7.gt.i6_3) then
            print*,'Error in sgpmtnum ',n5,n6,n7
            stop
         endif
         NOPMT = ibot6(n5,n6,n7)

      else if(NAMES.eq.'PXB8') then
         if(n5.gt.i7_bot) then
            print*,'Error in sgpmtnum ',n5
            stop
         endif
         NOPMT = ibot7(n5)

      else
         print*,'Error in sgpmtnum ',names
         stop
      endif

      if(NOPMT.eq.3.or.NOPMT.eq.27.or.NOPMT.eq.633.or.NOPMT.eq.1185.or.
     &     NOPMT.eq.1263.or.NOPMT.eq.1383.or.NOPMT.eq.1407.or.
     &     NOPMT.eq.2013.or.NOPMT.eq.2565.or.NOPMT.eq.2643.or.
     &     NOPMT.eq.3393.or.NOPMT.eq.3945.or.NOPMT.eq.4143.or.
     &     NOPMT.eq.4167.or.NOPMT.eq.4695.or.NOPMT.eq.4773.or.
     &     NOPMT.eq.5325.or.NOPMT.eq.5403.or.NOPMT.eq.5523.or.
     &     NOPMT.eq.5547.or.NOPMT.eq.6075.or.NOPMT.eq.6153.or.
     &     NOPMT.eq.6705.or.NOPMT.eq.6903.or.NOPMT.eq.6927.or.
     &     NOPMT.eq.7455.or.NOPMT.eq.7533.or.
     &
     &     NOPMT.eq.9494.or.NOPMT.eq.9506.or.NOPMT.eq.9733.or.
     &     NOPMT.eq.9954.or.NOPMT.eq.9990.or.NOPMT.eq.10510.or.
     &     NOPMT.eq.10546.or.NOPMT.eq.10812.or.NOPMT.eq.11009.or.
     &     NOPMT.eq.11021) then
         iflag = 0
      else
         iflag = 1
      endif

      ELSE !SK-III not made yet
         write(*,*) 'Error SGPMTNUM wrong SK_GEOMETRY flag',SK_GEOMETRY
         stop
      ENDIF

      RETURN
      END
