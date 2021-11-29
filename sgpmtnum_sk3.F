************************************************************************
*     -------------------------------------------------
      SUBROUTINE SGPMTNUM_SK3(NAMES,n5,n6,n7,n8,NOPMT,iflag)
*     -------------------------------------------------
*     
*     (Purpose)
*     This routine is used to determine which pmt site a GEANT 
*     volume corresponds to.
*     This routine determines the pmt cable number from the GEANT volume's 
*     position in the geometrical heirarchy, as defined by NAMES 
*     and NUMBER(*).  As the Acrylic cover, Water gap, and inner pmt at a 
*     location are all daughters of the same mother, this routine returns the
*     same cable number for all three volumes at a pmt site.  
*     iflag notes if this volume corresponds to an irregularly 
*     removed pmt, i.e. if this pmt is defined in GEANT but is not in SK-III.
*     
*     (Input)
*     NAMES ; name of geometry (5th array)
*     n5    ; number(5) in gcvol
*     n6    ; number(6) in gcvol
*     n7    ; number(7) in gcvol
*     n8    ; number(8) in gcvol
*     
*     (Output)
*     NOPMT ; cable number
*     iflag ; 1:pmt exist, 0:pmt missing
*     
*     (Creation Date and Author)
*     2002.10.23 ; First version by Y.Koshio
*     2006.04.15 ; Split sgpmtnum.F into 3 routines, i.e. into SK-I,II,III
*                     by G.Mitsuka
*
************************************************************************

      implicit none
#include "skdetsim_parm.h"
#include "skparm.h"
#include "geopmt.h"
#include "skhead.h"


      character*4 names
      integer n5,n6,n7,n8,nopmt,iflag

      IF(NAMES.NE.'PMPW'.and.NAMES.NE.'PMP1'.and.
     &     NAMES.NE.'PMP3'.and.
     &     NAMES.NE.'PMP4'.and.NAMES.NE.'PMP5'.and.
     &     NAMES.NE.'PMPT'.and.NAMES.NE.'PMPB') then
         print*,'Error in sgpmtnum_sk3 wrong NAMES',names
         stop
      endif

** for barrel
      if(NAMES.eq.'PMPW') then
         if(n7.ge.1.and.n7.lt.10) then ! from SM30 to SM38
            NOPMT = (n7-1)*204 + (n8-1)/3*51
     &           + ((n6-1)*3+1) + mod((n8-1),3) + 1 + 5814
         else if(n7.ge.10) then ! from SM1 to SM28
            NOPMT = 204 + (n7-1)*204 + (n8-1)/3*51
     &           + ((n6-1)*3+1) + mod((n8-1),3) - 2039
         end if
      else if(NAMES.eq.'PMP1') then ! SM29
         NOPMT = 5712 + (n7-1001)/3*51
     &        + ((n6-1)*3+1) + mod((n7-1001),3) + 1
      else if(NAMES.eq.'PMP3') then
** for edge correction on bottom part
         if(n6.eq.1) then
            NOPMT = 10320
         else if(n6.eq.2) then
            NOPMT = 10596
         else if(n6.eq.3) then
            NOPMT = 10812
         else if(n6.eq.4) then
            NOPMT = 10928
         else if(n6.eq.5) then
            NOPMT = 10964
         else if(n6.eq.6) then
            NOPMT = 10998
         else if(n6.eq.7) then
            NOPMT = 11030
         else if(n6.eq.8) then
            NOPMT = 11130
         else if(n6.eq.9) then
            NOPMT = 11144
         else if(n6.eq.10) then
            NOPMT = 11146
         else if(n6.eq.11) then
            NOPMT = 11145
         else if(n6.eq.12) then
            NOPMT = 11131
         else if(n6.eq.13) then
            NOPMT = 11111
         else if(n6.eq.14) then
            NOPMT = 10999
         else if(n6.eq.15) then
            NOPMT = 10965
         else if(n6.eq.16) then
            NOPMT = 10929
         else if(n6.eq.17) then
            NOPMT = 10891
         else if(n6.eq.18) then
            NOPMT = 10771
         else if(n6.eq.19) then
            NOPMT = 10551
         else if(n6.eq.20) then
            NOPMT = 10273
         else if(n6.eq.21) then
            NOPMT = 10225
         else if(n6.eq.22) then
            NOPMT = 9949
         else if(n6.eq.23) then
            NOPMT = 9733
         else if(n6.eq.24) then
            NOPMT = 9617
         else if(n6.eq.25) then
            NOPMT = 9581
         else if(n6.eq.26) then
            NOPMT = 9547
         else if(n6.eq.27) then
            NOPMT = 9515
         else if(n6.eq.28) then
            NOPMT = 9485
         else if(n6.eq.29) then
            NOPMT = 9435
         else if(n6.eq.30) then
            NOPMT = 9401
         else if(n6.eq.31) then
            NOPMT = 9399
         else if(n6.eq.32) then
            NOPMT = 9400
         else if(n6.eq.33) then
            NOPMT = 9414
         else if(n6.eq.34) then
            NOPMT = 9458
         else if(n6.eq.35) then
            NOPMT = 9514
         else if(n6.eq.36) then
            NOPMT = 9546
         else if(n6.eq.37) then
            NOPMT = 9580
         else if(n6.eq.38) then
            NOPMT = 9616
         else if(n6.eq.39) then
            NOPMT = 9654
         else if(n6.eq.40) then
            NOPMT = 9774
         else if(n6.eq.41) then
            NOPMT = 9860
         else if(n6.eq.42) then
            NOPMT = 9994
         else if(n6.eq.43) then
            NOPMT = 10272
         else
            NOPMT = (n6-3001)*51 + 1
         endif
      else if(NAMES.eq.'PMP4') then
         NOPMT = (n6-4000)*51 - 1
      else if(NAMES.eq.'PMP5') then
** for edge correction on top part
         if(n6.eq.1) then
            NOPMT = 8572
         else if(n6.eq.2) then
            NOPMT = 8848
         else if(n6.eq.3) then
            NOPMT = 9064
         else if(n6.eq.4) then
            NOPMT = 9180
         else if(n6.eq.5) then
            NOPMT = 9216
         else if(n6.eq.6) then
            NOPMT = 9250
         else if(n6.eq.7) then
            NOPMT = 9282
         else if(n6.eq.8) then
            NOPMT = 9382
         else if(n6.eq.9) then
            NOPMT = 9396
         else if(n6.eq.10) then
            NOPMT = 9398
         else if(n6.eq.11) then
            NOPMT = 9397
         else if(n6.eq.12) then
            NOPMT = 9383
         else if(n6.eq.13) then
            NOPMT = 9363
         else if(n6.eq.14) then
            NOPMT = 9251
         else if(n6.eq.15) then
            NOPMT = 9217
         else if(n6.eq.16) then
            NOPMT = 9181
         else if(n6.eq.17) then
            NOPMT = 9143
         else if(n6.eq.18) then
            NOPMT = 9023
         else if(n6.eq.19) then
            NOPMT = 8803
         else if(n6.eq.20) then
            NOPMT = 8525
         else if(n6.eq.21) then
            NOPMT = 8477
         else if(n6.eq.22) then
            NOPMT = 8201
         else if(n6.eq.23) then
            NOPMT = 8069
         else if(n6.eq.24) then
            NOPMT = 7985
         else if(n6.eq.25) then
            NOPMT = 7869
         else if(n6.eq.26) then
            NOPMT = 7833
         else if(n6.eq.27) then
            NOPMT = 7799
         else if(n6.eq.28) then
            NOPMT = 7767
         else if(n6.eq.29) then
            NOPMT = 7737
         else if(n6.eq.30) then
            NOPMT = 7687
         else if(n6.eq.31) then
            NOPMT = 7653
         else if(n6.eq.32) then
            NOPMT = 7651
         else if(n6.eq.33) then
            NOPMT = 7652
         else if(n6.eq.34) then
            NOPMT = 7666
         else if(n6.eq.35) then
            NOPMT = 7710
         else if(n6.eq.36) then
            NOPMT = 7766
         else if(n6.eq.37) then
            NOPMT = 7798
         else if(n6.eq.38) then
            NOPMT = 7832
         else if(n6.eq.39) then
            NOPMT = 7868
         else if(n6.eq.40) then
            NOPMT = 7906
         else if(n6.eq.41) then
            NOPMT = 8026
         else if(n6.eq.42) then
            NOPMT = 8112
         else if(n6.eq.43) then
            NOPMT = 8246
         else if(n6.eq.44) then
            NOPMT = 8524
         else
            NOPMT = (n6-5000)*51
         endif
** for top
      else if(NAMES.eq.'PMPT') then
         NOPMT = ISEQT(N7,N8)
** for bottom
      else if(NAMES.eq.'PMPB') then
         NOPMT = ISEQB(N7,N8)
      else
         print*,'Error in sgpmtnum ',names
         stop
      endif
      iflag = 1
      IF(NOPMT.EQ.0 ! PMTs not included in PMPT ,but PBXT
     &     .OR.NOPMT.EQ.7667.OR.NOPMT.EQ.7686.OR.NOPMT.EQ.8937 !not installed 16PMTs(Top)
     &     .OR.NOPMT.EQ.8980.OR.NOPMT.EQ.9283.OR.NOPMT.EQ.9312
     &     .OR.NOPMT.EQ.9339.OR.NOPMT.EQ.9362
     &     .OR.NOPMT.EQ.9415.OR.NOPMT.EQ.9434.OR.NOPMT.EQ.9817 !not installed 17PMTs(Bot)
     &     .OR.NOPMT.EQ.10685.OR.NOPMT.EQ.10728.OR.NOPMT.EQ.11031
     &     .OR.NOPMT.EQ.11060.OR.NOPMT.EQ.11087.OR.NOPMT.EQ.11110) THEN
         IFLAG=0
      ENDIF
      RETURN
      END
