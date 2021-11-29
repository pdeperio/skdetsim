*********************************************************************
*
*     -----------------
      SUBROUTINE SGGEOM_SK2
*     -----------------
*
*     (Purpose)
*       Define the detector geometry and tracking media for SK-II
*
*     (Input)
*       none
*
*     (Output)
*       none
*
*     (Creation Date and Author)
*       1992.04.24 ; First version by J.Kanzaki
*       1993.02.01 ; Modified by Y.Koshio
*            02.13 ; Change the fastest version (rotation)
*       1997.02.04 ; Modified to correct the dead space by R.Doyle
*
************************************************************************

C     COMMONs for geometries

#include "skparm.h"
#include "geotnk.h"
#include "geopmt.h"
#include "geopma.h"

C     COMMONs for geometries (sksim local)

#include "./guskge.h"

C     COMMONs for materials

#include "./gumate.h"
#ifdef UPMU
#include "./gumed_u.h"
#else
#include "./gumed.h"
#endif
      REAL*4 PAR(3),PRTANK(3),PRINVO(3),PRFIDU(3)
c      REAL*4 PRPMTW(3),PRBOX(3),PRPMTT(5),PRPMTB(5)
      REAL*4 PRPMTW(5),PRBOX(3)
      REAL*4 PRPMTU(6),PRPMTL(6)
      REAL*4 PRXMED(3),PRWMED(3),PRXANT(3)
      REAL*4 PRACR(6)
      real ppckov(2),absco(2),effic(2),rindex(2)

      integer iyadd(48)
      real zeroz, zerot
      parameter(zeroz = 3.7, zerot = 2./(DIINTK*3.1416) *360.)

#ifdef NICKEL
c     COMMONs for source position

#include "nicpara.h"

      real*4 parsor(3)
      dimension asor(3),zsor(3),wmsor(3)
#endif

*
*     Set black sheet position:
*     Note: this is a common block parameter also used by tracking.
      ZBLST = ZBLST_SK_II
*
*     Set geometry parameters: 
      CALL GEOSET
      CALL GXYZSET
*     Stores the constants for the material in the data structure JMATE

c     CALL GSMATE(IMTWTR,'WATER$',
c    &            AWTR,ZWTR,DNWTR,RLWTR,ALWTR,UBWTR,NBWTR)

*    water
      AWTR(1)  =  AWTR1
      AWTR(2)  =  AWTR2
      ZWTR(1)  =  ZWTR1
      ZWTR(2)  =  ZWTR2
      WMWTR(1) =  WMWTR1
      WMWTR(2) =  WMWTR2
      CALL GSMIXT(IMTWTR,'WATER$',
     &            AWTR,ZWTR,DNWTR,NMWTR,WMWTR)

*    acryl

      AACR(1)  =  AACR1
      AACR(2)  =  AACR2
      AACR(3)  =  AACR3
      ZACR(1)  =  ZACR1
      ZACR(2)  =  ZACR2
      ZACR(3)  =  ZACR3
      WMACR(1) =  WMACR1
      WMACR(2) =  WMACR2
      WMACR(3) =  WMACR3
      CALL GSMIXT(IMTACR,'ACRYL$',
     &     AACR,ZACR,DNACR,NMACR,WMACR)

#ifdef NICKEL
c        material for Ni source
      asor(1) = awtr1
      asor(2) = awtr2
      asor(3) = 58.7
      zsor(1) = zwtr1
      zsor(2) = zwtr2
      zsor(3) = 28.
      wmsor(1)= 0.056
      wmsor(2)= 0.500
      wmsor(3)= 0.444
      dens    = 1.65
      call gsmixt(99,'SOURCE',asor,zsor,dens,3,wmsor)
      call gstmed(50,'Ni source',99,0,0,0.,1.,0.5,
     &     0.05,0.5,0.5,ub,nwb)
      call gstmed(51,'Ni source',99,0,0,0.,1.,0.5,
     &     0.05,0.5,0.5,ub,nwb)
      call gstmed(52,'air',15,0,0,0.,1.,0.5,
     &     0.05,0.5,0.5,ub,nwb)
#endif

C     CALL GPMATE(0)

*     Stores the parameters of the tracking medium in the data structure
*     JTMED

*     Outer water volume

      CALL GSTMED(IMDWTO,'Outer Water volume$',IMTWTR,
     &            ISWTO,IFWTO,FLDWTO,TMXWTO,DMXWTO,DEXWTO,EPSWTO,STMWTO,
     &            UBWTO,NBWTO)

*     following value is DUMMY
      ppckov(1) = 1.768E-09
      ppckov(2) = 4.126E-09
      absco(1)  = 6000.
      absco(2)  = 6000.
      effic(1)  = 0.05
      effic(2)  = 0.05
      rindex(1) = 1.344
      rindex(2) = 1.344
**
      call gsckov(imdwto,2,ppckov,absco,effic,rindex)

*     Inner water volume

      CALL GSTMED(IMDWTI,'Inner Water volume$',IMTWTR,
     &            ISWTI,IFWTI,FLDWTI,TMXWTI,DMXWTI,DEXWTI,EPSWTI,STMWTI,
     &            UBWTI,NBWTI)

      call gsckov(imdwti,2,ppckov,absco,effic,rindex)

*     PMT cylinder volume

      CALL GSTMED(IMDPMC,'PMT Cylinder volume$',IMTVAC,
     &            ISPMC,IFPMC,FLDPMC,TMXPMC,DMXPMC,DEXPMC,EPSPMC,STMPMC,
     &            UBPMC,NBPMC)

*     The place between inner and anti volume (not sensitive)

      CALL GSTMED(IMDWTM,'mediam volume$',IMTWTR,
     &            ISWTM,IFWTM,FLDWTM,TMXWTM,DMXWTM,DEXWTM,EPSWTM,STMWTM,
     &            UBWTM,NBWTM)

*     Anti PMT and wave length shifter

      CALL GSTMED(IMDPMA,'anti PMA & W.S.$',IMTVAC,
     &            ISPMA,IFPMA,FLDPMA,TMXPMA,DMXPMA,DEXPMA,EPSPMA,STMPMA,
     &            UBPMA,NBPMA)

*     Out of the tank
#ifdef UPMU
*    Add rock components 2007 4/20 T.Tanaka
*--- Modified by Yoshida
*    Define Rock out of the tank
      ARCK(1)  =  ARCK1
      ARCK(2)  =  ARCK2
      ZRCK(1)  =  ZRCK1
      ZRCK(2)  =  ZRCK2
      WMRCK(1) =  WMRCK1
      WMRCK(2) =  WMRCK2
      CALL GSMIXT(IMTRCK,'ROCK$',
     &            ARCK,ZRCK,DNRCK,NMRCK,WMRCK)
      CALL GSTMED(IMDMST,'Out of the tank$',IMTRCK,

     &            ISMST,IFMST,FLDMST,TMXMST,DMXMST,DEXMST,EPSMST,STMMST,
     &            UBMST,NBMST)
*---> add by saji
      CALL GSTPAR(IMDMST,'LOSS',2.)
*      CALL GSTPAR(IMDMST,'DRAY',0.)
#else

      CALL GSTMED(IMDMST,'Out of the tank$',IMTVAC,
     &            ISMST,IFMST,FLDMST,TMXMST,DMXMST,DEXMST,EPSMST,STMMST,
     &            UBMST,NBMST)
#endif
*     Acrylic case

      CALL GSTMED(IMDACR,'Acrylic case$',IMTACR,
     &            ISACR,IFACR,FLDACR,TMXACR,DMXACR,DEXACR,EPSACR,STMACR,
     &            UBACR,NBACR)

      ppckov(1) = 1.768E-09
      ppckov(2) = 4.126E-09
      absco(1)  = 6000.
      absco(2)  = 6000.
      effic(1)  = 0.05
      effic(2)  = 0.05
      rindex(1) = 1.49
      rindex(2) = 1.49

      call gsckov(IMDACR,2,ppckov,absco,effic,rindex)

*     Water Between PMT and Acrylic case

      CALL GSTMED(IMDWTS,'Water between PMT and Acrylic case$',IMTWTR,
     &            ISWTS,IFWTS,FLDWTS,TMXWTS,DMXWTS,DEXWTS,EPSWTS,STMWTS,
     &            UBWTS,NBWTS)

      call gsckov(imdWTS,2,ppckov,absco,effic,rindex)

c     CALL GPTMED(0)

*     Define detector volumes

*     First define an imaginary box which includes the whole detector
*     and defines MAster Referecne System (MARS).
#ifdef UPMU
*----> saji modified 8km rock box
      PAR(1)= 800000.
      PAR(2)= 800000.
      PAR(3)= 800000.
#else
      PAR(1)= DXMAST
      PAR(2)= DYMAST
      PAR(3)= DZMAST
#endif
      CALL GSVOLU('MAST','BOX ', IMDMST, PAR, 3, IVOLU)

*     Define a rotation matrix which is not rotated at all.

      IRMAST = 1
      CALL GSROTM(IRMAST, 90., 0., 90., 90., 0., 0.)

*     Then define the tank in this box

      NPTANK    = 3
      PRTANK(1) = 0.
      PRTANK(2) = RTKTK
      PRTANK(3) = 0.5*HITKTK
      CALL GSVOLU('TANK','TUBE', IMDWTO, PRTANK, NPTANK, IVTANK)
      XTANK = 0.
      YTANK = 0.
      ZTANK = 0.5*HITKTK - ZCNTTK
      IRTANK = IRMAST
      CALL GSPOS('TANK',1,'MAST',XTANK,YTANK,ZTANK,IRTANK,'ONLY')

*** MEDIUM
*     the place between inner and anti volume (not sensitive)

*     Barrel
      NPWMED = 3
      PRWMED(1) = RINTK-ZBLST
      PRWMED(2) = RINTK+RMED
******Modified to shorten the dead space which extends too high*
      PRWMED(3) = ZPINTK+ZMED
****************************************************************
      CALL GSVOLU('WMED','TUBE', IMDWTM, PRWMED, NPWMED, IVWMED)

      XWMED = 0.
      YWMED = 0.
      ZWMED = 0.
      IRWMED = IRMAST
      CALL GSPOS('WMED',1,'TANK',XWMED,YWMED,ZWMED,IRWMED,'ONLY')

*     Top & Bottom
      NPXMED = 3
      PRXMED(1) = 0.
      PRXMED(2) = RINTK-ZBLST
      PRXMED(3) = ZPINTK+ZMED
      CALL GSVOLU('XMED','TUBE', IMDWTM, PRXMED, NPXMED, IVXMED)

      XXMED = 0.
      YXMED = 0.
      ZXMED = 0.
      IRXMED = IRMAST
      CALL GSPOS('XMED',1,'TANK',XXMED,YXMED,ZXMED,IRXMED,'ONLY')

*** ANTI

*     Barrel

      NPXANT = 3
      PRXANT(1) = RINTK + RMED - ZPMCYA
      PRXANT(2) = RINTK + RMED
      PRXANT(3) = 0.5*HIWAL
      CALL GSVOLU('WANT','TUBE', IMDPMA, PRXANT, NPXANT, IVXANT)

      XXANT = 0.
      YXANT = 0.
      ZXANT = 0.
      IRXANT = 2
      PHII   = 360./FLOAT(MXWTPM)/2.
      CALL GSROTM(IRXANT, 90.,-PHII, 90., 90.-PHII, 0., 0.)
      CALL GSPOS('WANT',1,'WMED',XXANT,YXANT,ZXANT,IRXANT,'ONLY')

      CALL GSDVN('PMAZ','WANT',MXWZPMA,3)
      CALL GSDVN('PMAW','PMAZ',MXWTPMA,2)

*     Top & Bottom

*     plate
      NPXANT = 3
      PRXANT(1) = 0.
      PRXANT(2) = RINTK
      PRXANT(3) = ZPMCYA/2.
      CALL GSVOLU('ANTB','TUBE', IMDPMA, PRXANT, NPXANT, IVXANT)

      XXANT = 0.
      YXANT = 0.
      ZXANT = ZPINTK + ZMED
      IRXANT = IRMAST
*     for top
      CALL GSPOS('ANTB',1,'XMED',XXANT,YXANT,ZXANT,IRXANT,'ONLY')
*     for bottom
      ZXANT =-ZXANT
      CALL GSPOS('ANTB',2,'XMED',XXANT,YXANT,ZXANT,IRXANT,'ONLY')

*       center
*     anti PMT volume
      NPXANT = 3
      PRXANT(1) = 1272.6
      PRXANT(2) = 1131.2
      PRXANT(3) = ZPMCYA/2.
      CALL GSVOLU('ANTC','BOX ', IMDPMA, PRXANT, NPXANT, IVXANT)

      XXANT = 0.
      YXANT = 0.
      ZXANT = 0.
      IRXANT = IRMAST
      CALL GSPOS('ANTC',1,'ANTB',XXANT,YXANT,ZXANT,IRXANT,'ONLY')

      CALL GSDVN('ANDV','ANTC',12,1)
      CALL GSDVN('ANPM','ANDV',16,2)

*       side 1
      NPXANT = 3
      PRXANT(1) = 106.05
C      PRXANT(2) = 848.4
      PRXANT(2) = 989.8
      PRXANT(3) = ZPMCYA/2.
      CALL GSVOLU('ANT1','BOX ', IMDPMA, PRXANT, NPXANT, IVXANT)

      XXANT = 1272.6 + 106.05
      YXANT = 0.
      ZXANT = 0.
      IRXANT = IRMAST
      CALL GSPOS('ANT1',1,'ANTB',XXANT,YXANT,ZXANT,IRXANT,'ONLY')
      XXANT = -XXANT
      CALL GSPOS('ANT1',2,'ANTB',XXANT,YXANT,ZXANT,IRXANT,'ONLY')

C      CALL GSDVN('ANP1','ANT1',12,2)
      CALL GSDVN('ANP1','ANT1',14,2)

*       side 2
      NPXANT = 3
C      PRXANT(1) = 1060.5
      PRXANT(1) = 1272.6
      PRXANT(2) = 141.4
      PRXANT(3) = ZPMCYA/2.
      CALL GSVOLU('ANT2','BOX ', IMDPMA, PRXANT, NPXANT, IVXANT)

      XXANT = 0.
      YXANT = 1131.2 + 141.4
      ZXANT = 0.
      IRXANT = IRMAST
      CALL GSPOS('ANT2',1,'ANTB',XXANT,YXANT,ZXANT,IRXANT,'ONLY')
      YXANT = -YXANT
      CALL GSPOS('ANT2',2,'ANTB',XXANT,YXANT,ZXANT,IRXANT,'ONLY')

C      CALL GSDVN('AND2','ANT2',10,1)
      CALL GSDVN('AND2','ANT2',12,1)
      CALL GSDVN('ANP2','AND2',2,2)

*       side 3
      NPXANT = 3
      PRXANT(1) = 106.05
C      PRXANT(2) = 565.6
      PRXANT(2) = 707.9
      PRXANT(3) = ZPMCYA/2.
      CALL GSVOLU('ANT3','BOX ', IMDPMA, PRXANT, NPXANT, IVXANT)

      XXANT = 1272.6 + 212.1 + 106.05
      YXANT = 0.
      ZXANT = 0.
      IRXANT = IRMAST
      CALL GSPOS('ANT3',1,'ANTB',XXANT,YXANT,ZXANT,IRXANT,'ONLY')
      XXANT = -XXANT
      CALL GSPOS('ANT3',2,'ANTB',XXANT,YXANT,ZXANT,IRXANT,'ONLY')

C      CALL GSDVN('ANP3','ANT3',8,2)
      CALL GSDVN('ANP3','ANT3',10,2)

*       side 4
      NPXANT = 3
C      PRXANT(1) = 636.3
      PRXANT(1) = 848.8
C      PRXANT(2) = 141.4
      PRXANT(2) = 110.0
      PRXANT(3) = ZPMCYA/2.
      CALL GSVOLU('ANT4','BOX ', IMDPMA, PRXANT, NPXANT, IVXANT)

      XXANT = 0.
C      YXANT = 1131.2 + 282.8 + 141.4
      YXANT = 1131.2 + 282.8 + 110.0
      ZXANT = 0.
      IRXANT = IRMAST
      CALL GSPOS('ANT4',1,'ANTB',XXANT,YXANT,ZXANT,IRXANT,'ONLY')
      YXANT = -YXANT
      CALL GSPOS('ANT4',2,'ANTB',XXANT,YXANT,ZXANT,IRXANT,'ONLY')

C      CALL GSDVN('AND4','ANT4',6,1)
      CALL GSDVN('AND4','ANT4',8,1)
      CALL GSDVN('ANP4','AND4',2,2)

*** INNER

*     Define the inner volume in the tank.

      NPINVO = 3
      PRINVO(1) = 0.
      PRINVO(2) = RINTK-ZBLST
      PRINVO(3) = 0.5*HIINTK-ZBLST
      CALL GSVOLU('INVO','TUBE', IMDWTI, PRINVO, NPINVO, IVINVO)
      XINVO = 0.
      YINVO = 0.
C     ZINVO = 0.
      ZINVO = -ZTANK
      IRINVO = IRMAST
      CALL GSPOS('INVO',1,'XMED',XINVO,YINVO,ZINVO,IRINVO,'ONLY')

*     Define the pseudo-fiducial volume in the tank.
c
c      NPFIDU = 3
c      PRFIDU(1) = 0.
cc      PRFIDU(2) = RINTK - 2.*ZPMCY
cc      PRFIDU(3) = 0.5*HIINTK - 2.*ZPMCY
c      PRFIDU(2) = RINTK - ZACCY - WIDACR/2.0
c      PRFIDU(3) = 0.5*HIINTK - ZACCY - WIDACR/2.0
c      CALL GSVOLU('FIDU','TUBE', IMDWTI, PRFIDU, NPFIDU, IVFIDU)
c      XFIDU = XINVO
c      YFIDU = YINVO
c      ZFIDU = ZINVO
c      IRFIDU = IRMAST
c      CALL GSPOS('FIDU',1,'INVO',XFIDU,YFIDU,ZFIDU,IRFIDU,'ONLY')

#ifdef LINAC
      call lingeom1
#endif

#ifdef NICKEL
      parsor(1)=0.
      parsor(2)=9.
      parsor(3)=2.5
      call gsvolu('sorc','TUBE',50,parsor,3,ivolu)
      call gspos('sorc',1,'INVO',sorx,sory,sorz-6.,irmast,'only')
      parsor(1)=2.5
      parsor(2)=9.
      parsor(3)=6.
      call gsvolu('sor2','TUBE',51,parsor,3,ivolu)
      call gspos('sor2',1,'INVO',sorx,sory,sorz+2.5,irmast,'only')
      parsor(1)=0.
      parsor(2)=2.5
      parsor(3)=6.
      call gsvolu('airr','TUBE',52,parsor,3,ivolu)
      call gspos('airr',1,'INVO',sorx,sory,sorz+2.5,irmast,'only')
#endif

************
*     Define the PMT volume on the wall
************

c      NPPMTW = 3
c      PRPMTW(1) = RINTK - ZPMCY*2.
c      PRPMTW(2) = RINTK
c      PRPMTW(3) = 0.5*HIWAL
c      CALL GSVOLU('PMPW','TUBE', IMDWTI, PRPMTW, NPPMTW, IVPMTW)
c
c      XPMTW = 0.
c      YPMTW = 0.
cC     ZPMTW = 0.
c      ZPMTW = -ZTANK
c      IRPMTW = 2
c      CALL GSPOS('PMPW',1,'INVO',XPMTW,YPMTW,ZPMTW,IRPMTW,'ONLY')
c
c*     Divide this volume into cells
c
c      CALL GSDVN('PMCZ','PMPW',MXWZPM,3)
cC     CALL GSORD('PMCZ',3)
c
c      CALL GSDVN('PMCW','PMCZ',MXWTPM,2)
cC     CALL GSORD('PMCW',6)

*     Define a water hemisphere between acrylic case and PMT surface
*     medium number is IMDWTS ( = 8)

      NPACR = 6
      PRACR(1) = 0.
      PRACR(2) = RACR+WIDACR/2.
      PRACR(3) = 0.
      PRACR(4) = 90.
      PRACR(5) = 0.
      PRACR(6) = 360.
      CALL GSVOLU('SPCW','SPHE',IMDWTS,PRACR,NPACR,IVACR)

      NPACR = 6
      PRACR(1) = 0
      PRACR(2) = RACR+WIDACR/2.
      PRACR(3) = 90.
      PRACR(4) = 180.
      PRACR(5) = 0.
      PRACR(6) = 360.
      CALL GSVOLU('SPCT','SPHE',IMDWTS,PRACR,NPACR,IVACR)

      PRACR(1) = 0
      PRACR(2) = RACR+WIDACR/2.
      PRACR(3) = 0.
      PRACR(4) = 90.
      PRACR(5) = 0.
      PRACR(6) = 360.
      CALL GSVOLU('SPCB','SPHE',IMDWTS,PRACR,NPACR,IVACR)

* main position

      NPPMTW = 5
      PRPMTW(1) = RINTK - ZACCY - WIDACR/2.0
      PRPMTW(2) = RINTK - ZBLST
      PRPMTW(3) = 0.5*(HIWAL-212.1)
      PRPMTW(4) = 9.6
      PRPMTW(5) = 9.6 + 259.2   ! = offset + 360.*(108/150)
      CALL GSVOLU('PMPW','TUBS', IMDWTI, PRPMTW, NPPMTW, IVPMTW)

      XPMTW = 0.
      YPMTW = 0.
c      ZPMTW = -ZTANK
      ZPMTW = -35.35
      IRPMTW = 2
      CALL GSPOS('PMPW',1,'INVO',XPMTW,YPMTW,ZPMTW,IRPMTW,'ONLY')

      CALL GSDVN('PMCZ','PMPW',8,3)
      CALL GSDVN('PMCW','PMCZ',3,2)

      do i=1,108
         if(i.ne.7.and.i.ne.26.and.i.ne.43.and.
     &        i.ne.63.and.i.ne.80.and.i.ne.99) then
            if(mod(i,3).eq.1) then
               if(mod((i-1)/3,2).eq.0) then
                  ZPMTU = -106.05 - zeroz
               else
                  ZPMTU = -176.75
               endif
            else if(mod(i,3).eq.2) then
               if(mod((i-1)/3,2).eq.0) then
                  ZPMTU = 35.35
               else
                  ZPMTU = -35.35 + zeroz
               endif
            else if(mod(i,3).eq.0) then
               if(mod((i-1)/3,2).eq.0) then
                  ZPMTU = 176.75 + zeroz
               else
                  ZPMTU = 106.05 - zeroz
               endif
            endif

            degtub = float((i-1)/3)*2.4-42. ! 2.4=360./150
            if(mod(i,12).eq.1.or.mod(i,12).eq.2.or.mod(i,12).eq.3) then
               degtub = degtub + zerot
            else if(mod(i,12).eq.10.or.mod(i,12).eq.11.or.mod(i,12).eq.0) then
               degtub = degtub - zerot
            endif
            radtub = degtub*3.1416/180.

            IRPMT1 = 500+i
            rottubx = 90.+degtub
            rottubz = 180.+degtub
            CALL GSROTM(IRPMT1, 90., rottubx, 180., 0., 90., rottubz)

            XACRW = (RINTK + (RACR - ZACCY))*cos(radtub)
            YACRW = (RINTK + (RACR - ZACCY))*sin(radtub)
            ZACRW = ZPMTU
            IRACRW= IRPMT1
            CALL GSPOS('SPCW',i,'PMCW',XACRW,YACRW,ZACRW,IRACRW,'MANY')

         endif
      enddo

* extra position 1 (29th SM)

      NPPMTW = 5
      PRPMTW(1) = RINTK - ZACCY - WIDACR/2.0
      PRPMTW(2) = RINTK - ZBLST
      PRPMTW(3) = 0.5*(HIWAL-212.1)
      PRPMTW(4) = 268.8         ! 360.*(112/150)
      PRPMTW(5) = 268.8 + 4.8   ! + 360.*(2/150)
      CALL GSVOLU('PMP1','TUBS', IMDWTI, PRPMTW, NPPMTW, IVPMTW)

      XPMTW = 0.
      YPMTW = 0.
c      ZPMTW = -ZTANK
      ZPMTW = -35.35
      IRPMTW = 2
      CALL GSPOS('PMP1',1,'INVO',XPMTW,YPMTW,ZPMTW,IRPMTW,'ONLY')

      CALL GSDVN('PMC1','PMP1',8,3)

      do i=1,6
         if(i.eq.1) then
            ZPMTU = -106.05 - zeroz
            degtub = 270.-360. + zerot
         else if(i.eq.2) then
            ZPMTU = 35.35
            degtub = 270.-360. + zerot
         else if(i.eq.3) then
            ZPMTU = 176.75 + zeroz
            degtub = 270.-360. + zerot
         else if(i.eq.4) then
            ZPMTU = -176.75
            degtub = 272.4-360.-zerot
         else if(i.eq.5) then
            ZPMTU = -35.35 + zeroz
            degtub = 272.4-360.-zerot
         else if(i.eq.6) then
            ZPMTU = 106.05 - zeroz
            degtub = 272.4-360.-zerot
         endif
         radtub = degtub*3.1416/180.

         inum = i+1000
         IRPMT1 = inum
         rottubx = 90.+degtub
         rottubz = 180.+degtub
         CALL GSROTM(IRPMT1, 90., rottubx, 180., 0., 90., rottubz)

         XACRW = (RINTK + (RACR - ZACCY))*cos(radtub)
         YACRW = (RINTK + (RACR - ZACCY))*sin(radtub)
         ZACRW = ZPMTU
         IRACRW= IRPMT1
         CALL GSPOS('SPCW',inum,'PMC1',XACRW,YACRW,ZACRW,IRACRW,'MANY')
      enddo

* extra position 2

      NPPMTW = 5
      PRPMTW(1) = RINTK - ZACCY - WIDACR/2.0
      PRPMTW(2) = RINTK - ZBLST
      PRPMTW(3) = 0.5*(HIWAL-212.1)
      PRPMTW(4) = 273.6         ! 360.*(114/150) 
      PRPMTW(5) = 360. + 9.6    ! + 360.*(4/150) 
      CALL GSVOLU('PMP2','TUBS', IMDWTI, PRPMTW, NPPMTW, IVPMTW)

      XPMTW = 0.
      YPMTW = 0.
      ZPMTW = -35.35
      IRPMTW = 2
      CALL GSPOS('PMP2',1,'INVO',XPMTW,YPMTW,ZPMTW,IRPMTW,'ONLY')

      CALL GSDVN('PMC2','PMP2',8,3)

      do i=1,120
         if(i.ne.1.and.i.ne.20.and.i.ne.37.and.i.ne.57.and.
     &        i.ne.74.and.i.ne.93.and.i.ne.111) then

            if(mod(i,3).eq.1) then
               if(mod((i-1)/3,2).eq.0) then
                  ZPMTU = -106.05 - zeroz
               else
                  ZPMTU = -176.75
               endif
            else if(mod(i,3).eq.2) then
               if(mod((i-1)/3,2).eq.0) then
                  ZPMTU = 35.35
               else
                  ZPMTU = -35.35 + zeroz
               endif
            else if(mod(i,3).eq.0) then
               if(mod((i-1)/3,2).eq.0) then
                  ZPMTU = 176.75 + zeroz
               else
                  ZPMTU = 106.05 - zeroz
               endif
            endif

            degtub = 274.8 + float((i-1)/3)*2.4 - 360.
            if(mod(i,12).eq.1.or.mod(i,12).eq.2.or.mod(i,12).eq.3) then
               degtub = degtub + zerot
            else if(mod(i,12).eq.10.or.mod(i,12).eq.11.or.mod(i,12).eq.0) then
               degtub = degtub - zerot
            endif
            radtub = degtub*3.1416/180.

            inum=i+2000
            IRPMT1 = inum
            rottubx = 90.+degtub
            rottubz = 180.+degtub
            CALL GSROTM(IRPMT1, 90., rottubx, 180., 0., 90., rottubz)

            XACRW = (RINTK + (RACR - ZACCY))*cos(radtub)
            YACRW = (RINTK + (RACR - ZACCY))*sin(radtub)
            ZACRW = ZPMTU
            IRACRW= IRPMT1
            CALL GSPOS('SPCW',inum,'PMC2',XACRW,YACRW,ZACRW,IRACRW,'MANY')

         endif
      enddo

* extra position 3

      NPPMTW = 3
      PRPMTW(1) = RINTK - ZACCY - WIDACR/2.0
      PRPMTW(2) = RINTK - ZBLST
      PRPMTW(3) = 0.5*70.7
      CALL GSVOLU('PMP3','TUBE', IMDWTI, PRPMTW, NPPMTW, IVPMTW)

      XPMTW = 0.
      YPMTW = 0.
c      ZPMTW = -ZTANK
      ZPMTW = -1767.5
      IRPMTW = 1
      CALL GSPOS('PMP3',1,'INVO',XPMTW,YPMTW,ZPMTW,IRPMTW,'ONLY')

      do i=1,150
         if(mod(i,2).ne.0) then
            if(i.ne.1.and.i.ne.13.and.i.ne.25.and.i.ne.37.and.
     &           i.ne.61.and.i.ne.73.and.i.ne.97.and.i.ne.109.and.
     &           i.ne.113.and.i.ne.133.and.i.ne.135.and.i.ne.145) then

               ZPMTU = 0. + zeroz
               degtub = float(i-1)*2.4
               if(i.le.112) then
                  if(mod(i,4).eq.1) then
                     degtub = degtub + zerot
                  else if(mod(i,4).eq.0) then
                     degtub = degtub - zerot
                  endif
               else if(i.eq.113) then ! 29th
                  degtub = degtub + zerot
               else if(i.eq.114) then ! 29th
                  degtub = degtub - zerot
               else
                  if(mod(i,4).eq.3) then
                     degtub = degtub + zerot
                  else if(mod(i,4).eq.2) then
                     degtub = degtub - zerot
                  endif
               endif
               radtub = degtub*3.1416/180.

               inum=i+3000
               IRPMT1 = inum
               rottubx = 90.+degtub
               rottubz = 180.+degtub
               CALL GSROTM(IRPMT1, 90., rottubx, 180., 0., 90., rottubz)

               XACRW = (RINTK + (RACR - ZACCY))*cos(radtub)
               YACRW = (RINTK + (RACR - ZACCY))*sin(radtub)
               ZACRW = ZPMTU
               CALL GSPOS('SPCW',inum,'PMP3',XACRW,YACRW,ZACRW,IRPMT1,'MANY')
            endif
         endif
      enddo

*  edge correction
      zacrb = -52.9
*    #9580
      xacrb = 1166.55
      yacrb = -1166.55
      CALL GSPOS('SPCB',1,'PMP3',XACRB,YACRB,ZACRB,0,'MANY')
*    #9733
      xacrb = -1449.35
      yacrb = -813.05
      CALL GSPOS('SPCB',2,'PMP3',XACRB,YACRB,ZACRB,0,'MANY')
*    #9994
      xacrb = 1590.75
      yacrb = -459.55
      CALL GSPOS('SPCB',3,'PMP3',XACRB,YACRB,ZACRB,0,'MANY')
*    #10225
      xacrb = -1661.45
      yacrb = -35.35
      CALL GSPOS('SPCB',4,'PMP3',XACRB,YACRB,ZACRB,0,'MANY')
*    #10812
      xacrb = 1449.35
      yacrb = 813.05
      CALL GSPOS('SPCB',5,'PMP3',XACRB,YACRB,ZACRB,0,'MANY')

* extra position 4

      NPPMTW = 3
      PRPMTW(1) = RINTK - ZACCY - WIDACR/2.0
      PRPMTW(2) = RINTK - ZBLST
      PRPMTW(3) = 0.5*70.7
      CALL GSVOLU('PMP4','TUBE', IMDWTI, PRPMTW, NPPMTW, IVPMTW)

      XPMTW = 0.
      YPMTW = 0.
c      ZPMTW = -ZTANK
      ZPMTW = 1696.8
      IRPMTW = 1
      CALL GSPOS('PMP4',1,'INVO',XPMTW,YPMTW,ZPMTW,IRPMTW,'ONLY')

      do i=1,150
         if(mod(i,2).eq.0) then
            ZPMTU = 0.
            degtub = float(i-1)*2.4
            if(i.le.112) then
               if(mod(i,4).eq.1) then
                  degtub = degtub + zerot
               else if(mod(i,4).eq.0) then
                  degtub = degtub - zerot
               endif
            else if(i.eq.113) then ! 29th
               degtub = degtub + zerot
            else if(i.eq.114) then ! 29th
               degtub = degtub - zerot
            else
               if(mod(i,4).eq.3) then
                  degtub = degtub + zerot
               else if(mod(i,4).eq.2) then
                  degtub = degtub - zerot
               endif
            endif
            radtub = degtub*3.1416/180.

            inum=i+4000
            IRPMT1 = inum
            rottubx = 90.+degtub
            rottubz = 180.+degtub
            CALL GSROTM(IRPMT1, 90., rottubx, 180., 0., 90., rottubz)

            XACRW = (RINTK + (RACR - ZACCY))*cos(radtub)
            YACRW = (RINTK + (RACR - ZACCY))*sin(radtub)
            ZACRW = ZPMTU
            CALL GSPOS('SPCW',inum,'PMP4',XACRW,YACRW,ZACRW,IRPMT1,'MANY')
         endif
      enddo

* extra position 5

      NPPMTW = 3
      PRPMTW(1) = RINTK - ZACCY - WIDACR/2.0
      PRPMTW(2) = RINTK - ZBLST
      PRPMTW(3) = 0.5*70.7
      CALL GSVOLU('PMP5','TUBE', IMDWTI, PRPMTW, NPPMTW, IVPMTW)

      XPMTW = 0.
      YPMTW = 0.
c      ZPMTW = -ZTANK
      ZPMTW = 1767.5
      IRPMTW = 1
      CALL GSPOS('PMP5',1,'INVO',XPMTW,YPMTW,ZPMTW,IRPMTW,'ONLY')

      do i=1,150
         if(mod(i,2).ne.0) then
            if(i.ne.7.and.i.ne.17.and.i.ne.19.and.i.ne.39.and.
     &           i.ne.43.and.i.ne.55.and.i.ne.79.and.i.ne.91.and.
     &           i.ne.115.and.i.ne.127.and.i.ne.139) then

               ZPMTU = 0. - zeroz
               degtub = float(i-1)*2.4
               if(i.le.112) then
                  if(mod(i,4).eq.1) then
                     degtub = degtub + zerot
                  else if(mod(i,4).eq.0) then
                     degtub = degtub - zerot
                  endif
               else if(i.eq.113) then ! 29th
                  degtub = degtub + zerot
               else if(i.eq.114) then ! 29th
                  degtub = degtub - zerot
               else
                  if(mod(i,4).eq.3) then
                     degtub = degtub + zerot
                  else if(mod(i,4).eq.2) then
                     degtub = degtub - zerot
                  endif
               endif
               radtub = degtub*3.1416/180.

               inum=i+5000
               IRPMT1 = inum
               rottubx = 90.+degtub
               rottubz = 180.+degtub
               CALL GSROTM(IRPMT1, 90., rottubx, 180., 0., 90., rottubz)

               XACRW = (RINTK + (RACR - ZACCY))*cos(radtub)
               YACRW = (RINTK + (RACR - ZACCY))*sin(radtub)
               ZACRW = ZPMTU
               CALL GSPOS('SPCW',inum,'PMP5',XACRW,YACRW,ZACRW,IRPMT1,'MANY')
            endif
         endif
      enddo

*  edge correction
      zacrb = 52.9
*    #7652
      xacrb = 35.35
      yacrb = -1661.45
      CALL GSPOS('SPCT',1,'PMP5',XACRB,YACRB,ZACRB,0,'MANY')
*    #7710
      xacrb = 813.05
      yacrb = -1449.35
      CALL GSPOS('SPCT',2,'PMP5',XACRB,YACRB,ZACRB,0,'MANY')
*    #7737
      xacrb = -1025.15
      yacrb = -1307.95
      CALL GSPOS('SPCT',3,'PMP5',XACRB,YACRB,ZACRB,0,'MANY')
*    #7767
      xacrb = -1095.85
      yacrb = -1237.25
      CALL GSPOS('SPCT',4,'PMP5',XACRB,YACRB,ZACRB,0,'MANY')
*    #7869
      xacrb = -1307.95
      yacrb = -1025.15
      CALL GSPOS('SPCT',5,'PMP5',XACRB,YACRB,ZACRB,0,'MANY')
*    #8026
      xacrb = 1449.35
      yacrb = -813.05
      CALL GSPOS('SPCT',6,'PMP5',XACRB,YACRB,ZACRB,0,'MANY')
*    #8069
      xacrb = -1520.05
      yacrb = -671.65
      CALL GSPOS('SPCT',7,'PMP5',XACRB,YACRB,ZACRB,0,'MANY')
*    #8524
      xacrb = 1661.45
      yacrb = -35.35
      CALL GSPOS('SPCT',8,'PMP5',XACRB,YACRB,ZACRB,0,'MANY')
*    #8525
      xacrb = -1661.45
      yacrb = 35.35
      CALL GSPOS('SPCT',9,'PMP5',XACRB,YACRB,ZACRB,0,'MANY')
*    #8848
      xacrb = 1590.75
      yacrb = 459.55
      CALL GSPOS('SPCT',10,'PMP5',XACRB,YACRB,ZACRB,0,'MANY')
*    #9023
      xacrb = -1449.35
      yacrb = 813.05
      CALL GSPOS('SPCT',11,'PMP5',XACRB,YACRB,ZACRB,0,'MANY')
*    #9180
      xacrb = 1307.95
      yacrb = 1025.15
      CALL GSPOS('SPCT',12,'PMP5',XACRB,YACRB,ZACRB,0,'MANY')
*    #9216
      xacrb = 1237.25
      yacrb = 1095.85
      CALL GSPOS('SPCT',13,'PMP5',XACRB,YACRB,ZACRB,0,'MANY')
*    #9250
      xacrb = 1166.55
      yacrb = 1166.55
      CALL GSPOS('SPCT',14,'PMP5',XACRB,YACRB,ZACRB,0,'MANY')
*    #9396
      xacrb = 459.55
      yacrb = 1590.75
      CALL GSPOS('SPCT',15,'PMP5',XACRB,YACRB,ZACRB,0,'MANY')
*    #9397
      xacrb = -35.35
      yacrb = 1661.45
      CALL GSPOS('SPCT',16,'PMP5',XACRB,YACRB,ZACRB,0,'MANY')

*     Define a Acrylic case
*     medium number is IMDACR ( = 7)

      NPACR = 6
      PRACR(1) = RACR-WIDACR/2.
      PRACR(2) = RACR+WIDACR/2.
      PRACR(3) = 0.
      PRACR(4) = 90.
      PRACR(5) = 0.
      PRACR(6) = 360.
      CALL GSVOLU('ACRW','SPHE',IMDACR,PRACR,NPACR,IVACR)

      XACRW = 0.
      YACRW = 0.
      ZACRW = 0.
c      IRACRW= IRPMT1
      IRACRW= 0
      CALL GSPOS('ACRW',1,'SPCW',XACRW,YACRW,ZACRW,IRACRW,'MANY')

*     Define a sphere which is surface of PMT.
*     two kinds of SPHERE

* 1st
      NPPMTU = 6
      PRPMTU(1) = 0.
      PRPMTU(2) = RPMTU
      PRPMTU(3) = 0.
      PRPMTU(4) = 90.
      PRPMTU(5) = 0.
      PRPMTU(6) = 360.
c      CALL GSVOLU('PMW1','SPHE',IMDWTI,PRPMTU,NPPMTU,IVPMTU)
      CALL GSVOLU('PMW1','SPHE',IMDWTS,PRPMTU,NPPMTU,IVPMTU)

      XPMTU = 0.
      YPMTU = 0.
      ZPMTU = (RACR - ZACCY)-(RPMTU - ZPMCY*2.) 
      IRPMT1 = 0
      CALL GSPOS('PMW1',1,'SPCW',XPMTU,YPMTU,ZPMTU,IRPMT1,'MANY')

* 2nd
      NPPMTL = 6
      RPMTL  = SQRT((RPMT)**2+(APMTL)**2)
      PRPMTL(1) = 0.
      PRPMTL(2) = RPMTL
      PRPMTL(3) = 0.
      PRPMTL(4) = 90.
      PRPMTL(5) = 0.
      PRPMTL(6) = 360.
      CALL GSVOLU('PMTW','SPHE',IMDPMC,PRPMTL,NPPMTL,IVPMTL)

      XPMT2 = 0.
      YPMT2 = 0.
      ZPMT2 = (RPMTU - ZPMCY*2.) - APMTL
      IRPMT2= 0
      CALL GSPOS('PMTW',1,'PMW1',XPMT2,YPMT2,ZPMT2,IRPMT2,'MANY')

************
*    Define the PMT volume on the top
************

c      NPPMTT = 3
c      PRPMTT(1) = 0.
c      PRPMTT(2) = RINTK
c      PRPMTT(3) = ZPMCY
c      CALL GSVOLU('PMPT','TUBE', IMDWTI, PRPMTT, NPPMTT, IVPMTP)
c
c*     Position it at the top
c
c      XPMTT = 0.
c      YPMTT = 0.
cC     ZPMTT = 0.
c      ZPMTT = -ZTANK + 0.5*HIINTK - ZPMCY
c
c      CALL GSPOS('PMPT',1,'INVO',XPMTT,YPMTT,ZPMTT,0,'ONLY')

*   define top plate

*   Z size and position of the plate (common)
C      PRBOX(3) = ZPMCY
      PRBOX(3) = (ZACCY+WIDACR/2.0-ZBLST)/2.0
C      ZPMTT = 0.
C      ZPMTT = -ZTANK + 0.5*HIINTK - ZPMCY
      ZPMTT = -ZTANK + 0.5*HIINTK - ZBLST - (ZACCY+WIDACR/2.0-ZBLST)/2.0

      NPBOX = 3
      PRBOX(1) = 24.*STEPX/2.
      PRBOX(2) = 36.*STEPY/2.
      CALL GSVOLU('PBXT','BOX ', IMDWTI, PRBOX, NPBOX, IVPBXT)

      XPMTT = 0.
      YPMTT = 0.
      CALL GSPOS('PBXT',1,'INVO',xpmtt,ypmtt,zpmtt,0,'ONLY')

      CALL GSDVN('PXT1','PBXT',4,1)
      CALL GSDVN('PXT2','PXT1',6,2)

      NPBOX = 3
      PRBOX(1) = 6.*STEPX/2.
      PRBOX(2) = 24.*STEPY/2.
      CALL GSVOLU('PBX3','BOX ', IMDWTI, PRBOX, NPBOX, IVPBXT)

      XPMTT = -15.*STEPX
      YPMTT = 0.
      CALL GSPOS('PBX3',1,'INVO',xpmtt,ypmtt,zpmtt,0,'ONLY')

      CALL GSDVN('PXT3','PBX3',4,2)

      NPBOX = 3
      PRBOX(1) = 6.*STEPX/2.
      PRBOX(2) = 12.*STEPY/2.
      CALL GSVOLU('PBX9','BOX ', IMDWTI, PRBOX, NPBOX, IVPBXT)

      XPMTT = 15.*STEPX
      YPMTT = -6.*STEPY
      CALL GSPOS('PBX9',1,'INVO',xpmtt,ypmtt,zpmtt,0,'ONLY')

      CALL GSDVN('PXT9','PBX9',2,2)

      NPBOX = 3
      PRBOX(1) = 6.*STEPX/2.
      PRBOX(2) = 6.*STEPY/2.
      CALL GSVOLU('PXT0','BOX ', IMDWTI, PRBOX, NPBOX, IVPBXT)

      XPMTT = 15.*STEPX
      YPMTT = 9.*STEPY
      CALL GSPOS('PXT0',1,'INVO',xpmtt,ypmtt,zpmtt,0,'ONLY')

      NPBOX = 3
      PRBOX(1) = 6.*STEPX/2.
      PRBOX(2) = 6.*STEPY/2.
      CALL GSVOLU('PXTT','BOX ', IMDWTI, PRBOX, NPBOX, IVPBXT)

      XPMTT = 15.*STEPX
      YPMTT = 3.*STEPY
      CALL GSPOS('PXTT',1,'INVO',xpmtt,ypmtt,zpmtt,0,'ONLY')

      NPBOX = 3
      PRBOX(1) = 3.*STEPX/2.
      PRBOX(2) = 24.*STEPY/2.
      CALL GSVOLU('PBX4','BOX ', IMDWTI, PRBOX, NPBOX, IVPBXT)
      CALL GSVOLU('PBX5','BOX ', IMDWTI, PRBOX, NPBOX, IVPBXT)

      XPMTT = -19.5*STEPX
      YPMTT = 0.
      CALL GSPOS('PBX4',1,'INVO',xpmtt,ypmtt,zpmtt,0,'ONLY')

      CALL GSDVN('PXT4','PBX4',4,2)

      XPMTT = 19.5*STEPX
      YPMTT = 0.
      CALL GSPOS('PBX5',1,'INVO',xpmtt,ypmtt,zpmtt,0,'ONLY')

      CALL GSDVN('PXT5','PBX5',12,2)

      NPBOX = 3
      PRBOX(1) = 3.*STEPX/2.
      PRBOX(2) = 6.*STEPY/2.
      CALL GSVOLU('PBX6','BOX ', IMDWTI, PRBOX, NPBOX, IVPBXT)

      XPMTT = 13.5*STEPX
      YPMTT = 15.*STEPY
      CALL GSPOS('PBX6',1,'INVO',xpmtt,ypmtt,zpmtt,0,'ONLY')
      YPMTT = -15.*STEPY
      CALL GSPOS('PBX6',2,'INVO',xpmtt,ypmtt,zpmtt,0,'ONLY')

      CALL GSDVN('PXT6','PBX6',3,2)

      NPBOX = 3
      PRBOX(1) = 24.*STEPX/2.
      PRBOX(2) = 2.*STEPY/2.
      CALL GSVOLU('PBX7','BOX ', IMDWTI, PRBOX, NPBOX, IVPBXT)

      XPMTT = 0.
      YPMTT = 19.*STEPY
      CALL GSPOS('PBX7',1,'INVO',xpmtt,ypmtt,zpmtt,0,'ONLY')
      YPMTT = -19.*STEPY
      CALL GSPOS('PBX7',2,'INVO',xpmtt,ypmtt,zpmtt,0,'ONLY')

      CALL GSDVN('PXT7','PBX7',4,1)

      NPBOX = 3
      PRBOX(1) = STEPX/2.
      PRBOX(2) = STEPY/2.
      CALL GSVOLU('PXT8','BOX ', IMDWTI, PRBOX, NPBOX, IVPBXT)

      inum=0
      do ix=-24,-1
         XPMTT = float(ix)*STEPX + STEPX/2.
         call topbot(2,ix,num,iyadd)
         do iy=1,num
            inum=inum+1
            if(iyadd(iy).le.0) then
               YPMTT = float(iyadd(iy))*STEPY + STEPY/2.
            else
               YPMTT = float(iyadd(iy))*STEPY - STEPY/2.
            endif
            CALL GSPOS('PXT8',inum,'INVO',xpmtt,ypmtt,zpmtt,0,'ONLY')
         enddo
      enddo
      do ix=1,24
         XPMTT = float(ix)*STEPX - STEPX/2.
         call topbot(2,ix,num,iyadd)
         do iy=1,num
            inum=inum+1
            if(iyadd(iy).le.0) then
               YPMTT = float(iyadd(iy))*STEPY + STEPY/2.
            else
               YPMTT = float(iyadd(iy))*STEPY - STEPY/2.
            endif
            CALL GSPOS('PXT8',inum,'INVO',xpmtt,ypmtt,zpmtt,0,'ONLY')
         enddo
      enddo

*     Define a water hemisphere between acrylic case and PMT surface.
*     Medium number is IMDWTS ( = 8).

      ZACRT = RACR - ( (ZACCY+WIDACR/2-ZBLST)/2 - WIDACR/2 )

      do i=1,18
         if(i.ne.11) then
            if(i.eq.1.or.i.eq.7.or.i.eq.13) then
               XACRT = -176.75
            else if(i.eq.2.or.i.eq.8.or.i.eq.14) then
               XACRT = -35.35
            else if(i.eq.3.or.i.eq.9.or.i.eq.15) then
               XACRT = 106.05
            else if(i.eq.4.or.i.eq.10.or.i.eq.16) then
               XACRT = -106.05
            else if(i.eq.5.or.i.eq.17) then
               XACRT = 35.35
            else if(i.eq.6.or.i.eq.12.or.i.eq.18) then
               XACRT = 176.75
            endif

            if(i.eq.1.or.i.eq.2.or.i.eq.3) then
               YACRT = -176.75
            else if(i.eq.4.or.i.eq.5.or.i.eq.6) then
               YACRT = -106.05
            else if(i.eq.7.or.i.eq.8.or.i.eq.9) then
               YACRT = -35.35
            else if(i.eq.10.or.i.eq.12) then
               YACRT = 35.35
            else if(i.eq.13.or.i.eq.14.or.i.eq.15) then
               YACRT = 106.05
            else if(i.eq.16.or.i.eq.17.or.i.eq.18) then
               YACRT = 176.75
            endif

            CALL GSPOS('SPCT',i,'PXT2',XACRT,YACRT,ZACRT,0,'MANY')
            CALL GSPOS('SPCT',i,'PXT3',XACRT,YACRT,ZACRT,0,'MANY')
            CALL GSPOS('SPCT',i,'PXT9',XACRT,YACRT,ZACRT,0,'MANY')
            CALL GSPOS('SPCT',i,'PXT0',XACRT,YACRT,ZACRT,0,'MANY')
         endif
      enddo

* tenken
      do i=1,18
         if(i.ne.5.and.i.ne.9) then
            if(i.eq.1.or.i.eq.7.or.i.eq.13) then
               XACRT = -176.75
            else if(i.eq.2.or.i.eq.8.or.i.eq.14) then
               XACRT = -35.35
            else if(i.eq.3.or.i.eq.15) then
               XACRT = 106.05
            else if(i.eq.4.or.i.eq.10.or.i.eq.16) then
               XACRT = -106.05
            else if(i.eq.11.or.i.eq.17) then
               XACRT = 35.35
            else if(i.eq.6.or.i.eq.12.or.i.eq.18) then
               XACRT = 176.75
            endif

            if(i.eq.1.or.i.eq.2.or.i.eq.3) then
               YACRT = -176.75
            else if(i.eq.4.or.i.eq.6) then
               YACRT = -106.05
            else if(i.eq.7.or.i.eq.8) then
               YACRT = -35.35
            else if(i.eq.10.or.i.eq.11.or.i.eq.12) then
               YACRT = 35.35
            else if(i.eq.13.or.i.eq.14.or.i.eq.15) then
               YACRT = 106.05
            else if(i.eq.16.or.i.eq.17.or.i.eq.18) then
               YACRT = 176.75
            endif

            CALL GSPOS('SPCT',i,'PXTT',XACRT,YACRT,ZACRT,0,'MANY')
         endif
      enddo

      do i=1,9
         if(i.ne.5) then
            if(i.eq.1.or.i.eq.4.or.i.eq.7) then
               XACRT = 0.
            else if(i.eq.2.or.i.eq.8) then
               XACRT = -70.7
            else if(i.eq.3.or.i.eq.6.or.i.eq.9) then
               XACRT = 70.7
            endif

            if(i.eq.1) then
               YACRT = -176.75
            else if(i.eq.2.or.i.eq.3) then
               YACRT = -106.05
            else if(i.eq.4) then
               YACRT = -35.35
            else if(i.eq.6) then
               YACRT = 35.35
            else if(i.eq.7) then
               YACRT = 106.05
            else if(i.eq.8.or.i.eq.9) then
               YACRT = 176.75
            endif

            CALL GSPOS('SPCT',i,'PXT4',XACRT,YACRT,ZACRT,0,'MANY')
         endif
      enddo

      do i=1,3
         if(i.eq.1) then
            XACRT = -70.7
            YACRT = -35.35
         else if(i.eq.2) then
            XACRT = 70.7
            YACRT = -35.35
         else if(i.eq.3) then
            XACRT = 0.
            YACRT = 35.35
         endif
         CALL GSPOS('SPCT',i,'PXT5',XACRT,YACRT,ZACRT,0,'MANY')
         CALL GSPOS('SPCT',i,'PXT6',XACRT,YACRT,ZACRT,0,'MANY')
      enddo

      do i=1,6
         if(i.eq.1.or.i.eq.2.or.i.eq.3) then
            if(i.eq.1) then
               XACRT = -176.75
            else if(i.eq.2) then
               XACRT = -35.35
            else if(i.eq.3) then
               XACRT = 106.05
            endif
            YACRT = -35.35
         else if(i.eq.4.or.i.eq.5.or.i.eq.6) then
            if(i.eq.4) then
               XACRT = -106.05
            else if(i.eq.5) then
               XACRT = 35.35
            else if(i.eq.6) then
               XACRT = 176.75
            endif
            YACRT = 35.35
         endif
         CALL GSPOS('SPCT',i,'PXT7',XACRT,YACRT,ZACRT,0,'MANY')
      enddo

      XACRT = 0.
      YACRT = 0.
      CALL GSPOS('SPCT',1,'PXT8',XACRT,YACRT,ZACRT,0,'MANY')

*     Define Acrylic case on the top.
*     Medium number is IMDACR ( = 7).

      NPACR = 6
      PRACR(1) = RACR-WIDACR/2.
      PRACR(2) = RACR+WIDACR/2.
      PRACR(3) = 90.
      PRACR(4) = 180.
      PRACR(5) = 0.
      PRACR(6) = 360.
      CALL GSVOLU('ACRT','SPHE',IMDACR,PRACR,NPACR,IVACR)

      XACRT = 0.
      YACRT = 0.
      ZACRT = 0.
      CALL GSPOS('ACRT',1,'SPCT',XACRT,YACRT,ZACRT,0,'MANY')

*     Define a sphere which is surface of top PMT.
*     two kinds of SPHERE

* 1st
      NPPMTU = 6
      PRPMTU(1) = 0.
      PRPMTU(2) = RPMTU
      PRPMTU(3) = 90.
      PRPMTU(4) = 180.
      PRPMTU(5) = 0.
      PRPMTU(6) = 360.
      CALL GSVOLU('PMT1','SPHE',IMDWTS,PRPMTU,NPPMTU,IVPMTU)

      XPMTU = 0.
      YPMTU = 0.
      ZPMTU = RPMTU - 2*ZPMCY + (ZACCY+WIDACR/2)/2
     &        -RACR + ( (ZACCY+WIDACR/2)/2 - WIDACR/2 )

      CALL GSPOS('PMT1',1,'SPCT',XPMTU,YPMTU,ZPMTU,0,'MANY')

* 2nd
      NPPMTL = 6
      RPMTL  = SQRT((RPMT)**2+(APMTL)**2)
      PRPMTL(1) = 0.
      PRPMTL(2) = RPMTL
      PRPMTL(3) = 90.
      PRPMTL(4) = 180.
      PRPMTL(5) = 0.
      PRPMTL(6) = 360.
      CALL GSVOLU('PMTT','SPHE',IMDPMC,PRPMTL,NPPMTL,IVPMTL)

      XPMT2 = 0.
      YPMT2 = 0.
      ZPMT2 = - ((RPMTU - ZPMCY*2.) - APMTL)
      IRPMT2= 0
      CALL GSPOS('PMTT',1,'PMT1',XPMT2,YPMT2,ZPMT2,IRPMT2,'MANY')

************
**    Define the PMT volume on the bottom
************
c
c      NPPMTB = 3
c      PRPMTB(1) = 0.
c      PRPMTB(2) = RINTK
c      PRPMTB(3) = ZPMCY
c      CALL GSVOLU('PMPB','TUBE', IMDWTI, PRPMTB, NPPMTB, IVPMBP)
c
c*     Position it at the bottom
c
c      XPMTB = 0.
c      YPMTB = 0.
cC     ZPMTB = 0.
c      ZPMTB = -ZTANK - 0.5*HIINTK + ZPMCY
c
c      CALL GSPOS('PMPB',1,'INVO',XPMTB,YPMTB,ZPMTB,0,'ONLY')
c

*    define bottom plate

*   Z size and position of the plate (common)
C      PRBOX(3) = ZPMCY
      PRBOX(3) = (ZACCY+WIDACR/2.0-ZBLST)/2.0
C      ZPMTB = 0.
C      ZPMTB = -ZTANK - 0.5*HIINTK + ZPMCY
      ZPMTB = -ZTANK - 0.5*HIINTK +ZBLST + (ZACCY+WIDACR/2.0-ZBLST)/2.0

      NPBOX = 3
      PRBOX(1) = 24.*STEPX/2.
      PRBOX(2) = 36.*STEPY/2.
      CALL GSVOLU('PBXB','BOX ', IMDWTI, PRBOX, NPBOX, IVPBXT)

      XPMTB = 0.
      YPMTB = 0.
      CALL GSPOS('PBXB',1,'INVO',xpmtb,ypmtb,zpmtb,0,'MANY')

      CALL GSDVN('PXB1','PBXB',4,1)
      CALL GSDVN('PXB2','PXB1',6,2)

      NPBOX = 3
      PRBOX(1) = 6.*STEPX/2.
      PRBOX(2) = 24.*STEPY/2.
      CALL GSVOLU('PBB3','BOX ', IMDWTI, PRBOX, NPBOX, IVPBXT)

      XPMTB = 15.*STEPX
      YPMTB = 0.
      CALL GSPOS('PBB3',1,'INVO',xpmtb,ypmtb,zpmtb,0,'MANY')
      XPMTB = -15.*STEPX
      CALL GSPOS('PBB3',2,'INVO',xpmtb,ypmtb,zpmtb,0,'MANY')

      CALL GSDVN('PXB3','PBB3',4,2)

      NPBOX = 3
      PRBOX(1) = 3.*STEPX/2.
      PRBOX(2) = 24.*STEPY/2.
      CALL GSVOLU('PBB4','BOX ', IMDWTI, PRBOX, NPBOX, IVPBXT)
      CALL GSVOLU('PBB5','BOX ', IMDWTI, PRBOX, NPBOX, IVPBXT)

      XPMTB = -19.5*STEPX
      YPMTB = 0.
      CALL GSPOS('PBB4',1,'INVO',xpmtb,ypmtb,zpmtb,0,'MANY')

      CALL GSDVN('PXB4','PBB4',4,2)

      XPMTB = 19.5*STEPX
      YPMTB = 0.
      CALL GSPOS('PBB5',1,'INVO',xpmtb,ypmtb,zpmtb,0,'MANY')

      CALL GSDVN('PXB5','PBB5',12,2)

      NPBOX = 3
      PRBOX(1) = 3.*STEPX/2.
      PRBOX(2) = 6.*STEPY/2.
      CALL GSVOLU('PBB6','BOX ', IMDWTI, PRBOX, NPBOX, IVPBXT)

      XPMTB = 13.5*STEPX
      YPMTB = 15.*STEPY
      CALL GSPOS('PBB6',1,'INVO',xpmtb,ypmtb,zpmtb,0,'MANY')
      YPMTB = -15.*STEPY
      CALL GSPOS('PBB6',2,'INVO',xpmtb,ypmtb,zpmtb,0,'MANY')

      CALL GSDVN('PXB6','PBB6',3,2)

      NPBOX = 3
      PRBOX(1) = 24.*STEPX/2.
      PRBOX(2) = 2.*STEPY/2.
      CALL GSVOLU('PBB7','BOX ', IMDWTI, PRBOX, NPBOX, IVPBXT)

      XPMTB = 0.
      YPMTB = 19.*STEPY
      CALL GSPOS('PBB7',1,'INVO',xpmtb,ypmtb,zpmtb,0,'MANY')
      YPMTB = -19.*STEPY
      CALL GSPOS('PBB7',2,'INVO',xpmtb,ypmtb,zpmtb,0,'MANY')

      CALL GSDVN('PXB7','PBB7',4,1)

      NPBOX = 3
      PRBOX(1) = STEPX/2.
      PRBOX(2) = STEPY/2.
      CALL GSVOLU('PXB8','BOX ', IMDWTI, PRBOX, NPBOX, IVPBXT)

      inum=0
      do ix=-24,-1
         XPMTB = float(ix)*STEPX + STEPX/2.
         call topbot(3,ix,num,iyadd)
         do iy=1,num
            inum=inum+1
            if(iyadd(iy).le.0) then
               YPMTB = float(iyadd(iy))*STEPY + STEPY/2.
            else
               YPMTB = float(iyadd(iy))*STEPY - STEPY/2.
            endif
            CALL GSPOS('PXB8',inum,'INVO',xpmtb,ypmtb,zpmtb,0,'MANY')
         enddo
      enddo
      do ix=1,24
         XPMTB = float(ix)*STEPX - STEPX/2.
         call topbot(3,ix,num,iyadd)
         do iy=1,num
            inum=inum+1
            if(iyadd(iy).le.0) then
               YPMTB = float(iyadd(iy))*STEPY + STEPY/2.
            else
               YPMTB = float(iyadd(iy))*STEPY - STEPY/2.
            endif
            CALL GSPOS('PXB8',inum,'INVO',xpmtb,ypmtb,zpmtb,0,'MANY')
         enddo
      enddo

*     Define a water hemisphere between acrylic case and PMT surface.
*     Medium number is IMDWTS ( = 8).

      ZACRB = - ( RACR - ( (ZACCY+WIDACR/2-ZBLST)/2 - WIDACR/2 ))

      do i=1,18
         if(i.ne.8) then
            if(i.eq.1.or.i.eq.7.or.i.eq.13) then
               XACRB = -106.05
            else if(i.eq.2.or.i.eq.14) then
               XACRB = 35.35
            else if(i.eq.3.or.i.eq.9.or.i.eq.15) then
               XACRB = 176.75
            else if(i.eq.4.or.i.eq.10.or.i.eq.16) then
               XACRB = -176.75
            else if(i.eq.5.or.i.eq.11.or.i.eq.17) then
               XACRB = -35.35
            else if(i.eq.6.or.i.eq.12.or.i.eq.18) then
               XACRB = 106.05
            endif

            if(i.eq.1.or.i.eq.2.or.i.eq.3) then
               YACRB = -176.75
            else if(i.eq.4.or.i.eq.5.or.i.eq.6) then
               YACRB = -106.05
            else if(i.eq.7.or.i.eq.9) then
               YACRB = -35.35
            else if(i.eq.10.or.i.eq.11.or.i.eq.12) then
               YACRB = 35.35
            else if(i.eq.13.or.i.eq.14.or.i.eq.15) then
               YACRB = 106.05
            else if(i.eq.16.or.i.eq.17.or.i.eq.18) then
               YACRB = 176.75
            endif

            CALL GSPOS('SPCB',i,'PXB2',XACRB,YACRB,ZACRB,0,'MANY')
            CALL GSPOS('SPCB',i,'PXB3',XACRB,YACRB,ZACRB,0,'MANY')
         endif
      enddo

      do i=1,9
         if(i.ne.4) then
            if(i.eq.1.or.i.eq.7) then
               XACRB = -70.7
            else if(i.eq.2.or.i.eq.5.or.i.eq.8) then
               XACRB = 70.7
            else if(i.eq.3.or.i.eq.6.or.i.eq.9) then
               XACRB = 0.
            endif

            if(i.eq.1.or.i.eq.2) then
               YACRB = -176.75
            else if(i.eq.3) then
               YACRB = -106.05
            else if(i.eq.5) then
               YACRB = -35.35
            else if(i.eq.6) then
               YACRB = 35.35
            else if(i.eq.7.or.i.eq.8) then
               YACRB = 106.05
            else if(i.eq.9) then
               YACRB = 176.75
            endif
            CALL GSPOS('SPCB',i,'PXB4',XACRB,YACRB,ZACRB,0,'MANY')
         endif
      enddo

      do i=1,3
         if(i.eq.1) then
            XACRB = 0.
            YACRB = -35.35
         else if(i.eq.2) then
            XACRB = -70.7
            YACRB = 35.35
         else if(i.eq.3) then
            XACRB = 70.7
            YACRB = 35.35
         endif

         CALL GSPOS('SPCB',i,'PXB5',XACRB,YACRB,ZACRB,0,'MANY')
         CALL GSPOS('SPCB',i,'PXB6',XACRB,YACRB,ZACRB,0,'MANY')
      enddo

      do i=1,6
         if(i.eq.1.or.i.eq.2.or.i.eq.3) then
            if(i.eq.1) then
               XACRB = -106.05
            else if(i.eq.2) then
               XACRB = 35.35
            else if(i.eq.3) then
               XACRB = 176.75
            endif
            YACRB = -35.35
         else if(i.eq.4.or.i.eq.5.or.i.eq.6) then
            if(i.eq.4) then
               XACRB = -176.75
            else if(i.eq.5) then
               XACRB = -35.35
            else if(i.eq.6) then
               XACRB = 106.05
            endif
            YACRB = 35.35
         endif

         CALL GSPOS('SPCB',i,'PXB7',XACRB,YACRB,ZACRB,0,'MANY')
      enddo

      XACRB = 0.
      YACRB = 0.
      CALL GSPOS('SPCB',i,'PXB8',XACRB,YACRB,ZACRB,0,'MANY')

*     Define Acrylic case on the bottom
*     medium number is IMDACR ( = 7)

      NPACR    = 6
      PRACR(1) = RACR-WIDACR/2.
      PRACR(2) = RACR+WIDACR/2.
      PRACR(3) = 0.
      PRACR(4) = 90.
      PRACR(5) = 0.
      PRACR(6) = 360.
      CALL GSVOLU('ACRB','SPHE',IMDACR,PRACR,NPACR,IVACR)

      XACRB = 0.
      YACRB = 0.
      ZACRB = 0.
c      ZACRB = - ( RACR - ( (ZACCY+WIDACR/2)/2 - WIDACR/2 ))

      CALL GSPOS('ACRB',1,'SPCB',XACRB,YACRB,ZACRB,0,'MANY')


*     Define a sphere which is surface of bottom PMT.
*     two kinds of SPHERE

* 1st
      NPPMTU = 6
      PRPMTU(1) = 0.
      PRPMTU(2) = RPMTU
      PRPMTU(3) = 0.0
      PRPMTU(4) = 90.
      PRPMTU(5) = 0.
      PRPMTU(6) = 360.
      CALL GSVOLU('PMB1','SPHE',IMDWTS,PRPMTU,NPPMTU,IVPMTU)

      XPMTU = 0.
      YPMTU = 0.
c      ZPMTU = -( RPMTU - 2*ZPMCY + (ZACCY+WIDACR/2)/2)
      ZPMTU = -( RPMTU - 2*ZPMCY + (ZACCY+WIDACR/2)/2)
     &       + ( RACR - ( (ZACCY+WIDACR/2)/2 - WIDACR/2 ) )
      CALL GSPOS('PMB1',1,'SPCB',XPMTU,YPMTU,ZPMTU,0,'MANY')

* 2nd

      NPPMTL = 6
      RPMTL  = SQRT((RPMT)**2+(APMTL)**2)
      PRPMTL(1) = 0.
      PRPMTL(2) = RPMTL
      PRPMTL(3) = 0.
      PRPMTL(4) = 90.
      PRPMTL(5) = 0.
      PRPMTL(6) = 360.
      CALL GSVOLU('PMTB','SPHE',IMDPMC,PRPMTL,NPPMTL,IVPMTL)

      XPMT2 = 0.
      YPMT2 = 0.
      ZPMT2 = (RPMTU - ZPMCY*2.) - APMTL
      IRPMT2= 0
      CALL GSPOS('PMTB',1,'PMB1',XPMT2,YPMT2,ZPMT2,IRPMT2,'MANY')

**
CCCCC CALL GPROTM(0)

C     CALL GSORD('INVO',4)

c     CALL GPVOLU(0)
c     CALL GPSETS(0,0)

*     End of Geometry Initialization

c     CALL GSORD('PMTT',1)
      CALL GGCLOS

      RETURN
      END

      subroutine gxyzset

      implicit none

#include "skcnst.h"
#include "skparm.h"
#include "geotnk.h"
#include "geopmt.h"
#include "geopma.h"
#include "guxyzpm.h"

      integer ipmt, ixyz, ip, it, iz, icount, iy, ix, ixx, iyy, icol
      real stepz, zeroz, stept, zerot, zerox, zeroy, zz, tt, ttrad,
     &     xx, yy

      do ipmt=1,maxpma
         do ixyz=1,3
            gxyzpma(ixyz,ipmt) = xyzpma(ixyz,ipmt)
         enddo
      enddo

C SET BOUNDARYS FOR SERIAL NUMBER
C FIRST PMT # FOR WALL
      NFWAL=1
C LAST PMT # FOR WAL
      NLWAL=NFWAL+NTWAL-1
C FIRST PMT # FOR TOP
      NFTOP=NLWAL+1
C LAST PMT # FOR TOP
      NLTOP=NFTOP+NTTOP-1
C FIRST PMT # FOR BOTTOM
      NFBOT=NLTOP+1
C LAST PMT # FOR BOTTOM
      NLBOT=NFBOT+NTBOT-1

C BASIC SPACING AND ZERO BIASES

      STEPZ=HIWAL/FLOAT(MXWZPM)
      ZEROZ=0.
      STEPT=360./FLOAT(MXWTPM)
C     ZEROT=STEPT/2.
      ZEROT=0.
      ZEROX=STEPX/2.
      ZEROY=STEPY/2.

      DO 400 IP=1,MAXPM
      IF(IP.LE.NLWAL)THEN
C WALL K=1
         IT=(IP-1)/MXWZPM + 1
	 IF(IT.LE.0.OR.IT.GT.MXWTPM) THEN
	   CALL MSGOUT('GESTPM','Anomalous PMT number ', IP)
	   STOP
         END IF
	 IZ=IP-MXWZPM*(IT-1) - (MXWZPM/2 + 1)

** Change for SK-II (by Koshio)
         if(mod(ip,3).eq.1) then
            ZEROZ = 3.7
         else if(mod(ip,3).eq.2) then
            ZEROZ = 0.
         else if(mod(ip,3).eq.0) then
            ZEROZ = -3.7
         endif

         icol = (ip-1)/51 + 1
         if(icol.le.28*4) then
            if(mod(icol,4).eq.1) then
               ZEROT = 2./(DIINTK*PI) * 360.
            else if(mod(icol,4).eq.0) then
               ZEROT = -2./(DIINTK*PI) * 360.
            else
               ZEROT = 0.
            endif
         else if(icol.eq.28*4+1) then
            ZEROT = 2./(DIINTK*PI) * 360.
         else if(icol.eq.28*4+2) then
            ZEROT = -2./(DIINTK*PI) * 360.
         else
            if(mod(icol,4).eq.3) then
               ZEROT = 2./(DIINTK*PI) * 360.
            else if(mod(icol,4).eq.2) then
               ZEROT = -2./(DIINTK*PI) * 360.
            else
               ZEROT = 0.
            endif
         endif
*****

         ZZ=ZEROZ+IZ*STEPZ
         TT=ZEROT+(IT-1)*STEPT
         TTRAD=TT*3.141592/180
         XX=COS(TTRAD)*RINTK
         YY=SIN(TTRAD)*RINTK
         gXYZPM(1,IP)=XX
         gXYZPM(2,IP)=YY
         gXYZPM(3,IP)=ZZ
         gDXYZPM(1,IP)=-COS(TTRAD)
         gDXYZPM(2,IP)=-SIN(TTRAD)
         gDXYZPM(3,IP)=0.
      ELSE
         IF(IP.GE.NFTOP .AND. IP.LE.NLTOP) THEN
C TOP K=2
            IPMT=IP-NFTOP+1
            IF(IPMT.GT.NLTOP) THEN
               CALL MSGOUT('GESTPM','Anomalous PMT number 222',IPMT)
               STOP
            END IF

            ICOUNT=0

            DO 500 IY=-MXTYPM/2,MXTYPM/2
            IF(IY.NE.0) THEN
               DO 510 IX=-NMTYPM(IY)/2,NMTYPM(IY)/2
               IF(IX.NE.0) THEN
                  ICOUNT=ICOUNT+1
                  IF(IPMT.EQ.ICOUNT)THEN
C SELECTE IX,IY
                     IXX=IX
                     IYY=IY
                     IF(IX.LE.-1)IXX=IX+1
                     IF(IY.LE.-1)IYY=IY+1
                     gXYZPM(1,IP)=(IXX-1)*STEPX+ZEROX
                     gXYZPM(2,IP)=(IYY-1)*STEPY+ZEROY
                     gXYZPM(3,IP)=ZPINTK
                     gDXYZPM(1,IP)=0.
                     gDXYZPM(2,IP)=0.
                     gDXYZPM(3,IP)=-1.
                     GO TO 501
                  END IF
               END IF
 510           CONTINUE
            ENDIF
500         CONTINUE
501         CONTINUE
         ELSE
            IF(IP.GE.NFBOT)THEN
C BOTTOM K=3
               IPMT=IP-NFBOT+1
               IF(IPMT.GT.NLBOT) THEN
                  CALL MSGOUT('GESTPM','Anomalous PMT number 333',IPMT)
                  STOP
               END IF
               ICOUNT=0
               DO 600 IY=-MXBYPM/2,MXBYPM/2
               IF(IY.NE.0) THEN
                  DO 610 IX=-NMBYPM(IY)/2,NMBYPM(IY)/2
                  IF(IX.NE.0) THEN
                     ICOUNT=ICOUNT+1
                     IF(IPMT.EQ.ICOUNT)THEN
C SELECT IX,IY
                        IXX=IX
                        IYY=IY
                        IF(IX.LE.-1)IXX=IX+1
                        IF(IY.LE.-1)IYY=IY+1
                        gXYZPM(1,IP)=(IXX-1)*STEPX+ZEROX
                        gXYZPM(2,IP)=(IYY-1)*STEPY+ZEROY
                        gXYZPM(3,IP)=ZMINTK
                        gDXYZPM(1,IP)=0.
                        gDXYZPM(2,IP)=0.
                        gDXYZPM(3,IP)=1.
                        GO TO 601
                     ENDIF
                  END IF
  610             CONTINUE
               END IF
  600          CONTINUE
  601          CONTINUE
            ELSE
               CALL MSGOUT('GESTPM','Anomalous PMT number at loop 400',
     &                     IPMT)
               STOP
            ENDIF
         ENDIF
      ENDIF
400   CONTINUE

      return
      end

      subroutine topbot(i,ix,num,iyadd)
c
c  input
c     i     : 2 for top, or 3 for bottom
c     ix    :
c  output
c     num   :
c     iyadd :
c

      implicit none

      integer i,ix,num,iyadd(48)

      if(ix.gt.24.or.ix.lt.-24.or.ix.eq.0) then
         print*,'ix num is wrong in topbot'
         stop
      endif

      if(i.eq.2) then           ! for top
         if(ix.eq.-24) then
            num = 1
            iyadd(1) = 1
         else if(ix.eq.-23) then
            num = 6
            iyadd(1) = -5
            iyadd(2) = -3
            iyadd(3) = -1
            iyadd(4) = 2
            iyadd(5) = 4
            iyadd(6) = 6
         else if(ix.eq.-22) then
            num = 10
            iyadd(1) = -10
            iyadd(2) = -8
            iyadd(3) = -6
            iyadd(4) = -4
            iyadd(5) = -2
            iyadd(6) = 1
            iyadd(7) = 3
            iyadd(8) = 5
            iyadd(9) = 7
            iyadd(10)= 9
         else if(ix.eq.-21) then
            num = 0
         else if(ix.eq.-20) then
            num = 1
            iyadd(1) = 13
         else if(ix.eq.-19) then
            num = 3
            iyadd(1) = -15
            iyadd(2) = -13
            iyadd(3) = 14
         else if(ix.eq.-18) then
            num = 3
            iyadd(1) = -14
            iyadd(2) = 13
            iyadd(3) = 15
         else if(ix.eq.-17) then
            num = 4
            iyadd(1) = -15
            iyadd(2) = -13
            iyadd(3) = 14
            iyadd(4) = 16
         else if(ix.eq.-16) then
            num = 6
            iyadd(1) = -18
            iyadd(2) = -16
            iyadd(3) = -14
            iyadd(4) = 13
            iyadd(5) = 15
            iyadd(6) = 17
         else if(ix.eq.-15) then
            num = 5
            iyadd(1) = -19
            iyadd(2) = -17
            iyadd(3) = -13
            iyadd(4) = 14
            iyadd(5) = 18
         else if(ix.eq.-14) then
            num = 7
            iyadd(1) = -18
            iyadd(2) = -16
            iyadd(3) = -14
            iyadd(4) = 13
            iyadd(5) = 15
            iyadd(6) = 17
            iyadd(7) = 19
         else if(ix.eq.-13) then
            num = 8
            iyadd(1) = -19
            iyadd(2) = -17
            iyadd(3) = -15
            iyadd(4) = -13
            iyadd(5) = 14
            iyadd(6) = 16
            iyadd(7) = 18
            iyadd(8) = 20
         else if(ix.eq.-12) then
            num = 0
         else if(ix.eq.-11) then
            num = 1
            iyadd(1) = -21
         else if(ix.eq.-10) then
            num = 1
            iyadd(1) = 21
         else if(ix.eq.-9) then
            num = 1
            iyadd(1) = 22
         else if(ix.eq.-8) then
            num = 2
            iyadd(1) = -22
            iyadd(2) = 21
         else if(ix.eq.-7) then
            num = 2
            iyadd(1) = -21
            iyadd(2) = 22
         else if(ix.eq.-6) then
            num = 3
            iyadd(1) = -22
            iyadd(2) = 21
            iyadd(3) = 23
         else if(ix.eq.-5) then
            num = 3
            iyadd(1) = -23
            iyadd(2) = -21
            iyadd(3) = 22
         else if(ix.eq.-4) then
            num = 3
            iyadd(1) = -22
            iyadd(2) = 21
            iyadd(3) = 23
         else if(ix.eq.-3) then
            num = 1
            iyadd(1) = -23
         else if(ix.eq.-2) then
            num = 3
            iyadd(1) = -22
            iyadd(2) = 21
            iyadd(3) = 23
         else if(ix.eq.-1) then
            num = 4
            iyadd(1) = -23
            iyadd(2) = -21
            iyadd(3) = 22
            iyadd(4) = 24
         else if(ix.eq.1)  then
            num = 4
            iyadd(1) = -24
            iyadd(2) = -22
            iyadd(3) = 21
            iyadd(4) = 23
         else if(ix.eq.2)  then
            num = 3
            iyadd(1) = -23
            iyadd(2) = -21
            iyadd(3) = 22
         else if(ix.eq.3)  then
            num = 3
            iyadd(1) = -22
            iyadd(2) = 21
            iyadd(3) = 23
         else if(ix.eq.4)  then
            num = 1
            iyadd(1) = -23
         else if(ix.eq.5)  then
            num = 3
            iyadd(1) = -22
            iyadd(2) = 21
            iyadd(3) = 23
         else if(ix.eq.6)  then
            num = 3
            iyadd(1) = -23
            iyadd(2) = -21
            iyadd(3) = 22
         else if(ix.eq.7)  then
            num = 3
            iyadd(1) = -22
            iyadd(2) = 21
            iyadd(3) = 23
         else if(ix.eq.8)  then
            num = 2
            iyadd(1) = -21
            iyadd(2) = 22
         else if(ix.eq.9)  then
            num = 2
            iyadd(1) = -22
            iyadd(2) = 21
         else if(ix.eq.10) then
            num = 0
         else if(ix.eq.11) then
            num = 1
            iyadd(1) = 21
         else if(ix.eq.12) then
            num = 1
            iyadd(1) = -21
         else if(ix.eq.13) then
            num = 2
            iyadd(1) = -20
            iyadd(2) = 19
         else if(ix.eq.14) then
            num = 1
            iyadd(1) = -19
         else if(ix.eq.15) then
            num = 0
         else if(ix.eq.16) then
            num = 3
            iyadd(1) = -17
            iyadd(2) = -13
            iyadd(3) = 14
         else if(ix.eq.17) then
            num = 5
            iyadd(1) = -16
            iyadd(2) = -14
            iyadd(3) = 13
            iyadd(4) = 15
            iyadd(5) = 17
         else if(ix.eq.18) then
            num = 4
            iyadd(1) = -15
            iyadd(2) = -13
            iyadd(3) = 14
            iyadd(4) = 16
         else if(ix.eq.19) then
            num = 3
            iyadd(1) = -14
            iyadd(2) = 13
            iyadd(3) = 15
         else if(ix.eq.20) then
            num = 1
            iyadd(1) = -13
         else if(ix.eq.21) then
            num = 0
         else if(ix.eq.22) then
            num = 6
            iyadd(1) = -7
            iyadd(2) = -5
            iyadd(3) = -1
            iyadd(4) = 2
            iyadd(5) = 6
            iyadd(6) = 8
         else if(ix.eq.23) then
            num = 7
            iyadd(1) = -6
            iyadd(2) = -4
            iyadd(3) = -2
            iyadd(4) = 1
            iyadd(5) = 3
            iyadd(6) = 5
            iyadd(7) = 7
         else if(ix.eq.24) then
            num = 1
            iyadd(1) = -1
         endif

      else if(i.eq.3) then      ! for bottom
         if(ix.eq.-24) then
            num = 1
            iyadd(1) = -1
         else if(ix.eq.-23) then
            num = 6
            iyadd(1) = -6
            iyadd(2) = -4
            iyadd(3) = -2
            iyadd(4) = 1
            iyadd(5) = 3
            iyadd(6) = 5
         else if(ix.eq.-22) then
            num = 9
            iyadd(1) = -9
            iyadd(2) = -7
            iyadd(3) = -5
            iyadd(4) = -3
            iyadd(5) = -1
            iyadd(6) = 2
            iyadd(7) = 4
            iyadd(8) = 6
            iyadd(9) = 8
         else if(ix.eq.-21) then
            num = 0
         else if(ix.eq.-20) then
            num = 1
            iyadd(1) = -13
         else if(ix.eq.-19) then
            num = 2
            iyadd(1) = -14
            iyadd(2) = 13
         else if(ix.eq.-18) then
            num = 3
            iyadd(1) = -15
            iyadd(2) = -13
            iyadd(3) = 14
         else if(ix.eq.-17) then
            num = 4
            iyadd(1) = -16
            iyadd(2) = -14
            iyadd(3) = 13
            iyadd(4) = 15
         else if(ix.eq.-16) then
            num = 5
            iyadd(1) = -17
            iyadd(2) = -15
            iyadd(3) = -13
            iyadd(4) = 14
            iyadd(5) = 16
         else if(ix.eq.-15) then
            num = 4
            iyadd(1) = -18
            iyadd(2) = -14
            iyadd(3) = 13
            iyadd(4) = 17
         else if(ix.eq.-14) then
            num = 7
            iyadd(1) = -19
            iyadd(2) = -17
            iyadd(3) = -15
            iyadd(4) = -13
            iyadd(5) = 14
            iyadd(6) = 16
            iyadd(7) = 18
         else if(ix.eq.-13) then
            num = 8
            iyadd(1) = -20
            iyadd(2) = -18
            iyadd(3) = -16
            iyadd(4) = -14
            iyadd(5) = 13
            iyadd(6) = 15
            iyadd(7) = 17
            iyadd(8) = 19
         else if(ix.eq.-12) then
            num = 0
         else if(ix.eq.-11) then
            num = 1
            iyadd(1) = 21
         else if(ix.eq.-10) then
            num = 1
            iyadd(1) = -21
         else if(ix.eq.-9) then
            num = 0
         else if(ix.eq.-8) then
            num = 2
            iyadd(1) = -21
            iyadd(2) = 22
         else if(ix.eq.-7) then
            num = 2
            iyadd(1) = -22
            iyadd(2) = 21
         else if(ix.eq.-6) then
            num = 3
            iyadd(1) = -23
            iyadd(2) = -21
            iyadd(3) = 22
         else if(ix.eq.-5) then
            num = 3
            iyadd(1) = -22
            iyadd(2) = 21
            iyadd(3) = 23
         else if(ix.eq.-4) then
            num = 3
            iyadd(1) = -23
            iyadd(2) = -21
            iyadd(3) = 22
         else if(ix.eq.-3) then
            num = 1
            iyadd(1) = 23
         else if(ix.eq.-2) then
            num = 3
            iyadd(1) = -23
            iyadd(2) = -21
            iyadd(3) = 22
         else if(ix.eq.-1) then
            num = 3
            iyadd(1) = -22
            iyadd(2) = 21
            iyadd(3) = 23
         else if(ix.eq.1)  then
            num = 3
            iyadd(1) = -23
            iyadd(2) = -21
            iyadd(3) = 22
         else if(ix.eq.2)  then
            num = 3
            iyadd(1) = -22
            iyadd(2) = 21
            iyadd(3) = 23
         else if(ix.eq.3)  then
            num = 3
            iyadd(1) = -23
            iyadd(2) = -21
            iyadd(3) = 22
         else if(ix.eq.4)  then
            num = 1
            iyadd(1) = 23
         else if(ix.eq.5)  then
            num = 3
            iyadd(1) = -23
            iyadd(2) = -21
            iyadd(3) = 22
         else if(ix.eq.6)  then
            num = 3
            iyadd(1) = -22
            iyadd(2) = 21
            iyadd(3) = 23
         else if(ix.eq.7)  then
            num = 2
            iyadd(1) = -21
            iyadd(2) = 22
         else if(ix.eq.8)  then
            num = 2
            iyadd(1) = -22
            iyadd(2) = 21
         else if(ix.eq.9)  then
            num = 2
            iyadd(1) = -21
            iyadd(2) = 22
         else if(ix.eq.10) then
            num = 0
         else if(ix.eq.11) then
            num = 1
            iyadd(1) = -21
         else if(ix.eq.12) then
            num = 0
         else if(ix.eq.13) then
            num = 2
            iyadd(1) = -19
            iyadd(2) = 20
         else if(ix.eq.14) then
            num = 1
            iyadd(1) = 19
         else if(ix.eq.15) then
            num = 0
         else if(ix.eq.16) then
            num = 3
            iyadd(1) = -14
            iyadd(2) = 13
            iyadd(3) = 17
         else if(ix.eq.17) then
            num = 5
            iyadd(1) = -17
            iyadd(2) = -15
            iyadd(3) = -13
            iyadd(4) = 14
            iyadd(5) = 16
         else if(ix.eq.18) then
            num = 3
            iyadd(1) = -14
            iyadd(2) = 13
            iyadd(3) = 15
         else if(ix.eq.19) then
            num = 2
            iyadd(1) = -13
            iyadd(2) = 14
         else if(ix.eq.20) then
            num = 1
            iyadd(1) = 13
         else if(ix.eq.21) then
            num = 0
         else if(ix.eq.22) then
            num = 6
            iyadd(1) = -8
            iyadd(2) = -6
            iyadd(3) = -2
            iyadd(4) = 1
            iyadd(5) = 5
            iyadd(6) = 7
         else if(ix.eq.23) then
            num = 7
            iyadd(1) = -7
            iyadd(2) = -5
            iyadd(3) = -3
            iyadd(4) = -1
            iyadd(5) = 2
            iyadd(6) = 4
            iyadd(7) = 6
         else if(ix.eq.24) then
            num = 0
         endif
      else
         print*,'i num is wrong in topbot'
         stop
      endif

      return
      end
