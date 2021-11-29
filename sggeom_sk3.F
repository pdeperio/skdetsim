*********************************************************************
*     
*     -----------------
      SUBROUTINE SGGEOM_SK3
*     -----------------
*     
*     (Purpose)
*     Define the detector geometry and tracking media for SK-III
*     
*     (Input)
*     none
*     
*     (Output)
*     none
*     
*     (Creation Date and Author)
*     1992.04.24 ; First version by J.Kanzaki
*     1993.02.01 ; Modified by Y.Koshio
*     02.13 ; Change the fastest version (rotation)
*     1997.02.04 ; Modified to correct the dead space by R.Doyle
*     2006.04.15 ; SK-III version by G.Mitsuka
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
      REAL*4 PRPMTW(5),PRBOX(3),PRPMTT(5),PRPMTB(5)
      REAL*4 PRENDCAP(3)

C     REAL*4 PRPMTW(5),PRBOX(3)
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
      ZBLST = ZBLST_SK_III

      CALL GEOSET
      CALL GXYZSET_SK3
*     Stores the constants for the material in the data structure JMATE

c     CALL GSMATE(IMTWTR,'WATER$',
c     &            AWTR,ZWTR,DNWTR,RLWTR,ALWTR,UBWTR,NBWTR)

*     water
      AWTR(1)  =  AWTR1
      AWTR(2)  =  AWTR2
      ZWTR(1)  =  ZWTR1
      ZWTR(2)  =  ZWTR2
      WMWTR(1) =  WMWTR1
      WMWTR(2) =  WMWTR2
      CALL GSMIXT(IMTWTR,'WATER$',
     &     AWTR,ZWTR,DNWTR,NMWTR,WMWTR)

*     acryl

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
c     material for Ni source
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
     &     ISWTO,IFWTO,FLDWTO,TMXWTO,DMXWTO,DEXWTO,EPSWTO,STMWTO,
     &     UBWTO,NBWTO)

      CALL GSTMED(IMDWOE,'Outer Water endcap$',IMTWTR,
     &     ISWTO,IFWTO,FLDWTO,TMXWTO,DMXWTO,DEXWTO,EPSWTO,STMWTO,
     &     UBWTO,NBWTO)

*     following value is DUMMY
      ppckov(1) = 1.768E-09
      ppckov(2) = 4.126E-09
      absco(1)  = 6000.
      absco(2)  = 6000.
      effic(1)  = 0.05
      effic(2)  = 0.05
      rindex(1) = 1.344
      rindex(2) = 1.344
*     *
      call gsckov(imdwto,2,ppckov,absco,effic,rindex)
      call gsckov(imdwoe,2,ppckov,absco,effic,rindex)

*     Inner water volume

      CALL GSTMED(IMDWTI,'Inner Water volume$',IMTWTR,
     &     ISWTI,IFWTI,FLDWTI,TMXWTI,DMXWTI,DEXWTI,EPSWTI,STMWTI,
     &     UBWTI,NBWTI)

      call gsckov(imdwti,2,ppckov,absco,effic,rindex)

*     PMT cylinder volume

      CALL GSTMED(IMDPMC,'PMT Cylinder volume$',IMTVAC,
     &     ISPMC,IFPMC,FLDPMC,TMXPMC,DMXPMC,DEXPMC,EPSPMC,STMPMC,
     &     UBPMC,NBPMC)

*     The place between inner and anti volume (not sensitive)

      CALL GSTMED(IMDWTM,'mediam volume$',IMTWTR,
     &     ISWTM,IFWTM,FLDWTM,TMXWTM,DMXWTM,DEXWTM,EPSWTM,STMWTM,
     &     UBWTM,NBWTM)

*     Anti PMT and wave length shifter

      CALL GSTMED(IMDPMA,'anti PMA & W.S.$',IMTVAC,
     &     ISPMA,IFPMA,FLDPMA,TMXPMA,DMXPMA,DEXPMA,EPSPMA,STMPMA,
     &     UBPMA,NBPMA)

*     Out of the tank
#ifdef UPMU
*    Add rock component for upmu 2007 4/20 T.tanaka
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
     &     ISMST,IFMST,FLDMST,TMXMST,DMXMST,DEXMST,EPSMST,STMMST,
     &     UBMST,NBMST)
#endif
*     Acrylic case

      CALL GSTMED(IMDACR,'Acrylic case$',IMTACR,
     &     ISACR,IFACR,FLDACR,TMXACR,DMXACR,DEXACR,EPSACR,STMACR,
     &     UBACR,NBACR)

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
     &     ISWTS,IFWTS,FLDWTS,TMXWTS,DMXWTS,DEXWTS,EPSWTS,STMWTS,
     &     UBWTS,NBWTS)

      call gsckov(imdWTS,2,ppckov,absco,effic,rindex)

c     CALL GPTMED(0)

*     Define detector volumes

*     First define an imaginary box which includes the whole detector
*     and defines MAster Referecne System (MARS).
#ifdef UPMU
*----> saji modified 8km rock box (for upmu)
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

#ifdef UPMU
      CALL SGGEOM_DOME          ! Add SK dome to geometry
#endif

***   Separated top and bottom tank volumes for OD endcaps

      NPENDCAP    = 3
      PRENDCAP(1) = 0.
      PRENDCAP(2) = RINTK+RMED
      PRENDCAP(3) = 0.5*(0.5*HITKTK-(ZPINTK+ZMED))
      CALL GSVOLU('ODEC','TUBE', IMDWOE, PRENDCAP, NPENDCAP, IVENDCAP)
      XENDCAP = 0.
      YENDCAP = 0.
      ZENDCAP = ZPINTK + ZMED  + PRENDCAP(3)
      write(6,*)'prendcap',prendcap
      write(6,*)'x,y,zendcap',xendcap,yendcap,zendcap
      IRENDCAP = IRMAST
      CALL GSPOS('ODEC',1,'TANK',XENDCAP,YENDCAP,ZENDCAP,IRTANK,'ONLY')
      XENDCAP = 0.
      YENDCAP = 0.
      ZENDCAP = -ZENDCAP
      IRENDCAP = IRMAST
      CALL GSPOS('ODEC',2,'TANK',XENDCAP,YENDCAP,ZENDCAP,IRTANK,'ONLY')

***   MEDIUM
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

***   ANTI

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

*     center
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

*     side 1
      NPXANT = 3
      PRXANT(1) = 106.05
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

      CALL GSDVN('ANP1','ANT1',14,2)

*     side 2
      NPXANT = 3
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

      CALL GSDVN('AND2','ANT2',12,1)
      CALL GSDVN('ANP2','AND2',2,2)

*     side 3
      NPXANT = 3
      PRXANT(1) = 106.05
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

      CALL GSDVN('ANP3','ANT3',10,2)

*     side 4
      NPXANT = 3
      PRXANT(1) = 848.8
      PRXANT(2) = 110.0
      PRXANT(3) = ZPMCYA/2.
      CALL GSVOLU('ANT4','BOX ', IMDPMA, PRXANT, NPXANT, IVXANT)

      XXANT = 0.
      YXANT = 1131.2 + 282.8 + 110.0
      ZXANT = 0.
      IRXANT = IRMAST
      CALL GSPOS('ANT4',1,'ANTB',XXANT,YXANT,ZXANT,IRXANT,'ONLY')
      YXANT = -YXANT
      CALL GSPOS('ANT4',2,'ANTB',XXANT,YXANT,ZXANT,IRXANT,'ONLY')

      CALL GSDVN('AND4','ANT4',8,1)
      CALL GSDVN('ANP4','AND4',2,2)

***   INNER
*     Define the inner volume in the tank.
      NPINVO = 3
      PRINVO(1) = 0.
      PRINVO(2) = RINTK-ZBLST
      PRINVO(3) = 0.5*HIINTK-ZBLST
      CALL GSVOLU('INVO','TUBE', IMDWTI, PRINVO, NPINVO, IVINVO)
      XINVO = 0.
      YINVO = 0.
      ZINVO = -ZTANK
      IRINVO = IRMAST
      CALL GSPOS('INVO',1,'XMED',XINVO,YINVO,ZINVO,IRINVO,'ONLY')

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

*     main position

      NPPMTW = 5
      PRPMTW(1) = RINTK - ZACCY - WIDACR/2.0
      PRPMTW(2) = RINTK - ZBLST
      PRPMTW(3) = 0.5*(HIWAL-212.1)
      PRPMTW(4) = 268.8 + 4.8 ! start of SM30 i.e. end of SM29
      PRPMTW(5) = 360.0 + 268.8 ! end of SM28 i.e. start of SM29 + 2Pi
      CALL GSVOLU('PMPW','TUBS', IMDWTI, PRPMTW, NPPMTW, IVPMTW)

      XPMTW = 0.
      YPMTW = 0.
      ZPMTW = -35.35
      IRPMTW = 2
      CALL GSPOS('PMPW',1,'INVO',XPMTW,YPMTW,ZPMTW,IRPMTW,'ONLY')

      CALL GSDVN('PMCZ','PMPW',16,3)
      CALL GSDVN('PMCW','PMCZ',37,2)

      do i=1,12
         if(mod(i,3).eq.1) then
            ZPMTU = -70.7
         else if(mod(i,3).eq.2) then
            ZPMTU = 0.0 - zeroz
         else if(mod(i,3).eq.0) then
            ZPMTU = 70.7 + zeroz
         endif
         
         if(i.ge.1.and.i.lt.4) then
            degtub = -3.6 + zerot
         else if(i.ge.4.and.i.lt.7) then
            degtub = -1.2
         else if(i.ge.7.and.i.lt.10) then
            degtub = 1.2
         else if(i.ge.10.and.i.lt.13) then
            degtub = 3.6 - zerot
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

      enddo

*     extra position 1 (29th SM)

      NPPMTW = 5
      PRPMTW(1) = RINTK - ZACCY - WIDACR/2.0
      PRPMTW(2) = RINTK - ZBLST
      PRPMTW(3) = 0.5*(HIWAL-212.1)
      PRPMTW(4) = 268.8         ! 360.*((4*28)/150)
      PRPMTW(5) = 268.8 + 4.8   ! + 360.*(2/150)
      CALL GSVOLU('PMP1','TUBS', IMDWTI, PRPMTW, NPPMTW, IVPMTW)

      XPMTW = 0.
      YPMTW = 0.
      ZPMTW = -35.35
      IRPMTW = 2
      CALL GSPOS('PMP1',1,'INVO',XPMTW,YPMTW,ZPMTW,IRPMTW,'ONLY')

      CALL GSDVN('PMC1','PMP1',16,3)

      do i=1,6
         if(i.eq.1) then
            ZPMTU = -70.7
            degtub = 270.-360. + zerot
         else if(i.eq.2) then
            ZPMTU = 0.0 - zeroz
            degtub = 270.-360. + zerot
         else if(i.eq.3) then
            ZPMTU = 70.7 + zeroz
            degtub = 270.-360. + zerot
         else if(i.eq.4) then
            ZPMTU = -70.7
            degtub = 272.4-360.-zerot
         else if(i.eq.5) then
            ZPMTU = 0.0 - zeroz
            degtub = 272.4-360.-zerot
         else if(i.eq.6) then
            ZPMTU = 70.7 + zeroz
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

*     extra position 3

      NPPMTW = 3
      PRPMTW(1) = RINTK - ZACCY - WIDACR/2.0
      PRPMTW(2) = RINTK - ZBLST
C      PRPMTW(3) = 0.5*70.7
      PRPMTW(3) = 0.5*(0.5*HIINTK - ZBLST - (1767.5-0.5*70.7))
      CALL GSVOLU('PMP3','TUBE', IMDWTI, PRPMTW, NPPMTW, IVPMTW)

      XPMTW = 0.
      YPMTW = 0.
C      ZPMTW = -1767.5
      ZPMTW = -(1767.5-0.5*70.7) - PRPMTW(3)
      IRPMTW = 1
      CALL GSPOS('PMP3',1,'INVO',XPMTW,YPMTW,ZPMTW,IRPMTW,'ONLY')

      do i=1,150
         ZPMTU = (PRPMTW(3) - 70.7/2) + zeroz
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
      enddo

*     edge correction(Hut1 -> 2 -> 3 -> 4, Total 43 PMTs)
C      zacrb = -52.9
      zacrb = (PRPMTW(3) - 70.7/2) - 52.9
*     #10320
      xacrb = 1661.45
      yacrb = 35.35
      CALL GSPOS('SPCB',1,'PMP3',XACRB,YACRB,ZACRB,0,'MANY')
*     #10596
      xacrb = 1590.75
      yacrb = 459.55
      CALL GSPOS('SPCB',2,'PMP3',XACRB,YACRB,ZACRB,0,'MANY')
*     #10812
      xacrb = 1449.35
      yacrb = 813.05
      CALL GSPOS('SPCB',3,'PMP3',XACRB,YACRB,ZACRB,0,'MANY')
*     #10928
      xacrb = 1307.95
      yacrb = 1025.15
      CALL GSPOS('SPCB',4,'PMP3',XACRB,YACRB,ZACRB,0,'MANY')
*     #10964
      xacrb = 1237.25
      yacrb = 1095.85
      CALL GSPOS('SPCB',5,'PMP3',XACRB,YACRB,ZACRB,0,'MANY')
*     #10998
      xacrb = 1166.55
      yacrb = 1166.55
      CALL GSPOS('SPCB',6,'PMP3',XACRB,YACRB,ZACRB,0,'MANY')
*     #11030
      xacrb = 1095.85
      yacrb = 1237.25
      CALL GSPOS('SPCB',7,'PMP3',XACRB,YACRB,ZACRB,0,'MANY')
*     #11130
      xacrb = 671.65
      yacrb = 1520.05
      CALL GSPOS('SPCB',8,'PMP3',XACRB,YACRB,ZACRB,0,'MANY')
*     #11144
      xacrb = 459.55
      yacrb = 1590.75
      CALL GSPOS('SPCB',9,'PMP3',XACRB,YACRB,ZACRB,0,'MANY')
*     #11146
      xacrb = 35.35
      yacrb = 1661.45
      CALL GSPOS('SPCB',10,'PMP3',XACRB,YACRB,ZACRB,0,'MANY')
*     #11145
      xacrb = -35.35
      yacrb = 1661.45
      CALL GSPOS('SPCB',11,'PMP3',XACRB,YACRB,ZACRB,0,'MANY')
*     #11131
      xacrb = -459.55
      yacrb = 1590.75
      CALL GSPOS('SPCB',12,'PMP3',XACRB,YACRB,ZACRB,0,'MANY')
*     #11111
      xacrb = -671.65
      yacrb = 1520.05
      CALL GSPOS('SPCB',13,'PMP3',XACRB,YACRB,ZACRB,0,'MANY')
*     #10999     
      xacrb = -1095.85
      yacrb = 1237.25
      CALL GSPOS('SPCB',14,'PMP3',XACRB,YACRB,ZACRB,0,'MANY')
*     #10965
      xacrb = -1166.55
      yacrb = 1166.55
      CALL GSPOS('SPCB',15,'PMP3',XACRB,YACRB,ZACRB,0,'MANY')
*     #10929
      xacrb = -1237.25
      yacrb = 1095.85
      CALL GSPOS('SPCB',16,'PMP3',XACRB,YACRB,ZACRB,0,'MANY')
*     #10891
      xacrb = -1307.95
      yacrb = 1025.15
      CALL GSPOS('SPCB',17,'PMP3',XACRB,YACRB,ZACRB,0,'MANY')
*     #10771
      xacrb = -1449.35
      yacrb = 813.05
      CALL GSPOS('SPCB',18,'PMP3',XACRB,YACRB,ZACRB,0,'MANY')
*     #10551
      xacrb = -1590.75
      yacrb = 459.55
      CALL GSPOS('SPCB',19,'PMP3',XACRB,YACRB,ZACRB,0,'MANY')
*     #10273
      xacrb = -1661.45
      yacrb = 35.35
      CALL GSPOS('SPCB',20,'PMP3',XACRB,YACRB,ZACRB,0,'MANY')
*     #10225
      xacrb = -1661.45
      yacrb = -35.35
      CALL GSPOS('SPCB',21,'PMP3',XACRB,YACRB,ZACRB,0,'MANY')
*     #9949
      xacrb = -1590.75
      yacrb = -459.55
      CALL GSPOS('SPCB',22,'PMP3',XACRB,YACRB,ZACRB,0,'MANY')
*     #9733
      xacrb = -1449.35
      yacrb = -813.05
      CALL GSPOS('SPCB',23,'PMP3',XACRB,YACRB,ZACRB,0,'MANY')
*     #9617
      xacrb = -1307.95
      yacrb = -1025.15
      CALL GSPOS('SPCB',24,'PMP3',XACRB,YACRB,ZACRB,0,'MANY')
*     #9581
      xacrb = -1237.25
      yacrb = -1095.85
      CALL GSPOS('SPCB',25,'PMP3',XACRB,YACRB,ZACRB,0,'MANY')
*     #9547
      xacrb = -1166.55
      yacrb = -1166.55
      CALL GSPOS('SPCB',26,'PMP3',XACRB,YACRB,ZACRB,0,'MANY')
*     #9515
      xacrb = -1095.85
      yacrb = -1237.25
      CALL GSPOS('SPCB',27,'PMP3',XACRB,YACRB,ZACRB,0,'MANY')
*     #9485
      xacrb = -1025.15
      yacrb = -1307.95
      CALL GSPOS('SPCB',28,'PMP3',XACRB,YACRB,ZACRB,0,'MANY')
*     #9435
      xacrb = -813.05
      yacrb = -1449.35
      CALL GSPOS('SPCB',29,'PMP3',XACRB,YACRB,ZACRB,0,'MANY')
*     #9401
      xacrb = -459.55
      yacrb = -1590.75
      CALL GSPOS('SPCB',30,'PMP3',XACRB,YACRB,ZACRB,0,'MANY')
*     #9399
      xacrb = -35.35
      yacrb = -1661.45
      CALL GSPOS('SPCB',31,'PMP3',XACRB,YACRB,ZACRB,0,'MANY')
*     #9400
      xacrb = 35.35
      yacrb = -1661.45
      CALL GSPOS('SPCB',32,'PMP3',XACRB,YACRB,ZACRB,0,'MANY')
*     #9414
      xacrb = 459.55
      yacrb = -1590.75
      CALL GSPOS('SPCB',33,'PMP3',XACRB,YACRB,ZACRB,0,'MANY')
*     #9458
      xacrb = 813.05
      yacrb = -1449.35
      CALL GSPOS('SPCB',34,'PMP3',XACRB,YACRB,ZACRB,0,'MANY')
*     #9514
      xacrb = 1025.15
      yacrb = -1307.95
      CALL GSPOS('SPCB',35,'PMP3',XACRB,YACRB,ZACRB,0,'MANY')
*     #9546
      xacrb = 1095.85
      yacrb = -1237.25
      CALL GSPOS('SPCB',36,'PMP3',XACRB,YACRB,ZACRB,0,'MANY')
*     #9580
      xacrb = 1166.55
      yacrb = -1166.55
      CALL GSPOS('SPCB',37,'PMP3',XACRB,YACRB,ZACRB,0,'MANY')
*     #9616
      xacrb = 1237.25
      yacrb = -1095.85
      CALL GSPOS('SPCB',38,'PMP3',XACRB,YACRB,ZACRB,0,'MANY')
*     #9654
      xacrb = 1307.95
      yacrb = -1025.15
      CALL GSPOS('SPCB',39,'PMP3',XACRB,YACRB,ZACRB,0,'MANY')
*     #9774
      xacrb = 1449.35
      yacrb = -813.05
      CALL GSPOS('SPCB',40,'PMP3',XACRB,YACRB,ZACRB,0,'MANY')
*     #9860
      xacrb = 1520.05
      yacrb = -671.65
      CALL GSPOS('SPCB',41,'PMP3',XACRB,YACRB,ZACRB,0,'MANY')
*     #9994
      xacrb = 1590.75
      yacrb = -459.55
      CALL GSPOS('SPCB',42,'PMP3',XACRB,YACRB,ZACRB,0,'MANY')
*     #10272
      xacrb = 1661.45
      yacrb = -35.35
      CALL GSPOS('SPCB',43,'PMP3',XACRB,YACRB,ZACRB,0,'MANY')
      
*     extra position 4

      NPPMTW = 3
      PRPMTW(1) = RINTK - ZACCY - WIDACR/2.0
      PRPMTW(2) = RINTK - ZBLST
      PRPMTW(3) = 0.5*70.7
      CALL GSVOLU('PMP4','TUBE', IMDWTI, PRPMTW, NPPMTW, IVPMTW)

      XPMTW = 0.
      YPMTW = 0.
c     ZPMTW = -ZTANK
      ZPMTW = 1696.8
      IRPMTW = 1
      CALL GSPOS('PMP4',1,'INVO',XPMTW,YPMTW,ZPMTW,IRPMTW,'ONLY')

      do i=1,150
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
      enddo

*     extra position 5

      NPPMTW = 3
      PRPMTW(1) = RINTK - ZACCY - WIDACR/2.0
      PRPMTW(2) = RINTK - ZBLST
C      PRPMTW(3) = 0.5*70.7
      PRPMTW(3) = 0.5*(0.5*HIINTK - ZBLST - (1767.5-0.5*70.7))
      CALL GSVOLU('PMP5','TUBE', IMDWTI, PRPMTW, NPPMTW, IVPMTW)

      XPMTW = 0.
      YPMTW = 0.
C      ZPMTW = 1767.5
      ZPMTW = (1767.5-0.5*70.7) + PRPMTW(3)
      IRPMTW = 1
      CALL GSPOS('PMP5',1,'INVO',XPMTW,YPMTW,ZPMTW,IRPMTW,'ONLY')

      do i=1,150
C         ZPMTU = 0. - zeroz
         ZPMTU = -(PRPMTW(3) - 70.7/2) - zeroz
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
      enddo

*     edge correction(Hut1 -> 2 -> 3 -> 4, Total 44 PMTs)
C      zacrb = 52.9
      zacrb = -(PRPMTW(3) - 70.7/2) + 52.9
*     #8572
      xacrb = 1661.45
      yacrb = 35.35
      CALL GSPOS('SPCT',1,'PMP5',XACRB,YACRB,ZACRB,0,'MANY')
*     #8848
      xacrb = 1590.75
      yacrb = 459.55
      CALL GSPOS('SPCT',2,'PMP5',XACRB,YACRB,ZACRB,0,'MANY')
*     #9064
      xacrb = 1449.35
      yacrb = 813.05
      CALL GSPOS('SPCT',3,'PMP5',XACRB,YACRB,ZACRB,0,'MANY')
*     #9180
      xacrb = 1307.95
      yacrb = 1025.15
      CALL GSPOS('SPCT',4,'PMP5',XACRB,YACRB,ZACRB,0,'MANY')
*     #9216
      xacrb = 1237.25
      yacrb = 1095.85
      CALL GSPOS('SPCT',5,'PMP5',XACRB,YACRB,ZACRB,0,'MANY')
*     #9250
      xacrb = 1166.55
      yacrb = 1166.55
      CALL GSPOS('SPCT',6,'PMP5',XACRB,YACRB,ZACRB,0,'MANY')
*     #9282
      xacrb = 1095.85
      yacrb = 1237.25
      CALL GSPOS('SPCT',7,'PMP5',XACRB,YACRB,ZACRB,0,'MANY')
*     #9382
      xacrb = 671.65
      yacrb = 1520.05
      CALL GSPOS('SPCT',8,'PMP5',XACRB,YACRB,ZACRB,0,'MANY')
*     #9396
      xacrb = 459.55
      yacrb = 1590.75
      CALL GSPOS('SPCT',9,'PMP5',XACRB,YACRB,ZACRB,0,'MANY')
*     #9398
      xacrb = 35.35
      yacrb = 1661.45
      CALL GSPOS('SPCT',10,'PMP5',XACRB,YACRB,ZACRB,0,'MANY')
*     #9397
      xacrb = -35.35
      yacrb = 1661.45
      CALL GSPOS('SPCT',11,'PMP5',XACRB,YACRB,ZACRB,0,'MANY')
*     #9383
      xacrb = -459.55
      yacrb = 1590.75
      CALL GSPOS('SPCT',12,'PMP5',XACRB,YACRB,ZACRB,0,'MANY')
*     #9363
      xacrb = -671.65
      yacrb = 1520.05
      CALL GSPOS('SPCT',13,'PMP5',XACRB,YACRB,ZACRB,0,'MANY')
*     #9251
      xacrb = -1095.85
      yacrb = 1237.25
      CALL GSPOS('SPCT',14,'PMP5',XACRB,YACRB,ZACRB,0,'MANY')
*     #9217
      xacrb = -1166.55
      yacrb = 1166.55
      CALL GSPOS('SPCT',15,'PMP5',XACRB,YACRB,ZACRB,0,'MANY')
*     #9181
      xacrb = -1237.25
      yacrb = 1095.85
      CALL GSPOS('SPCT',16,'PMP5',XACRB,YACRB,ZACRB,0,'MANY')
*     #9143
      xacrb = -1307.95
      yacrb = 1025.15
      CALL GSPOS('SPCT',17,'PMP5',XACRB,YACRB,ZACRB,0,'MANY')
*     #9023
      xacrb = -1449.35
      yacrb = 813.05
      CALL GSPOS('SPCT',18,'PMP5',XACRB,YACRB,ZACRB,0,'MANY')
*     #8803
      xacrb = -1590.75
      yacrb = 459.55
      CALL GSPOS('SPCT',19,'PMP5',XACRB,YACRB,ZACRB,0,'MANY')
*     #8525
      xacrb = -1661.45
      yacrb = 35.35
      CALL GSPOS('SPCT',20,'PMP5',XACRB,YACRB,ZACRB,0,'MANY')
*     #8477
      xacrb = -1661.45
      yacrb = -35.35
      CALL GSPOS('SPCT',21,'PMP5',XACRB,YACRB,ZACRB,0,'MANY')
*     #8201
      xacrb = -1590.75
      yacrb = -459.55
      CALL GSPOS('SPCT',22,'PMP5',XACRB,YACRB,ZACRB,0,'MANY')
*     #8069
      xacrb = -1520.05
      yacrb = -671.65
      CALL GSPOS('SPCT',23,'PMP5',XACRB,YACRB,ZACRB,0,'MANY')
*     #7985
      xacrb = -1449.35
      yacrb = -813.05
      CALL GSPOS('SPCT',24,'PMP5',XACRB,YACRB,ZACRB,0,'MANY')
*     #7869
      xacrb = -1307.95
      yacrb = -1025.15
      CALL GSPOS('SPCT',25,'PMP5',XACRB,YACRB,ZACRB,0,'MANY')
*     #7833
      xacrb = -1237.25
      yacrb = -1095.85
      CALL GSPOS('SPCT',26,'PMP5',XACRB,YACRB,ZACRB,0,'MANY')
*     #7799
      xacrb = -1166.55
      yacrb = -1166.55
      CALL GSPOS('SPCT',27,'PMP5',XACRB,YACRB,ZACRB,0,'MANY')
*     #7767
      xacrb = -1095.85
      yacrb = -1237.25
      CALL GSPOS('SPCT',28,'PMP5',XACRB,YACRB,ZACRB,0,'MANY')
*     #7737
      xacrb = -1025.15
      yacrb = -1307.95
      CALL GSPOS('SPCT',29,'PMP5',XACRB,YACRB,ZACRB,0,'MANY')
*     #7687
      xacrb = -813.05
      yacrb = -1449.35
      CALL GSPOS('SPCT',30,'PMP5',XACRB,YACRB,ZACRB,0,'MANY')
*     #7653
      xacrb = -459.55
      yacrb = -1590.75
      CALL GSPOS('SPCT',31,'PMP5',XACRB,YACRB,ZACRB,0,'MANY')
*     #7651
      xacrb = -35.35
      yacrb = -1661.45
      CALL GSPOS('SPCT',32,'PMP5',XACRB,YACRB,ZACRB,0,'MANY')
*     #7652
      xacrb = 35.35
      yacrb = -1661.45
      CALL GSPOS('SPCT',33,'PMP5',XACRB,YACRB,ZACRB,0,'MANY')
*     #7666
      xacrb = 459.55
      yacrb = -1590.75
      CALL GSPOS('SPCT',34,'PMP5',XACRB,YACRB,ZACRB,0,'MANY')
*     #7710
      xacrb = 813.05
      yacrb = -1449.35
      CALL GSPOS('SPCT',35,'PMP5',XACRB,YACRB,ZACRB,0,'MANY')
*     #7766
      xacrb = 1025.15
      yacrb = -1307.95
      CALL GSPOS('SPCT',36,'PMP5',XACRB,YACRB,ZACRB,0,'MANY')
*     #7798
      xacrb = 1095.85
      yacrb = -1237.25
      CALL GSPOS('SPCT',37,'PMP5',XACRB,YACRB,ZACRB,0,'MANY')
*     #7832
      xacrb = 1166.55
      yacrb = -1166.55
      CALL GSPOS('SPCT',38,'PMP5',XACRB,YACRB,ZACRB,0,'MANY')
*     #7868
      xacrb = 1237.25
      yacrb = -1095.85
      CALL GSPOS('SPCT',39,'PMP5',XACRB,YACRB,ZACRB,0,'MANY')
*     #7906
      xacrb = 1307.95
      yacrb = -1025.15
      CALL GSPOS('SPCT',40,'PMP5',XACRB,YACRB,ZACRB,0,'MANY')
*     #8026
      xacrb = 1449.35
      yacrb = -813.05
      CALL GSPOS('SPCT',41,'PMP5',XACRB,YACRB,ZACRB,0,'MANY')
*     #8112
      xacrb = 1520.05
      yacrb = -671.65
      CALL GSPOS('SPCT',42,'PMP5',XACRB,YACRB,ZACRB,0,'MANY')
*     #8246
      xacrb = 1590.75
      yacrb = -459.55
      CALL GSPOS('SPCT',43,'PMP5',XACRB,YACRB,ZACRB,0,'MANY')
*     #8524
      xacrb = 1661.45
      yacrb = -35.35
      CALL GSPOS('SPCT',44,'PMP5',XACRB,YACRB,ZACRB,0,'MANY')

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
      IRACRW= 0
      CALL GSPOS('ACRW',1,'SPCW',XACRW,YACRW,ZACRW,IRACRW,'MANY')

*     Define a sphere which is surface of PMT.
*     two kinds of SPHERE
*     1st
      NPPMTU = 6
      PRPMTU(1) = 0.
      PRPMTU(2) = RPMTU
      PRPMTU(3) = 0.
      PRPMTU(4) = 90.
      PRPMTU(5) = 0.
      PRPMTU(6) = 360.
      CALL GSVOLU('PMW1','SPHE',IMDWTS,PRPMTU,NPPMTU,IVPMTU)

      XPMTU = 0.
      YPMTU = 0.
      ZPMTU = (RACR - ZACCY)-(RPMTU - ZPMCY*2.) 
      IRPMT1 = 0
      CALL GSPOS('PMW1',1,'SPCW',XPMTU,YPMTU,ZPMTU,IRPMT1,'MANY')

*     2nd
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
**     Define the PMT volume on the top
************
*     First, make the mother volume "PMPT"

      NPPMTT = 3
      PRPMTT(1) = 0.
      PRPMTT(2) = RINTK - ZACCY - WIDACR/2.0
      PRPMTT(3) = (ZACCY+WIDACR/2.0-ZBLST)/2.0
      CALL GSVOLU('PMPT','TUBE', IMDWTI, PRPMTT, NPPMTT, IVPMTP)

*     Position it at the top

      XPMTT = 0.
      YPMTT = 0.
      ZPMTT = -ZTANK + 0.5*HIINTK - ZBLST - (ZACCY+WIDACR/2.0-ZBLST)/2.0
      CALL GSPOS('PMPT',1,'INVO',XPMTT,YPMTT,ZPMTT,0,'ONLY')

*     Put top plate "PBXT" on PMPT and devide PBXT into MXTXPM*MXTYPM

      NPBOX = 3
      PRBOX(1) = FLOAT(MXTXPM)*STEPX/2.
      PRBOX(2) = FLOAT(MXTYPM)*STEPY/2.
      PRBOX(3) = (ZACCY+WIDACR/2.0-ZBLST)/2.0
      CALL GSVOLU('PBXT','BOX ', IMDWTI, PRBOX, NPBOX, IVPBXT)

      CALL GSPOS('PBXT',1,'PMPT',0.,0.,0.,0,'MANY')

      CALL GSDVN('PXT1','PBXT',MXTXPM,1)
      CALL GSDVN('PXT2','PXT1',MXTYPM,2)

*     Define a sphere which is surface of top PMT.

      XACRT = 0.
      YACRT = 0.
      ZACRT = RACR - ( (ZACCY+WIDACR/2-ZBLST)/2 - WIDACR/2 )
      CALL GSPOS('SPCT',1,'PXT2',XACRT,YACRT,ZACRT,0,'MANY')

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
*     1st
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
     &     -RACR + ( (ZACCY+WIDACR/2)/2 - WIDACR/2 )

      CALL GSPOS('PMT1',1,'SPCT',XPMTU,YPMTU,ZPMTU,0,'MANY')

*     2nd
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
*     First, make the mother volume "PMPB"

      NPPMTB = 3
      PRPMTB(1) = 0.
      PRPMTB(2) = RINTK - ZACCY - WIDACR/2.0
      PRPMTB(3) = (ZACCY+WIDACR/2.0-ZBLST)/2.0
      CALL GSVOLU('PMPB','TUBE', IMDWTI, PRPMTB, NPPMTB, IVPMBP)

*     Position it at the bottom

      XPMTB = 0.
      YPMTB = 0.
      ZPMTB = -ZTANK - 0.5*HIINTK +ZBLST + (ZACCY+WIDACR/2.0-ZBLST)/2.0
      CALL GSPOS('PMPB',1,'INVO',XPMTB,YPMTB,ZPMTB,0,'ONLY')

*     Put a top plate "PBXB" on PMPB and devide PBXB into MXTXPM*MXTYPM

      NPBOX = 3
      PRBOX(1) = FLOAT(MXTXPM)*STEPX/2.
      PRBOX(2) = FLOAT(MXTYPM)*STEPY/2.
      PRBOX(3) = (ZACCY+WIDACR/2.0-ZBLST)/2.0
      CALL GSVOLU('PBXB','BOX ', IMDWTI, PRBOX, NPBOX, IVPBXT)

      CALL GSPOS('PBXB',1,'PMPB',0.,0.,0.,0,'MANY')

      CALL GSDVN('PXB1','PBXB',MXTXPM,1)
      CALL GSDVN('PXB2','PXB1',MXTYPM,2)

*     Define a sphere which is surface of bottom PMT.

      XACRB = 0.
      YACRB = 0.
      ZACRB = - ( RACR - ( (ZACCY+WIDACR/2-ZBLST)/2 - WIDACR/2 ))
      CALL GSPOS('SPCB',1,'PXB2',XACRB,YACRB,ZACRB,0,'MANY')

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

      CALL GSPOS('ACRB',1,'SPCB',XACRB,YACRB,ZACRB,0,'MANY')


*     Define a sphere which is surface of bottom PMT.
*     two kinds of SPHERE
*     1st
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
      ZPMTU = -( RPMTU - 2*ZPMCY + (ZACCY+WIDACR/2)/2)
     &     + ( RACR - ( (ZACCY+WIDACR/2)/2 - WIDACR/2 ) )
      CALL GSPOS('PMB1',1,'SPCB',XPMTU,YPMTU,ZPMTU,0,'MANY')

*     2nd

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

*     *
CCCCC CALL GPROTM(0)

C     CALL GSORD('INVO',4)

c     CALL GPVOLU(0)
c     CALL GPSETS(0,0)

*     End of Geometry Initialization

c     CALL GSORD('PMTT',1)

      CALL GGCLOS

c      call grfile(21,'detector.geom','ON')

      RETURN
      END

**************************************************
      subroutine gxyzset_sk3

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

C     SET BOUNDARYS FOR SERIAL NUMBER
C     FIRST PMT # FOR WALL
      NFWAL=1
C     LAST PMT # FOR WAL
      NLWAL=NFWAL+NTWAL-1
C     FIRST PMT # FOR TOP
      NFTOP=NLWAL+1
C     LAST PMT # FOR TOP
      NLTOP=NFTOP+NTTOP-1
C     FIRST PMT # FOR BOTTOM
      NFBOT=NLTOP+1
C     LAST PMT # FOR BOTTOM
      NLBOT=NFBOT+NTBOT-1

C     BASIC SPACING AND ZERO BIASES

      STEPZ=HIWAL/FLOAT(MXWZPM)
      ZEROZ=0.
      STEPT=360./FLOAT(MXWTPM)
C     ZEROT=STEPT/2.
      ZEROT=0.
      ZEROX=STEPX/2.
      ZEROY=STEPY/2.

      DO 400 IP=1,MAXPM
         IF(IP.LE.NLWAL)THEN
C     WALL K=1
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
C     TOP K=2
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
C     SELECTE IX,IY
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
 510                 CONTINUE
                  ENDIF
 500           CONTINUE
 501           CONTINUE
            ELSE
               IF(IP.GE.NFBOT)THEN
C     BOTTOM K=3
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
C     SELECT IX,IY
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
 610                    CONTINUE
                     END IF
 600              CONTINUE
 601              CONTINUE
               ELSE
                  CALL MSGOUT('GESTPM','Anomalous PMT number at loop 400',
     &                 IPMT)
                  STOP
               ENDIF
            ENDIF
         ENDIF
 400  CONTINUE

      return
      end
