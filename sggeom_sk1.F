*********************************************************************
*
*     -----------------
      SUBROUTINE SGGEOM_SK1
*     -----------------
*
*     (Purpose)
*       Define the detector geometry and tracking media for SK-I
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
*
************************************************************************

C     COMMONs for geometries

#include "skparm.h"
#include "geotnk.h"
#include "geopmt.h"
#include "geopma.h"

C     COMMONs for geometries (sksim local)

#include "guskge.h"

C     COMMONs for materials

#include "gumate.h"
#ifdef UPMU
#include "gumed_u.h"
#else
#include "gumed.h"
#endif
#include "guxyzpm.h"
      integer i


      REAL*4 PAR(3),PRTANK(3),PRINVO(3),PRFIDU(3)
      REAL*4 PRPMTW(3),PRBOX(3),PRPMTT(5),PRPMTB(5)
      REAL*4 PRPMTU(6),PRPMTL(6)
      REAL*4 PRXMED(3),PRWMED(3),PRXANT(3)

      real ppckov(2),absco(2),effic(2),rindex(2)

**  set the value of Black-sheet 
**
      ZBLST = 0.0

      CALL GEOSET
      CALL GXYZSET_SK1

*     Stores the constants for the material in the data structure JMATE

c     CALL GSMATE(IMTWTR,'WATER$',
c    &            AWTR,ZWTR,DNWTR,RLWTR,ALWTR,UBWTR,NBWTR)

      AWTR(1)  =  AWTR1
      AWTR(2)  =  AWTR2
      ZWTR(1)  =  ZWTR1
      ZWTR(2)  =  ZWTR2
      WMWTR(1) =  WMWTR1
      WMWTR(2) =  WMWTR2
      CALL GSMIXT(IMTWTR,'WATER$',
     &            AWTR,ZWTR,DNWTR,NMWTR,WMWTR)

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

c     CALL GPTMED(0)

*     Define detector volumes

*     First define an imaginary box which includes the whole detector
*     and defines MAster Referecne System (MARS).
#ifdef UPMU
! *----> saji modified 8km rock box
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
      PRWMED(1) = RINTK
      PRWMED(2) = RINTK+RMED
C      PRWMED(3) = 0.5*HITKTK
      PRWMED(3) = ZPINTK+ZMED ! debug
      CALL GSVOLU('WMED','TUBE', IMDWTM, PRWMED, NPWMED, IVWMED)

      XWMED = 0.
      YWMED = 0.
      ZWMED = 0.
      IRWMED = IRMAST
      CALL GSPOS('WMED',1,'TANK',XWMED,YWMED,ZWMED,IRWMED,'ONLY')

*     Top & Bottom
      NPXMED = 3
      PRXMED(1) = 0.
      PRXMED(2) = RINTK
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
      PRXANT(2) = 848.4
      PRXANT(3) = ZPMCYA/2.
      CALL GSVOLU('ANT1','BOX ', IMDPMA, PRXANT, NPXANT, IVXANT)

      XXANT = 1272.6 + 106.05
      YXANT = 0.
      ZXANT = 0.
      IRXANT = IRMAST
      CALL GSPOS('ANT1',1,'ANTB',XXANT,YXANT,ZXANT,IRXANT,'ONLY')
      XXANT = -XXANT
      CALL GSPOS('ANT1',2,'ANTB',XXANT,YXANT,ZXANT,IRXANT,'ONLY')

      CALL GSDVN('ANP1','ANT1',12,2)

*       side 2
      NPXANT = 3
      PRXANT(1) = 1060.5
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

      CALL GSDVN('AND2','ANT2',10,1)
      CALL GSDVN('ANP2','AND2',2,2)

*       side 3
      NPXANT = 3
      PRXANT(1) = 106.05
      PRXANT(2) = 565.6
      PRXANT(3) = ZPMCYA/2.
      CALL GSVOLU('ANT3','BOX ', IMDPMA, PRXANT, NPXANT, IVXANT)

      XXANT = 1272.6 + 212.1 + 106.05
      YXANT = 0.
      ZXANT = 0.
      IRXANT = IRMAST
      CALL GSPOS('ANT3',1,'ANTB',XXANT,YXANT,ZXANT,IRXANT,'ONLY')
      XXANT = -XXANT
      CALL GSPOS('ANT3',2,'ANTB',XXANT,YXANT,ZXANT,IRXANT,'ONLY')

      CALL GSDVN('ANP3','ANT3',8,2)

*       side 4
      NPXANT = 3
      PRXANT(1) = 636.3
      PRXANT(2) = 141.4
      PRXANT(3) = ZPMCYA/2.
      CALL GSVOLU('ANT4','BOX ', IMDPMA, PRXANT, NPXANT, IVXANT)

      XXANT = 0.
      YXANT = 1131.2 + 282.8 + 141.4
      ZXANT = 0.
      IRXANT = IRMAST
      CALL GSPOS('ANT4',1,'ANTB',XXANT,YXANT,ZXANT,IRXANT,'ONLY')
      YXANT = -YXANT
      CALL GSPOS('ANT4',2,'ANTB',XXANT,YXANT,ZXANT,IRXANT,'ONLY')

      CALL GSDVN('AND4','ANT4',6,1)
      CALL GSDVN('ANP4','AND4',2,2)

*** INNER

*     Define the inner volume in the tank.

      NPINVO = 3
      PRINVO(1) = 0.
      PRINVO(2) = RINTK
      PRINVO(3) = 0.5*HIINTK
      CALL GSVOLU('INVO','TUBE', IMDWTI, PRINVO, NPINVO, IVINVO)
      XINVO = 0.
      YINVO = 0.
C     ZINVO = 0.
      ZINVO = -ZTANK
      IRINVO = IRMAST
      CALL GSPOS('INVO',1,'XMED',XINVO,YINVO,ZINVO,IRINVO,'ONLY')

*     Define the pseudo-fiducial volume in the tank.

      NPFIDU = 3
      PRFIDU(1) = 0.
      PRFIDU(2) = RINTK - 2.*ZPMCY
      PRFIDU(3) = 0.5*HIINTK - 2.*ZPMCY
      CALL GSVOLU('FIDU','TUBE', IMDWTI, PRFIDU, NPFIDU, IVFIDU)
      XFIDU = XINVO
      YFIDU = YINVO
      ZFIDU = ZINVO
      IRFIDU = IRMAST
      CALL GSPOS('FIDU',1,'INVO',XFIDU,YFIDU,ZFIDU,IRFIDU,'ONLY')

*     Define the PMT volume on the wall

      NPPMTW = 3
      PRPMTW(1) = RINTK - ZPMCY*2.
      PRPMTW(2) = RINTK
      PRPMTW(3) = 0.5*HIWAL
      CALL GSVOLU('PMPW','TUBE', IMDWTI, PRPMTW, NPPMTW, IVPMTW)

      XPMTW = 0.
      YPMTW = 0.
C     ZPMTW = 0.
      ZPMTW = -ZTANK
      IRPMTW = 2
      CALL GSPOS('PMPW',1,'INVO',XPMTW,YPMTW,ZPMTW,IRPMTW,'ONLY')

*     Divide this volume into cells

      CALL GSDVN('PMCZ','PMPW',MXWZPM,3)
C     CALL GSORD('PMCZ',3)

      CALL GSDVN('PMCW','PMCZ',MXWTPM,2)
C     CALL GSORD('PMCW',6)

*     Rotation matrices for PMT's on the wall (IRWALL)

      IRWALL = 5
      CALL GSROTM(IRWALL, 90., 90., 180., 0., 90., 180.)

*     Define a sphere which is surface of PMT.

      NPPMTU = 6
      PRPMTU(1) = 0.
      PRPMTU(2) = RPMTU
      PRPMTU(3) = 0.
      PRPMTU(4) = 90.
      PRPMTU(5) = 0.
      PRPMTU(6) = 360.
      CALL GSVOLU('PMW1','SPHE',IMDWTI,PRPMTU,NPPMTU,IVPMTU)

      XPMTU = RINTK + (RPMTU - ZPMCY*2.)
      YPMTU = 0.
      ZPMTU = 0.
      IRPMT1= IRWALL
      CALL GSPOS('PMW1',1,'PMCW',XPMTU,YPMTU,ZPMTU,IRPMT1,'MANY')

*

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

**    Define the PMT volume on the top

      NPPMTT = 3
      PRPMTT(1) = 0.
      PRPMTT(2) = RINTK
      PRPMTT(3) = ZPMCY
      CALL GSVOLU('PMPT','TUBE', IMDWTI, PRPMTT, NPPMTT, IVPMTP)

*     Position it at the top

      XPMTT = 0.
      YPMTT = 0.
C     ZPMTT = 0.
      ZPMTT = -ZTANK + 0.5*HIINTK - ZPMCY

      CALL GSPOS('PMPT',1,'INVO',XPMTT,YPMTT,ZPMTT,0,'ONLY')

*     Devide top plate

      NPBOX = 3
      PRBOX(1) = FLOAT(MXTXPM)*STEPX/2.
      PRBOX(2) = FLOAT(MXTYPM)*STEPY/2.
      PRBOX(3) = ZPMCY
      CALL GSVOLU('PBXT','BOX ', IMDWTI, PRBOX, NPBOX, IVPBXT)

      CALL GSPOS('PBXT',1,'PMPT',0.,0.,0.,0,'MANY')

      CALL GSDVN('PXT1','PBXT',MXTXPM,1)
      CALL GSDVN('PXT2','PXT1',MXTYPM,2)

*     Define a sphere which is surface of top PMT.

      NPPMTU = 6
      PRPMTU(1) = 0.
      PRPMTU(2) = RPMTU
      PRPMTU(3) = 90.
      PRPMTU(4) = 180.
      PRPMTU(5) = 0.
      PRPMTU(6) = 360.
      CALL GSVOLU('PMT1','SPHE',IMDWTI,PRPMTU,NPPMTU,IVPMTU)

      XPMTU = 0.
      YPMTU = 0.
      ZPMTU = RPMTU - ZPMCY
      CALL GSPOS('PMT1',1,'PXT2',XPMTU,YPMTU,ZPMTU,0,'MANY')

*

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

**    Define the PMT volume on the bottom

      NPPMTB = 3
      PRPMTB(1) = 0.
      PRPMTB(2) = RINTK
      PRPMTB(3) = ZPMCY
      CALL GSVOLU('PMPB','TUBE', IMDWTI, PRPMTB, NPPMTB, IVPMBP)

*     Position it at the bottom

      XPMTB = 0.
      YPMTB = 0.
C     ZPMTB = 0.
      ZPMTB = -ZTANK - 0.5*HIINTK + ZPMCY

      CALL GSPOS('PMPB',1,'INVO',XPMTB,YPMTB,ZPMTB,0,'ONLY')

*     Devide bottom plate

      NPBOX = 3
      PRBOX(1) = FLOAT(MXTXPM)*STEPX/2.
      PRBOX(2) = FLOAT(MXTYPM)*STEPY/2.
      PRBOX(3) = ZPMCY
      CALL GSVOLU('PBXB','BOX ', IMDWTI, PRBOX, NPBOX, IVPBXT)

      CALL GSPOS('PBXB',1,'PMPB',0.,0.,0.,0,'MANY')

      CALL GSDVN('PXB1','PBXB',MXTXPM,1)
      CALL GSDVN('PXB2','PXB1',MXTYPM,2)

*     Define a sphere which is surface of bottom PMT.

      NPPMTU = 6
      PRPMTU(1) = 0.
      PRPMTU(2) = RPMTU
      PRPMTU(3) = 0.
      PRPMTU(4) = 90.
      PRPMTU(5) = 0.
      PRPMTU(6) = 360.
      CALL GSVOLU('PMB1','SPHE',IMDWTI,PRPMTU,NPPMTU,IVPMTU)

      XPMTU = 0.
      YPMTU = 0.
      ZPMTU = - (RPMTU - ZPMCY)
      CALL GSPOS('PMB1',1,'PXB2',XPMTU,YPMTU,ZPMTU,0,'MANY')

*

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
 
C -----------------------------------------------     
      Subroutine GXYZSET_SK1
C
C     Copies geoset values xyzpm, dxyzpm, etc. to local geant variables
C
C ------------------------------------------------
      implicit none
#include "skparm.h"
#include "guxyzpm.h"
#include "geopmt.h"
#include "geopma.h"

      integer I, J
      
      DO I=1,MAXPM
         DO J=1,3
            GXYZPM(J,I) = XYZPM(J,I)
            GDXYZPM(J,I) = DXYZPM(J,I)
         ENDDO
      ENDDO

      DO I=1,MAXPMA
         DO J=1,3
            GXYZPMA(J,I) = XYZPMA(J,I)
         ENDDO
      ENDDO

      RETURN
      END
