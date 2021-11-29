************************************************************************
*     -----------------------
      SUBROUTINE SGPMT(NOPMT)
*     -----------------------
*
*     (Purpose)
*     ID PMT
*     (Input)
*       NOPMT::  PMT cable number
*
*     (Output)
*       none
*
*
*     (Creation Date and Author)
*       1992.12.02 ; First version by Y.Koshio
*
C     Feb.-5-97  poralization is considered by M.Shiozawa
C     Mar-9-98 BUG fixed / add #include "gconst.h" for clight by y.Itow
C     APR-18-2006 Added SK-III case, still references SK-II QE table by K. Bays
************************************************************************
      implicit none
#include "skdetsim_parm.h"
#include "skhead.h"
#include "geant321/gctrak.inc"
#include "geant321/gcvolu.inc"
#include "geant321/gckine.inc"
#include "geant321/gcking.inc"
#include "skparm.h"
#include "guxyzpm.h"
#include "skmc.h"
#include "gupht.h"
#include "guhigh.h"
#include "dsbnkf.h"
#include "guflag.h"
#include "dsparm.h"
#include "guref.h"
#include "geant321/gconst.inc"
#include "skdetsim_test.h"

#ifdef FADC      
#include "fadcwrbnk.h"
#endif

      INTEGER NOPMT
      integer IDATA(2)
      real RDATA(2),VECNEW(3),POLARNEW(3)
      EQUIVALENCE (IDATA(1),RDATA(1))
      real pos2(3),dir2(3),pos3(3),dirran(3),posran(3)
      real FIXSTEP
      PARAMETER(FIXSTEP = .1)
      real blamda,refwtr,effnsg,costh,rdummy,PMTQESG,TOFGSK,del
      integer ihit,IFLGREF, I, isegm, ndata
      character*4 cname
      real PRABSL0

      character*128 fname, fqetab
      integer icall,ipmt,icab
      real    factor(maxpm),qetable
      real    collpmt(maxpm)
      data    icall/0/
      save    icall,collpmt,factor


      integer ikeksk
      real qe, norm, abc(3), pol

      CALL UHTOC(NAMES(5),4,CNAME,4)
      IF(NOPMT.LE.0.OR.NOPMT.GT.MAXPM) THEN
         WRITE(6,1234) CNAME,NUMBER(7),NUMBER(8)
 1234    FORMAT(' WARNING in SGPMT',A5,2I3)
         RETURN
      ENDIF     


***  correction of actual PMT shape 
      if(SK_GEOMETRY.eq.SK_I .and. REFLECTION_TUNE.eq.LE_TUNE) then
         call skdonuts_sk1(gxyzpm(1,nopmt),gdxyzpm(1,nopmt),
     &        vect,vect(4),pos2,dir2,pos3,ihit)
      else
         call skdonuts(gxyzpm(1,nopmt),gdxyzpm(1,nopmt),
     &        vect,vect(4),pos2,dir2,pos3,dirran,posran,ihit)
      endif

c      write(*,*) 'PMT Number =',NOPMT
c      write(*,*) 'Pos of PMT=',xyzpm(1,nopmt),xyzpm(2,nopmt),xyzpm(3,nopmt)
c      write(*,*) 'IHIT=',ihit

      if(ihit.eq.0)then
        ngphot=1
        step=sqrt((pos2(1)-vect(1))**2+
     &            (pos2(2)-vect(2))**2+
     &            (pos2(3)-vect(3))**2)+fixstep
        xphot(1,1)=vect(1)+vect(4)*step
        xphot(2,1)=vect(2)+vect(5)*step
        xphot(3,1)=vect(3)+vect(6)*step
        xphot(4,1)=vect(4)
        xphot(5,1)=vect(5)
        xphot(6,1)=vect(6)
        BLAMDA = 1237.8E-09 / GEKIN
        refwtr=effnsg(blamda)
c       write(6,*)'step,blamda,refwtr=',step,blamda,refwtr
        xphot(7,1)=vect(7)
        if(SK_GEOMETRY .ge. SK_V) then ! add polarization from SK-V
           xphot(8,1)=polar(1)
           xphot(9,1)=polar(2)
           xphot(10,1)=polar(3)
        endif
        xphot(11,1)=tofg+step/clight*refwtr

        if (OUTPUT_SCATTABLE.eq.1) then
           call gskpho(-ISTORY) ! track photon scattering and reflections
        else
           call gskpho(0)
        endif

        ngphot=0
c       write(6,*)nopmt,vect
        return
      else if ( ihit.eq.2 ) then    
**  absorption by PMT band
         istop=1
         return
      endif

      if(SK_GEOMETRY.eq.SK_I .and. REFLECTION_TUNE.eq.LE_TUNE) then
**  50% reflection by PMT band in ledetsim SK-1
         CALL UHTOC(NAMES(5),4,CNAME,4)
         IF(CNAME.EQ."PMPW" .and. 
     &        (gdxyzpm(1,nopmt)*(pos3(1)-gxyzpm(1,nopmt))
     &        +gdxyzpm(2,nopmt)*(pos3(2)-gxyzpm(2,nopmt))
     &        +gdxyzpm(3,nopmt)*(pos3(3)-gxyzpm(3,nopmt))).le.2.) then
c turn on/off stainless band
            CALL GRNDM(RDUMMY,1)
            if(RDUMMY.lt.0.5) then ! on the stainless band 50% of photons are
               NGPHOT = 1       ! reflected
               step=sqrt((pos2(1)-pos3(1))**2+
     &              (pos2(2)-pos3(2))**2+
     &              (pos2(3)-pos3(3))**2)+
     &              sqrt((pos2(1)-vect(1))**2+
     &              (pos2(2)-vect(2))**2+
     &              (pos2(3)-vect(3))**2)+fixstep
               norm = sqrt(dir2(1)**2+dir2(2)**2+dir2(3)**2)
               if(norm.gt.0)then
                  dir2(1) = dir2(1)/norm
                  dir2(2) = dir2(2)/norm
                  dir2(3) = dir2(3)/norm
               else
                  write(*,*)'ERROR in sgpmt dir2 = 0'
                  stop
               end if
               DO I = 1,3
                  XPHOT(I+3,1) = dir2(I)
                  XPHOT(I,1) = pos3(I) + fixstep*dir2(I)
               ENDDO

               BLAMDA = 1237.8E-09 / GEKIN
               refwtr=effnsg(blamda)
               xphot(7,1)=VECT(7)
               XPHOT(11,1) = TOFG+step/clight*refwtr
**** porlarization is not yet
               if (OUTPUT_SCATTABLE.eq.1) then
                  CALL GSKPHO(-ISTORY-20000000) ! track photon scattering and reflections
               else
                  CALL GSKPHO(0)
               endif

               NGPHOT=0
            else
               istop=1
            endif
c ---end of steel band
            return
         endif
      endif

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C   >=SK-III && ! TRACK_QE_PHOTON && ! AP_REFLECTION
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      IF ( 
     &   (REFLECTION_TUNE .eq. UNI_TUNE) .or.
     &  ((REFLECTION_TUNE .eq. LE_TUNE) .and.
     &  ((SK_GEOMETRY .eq. SK_III) .or. (SK_GEOMETRY .eq. SK_IV) .or. (SK_GEOMETRY .eq. SK_V) )
     &  )) then
 
*----------------------------------------
* apply QE table
*----------------------------------------
 
         if(icall.eq.0) then
            icall = 1
            ipmt=1

          
           if(SK_GEOMETRY .ge. SK_III) then
            fqetab='qetable'//cqetab(1:3)//'.dat'
            call findconsts(fname,'QE table',
c     $           'qetable3_0.dat',0)
     $           fqetab,0)
            open(unit=91,file=fname,status='old',
     $           form='formatted')
 410        continue
            read(91,*,end=411) icab, qetable
            factor(ipmt)=qetable
            if (OUTPUT_SCATTABLE.eq.1) then
               factor(ipmt) = 1.0
            endif
            ipmt=ipmt+1
            goto 410
 
 411        continue
            close(91)
            if(ipmt.ne.maxpm+1) then
               write(6,*) 'Entry of QE table is wrong !'
               stop
            endif

           else if(SK_GEOMETRY.eq.SK_II) then
               if(icall.eq.0) then
                  icall = 1
                  ipmt=1

                  call findconsts(fname,'QE table',
     $                 'qetable_1.dat',0)
                  open(unit=91,file=fname,status='old',
     $                 form='formatted')
 220              continue
                  read(91,*,end=221) icab, qetable
                  factor(ipmt)=qetable
                  if (OUTPUT_SCATTABLE.eq.1) then
                     factor(ipmt)=1.0
                  endif
                  ipmt=ipmt+1
                  goto 220

 221              continue
                  close(91)
                  if(ipmt.ne.maxpm+1) then
                     write(6,*) 'Entry of QE table is wrong !'
                     stop
                  endif
               endif
           else if(SK_GEOMETRY.eq.SK_I) then
              do ipmt = 1, MAXPM
                  factor(ipmt)=1.0
              enddo
           endif

         endif

         BLAMDA = 1237.8E-09 / GEKIN

         if(TRACK_PHOTON.eq.TRACK_ALL_PHOTON) then
            QE=pmtqesg(blamda)*factor(NOPMT)*COREPMT/PRABSL0(blamda)
         elseif (TRACK_PHOTON.eq.TRACK_QE_OVER_ABS_PHOTON) then ! TRACK_QE_OVER_ABS_PHOTON
            QE=factor(NOPMT)/maxfct
         elseif (TRACK_PHOTON.eq.TRACK_QE_OVER_ABS_SAVE_MEMORY) then 
            QE=factor(NOPMT)/maxfct
         elseif  (TRACK_PHOTON.eq.TRACK_QE_PHOTON) then 
            QE=factor(NOPMT)/maxfct/prabsl0(blamda)*0.75
         else
            write(*,*) 'SGPMT: TRACK_PHOTON is irregular!!: ',TRACK_PHOTON
            stop
         endif

         if (qe .gt. 1.0) then
            write(*,*) 'SGPMT: QE is greater than 1. Irregal.'
         endif

         CALL SGREFPMT(NOPMT,VECT,QE,dir2,costh,IFLGREF)
    
*  If photon is absorbed w/o p.e.(5) or transmitted(4), do nothing.
* Just to be fast, this condition is here. But this effect is very small.
C         IF(IFLGREF.eq.5 ) THEN
         IF(IFLGREF.eq.5 .or. IFLGREF.eq.4 ) THEN

*  If photon is reflected on PMT (specular) (1)
         ELSE IF(IFLGREF.EQ.1) THEN
**** porlarization is not yet

            NGPHOT = 1

            step=sqrt((pos2(1)-pos3(1))**2+
     &           (pos2(2)-pos3(2))**2+
     &           (pos2(3)-pos3(3))**2)+
     &           sqrt((pos2(1)-vect(1))**2+
     &           (pos2(2)-vect(2))**2+
     &           (pos2(3)-vect(3))**2)+fixstep
            norm = sqrt(dir2(1)**2+dir2(2)**2+dir2(3)**2)
            if(norm.gt.0)then
               dir2(1) = dir2(1)/norm
               dir2(2) = dir2(2)/norm
               dir2(3) = dir2(3)/norm
            else
               write(*,*)'ERROR in sgpmt dir2 = 0'
               stop
            end if
            DO I = 1,3
               XPHOT(I+3,1) = dir2(I)
               XPHOT(I,1) = pos3(I) + fixstep*dir2(I)
            ENDDO

            refwtr=effnsg(blamda)
            xphot(7,1)=VECT(7)
C
C polarization (refer from random reflection formular in Geant4) from SK-V
C
            if(SK_GEOMETRY .ge. SK_V) then
               abc(1) = -vect(4)+dir2(1)
               abc(2) = -vect(5)+dir2(2)
               abc(3) = -vect(6)+dir2(3)
               norm = sqrt(abc(1)**2+abc(2)**2+abc(3)**2)
               abc(1) = abc(1)/norm
               abc(2) = abc(2)/norm
               abc(3) = abc(3)/norm
               pol = polar(1)*abc(1) + polar(2)*abc(2) + polar(3)*abc(3)
               xphot(8,1)  = -polar(1) + 2.*pol*abc(1)
               xphot(9,1)  = -polar(2) + 2.*pol*abc(2)
               xphot(10,1) = -polar(3) + 2.*pol*abc(3)
               norm = sqrt(xphot(8,1)**2 + xphot(9,1)**2 + xphot(10,1)**2)
               xphot(8,1)  = xphot(8,1)/norm
               xphot(9,1)  = xphot(9,1)/norm
               xphot(10,1) = xphot(10,1)/norm
            endif
cccccc
            XPHOT(11,1) = TOFG+step/clight*refwtr

            if (OUTPUT_SCATTABLE.eq.1) then
               CALL GSKPHO(-ISTORY-20000000) ! track photon scattering and reflections
            else
               CALL GSKPHO(0)
            endif

            NGPHOT=0
     
*  If p.e. is produced(3)
         ELSE IF(IFLGREF.eq.3) THEN

            if (OUTPUT_SCATTABLE.eq.1) then
               call addtoscattable(nopmt,ISTORY)
            endif

            NOPH(ITRA)=NOPH(ITRA)+1
            TOFGSK = TOFG * (1.E+09)
            CALL GRNDM(RDUMMY,1)
            if((SK_GEOMETRY .lt.3).and.(WATER_TUNE .eq. UNI_TUNE))then
               if(RDUMMY.lt.0.0115) then
                  call rngausk2(107.,6.,del)
                  tofgsk=tofgsk+del
               elseif(RDUMMY.lt.0.019) then
c     call rngausk2(37.5,2.,del)
                  tofgsk=tofgsk+36.5
               elseif(RDUMMY.lt.0.04275) then
                  call rngausk2(87.,22.,del)
                  tofgsk=tofgsk+del
               elseif(RDUMMY.lt.0.043) then
                  call rngausk2(-47.,8.,del)
                  tofgsk=tofgsk+del
               elseif(RDUMMY.lt.0.05) then
                  call rngausk(15.,4.5,del)
                  tofgsk=tofgsk+del
               endif
            else
               if(RDUMMY.lt.0.014) then
                  call rngausk2( 111.5,7.0,del)
                  tofgsk=tofgsk+del
               elseif(RDUMMY.lt.0.027) then
                  call rngausk2(39.,5.0,del)
                  tofgsk=tofgsk+del
               elseif(RDUMMY.lt.0.05075) then
                  call rngausk2(80.,25.,del)
                  tofgsk=tofgsk+del
               elseif(RDUMMY.lt.0.051) then
                  call rngausk2(-47.,8.,del)
                  tofgsk=tofgsk+del
               elseif(RDUMMY.lt.0.061) then
                  call rngausk(15.,4.5,del)
                  tofgsk=tofgsk+del
               endif
            endif
            if ( 
     &  (TRACK_PHOTON .eq. TRACK_ALL_PHOTON) .or.
     &  (TRACK_PHOTON .eq. TRACK_QE_OVER_ABS_PHOTON) 
     &  ) then
               ISEGM = ITRA
               IDATA(1) = NOPMT
               RDATA(2) = TOFGSK
               NDATA = 2
               CALL KZREP2('MCPHOTON',ISEGM,NOPH(ISEGM),NDATA,IDATA)
            ENDIF

           if (
     &         (TRACK_PHOTON .eq. TRACK_ALL_PHOTON) .or.
     &         (TRACK_PHOTON .eq. TRACK_QE_PHOTON)
     &         ) then

            IF(NQGESK(NOPMT,ITRA).EQ.0) THEN
               TGESK(NOPMT,ITRA)=TOFGSK
            ELSE IF(TOFGSK.LE.TGESK(NOPMT,ITRA)) THEN
                    TGESK(NOPMT,ITRA)=TOFGSK
            ENDIF

            IF (TOFGSK.le.TGESK(NOPMT,ITRA)+SGGATE) THEN
               NQGESK(NOPMT,ITRA)=NQGESK(NOPMT,ITRA)+1
            ENDIF

           elseif (
     &         (TRACK_PHOTON .eq. TRACK_QE_OVER_ABS_PHOTON) .or.
     &         (TRACK_PHOTON .eq. TRACK_QE_OVER_ABS_SAVE_MEMORY)
     &         ) then
              call sgcnt_photon(NOPMT)
           else
              write(*,*) 'WARINING: SGPMT: TRACK_PHOTON flag is irregular !!'
           endif            

C
C+ Standard Diagnostic Test Histograms
C
            if (OUTPUT_HBOOK .eq. OUT_TEST) then
               TEST_Coll_phot=TEST_Coll_phot+1
               TEST_Hitmap(NOPMT)=TEST_Hitmap(NOPMT)+1.
            endif

* diffused reflection(2) is not supported yet. And another is strange
*         ELSE IF( IFLGREF.ne.2 ) THEN
         ELSE
               WRITE(6,*) 'Warning in SGPMT not suported IFLGREF:',IFLGREF
         ENDIF

      ELSE      ! SK_I, SK_II, SK_III&&TRACK_QE_PHOTON

******
* Original apdetsim order (reflection -> correction eff.)
******

        if (REFLECTION_TUNE.eq.AP_TUNE) then

               BLAMDA = 1237.8E-09 / GEKIN

*  Calculate abs. or reflection

*    costh is incident angle from skdonuts : cos2th = (vect. dir2)
         costh=-1.*(vect(4)*dir2(1)+vect(5)*dir2(2)+vect(6)*dir2(3))   
** 97/3/10 temporary modified by Y.Itow
         if(costh.lt.-1)then
            costh = costh + 0.000001
         elseif(costh.gt.1)then
            costh = costh - 0.000001
         endif
**      
         costh=sqrt((costh+1.)/2.)
         CALL SGREFPMT(NOPMT,VECT,QE,VECNEW,COSTH,IFLGREF)
C        CALL SGREFPMT(NOPMT,VECT,IFLGREF,VECNEW,COSTH)
CKO      CALL SGREFPMT(NOPMT,VECT,IFLGREF,VECNEW)
CMS      CALL SGREFPMT(NOPMT,VECT,IFLGREF,VECNEW,POLARNEW)

c      write(*,*) 'Ref or Caught ?=',IFLGREF

*  If photon is reflected by PMT 

         IF(IFLGREF.EQ.1.OR.IFLGREF.EQ.2) THEN

cc            CALL GNEXT(VECT,SNEXT,SAFETY)
cc            STEP = SNEXT + PREC
            STEP = FIXSTEP
         
            NGPHOT = 1
            DO I = 1,3
**** use position and direction calculated by skdonuts
C            XPHOT(I+3,1) = VECNEW(I)
C            XPHOT(I,1) = VECT(I) + STEP*VECNEW(I)
               if (IFLGREF.EQ.1) then
                  XPHOT(I+3,1) = dir2(i)
                  XPHOT(I,1) = pos3(i) + STEP*dir2(i)
               else
                  XPHOT(I+3,1) = dirran(i)
                  XPHOT(I,1) = posran(i) + STEP*dirran(i)
               endif
****  porlarization is not yet suported
C     XPHOT(I+7,1) = POLARNEW(I)
            ENDDO
         
            XPHOT(7,1)  = VECT(7)
c     XPHOT(11,1) = TOFG + STEP/CLIGHT*refwtr
            XPHOT(11,1) = TOFG

            if (OUTPUT_SCATTABLE.eq.1) then
               CALL GSKPHO(-ISTORY-20000000) ! track photon scattering and reflections
            else
               CALL GSKPHO(0)
            endif

            NGPHOT=0

*     Another flag is not supported yet
         ELSE if(IFLGREF.ne.0) then
            WRITE(6,*) 'Warning in SGPMT'

         Else

* Q.E.
            IF(TRACK_PHOTON.eq.TRACK_ALL_PHOTON) then
               CALL GRNDM(RDUMMY,1)
C               IF(RDUMMY.GT.(PMTQESG(BLAMDA)/(1.0-0.25))) THEN
               IF(RDUMMY.GT.(PMTQESG(BLAMDA)*REDUCFCT/(1.0-0.25))) THEN
                  RETURN
               ENDIF
            ELSEIF ((TRACK_PHOTON.eq.TRACK_QE_OVER_ABS_PHOTON) .or. 
     &              (TRACK_PHOTON.eq.TRACK_QE_OVER_ABS_SAVE_MEMORY)) then
               CALL GRNDM(RDUMMY,1)
               IF(RDUMMY .GT. (PRABSL0(BLAMDA)/0.75)) then
                  RETURN
               ENDIF
            ENDIF

*  Correction eff.
*----------------------------------------
* apply QE table (only for SK-II or SK-III)
*----------------------------------------

            if(icall.eq.0) then
               icall = 1

*=== SK-III or SK-IV ===

               if(SK_GEOMETRY.ge.SK_III) then
                  ipmt=1
CCCC moved to dsinip.F 2008.04.24 Byeongsu Yang
CC                  maxfct=1.3

                  fqetab='qetable'//cqetab(1:3)//'.dat'
                  call findconsts(fname,'QE table',
c     $                 'qetable3_0.dat',0)
     $                 fqetab,0)
                  open(unit=91,file=fname,status='old',
     $                 form='formatted')
 310              continue
                  read(91,*,end=311) icab,qetable
                  factor(icab)=qetable/maxfct
                  IF(factor(icab).gt.1.0) factor(icab)=1.0
                  ipmt=ipmt+1
                  goto 310

 311              continue
                  close(91)
                  if(ipmt.ne.maxpm+1) then
                     write(6,*) 'Entry of QE table is wrong !'
                     stop
                  endif
CCCC moved to dsinip.F 2008.04.24 Byeongsu Yang
CC* tuned by cosmic ray muon
CCC                  COREPMT = COREPMT*0.929
CCC     SK-III charge scale is tuned by G.Mitsuka(Apr11,2007)
CCC                  COREPMT = COREPMT*0.978
CCC     SK-III charge scale is re-tuned after BLST bug-fix by G.Mitsuka(Oct,2007)
CC                  COREPMT = COREPMT*0.989

                  do ipmt=1,maxpm
                     collpmt(ipmt)=COREPMT*factor(ipmt)
                  enddo
                  

*=== SK-II ===
              else if(SK_GEOMETRY.eq.SK_II) then
                  ipmt=1
CCC moved to dsinip.F 2008.04.24 Byeongsu Yang
C                  maxfct=1.3

                  call findconsts(fname,'QE table',
     $                 'qetable_1.dat',0)
                  open(unit=91,file=fname,status='old',
     $                 form='formatted')
 110              continue
                  read(91,*,end=111) icab,qetable
                  factor(icab)=qetable/maxfct
                  IF(factor(icab).gt.1.0) factor(icab)=1.0
                  ipmt=ipmt+1
                  goto 110

 111              continue
                  close(91)
                  if(ipmt.ne.maxpm+1) then
                     write(6,*) 'Entry of QE table is wrong !'
                     stop
                  endif
CCC moved to dsinip.F 2008.04.24 Byeongsu Yang
* tuned by cosmic ray muon
C                  COREPMT = COREPMT*0.929
                  do ipmt=1,maxpm
                     collpmt(ipmt)=COREPMT*factor(ipmt)
                  enddo
*=== SK-I ===
               else if(SK_GEOMETRY.eq.SK_I) then
                  do ipmt=1,maxpm
                     collpmt(ipmt)=COREPMT
                  enddo
               endif
            endif

* apply QE table (only for SK-II or SK-III)
*to here----------------------------------------

            CALL GRNDM(RDUMMY,1)
C            IF(RDUMMY.LE.COREPMT) THEN
            IF(RDUMMY.LE.collpmt(NOPMT)) THEN

         
*     If photon is caught by PMT
            
               NOPH(ITRA)=NOPH(ITRA)+1
            
               TOFGSK = TOFG * (1.E+09)

               call rngausk(100.,10.,del)               
C     if(rndm().lt.0.02) tofgsk=tofgsk+del
               CALL GRNDM(RDUMMY,1)
               if(RDUMMY.lt.0.02) tofgsk=tofgsk+del
            
C--   CHANGE-Y.H/96/11/16            IF(TOFGSK.LE.SGGATE) THEN
C     IF(NQGESK(NOPMT,ITRA).EQ.0)THEN
C     TGESK(NOPMT,ITRA)=TOFGSK
C     ELSE
C     IF(TOFGSK.LE.TGESK(NOPMT,ITRA)) 
C     &              TGESK(NOPMT,ITRA)=TOFGSK
C     ENDIF
C     NQGESK(NOPMT,ITRA)=NQGESK(NOPMT,ITRA)+1
C     
C--   CHANGE-Y.H/96/11/16            ENDIF

            if ( 
     &  (TRACK_PHOTON .eq. TRACK_ALL_PHOTON) .or.
     &  (TRACK_PHOTON .eq. TRACK_QE_OVER_ABS_PHOTON) 
     &  ) then
                  ISEGM = ITRA
                  IDATA(1) = NOPMT
                  RDATA(2) = TOFGSK
                  NDATA = 2
                  CALL KZREP2('MCPHOTON',ISEGM,NOPH(ISEGM),NDATA,IDATA)
            endif

           if ((TRACK_PHOTON .eq. TRACK_ALL_PHOTON) .or.
     &         (TRACK_PHOTON .eq. TRACK_QE_PHOTON)
     &         ) then
               IF(NQGESK(NOPMT,ITRA).EQ.0)THEN
                  TGESK(NOPMT,ITRA)=TOFGSK
               ELSE
                  IF(TOFGSK.LE.TGESK(NOPMT,ITRA)) 
     &                 TGESK(NOPMT,ITRA)=TOFGSK
               ENDIF
            
               IF (TOFGSK.le.TGESK(NOPMT,ITRA)+SGGATE) THEN
                  NQGESK(NOPMT,ITRA)=NQGESK(NOPMT,ITRA)+1
               ENDIF
            
           elseif (
     &         (TRACK_PHOTON .eq. TRACK_QE_OVER_ABS_PHOTON) .or.
     &         (TRACK_PHOTON .eq. TRACK_QE_OVER_ABS_SAVE_MEMORY)
     &         ) then
              call sgcnt_photon(NOPMT)
           else
              write(*,*) 'WARINING: SGPMT: TRACK_PHOTON flag is irregular !!'
           endif
C
C+ Standard Diagnostic Test Histograms
C
               if (OUTPUT_HBOOK .eq. OUT_TEST) then
                  TEST_Coll_phot=TEST_Coll_phot+1
                  TEST_Hitmap(NOPMT)=TEST_Hitmap(NOPMT)+1.
               endif

            ENDIF     ! collpmt(NOPMT) condition
         ENDIF        ! IFLGREF condition

******
* Original ledetsim order (QE -> correction eff.-> reflection)
******

        else if (REFLECTION_TUNE.eq.LE_TUNE) then

*  Q.E.

         BLAMDA = 1237.8E-09 / GEKIN
c         IF(TRACK_PHOTON.eq.TRACK_ALL_PHOTON) then
c            qe=pmtqesg(blamda)
            IF(SK_GEOMETRY.eq.SK_I) then !Good Old SK-I

             if (TRACK_PHOTON .eq. TRACK_ALL_PHOTON) then
              QE = pmtqesg(blamda) 
             elseif (TRACK_PHOTON .eq. TRACK_QE_PHOTON) then
              QE = 0.75/REDUCFCT
             elseif (TRACK_PHOTON .eq. TRACK_QE_OVER_ABS_PHOTON) then
              QE = PRABSL0(blamda)/REDUCFCT
             elseif (TRACK_PHOTON .eq. TRACK_QE_OVER_ABS_SAVE_MEMORY) then
              QE = PRABSL0(blamda)/REDUCFCT
             endif

             if(ikeksk(nopmt).eq.1) qe=qe*1.2

*----------------------------------------
* apply QE table (only for SK-II or SK-III)
*----------------------------------------

C*=== SK-III ===
C            else if(SK_GEOMETRY.eq.SK_III) then
C               if(icall.eq.0) then
C                  icall = 1
C                  ipmt=1
C
C                  call findconsts(fname,'QE table',
C     $                 'qetable3_0.dat',0)
C                  open(unit=91,file=fname,status='old',
C     $                 form='formatted')
C 410              continue
C                  read(91,*,end=411) icab, qetable
C                  factor(ipmt)=qetable
C                  ipmt=ipmt+1
C                  goto 410
C
C 411              continue
C                  close(91)
C                  if(ipmt.ne.maxpm+1) then
C                     write(6,*) 'Entry of QE table is wrong !'
C                     stop
C                  endif
C               endif
C
C               QE = QE * factor(NOPMT)
C
C               GOTO 3000	! do PMT simulation(p.e. reflection and so on.)
*=== SK-II ===
            else if(SK_GEOMETRY.eq.SK_II) then
               if(icall.eq.0) then
                  icall = 1
                  ipmt=1

                  call findconsts(fname,'QE table',
     $                 'qetable_1.dat',0)
                  open(unit=91,file=fname,status='old',
     $                 form='formatted')
 210              continue
                  read(91,*,end=211) icab, qetable
                  factor(ipmt)=qetable
                  ipmt=ipmt+1
                  goto 210

 211              continue
                  close(91)
                  if(ipmt.ne.maxpm+1) then
                     write(6,*) 'Entry of QE table is wrong !'
                     stop
                  endif
               endif

c               QE = QE * factor(NOPMT)

               if (TRACK_PHOTON .eq. TRACK_ALL_PHOTON) then
                QE = pmtqesg(blamda) * factor(NOPMT)
               elseif (TRACK_PHOTON .eq. TRACK_QE_PHOTON) then
                QE = factor(NOPMT)*0.75/REDUCFCT
               elseif (TRACK_PHOTON .eq. TRACK_QE_OVER_ABS_PHOTON) then
                QE = factor(NOPMT)*PRABSL0(blamda)/REDUCFCT
               elseif (TRACK_PHOTON .eq. TRACK_QE_OVER_ABS_SAVE_MEMORY) then
                QE = factor(NOPMT)*PRABSL0(blamda)/REDUCFCT
               endif
*----------------------------------------

            Endif
            if(qe.gt.1.) write(*,*) 'error in sgpmt',qe

cc         ELSE IF(TRACK_PHOTON.eq.TRACK_QE_PHOTON) THEN !Q.E. cut already made
cc            QE=2. !always pass this Q.E. cut
cc         ENDIF

         CALL GRNDM(RDUMMY,1)
         IF(RDUMMY.LE.QE) THEN

            CALL GRNDM(RDUMMY,1)
            IF(RDUMMY.LE.COREPMT) THEN

*  If photon is caught by PMT

               NOPH(ITRA)=NOPH(ITRA)+1

               TOFGSK = TOFG * (1.E+09)

               CALL GRNDM(RDUMMY,1)
               if(RDUMMY.lt.0.0115) then
                  call rngausk2(107.,6.,del)               
                  tofgsk=tofgsk+del
               elseif(RDUMMY.lt.0.019) then
c                  call rngausk2(37.5,2.,del)               
                  tofgsk=tofgsk+36.5
               elseif(RDUMMY.lt.0.04275) then
                  call rngausk2(87.,22.,del)               
                  tofgsk=tofgsk+del
               elseif(RDUMMY.lt.0.043) then
                  call rngausk2(-47.,8.,del)               
                  tofgsk=tofgsk+del
               elseif(RDUMMY.lt.0.05) then
                  call rngausk(15.,4.5,del)
                  tofgsk=tofgsk+del
               endif

            if ( 
     &  (TRACK_PHOTON .eq. TRACK_ALL_PHOTON) .or.
     &  (TRACK_PHOTON .eq. TRACK_QE_OVER_ABS_PHOTON) 
     &  ) then
                  ISEGM = ITRA
                  IDATA(1) = NOPMT
                  RDATA(2) = TOFGSK
                  NDATA = 2
                  CALL KZREP2('MCPHOTON',ISEGM,NOPH(ISEGM),NDATA,IDATA)
            endif

           if ((TRACK_PHOTON .eq. TRACK_ALL_PHOTON) .or.
     &         (TRACK_PHOTON .eq. TRACK_QE_PHOTON)
     &         ) then

               IF(NQGESK(NOPMT,ITRA).EQ.0)THEN
                  TGESK(NOPMT,ITRA)=TOFGSK
               ELSE
                  IF(TOFGSK.LE.TGESK(NOPMT,ITRA))
     &                 TGESK(NOPMT,ITRA)=TOFGSK
               ENDIF
            
               IF (TOFGSK.le.TGESK(NOPMT,ITRA)+SGGATE) THEN
                  NQGESK(NOPMT,ITRA)=NQGESK(NOPMT,ITRA)+1
               ENDIF
           elseif (
     &         (TRACK_PHOTON .eq. TRACK_QE_OVER_ABS_PHOTON) .or.
     &         (TRACK_PHOTON .eq. TRACK_QE_OVER_ABS_SAVE_MEMORY)
     &         ) then
              call sgcnt_photon(NOPMT)
           else
              write(*,*) 'WARINING: SGPMT: TRACK_PHOTON flag is irregular !!'
           endif
C
C+ Standard Diagnostic Test Histograms
C
               if (OUTPUT_HBOOK .eq. OUT_TEST) then
                  TEST_Coll_phot=TEST_Coll_phot+1
                  TEST_Hitmap(NOPMT)=TEST_Hitmap(NOPMT)+1.
               endif

            ENDIF

         ELSE   ! QE condition

*  Calculate abs. or reflection

            CALL SGREFPMT(NOPMT,VECT,QE,dir2,costh,IFLGREF)

*  If photon is reflected by PMT 
            IF(IFLGREF.EQ.1.OR.IFLGREF.EQ.2) THEN

               NGPHOT = 1
               step=sqrt((pos2(1)-pos3(1))**2+
     &              (pos2(2)-pos3(2))**2+
     &              (pos2(3)-pos3(3))**2)+
     &              sqrt((pos2(1)-vect(1))**2+
     &              (pos2(2)-vect(2))**2+
     &              (pos2(3)-vect(3))**2)+fixstep
               norm = sqrt(dir2(1)**2+dir2(2)**2+dir2(3)**2)
               if(norm.gt.0)then
                  dir2(1) = dir2(1)/norm
                  dir2(2) = dir2(2)/norm
                  dir2(3) = dir2(3)/norm
               else
                  write(*,*)'ERROR in sgpmt dir2 = 0'
                  stop
               end if
               DO I = 1,3
                  XPHOT(I+3,1) = dir2(I)
                  XPHOT(I,1) = pos3(I) + fixstep*dir2(I)
               ENDDO

               refwtr=effnsg(blamda)
               xphot(7,1)=VECT(7)
               XPHOT(11,1) = TOFG+step/clight*refwtr

**** porlarization is not yet

               if (OUTPUT_SCATTABLE.eq.1) then
                  CALL GSKPHO(-ISTORY-20000000) ! track photon scattering and reflections
               else
                  CALL GSKPHO(0)
               endif

               NGPHOT=0

*  Another flag is not supported yet
            ELSE if(IFLGREF.ne.0) then
               WRITE(6,*) 'Warning in SGPMT'
            ENDIF

         ENDIF  ! QE condition
        endif   ! REFLECTION_TUNE condition

      Endif     ! SK_GEOMETRY condition

c      print*,'pmtout/ca--tch',gekin

#ifdef FADC
      if (IFADCI) then
c     WRITE(6,*)'IFADCI = ',IFADCI,' --> Call fadcrecord '
         call fadcrecord(NOPMT,TOFGSK)
      endif
#endif

      RETURN
      END
