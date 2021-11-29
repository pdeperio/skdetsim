************************************************************************
*
*     -----------------
      SUBROUTINE GUTREV
*     -----------------
*
*     (Purpose)
*       Steering routine to trace all tracks in the current event.
*
*     (Parameter)
*       NPC   : Particle Code
*       NOPH  : Number of Photon for 1 event
*       NCBL  : Cable Number
*       TOFP  : Hit Time
*
*     (Creation Date and Author)
*       1992.11.30 : First version by Y.Koshio
*
*       1998.03.16 ; Y.Hayato : add new secondary particle storage code
*       2011.08.29 ; P.de Perio: add NEUT cascade history bank
*
************************************************************************
      implicit none
#include "geant321/gckine.inc"
#include "skmc.h"
#include "skhead.h"
#include "skparm.h"
#include "gupht.h"

#include "geant321/gcflag.inc"
#include "guanti.h"

#include "guflag.h"
#include "dsbnkf.h"
#include "guhigh.h"

#include "dsdata.h"
#include "dsmulti.h"
#include "skdetsim_parm.h"

*--  temporary common
#include "apmultmp.h"

*** 1998/03/18 Add Y.H. new common for secondary particle
#include "vcwork.h"
#include "apscndry.h"

*** 2011/08/29 P. de Perio - Common block of NEUT cascade option 
#ifdef SITE_HAS_NEUT
#include "apcard.h"
#include "fsihist.h"

      real*4    EVPIPROB
      external  EVPIPROB
#endif

#include "aptrkseg.h"
      
#ifdef FADC
#include "fadcwrbnk.h"
#endif
      integer idata,ibuf,nwmast,nw,ierr,nsegm,lsegm,isegm,isseg,ndata,i
      integer iseq,iseed1,iseed2,j,numanti
      real rdata
      dimension idata(27),rdata(27)

      DIMENSION IBUF(2*MAXRNG+2)
      EQUIVALENCE (IDATA(1),RDATA(1))


      nevt = nevt + 1

*--  initialize 

      call vzero(noph,maxrng)
      
*--  create banks

*--  MCPHOTON

      if ((TRACK_PHOTON .ne. NOT_TRACK_PHOTON).and.(OUTPUT_ZBS)) then
*      IF(IFHIGH.EQ.0) THEN
c      IF(TRACK_PHOTON .eq. TRACK_ALL_PHOTON) THEN
      IF( 
     & (TRACK_PHOTON .eq. TRACK_ALL_PHOTON) .or.
     & (TRACK_PHOTON .eq. TRACK_QE_OVER_ABS_PHOTON) 
     & ) THEN

       if ((TRACK_PHOTON .eq. TRACK_QE_OVER_ABS_PHOTON) .or. 
     & ((TRACK_PHOTON .eq. TRACK_ALL_PHOTON).and.(SK_GEOMETRY.ge.SK_IV))
     &)then
         NWMAST=(NRING+1)*2+3  ! dsnois_multi doesn't use MCTQ bank.
         NW=2
         CALL KZBCR2('MCPHOTON',NWMAST,'3I',NRING+1,NW,'1I 1F',IERR)
       else
         NWMAST=NRING*2+3
         NW=2
         CALL KZBCR2('MCPHOTON',NWMAST,'3I',NRING,NW,'1I 1F',IERR)
       endif

      ENDIF
      

*--  MCRAWTQ

      NWMAST=NRING*2+1
      NW=3
      if( sk_geometry .lt. sk_iv )then
        CALL KZBCR2('MCRAWTQ',NWMAST,'3I',NRING,NW,'1I 1F 1I',IERR)
      else
        CALL KZBCR2('MCRAWTQ',3,'3I',1,NW,'1I 1F 1I',IERR)
      endif

*--  MCRAWTQA
      NWMAST=NRING*2+1
      NW=3
      if( sk_geometry .lt. sk_iv )then
        CALL KZBCR2('MCRAWTQA',NWMAST,'3I',NRING,NW,'1I 1F 1I',IERR)
      else if( sk_geometry .ge. sk_iv .and.  od_photon .eq. count_all_od_photon )then 
        CALL KZBCR2('MCPHOTAN',3,'3I',1,2 ,'1I 1F'   ,IERR)
      else
        CALL KZBCR2('MCRAWTQA',3,'3I',1,NW,'1I 1F 1I',IERR)
      endif

      endif ! end of track_photon check


*--  CONVVECT (temporary)
      flag_warning_convvect = 0
      npi0=0
      nwmast=1
      nsegm =1
      lsegm =11
      CALL KZBCR1('CONVVECT',nwmast,'1I',nsegm,lsegm,'2I 7F 2I',IERR)

*** 1998/03/18 Add Y.H. new common for secondary particle
      nscndprt=0
      nwmast = 1
      nsegm  = 1
      lsegm  = 27
      CALL KZBCR1('SCNDPRT',nwmast,'1I',nsegm,lsegm,'3I 16F 8I',IERR)

cc   control flag for warning of too many scondary particles.
      flag_warning_scndprt = 0

      do i=1,secmaxrng
         itrkscnd(i)=0
         istakscnd(i)=0
         iprtscnd(i)=0
         do j=1,3
            vtxscnd(j,i)=0
            pscnd(j,i)=0
            pprnt(j,i)=0
            pprntinit(j,i)=0
            vtxprnt(j,i)=0
         end do
         tscnd(i)=0
         iprntprt(i)=0
         lmecscnd(i)=0
         iprnttrk(i)=0
         iorgprt(i)=0
         iflgscnd(i)=0
         iprntidx(i)=0
         nchilds(i)=0
         ichildidx(i)=0
      end do

c     New bank for linking input particles (MCVECT) to secondaries stack
C     nvc = shared from vcwork.h
      nwmast = 1
      nsegm  = 1
      lsegm  = 2
      CALL KZBCR1('SCNDPRTVC',nwmast,'1I',nsegm,lsegm,'2I',IERR)

      npvcscnd=nvc
      do i=1,MAXVC
         nchildsvc(i) = 0
         ichildidxvc(i) = 0         
      end do

*** 2011/08/29 P. de Perio - Common block of NEUT cascade option 
C     for hadronic interactions
#ifdef SITE_HAS_NEUT
      if (IAPIMOD.eq.2) then
         call EFCLFSI
      end if
#endif

*---- Traces all tracks

      CALL GTREVE
      

*---- Fill the data into bank

ccc      print*,'IEVENT=',ievent,' ; NOPH=',noph(1),noph(2),noph(3),
ccc     &               noph(4),noph(5),noph(6),noph(7)

*     INNER

      if ((TRACK_PHOTON .ne. NOT_TRACK_PHOTON) .and. (OUTPUT_ZBS)) then
*      IF(IFHIGH.EQ.0)THEN
c      IF(TRACK_PHOTON .eq. TRACK_ALL_PHOTON) THEN

      IF (
     & (TRACK_PHOTON .eq. TRACK_ALL_PHOTON) .or.
     & (TRACK_PHOTON .eq. TRACK_QE_OVER_ABS_PHOTON)
     &  ) THEN

         ISEGM = 0
         ISSEG = 0
c         NDATA = 2*NRING+3
         IBUF(1)=NRING	
         DO 200 I=1,NRING
            IBUF(2*I)   = IPRT(I)
            IBUF(2*I+1) = NOPH(I)
            NOPH(I)=0
 200     CONTINUE

         iseq = -1
         call grndmq(iseed1,iseed2,iseq,'G')

       if ((TRACK_PHOTON .eq. TRACK_QE_OVER_ABS_PHOTON) .or. 
     & ((TRACK_PHOTON .eq. TRACK_ALL_PHOTON).and.(SK_GEOMETRY.ge.SK_IV))
     &)then
         NDATA = 2*(NRING+1)+3
         IBUF(2+2*(NRING+1))=ISEED1
         IBUF(3+2*(NRING+1))=ISEED2
       else
         NDATA = 2*NRING+3
         IBUF(2+2*NRING)=ISEED1
         IBUF(3+2*NRING)=ISEED2
       endif

ccc         CALL RDMOUT(ISEED)
ccc         IBUF(2+2*NRING)=ISEED
         CALL KZREP2('MCPHOTON',ISEGM,ISSEG,NDATA,IBUF)
      ENDIF

      IF (
     & (TRACK_PHOTON .eq. TRACK_ALL_PHOTON) .or.
     & (TRACK_PHOTON .eq. TRACK_QE_PHOTON)
     &  ) THEN

      DO 300 J=1,NRING
         ISEGM = J
         ISSEG = 0
         NOPH(J)=0
         DO 400 I=1,MAXPM
            IF(NQGESK(I,J).NE.0)THEN
               IDATA(1) = I
               RDATA(2) = TGESK(I,J)
               IDATA(3) = NQGESK(I,J)
               TGESK(I,J) = 0.
               NQGESK(I,J)= 0
               NDATA = 3
               ISSEG=ISSEG+1
               CALL KZREP2('MCRAWTQ',ISEGM,ISSEG,NDATA,IDATA)
               NOPH(J)=NOPH(J)+1
            ENDIF
 400     CONTINUE
 300  CONTINUE

      ISEGM = 0
      ISSEG = 0
      NDATA = 1+2*NRING
      IBUF(1)=NRING	
      DO 500 I=1,NRING
         IBUF(2*I)   = IPRT(I)
         IBUF(2*I+1) = NOPH(I)
         NOPH(I)=0
 500  CONTINUE
      CALL KZREP2('MCRAWTQ',ISEGM,ISSEG,NDATA,IBUF)

      else if (
     & (TRACK_PHOTON .eq. TRACK_QE_OVER_ABS_PHOTON) .or.
     & (TRACK_PHOTON .eq. TRACK_QE_OVER_ABS_SAVE_MEMORY)
     &  ) THEN

      DO 310 J=1,1
         ISEGM = J
         ISSEG = 0
         NOPH(J)=0
         DO 410 I=1,num_gate_all
            IF(NPHOTO_GATE(I).NE.0)THEN
               IDATA(1) = Pmt_id(I)
               RDATA(2) = Timing_gate(I)
               IDATA(3) = NPHOTO_GATE(I)
               NDATA = 3
               ISSEG=ISSEG+1
               CALL KZREP2('MCRAWTQ',ISEGM,ISSEG,NDATA,IDATA)
c       write(88,*) ievent,pmt_id(i),timing_gate(i),nphoto_gate(i),width_photons(i)
               NOPH(J)=NOPH(J)+1
            ENDIF
 410     CONTINUE
 310  CONTINUE

      ISEGM = 0
      ISSEG = 0
      NDATA = 1+2*1
      IBUF(1)= 1
      DO 510 I=1,1
         IBUF(2*I)   = IPRT(I)
         IBUF(2*I+1) = NOPH(I)
         NOPH(I)=0
 510  CONTINUE
      CALL KZREP2('MCRAWTQ',ISEGM,ISSEG,NDATA,IBUF)

      call vzero(pmt_id,MXGATE)
      call vzero(nphoto_gate,MXGATE)
      call vzero(timing_gate,MXGATE)
      call vzero(width_photons,MXGATE)
      call vzero(list_gate,MAXPM*50)
      call vzero(num_gate,MAXPM)
      call vzero(noph,MAXRNG)

      num_gate_all = 0

      endif

      endif

*    CONVVECT (temporary)

      idata(1)=npi0
      ndata   = 1
      isegm   = 0
      CALL KZREP1('CONVVECT',isegm,ndata,idata)
c      print*,npi0,(ippi0(i),i=1,npi0)
c      do i=1,npi0
c         print*,(vpi0(j,i),j=1,3),(ppi0(j,i),j=1,3)
c      enddo
      do 501 i=1,npi0
         idata(1)=itrcon(i)
         idata(2)=iprtpi0(i)
         do 502 j=1,3
            rdata(j+2)=vpi0(j,i)
            rdata(j+5)=ppi0(j,i)
 502     continue
         rdata(9)=tpi0(i)
C-ADD Y.H 96/11/18
         idata(10)=ippi0(i)
         idata(11)=imeca(i)
         lsegm    = 11
         CALL KZREP1('CONVVECT',i,lsegm,idata)
 501  continue

*** 1998/03/18 Add Y.H. new common for secondary particle
*    SCNDPRT

      idata(1)=nscndprt
      ndata   = 1
      isegm   = 0
      CALL KZREP1('SCNDPRT',isegm,ndata,idata)

      do i=1,nscndprt
         idata(1)=itrkscnd(i)
         idata(2)=istakscnd(i)
         idata(3)=iprtscnd(i)
         do j=1,3
            rdata(j+3)=vtxscnd(j,i)
            rdata(j+6)=pscnd(j,i)
            rdata(j+9)=pprnt(j,i)            
            rdata(j+12)=pprntinit(j,i)            
            rdata(j+15)=vtxprnt(j,i)            
         enddo
         rdata(19)=tscnd(i)
         idata(20)=iprntprt(i)
         idata(21)=lmecscnd(i)
         idata(22)=iprnttrk(i)
         idata(23)=iorgprt(i)
         idata(24)=iflgscnd(i)
         idata(25)=iprntidx(i)
         idata(26)=nchilds(i)
         idata(27)=ichildidx(i)

         lsegm    = 27
         CALL KZREP1('SCNDPRT',i,lsegm,idata)
      enddo

*     SCNDPRTVC
      idata(1)=nvc
      ndata   = 1
      isegm   = 0
      CALL KZREP1('SCNDPRTVC',isegm,ndata,idata)

      do i=1,nvc
         idata(1)=nchildsvc(i)
         idata(2)=ichildidxvc(i)
         lsegm   = 2
         CALL KZREP1('SCNDPRTVC',i,lsegm,idata)
      enddo

*** 2011/08/29 P. de Perio - Common block of NEUT cascade option 
C     for hadronic interactions
#ifdef SITE_HAS_NEUT
      if (IAPIMOD.eq.2) then

C -- Calculate FSIPROB after all possible instances of EFTRACE
C     (This is for the NEUT cascade reweighting)
         FSIPROB=EVPIPROB()
    
         IDATA(1)=NVERT

         nwmast = 1
         nsegm  = NVERT
         lsegm  = 4
         CALL KZBCR1('NEUTVERTSK',nwmast,'1I',nsegm,lsegm,'3F 1I',IERR)

         CALL KZREP1('NEUTVERTSK',0,nwmast,idata(1))

C         write (*,*) "NVERT = ", idata(1)
         
         do i=1,nvert

C            write(*,*) i,  " posvert = "

            do j=1,3
               rdata(j)=posvert(j,i)
C               write(*,*) rdata(j), " "
            end do
            idata(4)=iflgvert(i)

C            write(*,*) "iflgvert = ", idata(4)

            CALL KZREP1('NEUTVERTSK',i,lsegm,rdata)
         end do
         

         idata(1)= nvcvert
         rdata(2)= fsiprob

         nwmast = 2
         nsegm  = nvcvert
         lsegm  = 8

         CALL KZBCR1('NEUTPARTSK',nwmast,'1I 1F',nsegm,lsegm,'5F 3I',IERR)
         CALL KZREP1('NEUTPARTSK',0,nwmast,idata(1))

C         write(*,*) "Nvcvert = ",idata(1) , " fsiprob = ", rdata(2)

         do i=1,nvcvert
            
C            write(*,*) i,  " dirvert = "

            do j=1,3
               rdata(j)=dirvert(j,i)
C               write(*,*) rdata(j), " "
           end do
            rdata(4)=abspvert(i)
            rdata(5)=abstpvert(i)
            idata(6)=ipvert(i)
            idata(7)=iverti(i)
            idata(8)=ivertf(i)

C            write(*,*) rdata(4),rdata(5),idata(6),idata(7),idata(8)

            CALL KZREP1('NEUTPARTSK',i,lsegm,rdata)
         end do

      end if
#endif

!     ANTI
! ---------------------  Start of MCRAWTQA creation

      if ((TRACK_PHOTON .ne. NOT_TRACK_PHOTON).and.(OUTPUT_ZBS)) then
      ! atm behavior is simulated for SK-I,II,III
      if (SK_GEOMETRY .le. SK_III) then
      DO 900 J=1,NRING

         ISEGM = J
         ISSEG = 0
         NOPH(J)=0

         numanti = 0
         DO 1000 I=1,MAXPMA

            IF(NQGESKA(I,J).NE.0)THEN

               IDATA(1) = I
               RDATA(2) = TGESKA(I,J)
               IDATA(3) = NQGESKA(I,J)
               numanti=numanti+nqgeska(i,j)
               TGESKA(I,J) = 0.
               NQGESKA(I,J)= 0
               NDATA = 3
               ISSEG=ISSEG+1
               CALL KZREP2('MCRAWTQA',ISEGM,ISSEG,NDATA,IDATA)
               NOPH(J)=NOPH(J)+1
            ENDIF
 1000    CONTINUE
 900  CONTINUE

      ISEGM = 0
      ISSEG = 0
      NDATA = 1+2*NRING
      IBUF(1)=NRING	
      DO 1100 I=1,NRING
         IBUF(2*I)   = IPRT(I)
         IBUF(2*I+1) = NOPH(I)
         NOPH(I)=0
 1100 CONTINUE
      CALL KZREP2('MCRAWTQA',ISEGM,ISSEG,NDATA,IBUF)


!     QB based MCRAWTQA
!     N.B. there is no looping over rings here
      else if (SK_GEOMETRY .ge. SK_IV) then

       if( od_photon .ne. count_all_od_photon )then 
         isegm = 1
         isseg = 0
         noph(1)=0
         DO 1010 i=1,num_gate_all_anti
            if(nphoto_gate_anti(i).ne.0)then
               idata(1) = pmt_id_anti(i)
               rdata(2) = timing_gate_anti(i)
               idata(3) = nphoto_gate_anti(i)
               ndata = 3
               isseg=isseg+1
               call kzrep2('MCRAWTQA',isegm  ,isseg,ndata,idata)
               noph(1)=noph(1)+1
            endif
 1010    continue
       endif  ! if we apply gating within sgpmt

      isegm   = 0
      isseg   = 0
      ndata   = 3
      ibuf(1) = 1
      ibuf(2) = iprt(1)
      ibuf(3) = noph(1)

      ! has gating been applying in sgpmta?
      if( od_photon .ne. count_all_od_photon )then 
        call kzrep2('MCRAWTQA',0,0,3,ibuf)
      else  
      ! otherwise use all photons inside of dsanti_sk4
        ibuf(3) = nqgeska(1,1)
        call kzrep2('MCPHOTAN',0,0,3,ibuf)
      endif 

      noph(1) = 0
      nqgeska(1,1) = 0

      num_gate_all_anti = 0
      call vzero(pmt_id_anti,MXGATE)
      call vzero(nphoto_gate_anti,MXGATE)
      call vzero(timing_gate_anti,MXGATE)
      call vzero(width_photons_anti,MXGATE)
      call vzero(list_gate_anti,maxpma*50)
      call vzero(num_gate_anti,maxpma)

      endif
! ---------------------  End of MCRAWTQA creation




*     Trigger Simulation
      
      CALL DS
      endif

#ifdef FADC
      if (IFADCI) then
         if (flag_fadc_sube .eq. 0) then
*     Compress and Fill Bank
            write(6,'(A40)')'====== WRITING FADC DATA TO BANK ======='
            Call fadcwrite
            write(6,'(A40)')'========================================'
*     FADC initialization per event
            Call fadcwrinit
         endif
         flag_fadc_sube = 0
      endif
#endif

      if (OUTPUT_SCATTABLE.eq.1) then
         call fillscattable
      endif

      RETURN
      END
