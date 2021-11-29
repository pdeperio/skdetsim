************************************************************************
*     -------------------------------------
      SUBROUTINE DSANTI_SK4
*     -------------------------------------
*
*     (Purpose)
*        Calculate TQ of the ANTI PMT
*
*     (Input)
*
*     (Output)
*        'TQA' bank
*
*     (Creation Date and Author)
*        2008.12.25 : update for sk-iv 
*                   : R. Wendell (rvw)
*
*        2019.07.12 : update for sk-V
*                   : Baran Bodur
*
*
************************************************************************
      IMPLICIT NONE

#include "skparm.h"
#include "dsparma.h"
#include "skmc.h"
#include "dsoffset.h"
#include "skdetsim_parm.h"
#include "guflag.h"
#include "skhead.h"
#include "geopma.h"
#include "geotnk.h"

#include "dsdata.h"
#include "dsparm.h"
! so that we may interface with the software trigger
C#include "sktq.h"
#include "sktq_raw.h"


      INTEGER IBUF(MAXRNG*2+1)
      REAL RBUF(MAXRNG*2+1)
      EQUIVALENCE (IBUF(1),RBUF(1))
      REAL TLIST(MAXPMA), QLIST(MAXPMA), q_tmp
      INTEGER NUM(MAXPMA)
      INTEGER IOLD1, IOLD2
      INTEGER I, IERR, NW, NDATA,j  
      INTEGER NRING, IRING, ICABLE, NHIT, IHTMX, IHITS
      REAL RDUMMY, AOFFST
!     3 = top, 2 = bottom, 1 = wall
      REAL RDARK(3), DARKT, PE1
      REAL SATPARAM(MAXPMA)
      REAL EFF_FAC(MAXPMA)
      REAL PED(MAXPMA), SIGPED
      REAL SMEAR, SIGSMEAR
      LOGICAL SET_GAIN, EXP_RSP, PLTEXP_RSP, GAUS_RSP
      LOGICAL SET_SAT, SET_EFF
      DATA SIGPED, SIGSMEAR/0.12,0.2/
      DATA SET_GAIN/.TRUE./
      DATA SET_SAT/.TRUE./
      DATA SET_EFF/.TRUE./
      DATA EXP_RSP, PLTEXP_RSP, GAUS_RSP/.FALSE.,.TRUE., .FALSE./

      external GRNDM

! following dsres_sk4.F, MXPHDS is found in dsdata.h
      integer nphoton
      integer nphoton_tmp
      integer photon_pmt_list(MXPHDS)
      real    photon_t_list(MXPHDS)
      real    q

      real    gate_tlist(MXGATE), dummy_list(MXGATE)
      integer gate_clist(MXGATE),gate_qlist(MXGATE)
      integer ngates
      integer gate
      integer n_ingate_by_track
      integer ndark, dk_i 
      integer locale
      integer skGeom

      integer iq_small, iq_medium, iq_large, flag_gain
      real pc2pe_mc, cnt2pc_s, cnt2pc_m, cnt2pc_l, pc2pe_new, pc2pe_old
      parameter(pc2pe_new=2.465)
      parameter(pc2pe_old=4.2)
      parameter(cnt2pc_s=0.100, cnt2pc_m=0.732, cnt2pc_l=5.141)
      real tmin, tmax
      integer darkRun
    
      ! prevent OD simulation when requested -rvw
      if( ifanti.eq.0 ) return
      
      !get sk geometry
      skGeom = SK_GEOMETRY ! use this when skdetsim is updated
      !skGeom = 5 ! this was to test skdetsim sk-5 update prior to 

      ! darkRun should be set by card file eventually
      darkRun = ODDARKRUN

      ! open the root single pe histograms 
      call open_root(skGeom,darkRun) !  
      print *, "od info loaded: ", skGeom , "\n"
      
      ! if sk5 then read pedestal and imb smearing from card
      if(skGeom.ge.5) then
        SIGPED = ODPEDSMEAR
        SIGSMEAR = ODIMBSMEAR
      endif
      
      !print *, "ped semaring: ", sigped
      !print *, "imb semaring: ", sigsmear
      
      ! random pedestals and store seeds 
      IF (SET_GAIN) THEN
         CALL GRNDMQ(IOLD1,IOLD2,0,'G')     ! store old seed
         CALL GRNDMQ(20030412,92612,1,'S')  ! initialize new seed
         CALL RNORMX(PED,MAXPMA,GRNDM)      ! generate random pedestals
         DO I = 1, MAXPMA
            PED(I) = PED(I)*SIGPED
         END DO
         CALL GRNDMQ(IOLD1,IOLD2,1,'S')     ! restore old seed
         SET_GAIN = .FALSE.
      END IF


      IF (SET_SAT) THEN
         CALL GRNDMQ(IOLD1,IOLD2,0,'G')     ! store old seed
         CALL GRNDMQ(20030442,921612,1,'S')  ! initialize new seed
         DO I = 1, MAXPMA
         ! tube-by-tube sat parameters
         ! computed by looking at maximum charge 
         ! distributions in the data
         
         ! New in 2019,sk5: computed from laser data for 538 pmts,
         ! other were selected randomly (but based on type) from 
         !the distribution of those 538 pmts saturation...
            call od_pmt_sat( satparam(i), i )
         ENDDO

         CALL GRNDMQ(IOLD1,IOLD2,1,'S')     ! restore old seed
         SET_SAT = .FALSE.
      END IF

      IF (SET_EFF) THEN
         CALL GRNDMQ(IOLD1,IOLD2,0,'G')     ! store old seed
         CALL GRNDMQ(20050442,921612,1,'S')  ! initialize new seed
         DO I = 1, MAXPMA
 
! Apply extra efficiency factor for old tubes.  Note this is only
!   correct if old tubes have smaller efficiency, apmtefr<1
!  (if not, should modify sgpmta.F)
            IF (OD_TUBE_TYPE(I).EQ.2) THEN
               EFF_FAC(I) = APMTEFR
            ELSE
               EFF_FAC(I) = 1.
            ENDIF

! Now apply efficiency factor for top/bottom/wall

            IF (XYZPMA(3,I).GE.ZPINTK) THEN
               EFF_FAC(I) = EFF_FAC(I)*APMTEFT
            ELSE IF (XYZPMA(3,I).LE.-ZPINTK) THEN
               EFF_FAC(I) = EFF_FAC(I)*APMTEFB
            ELSE
               EFF_FAC(I) = EFF_FAC(I)*APMTEFW
            ENDIF

         ENDDO
         CALL GRNDMQ(IOLD1,IOLD2,1,'S')     ! restore old seed
         SET_EFF = .FALSE.

      ENDIF

!
!     Initialize arrays
!
      CALL VZERO(TLIST,MAXPMA)  ! list of time for each photon hit
      CALL VZERO(QLIST,MAXPMA)  ! list of charge of each photon hit

      ! zero out arrays that will be used for the software trigger
      nhitaz_raw = 0
      call vzero( itaskz_raw , 30*maxpma ) 
      call vzero( iqaskz_raw , 30*maxpma ) 
      call vzero( taskz_raw  , 30*maxpma ) 
      call vzero( qaskz_raw  , 30*maxpma ) 
      call vzero( icabaz_raw , 30*maxpma ) 


!     Read MCPHOTAN bank
      if( od_photon .eq. count_all_od_photon )then 
        call kzget2('MCPHOTAN', 0 ,0, NW, IBUF )
        nphoton = ibuf(3)
        do i= 1, nphoton
           nw = 2
           call kzget2('MCPHOTAN', 1 ,i, NW, IBUF )
           photon_pmt_list(i) = IBUF(1)
           photon_t_list(i)   = RBUF(2)*1.e9 + AOFFST
        enddo
      endif ! end of counting hits with all photons

!
!     Loop over M.C. particles
!

!     Read MCRAWTQA bank
      if( od_photon .ne. count_all_od_photon )then 
         call kzget2('MCRAWTQA',0,0,NW,IBUF)
         NRING=IBUF(1)
         DO I=1,NRING
            NUM(I) = IBUF(2*I+1)
         END DO
         CALL KZGET2('MCRAWTQA',1,1,NW,IBUF)

!     in sk4 we want to store all of the times then count the in gate
!     hits using dscnt.F
         nphoton = 0
         DO IRING=1,NRING
            DO I=1,NUM(IRING)
 
               CALL KZGET2('MCRAWTQA',IRING,I,NW,IBUF)
 
               ICABLE = IBUF(1)
 
               AOFFST=COMOFST
               ! convert to shifted [ns]
               RBUF(2) = RBUF(2)*1.e9 + AOFFST
 
               n_ingate_by_track   = ibuf(3)
               ! this is done to provide all of the hits to dscnt
               do j = 1 , n_ingate_by_track
                  nphoton = nphoton + 1
                  photon_pmt_list(nphoton) = IBUF(1)
                  photon_t_list(nphoton)   = RBUF(2)
               enddo
            enddo ! loop on hits for this particle
         enddo ! loop on particle
      endif ! end of read from mcrawtqa


      NHIT = 0

      ! use the same time window as ID (see dsnois_sk4.F, dsdata.h)
      tmin = -trg_offset - twid_block*2
#ifdef NEUTRON
      tmax = twid_block
#else
      tmax = (twid_block - trg_offset) + twid_block*3
#endif

      ! time in [ns]
      FULWIDA = tmax - tmin
      RDARK(1)   = darkdsa(1)*FULWIDA*1.e-9
      RDARK(2)   = darkdsa(2)*FULWIDA*1.e-9
      RDARK(3)   = darkdsa(3)*FULWIDA*1.e-9

!     Loop over all tubes/gates, add up hits for each tube
      DO I=1,MAXPMA

         ! select dark rate by tank geometry
         if (xyzpma(3,i).ge. zpintk) then
            locale = 3 ! top
         else if ( xyzpma(3,i).le. -zpintk) then
            locale = 2 ! bottom
         else
            locale = 1 ! barrel
         endif 

! rvw 20091020 odtuning, try individual dark rates
        ! returns value in Hz
        
        ! call od_pmt_dark( rdark(locale), I, sk_geometry )          
        
        call od_pmt_dark( rdark(locale), I, skGeom, darkRun ) !can select darkrates from run         
           
        ! 20091026 
        ! since the measured dark rates are taken after tube efficiency 
        ! effects and input to the simulation we 
        ! correct for the efficiency factor 
        ! to make the dark rates correct in the output. 
         rdark(locale) = rdark(locale)*FULWIDA*1.e-9/eff_fac(I)

         call gpoiss(rdark(locale),ndark,1)

         do dk_i = 1, ndark 
!           Make dark hits: sk4 allows for multiple hits so 
!           disperse the darknoise throughout the time window
!           allow noise in blocks [N-2, N+3] 

            call grndm(rdummy,1)
            darkt=rdummy*fulwida+tmin-DSANTT0+AOFFST

            nphoton = nphoton + 1
            photon_pmt_list(nphoton) = I 
            ! darkt is in [ns]
            photon_t_list(nphoton)   = darkt 
         enddo
      enddo   ! end of loop over PMTs


!     routine returns lists indexed by pmt gate number
!     thus allowing for multiple hits
      call dscnt(nphoton , photon_t_list, photon_pmt_list, 
     &           ngates  , gate_clist   , gate_tlist     , gate_qlist, dummy_list)


      ! gate_tlist is in [ns]
      do 400 gate = 1, ngates

         i = gate_clist( gate )   ! store the cable number

         ihtmx = gate_qlist(gate)
         QLIST(I)=0.

!        Loop over hits for this tube/gate
         DO 401 IHITS=1,IHTMX
!           Apply extra efficiency factor to this hit
            CALL GRNDM(RDUMMY,1)

            if (rdummy.gt.EFF_FAC(I)) goto 401

            !!! NEXT PART is al TEMPPP !!! CHANGE BACK

            if(skGeom.lt.5) then
              if (od_tube_type(i).eq.2) then
                call od_old1pe( pe1 )
              else
                call od_new1pe( pe1 )
              endif
            elseif(skGeom.ge.5) then
              if(ODINDSPE) then
                !print *, "Using Individual PMT SPE"
                call od_pmt1pe(pe1, i) !spe by individual pmt
              else
                !print *, "Using SPE dist by PMT type"
                call od_type1pe(pe1, i) !spe by type
             endif
           endif
  

            qlist(i)=qlist(i)+pe1
 401     CONTINUE ! end of loop on gates


!        Apply discriminator threshold
         IF (OD_TUBE_TYPE(I).EQ.1) THEN
            IF (QLIST(I) .LE. ADSCTHR*ADSTHRF) THEN
               QLIST(I) = 0.
               TLIST(I) = 0.
               GO TO 400               ! skip to the next gate
            ENDIF
         ELSE
            IF (QLIST(I) .LE. ADSCTHR) THEN
               QLIST(I) = 0.
               TLIST(I) = 0.
               GO TO 400               ! skip to the next gate
            ENDIF
         ENDIF



!        Final noise smearing after threshold for old tubes

         IF (OD_TUBE_TYPE(I).EQ.2) THEN
            CALL RNORMX(RDUMMY,1,GRNDM)
            SMEAR = 1. + RDUMMY*SIGSMEAR
            QLIST(I) = QLIST(I)*SMEAR
         endif

!
!        Pedestal smearing
!
         QLIST(I) = QLIST(I) + PED(I)
         

!
! tube-by-tube saturation curves
! 
         q = qlist(i)
         q = q/(1+SATPARAM(I)*q) 
         !if(q.ge.2) then
         !  print *, "PMT: ", I, " Sat param: ", SATPARAM(I)
         !  print *, "qReal: ", qlist(i), " qSaturated: ", q
         !endif
         qlist(i) = q


!        Suppress negative values (need to understand real data behavior- ETK)
         IF (QLIST(I) .LT. 0.0) goto 400 
         nhit = nhit + 1

         if(skGeom.lt.5) then
           if (od_tube_type(i).eq.2) then
             pc2pe_mc = pc2pe_old
           else
             pc2pe_mc = pc2pe_new
           endif
         elseif(skGeom.ge.5) then ! use actual pc2pe value
           call od_pc2pe(pc2pe_mc,i)
         endif

         !print *, "PMT: ", i, " od_pc2pe: ", pc2pe_mc
         
         iq_small  = int( qlist(i)/pc2pe_mc/cnt2pc_s ) + 961
         iq_medium = int( qlist(i)/pc2pe_mc/cnt2pc_m ) + 961
         iq_large  = int( qlist(i)/pc2pe_mc/cnt2pc_l ) + 961


         ! originally was a single channel conversion
         ! qcnvsk = 13.34 Tcounts / pe
         ! qlist(i)=qlist(i) * qcnvsk
                                   q_tmp    = qlist(i)/pc2pe_mc/cnt2pc_s 
         if( iq_small  .ge. 1350 ) q_tmp    = qlist(i)/pc2pe_mc/cnt2pc_m 
         if( iq_medium .ge. 1350 ) q_tmp    = qlist(i)/pc2pe_mc/cnt2pc_l

!       set the ADC gain level
                                   flag_gain = 0
        if ( iq_small  .ge. 1350 ) flag_gain = 1
        if ( iq_medium .ge. 1350 ) flag_gain = 2

!       Saturation of the QB ADC
        if ( q_tmp .ge. 2**11-1) q_tmp = 2**11-1

!        i  is the cable number
!        tdc_unit_sk4 = 0.52 [ cnt / ns ], (gate_tlist is in ns)
!        units are taskz: [ns]  , qaskz: [ pe ]
         taskz_raw( nhit ) = gate_tlist(gate)
         qaskz_raw( nhit ) = qlist(i)

         ! qlist is now in counts 
         qlist(i) = q_tmp +  ( flag_gain * 2**14 );


!        i is the cable for this gate, it is update for every hit
!        fill variables that are in counts
!        this information will be used to fill the TQAZ banks, and run the 
!        software trigger
         itaskz_raw( nhit ) = int( (taskz_raw(nhit)+trg_offset) / tdc_unit_sk4 ) +(ntrg+2)*2**15
         iqaskz_raw( nhit ) = qlist( i )

C-- Add 20000 here! ( 2009/05/06 Y.Hayato )
         icabaz_raw( nhit ) = i + 20000
                                 ! software trigger needs normal "c" indexing
                  
 400  CONTINUE    !! end of loop on gates
       nhitaz_raw = nhit
      ! commons to access the software trigger should be filled by now
      ! close od single pe histograms
      
      call close_root(skGeom,darkRun)

! ---
      RETURN
      END
