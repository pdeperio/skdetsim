***********************************************************************
*     -----------------
      SUBROUTINE GUSTEP
*     -----------------
*
*     (Purpose)
*       This routine is called every one step of tracing is finished.
*     
*     (Input)
*       none
*
*     (Output)
*       none
*
*     (Creation Date and Author)
*       1992.04.24 ; First version by J.Kanzaki
*
*       1992.11.27 ; Modified for SK : Geantino for Cherencov light
*                           by Y.Suzuki/Y.Koshio
*       1998.03.16 ; Y.Hayato : add new secondary particle storage code
*
*       2002.08    ; C. Sterner : Modified for SK-2, add acrylic interfaces 
*       2006.04.15 ; SK-III version by G.Mitsuka
*
*       2006.05.10 ; Corrected a bug for OD segmentation. F.Dufour
*
*       2013.02.03 ; Clean up filling of secondaries bank and add more
*                    information. P. de Perio
*     
************************************************************************
      implicit none
#include "skdetsim_parm.h"
#include "geant321/gcflag.inc"
#include "geant321/gckine.inc"
#include "geant321/gcking.inc"
#include "geant321/gctrak.inc"
#include "geant321/gctmed.inc"
#include "geant321/gcvolu.inc"
#include "geant321/gcnum.inc"
#include "gumate.h"
#ifdef UPMU
#include "gumed_u.h"
#else
#include "gumed.h"
#endif
#include "skparm.h"
#include "skhead.h"
#include "skdetsim_test.h"
#include "guxyzpm.h"
#include "guskge.h"

**** Temporary
#include "gt2pd.h"
#include "skmc.h"
#include "apmultmp.h"

*** 1998/03/18 Add Y.H. new common for secondary particle
#include "vcwork.h"
#include "apscndry.h"
*** Muon kinematics common
#include "muvtx.h"
*** Need to read transmission probability from card file (F.Dufour)
#include "dsparma.h"
      

#include "geant321/gcphys.inc"
#include "geant321/gconsp.inc"
#include "geant321/gccuts.inc"

c      dimension rdata(3)
      integer iii, j, i, mecs, itmptra, isrchprt, onstack, samepart
      integer nphots_nucint, npions_nucint
      real ACRTRAN(3), distance, dummy, vert_temp(3), epsilon, absmom
      integer numedold, NOPMT, existflag, ingap, nearpmt
      character*4 cname

      parameter (epsilon=0.0001)

#ifdef LINAC
      integer linac_called /0/
      save linac_called
#endif

#ifdef UPMU
      real particlepos(3)
      real distancetodetector
#endif

***   Check for irregular removed pmt/acr/gap in SK-II
***   Find NOPMT
         !In SK-II some missing pmts could not be fit into the regular 
         !patterns defined by sggeom, they are removed individually here
      IF(INWVOL.EQ.1) THEN
         IF(NUMED.EQ.IMDPMC.OR.NUMED.EQ.IMDACR.OR.NUMED.EQ.IMDWTS) THEN
            CALL UHTOC(NAMES(5),4,CNAME,4)
C            CALL SGPMTNUM(CNAME,NUMBER(5),NUMBER(6),NUMBER(7)
C     &           ,NUMBER(8),NOPMT,existflag)
C     SGPMTNUM is splited into 3 routines (G.Mitsuka)
            IF(SK_GEOMETRY.eq.SK_I) then !Good Old SK-I
               CALL SGPMTNUM_SK1(CNAME,NUMBER(5),NUMBER(6),NUMBER(7)
     &              ,NUMBER(8),NOPMT,existflag)
            ELSE IF(SK_GEOMETRY.eq.SK_II) then !SK-II
               CALL SGPMTNUM_SK2(CNAME,NUMBER(5),NUMBER(6),NUMBER(7)
     &              ,NUMBER(8),NOPMT,existflag)
            ELSE IF(SK_GEOMETRY.eq.SK_III) then !SK-III
               CALL SGPMTNUM_SK3(CNAME,NUMBER(5),NUMBER(6),NUMBER(7)
     &              ,NUMBER(8),NOPMT,existflag)
            ELSE IF(SK_GEOMETRY.eq.SK_IV) then !SK-IV
               CALL SGPMTNUM_SK3(CNAME,NUMBER(5),NUMBER(6),NUMBER(7)
     &              ,NUMBER(8),NOPMT,existflag)
            ELSE IF(SK_GEOMETRY.eq.SK_V) then !SK-V
               CALL SGPMTNUM_SK3(CNAME,NUMBER(5),NUMBER(6),NUMBER(7)
     &              ,NUMBER(8),NOPMT,existflag)               
            ENDIF
            IF(existflag.eq.0) then !PMT and cover missing
               NUMED=IMDWTI !replace with ID water to remove
            ENDIF
            IF(NOPMT.gt. 7650 ) then
               if(ipart.eq.50) then
                  call grndm(dummy,1)
                  if(dummy.lt. 0.0001 ) then
                     istop=1
                     numedold=0
                     return
                  endif
               endif
            endif
            IF(OUTPUT_HBOOK.eq.OUT_TEST) then !store pmt existance
               if(existflag.eq.1) then
                  if(TEST_pmtmap(NOPMT).ne.1.) then
                     TEST_pmtmap(NOPMT)=1.
                     TEST_exist=1
                     TEST_gx=gxyzpm(1,nopmt)
                     TEST_gy=gxyzpm(2,nopmt)
                     TEST_gz=gxyzpm(3,nopmt)
                     TEST_NOPMT=NOPMT
                     call hfnt(23)
                  endif
               endif
            endif
         ENDIF
      ENDIF

***   "Cherenkov light"
      IF(IPART.EQ.50) THEN
**    Check if a Cherehkov photon goes into the same volume
         IF(INWVOL.EQ.0) THEN   !Same volume transport  

*     Check if the media is Water or Acrylic

c            IF(NUMED.EQ.IMDWTO.OR.NUMED.EQ.IMDWTI.OR.NUMED.EQ.IMDACR
c     &           .OR.NUMED.EQ.IMDWTS) then


* Use new defined volume IMDWOE (= outer endcap)
            IF(NUMED.EQ.IMDWTO.OR.NUMED.EQ.IMDWOE
     &           .OR.NUMED.EQ.IMDWTI.OR.NUMED.EQ.IMDACR
     &           .OR.NUMED.EQ.IMDWTS) then

               if(ISTOP .NE. 0) then
                  numedold=0
               ENDIF
               goto 999
#ifdef UPMU
*     kill cherenkov photon track in the dome
            ELSE IF ( NUMED.EQ.IMD_DOME_AIR) THEN
               istop = 1        
#endif
*     otherwise (error)
            ELSE
               WRITE(6,*) "ERROR in GUSTEP 1"
               istop = 1
            ENDIF


         ELSE IF(INWVOL.EQ.1) THEN !ENTER NEW VOLUME

            if(GEKIN.le.1.e-11) istop=1

            IF(numed.eq.imdwto.or.numed.eq.imdwoe) then !ENTER OD WATER
* All this part deals with OD segmentation. (F.Dufour, May 10th 2006)
               
               if ((numedold.eq.imdwto.and.numed.eq.imdwoe).or.
     &              (numed.eq.imdwto.and.numedold.eq.imdwoe)) then

                  call sgtyvk_odsegm
                  
                  istop=1
               else

                  goto 999
                  
               endif

            else IF(numed.eq.imdwti) then !ENTER ID WATER
               if(numedold.eq.IMDACR) then !ENTER ID WATER FROM ACRYLIC
                  call sgacrylic(ACRTRAN) !acrylic/water interface
                  istop=1

               else if(numedold.eq.IMDWTS) then
c                  !ENTER ID WATER FROM GAP WATER directly, rare error
                  goto 999
               else
                  goto 999   
               endif

            else IF(numed.eq.imdpmc) then !ENTER PMT
               if(numedold.eq.imdwti.OR.numedold.eq.IMDWTS) then
                                           !FROM ID WATER OR WATER GAP
                 call sgpmt(NOPMT)
                  istop=1
               else
c                  write(6,*) "C is generated in PMT: numedold=", numedold
                  istop=1
               endif
               
            else IF(numed.eq.IMDACR) then !ENTER ACRYLIC
c               Store acrylic medium position for later normal calculation
               DO III = 1,3
                  ACRTRAN(III) = GTRAN(III,NLEVEL)
               ENDDO

               if(numedold.eq.imdwti.OR.numedold.eq.IMDWTS) then 
c                  WRITE(*,*) 'ENTER ACRYLIC from water', numedold
                  call sgacrylic(ACRTRAN)  !acrylic/water interface
                  istop=1                
               else
                  goto 999
               endif 

               
            else IF(numed.eq.IMDWTS) then !ENTER WATER GAP BETWEEN ACR, PMT
               if(numedold.eq.IMDACR) then !FROM ACRYLIC
c                  WRITE(*,*) 'EXIT ACRYLIC INTO gap'
                  call sgacrylic(ACRTRAN)  !acrylic/water interface
                  istop=1

               else if(numedold.eq.IMDWTI) then
C               ENTER WATER GAP FROM ID WATER, rare error
                  goto 999
               else
                  goto 999
               endif 

            else IF(numed.eq.imdwtm) then !ENTER NON-SENSITIVE REGION

               if(numedold.eq.imdwto.or.numedold.eq.imdwoe) ! (Correction F.Dufour)
     &              then        !FROM OD (HIT TYVEK)
                  call sgtyvk
                  istop=1

               else if(numedold.eq.imdwti) then !FROM ID WATER (HIT BST)
                  call sgholebst(nearpmt,ingap,distance)
                  if(ingap.eq.0) then
                     call sgbst
                  else if(ingap.eq.1) then
                     if(distance.gt.sqrt((RACR-WIDACR/2.)**2
     &                    -(ZBLST+(RACR - ZACCY))**2)) then
                        call sgacredge(nearpmt,0)
                     else
                        call sggapedge(nearpmt)
                     endif
                     istop=1
                  endif

               else if(numedold.eq.IMDACR) then !FROM ACR (HIT acrylic edge)
c                  WRITE(*,*) 'EXIT ACRYLIC INTO edge', nstep
                  call sgacredge(NOPMT,1)  !edge of acrylic medium
                  istop=1
               else if(numedold.eq.IMDWTS) then 
                                    !FROM WATER GAP (between pmt, cover)
c                  WRITE(*,*) 'EXIT GAP INTO edge'
                  call sggapedge(nopmt)  !edge of water gap between acr and pmt
                  istop=1
               else if(numedold.eq.0) then
c                  write(6,*) "C is generated in nonsensitive region"
                  istop=1
               else
C                  write(6,*) "ERROR in GUSTEP 2"
                  istop = 1
               endif

            else IF(numed.eq.imdpma) then !ENTER OD PMTA
               if(numedold.eq.imdwto.or.numedold.eq.imdwoe) then !FROM OD WATER
                  call sgpmta
                  istop=1

               else if(numedold.eq.0) then
c                  write(6,*) "C is generated in PMA"
                  istop=1
               else

C                  write(6,*) "ERROR in GUSTEP 3"
                  istop = 1
               endif
            else IF(numed.eq.imdmst) then !EXIT TANK
               if(numedold.eq.imdwto.or.numedold.eq.imdwoe) then !FROM OD WATER
                  call sgtyvk
                  istop=1
               else if(numedold.eq.0) then
c     write(6,*) "C is generated out of TANK"
                  istop=1
               else
C     write(6,*) "ERROR in GUSTEP 4 Photon Exit Tank from",numedold
                  istop = 1
               endif
#ifdef UPMU
*     Exits into dome from the OD top cap
            else IF(numed.eq.imd_dome_air) then !EXIT TANK
               if(numedold.eq.imdwto.or.numedold.eq.imdwoe) then !FROM OD WATER
                  call sgtyvk
                  istop=1
               else if(numedold.eq.0) then
c     write(6,*) "C is generated out of TANK"
                  istop=1
               else
C     write(6,*) "ERROR in GUSTEP 4 Photon Exit Tank from",numedold
                  istop = 1
               endif
#endif               
               
#ifdef NICKEL
            else if(numed.ge.50) then ! hit nickel sources
               istop=1
#endif

#ifdef LINAC
            elseif(numed.ge.50)then ! hit linac geometry
               if(linac_called.eq.0) then
                  write(*,*) ' *************************************'
                  write(*,*) ' *** GUSTEP if for LINAC           ***'
                  write(*,*) ' *************************************'
                  call sleep(2)
                  linac_called=1
               endif
               if(numed.eq.50)then
                  call rfpip1
c                  call rfpip1_5
                  istop=1               
               else if(numed.eq.54)then
                  call rfpip1_6
                  istop=1               
               else if(numed.eq.55)then
                  call rfpip1_7
                  istop=1               
               else if(numed.eq.51)then
                 call rfpip2
                  istop=1               
               else if(numed.eq.52)then
                 call rffra
                  istop=1               
               else if(numed.eq.60)then
                  call rfwin
                  istop=1               
C--REP MN
C               else if(numed.eq.73)then
               else if(numed.eq.53)then
C--
                  call rfpip
                  istop=1   
CCS-03-05
               else if(numed.eq.73) then
                  istop=1
                  write(*,*) 'GUSTEP: warning rare LINAC lost photon in',numed
               else if(numed.eq.74) then
                  istop=1
                  write(*,*) 'GUSTEP: warning rare LINAC lost photon in',numed
               endif
#endif
            endif

         ELSE IF(INWVOL.EQ.2) THEN !EXIT VOLUME
            numedold=numed
            goto 999

         ELSE
            ISTOP=2

         ENDIF
         

         IF(istop.ne.0) THEN !Photon stopped
            numedold=0
         ENDIF

***   For other (non-Cherenkov photon) particles
      ELSE


**    Store secondery particle into stacks

C--
C--RUNNING PARTICLE / Charge EXCHANGE / Decay / etc... in the GEANT ---
C--

*** HAT 12/2011: check of NMEC from Patrick to bail if we are on 
*** the first event with nmec uninitialized
         if(nmec .le. 0)then
            if(ngkine.gt. 0)then
               write (*,*)"GUSTEP Warning: ngkine=", ngkine, 
     +                    ",but nmec=", nmec
            endif
            goto 104
         end if

C        Check that there's still space in secondary bank
         if(nscndprt.ge.secmaxrng) then
            if (flag_warning_scndprt .eq. 0) then
               flag_warning_scndprt = 1
               write(6,*) '***************************************************'
               write(6,*) 'warning; # of secondary particles >',secmaxrng
               write(6,*) 'Information of secondary particles is no longer stored'
               write(6,*) 'in CONVVECT bank hereafter in this event.'
               write(6,*) '***************************************************'
            endif
            goto 104
         endif

C        Select parent particle types to store (and decays)
C        (WARNING: Make sure to modify conditional below at "C 107" to keep
C         previous behaviour of CONVVECT bank, if modifying this conditional)

C         if( ( ((ipart.eq.1).or.(ipart.ge.5.and.ipart.le.12)
C     $        .or.(ipart.ge.14.and.ipart.le.24)
C     $        .or.(ipart.ge.26.and.ipart.le.32))

         if( ( ( ipart.eq.1 .or. (ipart.ge.5 .and. ipart.le.32) )
     $        .and.istop.ne.0.and.ngkine.ge.1 )
     $        .or. (lmec(nmec).eq.5 .and. ngkine.ge.1)
     $        ) then

            mecs=lmec(nmec)

C           Skip undesired secondary production mechanisms.
C           Don't skip KLKS regeneration mecs==42 (cregis), 
C           nor photonuclear absorption mecs==23,
C           nor Rayleigh effect mecs==25, nor below energy threshold mecs==30
C           (WARNING: Make sure to modify conditional below at "C 107" to keep
C            previous behaviour of CONVVECT bank, if modifying this conditional)
            if (mecs.gt.21 .and. mecs.ne.23 .and. mecs.ne.25 .and. mecs.ne.30 .and. mecs.ne.42) GOTO 104
            if (mecs.lt.5)  GOTO 104
            if (mecs.eq.10)  GOTO 104 ! Delta-ray 
C            if (mecs.ge.7 .and. mecs.le.10) goto 104
            if (mecs.ge.14.and.mecs.le.17) GOTO 104

C     Comment out to keep scattered pions
C            if (ngkine.eq.1) THEN
C               IF (int(gkin(5,1)+0.5).eq.ipart) GOTO 104
C            ENDIF

C           Count number of photons produced from proton interaction,
C           where the nucleon was originally above Cherenkov threshold,
C           for ignoring them in indexing and storage later
            nphots_nucint=0
            npions_nucint=0
            absmom = sqrt(pvert(1)**2+pvert(2)**2+pvert(3)**2)

C           For protons and neutrons
C            if ( (13.le.ipart .and. ipart.le.15) .or. ipart.eq.25 ) then

C           For hadronic interactions
            if (mecs.eq.12 .or. mecs.eq.13 .or. mecs.eq.19 .or. mecs.eq.20) then

               do j=1,ngkine
                     
C                 Count photons
C 106             (This 'if' conditional, including above, must be exactly the same as the one below at "C 106")
                  if (gkin(5,j).eq.1) nphots_nucint = nphots_nucint + 1

C                 Count pions
                  if (7.le.gkin(5,j).and.gkin(5,j).le.9) npions_nucint = npions_nucint + 1

               end do
               
C              Skip parent neutrons and below Cherenkov threshold parent protons 
C              if no pions nor neutrons were produced
               if ( (npions_nucint.eq.0) .and.
     $              ( (ipart.eq.13 .or. ipart.eq.25) .or.
     $              ( (ipart.eq.14 .or. ipart.eq.15) .and. absmom.lt.1.3) ) ) then 
                  goto 104
                  
               end if
               
            end if
            
            samepart = 0
            onstack = 0

C           Loop over secondaries generated in this step
            do 105 j=1,ngkine
                           
C              Skip photons from hadronic interactions
C 106          (This 'if' conditional must be exactly the same as the one above at "C 106")
               if ( gkin(5,j).eq.1 .and. 
     $              (mecs.eq.12 .or. mecs.eq.13 .or. mecs.eq.19 .or. mecs.eq.20) ) GOTO 105
C     $              ((13.le.ipart .and. ipart.le.15).or.ipart.eq.25) ) GOTO 105

               iprntid_temp  = 0
               ichildid_temp = 0

C     Only store this secondary if parent vertex exists on the primary or secondary stack
c     (e.g. to prevent storage of entire photon/electron shower)

C              Check initial track vertex (e.g. from NEUT)
               CALL GFVERT(1,vert_temp,dummy,dummy,dummy,dummy,dummy)
               do i=1,3

                  if (abs(vert_temp(i)-vert(i)).lt.epsilon) then
                     onstack=1
                  else
                     onstack=0
                  endif
                  
                  if (onstack.eq.0) exit
               enddo

C              Vertex matched initial vertex
               if (onstack.eq.1) then

C                 Look for the parent particle
                  samepart = 0
                  do isrchprt=1,NVC

C                    Check PID
                     if (ipdg(ipart).ne.IPVC(isrchprt)) cycle

C                    Check momentum
                     do i=1,3

                        if ( abs(pvert(i)*1000-pvc(i,isrchprt)) .lt. epsilon ) then
                           samepart = 1
                        else
                           samepart = 0
                        end if
                        if (samepart.eq.0) exit
                     enddo

C                    Not the same particle, check the next one
                     if (samepart.eq.0) then
                        cycle

C                    Found a matching particle!
                     else 
                        exit
                        
                     end if

                  enddo

C                 Assign negative index for particles found in MCVECT (NEUT) stack
                  if (samepart.eq.1) then
                     iprntid_temp = -isrchprt
                     nchildsvc(isrchprt) = ngkine

C                    Reverse child index since GEANT likes to track GCKINE stack backwards
                     ichildidxvc(isrchprt) = nscndprt+2-ngkine

                     if (nphots_nucint.ne.0) then
                        ichildidxvc(isrchprt) = ichildidxvc(isrchprt) + nphots_nucint
                        nchildsvc(isrchprt) = ngkine - nphots_nucint
                     end if
 
C                    Found a matching parent particle in the stack so continue
C                    filling this secondary
                     goto 103

C                 Could not find a matching particle in the MCVECT (NEUT) stack
C                 so search the secondaries stack
C                  else 
C                     goto 102
                     
                  end if
                  
               end if

C              Check secondary particles
               if(nscndprt.gt.secmaxrng) then
                  if (flag_warning_scndprt .eq. 0) then
                     flag_warning_scndprt = 1
                     write(6,*) '***************************************************'
                     write(6,*) 'warning; # of secondary particles >',secmaxrng
                     write(6,*) 'Information of secondary particles is no longer stored'
                     write(6,*) 'in CONVVECT bank hereafter in this event.'
                     write(6,*) '***************************************************'
                  endif
                  goto 104
               endif

               do iii=1, nscndprt

C                 Look for the parent particle
                  samepart = 0
                 
C                 Check if secondary vertex matches partent vertex
                  do i=1,3
                     if (abs(vtxscnd(i,iii)-vert(i)).lt.epsilon) then
                        onstack=1
                     else
                        onstack=0
                     endif
                     
                     if (onstack.eq.0) exit
                  enddo

                  if (onstack.eq.1) then
                     
C                    Check PID
                     if (ipdg(ipart).ne.iprtscnd(iii)) cycle
                                          
C                    Check momentum
                     do i=1,3
                        if ( abs(pvert(i)*1000-pscnd(i,iii)) .lt. epsilon ) then
                           samepart = 1
                        else
                           samepart = 0
                        end if
                        if (samepart.eq.0) exit
                     enddo

C                    Not the same particle, check the next one
                     if (samepart.eq.0) then
                        cycle

C                    Found a matching particle!
                     else 
                        exit

                     end if                     
                     
                  endif
                  
               enddo

C              Assign positive index for particles found in secondaries stack
               if (samepart.eq.1) then
                  iprntid_temp = iii
                  ichildid_temp = nscndprt+1

C                 Only fill child information for first child 
                  if (nchilds(iii).eq.0) then
                     ichildidx(iii) = ichildid_temp
                     nchilds(iii) = ngkine
                  end if

                  if (nphots_nucint.ne.0) then
                     nchilds(iii) = ngkine - nphots_nucint
                  end if


C                 Reset so that further secondaries don't have incorrect child index
                  ichildid_temp = 0
                     
C              Parent particle is not on any stack so skip filling this secondary
               else

C                 Fill all decays, hadronic interactions, bremsstrahlung, etc. for high Cherenkov
C                 yield particles even without matching vertex (since these tracks are killed at 
C                 each step, making tracking by parent vertex impossible)
                  if (mecs.ne.6 .and. mecs.ne.7 .and. mecs.ne.8 .and.
     $                mecs.ne.23 .and. mecs.ne.25 .and. mecs.ne.30) then
                     goto 103
                  end if

C                 Skip the rest (e.g. kills delta-ray photon interaction chains)
                  goto 104
                  
               end if

 103           continue

C 107          Keep old behaviour of VECT2 stack (especially for superscan)
               if ( (ipart.ne.1.and.ipart.ne.13.and.ipart.ne.14.and.ipart.ne.15.and.ipart.ne.25)
     &              .and. mecs.ne.23 .and. mecs.ne.25 .and. mecs.ne.30 .and.
     &              (mecs.lt.7 .or. mecs.gt.10) 
     &              ) then

                  if(npi0.ge.MAX2NDP) then
                     if (flag_warning_convvect .eq. 0) then
                        flag_warning_convvect = 1
                        write(6,*) '***************************************************'
                        write(6,*) 'warning; # of secondary particles >',MAX2NDP
                        write(6,*) 'Information of secondary particles is no longer stored'
                        write(6,*) 'in CONVVECT bank hereafter in this event.'
                        write(6,*) '***************************************************'
                     endif
                     goto 104
                  endif

                  npi0=npi0+1
                  
                  ippi0(npi0)       =ipdg(ipart)
                  imeca(npi0)       =mecs
                  itrcon(npi0)      =itra
                  iprtpi0(npi0)     =ipdg(int(gkin(5,j)+0.5))
                  tpi0(npi0)        =tofg*1.e+9

                  do i=1,3
                     vpi0(i,npi0)   =gpos(i,j)
                     ppi0(i,npi0)   =gkin(i,j)*1000.
                  enddo

               end if
C              End of filling VECT2 (CONVVECT) stack
               
C               Check that current track ID is the same as original parent ID
C               (Should always be true)
C               if (ipart.ne.ipaold) then
C                  write(*,*) ipart, " != ", ipaold
C                  stop
C               endif

C              Start filling full information secondaries stack

C              Store this secondary track and vertex permanently
C              for matching later
C               IFLGK(J) = 1

C              Check that there's still space in secondary bank
               if(nscndprt.ge.secmaxrng) then
                  if (flag_warning_scndprt .eq. 0) then
                     flag_warning_scndprt = 1
                     write(6,*) '***************************************************'
                     write(6,*) 'warning; # of secondary particles >',secmaxrng
                     write(6,*) 'Information of secondary particles is no longer stored'
                     write(6,*) 'in CONVVECT bank hereafter in this event.'
                     write(6,*) '***************************************************'
                  endif
                  goto 104
               endif

               nscndprt=nscndprt+1

               itrkscnd(nscndprt)=itra
               istakscnd(nscndprt)=istak
               iprntprt(nscndprt)=ipdg(ipart)
               lmecscnd(nscndprt)=mecs
               iprtscnd(nscndprt)=ipdg(int(gkin(5,j)+0.5))
               tscnd(nscndprt)   =tofg*1.e+9
               iflgscnd(nscndprt)=iflgscnd_temp
               iprntidx(nscndprt)=iprntid_temp
               nchilds(nscndprt) =0
               ichildidx(nscndprt)=ichildid_temp

               do i=1,3
                  vtxscnd(i,nscndprt)=gpos(i,j)
                  pscnd(i,nscndprt)  =gkin(i,j)*1000.
                  pprnt(i,nscndprt)  =VECT(7)*VECT(I+3)*1000.  
                  pprntinit(i,nscndprt)=pvert(i)*1000
                  vtxprnt(i,nscndprt)=vert(i)
               enddo               

               itmptra=0
               do isrchprt=1,NVC
                  if (ICRNVC(isrchprt).eq.1) then
                     itmptra=itmptra+1
                  endif
                  if (itmptra.eq.itra) goto 505
               enddo
 505           iprnttrk(nscndprt)=isrchprt
               iorgprt(nscndprt) =ipvc(isrchprt)

 105        continue

         endif

 104     continue

         iflgscnd_temp=0        ! reset
            
C Alec, Sept. 2002
C Save muon kinematics if a muon enters or leaves ID or OD
C 
C primary check - was a volume entered by a primary muon?
         if (inwvol.eq.1.and.(ipart.eq.5.or.ipart.eq.6)
     &        .and.(itra.eq.first_id_mu)) 
     &        call save_muon_path
         
C If we have a zero momentum muon, save the stopping info.
C Note that this is an independant way of figuring stopping that the more
C complicated way in skmkapmcstop.F
         if (vect(7).eq.0.0.and.(ipart.eq.5.or.ipart.eq.6)) then
C is this the last primary muon?
            if (first_id_mu.eq.itra) then
c     print*,'*** stopping point ',
c     &              ipart,itra,(vect(iii),iii=1,7),idmu
               do iii=1,3
                  mupos(3,iii)=vect(iii)
                  mudir(3,iii)=vect(iii+3)
               enddo
               mumom(3)=vect(7)
            endif ! first muon in ID
         endif                  ! check for stops

#ifdef UPMU
C Add rock component for upmu 2007 4/20 T.Tanaka
 
C Check for electrons or gammas in the rock, and kill them
C Or should we put this back up at the child particle level?  That would
C be more efficient, but I don't trust it.
         if ((numed.eq.imdmst).and.
     &        (ipart.eq.2.or.ipart.eq.3.or.ipart.eq.1)) then
c     because this is for OD simulations, we allow some electromagnetic particles to propagate near the outside of the detector
c     we keep gammas and electrons which are within 150 cm of the detector. This is about 10 radiation lengths in rock.
            do iii=1,3
               particlepos(iii) = vect(iii)
            enddo
            call disttodet( particlepos, distancetodetector )
            if (distancetodetector.gt.150) then
               istop = 1
            endif
         endif
c         if ((numed.eq.imdmst.or.numed.eq.IMD_DOME_AIR).and.
c     &        (ipart.eq.2.or.ipart.eq.3.or.ipart.eq.1)) then
c            istop = 1
c         endif
#endif

         IF(NGKINE.GT.0) then
            CALL GSKING(0)
            do I=1,NGKINE
               IFLGK(I) = 0
            end do
         end if
         NGKINE=0
         
         IF(NGPHOT.GT.0) CALL GSKPHO(0)
         NGPHOT=0

C Alec, 9/02 -
C The commented section below is ancient cruft, I have no idea where it
C came from, but it's not part of the muon vertex changes I made above
**    When the particle goes into inner tank from anti tank,
**    store the kinematics
c         if(inwvol.eq.1.and.numed.eq.imdwti.and.numedoldp.eq.imdwtm)then
c            if(ipart.eq.6)then
c               rdata(1)=vect(1)
c               rdata(2)=vect(2)
c               rdata(3)=vect(3)
c               call kzput1('MCVERTEX',1,0,3,rdata)
c               rdata(1)=vect(4)*vect(7)*1000.
c               rdata(2)=vect(5)*vect(7)*1000.
c               rdata(3)=vect(6)*vect(7)*1000.
c               call kzput1('MCVECT',1,2,3,rdata)
c               print*,ipart,(vect(iii),iii=1,6)
c               numedoldp=0
c            endif
c         endif
c         if(inwvol.eq.2)numedoldp=numed
c         if(istop.ne.0)numedoldp=0
      END IF
      
 999  RETURN
      END
