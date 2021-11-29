***********************************************************************
*     -----------------
      SUBROUTINE SGNUCEFF
*     -----------------
*
*     (Purpose)
*       GET THE KIND OF INTERACTION(INTKND) FOR THE GIVEN MOMENTUM(P)
*       USING NEUT LIBRARIES
*
*     (Input)
*       none
*
*     (Output)
*       none
*
*     (Creation Date and Author)
*       2011.06.10 ; First version   by  P.de Perio
*
************************************************************************
      implicit none

#include "geant321/gccuts.inc"
#include "geant321/gckine.inc"
#include "geant321/gcking.inc"
#include "geant321/gctrak.inc"
#include "geant321/gctmed.inc"

#include "skcnst.h"
#include "gt2pd.h"
#include "apscndry.h"

#ifdef SITE_HAS_NEUT
#include "vcwork.h"
#include "posinnuc.h"
#include "efpion.h"
#include "fsihist.h"
#include "necard.h"

#include "skdetsim_parm.h"
C Copied from vcwork.h to temporarily store initial particle stack
C before it is overwritten by NEUT cascade in this function
      
      INTEGER NVC_temp, IPVC_temp, IORGVC_temp, IFLGVC_temp
      INTEGER ICRNVC_temp, IVTIVC_temp, IVTFVC_temp
      REAL    POSVC_temp, AMASVC_temp, PVC_temp
      REAL TIMVC_temp, POSIVC_temp, POSFVC_temp
      COMMON /VCWORK_temp/ NVC_temp,POSVC_temp(3),
     $    IPVC_temp(MAXVC),AMASVC_temp(MAXVC),PVC_temp(3,MAXVC),
     $    IORGVC_temp(MAXVC),IFLGVC_temp(MAXVC),ICRNVC_temp(MAXVC),
     $    TIMVC_temp(MAXVC),POSIVC_temp(3,MAXVC),IVTIVC_temp(MAXVC),
     $    POSFVC_temp(3,MAXVC),IVTFVC_temp(MAXVC)

#endif

C#define DEBUG

      REAL P, T, COST, CONE, R
      REAL PVECT(3)
      INTEGER I, J, IVC, NLOOPS
      INTEGER GEANTID,PDGID, ireaction

      REAL tot,absb,inel,scx,dcx,elas,rxn
      REAL PABS,PINEL,PSCX,PDCX,PELAS,PRXN,PROBSUM

      REAL ROH,ROO
      DATA ROH,ROO/6.689E22,3.344E22/
      REAL totcrswater, crshydrogen

      REAL SIGEL,SIGINEL,SIGCX,SIGTOT
      INTEGER PICHARGE,NUCCHARGE,IERROR

      REAL RMAX, RAD, posv(3)

      REAL PIPINT
      EXTERNAL PIPINT

      REAL PIOINTNEUT
      EXTERNAL PIOINTNEUT

      REAL prob_radcapture
      EXTERNAL prob_radcapture

      REAL DIR(3)

C new variable for pi- capture
      REAL pcm2(2,4), pcm3(3,4)

C      real*4    EVPIPROB
C      external  EVPIPROB
c---------for calculating impact position on nucleus      
      real PMOM_MAG,SINPHI,COSPHI,COSTHETA,SINTHETA
      REAL PVECTOUT(3)

C For FSI history common block
      INTEGER NVERT_TMP, NVCVERT_TMP, B_OXYGEN, B_STORE_FREE

C For pi- capture
      REAL DIRDUMMY(3), g_cap_ene, PMOM1(3), XMASS1, UV(3)
      REAL NMASS, PIMASS, GAM, W, AMAS(10), PMOMD(3,10), beta(4)
      integer ID, NDECAY, IERR

#ifdef WithG4
C----- Geant4 pi- capture: Added by RA Sep. 2018
      REAL*8 secVec(9,250)
      INTEGER Nscnd,iscnd
#endif


#ifdef SITE_HAS_NEUT
      IF(GEKIN.LT.CUTHAD)THEN
         ISTOP = 1
         NMEC = NMEC +1
         LMEC(NMEC) = 5
         
         IF(IPART.EQ.8) then
            CALL GDECAY

         ELSE

#ifdef WithG4
C         pi- capture by Geant4
C         RA 2018
          if (USEG4.eq.1) then
            CALL g4_pimicap(VECT(4)*VECT(7), VECT(5)*VECT(7),
     +                      VECT(6)*VECT(7), Nscnd, secVec)
            NGKINE=0
            DO iscnd=1, Nscnd
                NGKINE=NGKINE+1
                GKIN(1,NGKINE) =secVec(6,iscnd)
                GKIN(2,NGKINE) =secVec(7,iscnd)
                GKIN(3,NGKINE) =secVec(8,iscnd)
                GKIN(4,NGKINE) =secVec(9,iscnd)
                GKIN(5,NGKINE) =secVec(1,iscnd)
                GPOS(1,NGKINE) =VECT(1)
                GPOS(2,NGKINE) =VECT(2)
                GPOS(3,NGKINE) =VECT(3)
                TOFD(NGKINE) = 0
            END DO
          else 
#endif
C        pi- capture written by Y.SONODA, 2017.12.01
C         NOTICE : Not calculating any nuclear effects. ( fermi motion. etc...)
            ID = 2112
            call mcmass(ID, NMASS)
            ID = -211
            call mcmass(ID, PIMASS)

            call grndm(r,1)
C           2% -> radiative capture
            if ( r .lt. 0.02 ) then !1
              call grndm(r,1)
              g_cap_ene = prob_radcapture(r)

C             gamma energy greater than 100MeV -> 1 gammma + 1 neutron
              if ( g_cap_ene .gt. 100. ) then
                NGKINE=2
                GKIN(5,1) = 1
                GKIN(5,2) = 13
                GKIN(4,1) = g_cap_ene * 0.001

                call rndir(dirdummy)
                DO I=1,3
                  GKIN(I,1) = GKIN(4,1) * dirdummy(I)
                  GKIN(I,2) = -GKIN(I,1) 
                  GPOS(I,1)=VECT(I)
                  GPOS(I,2)=VECT(I)
                end do
                GKIN(4,2) = sqrt(GKIN(1,1)**2+GKIN(2,1)**2+GKIN(3,1)**2+(NMASS*0.001)**2)

                IFLGK(1) = 0
                IFLGK(2) = 0
                TOFD(1) = 0.
                TOFD(2) = 0.

                NMEC=NMEC+1
                LMEC(NMEC) = 12

C             gamma energy greater than 100MeV -> 1 gammma + 2 neutron
              else
                NGKINE=3
                GKIN(5,1) = 1
                GKIN(5,2) = 13
                GKIN(5,3) = 13
                GKIN(4,1) = g_cap_ene * 0.001
                call rndir(dirdummy)
                DO I=1,3
                  GKIN(I,1) = GKIN(4,1) * dirdummy(I)
                end do

                W = PIMASS + (2*NMASS) - g_cap_ene
                XMASS1 = sqrt(W**2 - g_cap_ene**2)
                do i=1,3
                  PMOM1(I) = -GKIN(I,1) * 1000
                enddo
                BETA(1) = - PMOM1(1) / W
                BETA(2) = - PMOM1(2) / W
                BETA(3) = - PMOM1(3) / W
                BETA(4) = W / XMASS1
                UV(1) = BETA(1) / sqrt(BETA(1)**2+BETA(2)**2+BETA(3)**2)
                UV(2) = BETA(2) / sqrt(BETA(1)**2+BETA(2)**2+BETA(3)**2)
                UV(3) = BETA(3) / sqrt(BETA(1)**2+BETA(2)**2+BETA(3)**2)
                GAM = BETA(4)
                call mcvecbst(PMOM1, XMASS1, UV, GAM)

                W = sqrt(PMOM1(1)**2+PMOM1(2)**2+PMOM1(3)**2+XMASS1**2)
                AMAS(1) = NMASS
                AMAS(2) = NMASS
                NDECAY = 2
                call VCPHSVCT( W, AMAS, NDECAY, PMOMD, IERR)

                do i=1,3
                  UV(I) = -UV(I)
                  PMOM1(I) = PMOMD(I,1)
                enddo
                call mcvecbst(PMOM1, NMASS, UV, GAM)
                GKIN(1,2) = PMOM1(1) * 0.001
                GKIN(2,2) = PMOM1(2) * 0.001
                GKIN(3,2) = PMOM1(3) * 0.001
                GKIN(4,2) = sqrt(GKIN(1,2)**2+GKIN(2,2)**2+GKIN(3,2)**2+(NMASS*0.001)**2)

                do i=1,3
                  PMOM1(I) = PMOMD(I,2)
                enddo
                call mcvecbst(PMOM1, NMASS, UV, GAM)
                GKIN(1,3) = PMOM1(1) * 0.001
                GKIN(2,3) = PMOM1(2) * 0.001
                GKIN(3,3) = PMOM1(3) * 0.001
                GKIN(4,3) = sqrt(GKIN(1,3)**2+GKIN(2,3)**2+GKIN(3,3)**2+(NMASS*0.001)**2)

                DO I=1,3
                  GPOS(I,1)=VECT(I)
                  GPOS(I,2)=VECT(I)
                  GPOS(I,3)=VECT(I)
                  IFLGK(I) = 0
                  TOFD(I) = 0.
                end do

                NMEC=NMEC+1
                LMEC(NMEC) = 12
              end if 

C           98% -> two neutron emission
C           pi- p n -> n n : 73.36%, pi- p p -> p n : 17.64%
C           this ratio is tentative. it should be modified by reference.
            else
              call grndm(r,1)
              if ( r .lt. 0.17647 ) then
                NGKINE=1
                GKIN(5,1) = 13
                GKIN(4,1) = (NMASS + PIMASS) * 0.001
                call rndir(dirdummy)
                DO I=1,3
                  GKIN(I,1) = sqrt(GKIN(4,1)**2 - (NMASS*0.001)**2)*DIRDUMMY(I)
                  GPOS(I,1) = VECT(I)
                enddo
                IFLGK(1) = 0
                TOFD(1) = 0.
                NMEC=NMEC+1
                LMEC(NMEC) = 12

              else
                NGKINE=2
                GKIN(5,1) = 13
                GKIN(5,2) = 13
                GKIN(4,1) = (NMASS + PIMASS / 2.0) * 0.001
                GKIN(4,2) = (NMASS + PIMASS / 2.0) * 0.001

                call rndir(dirdummy)
                DO I=1,3
                  GKIN(I,1) =  sqrt(GKIN(4,1)**2 - (NMASS*0.001)**2)*DIRDUMMY(I)
                  GKIN(I,2) = -sqrt(GKIN(4,2)**2 - (NMASS*0.001)**2)*DIRDUMMY(I)
                  GPOS(I,1) = VECT(I)
                  GPOS(I,2) = VECT(I)
                end do
                IFLGK(1) = 0
                IFLGK(2) = 0
                TOFD(1) = 0.
                TOFD(2) = 0.

                NMEC=NMEC+1
                LMEC(NMEC) = 12
              end if
            end if
#ifdef WithG4
          endif ! G4 or handmade pi- capture
#endif
         END IF
         GOTO 999
      ENDIF

      P=VECT(7)*1000.
      T=GEKIN*1000.
      IF(T.GT.380.) then
         write (*,*) 'SGNUCEFF: Error NEUT cross sections not valid for 
     &                tpi (',T,') > 380, setting TPI=380 MeV'
         stop
      end if

C     Get momentum vector
      DO I=1,3
         PVECT(I)=P*VECT(I+3)
      end do

      NGKINE=0

C-----Determine if interaction on Oxygen or Free Hydrogen
      
      B_OXYGEN=0
      
      NUCCHARGE=1
      if (IPART.eq.8) then
         PICHARGE=1
      else if (IPART.eq.9) then
         PICHARGE=-1
      end if

      call EFPNCRSH_PNSEP(P,PICHARGE,NUCCHARGE,SIGEL,SIGINEL,SIGCX,IERROR)

      SIGTOT=(SIGEL+SIGINEL+SIGCX)

      crshydrogen = SIGTOT*ROH/1000
      totcrswater = crshydrogen+PIOINTNEUT(P,IPART)*ROO

#ifdef DEBUG
      write(*,*) "pion P      = ", P
      write(*,*) "crshydrogen = ", crshydrogen
      write(*,*) "totcrswater = ", totcrswater
      write(*,*) "ratio       = ", crshydrogen/totcrswater
#endif 

      CALL GRNDM(R,1)
      if (R.lt.crshydrogen/totcrswater) then
         NUMBNDP=0
         NUMBNDN=0
         NUMFREP=1
         NUMATOM=1
         goto 200

      end if

#ifdef DEBUG
      write(*,*) "SGNUCEFF: Pion interaction on oxygen"
#endif
      B_OXYGEN=1

C-----Oxygen only for now (This will need to be modified if 
C-----calling this routine on rock or other materials)
      NUMBNDP=8
      NUMBNDN=8
      NUMFREP=2
      NUMATOM=16
      call nesetfgparams

C-----Get interaction cross sections
      Call getneutxs(P,tot,absb,inel,scx,dcx,elas)

      rxn=absb+inel+scx+dcx

C-----Calculate THE RATIO OF EACH INTERACTION
  
      PRXN=rxn/tot
      PELAS=elas/tot

#ifdef DEBUG
      write(*,*)"rxn, tot, PRXN, PELAS = ",rxn, tot, PRXN, PELAS
#endif

      PROBSUM=PRXN+PELAS
      IF (PROBSUM.ne.1.) then
         PRXN=PRXN/PROBSUM
         PELAS=PELAS/PROBSUM
      ENDIF

      PELAS=PRXN+PELAS

C-----GENERATE RANDOM NUMBER
C      R=RNDM()
      CALL GRNDM(R,1)
      IF(R.GE.0.AND.R.LT.PRXN) then
         GOTO 100
      ELSE IF(R.GE.PRXN) then
         GOTO 400
      end if

#ifdef DEBUG
      write(*,*) "SGNUCEFF Warning: Did not pick an interaction type
     $ within nucleus...assuming nothing happened."
#endif

      RETURN
C     
C---Reaction
C     
 100  ISTOP=1

C     Temporarily store VCWORK variables
      NVC_temp = NVC 
      do i=1,3
         POSVC_temp(i)    = POSVC(i)    
      enddo
      do IVC=1,NVC
         IPVC_temp(ivc)   = IPVC(ivc)    
         AMASVC_temp(ivc) = AMASVC(ivc)
         IORGVC_temp(ivc) = IORGVC(ivc)  
         IFLGVC_temp(ivc) = IFLGVC(ivc)  
         ICRNVC_temp(ivc) = ICRNVC(ivc)  
         TIMVC_temp(ivc)  = TIMVC(ivc)   
         IVTIVC_temp(ivc) = IVTIVC(ivc)  
         IVTFVC_temp(ivc) = IVTFVC(ivc)    
         
         do i=1,3   
            POSIVC_temp(i,IVC) = POSIVC(i,IVC) 
            POSFVC_temp(i,IVC) = POSFVC(i,IVC) 
            PVC_temp(i,IVC)    = PVC(i,IVC)    
         enddo
      enddo


C     Set number of VCWORK particle
      NVC=5

C     Bound nucleus for NEUT
      IBOUND=1

C     Maximum nuclear radius based on Woods-Saxon nuclear radius parameter
      RMAX = C * 2.5
 
C     Generate random position incident on nucleus
      NLOOPS=0

C     Rotate momentum vector to the z-axis
      PMOM_MAG=SQRT(PVECT(1)**2+PVECT(2)**2+PVECT(3)**2)
      COSTHETA=PVECT(3)/PMOM_MAG

      IF(ABS(COSTHETA).EQ.1)THEN
         SINTHETA=0
      ELSE
         SINTHETA=SQRT(1-COSTHETA**2)
      ENDIF

      IF(SINTHETA.GT.0.00001)THEN
         SINPHI=PVECT(2)/(PMOM_MAG*SINTHETA)
         COSPHI=PVECT(1)/(PMOM_MAG*SINTHETA)
      ELSE
         COSPHI=1.0
         SINPHI=0.0
      ENDIF

#ifdef DEBUG
      WRITE(*,*)"Initial PVECT = ",PVECT(1), PVECT(2), PVECT(3)
      write(*,*)"Rotating by angles = ", COSPHI, -SINPHI, COSTHETA, -SINTHETA
#endif

      CALL ROTATEVECT(PVECT,PVECTOUT,COSPHI,-SINPHI,COSTHETA,-SINTHETA)

#ifdef DEBUG
      write(*,*)"Rotated vector = ",PVECTOUT(1), PVECTOUT(2), PVECTOUT(3)
      write(*,*)"NVERT, NVCVERT prior to NEFFECT loop: ",NVERT, NVCVERT
#endif

C     Store FSI history variables in case NEFFECT needs to be re-done
      NVERT_TMP = NVERT
      NVCVERT_TMP = NVCVERT

 110  continue
C     Select random incident position in nucleus frame
C     (Should be same as in NEUT piscat simulation)
      CALL GRNDM(R,1)
      posv(1)=-RMAX+R*RMAX*2.
      CALL GRNDM(R,1)
      posv(2)=-RMAX+R*RMAX*2.
      RAD   =SQRT(posv(1)**2+posv(2)**2)
      IF(RAD.GE.RMAX)GO TO 110
      posv(3)=-SQRT(RMAX**2-RAD**2)+0.0001

C     Setup NEUT common blocks      
      DO I=1,3
         posnuc(I,5) = posv(I)
      end do

      IPVC(5)=IPDG(IPART)

      CALL UCOPY(PVECTOUT,PVC(1,5),3)
      IORGVC(5)=0
      IFLGVC(5)=0
      ICRNVC(5)=1
      CALL UCOPY(posv,POSIVC(1,5),3)
      IVTIVC(5)=0
      CALL MCMASS(IPVC(5),AMASVC(5))

#ifdef DEBUG
      write(*,*) "Calling NEFFECT with inputs..."
      write(*,*) "PDG = ",IPVC(5), " GEANT ID = ", IPART
      write(*,*) "SGNUCEFF pos: ", POSIVC(1,5), POSIVC(2,5), POSIVC(3,5)
      write(*,*) "SGNUCEFF pvc: ", PVC(1,5), PVC(2,5), PVC(3,5)
#endif

C     The NEUT cascade
C      write(*,*) "SGNUCEFF: Before NEFFECT, NVERT = ",NVERT
      CALL NEFFECT

#ifdef DEBUG
      write(*,*) "SGNUCEFF: Interaction flag = ", IFLGVC(5)
#endif
      
C     Loop until an interaction occurs
      if (IFLGVC(5).eq.0) then
         NLOOPS=NLOOPS+1
         
         NVERT = NVERT_TMP
         NVCVERT = NVCVERT_TMP

#ifdef DEBUG
         write(*,*) "Last position:",posnuc(1,5),posnuc(2,5),posnuc(3,5)
         write(*,*) "PDG  = ",IPVC(5)
         write(*,*) "P(3) = ",PVC(1,5),PVC(2,5),PVC(3,5)         
         write(*,*) "RMAX=",RMAX, ", NLOOPS = ",NLOOPS
         write(*,*)
#endif

C        Check to prevent infinite loop (seems to be unecessary)
c         if (NLOOPS.ge.500) then
C           write(*,*) "Error SGNUCEFF: Stuck in infinite loop"
C           write(*,*) "PDG = ",IPVC(5), "P = ",P, ", not interacting"
C           write(*,*) "Do nothing with this particle." 
C           RETURN
C         end if

C         write(*,*) "SGNUCEFF: No interaction, redoing NEFFECT"

         goto 110
      end if

#ifdef DEBUG
      write(*,*) "SGNUCEFF: Interaction flag = ", IFLGVC(5)
#endif
C      write(*,*) "SGNUCEFF: After NEFFECT, NVERT = ",NVERT111

C
C --  CONSIDER NULCEAR EXCITATION OF GAMMA
C         (Cannot interface at the moment)
C      CALL NEMKGAM ! only 16O
C
C --  nucleon rescattering
C
c      call nrintr
C      
C --  nucleon ejection following pion/delta absorption
C
C     Comment out for now since GEANT seems to enter an infinite
C     tracking loop with these nucleons (needs diagnosis and fix)
C      call absneject

C     Store interaction flag in secondaries stack
      iflgscnd_temp = IFLGVC(5)

C     Add placeholder proton for absorption events instead
      if (IFLGVC(5).eq.3) then
         
         NGKINE=1
         
         DO I=1,3
            GKIN(I,NGKINE)=0
            GPOS(I,NGKINE)=VECT(I)
         end do
         GKIN(5,NGKINE) = 14
         TOFD(NGKINE) = 0.
         IFLGK(NGKINE) = -1

         KCASE = NAMEC(12)
         
         NMEC = NMEC + 1
         LMEC(NMEC) = 12

      else
         
         KCASE = NAMEC(12)
         NMEC = NMEC + 1
         LMEC(NMEC) = 12

      end if

C -- Translate VCWORK to GEKIN (GEANT secondaries)
C
      do ivc=6,NVC

C        Skip particles that died within the nucleus
         if (ICRNVC(ivc).ne.1) cycle

C        Convert PDG to GEANT ID (see gt2pd.h)
         do GEANTID=1,49
            if (IPDG(GEANTID).eq.IPVC(ivc)) exit
         end do
         if (GEANTID.gt.49) then
            write(*,*) "SGNUCEFF Warning: Unknown PDG code = ",IPVC(ivc)
            write(*,*) "Ignoring this secondary particle."
            cycle
         end if

         NGKINE=NGKINE+1

C        Rotate back to incoming pion direction frame       
         CALL ROTATEVECT(PVC(1,ivc),PVECTOUT,COSPHI,SINPHI,
     $        COSTHETA,SINTHETA)

         do I=1,3
            GKIN(I,NGKINE)=PVECTOUT(I)/1000.
            GPOS(I,NGKINE)=VECT(I)
         end do

         GKIN(4,NGKINE) = sqrt(GKIN(1,NGKINE)**2+GKIN(2,NGKINE)**2+
     $        GKIN(3,NGKINE)**2+(AMASVC(IVC)/1000)**2)

         GKIN(5,NGKINE) = GEANTID
 
         TOFD(NGKINE) = 0.

#ifdef DEBUG
         write(*,*) "Secondary #",NGKINE
         write(*,*) "Geant ID = ",GKIN(5,NGKINE), ", PDG = ", IPVC(ivc)
         write(*,*) "P(3) GeV = ", GKIN(1,NGKINE), GKIN(2,NGKINE), GKIN(3,NGKINE)
         write(*,*) "E (GeV)  = ", GKIN(4,NGKINE)
         write(*,*) ""
#endif

      end do

      GOTO 998

C
C---Interaction on Free Hydrogen (or hadronic elastic scatter at 400)
C     

C     Pi-minus may charge exchange
 200  continue

#ifdef DEBUG
      write(*,*) "SGNUCEFF: Pion interaction on free hydrogen"

      write(*,*) "In free proton handling..."
#endif

      if (IPART.eq.9) then
         CALL GRNDM(R,1)
         
C     Charge exchange based on predicted CX cross section by SAID
C     fit of pi-N data (see Figure 2.5 in T2K-TN-033).
C     Beware: We are modelling inelastic scatter (pion production)
C             with SCX only, but should be ok since p_max=500 MeV/c
         if (R.lt.(SIGCX+SIGINEL)/SIGTOT) then

            ISTOP=1

C     Fill common block

C           Store a single vertex
            B_STORE_FREE = 0
            if (NVERT.lt.MAXVERT.and.NVCVERT.lt.MAXVCVERT-1) then

               B_STORE_FREE = 1
               
               NVERT = NVERT+1
               do I=1,3
                  POSVERT(I,NVERT) = 0
               end do

C     See fsihist.h for flag convention            
               IFLGVERT(NVERT) = 400
               iflgscnd_temp = IFLGVERT(NVERT)

C     Incoming particle information
               NVCVERT = NVCVERT + 1
               do I=1,3
                  DIRVERT(I,NVCVERT) = VECT(I+3)
               end do
               ABSPVERT(NVCVERT) = P
               ABSTPVERT(NVCVERT) = P
               IPVERT(NVCVERT) = IPDG(IPART)
               IVERTI(NVCVERT) = NVERT
               IVERTF(NVCVERT) = NVERT
            end if

C           Old SGPINT routine which may be too forward going               
C            CALL SGPICH(P,COST)
C            CONE = ACOS(COST)*180./PI
C            CALL RNCONE(VECT(4),CONE, DIR)

            if (NUMATOM.ne.1) then
               write(*,*) "SGNUCEFF Error: NUMATOM != 1 for pim free proton CX"
               stop
            end if

            ireaction=3
            CALL EFHINELS_PS(IPDG(IPART),PVECT,PDGID,PVECTOUT,1,0,ireaction)
C            call pd2gt(PDGID,GEANTID)
C            if (GEANTID.ne.7) then
C               write(*,*) "SGNUCEFF Error: Expected outgoing pi0 but returned PDGID=",PDGID
C               stop
C            end if

            NGKINE=1

            P = 0
            DO  I=1,3
               GKIN(I,NGKINE)=PVECTOUT(I)/1000.
               GPOS(I,NGKINE)=VECT(I)

               P = P + PVECTOUT(I)**2
            end do
            P = sqrt(P)
            
            DO  I=1,3
               DIR(I) = PVECTOUT(I)/P
            end do

            GKIN(5,NGKINE) = 7.
            TOFD(NGKINE) = 0.

            GKIN(4,NGKINE) = sqrt(GKIN(1,NGKINE)**2+GKIN(2,NGKINE)**2+
     $           GKIN(3,NGKINE)**2+(0.1349766)**2)

            KCASE = NAMEC(13)
            NMEC = NMEC + 1
            LMEC(NMEC) = 13

C     Fill common block

C           Outgoing particle information
            if (B_STORE_FREE.eq.1.and.NVCVERT.lt.MAXVCVERT) then

               NVCVERT = NVCVERT + 1
               do I=1,3
                  DIRVERT(I,NVCVERT) = DIR(I)
               end do
               ABSPVERT(NVCVERT) = P
               ABSTPVERT(NVCVERT) = P
               IPVERT(NVCVERT) = IPDG(GKIN(5,NGKINE))
               IVERTI(NVCVERT) = NVERT
               IVERTF(NVCVERT) = NVERT
            end if

#ifdef DEBUG
            write(*,*) "SGNUCEFF: pi- charge exchange"

            write(*,*) "Secondary #",NGKINE
            write(*,*) "Geant ID = ",GKIN(5,NGKINE), ", PDG = ", IPDG(GKIN(5,NGKINE))
            write(*,*) "P(3) GeV = ", GKIN(1,NGKINE), GKIN(2,NGKINE), GKIN(3,NGKINE)
            write(*,*) "E (GeV)  = ", GKIN(4,NGKINE)
            write(*,*) ""
#endif

            goto 999
            
         end if

      end if

C     If no pi-minus CX, then it scatters elastically
C     Come here for elastic scattering on oxygen too
 400  continue

#ifdef DEBUG
      write(*,*) "Calling SGPIEL..."
#endif

      ISTOP=1

C     Fill common block

C     Store a single vertex
      B_STORE_FREE = 0
      if (NVERT.lt.MAXVERT.and.NVCVERT.lt.MAXVCVERT-1) then
         B_STORE_FREE = 1

         NVERT = NVERT+1
         do I=1,3
            POSVERT(I,NVERT) = 0
         end do
         
C     See fsihist.h for flag convention            
         IFLGVERT(NVERT) = 900 + B_OXYGEN
         iflgscnd_temp = IFLGVERT(NVERT)
              
C     Incoming particle information
         NVCVERT = NVCVERT + 1
         do I=1,3
            DIRVERT(I,NVCVERT) = VECT(I+3)
         end do
         ABSPVERT(NVCVERT) = P
         ABSTPVERT(NVCVERT) = P
         IPVERT(NVCVERT) = IPDG(IPART)
         IVERTI(NVCVERT) = NVERT
         IVERTF(NVCVERT) = NVERT
      end if

      NGKINE=1
C     (Quasi-)elastic scatter on free proton
      if (NUMATOM.eq.1) then
         ireaction=1
C         write(*,*) "Before EFHINELS: ", PVECT(1), PVECT(2), PVECT(3)
         CALL EFHINELS_PS(IPDG(IPART),PVECT,PDGID,PVECTOUT,1,0,ireaction)

         P = 0
         DO  I=1,3
            GKIN(I,NGKINE)=PVECTOUT(I)/1000.

            P = P + PVECTOUT(I)**2
         end do
         P = sqrt(P)
         
         DO  I=1,3
            DIR(I) = PVECTOUT(I)/P
         end do

C         write(*,*) "After EFHINELS: ", GKIN(1,NGKINE)*1000, GKIN(2,NGKINE)*1000, GKIN(3,NGKINE)*1000

C     Elastic scatter on oxygen
      else
         
         CALL SGPIEL(P,COST)
      
C        Correct pion for energy lost and scattered angle
         CONE = ACOS(COST)*180./PI
         CALL RNCONE(VECT(4),CONE, DIR)
         DO  I=1,3
            GKIN(I,NGKINE)=DIR(I)*P/1000.
         end do
         
      end if

C     Set starting position of secondary
      DO  I=1,3
         GPOS(I,NGKINE)=VECT(I)
      end do

      GKIN(5,NGKINE) = IPART
      TOFD(NGKINE) = 0.

      GKIN(4,NGKINE) = sqrt(GKIN(1,NGKINE)**2+GKIN(2,NGKINE)**2+
     $     GKIN(3,NGKINE)**2+(0.13957018)**2)

      KCASE = NAMEC(13)
      NMEC = NMEC + 1
      LMEC(NMEC) = 13

C     Fill common block

C     Outgoing particle information
      if (B_STORE_FREE.eq.1.and.NVCVERT.lt.MAXVCVERT) then
         NVCVERT = NVCVERT + 1
         do I=1,3
            DIRVERT(I,NVCVERT) = DIR(I)
         end do
         ABSPVERT(NVCVERT) = P
         ABSTPVERT(NVCVERT) = P
         IPVERT(NVCVERT) = IPDG(GKIN(5,NGKINE))
         IVERTI(NVCVERT) = NVERT
         IVERTF(NVCVERT) = NVERT            
      end if


#ifdef DEBUG
      write(*,*) "SGNUCEFF: Pion elastic scatter"

      write(*,*) "Secondary #",NGKINE
      write(*,*) "Geant ID = ",GKIN(5,NGKINE), ", PDG = ", IPDG(GKIN(5,NGKINE))
      write(*,*) "P(3) GeV = ", GKIN(1,NGKINE), GKIN(2,NGKINE), GKIN(3,NGKINE)
      write(*,*) "E (GeV)  = ", GKIN(4,NGKINE)
      write(*,*) ""
#endif

      goto 999

C     Restore VCWORK bank after reactive type reaction (NEUT Cascade)
 998  NVC = NVC_temp 
      do i=1,3
         POSVC(i)    = POSVC_temp(i)    
      enddo
      do IVC=1,NVC
         IPVC(ivc)   = IPVC_temp(ivc)    
         AMASVC(ivc) = AMASVC_temp(ivc)
         IORGVC(ivc) = IORGVC_temp(ivc)  
         IFLGVC(ivc) = IFLGVC_temp(ivc)  
         ICRNVC(ivc) = ICRNVC_temp(ivc)  
         TIMVC(ivc)  = TIMVC_temp(ivc)   
         IVTIVC(ivc) = IVTIVC_temp(ivc)  
         IVTFVC(ivc) = IVTFVC_temp(ivc)    
         
         do i=1,3   
            POSIVC(i,IVC) = POSIVC_temp(i,IVC) 
            POSFVC(i,IVC) = POSFVC_temp(i,IVC) 
            PVC(i,IVC)    = PVC_temp(i,IVC)    
         enddo
      enddo
      
 999  RETURN
      end
#else

      write(*,*) 
     $     "Error: Trying to use SGNUCEFF with NEUT_ROOT undefined"
      stop
      end
#endif

C   this fuction is used in pi- capture written by Y.SONODA
C   the values are picked up in Fig.24 (purple hist), arXiv:1701.03628v2
C   contains 

      real function prob_radcapture (r)
        real r, m
        real g_energy(73), g_cdf(73)
        data g_energy/60.000000, 61.448800, 62.419900, 63.459600,
     &   64.494800, 65.417500, 66.100300, 67.009000, 67.987400,
     &   69.016500, 69.668900, 70.530500, 71.502900, 72.565400,
     &   73.633100, 74.490100, 75.462200, 76.552800, 77.477600,
     &   78.325200, 79.765000, 80.801700, 81.917300, 82.991700,
     &   84.269400, 85.625300, 86.593900, 87.486300, 88.295700,
     &   89.386500, 90.491200, 91.539200, 92.506700, 93.425000,
     &   94.544800, 95.440900, 96.438600, 97.503600, 98.375900,
     &   99.336400, 100.348000, 101.342000, 102.413000, 103.355000,
     &   104.385000, 105.388000, 106.270000, 107.254000, 108.067000,
     &   109.297000, 110.252000, 111.297000, 112.146000, 113.196000,
     &   114.185000, 115.094000, 116.164000, 117.207000, 118.108000,
     &   119.182000, 120.190000, 121.231000, 122.287000, 123.373000,
     &   124.325000, 125.273000, 126.247000, 127.195000, 128.340000,
     &   129.649000, 131.086000, 132.288000, 133.000000/
        data g_cdf/0.000000, 0.000477, 0.002952, 0.004158, 0.005393,
     &   0.007238, 0.009068, 0.009668, 0.010273, 0.011837, 0.013355,
     &   0.014031, 0.015696, 0.018583, 0.021028, 0.025319, 0.027867,
     &   0.030919, 0.035989, 0.040136, 0.045550, 0.050917, 0.057188,
     &   0.063540, 0.071866, 0.079711, 0.085857, 0.093147, 0.104554,
     &   0.115483, 0.124896, 0.136360, 0.148639, 0.164711, 0.177353,
     &   0.192614, 0.209640, 0.224536, 0.240913, 0.256249, 0.275275,
     &   0.298582, 0.318703, 0.343045, 0.366869, 0.395568, 0.418914,
     &   0.438482, 0.464736, 0.490598, 0.517616, 0.544408, 0.573926,
     &   0.602792, 0.642383, 0.689892, 0.741823, 0.783810, 0.830299,
     &   0.885076, 0.922332, 0.935304, 0.947855, 0.968720, 0.985641,
     &   0.992852, 0.995085, 0.995393, 0.996262, 0.996649, 0.997494,
     &   0.997891, 1.000000/

        do I=1,72
          if ( r .ge. g_cdf(I) ) then
            m = (r - g_cdf(I)) / (g_cdf(I+1) - g_cdf(I))
            prob_radcapture = (m * g_energy(I+1) + (1 - m) * g_energy(I))
          end if
        end do
        return
      end function prob_radcapture

      END


