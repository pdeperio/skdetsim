***********************************************************************
*     -----------------
      SUBROUTINE SKTEST
*     -----------------
*
*     These routines are for the standard diagnostic histogram package
*     that has been added to skdetsim.  sktest_init initializes and books
*     histograms and ntuples.  sktest_close closes these HBOOK items and
*     prints a summarey of the run just completed for quick reference.
*
*     2003.2 ; C. Sterner First version 
*     2006.4.18 ; K Bays updated for SK-III
***********************************************************************

      implicit none
#include "skparm.h"
#include "skdetsim_parm.h"
#include "skdetsim_test.h"
#include "skhead.h"
#include "guxyzpm.h"

      INTEGER IMVTX,IMDIR,IPART,IMMOM
      REAL RVTX, RWALL, RDIR, RMOM
      INTEGER RNDCARD
      COMMON/VECTCARD/IMVTX,RVTX(3),RWALL,IMDIR,RDIR(3),
     &     IPART,IMMOM,RMOM(2),RNDCARD(5)
      Integer nhbk,ii, lunout, ierr, iiexist, ll, nn, order, ordere
      integer TEST_numpmt
      parameter (nhbk=1000000)
      integer h
      common /pawc/ h(nhbk)
      real hstati
      character*600 hbktag
      real distvtx

      entry sktest_init
         call hlimit(-nhbk)         
         lunout = 89
         ierr = 0
         call hropen(lunout,'MCTEST',TEST_fname,'N',1024,ierr)
         if (ierr.ne.0) then
            write(6,*)'cannot open hbook file'
            stop
         endif
         call hbnt(22,'MCTEST1','D')
         hbktag = 'event:I, Gen:I , Coll:I , Colla:I '
         call hbname(22,'ALL',TEST_EVENTS,hbktag)
         call hbnt(23,'MCTESTPMT','D')
         hbktag = 'nopmt:I, x:R, y:R, z:R, exist:I '
         call hbname(23,'TESTPMT',TEST_NOPMT,hbktag)

C     Book histograms (control bounds for some standard event types)
         if(IPART.eq.11.and.RMOM(1).eq.10.0) then !electron
               CALL HBOOK1(40,'Number of Photons Generated'
     &              ,100,0.,2500.,0)
               if (sk_geometry.eq.SK_I) then
               CALL HBOOK1(41,'Number of Photons Collected in ID'
     &              ,100,0.,200.,0)
               else if (sk_geometry.eq.SK_II) then
                  CALL HBOOK1(41,'Number of Photons Collected in ID'
     &              ,100,0.,100.,0)
               else if (sk_geometry.eq.SK_III) then
                  CALL HBOOK1(41,'Number of Photons Collected in ID'
     &              ,100,0.,200.,0)
               endif   
         else if(IPART.eq.13.and.RMOM(1).eq.1000.) then !muon
               CALL HBOOK1(40,'Number of Photons Generated',400
     &              ,20000.,40000.,0)
               if (sk_geometry.eq.SK_I) then
               CALL HBOOK1(41,'Number of Photons Collected in ID'
     &              ,400,5200.,11600.,0)
               else if (sk_geometry.eq.SK_II) then
                  CALL HBOOK1(41,'Number of Photons Collected in ID'
     &              ,400,2600.,5800.,0)
               else if (sk_geometry.eq.SK_III) then
                  CALL HBOOK1(41,'Number of Photons Collected in ID'
     &              ,400,5200.,11600.,0)
               endif   
             
         else !use automatic bounds and bins
            CALL HBOOK1(40,'Number of Photons Generated',100,0.,-20.,0)
            CALL HBOOK1(41,'Number of Photons Collected in ID',100,0.,-20.,0)
         endif
         CALL HBOOK1(42,'Number of Photons Collected in OD',100,0.,-20.,0)
         CALL HBOOK1(43,'Number of ID PMT Hit',100,0.,-20.,0)
C     Hitmaps
         CALL HBOOK1(50,'Hitmap of all ID PMTs',MAXPM,1.,MAXPM,0)
         CALL HBOOK1(51,'Hitmap of only existing ID PMTs',MAXPM,1.,MAXPM,0)
         CALL HBOOK1(60,'Hitmap of OD PMTs',MAXPMA,1.,MAXPMA,0)
         if(IMVTX.eq.1) then
         CALL HBOOK1(52,'Number of Photons Collected vs Distance from Vertex'
     &        ,MAXPM,0.,2500.0,0)
         CALL HBOOK1(53,'Number of ID PMT Hit vs Distance from Vertex'
     &        ,MAXPM,0.,2500.0,0)
         CALL HBOOK1(54,'Number of ID PMT That exist vs Distance from Vertex'
     &        ,MAXPM,0.,2500.0,0)
         endif
         CALL HBOOK1(55,'Hitmap in Horizontal Bands of All ID Wall',
     &        MAXPM,1.,MAXPM,0)
         CALL HBOOK1(56,'Hitmap in Horizontal Bands of Existing ID Wall',
     &        MAXPM,1.,MAXPM,0)

         DO ii=1,MAXPM
            TEST_pmtmap(ii)=0.
         ENDDO
      
         return
         
         entry sktest_close

         if(TEST_Events.ne.0) then
            iiexist=0
            TEST_numpmt =0
            DO ii=1,MAXPM
               IF(TEST_Hitmap(ii).gt.0.0) then
                  TEST_numpmt=TEST_numpmt+1
               endif
               TEST_sum_coll=TEST_sum_coll+TEST_Hitmap(ii)
               TEST_Hitmap(ii)=TEST_Hitmap(ii)/TEST_Events
               if(TEST_pmtmap(ii).eq.1) then
                  iiexist=iiexist+1
                  call HFILL(51,FLOAT(iiexist),0.,TEST_Hitmap(ii))
                  if(IMVTX.eq.1) then
                     distvtx=sqrt((gxyzpm(1,ii)-RVTX(1))**2
     &                    + (gxyzpm(2,ii)-RVTX(2))**2
     &                    + (gxyzpm(3,ii)-RVTX(3))**2)
                     call HFILL(52,distvtx,0.,TEST_Hitmap(ii))
                     TEST_NHitmap(ii)=TEST_NHitmap(ii)/TEST_Events
                     call HFILL(53,distvtx,0.,TEST_Nhitmap(ii))
                     call HFILL(54,distvtx,0.,1.0)
                  endif
               endif
            enddo
            DO ii=1,MAXPMA
               TEST_sum_colla=TEST_sum_colla+TEST_Hitmapa(ii)
               TEST_hitmapa(ii)=TEST_Hitmapa(ii)/TEST_Events
            ENDDO
C       Hitmaps with horizontal not vertical bands
            DO ll=0,50
               nn=ll
               DO ii=1,150
                  nn=nn+1
                  order=order+1.
                  CALL HFILL(55,FLOAT(order),0.,TEST_Hitmap(nn))
                  if(TEST_pmtmap(nn).eq.1) then
                     ordere=ordere+1.
                     CALL HFILL(56,FLOAT(ordere),0.,TEST_Hitmap(nn))
                  endif
                  nn=nn+50
               ENDDO
            ENDDO
            DO order=7651,MAXPM
               nn=nn+1
               CALL HFILL(55,FLOAT(order),0.,TEST_Hitmap(nn))
               if(TEST_pmtmap(nn).eq.1) then
                  ordere=ordere+1.
                  CALL HFILL(56,FLOAT(ordere),0.,TEST_Hitmap(nn))
               endif
            ENDDO

               
         endif
         call HPAK(50,TEST_Hitmap)
         call HPAK(60,TEST_Hitmapa)
         call hrout(0,ierr,' ')
C
C+ Output summary of Diagnostic
C     
         write(6,*) '*******************************************'
         write(6,*) '*    SKDETSIM internal diagnostic summary *'
         write(6,*) '*******************************************' 
         write(6,*) 'Diagnostic Histograms output to: ',TEST_fname
         write(6,*) 'Events Processed:', TEST_Events
         Write(6,*) 'Number of PMTs counted :by gustep =',iiexist
     &        ,' by sgpmt =',TEST_numpmt
         write(6,*) 'Mean Number of Cherenkov Photons Per Event:',
     &        HSTATI(40,1,' ',1)
         write(6,*) 'Mean Number of Photons Collected in ID per Event:',
     &        HSTATI(41,1,' ',1)
         write(6,*) 'Mean Number of Photons Collected in OD per Event:',
     &        HSTATI(42,1,' ',1)

         write(6,*) 'Average Number of Photons Collected Per ID PMT:',
     &       FLOAT( TEST_sum_coll)/FLOAT(iiexist)
         write(6,*) 'Average Number of Photons Collected Per OD PMT:',
     &         FLOAT(TEST_sum_colla)/FLOAT(MAXPMA)
         write(6,*) ' '
         write(6,*) '*******************************************' 


         call hrend('MCTEST')
         close(lunout)
         
         return
         end
