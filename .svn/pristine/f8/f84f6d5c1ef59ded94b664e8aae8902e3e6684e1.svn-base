************************************************************************
*     -------------------------------------
      SUBROUTINE DSRES
*     -------------------------------------
*
*     (Purpose)
*        Assemble the individual photons' information.
*        Consider ADC and TDC resolutions
*
*     (Input)
*        'MCTQ' bank
*
*     (Output)
*        'DATATMP' bank (deleted in "dstrig.F")
*       
*     (Creation Date and Author)
*       1992.05.21 ; First version by K.S.Hirata
*       1994.02.18 ; Modified for trigger simulation by K. S. Hirata
*                    (Ring informations are lost.)
*       1995.02.03 ; Modified by Koshio
*       2003.01.07 ; Modified for SK-II by Koshio
*       2007.01.23 ; Modified for SK-III by Koshio
*
************************************************************************

      implicit none

#include "skparm.h"
#include "skhead.h"
#include "skdetsim_parm.h"
#include "dsparm.h"
#include "skmc.h"
#include "dsdata.h"
#include "skpmt.h"
      
      integer IBUF(3)
      real RBUF(3)
      EQUIVALENCE (IBUF(1),RBUF(1))
      integer ICLIST(MXPHDS),IQLIST(MXPHDS),IFLAG(MAXPM)
      real QLIST(MXPHDS),TLIST(MXPHDS)

      integer icnt,iring,iph,ndata,ic,iq,ibef,ierr,i,ifthr,icab,
     &     nwmast,nsegm,lsegm,isseg,isegm
      real t,dt,absdt,pe1,pe,thr,tsig,trnd,qpc
      real tsigE, tsigL

C --- Consider T & Q resolutions

c      write(*,*) '*******************************'
c      write(*,*) '*******************************'
c      write(*,*) '******   DSRES is called ******'
c      write(*,*) '*******************************'
c      write(*,*) '*******************************'
      ICNT= 0
      CALL VZERO(IFLAG,MAXPM)	
      CALL VZERO(QLIST,MAXPM)	
      
      DO 10 IRING=1,NRNGDS+1
         
         DO 20 IPH=1,NPHRDS(IRING)
            
            CALL KZGET2('MCTQ',IRING,IPH,NDATA,IBUF)

            IC = IBUF(1)
            T  = RBUF(2)
            IQ = IBUF(3)

            IF(IFLAG(IC) .EQ. 0) THEN

               ICNT = ICNT+1
               ICLIST(ICNT) = IC
               TLIST(ICNT) = T
               IQLIST(ICNT) = IQ
               IFLAG(IC) = ICNT

            ELSE

               IBEF = IFLAG(IC)
               	
               DT = T - TLIST(IBEF)
               ABSDT = ABS(DT)

               IF(ABSDT .LT. SIGWDS) THEN
                  IQLIST(IBEF) = IQLIST(IBEF) + IQ
                  IF(DT .LT. 0.) TLIST(IBEF) = T

               ELSE IF(ABSDT.GE.SIGWDS.AND.ABSDT.LE.DEADDS) THEN
C                  PRINT *,' === ATM veto ===',IC,ABSDT 
                  IF(DT .LT. 0.) THEN
                     IQLIST(IBEF) = IQ
                     TLIST(IBEF) = T
                  ENDIF

*** Not supposed yet (for multi event)
               ELSE
C                  PRINT *,' === 2 buffers are full. === ',IC,ABSDT

               END IF
            END IF

 20      CONTINUE

 10   CONTINUE   

C --- Save data into a temporary bank

      NWMAST = 0
      NSEGM = 1
      LSEGM = 3
      CALL KZBCR2('DATATMP',NWMAST,'*I',NSEGM,LSEGM,'1I 2F',IERR)

      ISSEG = 0
      DO 30 I=1,ICNT	

         If (sk_geometry.eq.SK_I) then
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     single p.e. generation for SK-I
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
            DO 35 IQ=1,IQLIST(I)
               CALL SKRN1PE(PE1)
               QLIST(I)=QLIST(I)+PE1
 35         CONTINUE
            PE = QLIST(I)
            
            if((DS_TUNE .eq. LE_TUNE) .or. (DS_TUNE .eq. UNI_TUNE)) then
               thr=0.32
               CALL DSTHR(PE,IFTHR,thr)
            else if(DS_TUNE .eq. AP_TUNE) then
               CALL DSTHR(PE,IFTHR,thr)
            else
               write(*,*) 'DSRES: Not Supported SK_GEOMETRY',
     &              SK_GEOMETRY,' and DS_TUNE',DS_TUNE
               stop
            endif
            
            IF(IFTHR.EQ.1) GOTO 30
            if(
     & ((DS_TUNE.eq.LE_TUNE) .or. (DS_TUNE .eq. UNI_TUNE)) .and.
     &  (iqlist(i).eq.1) ) pe=pe*0.935
            QLIST(I) = PE
            
            if((DS_TUNE .eq. LE_TUNE) .or. (DS_TUNE .eq. UNI_TUNE)) then
               call sktsig_le(qlist(i),tsig)
               tlist(i)=tlist(i)+4.4*exp(-7*pe**2)+0.71
               CALL rntwogauss(tlist(i),tsig,trnd)
               TLIST(I) = TRND
               
            else if(DS_TUNE.eq.AP_TUNE) then
               qpc = qlist(i) * 2.0555
               call sktsig(qpc,tsig)
c     CALL SKTSIG(QLIST(I),TSIG)
               CALL RNGAUS(0.,TSIG,TRND)
               TLIST(I) = TLIST(I)+TRND
            else
               write(*,*) 'DSRES: Not Supported DS_TUNE',DS_TUNE
               stop
            endif
         Else if (sk_geometry.eq.SK_II) then
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     single p.e. generation for SK-II
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     only for lowe analyses (highe uses dshigh instead of this subroutine)
            DO 36 IQ=1,IQLIST(I)
               CALL SKRN1PE_SK2(PE1)
               QLIST(I)=QLIST(I)+PE1
 36         CONTINUE
            PE = QLIST(I)
            
            CALL DSTHR(PE,IFTHR,thr)
            
            IF(IFTHR.EQ.1) GOTO 30
            QLIST(I) = PE
            
            call sktsig_le(qlist(i),tsig)
            tlist(i)=tlist(i)+4.4*exp(-7*pe**2)+0.71
            CALL rntwogauss(tlist(i),tsig,trnd)
            TLIST(I) = TRND
               
         Else if (sk_geometry.eq.SK_III) then
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     single p.e. generation for SK-III
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     only for lowe analyses (highe uses dshigh instead of this subroutine)
            icab = iclist(i)

            do IQ=1,IQLIST(I)
               CALL SKRN1PE_SK3(icab,PE1)
               QLIST(I)=QLIST(I)+PE1
            enddo
            PE = QLIST(I)
            
            CALL dsthr_sk3(icab,PE,IFTHR)
            
            IF(IFTHR.EQ.1) GOTO 30
            QLIST(I) = PE

            call sktsig_le(qlist(i),tsig)
            tlist(i)=tlist(i)+4.4*exp(-7*pe**2)+0.71
            CALL rntwogauss(tlist(i),tsig,trnd)
            TLIST(I) = TRND
               
         Else if (sk_geometry.eq.SK_IV) then
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     single p.e. generation for SK-IV
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     only for lowe analyses (highe uses dshigh instead of this subroutine)
            icab = iclist(i)

            do IQ=1,IQLIST(I)
               CALL SKRN1PE_SK3(icab,PE1)
               QLIST(I)=QLIST(I)+PE1
            enddo
            PE = QLIST(I)
            
            CALL dsthr_sk3(icab,PE,IFTHR)
            
            IF(IFTHR.EQ.1) GOTO 30
            QLIST(I) = PE

            call sktsig_sk4_asy(qlist(i), tsigE, tsigL)
            tlist(i)=tlist(i)+4.4*exp(-7*pe**2)+0.71
            CALL rntwogauss_asy(tlist(i), tsigE, tsigL, trnd)
            TLIST(I) = TRND

         Else if (sk_geometry.eq.SK_V) then
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     single p.e. generation for SK-V (just copy from SK-IV) 
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     only for lowe analyses (highe uses dshigh instead of this subroutine)
            icab = iclist(i)

            do IQ=1,IQLIST(I)
               CALL SKRN1PE_SK3(icab,PE1)
               QLIST(I)=QLIST(I)+PE1
            enddo
            PE = QLIST(I)
            
            CALL dsthr_sk3(icab,PE,IFTHR)
            
            IF(IFTHR.EQ.1) GOTO 30
            QLIST(I) = PE

            if ( (nskpmt(icab) .eq. 2) .or. (nskpmt(icab) .eq. 3)) then
               call sktsig_sk5_asy(qlist(i), tsigE, tsigL)
            elseif (nskpmt(icab) .eq. 5)  then
               call sktsig_hk_asy(qlist(i), tsigE, tsigL)
            else
               write(6,*) 'DSRES: PMT category is strange. Stop.'
            endif
            
            tlist(i)=tlist(i)+4.4*exp(-7*pe**2)+0.71
            CALL rntwogauss_asy(tlist(i), tsigE, tsigL, trnd)
            TLIST(I) = TRND            
               
         Else
            write(6,*) 'supported only sk_geometry 1,2,3,4,5 in dsres'
            stop
         Endif

         IBUF(1) = ICLIST(I)
         RBUF(2) = TLIST(I)
         RBUF(3) = QLIST(I)

         ISEGM = 1
         NDATA = 3
         ISSEG = ISSEG + 1
         CALL KZREP2('DATATMP',NSEGM,ISSEG,NDATA,IBUF)

 30   CONTINUE   

C ---
      RETURN
      END
