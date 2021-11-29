************************************************************************
*     -------------------------------------
      SUBROUTINE DSHIGH
*     -------------------------------------
*
*     (Purpose)
*        Dark noise of PMT's
*        Channel conversion of T & Q
*
*     (Creation Date and Author)
*       1992.05.28 ; First version by K.S.Hirata
*       1994.02.19 ; Added for dark noise study by K.S.Hirata
*       1995.02.03 ; Modified by Koshio
*       1998.03.24 ; Modified by Ishihara
*       1998.10.18 ; Modified by Kaneyuki 
*       1999.05.23 ; remove random number reseeding (D.Casper)
*       2001. 2. 9 : RNDM -> GRNDM by K.K.
*       2001. 7. 2 : bug fix (RIDUMMY -> RDUMMY) by K.K.
*       2001. 8.23 : rn1pe -> aprn1pe
*       2003. 1.17 : aprn1pe -> skrn1pe
*       2003. 4.17 : skrn1pe_sk2 is added for SK-II MC by M.Shiozawa
*       2003. 7. 2 : Same seed for initalize always by C. Sterner
*       2004. 4.30 : Add FADC if defined 'FADC' by M.Ishitsuka
************************************************************************

#include "skparm.h"
#include "skhead.h"
#include "dsparm.h"
#include "dsparma.h"
#include "skmc.h"
#include "dsdata.h"
#include "dsoffset.h"
#include "skmisch.h"
#include "skdetsim_parm.h"

#include <zbs.h>

#ifdef FADC      
#include "fadcwrbnk.h"
#endif

C --  ADD SYSTEMATIC UNCERTAINTY K.K. 19-OCT-98
 
      REAL TSIGMA,QSIGMA
      PARAMETER(TSIGMA=0.7,QSIGMA=0.05)
      REAL TOFF(MAXPM),QOFF(MAXPM)
      INTEGER ICALL/0/


      DIMENSION JBUF(3),SBUF(3),KBUF(MAXRNG*2+3)
      EQUIVALENCE(JBUF(1),SBUF(1))
      DIMENSION IBUF(MAXPM*2),RBUF(MAXPM*2)
      EQUIVALENCE (IBUF(1),RBUF(1))
      DIMENSION ICTQ(3,MXPHDS),RCTQ(3,MXPHDS),IFLAG(MAXPM)
      EQUIVALENCE(ICTQ,RCTQ)

C-Y.H.-Support TQREAL BANK
      INTEGER*4 ITQBANK(3,MAXPM)
      REAL*4    RTQBANK(3,MAXPM)
      EQUIVALENCE(ITQBANK(1,1),RTQBANK(1,1))
      INTEGER*4 ITQRHD(3)
      REAL*4    RTQRHD(3)
      EQUIVALENCE (ITQRHD(1),RTQRHD(1))

C-Y.H.-Support TQREAL BANK

c     PARAMETER (OFFSET0 = 660.)
c     PARAMETER (OFFSET0 = 960.) ! it is moved to dsinip
      Integer offset, i, fasthit, start_j
      Integer tbuf(20000)

      Real satmean,satsigma
      Data satmean1/247.4/, satsigma1/7.662/
      Data satmean2/220.9/, satsigma2/14.16/

      real thr_dummy

C --  set saturation point
C     saturation charge of SK-II is changed mainly due to new pc2pe
C     changed by Shiozawa and data given by Nakayama
      If (SK_GEOMETRY.eq.SK_I) then
         satmean = satmean1
         satsigma = satsigma1
      Else If (SK_GEOMETRY.eq.SK_II) then
         satmean = satmean2
         satsigma = satsigma2
      Else
         satmean = satmean2
         satsigma = satsigma2
      Endif

C --  CALC SYSTEMATIC UNCERTAINTY 

      IF (ICALL.EQ.0) THEN
         ICALL=1
CCS      Store random number seed
         CALL RLUXAT(LUX,ISEED_RDM1,ISEED_RDM2,ISEED_RDM3)
CCS      Use same seeds for initialization every time, otherwise first
C           event differs and random seeds cannot be used to restart
         CALL RLUXGO(3,314159265,293928,0)
         DO 2 I=1,MAXPM
            Call rngausk(0.0,TSIGMA,TOFF(I))
            Call rngausk(1.0,QSIGMA,QOFF(I))
 2       CONTINUE
CCS      Restore random number seed
         CALL RLUXGO(LUX,ISEED_RDM1,ISEED_RDM2,ISEED_RDM3)
      ENDIF

      COMOFST=0.

C --- Read data in MCRAWTQ' bank

      isegm = 0
      issegm = 0
      CALL KZGET2('MCRAWTQ',isegm,issegm,NDATA,KBUF)

      NRNGDS = KBUF(1)
      NPHADS = 0
      DO 1 IRING=1,NRNGDS
         IPARDS(IRING) = KBUF(2*IRING)
         NPHRDS(IRING) = KBUF(2*IRING+1)
         NPHADS = NPHADS+NPHRDS(IRING)
 1    CONTINUE


C --- decide t0
      CALL VZERO(tbuf,20000)
      DO 7010 IRING=1,NRNGDS
         DO 7020 IPH=1,NPHRDS(IRING)
            CALL KZGET2('MCRAWTQ',IRING,IPH,NDATA,JBUF)
            T  = SBUF(2)
            If (t.lt.100000.0) then
               tbuf(int(t/5.0)+1) = tbuf(int(t/5.0)+1)+1
            Endif
 7020    Continue
 7010 Continue

      fasthit = 0
      Do 7030 i=1,20000
         start_j=i-40+1
         If (start_j.le.1) start_j=1
         fasthit = 0
         Do j=start_j,i
            fasthit = fasthit+tbuf(j)
         Enddo
         If (fasthit.gt.25.) then
            offset = offset0 - real(i)*5
            Goto 7040
         Endif
 7030 Continue
c      offset = offset0 - 500.  !-----1997-Jan-4 hayakawa
      offset = offset0
 7040 Continue

#ifdef WTSCTT
      offset = offset0
#endif

c      TMIN = offset - 300.
c      TMAX = offset + 1050.  !-----1997-Jan-4 hayakawa
C      TMIN = offset - 400.
C      TMAX = offset + 950.
C     changed by M.Shiozawa  May-8-97
      TMIN = offset0 - 400.
      TMAX = offset0 + 950.

      COMOFST=OFFSET
*     TRG 48clk. info.
      nt48sk(1) = 0
      nt48sk(2) = 0
      nt48sk(3) = ifix((offset0-offset)/20.)

      ICNT = 0
      CALL VZERO(IFLAG,MAXPM)

      DO 10 IRING=1,NRNGDS

         DO 20 IPH=1,NPHRDS(IRING)

            CALL KZGET2('MCRAWTQ',IRING,IPH,NDATA,JBUF)

            IC = JBUF(1)
ccc            T  = SBUF(2) + (GATEDS-GLTRDS)
            T  = SBUF(2) + OFFSET
            Q  = FLOAT(JBUF(3))

            IF(IFLAG(IC).EQ.0) THEN

C     if(T.gt.TMIN.or.T.lt.TMAX)then
C     changed by M.Shiozawa  May-8-97
               if(T.gt.TMIN.and.T.lt.TMAX)then
                  ICNT = ICNT+1
                  ICTQ(1,ICNT) = IC
                  RCTQ(2,ICNT) = T
                  RCTQ(3,ICNT) = Q
                  IFLAG(IC)    = ICNT
               end if

            ELSE

               IBEF = IFLAG(IC)
               TBEF = RCTQ(2,IBEF)
C     if(T.gt.TMIN.or.T.lt.TMAX)then
C     changed by M.Shiozawa  May-8-97
               if(T.gt.TMIN.and.T.lt.TMAX)then
                  IF(T.LT.TBEF) RCTQ(2,IBEF) = T
                  QBEF = RCTQ(3,IBEF)
                  RCTQ(3,IBEF) = QBEF + Q
               end if

            ENDIF

C            IF(ICNT.EQ.1) THEN
C               TMIN = T
C               TMAX = T
C            ELSE
C               IF(T.LE.TMIN) TMIN = T
C               IF(T.GE.TMAX) TMAX = T
C            ENDIF

 20      CONTINUE

 10   CONTINUE   
C
C     If bank is already present delete it before creating new one
C     
      CALL KZBLOC('TQ',IERR)
      IF (IERR.NE.0) CALL KZBDEL('TQ')
      CALL KZBCR0('TQ',IERR)
      NHIT = 0
      NWRD = 0

      IF(ICNT.GE.1) THEN
         CALL SORTI(RCTQ,3,ICNT,1)

      ELSE
         WRITE(6,*) 'NO HIT !!!'

      END IF

*     1998-Mar-24 comment out by K.Ishihara
*     ********************************************1997-Jan-4 Hayakawa *
*      DARKDS = 3.4E3 
*     ******************************************** 3.4 kHz ************
      RDARK = DARKDS * (TMAX-TMIN) * 1.E-09
      RDHIT = DARKDS * (GATEDS-GLTRDS-SIGWDS) * 1.E-09
      RNOHT = DARKDS * (DEADDS+GLTRDS-GATEDS) * 1.E-09

      J = 1
      DO 30 IC=1,MAXPM

*        skip not-existing PMTs

         if(imis(ic).ne.0) goto 30

         IF(IC.LT.ICTQ(1,J).OR.J.EQ.ICNT+1) THEN
*           -------------------------------no hit cable (Cal. dark noise)
C            IF(RNDM().GE.RDARK) GOTO 30
            CALL GRNDM(RDUMMY,1)
            IF(RDUMMY.GE.RDARK) GOTO 30
C            T = TMIN + RNDM()*(TMAX-TMIN)
            CALL GRNDM(RDUMMY,1)
            T = TMIN + RDUMMY*(TMAX-TMIN)

c           CALL APRN1PE(PE1)
            CALL RNDARK(PE1)
C    correction using dark (600-800) by K.Ishihara
            Q = PE1*1.13
         ELSE IF(IC.EQ.ICTQ(1,J)) THEN
*           ------------------------------- hit cable
            T = RCTQ(2,J)
            Q = RCTQ(3,J)
            J = J+1
C            IF(RNDM().LT.RNOHT) GOTO 30
            CALL GRNDM(RDUMMY,1)
            IF(RDUMMY.LT.RNOHT) GOTO 30
C            IF(RNDM().LT.RDHIT) THEN
            CALL GRNDM(RDUMMY,1)
            IF(RDUMMY.LT.RDHIT) THEN
*              ---------------------------noise hit is earlier than signal
               CALL RNDARK(PE1)
C    correction using dark (600-800) by K.Ishihara
               Q = PE1*1.13
C               T=T-SIGWDS-RNDM()*(GATEDS-GLTRDS-SIGWDS)
               CALL GRNDM(RDUMMY,1)
               T=T-SIGWDS-RDUMMY*(GATEDS-GLTRDS-SIGWDS)
            ELSE
               If (sk_geometry.eq.SK_I) then
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     single p.e. generation for SK-I
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
                  IQTMP = int(Q +.1)
                  Q=0.
                  do 35 ITMP=1,IQTMP
c     CALL APRN1PE(PE1)
                     CALL SKRN1PE(PE1)
                     Q=Q+PE1
 35               CONTINUE
C     IF(Q .LE. 0.) Q = 0.
                  
C     --           ADD K.K.
                  Q=Q*QOFF(IC)
                  
                  CALL SKTSIG(Q,TSIG)      
                  CALL RNGAUS(0.,TSIG,TRND)
                  T = T + TRND + TOFF(IC)
                  CALL DSTHR(Q,IFTHR,thr_dummy)
#ifdef TUNED_1PE
C     
C     charge dist. is modified for tuning MC of LINAC and decay-e
C     May-6-97  by M.Shiozawa
C     
                  q = q*0.95
#endif TUNED_1PE
                  IF(IFTHR.EQ.1) GOTO 30
               Elseif (sk_geometry.eq.SK_II) then
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     single p.e. generation for SK-II
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
                  IQTMP = int(Q +.1)
                  Q=0.
                  do 36 ITMP=1,IQTMP
                     CALL SKRN1PE_SK2(PE1)
                     Q=Q+PE1
 36               CONTINUE

                  Q=Q*QOFF(IC)
                  
                  CALL SKTSIG(Q,TSIG)      
                  CALL RNGAUS(0.,TSIG,TRND)
                  T = T + TRND + TOFF(IC)
                  CALL DSTHR(Q,IFTHR,thr_dummy)
                  IF(IFTHR.EQ.1) GOTO 30
               Else
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     single p.e. generation for SK-III
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
                  IQTMP = int(Q +.1)
                  Q=0.
                  do 37 ITMP=1,IQTMP
                     CALL SKRN1PE_SK3(IC,PE1)
                     Q=Q+PE1
 37               CONTINUE

                  Q=Q*QOFF(IC)
                  
                  CALL SKTSIG(Q,TSIG)      
                  CALL RNGAUS(0.,TSIG,TRND)
                  T = T + TRND + TOFF(IC)
                  CALL DSTHR_SK3(IC,Q,IFTHR)
                  IF(IFTHR.EQ.1) GOTO 30
               Endif
            ENDIF
         ELSE
            WRITE(6,*) 'ERROR in DSHIGH 2'
         ENDIF

** Fill TQbank

ccc         ITCNT = 1000+INT(TCNVSK*(GMAX - T)+0.5)
ccc         ITCNT = INT(T*TCNVSK+0.5)
ccc         IQCNT = INT(Q*QCNVSK+0.5)
            TCNVSK=2.0
         ITCNT = 4095 - INT(TCNVSK*T+0.5)
         IQCNT = INT(Q*QCNVSK+0.5) + 1000

ccc         IF(ITCNT.GE.4096) ITCNT=4096
ccc         IF(ITCNT.LE.0)    ITCNT=0
         
         IF(ITCNT.GE.4096.OR.ITCNT.LE.0) GOTO 30

         IF(IQCNT.GE.4096) IQCNT=4096
         IF(IQCNT.LE.0)    IQCNT=0

         NHIT = NHIT+1

C-Y.H.-Support TQREAL BANK
c-------------consider overflow
         if (Q.ge.230.) then
C            Call rdmin(real(IC)**3.1415)
            Call rngausk(satmean,satsigma,Q)
         Endif

         RTQBANK(2,NHIT)=T
         RTQBANK(3,NHIT)=Q
         ITQBANK(1,NHIT)=IC

         IF(IC .EQ. ICSV) THEN
             
            NWRD = NWRD+1
            IBUF(NWRD) = 2**29 + ITCNT*2**16 
     &                 + 2**14 + 2**13 + IQCNT

         ELSE IF(IC .EQ. (ICSV+1)) THEN

            NWRD = NWRD+1
            IBUF(NWRD) = 2**30 + ITCNT*2**16 + IQCNT

         ELSE

            NWRD = NWRD+1
            IBUF(NWRD) = IC
            NWRD = NWRD+1
            IBUF(NWRD) = ITCNT*2**16 + IQCNT

         END IF   
         ICSV = IC

 30   CONTINUE

      IF(NWRD.GT.MAXPM*2) THEN
         WRITE(6,*)'ERROR in DSHIGH '
         STOP
      ENDIF

      isegm = 1
      CALL KZREP0('TQ',isegm,'I',NWRD,IBUF)

      IBUF(1) = NHIT
      IBUF(2) = NWRD
      NDATA = 2
      isegm = 0
      CALL KZREP0('TQ',isegm,'I',NDATA,IBUF)
C
C     If bank is already present delete it before creating new one
C     
      CALL KZBLOC('TQREAL',IERR)
      IF (IERR.NE.0) CALL KZBDEL('TQREAL')
      nwmast = 3
      lsegm  = 3
      CALL KZBCR1('TQREAL',nwmast,'1I 1F 1I',NHIT,lsegm,'1I 2F',IERR)
      IF (IERR.ne.0) THEN
         write(6,*) 'ERROR in CREATING TQREAL BANK'
         RETURN
      ENDIF

      ITQRHD(1)=NHIT
      RTQRHD(2)=2.055
      ITQRHD(3)=1
      isegm = 0
      ndata = 3
      CALL KZREP1('TQREAL',isegm,ndata,ITQRHD)
      
      IF (NHIT.gt.0) THEN
C         write(*,*) "***************************  NOW COPYING"
         CALL KZLDAT('TQREAL',LDAT)
         NENT=3*NHIT
C         DO 1529 I=1,NHIT
C            write(*,*) "CABLE:",ITQBANK(1,I),"/T=",RTQBANK(2,I),
C     $                 "/Q=",RTQBANK(3,I)
C 1529    CONTINUE
         CALL LBFCOPY(ITQBANK,IZBS(LDAT+3+1),NENT)
      ENDIF

C dark noise for FADC only when you defined FADC
C Note: the dark noise for FADC is not identical with original.
#ifdef FADC
      if (IFADCI) then
         if (flag_fadc_sube .eq. 0) then
            call fadcdark
         end if
      endif
#endif

      RETURN
      END
