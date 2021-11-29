      SUBROUTINE SKFAKEHD

#include "skhead.h"
#include "geant321/gcflag.inc"
#include "mcevent.h"

CMS      IMPLICIT NONE

      INTEGER*4 NDATA,IDATA(1000),IDATA1(1000)
      INTEGER*4 ITMP(3)
      INTEGER*4 IERR,LBANK
	  INTEGER*4 ICOUNT, OSUB
      CHARACTER*8 CD,CT,DUM1,DUM2

C K.Ishihara 8/20/98
	  DATA ICCOUNT/0/
	  DATA OSUB/0/

      CALL KZBLOC('HEADER',LBANK)
      IF (LBANK.eq.0) GOTO 150

      isegm = 0
      CALL KZGET0('HEADER',isegm,NDATA,IDATA1)
      isegm = 1
      CALL KZGET0('HEADER',isegm,NDATA,IDATA)
      CALL KZBDEL('HEADER')

 150  CALL KZBCR0('HEADER',IERR)
      IF (IERR.ne.0) THEN
         write(*,*) "Could not create HEADER BANK"
      ENDIF

      IDATA1(1)=4
      IDATA1(2)=2
      IDATA1(3)=12
      IDATA1(4)=58
      isegm = 0
      ndata = 4

      CALL KZREP0('HEADER',isegm,'I',ndata,IDATA1)
      
      IDATA(1)=0
      nrunsk=999999
      IDATA(2)=nrunsk
C K.Ishihara 8/20/98
      IDATA(3)=NSUBSK
      IEV=IEV+1
      IDATA(4)=IEV


      IF (LBANK.eq.0) THEN
C     K.Ishihara 29/10/98
C     for cyclic subrun
C     IF (OSUB.NE.NSUBSK) ICOUNT=0
         OSUB=NSUBSK
         ICOUNT=ICOUNT+1
         IDATA(4)=ICOUNT
         
         CALL DATIMH(CD,CT)
         
         READ(CD,100) ITMP(3),DUM1,ITMP(2),DUM2,ITMP(1)
         IDATA(5)=ITMP(1)*10000+ITMP(2)*100+ITMP(3)
         
         READ(CT,100) ITMP(1),DUM1,ITMP(2),DUM2,ITMP(3)
         IDATA(6)=ITMP(1)*10000+ITMP(2)*100+ITMP(3)
      ENDIF    

*      IDATA(7)=0
*      IDATA(8)=0
*      IDATA(9)=0
      IDATA(7)=nt48sk(1)
      IDATA(8)=nt48sk(2)
      IDATA(9)=nt48sk(3)
      IDATA(10)=3
* Set bits 0, 1, 31 (logically rather than arithmetically to avoid overflow error)
      IDATA(12) = IBSET(IBSET(IBSET(0,0),1),31)
C      IDATA(12)=1+2+32768*65536
      isegm = 1
      lsegm = 12

      CALL KZREP0('HEADER',isegm,'I',lsegm,IDATA)

      RETURN
 100  FORMAT(I2,A1,I2,A1,I2)

      END


