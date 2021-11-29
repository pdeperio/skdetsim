      SUBROUTINE APFILL2ND

      IMPLICIT NONE
#include "nework.h"
#include "vcwork.h"
#include "vcvrtx.h"
#include"skhead.h"
#include"apscndprt.h"
#include"apscndry.h"      

      Real wallsk
      
      INTEGER*4 IDATA(100),NDATA
      REAL*4    RDATA(100)
      EQUIVALENCE (IDATA(1),RDATA(1))

      INTEGER*4 lbank, i, j

      Integer tbuf(30000)
      Integer nsub
      Real p2ndabs

      numsubev = 0
      Call vzero(nsublist,(maxnsubev*maxnsubnev))
      Call vzero(nevinsube,maxnsubev)
      Call vzero(tbuf,30000)
      Call vzero(nevstartt,maxnsubev)

C------
C-Locate BANK
C      call kzdtab(6)
      CALL KZBLOC('CONVVECT',LBANK)
      IF (LBANK.eq.0) THEN
         write(*,*) "Event #",nevsk,":Doesn't contain CONVVECT bank"
         N2PRT=0
         NUMSUBEV=-1
         RETURN
      ENDIF

C-READ number of segments
      CALL KZGET1('CONVVECT',0,NDATA,IDATA)

      N2PRT=IDATA(1)
      
      IF (N2PRT.lt.0) THEN
         write(*,'(A17,$)') "Error in rd2ndprt:"
         write(*,*) "Number of secondary particles is less than 0!(",
     $              N2PRT,")"
         NUMSUBEV=-1
         RETURN
      ELSE
         IF (N2PRT.eq.0) THEN
C            write(*,*) "Event #",nevsk,":secondary particle is 0"
            NUMSUBEV=0
            RETURN
         ENDIF
      ENDIF 

C      write(*,*) "# of 2ndary particles:",N2PRT

C-READ informations of each 2ndary particle
      If (n2prt.gt.secmaxrng) n2prt = secmaxrng

      DO 10 I=1,n2prt
         CALL KZGET1('CONVVECT',I,NDATA,IDATA)
         IF (NDATA.ne.11) THEN
            write(0,'(A17,$)') "Error in rd2ndprt:"
            write(0,*) "Number of DATA is wrong(expected=11 !=",
     $                 NDATA,")"
            NUMSUBEV=-1
            RETURN
         ENDIF
            
         itrk(i)=IDATA(1)
         iprt2nd(i)=IDATA(2)
         DO 20 J=1,3
            vtx2nd(J,I)=RDATA(2+J)
            p2nd(J,I)  =RDATA(5+j)
 20      CONTINUE
         t2nd(I)   =RDATA(9)
         iprtorg(I)=IDATA(10)
         lmec(i)=IDATA(11)
         if (t2nd(I).lt.tingate) Goto 10

         p2ndabs = sqrt(p2nd(1,i)**2.0+p2nd(2,i)**2.0+p2nd(3,i)**2.0)
         If ((iabs(iprt2nd(i)).eq.11.or.iprt2nd(i).eq.22)
     &        .and.t2nd(i).lt.30000.
     &        .and.p2ndabs.gt.5.0
     &        .and.wallsk(rdata(3)).gt.0.) then
            tbuf(ifix(t2nd(i))) = tbuf(ifix(t2nd(i))) + 1
         Endif

 10   CONTINUE

      DO 60 I=1,30000
         If (tbuf(i).gt.0) then
            If (numsubev.eq.0) then
               numsubev = numsubev+1
               If (numsubev.gt.maxnsubev) then
                  print *, 'sub event exceed maximum.'
     &                 ,numsubev,maxnsubev
                  numsubev = maxnsubev
                  Goto 70
               Endif
               nevstartt(numsubev) = real(i)
            Else If (real(i)-nevstartt(numsubev).gt.tingate) then
               numsubev = numsubev+1
               If (numsubev.gt.maxnsubev) then
                  print *, 'sub event exceed maximum.'
     &                 ,numsubev,maxnsubev
                  numsubev = maxnsubev
                  Goto 70
               Endif
               nevstartt(numsubev) = real(i)
            Else
*     no new sub event.  this is within gate with previous event.
            Endif
         Endif
 60   CONTINUE
 70   CONTINUE

*      print *, 'numsubev=', numsubev
      Do 80 I=1,n2prt
         Do 81 nsub=1,numsubev
            If (t2nd(i).ge.nevstartt(nsub)
     &           .and.t2nd(i).lt.nevstartt(nsub)+tingate) then
               nevinsube(nsub) = nevinsube(nsub) + 1
               If (nevinsube(nsub).gt.maxnsubnev) then
                  print *, '#event in sub event gate exceed maximum.'
     &                 ,nevinsube(nsub),maxnsubnev
                  nevinsube(nsub) = maxnsubnev
                  Goto 82
               Endif
               nsublist(nsub,nevinsube(nsub)) = i
               Goto 82
            Endif
 81      Continue
 82      Continue
 80   Continue


      Do 90 nsub=1,numsubev
         print *, ' SUBE',nsub,': nev=',nevinsube(nsub)
         Do 91,i=1,nevinsube(nsub)
            print *, '         PID,DT,MOM='
     &           ,iprt2nd(nsublist(nsub,i)),t2nd(nsublist(nsub,i))
     &           ,sqrt(p2nd(1,nsublist(nsub,i))**2.0
     &           +p2nd(2,nsublist(nsub,i))**2.0
     &           +p2nd(3,nsublist(nsub,i))**2.0)
 91      Continue
 90   Continue


C- 1998/03/16 Add Y.H  to read new bank

C-Locate BANK
      CALL KZBLOC('SCNDPRT',LBANK)
      IF (LBANK.eq.0) THEN
         write(*,*) "Event #",nevsk,":Doesn't contain SCNDPRT bank"
         RETURN
      ENDIF

C-READ number of segments
      CALL KZGET1('SCNDPRT',0,NDATA,IDATA)
      nscndprt = IDATA(1)

      DO 100 I=1,nscndprt
         CALL KZGET1('SCNDPRT',I,NDATA,IDATA)
         IF (NDATA.ne.13 .and. NDATA.ne.17) THEN
            write(0,'(A17,$)') "Error in apfill2nd:"
            write(0,*) "Number of DATA is wrong(expected=13 or 17 !=",
     $                 NDATA,")"
            RETURN
         ENDIF
         
         IF (NDATA.eq.13) THEN
            itrkscnd(i)=IDATA(1)
            iprtscnd(i)=IDATA(2)
            DO J=1,3
               vtxscnd(J,I)=RDATA(2+J)
               pscnd(J,I)  =RDATA(5+j)
               pprnt(J,I)  =0
            END DO
            tscnd(I)   =RDATA(9)
            iprntprt(I)=IDATA(10)
            lmecscnd(i)=IDATA(11)
            iprnttrk(i)=idata(12)
            iorgprt(i)=idata(13)
            iflgscnd(i)=0

         ELSE IF (NDATA.eq.17) THEN
            itrkscnd(i)=IDATA(1)
            iprtscnd(i)=IDATA(2)
            DO J=1,3
               vtxscnd(J,I)=RDATA(2+J)
               pscnd(J,I)  =RDATA(5+j)
               pprnt(J,I)  =RDATA(8+j)
            end do
            tscnd(I)   =RDATA(12)
            iprntprt(I)=IDATA(13)
            lmecscnd(i)=IDATA(14)
            iprnttrk(i)=idata(15)
            iorgprt(i)=idata(16)
            iflgscnd(i)=idata(17)
         end if
 100  continue

      return
      END

