**********************************************************************
*     --------------------
      SUBROUTINE FADCWRITE
*     --------------------
*     (Creation Date and Author)
*     2001.12.27
*     2002.09.02 ; check hit channels by threshold cut ... T.Ooyabu
*     2003.05.27 ; change the bank format to consistent with observed data
*                                                          M.Ishitsuka
*     2003.08.25 ; debug last some signal T.Harada
*     2004.04.30 ; Arranged and put on CVS by M.Ishitsuka
**********************************************************************
      IMPLICIT NONE

#include "skhead.h"
#include "skparm.h"
#include "fadcparm.h"
#include "fadcbnk.h"
#include "fadc1pe.h"
#include "fadcwrbnk.h"

      Integer threshold_flag
C     Compress
      Integer i,j,m
      Integer hit_flag
      Integer pulse_cnt,forpls(7)
C     static variables for bank 
      Integer ii, nbit, nqq
      Integer ierr
      integer charch(NFADC_RANGE*2)
      integer ninttot, nsizetot, nsizech

      Call VZERO(ibufh,nseg0) 
      Call VZERO(ibuf1,nseg1) 
      Call VZERO(ibuf2,nseg2) 
      Call VZERO(ibuf3,nseg3) 

C*********************************************************************
C     check hit channels, threshold cut and number of hits
C*********************************************************************
      ibuf1(1)=0
      ibuf2(1)=0
      nsizetot=0
      ninttot=0
      DO i=1,MAXFADC
         threshold_flag = 0
         do j=1,NFADC_RANGE 
            pulse_cnt = int(RFCHPLS(j,i)*count2mv)
            if (pulse_cnt.ge.nthreshold) then
               threshold_flag = 1
            endif
         end do
         if(threshold_flag.eq.0) then
            ibuf1(i+1)=0
            ibuf2(i+1)=0
         else
            nsizech=10          !head
            charch(1)=0
            charch(5)=ibset(charch(5),5)
            charch(7)=ibset(charch(7),5)
            do nbit=0,7
               if(btest(nevsk,nbit+8)) then
                  charch(6) = ibset(charch(6),nbit)
               end if
               if(btest(nevsk,nbit)) then
                  charch(8) = ibset(charch(8),nbit)
               end if
            end do
            charch(9) = ibset(charch(9),6)
            do nbit=0,7
               if(btest(npedestal,nbit)) then
                  charch(10) = ibset(charch(10),nbit)
               end if
            end do            
            do j=1,NFADC_RANGE 
               Call VZERO(forpls,7)
               do m = 1,7
                  forpls(m) = int(RFCHPLS(j+m-1,i)*count2mv)
               enddo
               pulse_cnt = int(RFCHPLS(j,i)*count2mv)
               If ( pulse_cnt .gt. npedestal ) then
                  if ( hit_flag .eq. 0 ) then
                     ibuf1(i+1)=ibuf1(i+1)+1                     
                     charch(nsizech+1) = ibset(charch(nsizech+1),5)
                     charch(nsizech+1) = ibset(charch(nsizech+1),6)
                     do nbit=0,7
                        if(btest(NFADC_RANGE-j,nbit+8)) then
                           charch(nsizech+1)=ibset(charch(nsizech+1),nbit)
                        end if
                        if(btest(NFADC_RANGE-j,nbit)) then
                           charch(nsizech+2)=ibset(charch(nsizech+2),nbit)
                        end if
                     end do
C 818 debug
                     If ( charch(nsizech+1).eq.96 .and.
     $                    charch(nsizech+2).eq.0 ) then
                        charch(nsizech+1)=0
                        charch(nsizech+2)=0
                        goto 33
                     end if
                     
                     If ( charch(nsizech+1).eq.96 .and.
     $                    charch(nsizech+2).eq.1 ) then
                        charch(nsizech+1)=0
                        charch(nsizech+2)=0
                        goto 33
                     end if
                     
                     If ( charch(nsizech+1).eq.96 .and.
     $                    charch(nsizech+2).eq.2 ) then
                        charch(nsizech+1)=0
                        charch(nsizech+2)=0
                        goto 33
                     end if
                           
                     If ( charch(nsizech+1).eq.96 .and.
     $                    charch(nsizech+2).eq.3 ) then
                        charch(nsizech+1)=0
                        charch(nsizech+2)=0
                        goto 33 
                     end if
                                 
C     818 debug end
                                 
                     nsizech=nsizech+2 !hit position
                     nqq = 0
                     do ii=1,3
                        nqq = nqq + 1 
                        charch(nsizech+2+nqq)
     $                       = int(RFCHPLS(j+nqq-4,i)*count2mv)
                     enddo
                  endif
                  nqq=nqq+1
                  charch(nsizech+2+nqq)
     $                 = int(RFCHPLS(j,i)*count2mv)
                  hit_flag=1
                  
C     820 debug start
                  if(j.eq.NFADC_RANGE) then
                     nsizech = nsizech - 2
                     
                     if(charch(nsizech-1).eq.0 .and.
     $                    charch(nsizech).eq.0  ) then
                        charch(nsizech+1)=0
                        charch(nsizech+2)=0
                        charch(nsizech+3)=0
                        charch(nsizech+4)=0
                        charch(nsizech+5)=0
                        charch(nsizech+6)=0
                        charch(nsizech+7)=0
                        goto 33
                     end if
 
                     if(charch(nsizech-1).eq.64 .and.
     $                    charch(nsizech).eq.0  ) then
                        do ii=1,18
                           charch(ii)=0
                        end do
                        nsizech=0
                        goto 55
                     end if
                     
                  end if
C     820 debug end
                  
               else
                  if ( hit_flag .eq. 1 ) then
                     if ( forpls(1).le.npedestal .and. 
     $                    forpls(2).le.npedestal .and. 
     $                    forpls(3).le.npedestal .and. 
     $                    forpls(4).le.npedestal .and. 
     $                    forpls(5).le.npedestal .and. 
     $                    forpls(6).le.npedestal .and.
     $                    forpls(7).le.npedestal ) then
                        do ii=1,3
                           nqq=nqq+1
                           charch(nsizech+2+nqq)
     $                          = int(RFCHPLS(j+ii-1,i)*count2mv)
                        enddo
                        charch(nsizech+1)=ibset(charch(nsizech+1),7)
                        do nbit=0,4 !number of q&Q
                           if(btest(nqq,nbit+8)) then
                              charch(nsizech+1)
     &                             = ibset(charch(nsizech+1),nbit)
                           end if
                        end do
                        do nbit=0,7 !number of q&Q
                           if(btest(nqq,nbit)) then
                              charch(nsizech+2)
     &                             = ibset(charch(nsizech+2),nbit)
                           end if
                        end do
                        if(mod(nqq,2).ne.0) then
                           nqq=nqq+1
                           charch(nsizech+2+nqq)=255 !P
                        end if
                        nsizech = nsizech + 2 + nqq ! # of q&Q + qQ
                        charch(nsizech+1)=0 !D
                        charch(nsizech+2)=0 !D
                        nsizech=nsizech+2 !D
                        hit_flag=0 ! end of hit

C     820 debug start
                        if(nsizech.gt.NFADC_RANGE*2) then
                           nsizech=nsizech-4-nqq-2 ! Hit pos -2
                           goto 33
                        end if
C     820 debug end
                        
C     820 debug start
                     else if(j.eq.NFADC_RANGE) then
                        nsizech = nsizech - 2
                        
                        if(charch(nsizech-1).eq.0 .and.
     $                       charch(nsizech).eq.0  ) then
                           charch(nsizech+1)=0
                           charch(nsizech+2)=0
                           charch(nsizech+3)=0
                           charch(nsizech+4)=0
                           charch(nsizech+5)=0
                           charch(nsizech+6)=0
                           charch(nsizech+7)=0
                           goto 33
                        end if
                        
                        if(charch(nsizech-1).eq.64 .and.
     $                       charch(nsizech).eq.0  ) then
                           do ii=1,18
                              charch(ii)=0
                           end do
                           nsizech=0
                           goto 55
                        end if
                        
C     820 debug end
                        
                     else
                        nqq=nqq+1
                        charch(nsizech+2+nqq)
     $                       = int(RFCHPLS(j,i)*count2mv)
                        hit_flag=1
                     endif
                  endif
               Endif ! end of pedestal check
            end do
C 818 debug
 33         Continue

            charch(2)=btest(nsizech,16)
            do nbit=0,7
               if(btest(nsizech,nbit+8)) then
                  charch(3) = ibset(charch(3),nbit)
               end if
               if(btest(nsizech,nbit)) then
                  charch(4) = ibset(charch(4),nbit)
               end if
            end do
            ibuf2(i+1)=ninttot+1
            do ii=1,int((nsizech+3)/4)
               do nbit=0,7
                  if(btest(charch((ii-1)*4+1),nbit)) then
                     ibuf3(ninttot+ii)= ibset(ibuf3(ninttot+ii),nbit+24)
                  end if
                  if(btest(charch((ii-1)*4+2),nbit)) then
                     ibuf3(ninttot+ii)= ibset(ibuf3(ninttot+ii),nbit+16)
                  end if
                  if(btest(charch((ii-1)*4+3),nbit)) then
                     ibuf3(ninttot+ii)= ibset(ibuf3(ninttot+ii),nbit+8)
                  end if
                  if(btest(charch((ii-1)*4+4),nbit)) then
                     ibuf3(ninttot+ii)= ibset(ibuf3(ninttot+ii),nbit)
                  end if
               end do
            end do

 55         Continue

            ninttot=ninttot+int((nsizech+3)/4)
            nsizetot=nsizetot+int((nsizech+3)/4)*4
            Call VZERO(charch,NFADC_RANGE*2) 
         end if
      ENDDO

C*********************************************************************
C     Compress Data & Fill Bank
C*********************************************************************
*-- Creat Bank 'FADC' 
      Call KZBCR0('FADC',IERR)
*-- Header
      ibufh(1) = 6
      ibufh(2) = 1
      ibufh(3) = 3
      ibufh(4) = 961
      ibufh(5) = 961
      ibufh(6) = ninttot
      Call kzrep0('FADC',0,'I',nseg0,ibufh)
      call kzrep0('FADC',1,'I',nseg1,ibuf1)
      call kzrep0('FADC',2,'I',nseg2,ibuf2)
      call kzrep0('FADC',3,'I',ninttot,ibuf3)

      RETURN
      END
