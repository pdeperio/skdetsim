************************************************************************
*     -------------------------------------
      SUBROUTINE  dsreso_qb(icab,pc,rpc)
*     -------------------------------------
*
*     (Purpose)
*       estimate resolution effect of QBee     
*       (used in dsthr_sk4)
*     (INPUT) 
*      icab   : cable number
*      pc     : charge value (pC unit)
*            
*      (OUTPUT)
*      rpc    : charge value after considering resolution(pC unit)
*        
*            
*     (Creation Date and Author)
*       2009.01.18 ; First version by T.Tanaka
*
************************************************************************

      implicit none
      real pc,rpc,tem,tem2,cont,pctem,onethr
 
      integer n,icab,lun,lun1,lun2,lun3,lun4,all
      parameter (lun=99)
      parameter (lun1=98)
      parameter (lun2=97)
      parameter (lun3=96) 
      parameter (lun4=95)

      integer bin1, bin2
      parameter (bin1=180)
      parameter (bin2=500)
      real xth2(bin1),yth2(bin1),xre2(bin2),yre2(bin2),tempy
      real xth3(bin1),yth3(bin1),xre3(bin2),yre3(bin2)
      real xth(bin1),yth(bin1),xre(bin2),yre(bin2)
      integer ii,dbin,rbin
  
      integer maxpm
      parameter (maxpm=11146)
      integer n0, n1, n2, isk23(maxpm)
      

      character*132 file_name
      data file_name/'/skam/const/pmtinf.dat'/
  
      character*132 file_name1,file_name2
      character*132 file_name3,file_name4

      data file_name1/'sk3one_sk2pmt.dat'/
      data file_name2/'sk4one_sk2pmt.dat'/
      data file_name3/'sk3one_sk3pmt.dat'/
      data file_name4/'sk4one_sk3pmt.dat'/


      integer iresave/0/
      save iresave,xth2,yth2,xre2,yre2
      save xth3,yth3,xre3,yre3
  
C       write(6,*) 'ds',pc,rpc
 
CCC       open dat file for resolution 
           
       if(iresave.eq.0) then
          iresave = 1
          call findconsts(file_name,'pmtinf','pmtinf.dat',0)
          open(LUN,file=file_name,status='old')
 3        continue
          read(LUN,*,END=4) n0, n1, n2
           isk23(n0) = n1
          goto 3
 4        continue
          close(LUN)

           all=0
           call findconsts(file_name1,'sk3one_sk2pmt.dat','sk3one_sk2pmt.dat',0)
           open(lun1,file=file_name1,status='old')
 110       read(lun1,*,end=220) pctem,tem,tem2,cont 
            all=all+1
            xth2(all)=pctem
            yth2(all)=cont
            goto 110
 220       close(lun1)

           all =0
           call findconsts(file_name2,'sk4one_sk2pmt.dat','sk4one_sk2pmt.dat',0)
           open(lun2,file=file_name2,status='old')
 111       read(lun2,*,end=221) pctem,tem,cont
            all=all+1
            xre2(all)=pctem
            yre2(all)=cont

            goto 111
 221       close(lun2)

           all=0
           call findconsts(file_name3,'sk3one_sk3pmt.dat','sk3one_sk3pmt.dat',0)
           open(lun3,file=file_name3,status='old')
 112       read(lun3,*,end=222) pctem,tem,tem2,cont
            all=all+1 
            xth3(all)=pctem 
            yth3(all)=cont

            goto 112
 222       close(lun3)
 
           all =0

           call findconsts(file_name4,'sk4one_sk3pmt.dat','sk4one_sk3pmt.dat',0)
           open(lun4,file=file_name4,status='old')
 113       read(lun4,*,end=223) pctem,tem,cont
            all=all+1
            xre3(all)=pctem
            yre3(all)=cont

            goto 113
 223       close(lun4)

           endif

         if(isk23(icab).eq.2) then
         do n=1,bin1
         xth(n)=xth2(n)
         yth(n)=yth2(n)
         enddo 
         do n=1,bin2
         xre(n)=xre2(n)
         yre(n)=yre2(n)
         enddo
         onethr=0.80
         else  !! sk3pmt and new pmt

         do n=1,bin1
         xth(n)=xth3(n)
         yth(n)=yth3(n)
         enddo
         do n=1,bin2
         xre(n)=xre3(n)
         yre(n)=yre3(n)
         enddo
         onethr=0.83   

         endif
   
         
         

CCC     apply resolution effect
      if (pc.ge.onethr) then
       rpc = pc 
      else
C  search which bin exceed 'pc'
       do ii=bin1,1, -1
C       write(6,*) pc,xth(ii)
         if (pc. ge. xth(ii)) then
          dbin=ii
          goto 2000
         endif
       enddo
 2000    continue
C       write(6,*) dbin,xth(dbin)

C using interpolating function to serach more near horizontal 
C value for integraled value
        tempy=(yth(dbin+1)-yth(dbin))/(xth(dbin+1)-xth(dbin))
     &            *(pc-xth(dbin)) +yth(dbin) 

        do ii=bin2,1, -1
         if (tempy. ge. yre(ii)) then
         rbin=ii
         goto 3000
         endif
         enddo

 3000    continue
        rpc =((xre(rbin+1)-xre(rbin))/(yre(rbin+1)-yre(rbin)))
     &        *(tempy-yre(rbin)) +xre(rbin)      
      endif

      return
      end


       



