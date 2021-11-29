      subroutine lingeom1

      implicit none

#include "linpara.h"
      
      integer iii, irotate, ivolu

      real parflm(6),parfra(6),ASC(2),ZSC(2),WSC(2)
      real TMAXMD, STEMAX, DEEMAX, STMIN, FLDWTI, EPSIL

      real DSC, RLSC, tmidfra, scinz, scinr, guidex, guidey, guidez
      DATA DSC/1.032/
      DATA RLSC/42.4/

      write(*,*) ' ******************************************'
      write(*,*) ' *** LINAC geometry routine is called. ****'
      write(*,*) ' ******************************************'
      call sleep(2)

C---------> Scintilator (Plastic)
C
      ASC(1)= 12.01
      ASC(2)= 1.01
      ZSC(1)= 6.
      ZSC(2)= 1.
      WSC(1)= 0.9224
      WSC(2)= 0.0776

      CALL GSMATE(50,'TITAN$',47.88,22.,4.54,3.56,1.,0,0)
      CALL GSMATE(51,'IRON$',55.850,26.,7.87,1.76,1.,0,0)
      CALL GSMIXT(53,'SC$',ASC,ZSC,DSC,2,WSC)
c
      TMAXMD = 0.
      STEMAX = 0.001 
      DEEMAX = 0.01
      STMIN = 0.000001
      FLDWTI = 0.
      EPSIL = 0.0001

      DO 200 III=50,53
         CALL GSTMED( III, 'IRON$', 51, 0, 0,    
     +        FLDWTI,tmaxmd,STEMAX,DEEMAX,EPSIL,STMIN,0,0)
 200  CONTINUE

      CALL GSTMED( 60, 'titanium', 50, 0, 0,
     + FLDWTI,tmaxmd,STEMAX,DEEMAX,EPSIL,STMIN,0,0)

      CALL GSTMED( 70, 'scinti', 53, 0, 0,
     + FLDWTI,tmaxmd,STEMAX,DEEMAX,EPSIL,STMIN,0,0)

      CALL GSTMED( 80, 'lgde', 53, 0, 0,
     + FLDWTI,tmaxmd,STEMAX,DEEMAX,EPSIL,STMIN,0,0)
      CALL GSTMED( 81, 'lgde', 53, 0, 0,
     + FLDWTI,tmaxmd,STEMAX,DEEMAX,EPSIL,STMIN,0,0)
      CALL GSTMED( 82, 'lgde', 53, 0, 0,
     + FLDWTI,tmaxmd,STEMAX,DEEMAX,EPSIL,STMIN,0,0)

      STEMAX = 1.0E-4
      CALL GSTMED( 71, 'almi', 9, 0, 0,
     + FLDWTI,tmaxmd,STEMAX,DEEMAX,EPSIL,STMIN,0,0)

      CALL GSTMED( 72, 'almi', 9, 0, 0,
     + FLDWTI,tmaxmd,STEMAX,DEEMAX,EPSIL,STMIN,0,0)

      TMAXMD = 0.
      STEMAX = 1.0 
      DEEMAX = 0.01
      STMIN = 0.1
      FLDWTI = 0.
      EPSIL = 0.001
C--
      CALL GSTMED( 73, 'vacuum', 16, 0, 0,
     + FLDWTI,tmaxmd,STEMAX,DEEMAX,EPSIL,STMIN,0,0)
C--
      CALL GSTMED( 74, 'vacuum', 16, 0, 0,
     + FLDWTI,tmaxmd,STEMAX,DEEMAX,EPSIL,STMIN,0,0)

C>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
C---------> print all tracking medium
C
      CALL GPTMED(0)
C<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
c
      rpip1i=2.99/2.
      rpip1o=4.27/2.
      pip1l =5.0/2.
      rpip2i=7.21/2.
      rpip2o=7.63/2.
      pip2l =14.2/2.
      rfra  =16./2.
      tfra  =2./2.
C--
      twin=0.01/2.
c      twtr=0.2
      twtr=0.36                 !98/12/14 Teisei(Isino-san's tune)

      IROTATE = 1
      tmidfra = 0.55/2.0
C--
c
c   vacuum
c
      PARFLM(1)=0.
      PARFLM(2)=rpip2i
      PARFLM(3)=pip2l
      CALL GSVOLU('vac1','TUBE',73,PARFLM,3,IVOLU)
      CALL GSPOS('vac1',1,'INVO',sorx,sory,sorz+2.*pip1l+pip2l
     +     ,IROTATE,'MANY')      
c
      PARFLM(1)=0.
      PARFLM(2)=rpip1i
      PARFLM(3)=pip1l-twtr/2.
      CALL GSVOLU('vac2','TUBE',74,PARFLM,3,IVOLU)
      CALL GSPOS('vac2',1,'INVO',sorx,sory,sorz+pip1l+twtr/2.
     +     ,IROTATE,'MANY')      
c
c  titanium window
c
      PARFLM(1)=0.
      PARFLM(2)=rpip1i
      PARFLM(3)=twin
      CALL GSVOLU('kapt','TUBE',60,PARFLM,3,IVOLU)
      CALL GSPOS('kapt',1,'INVO',sorx,sory,sorz+twin+twtr,
     +     IROTATE,'ONLY')
c
c   pipe, frange
c
      PARFLM(1)=rpip1i
      PARFLM(2)=rpip1o
      PARFLM(3)=pip1l
      CALL GSVOLU('pip1','TUBE',50,PARFLM,3,IVOLU)
      CALL GSPOS('pip1',1,'INVO',sorx,sory,sorz+pip1l
     +     ,IROTATE,'ONLY')
c      
      PARFLM(1)=rpip2i
      PARFLM(2)=rpip2o
      PARFLM(3)=pip2l
      CALL GSVOLU('pip2','TUBE',51,PARFLM,3,IVOLU)
      CALL GSPOS('pip2',1,'INVO',sorx,sory,sorz+pip1l*2.+pip2l
     +     ,IROTATE,'ONLY')

      PARFLM(1)=rpip1i
      PARFLM(2)=rpip2i
      PARFLM(3)=tmidfra
      CALL GSVOLU('fra1','TUBE',53,PARFLM,3,IVOLU)
      CALL GSPOS('fra1',1,'INVO',sorx,sory,sorz+pip1l*2+tmidfra
     +     ,IROTATE,'ONLY')

c
      PARFLM(1)=0.
      PARFLM(2)=rfra
      PARFLM(3)=tfra
      CALL GSVOLU('fran','TUBE',52,PARFLM,3,IVOLU)
      CALL GSPOS('fran',1,'INVO',sorx,sory,sorz+(pip1l+pip2l)*2.+tfra
     +     ,IROTATE,'ONLY')
c
c    scinti, almi foil
c
c
      scinz=1.74
C--
      scinr=1.2
      parfra(1)=0.
      parfra(2)=scinr
      parfra(3)=0.05
      CALL GSVOLU('scin','TUBE',70,parfra,3,IVOLU)
      CALL GSPOS('scin',1,'INVO',sorx,sory,sorz+scinz,IROTATE,'ONLY')
c
      parfra(1)=0.
      parfra(2)=1.2
C--
      parfra(3)=1.5e-3/2.
      call gsvolu('alm1','TUBE',71,parfra,3,ivolu)
      call gspos('alm1',1,'INVO',sorx,sory,sorz+scinz-0.055
     +     ,irotate,'ONLY')
      call gsvolu('alm2','TUBE',72,parfra,3,ivolu)
      call gspos('alm2',1,'INVO',sorx,sory,sorz+scinz+0.055
     +     ,irotate,'ONLY')

c
C light guides
C
      guidex = 1.6/2.0
      guidey = 0.2/2.0
      guidez = 6.0/2.0
      parfra(1)=guidex
      parfra(2)=guidey
      parfra(3)=guidez
      CALL GSVOLU('lgde','BOX ',80,parfra,3,IVOLU)
      CALL GSPOS('lgde',1,'INVO',sorx,sory+1.0,
     &            sorz+scinz+guidez+0.1,IROTATE,'ONLY')
c
      guidex = 1.6/2.0
      guidey = 0.2/2.0
      guidez = 0.5/2.0
      parfra(1)=guidex
      parfra(2)=guidey
      parfra(3)=guidez
      CALL GSVOLU('lgd1','BOX ',81,parfra,3,IVOLU)
      CALL GSPOS('lgd1',1,'INVO',sorx,sory+0.8,
     &            sorz+scinz+guidez+0.1,IROTATE,'ONLY')
c
      guidex = 1.6/2.0
      guidey = 0.2/2.0
      guidez = 0.3/2.0
      parfra(1)=guidex
      parfra(2)=guidey
      parfra(3)=guidez
      CALL GSVOLU('lgd2','BOX ',82,parfra,3,IVOLU)
      CALL GSPOS('lgd2',1,'INVO',sorx,sory+0.8,
     &            sorz+scinz-guidez-0.1,IROTATE,'ONLY')
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      return
      end


