************************************************************************
*     -----------------
      SUBROUTINE SGCARD
*     -----------------
*
*     (Purpose)
*       Read VECT card and setup common blocks.
*       Called by : SGINIT
*
*     (Input)
*       card file
*
*     (Output)
*
*     (Creation Date and Author)
*       1995.1.30 ; Koshio
*       1999.5.23 ; D.Casper   Allow user-initialization of random seeds
*       2003.7. 4 ; C. Sterner skdetsim random seed initialization
************************************************************************

      IMPLICIT NONE

#include "dsbnkf.h"
#include "geant321/gcflag.inc"
#include "guflag.h"

      INTEGER RNDCARD(5)
      REAL RVTX(3),RDIR(3),RMOM(2),RWALL
      INTEGER IMVTX,IMDIR,IPART,IMMOM
      INTEGER IRET, ISEQ, iseed_rdm1, iseed_rdm2, iseed_rdm3, iseed_rdm4, iseed_rdm5, LUX
      COMMON/VECTCARD/IMVTX,RVTX,RWALL,IMDIR,RDIR,
     &     IPART,IMMOM,RMOM,RNDCARD
      CHARACTER*40 TMPFIL, COMMAND

C -- initialize parameter
      DATA RVTX,RDIR,RMOM /3*0.0,-1.,2*0.0,10.0,10000./
      RWALL = 0.

      IFLGV = -1
      IFANTI = 0
      NEVENT = 0
      IMVTX = 0
      IMDIR = 0
      IPART = 0
      IMMOM = 0

C -- extract VECT card

      CALL GETUNI(TMPFIL)
      OPEN(UNIT=90, FILE=FCARD, STATUS='OLD')
      OPEN(UNIT=91, FILE=TMPFIL,STATUS='UNKNOWN')      
      CALL EXCARD('VECT', 90, 91, IRET)
      CLOSE(UNIT=91)
      CLOSE(UNIT=90)

      OPEN(UNIT=92, FILE=TMPFIL,STATUS='OLD')
#ifdef absoft
C     Ugly work-around of problem passing LUN from Absoft to f2c or g77.
      call system('/bin/ln -sf '//TMPFIL//' fort.92')
#endif
      CALL FFSET('LINP', 92)
      CALL FFKEY('FILE',IFLGV, 1,'INTEGER')
      CALL FFKEY('NEVT',NEVENT,1,'INTEGER')
      CALL FFKEY('RAND',RNDCARD,5,'INTEGER')
      CALL FFKEY('MVTX',IMVTX, 1,'INTEGER')
      CALL FFKEY('VTX', RVTX , 3,'REAL')
      CALL FFKEY('WALL',RWALL, 1,'REAL')
      CALL FFKEY('MDIR',IMDIR, 1,'INTEGER')
      CALL FFKEY('DIR', RDIR , 3,'REAL')
      CALL FFKEY('PART',IPART, 1,'INTEGER')
      CALL FFKEY('MMOM',IMMOM, 1,'INTEGER')
      CALL FFKEY('MOM' ,RMOM , 2,'REAL')
      CALL FFGO
      CLOSE(UNIT=92)

C
C  Initialize random number seeds
C
CCS Comment on random number seeds:  skdetsim has been modified to use
C     random number seeds to restart at a specific event as of July 2003.
C     Two random number generators are used GRNDM (a GEANT random number
C     generator) and RANLUX (another CERNLIB generator).  Caution must be taken
C     in chosing random number seeds initially.  The FFREAD routine above 
C     which reads the supersim.card file only reads integers with less than 
C     about ten digits.  Trying to use more than ten digits in a seed give 
C     unpredictable results.
C     Also, the RANLUX generator can only be set with specific seeds. Not all
C     integer combinations work for its three seeds.  An error is output:
C Error in RESTARTING with RLUXGO:
C The values      97822      17466          0 cannot occur at luxury level    3
C     This error is output when seeds are initialized.  After this error all 
C     subsequent RANLUX seeds reported are useless and meaningless.  Please
C     make sure that RANLUX does not do this.  Setting the last three seeds
C     to 0 0 0 will give you a default RANLUX initialization.  
C     Seeds have been tested to be workable but problematic.
C     C. Sterner July 2003

      write(*,*) '::::::::::::::::::::::::::::::::::::::::::::::::::::::::'
      write(*,*) 'Initalizing random number seeds:'
      write(*,*) rndcard(1),rndcard(2),rndcard(3),rndcard(4),rndcard(5)

      IF(rndcard(2).NE.0)THEN
         WRITE(*,*) 'grndmq call type 1'
         CALL GRNDMQ(rndcard(1),rndcard(2),1,'S')

      ELSEIF(rndcard(1).GT.0)THEN
         ISEQ=rndcard(1)
         WRITE(*,*) 'grndmq call type 2'
         CALL GRNDMQ(rndcard(1),rndcard(2),ISEQ,'Q')
         CALL GRNDMQ(rndcard(1),rndcard(2),ISEQ,'S')
      ENDIF
      call grndmq(iseed_rdm1,iseed_rdm2,0,'G')

      if(rndcard(4)+rndcard(3)+rndcard(5).gt.0.0) then
         write(*,*) 'RANLUXGO call'
         CALL RLUXGO(3,rndcard(3),rndcard(4),rndcard(5))
         CALL RLUXAT(LUX,ISEED_RDM3,ISEED_RDM4,ISEED_RDM5)
      endif
      
      write(6,*) 'sgcard: random seeds initialized: '
     &     ,iseed_rdm1,iseed_rdm2,iseed_rdm3
     &     ,iseed_rdm4,iseed_rdm5
      write(*,*) '::::::::::::::::::::::::::::::::::::::::::::::::::::::::'

      WRITE(6,*) ''
      IF(IFLGV.EQ.2) THEN
        WRITE(6,*) '* read VECTOR from text file'
      ELSEIF(IFLGV.EQ.1) THEN
        WRITE(6,*) '* read VECTOR from card'
      ELSEIF(IFLGV.EQ.0) THEN
        WRITE(6,*) '* read VECTOR from ZBS file'
      ELSEIF(IFLGV.EQ.-1) THEN
        WRITE(6,*) '* read VECTOR from SGGVUS'
      ELSEIF(IFLGV.EQ.-2) THEN
        WRITE(6,*) '* read CONVVECT from file'
      ELSE
        WRITE(0,*) 'SGCARD: Error : VECT-FILE =', IFLGV
        STOP 1
      ENDIF

      IF(NEVENT.GT.0) THEN
        WRITE(6,*) '* the number of generated event was set to ', NEVENT
      ELSE
        WRITE(0,*) 'SGCARD: Error : VECT-NEVT ', NEVENT
        STOP 1
      ENDIF

      IF(IFLGV.EQ.1) THEN

        IF(IMVTX.EQ.1) THEN
          WRITE(6,*) '* Vetrex positon(X,Y,Z)(cm) : ', RVTX
        ELSEIF(IMVTX.EQ.2) THEN
          WRITE(6,*) '* Random vertex inside ', RWALL,'(cm) from PMT plane'
        ELSE
          WRITE(0,*) 'SGCARD: Error : VECT-VTX =', IMVTX
          STOP 1
        ENDIF

        IF(IMDIR.EQ.1) THEN
          WRITE(6,*) '* Direction Cosines : ',RDIR
        ELSEIF(IMDIR.EQ.2) THEN
          WRITE(6,*) '* Random direction'
        ELSE
          WRITE(0,*) 'SGCARD: Error : VECT-DIR =', IMDIR
          STOP 1
        ENDIF

        IF(IPART.NE.0) THEN
          WRITE(6,*) '* Particle ID : ',IPART
        ELSE
          WRITE(0,*) 'SGCARD: Error : VECT-PART =', IPART
          STOP 1
        ENDIF

        IF(IMMOM.EQ.1) THEN
          IF(RMOM(1).GE.0.) THEN
            WRITE(6,*) '* Fix momentum(MeV) : ',RMOM(1)
          ELSE
            WRITE(0,*) 'SGCARD: Error in momentum(MeV) : VECT-MOM ', RMOM(1)
            STOP 1
          ENDIF
        ELSEIF(IMMOM.EQ.2) THEN
          IF(RMOM(1).GE.0. .AND. RMOM(2).GE.0.) THEN
            WRITE(6,*) '* Randomly set momentum(MeV) between ', RMOM(1), ' and ', RMOM(2)
          ELSE
            WRITE(0,*) 'SGCARD: Error in momentum(MeV) : VECT-MOM ', RMOM
            STOP 1
          ENDIF
        ELSE
          WRITE(0,*) 'SGCARD: Error : VECT-MMOM =', IMMOM
          STOP 1
        ENDIF
      ENDIF
      WRITE(6,*) ''

#ifdef absoft
C     Ugly work-around of problem passing LUN from Absoft to f2c or g77.
      call system('/bin/rm fort.92')
#endif
      COMMAND = '/bin/rm '//TMPFIL
      CALL SYSTEM(COMMAND)

      RETURN
      END


