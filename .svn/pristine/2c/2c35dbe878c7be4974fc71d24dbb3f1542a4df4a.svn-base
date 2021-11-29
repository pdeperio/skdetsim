************************************************************************
*     -----------------
      SUBROUTINE SGFFIN 
*     -----------------
*
*     (Purpose)
*       Define GEANT parameters
*
*     (Creation Date and Author)
*       1991.12.25 ; First version by J.Kanzaki
*       1993.12.07 ; Modified by M. Nakahata (extract card and set LUN)
      
*
************************************************************************
#include "geant321/gccuts.inc"
#include "dsbnkf.h"

      CHARACTER*40 TMPFIL, COMMAND

C -- EXTRACT GEANT CARD

      CALL GETUNI(TMPFIL)
      OPEN(UNIT=90, FILE=FCARD, STATUS='OLD')
      OPEN(UNIT=91, FILE=TMPFIL,STATUS='UNKNOWN')
      CALL EXCARD('GEANT', 90, 91, IRET)
      CLOSE(UNIT=91)
      CLOSE(UNIT=90)

      OPEN(UNIT=96, FILE=TMPFIL,STATUS='OLD')
#ifdef absoft
C     Ugly work-around of problem passing LUN from Absoft to f2c or g77.
      call system('/bin/ln -sf '//TMPFIL//' fort.96')
#endif
      CALL FFINIT(0)
      CALL FFSET('LINP', 96)
      CALL GFFGO
      CLOSE(UNIT=96)
#ifdef absoft
C     Ugly work-around of problem passing LUN from Absoft to f2c or g77.
      call system('/bin/rm fort.96')
#endif
      COMMAND = '/bin/rm '//TMPFIL
      CALL SYSTEM(COMMAND)

C-- to set the upper limit of the CUTNEU (cut off energy of neutral hadrons)
C-- cutneu should be smaller than 1MeV (kinetic energy)

#ifdef NEUTRON
c      if (CUTNEU.gt.1.e-14) CUTNEU = 1.e-14
      if (CUTNEU.gt.1.e-17) CUTNEU = 1.e-17
#else
      if (CUTNEU.gt.5.e-5) CUTNEU = 5.e-5
#endif

      RETURN
      END
