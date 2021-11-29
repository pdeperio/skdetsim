      SUBROUTINE SKIDENTIFY(LUN)
c ======================================================================
c This prints a set of useful information to the specified LUN.
c Usually, this is the log file.
c
c Original author:    Ed Kearns   September 8, 2007
c ======================================================================
      IMPLICIT NONE
#include "skdetsim_svn_version.h"
      CHARACTER*128 TMPSTRING
      INTEGER I, LUN, LENSTR, IDU, NPAR, IARGC
      INTEGER LENOCC
      EXTERNAL LENOCC
      CHARACTER*20 BANNER
      DATA BANNER/'####################'/
c
c     Write banner
c 
      WRITE(LUN,'(A)') ' '
      WRITE(LUN,'(A,A,A,A)') BANNER,BANNER,BANNER,BANNER 
      WRITE(LUN,'(A,A,A,A)') BANNER,' SKDETSIM environment:'
      WRITE(LUN,'(A,A,A,A)') BANNER,BANNER,BANNER,BANNER 
c
c     Username (if set as environment variable)
c
      CALL GETENV('LOGNAME',TMPSTRING)
      IF (LENOCC(TMPSTRING) .GT. 0) THEN
         WRITE(LUN,'(A,A)') 'LOGNAME: ',
     &        TMPSTRING(1:LENOCC(TMPSTRING))
      END IF
c
c     User ID (must be set in Unix)
c
      CALL GETUIDF(IDU)
      WRITE(LUN,'(A,I6)') 'User ID: ',IDU
c
c     Executable program and arguments
c
      CALL GETARG(0,TMPSTRING)
      IF (LENOCC(TMPSTRING) .GT. 0) THEN
         WRITE(LUN,'(A,A)') 'Executable program: ',
     &        TMPSTRING(1:LENOCC(TMPSTRING))
      END IF
      NPAR = IARGC()
      DO I = 1, NPAR
         CALL GETARG(I,TMPSTRING)
         WRITE(LUN,'(A,I3,A,A)') 'Argument ',I,': ',
     &        TMPSTRING(1:LENOCC(TMPSTRING))
      END DO
c
c     Working directory
c
      CALL GETWDF(TMPSTRING)
      WRITE(LUN,'(A,A)') 'Working directory: ',
     &     TMPSTRING(1:LENOCC(TMPSTRING))
c
c     CERNLIB
c
      CALL GETENV('CERN_ROOT',TMPSTRING)
      IF (LENOCC(TMPSTRING) .GT. 0) THEN
         WRITE(LUN,'(A,A)') 'CERN_ROOT: ',
     &        TMPSTRING(1:LENOCC(TMPSTRING))
      END IF
c
c     Super-K specific environments
c
      CALL GETENV('NEUT_ROOT',TMPSTRING)
      IF (LENOCC(TMPSTRING) .GT. 0) THEN
         WRITE(LUN,'(A,A)') 'NEUT_ROOT: ',
     &        TMPSTRING(1:LENOCC(TMPSTRING))
      END IF
      CALL GETENV('EXPERIMENT_ROOT',TMPSTRING)
      IF (LENOCC(TMPSTRING) .GT. 0) THEN
         WRITE(LUN,'(A,A)') 'EXPERIMENT_ROOT: ',
     &        TMPSTRING(1:LENOCC(TMPSTRING))
      END IF
      CALL GETENV('PACKAGE_ROOT',TMPSTRING)
      IF (LENOCC(TMPSTRING) .GT. 0) THEN
         WRITE(LUN,'(A,A)') 'PACKAGE_ROOT: ',
     &        TMPSTRING(1:LENOCC(TMPSTRING))
      END IF
c
c     Write version numbers:
c
      WRITE(LUN,'(A,F10.3)') 'APDETSIM_VERSION',APDETSIM_VERSION
      WRITE(LUN,'(A,F10.3)') 'SKDETSIM_VERSION',SKDETSIM_VERSION
c
c     Write CVS keywords
c     NOTE: Since we have switched to SVN this no longer is useful
c     ETK commented it out
c
cc      WRITE(LUN,'(A)') "Automatically expanded CVS keywords: "
cc      WRITE(LUN,'(A)') "$Date: 2007-11-27 10:04:28 $"
cc      WRITE(LUN,'(A)') 
cc     &     "$Source: /home/cvs/skrep/repository/skdetsim/skidentify.F,v $"
cc      WRITE(LUN,'(A)') "$Revision: 1.6 $"
cc      WRITE(LUN,'(A)') "$Name: not supported by cvs2svn $"
cc      WRITE(LUN,'(A)') "$State: Exp $"

c
c     Write svn version number (added 3/2013, ETK)
c
      WRITE(LUN,'(A,A)') 'SVN VERSION NUMBER = ',SKDETSIM_SVN_VERSION
c
c     Write banner
c 
      WRITE(LUN,'(A,A,A,A)') BANNER,BANNER,BANNER,BANNER 
      WRITE(LUN,'(A)') ' '

      RETURN
      END
