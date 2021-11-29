***********************************************************************
*     -----------------
      SUBROUTINE SGACREFR(dir,normal,normalconv,costh,theta2,newvect)
*     -----------------
*
*     (Purpose)
*       Determine direction of refracted photon at water/acrylic interface
*
*     (Input)
*     dir(3)::      direction of incident photon in MARS
*     normal(3)::   direction of normal to interface
*     normalconv::  convention denoting which side of the interface the 
*                            normal points to
*                -1= angle between incoming photon and normal is >90
*                             normal points into new material
*
*                   incoming \       
*                             \   old material
*                              \  
*                        --------------------- interface
*                               |
*                               |   new material
*                        normal |
*                               V
*
*                +1= angle between incoming photon and normal is <90
*                             normal points into old material
*
*                               ^
*                  incoming \   |     old material 
*                            \  |
*                             \ | normal       
*                              \|  
*                        --------------------- interface
*                         new material            
*              
*     costh::       cosine of angle between dir(incident) and 
*                       normal to interface (>0;<90deg)
*     theta2::      angle of refracted photon
*
*     (Output)
*     newvect::     new direction of refracted photon
*
*     (Creation Date and Author)
*       2002.07.16 ; First version by Chris Sterner
*
************************************************************************

      IMPLICIT NONE 
      
C     INPUT
      REAL dir, normal, costh, theta2
      INTEGER normalconv
      DIMENSION dir(3), normal(3)
C     OUTPUT
      REAL newvect
      DIMENSION newvect(3)
C     LOCAL
      REAL length, vvect(3)
      INTEGER I
      

C-----Determine new direction of photon, newvect
C   create vector, vvect, orthogonal to normal and co-planar to normal&dir
C     vvect is tangent to surface. 
C     vvect and normal make 2-dim coordiante system 
       DO 40 I=1,3
          vvect(I)=dir(I)-1.0*FLOAT(normalconv)*costh*normal(I)
 40    continue
C       normalize vvect
       length=sqrt(vvect(1)**2+vvect(2)**2+vvect(3)**2)
       if(length.ne.0.AND.length.ne.1.0) then
          DO 50 I=1,3
             vvect(I)=vvect(I)/length
 50       continue
       endif
      
C      create newvect out of vector sum of vvect and normal
       DO 60 I=1,3
          newvect(I)=FLOAT(normalconv)*COS(theta2)*normal(I)+
     &         SIN(theta2)*vvect(I)
 60    continue
C      normalize newvect
       length=sqrt(newvect(1)**2+newvect(2)**2+newvect(3)**2)
       if(length.ne.0.AND.length.ne.1.0) then
          DO 70 I=1,3
             newvect(I)=newvect(I)/length
 70       continue
       endif

c       WRITE(*,*) 'IN REFR',(vvect(I),ii=1,3)


      RETURN 
      END
