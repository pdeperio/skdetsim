***********************************************************************
*     -----------------
      subroutine sgacrefl(dir,normal,normalconv,costh,newvect)
*     -----------------
*
*     (Purpose)
*        Determine direction of cherenkov photon reflected from
*     water/acrylic interface.
*     
*     (Input)
*       dir::       direction of incident photon in MasterReferenceSystem
*       normal::    direction of normal to the interface 
*       normalconv::  convention denoting which side of the interface the 
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
*       costh::     cosine of angle between dir and normal (>0;<90deg)
*
*     (Output)
*       newvect::   new direction of reflected photon
*
*     (Creation Date and Author)
*       2002.07.11 ; First version by Chris Sterner
*
************************************************************************


      implicit none
C     INPUT
      REAL dir, normal, costh
      INTEGER normalconv
C     OUTPUT
      REAL newvect
      DIMENSION dir(3), normal(3), newvect(3)


                                !direction of simple reflection
      newvect(1)=-1.0*FLOAT(normalconv)*2*costh*normal(1)+dir(1)
      newvect(2)=-1.0*FLOAT(normalconv)*2*costh*normal(2)+dir(2)
      newvect(3)=-1.0*FLOAT(normalconv)*2*costh*normal(3)+dir(3)

   
*     Random reflection direction...not yet

      RETURN
      END
