      SUBROUTINE SGGTXT(IEND)
*********************************************
*  Create MC particle banks based on input from
*  a USMC/Nuance style (text) kinematic file
*
*  IEND returns 0 if an event was read, -1 if EOF
*
*  D.Casper  28-7-98
*
*  Modifications:
*     20210608 - Guillaume Pronost: add multi-vertex
*     support. See examples/multivertex_input.dat 
*     for example
*
**********************************************
      implicit none
#include "dsbnkf.h"
#include "nework.h"
#include "nesolact.h"
#include "vcwork.h"
#include "vcvrtx.h"
#include "mcevent.h"
      integer iend
      integer max_particles, max_vertex
      parameter (max_particles = 100)
      parameter (max_vertex = 100)
      integer oxygen_code, tau_code
      parameter (oxygen_code = 8016)
      parameter (tau_code = 15)
      real vertex(4,max_vertex), direction(3,max_particles),
     &     energy(max_particles), mass, momentum
      integer n_particles, particle(max_particles),
     &     class(max_particles), icnuance
      integer n_vertex, partvtx(max_particles)
      integer j
C
      character*132 line
      integer jr, jl
      real vmod
C
C+  Common block used by CHPACK (CERNLIB M432)
C
      integer nd,ne,nf,ng,num(2),dummy(34)
      common /slate/ nd,ne,nf,ng,num,dummy
      integer iclocl, icnext, icfnbl, lnblnk
C
      call kzeclr
      call vcclcm

      iend = -1
      line = ' '
      jl = 1
      jr = 132
      n_particles = 0
      n_vertex = 0
      do while (iclocl('end',3,line,jl,jr).eq.0 .and.
     &          iclocl('stop',4,line,jl,jr).eq.0)
         read(lin,'(A)',end=100,err=100) line
         jl = icfnbl(line,1,132)
         if (ng.ne.0) then   !Not a blank line
            jr = lnblnk(line)
            if (iclocl('vertex',6,line,jl,jr).ne.0) then
               n_vertex = n_vertex + 1
               jl = iclocl('vertex',6,line,jl,jr)
               jl = icnext(line,jl,jr)   ! Find next blank after 'vertex'
               jl = ne
               read(line(jl:jr),*) (vertex(j,n_vertex),j=1,4)
            else if (iclocl('track',5,line,jl,jr).ne.0) then
               if (n_particles .eq. max_particles) then
                  write(6,*) 'SGGTXT - Maximum number of particles exceeded'
                  stop
               endif
               n_particles = n_particles + 1
               jl = iclocl('track',5,line,jl,jr)
               jl = icnext(line,jl,jr)
               jl = ne
               read(line(jl:jr),*)particle(n_particles),
     &              energy(n_particles),
     &              (direction(j,n_particles),j=1,3),
     &              class(n_particles)
               call vunit(direction(1,n_particles),
     &                    direction(1,n_particles),3)
               partvtx(n_particles) = n_vertex ! 20210608 - Guillaume: set vertex id for the track
            else if (iclocl('nuance',6,line,jl,jr).ne.0) then !information
               jl = iclocl('nuance',6,line,jl,jr)
               jl = icnext(line,jl,jr)
               jl = ne
               read(line(jl:jr),*) icnuance
            else if (iclocl('event',5,line,jl,jr).ne.0) then !event number
               jl = iclocl('event',5,line,jl,jr)
               jl = icnext(line,jl,jr)
               jl = ne
               read(line(jl:jr),*) iev
               iev=iev-1
            endif
         endif
      enddo
      if (iclocl('stop',4,line,jl,jr).ne.0) return
C
C+  Make HEADER and other canonical banks
C
      call mcmkhd
      call mcmkmh
      call mcmkwt
C
C+  Fill the NEUT bank with particles the primary particles
C
      numne = 0
C      call nuaneu(icnuance,icneut)
C      modene = ichannel   !No equivalence between modes with NEUT
C      modene = icneut
      modene = icnuance
      isol = -1
      solact  = -1.
      do j = 1, n_particles
         if (class(j).lt.0) then
            numne = numne + 1
            ipne(numne) = particle(j)
            call mcmass(particle(j),mass)
            if (energy(j).gt.mass) then
               momentum = sqrt(energy(j)**2-mass**2)/1000.
               call vscale(direction(1,j),momentum,pne(1,numne),3)
            else
               if ((ipne(numne).eq.2112 .or. ipne(numne).eq.2212).and.
     &              (class(j).eq.-1).and.(vmod(direction(1,j),3).gt.0.5)) then
                  energy(j) = energy(j) + 35.
                  if (energy(j).ge.mass) then
                     momentum = sqrt(energy(j)**2-mass**2)/1000.
                     call vscale(direction(1,j),momentum,pne(1,numne),3)
                  else
                     call vzero(pne(1,numne),3)
                  endif
               else
                  call vzero(pne(1,numne),3)
               endif
            endif
            iorgne(numne) = 0
            if (class(j).eq.-1) then   ! Initial state particle
               iflgne(numne) = -1
               icrnne(numne) = 0
            else                       ! Primary outgoing particle
               iflgne(numne) = 0
               icrnne(numne) = 1
            endif
         endif
      enddo
      call nemknebk(vertex)
C
C+  Fill the MCVERTEX bank with the primary
C
c      print*,"found vtx: ",n_vertex
      do j = 1, n_vertex
         nvtxvc = j
c         call ucopy(vertex,pvtxvc(1,j),3)
c ---    20210608 - Guillaume: Set initial vertex to previous vertex in the fileFill vertex array:         
         pvtxvc(1,j) = vertex(1,j)
         pvtxvc(2,j) = vertex(2,j)
         pvtxvc(3,j) = vertex(3,j)
         iflvvc(j)   = 0
         iparvc(j)   = 0
         timvvc(j)   = vertex(4,j) ! ! 20210608 - Guillaume: Set time
      enddo
      call vcmkvx
         
C
C+  Fill the MCVECT bank with particles leaving the nucleus
C+  Tau is special case, since we want it in the outgoing particles
C+  even though it is class -2
C
      nvc = 0
      do j = 1, n_particles
         if (class(j).ge.-1 .or. 
     &       ((class(j).eq.-2).and.(abs(particle(j)).eq.tau_code))) then
            nvc = nvc + 1
            ipvc(nvc) = particle(j)
            call mcmass(particle(j),amasvc(nvc))
c            print*,"Add particle",particle(j),energy(j),amasvc(nvc)
            if (energy(j).gt.amasvc(nvc)) then
               momentum = sqrt(energy(j)**2-amasvc(nvc)**2)
               call vscale(direction(1,j),momentum,pvc(1,nvc),3)
            else
               call vzero(pvc(1,nvc),3)
            endif
            iorgvc(nvc) = 0
            ivtivc(nvc) = partvtx(j) ! 20210608 - Guillaume: Set initial vertex to previous vertex in the file
            if (class(j).lt.0 .or. ipvc(nvc).eq.oxygen_code) then
               icrnvc(nvc) = 0
               iflgvc(nvc) = -1
            else
               icrnvc(nvc) = 1
               iflgvc(nvc) = 0
            endif
            ivtfvc(nvc) = partvtx(j) ! 20210608 - Guillaume: Set final vertex to previous vertex in the file
         endif
      enddo
      call vcmkvc
C+
      iend = 0
C+
 100  continue
      return
      end
