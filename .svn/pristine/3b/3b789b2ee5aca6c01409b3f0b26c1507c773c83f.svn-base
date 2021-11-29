      subroutine skdonuts_sk1(pmtpos,pmtdir,pos,dir,pos2,dir2,pos3,ihit)
*
**** Only for ledetsim in SK-I
*
c input
c   pmtpos: position of pmt
c   pmtdir: direction of pmt
c   pos   : position of photon
c   dir   : direction of photon
c output
c   pos2  : hit position or next position of photon when no hit
c   dir2  : new reflection direction when photon hits
c   pos3  : new photon position for geant pmt when reflected
c   ihit  : =0 no hit
c         : =1 hit
c
      implicit none
      integer ihit,idist,igeant
      real pmtpos(3),pmtdir(3),pos(3),dir(3),pos2(3),dir2(3),pos3(3)
      real v(3),d(3)
      real pi,dir1,d1,x1,y1,z1,r1,r2,hthr,hthrg
      parameter(pi=3.1415927,d1=12.7,r1=31.5,r2=25.4)
      parameter(x1=0.,y1=0.,z1=-d1)
      parameter(hthr=15.*d1/(r1-15.))
      parameter(hthrg=(r1**2-r2**2-d1**2)/2./d1)
      real x0,y0,z0,a,b,c,p,q,r,dist,x2,y2,z2,v1,v2,v3,vv1,dd
      real x3,y3,z3
      real dx,dy,dz,ad,pd2,a1,b1,c1,dd1,aa1
      real skddonuts,skdball
      real th,phi
c
c  move pmtpos --> (0,0,0)
      v1=pos(1)-pmtpos(1)
      v2=pos(2)-pmtpos(2)
      v3=pos(3)-pmtpos(3)
c
c  rotate pmtdir --> (0,0,1)
      th=acos(pmtdir(3))
      if(pmtdir(1).gt.0.)then
        if(pmtdir(2).gt.0.)then
          phi=atan(pmtdir(2)/pmtdir(1))
        else
          phi=atan(pmtdir(2)/pmtdir(1))+2.*pi
        endif
      elseif(pmtdir(1).lt.0.)then
        phi=atan(pmtdir(2)/pmtdir(1))+pi
      else
        if(pmtdir(2).gt.0.)then
          phi=pi/2.
        else
          phi=pi*1.5
        endif
      endif
c  rotate photon direction
c
      dir1=dir(1)*cos(phi)+dir(2)*sin(phi)
      d(2)=-dir(1)*sin(phi)+dir(2)*cos(phi)
c
      d(3)=dir(3)*cos(th)+dir1*sin(th)
      d(1)=-dir(3)*sin(th)+dir1*cos(th)
c  rotate photon position
      vv1=v1*cos(phi)+v2*sin(phi)
      v(2)=-v1*sin(phi)+v2*cos(phi)
c
      v(3)=v3*cos(th)+vv1*sin(th)
      v(1)=-v3*sin(th)+vv1*cos(th)
c
c  now photon vertex is 'v', and direction is 'd'
c  check hit
      x0=v(1)
      y0=v(2)
      z0=v(3)
c  photon direction
      a=d(1)
      b=d(2)
      c=d(3)
      ihit=0
      x2=x0
      y2=y0
      z2=z0
c  look for cross point on large sphere
      call skdpoint(x0,y0,z0,a,b,c,x1,y1,z1,r1,p,q,r,dist)
      if(dist.gt.r1)goto 200
      if(r.lt.hthrg)then
        call skdpoint(x0,y0,z0,a,b,c,x1,y1,0.,r2,p,q,r,dist)
        if(dist.gt.r2 .or. r.ge.hthrg)goto 200
      endif
c  check direction
      if((a*(p-x0)+b*(q-y0)+c*(r-z0)).lt.0. .and.
     &   sqrt((p-x0)**2+(q-y0)**2+(r-z0)**2).gt.1.)goto 200
c
      if(dist.le.r1)then
        if(r.ge.hthr)then
          ihit=1
          x2=p
          y2=q
          z2=r
          goto 200
        else
c  loop for check hit on donuts, 1mm step
          igeant=0
          do idist=1,508
            dd=float(idist)/10.
            x2=p+a*dd
            y2=q+b*dd
            z2=r+c*dd
            if(z2.lt.0.)goto 100
            if(skddonuts(x2,y2,z2).le.0.)then
              ihit=1
              goto 100
            endif
            if(z2.ge.hthrg)then
c  out of large sphere
              if(skdball(x1,y1,z1,r1,x2,y2,z2).gt.0.)goto 100
              igeant=1
            else
c  check geant pmt
              if(skdball(0.,0.,0.,r2,x2,y2,z2).le.0.)then
                igeant=1
              else
                if(igeant.eq.1)goto 100
              endif
            endif
          enddo
  100     continue
        endif
      endif
  200 continue
c  obtain reflection direction, position (x2,y2,z2)
c                     incident direction (a,b,c)
      if(ihit.eq.1)then
        if(z2.ge.hthr)then
c  normal vector, hit on large sphere
          dx=x2
          dy=y2
          dz=z2+d1
        else
c  normal vector, hit on donuts
          pd2=1.-10.4/sqrt(x2**2+y2**2)
          dx=x2*pd2
          dy=y2*pd2
          dz=z2
        endif
        ad=sqrt(dx**2+dy**2+dz**2)
        dx=dx/ad
        dy=dy/ad
        dz=dz/ad
        pd2=2.*(dx*a+dy*b+dz*c)
        a1=a-dx*pd2
        b1=b-dy*pd2
        c1=c-dz*pd2
        aa1=sqrt(a1**2+b1**2+c1**2)
        a1=a1/aa1
        b1=b1/aa1
        c1=c1/aa1
c  reflection direction is (a1,b1,c1)
c
c  obtain new vertex on geant surface
        if(z2.ge.hthr)then
          x3=x2
          y3=y2
          z3=z2
        else
          do idist=1,508
            dd=float(idist)/10.
            x3=x2+a1*dd
            y3=y2+b1*dd
            z3=z2+c1*dd
            if(z3.lt.0.)goto 300
            if(z3.ge.hthrg)then
              if(skdball(x1,y1,z1,r1,x3,y3,z3).gt.0.)goto 300
            else
              if(skdball(0.,0.,0.,r2,x3,y3,z3).gt.0.)goto 300
            endif
          enddo
        endif
  300   continue
      endif
        
c  inverse rotation
c   for position
      v(3)=z2*cos(th)-x2*sin(th)
      vv1=z2*sin(th)+x2*cos(th)
c
      v(1)=vv1*cos(phi)-y2*sin(phi)
      v(2)=vv1*sin(phi)+y2*cos(phi)
c
      pos2(1)=v(1)+pmtpos(1)
      pos2(2)=v(2)+pmtpos(2)
      pos2(3)=v(3)+pmtpos(3)
c
      v(3)=z3*cos(th)-x3*sin(th)
      vv1=z3*sin(th)+x3*cos(th)
c
      v(1)=vv1*cos(phi)-y3*sin(phi)
      v(2)=vv1*sin(phi)+y3*cos(phi)
c
      pos3(1)=v(1)+pmtpos(1)
      pos3(2)=v(2)+pmtpos(2)
      pos3(3)=v(3)+pmtpos(3)
c   for direction
      if(ihit.eq.1)then
        dir2(3)=c1*cos(th)-a1*sin(th)
        dd1=c1*sin(th)+a1*cos(th)
c
        dir2(1)=dd1*cos(phi)-b1*sin(phi)
        dir2(2)=dd1*sin(phi)+b1*cos(phi)
      else
        dir2(1)=0.
        dir2(2)=0.
        dir2(3)=0.
      endif
      return
      end
