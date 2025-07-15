      parameter(im=240,jm=121,km=37,im1=im+1,lm1=41,lm2=12)
      real uu(im,jm),vv(im,jm),fx(im,jm),fy(im,jm)
      real ua(im,jm),va(im,jm),a(im,jm,km)
      real um(im,jm),vm(im,jm)     

      real u(im1,jm),v(im1,jm)
      real u1(im1,jm),v1(im1,jm)
      real u2(im1,jm),v2(im1,jm)
      real u3(im1,jm),v3(im1,jm)
      real vor(im1,jm),psi(im1,jm)
      real div(im1,jm),xai(im1,jm)
      real*8 f(jm,im1),coslat(jm)
      real cor(im1,jm),t(im1,jm)
      real lat(jm)
      integer lmonth(12)
      data lmonth/31,28,31,30,31,30,31,31,30,31,30,31/

      open (11,file='height.8020.daily.37level.seq',
     &form='unformatted',
     &status='unknown',access='sequential')


      open (15,file='term1.8020.monthly.37level.seq',
     &form='unformatted',
     &status='unknown',access='sequential')

C       set the latitudes of grid points from north to south
C       as some subroutines prefer this order

      do j=1,jm
      lat(j)=90.-(j-1)*1.5
      enddo

C       define the Coriolis parameter: cor  (north to south)

      do i=1,im1 
      do j=1,jm
      cor(i,j)=2.*sin(lat(j)/180.*3.1415926535)*7.2921*0.00001   
      enddo
      enddo


      do l1=1,lm1
      do l2=1,lm2

C       this is the output field for monthly term1 values
      
      do k=1,km
      do j=1,jm
      do i=1,im
      a(i,j,k)=0.
      enddo
      enddo
      enddo

     
      do l3=1,lmonth(l2)

      do k=1,km

C       read in daily geopotential height field at each vertical level: 
C       uu (from south to north) 

      read(11)((uu(i,j),i=1,im),j=1,jm)

C       calculate geostrophic winds (fx,fy) from geopotential height uu
C       (from south to north)

      call gh2uv(uu,im,jm,fx,fy)


C       reshape fx&fy to u&v ordered from north to south   
C       as some subroutines prefer this order
   
      do i=1,im
      do j=1,jm
      u(i,j)=fx(i,jm+1-j)
      v(i,j)=fy(i,jm+1-j)
      enddo
      enddo

      do j=1,jm
      u(im1,j)=fx(1,jm+1-j)
      v(im1,j)=fy(1,jm+1-j)
      enddo

C       calculate relative vorticity (vor) from geostrophic winds (u,v, north to
C       south)

      call uv2vor
     o                (vor,
     i                 u,v,im,jm)

C       calculate absolutue vorticity vor + cor (north to south)

      do i=1,im1
      do j=1,jm
      vor(i,j)=vor(i,j)+cor(i,j)
      enddo
      enddo

C       calculate the gradient (u1, v1) of absolute vorticity (north to
C       south)

      call gradient(vor,u1,v1,im,jm)

C       calculate term1 (t) on a daily basis ( north to south)

      do i=1,im
      do j=1,jm
      t(i,j)=u(i,j)*u1(i,j)+v(i,j)*v1(i,j)
      enddo
      enddo

C       calculate monthly mean value of term1: a (north to south)

      do j=1,jm
      do i=1,im
      a(i,j,k)=a(i,j,k)+t(i,j)/real(lmonth(l2))
      enddo
      enddo

      enddo   
   
      enddo   

C       output monthly term1 values at 37 vertical levels ( from south
C       to north)

      do k=1,km
      write(15)((a(i,j,k),i=1,im),j=jm,1,-1)
      enddo


      enddo
      enddo
      end

      subroutine gradient(dat,u,v,nlo,nla)
      real       u   ( nlo+1, nla ), v   ( nlo+1, nla ),
     &           dat ( nlo+1, nla )
      real       phi ( nla ), cosphi(nla), sinphi(nla),
     &           delx( nla ), dely( nla )

      call const
     i          (nlo   , nla   ,
     o           a     , pai   ,
     o           phi   , cosphi, sinphi,
     o           dlamda, dphi  , delx  , dely , delyc)

      do 3000 j=2,nla-1
      do 2000 i=2,nlo
          u(i,j)=(dat(i+1,j)-dat(i-1,j))/delx(j)
          v(i,j)=-(dat(i,j+1)-dat(i,j-1)
     &)/dely(j)*cosphi(j)
          
 2000 continue
 3000 continue

      do 3001 j=2,nla-1
         u(1,j)=(dat(2,j)-dat(nlo,j))/delx(j)
         u(nlo+1,j) = u(1,j)
         v(1,j)=-(dat(1,j+1)-dat(1,j-1)
     &)/dely(j)*cosphi(j)
         v(nlo+1,j)=v(1,j)
 3001 continue

      do 2001 i=1,nlo+1
         u(i,1)  = u(i,2)
         u(i,nla)= u(i,nla-1)
         v(i,1)=v(i,2)
         v(i,nla)=v(i,nla-1)
 2001 continue
      return
      end

      subroutine xai2uv(dat,u,v,nlo,nla,cor)
      real       u   ( nlo+1, nla ), v   ( nlo+1, nla ),
     &           dat ( nlo+1, nla ), cor( nlo+1, nla )
      real       phi ( nla ), cosphi(nla), sinphi(nla),
     &           delx( nla ), dely( nla )

      call const
     i          (nlo   , nla   ,
     o           a     , pai   ,
     o           phi   , cosphi, sinphi,
     o           dlamda, dphi  , delx  , dely , delyc)

      do i=1,nlo+1
      do j=1,nla
      cor(i,j)=2*7.292*0.00001*sinphi(j)
      enddo
      enddo

      do 3000 j=2,nla-1
      do 2000 i=2,nlo
          v(i,j)=-(dat(i,j+1)-dat(i,j-1)
     &)/dely(j)*cosphi(j)
          u(i,j)=(dat(i+1,j)-dat(i-1,j))/delx(j)
 2000 continue
 3000 continue

      do 3001 j=2,nla-1
         v(1,j)=-(dat(1,j+1)-dat(1,j-1)
     &)/dely(j)*cosphi(j)
         v(nlo+1,j)=v(1,j)
         u(1,j)=(dat(2,j)-dat(nlo,j))/delx(j)
         u(nlo+1,j) = u(1,j)
 3001 continue

      do 2001 i=1,nlo+1
         v(i,1)  = v(i,2)
         v(i,nla)= v(i,nla-1)
         u(i,1)=u(i,2)
         u(i,nla)=u(i,nla-1)
 2001 continue
      return
      end

	subroutine gh2uv(h,ix,iy,u,v)
	real h(ix,iy),u(ix,iy),v(ix,iy)
cc
cc
	a=6.37e6
	omega=7.292e-5
	pi=3.1415926575
	dx=1.5*pi/180.0
	dy=1.5*pi/180.0
cc u,v
        do j=1,iy
	rlat=float(j-(iy+1)/2)*dy
	psin=sin(rlat)
	pcos=cos(rlat)
	ps=2.0*omega*psin

	do i=1,ix
	u(i,j)=-1.0/a*(h(i,j+1)-h(i,j-1))/ps/2.0/dy
	v(i,j)=1.0/a/pcos*(h(i+1,j)-h(i-1,j))/ps/2.0/dx
	enddo
	v(1,j)=1.0/a/pcos*(h(2,j)-h(ix,j))/ps/2.0/dx
	v(ix,j)=1.0/a/pcos*(h(1,j)-h(ix-1,j))/ps/2.0/dx
        u(1,j)=u(2,j)
        u(ix,j)=u(ix-1,j)
	enddo
cc
	return
	end
