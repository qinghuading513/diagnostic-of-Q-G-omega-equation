      parameter(im=240,jm=121,km=37,lm1=12)
      real lat(jm),dy,dx
      real p(km)
      integer month1(12),month2(12),mday(12)
      real u(im,jm,km,124),v(im,jm,km,124),t(im,jm,km,124)
      real o(im,jm,km,124),h(im,jm,km,124)
      real ua(im,jm,km),va(im,jm,km),ta(im,jm,km)
      real oa(im,jm,km),ha(im,jm,km)
      real tb(im,jm,km)
      real ad1(im,jm,km),ad2(im,jm,km),ad3(im,jm,km)
      real td(im,jm,km),q(im,jm,km)
      real hd1(im,jm,km),hd2(im,jm,km)
      real vd1(im,jm,km),vd2(im,jm,km)

      data month1/31,28,31,30,31,30,31,31,30,31,30,31/
      data month2/31,29,31,30,31,30,31,31,30,31,30,31/
      data p/1,2,3,5,7,10,20,30,50,70,100,125,150,175,200,225,
     &250,300,350,400,450,500,550,600,650,700,750,775,800,825,
     &850,875,900,925,950,975,1000/

      open(10,file='uwind.8020.6hour.37level.dat',form='unformatted'
     &,access='direct',recl=im*jm)
      open(20,file='vwind.8020.6hour.37level.dat',form='unformatted'
     &,access='direct',recl=im*jm)
      open(30,file='tempt.8021.6hour.37level.dat',form='unformatted'
     &,access='direct',recl=im*jm)
      open(40,file='omega.8020.6hour.37level.dat',form='unformatted'
     &,access='direct',recl=im*jm)
      open(50,file='diabatic.month.1980-2020.seq',
     &form='unformatted',access='sequential')

C       set the latitudes of grid points from south to north
      
      do j=1,jm
      lat(j)=(-90.+(j-1)*1.5)/180.*3.1415926575
      enddo

      irec=1

      do ll=1980,2020

C       determine if it is a leap year

      if(mod(ll,4).ne.0) then
      do i=1,12
      mday(i)=month1(i)
      enddo
      else
      do i=1,12
      mday(i)=month2(i)
      enddo
      endif

      do l1=1,12

      do k=1,km
      do j=1,jm
      do i=1,im

      ua(i,j,k)=0.
      va(i,j,k)=0.
      ta(i,j,k)=0.
      oa(i,j,k)=0.
      ha(i,j,k)=0.
      tb(i,j,k)=0.
      ad1(i,j,k)=0.
      ad2(i,j,k)=0.
      ad3(i,j,k)=0.

      do l2=1,mday(l1)*4
      u(i,j,k,l2)=0.
      v(i,j,k,l2)=0.
      t(i,j,k,l2)=0.
      o(i,j,k,l2)=0.
      h(i,j,k,l2)=0.
      enddo

      enddo
      enddo
      enddo


      do l2=1,mday(l1)*4

C        read in all 37 levels of 6-h fields 
C        u: zonal winds   ( south to north)
C        v: meridional winds  ( south to north)
C        t: temperature  ( south to north)
C        o: omega  ( south to north)

      do k=1,km
      read(10,rec=irec)((u(i,j,k,l2),i=1,im),j=1,jm)
      read(20,rec=irec)((v(i,j,k,l2),i=1,im),j=1,jm)
      read(30,rec=irec)((t(i,j,k,l2),i=1,im),j=1,jm)
      read(40,rec=irec)((o(i,j,k,l2),i=1,im),j=1,jm)
      irec=irec+1
      enddo

      enddo

C       read in one additional time step for temperature ( the first 6-h T in the next
C       month) so that we can calculate delta T for the current month

      irec1=irec

      do k=1,km
      read(30,rec=irec1)((tb(i,j,k),i=1,im),j=1,jm)
      irec1=irec1+1
      enddo

C     calculate the temperature difference over the entire one month 

      do k=1,km
      do j=1,jm
      do i=1,im
      tb(i,j,k)=tb(i,j,k)-t(i,j,k,1)
      enddo
      enddo
      enddo      

C       calculate monthly means of all five fields (u,v,t,o,h)
C       h:potential temperature    

      do k=1,km
      do j=1,jm
      do i=1,im
      do l2=1,mday(l1)*4
      h(i,j,k,l2)=t(i,j,k,l2)*((1000./p(k))**0.286)
      ua(i,j,k)=ua(i,j,k)+u(i,j,k,l2)/real(mday(l1)*4)
      va(i,j,k)=va(i,j,k)+v(i,j,k,l2)/real(mday(l1)*4)
      ta(i,j,k)=ta(i,j,k)+t(i,j,k,l2)/real(mday(l1)*4)
      oa(i,j,k)=oa(i,j,k)+o(i,j,k,l2)/real(mday(l1)*4)
      ha(i,j,k)=ha(i,j,k)+h(i,j,k,l2)/real(mday(l1)*4)
      enddo
      enddo
      enddo
      enddo

C       calculate deviations from the monthly means for each 6h data
C       point

      do k=1,km
      do j=1,jm
      do i=1,im
      do l2=1,mday(l1)*4
      u(i,j,k,l2)=u(i,j,k,l2)-ua(i,j,k)
      v(i,j,k,l2)=v(i,j,k,l2)-va(i,j,k)
      t(i,j,k,l2)=t(i,j,k,l2)-ta(i,j,k)
      o(i,j,k,l2)=o(i,j,k,l2)-oa(i,j,k)
      h(i,j,k,l2)=h(i,j,k,l2)-ha(i,j,k)
      enddo
      enddo
      enddo
      enddo

C       calculate the overbar terms (monthly means) in the Eq. from
C       Nigam et al. 2000

      do k=1,km
      do j=1,jm
      do i=1,im
      do l2=1,mday(l1)*4
      ad1(i,j,k)=ad1(i,j,k)+u(i,j,k,l2)*h(i,j,k,l2)/real(mday(l1)*4)
      ad2(i,j,k)=ad2(i,j,k)+v(i,j,k,l2)*h(i,j,k,l2)/real(mday(l1)*4)
      ad3(i,j,k)=ad3(i,j,k)+o(i,j,k,l2)*h(i,j,k,l2)/real(mday(l1)*4)
      enddo
      enddo
      enddo
      enddo

C       calculate the horizontal and vertical gradients

      do k=1,km
      do i=1,im
      do j=1,jm

      kk1=k-1
      kk2=k+1

      ii1=i-1
      ii2=i+1

      jj1=j-1
      jj2=j+1

      if(i.eq.1) then
      ii1=im
      endif
      if(i.eq.im) then
      ii2=1
      endif
      if(j.eq.1) then
      jj1=1
      jj2=3
      endif
      if(j.eq.jm) then
      jj2=jm
      jj1=jm-2
      endif
      if(k.eq.km) then
      kk2=km
      kk1=km-1
      endif
      if(k.eq.1) then
      kk1=1
      kk2=2
      endif


      dx=6378000.*2.*3.1415926/360.*1.5*cos(lat((jj1+jj2)/2))
      dy=6378000.*2.*3.1415926/360.*1.5

C       term A in the eq. from Nigam et al. 2000

      td(i,j,k)= tb(i,j,k)/mday(l1)/86400.

C       term B in the eq. from Nigam et al. 2000

      hd1(i,j,k)= (ta(ii2,j,k)-ta(ii1,j,k))/(2*dx)*ua(i,j,k)+
     &(ta(i,jj2,k)-ta(i,jj1,k))/(2*dy)*va(i,j,k)

C       term C in the eq. from Nigam et al. 2000

      vd1(i,j,k)=oa(i,j,k)*(ha(i,j,kk2)-ha(i,j,kk1))/
     &(p(kk2)*100.-p(kk1)*100.)*
     &((p(k)/1000.)**0.286)

C       term D in the eq. from Nigam et al. 2000

      hd2(i,j,k)= ((ad1(ii2,j,k)-ad1(ii1,j,k))/(2*dx)+
     &(ad2(i,jj2,k)-ad2(i,jj1,k))/(2*dy))*((p(k)/1000.)**0.286)

C       term E in the eq. from Nigam et al. 2000

      vd2(i,j,k)=(( ad3(i,j,kk2)-ad3(i,j,kk1) )/
     &(p(kk2)*100.-p(kk1)*100.) )*((p(k)/1000.)**0.286)

C       sum all 5 terms ( A to E) to calculate the Q term

      q(i,j,k)=td(i,j,k)+hd1(i,j,k)+vd1(i,j,k)+hd2(i,j,k)+vd2(i,j,k)

      enddo
      enddo
      enddo

C       output monthly mean diabatic heating term  ( south to north)

      do k=1,km
      write(50)((q(i,j,k),i=1,im),j=1,jm)
      enddo

      enddo 
      enddo

      end

