c*************************************
      parameter(im=240,jm=121,iyear=41)
      real g1(im,jm,iyear),g2(im,jm,iyear),x(iyear),y(iyear),z
      real cor(im,jm)
      open(13,file='prec.8021.annual.dat',
     &form='unformatted',
     &access='direct',recl=im*jm)
           open(14,file='omega.8021.annual.ave.dat',form='unformatted',
     &access='direct',recl=im*jm)
      open(15,file='corr.prec-omega.dat',
     &form='unformatted',access='direct',recl=im*jm)


      irec=1
      do l=1,iyear
      read(13,rec=irec)((g1(i,j,l),i=1,im),j=1,jm)
      read(14,rec=irec)((g2(i,j,l),i=1,im),j=1,jm)
      irec=irec+1
      enddo



      do i=1,im
      do j=1,jm

      do l=1,iyear
      x(l)=g1(i,j,l)
      y(l)=g2(i,j,l)
      enddo

      call corr(x,y,iyear,z)
      cor(i,j)=z
      enddo
      enddo

      write(15,rec=1)((cor(i,j),i=1,im),j=1,jm)
      end

      subroutine corr(dex,u,lyear,cor)
      real dex(lyear),u(lyear)
      real a,b,up1,down1,down2,cor

      a=0.
      b=0.
      up1=0.
      down1=0.
      down2=0.

      do l=1,lyear
      a=a+dex(l)/real(lyear)
      b=b+u(l)/real(lyear)
      enddo

      do l=1,lyear
      dex(l)=dex(l)-a
      u(l)=u(l)-b
      enddo

      do l=1,lyear
      up1=up1+u(l)*dex(l)
      down1=down1+u(l)**2
      down2=down2+dex(l)**2
      enddo

      cor=up1/sqrt(down1)/sqrt(down2)
      end

