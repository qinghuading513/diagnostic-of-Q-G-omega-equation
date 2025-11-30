      parameter(im=240,jm=121,km=41,lm=2)
      real r(im,jm,km),land(im,jm)
      real p(km),u(km),q(km)
      real s1(im,jm),s2(im,jm),s3(im,jm),s4(im,jm)
      DIMENSION X(km),Y(km),A(lm)
      DOUBLE PRECISION X,Y,A,DT1,DT2,DT3,B,C
      

      open(10,file='z850.annual.8021.dat',form='unformatted',
     &access='direct',recl=im*jm)
      open(12,file='z850.slope.41yr.dat',form='unformatted',
     &access='direct',recl=im*jm)
   
      irec=1
      do k=1,km
      read(10,rec=irec)((r(i,j,k),i=1,im),j=1,jm)
      irec=irec+1
      enddo


      do j=1,jm
      do i=1,im

      B=1.0
      ccc=0.0
      DO  k=1,km
      X(k)=B+(k-1)*1.
      ccc=ccc+x(k)/real(km) 
      Y(k)=r(i,j,k)
      q(k)=y(k)
      enddo

      ssx=0.
      do k=1,km
      ssx=ssx+(x(k)-ccc)**2
      enddo
 
      N=km
      M=lm
      CALL HPIR1(X,Y,A,N,M,DT1,DT2,DT3)

      c=0.0
      do k=1,km
      c=c+x(k)/real(km)
      enddo

      do k=1,km
      x(k)=x(k)-c
      enddo

      do k=1,km
      p(k)=a(1)+a(2)*x(k)
      enddo

      do k=1,km
      u(k)=y(k)-p(k)
      enddo

      std=0.
      do k=1,km
      std=std+u(k)**2
      enddo
      std=sqrt(std/real(km-2))
      s1(i,j)=(a(2)/std)*sqrt(km*(km-1)*(km+1)/12.)
      s2(i,j)=(a(2)/std)*sqrt(ssx)
      bb=0.
      call corr(q,p,km,bb)
      s3(i,j)=bb*sqrt(real(km-2.))/sqrt(1-bb**2)
      s4(i,j)=a(2)

      enddo
      enddo


      write(12,rec=1)((s1(i,j),i=1,im),j=1,jm)   ! trend significance 1
      write(12,rec=2)((s2(i,j),i=1,im),j=1,jm)   ! trend significance 2
      write(12,rec=3)((s3(i,j),i=1,im),j=1,jm)   ! trend significance 3
      write(12,rec=4)((s4(i,j),i=1,im),j=1,jm)   ! linear trend
     
      
      end
	SUBROUTINE HPIR1(X,Y,A,N,M,DT1,DT2,DT3)
	DIMENSION X(N),Y(N),A(M),S(20),T(20),B(20)
	DOUBLE PRECISION X,Y,A,S,T,B,DT1,DT2,DT3,
     *                   Z,D1,P,C,D2,G,Q,DT
	DO 5 I=1,M
5	A(I)=0.0
	IF (M.GT.N) M=N
	IF (M.GT.20) M=20
	Z=0.0
	DO 10 I=1,N
10	Z=Z+X(I)/N
	B(1)=1.0
	D1=N
	P=0.0
	C=0.0
	DO 20 I=1,N
	  P=P+(X(I)-Z)
	  C=C+Y(I)
20	CONTINUE
	C=C/D1
	P=P/D1
	A(1)=C*B(1)
	IF (M.GT.1) THEN
	  T(2)=1.0
	  T(1)=-P
	  D2=0.0
	  C=0.0
	  G=0.0
	  DO 30 I=1,N
	    Q=X(I)-Z-P
	    D2=D2+Q*Q
	    C=Y(I)*Q+C
	    G=(X(I)-Z)*Q*Q+G
30	  CONTINUE

	  C=C/D2
	  P=G/D2
	  Q=D2/D1
	  D1=D2
	  A(2)=C*T(2)
	  A(1)=C*T(1)+A(1)
	END IF
	DO 100 J=3,M
	  S(J)=T(J-1)
	  S(J-1)=-P*T(J-1)+T(J-2)
	  IF (J.GE.4) THEN
	    DO 40 K=J-2,2,-1
40	    S(K)=-P*T(K)+T(K-1)-Q*B(K)
	  END IF
	  S(1)=-P*T(1)-Q*B(1)
	  D2=0.0
	  C=0.0
	  G=0.0
	  DO 70 I=1,N
	    Q=S(J)
	    DO 60 K=J-1,1,-1
60	    Q=Q*(X(I)-Z)+S(K)
	    D2=D2+Q*Q
	    C=Y(I)*Q+C
	    G=(X(I)-Z)*Q*Q+G
70	  CONTINUE
	  C=C/D2
	  P=G/D2
	  Q=D2/D1
	  D1=D2
	  A(J)=C*S(J)
	  T(J)=S(J)
	  DO 80 K=J-1,1,-1
	    A(K)=C*S(K)+A(K)
	    B(K)=T(K)
	    T(K)=S(K)
80	  CONTINUE
100	CONTINUE
	DT1=0.0
	DT2=0.0
	DT3=0.0
	DO 120 I=1,N
	  Q=A(M)
	  DO 110 K=M-1,1,-1
110	  Q=Q*(X(I)-Z)+A(K)
	  DT=Q-Y(I)
	  IF (ABS(DT).GT.DT3) DT3=ABS(DT)
	  DT1=DT1+DT*DT
	  DT2=DT2+ABS(DT)
120	CONTINUE
	RETURN
	END

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

