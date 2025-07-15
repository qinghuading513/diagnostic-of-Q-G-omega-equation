c===================================================================
      subroutine const 
     i		     (nlo   , nla   ,
     o                a     , pai   ,
     o		      phi   , cosphi, sinphi,
     o                dlamda, dphi  , delx  , dely , delyc)
c===================================================================
c     purpose
c          set parameters for the basic dynamical diagnostics of 
c          global atmospheric circulation 
c          NOTE: !!! from north to south pole
c     *************************************************************
c     arguments
c     nlo:number of grid points in longitudinal direction(input)
c     nla:number of grid points in latitudinal direction(input)
c     a:earth radius  
c     pai: 
c     phi(nla):latitudes of the grid points
c     cosphi(nla):cos(phi(nla))
c     sinphi(nla):sin(phi(nla)) 
c     dlamda:longitude interval
c     dphi:latitude interval
c     delx(nla):a*cos(phi(nla))*2*dlamda
c     dely(nla):a*cos(phi(nla))*2*dphi
c     delyc:a*2*dphi
c
c===================================================================
      real  phi ( nla ), cosphi(nla), sinphi(nla),
     &      delx( nla ), dely( nla )
c===================================================================
c
      a=6.37122e+6
      pai=4.0*atan(1.0)
      dlamda = 2.*pai/float(nlo)
      dphi   = pai/float(nla-1)
      delyc  = a*2.0*dphi
      do j=1,nla
         phi(j)   = pai/2.-float(j-1)*dphi
         cosphi(j)   = cos(phi(j))
         sinphi(j)   = sin(phi(j))
         delx(j)  = a*cosphi(j)*2.0*dlamda
         dely(j)  = a*cosphi(j)*2.0*dphi
      end do
c     
      return
      end

      SUBROUTINE COSGEN (N,IJUMP,FNUM,FDEN,A)
C     =================================================================
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION       A(*)
C
      PI = PIMACH(DUM)
      IF (N .EQ. 0) GO TO 105
      IF (IJUMP .EQ. 1) GO TO 103
      K3 = N/IJUMP+1
      K4 = K3-1
      PIBYN = PI/DFLOAT(N+IJUMP)
      DO 102 K=1,IJUMP
         K1 = (K-1)*K3
         K5 = (K-1)*K4
         DO 101 I=1,K4
            X = K1+I
            K2 = K5+I
            A(K2) = -2.D0*COS(X*PIBYN)
  101    CONTINUE
  102 CONTINUE
      GO TO 105
  103 CONTINUE
      NP1 = N+1
      Y = PI/(DFLOAT(N)+FDEN)
      DO 104 I=1,N
         X = DFLOAT(NP1-I)-FNUM
         A(I) = 2.D0*COS(X*Y)
  104 CONTINUE
  105 CONTINUE
      RETURN
      END
      FUNCTION EPMACH (DUM)
C     =================================================================
C     THIS PROGRAM COMPUTES AN APPROXIMATE MACHIINE EPSILON (ACCURACY)
C
      IMPLICIT REAL*8(A-H,O-Z)
      COMMON /VALUE/  V
      EPS = 1.D0
  101 CONTINUE
        EPS = EPS/10.D0
        CALL STORE (EPS+1.D0)
      IF (V-1. .GT. 0.D0) GO TO 101
      EPMACH = 100.D0*EPS
      RETURN
      END
      SUBROUTINE GENBUN (NPEROD,N,MPEROD,M,A,B,C,IDIMY,Y,IERROR,W)
C     =================================================================
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION       Y(IDIMY,*)
      DIMENSION       W(*)       ,B(*)       ,A(*)       ,C(*)
      IERROR = 0
      IF (M .LE. 2) IERROR = 1
      IF (N .LE. 2) IERROR = 2
      IF (IDIMY .LT. M) IERROR = 3
      IF (NPEROD.LT.0 .OR. NPEROD.GT.4) IERROR = 4
      IF (MPEROD.LT.0 .OR. MPEROD.GT.1) IERROR = 5
      IF (MPEROD .EQ. 1) GO TO 102
      DO 101 I=2,M
         IF (A(I) .NE. C(1)) GO TO 103
         IF (C(I) .NE. C(1)) GO TO 103
         IF (B(I) .NE. B(1)) GO TO 103
  101 CONTINUE
      GO TO 104
  102 IF (A(1).NE.0.D0 .OR. C(M).NE.0.D0) IERROR = 7
      GO TO 104
  103 IERROR = 6
  104 IF (IERROR .NE. 0) RETURN
      MP1 = M+1
      IWBA = MP1
      IWBB = IWBA+M
      IWBC = IWBB+M
      IWB2 = IWBC+M
      IWB3 = IWB2+M
      IWW1 = IWB3+M
      IWW2 = IWW1+M
      IWW3 = IWW2+M
      IWD = IWW3+M
      IWTCOS = IWD+M
      IWP = IWTCOS+4*N
      DO 106 I=1,M
         K = IWBA+I-1
         W(K) = -A(I)
         K = IWBC+I-1
         W(K) = -C(I)
         K = IWBB+I-1
         W(K) = 2.D0-B(I)
         DO 105 J=1,N
            Y(I,J) = -Y(I,J)
  105    CONTINUE
  106 CONTINUE
      MP = MPEROD+1
      NP = NPEROD+1
      GO TO (114,107),MP
  107 GO TO (108,109,110,111,123),NP
  108 CALL POISP2 (M,N,W(IWBA),W(IWBB),W(IWBC),Y,IDIMY,W,W(IWB2),
     1             W(IWB3),W(IWW1),W(IWW2),W(IWW3),W(IWD),W(IWTCOS),
     2             W(IWP))
      GO TO 112
  109 CALL POISD2 (M,N,1,W(IWBA),W(IWBB),W(IWBC),Y,IDIMY,W,W(IWW1),
     1             W(IWD),W(IWTCOS),W(IWP))
      GO TO 112
  110 CALL POISN2 (M,N,1,2,W(IWBA),W(IWBB),W(IWBC),Y,IDIMY,W,W(IWB2),
     1             W(IWB3),W(IWW1),W(IWW2),W(IWW3),W(IWD),W(IWTCOS),
     2             W(IWP))
      GO TO 112
  111 CALL POISN2 (M,N,1,1,W(IWBA),W(IWBB),W(IWBC),Y,IDIMY,W,W(IWB2),
     1             W(IWB3),W(IWW1),W(IWW2),W(IWW3),W(IWD),W(IWTCOS),
     2             W(IWP))
  112 IPSTOR = W(IWW1)
      IREV = 2
      IF (NPEROD .EQ. 4) GO TO 124
  113 GO TO (127,133),MP
  114 CONTINUE
C
C     REORDER UNKNOWNS WHEN MP =0
C
      MH = (M+1)/2
      MHM1 = MH-1
      MODD = 1
      IF (MH*2 .EQ. M) MODD = 2
      DO 119 J=1,N
         DO 115 I=1,MHM1
            MHPI = MH+I
            MHMI = MH-I
            W(I) = Y(MHMI,J)-Y(MHPI,J)
            W(MHPI) = Y(MHMI,J)+Y(MHPI,J)
  115    CONTINUE
         W(MH) = 2.D0*Y(MH,J)
         GO TO (117,116),MODD
  116    W(M) = 2.D0*Y(M,J)
  117    CONTINUE
         DO 118 I=1,M
            Y(I,J) = W(I)
  118    CONTINUE
  119 CONTINUE
      K = IWBC+MHM1-1
      I = IWBA+MHM1
      W(K) = 0.D0
      W(I) = 0.D0
      W(K+1) = 2.D0*W(K+1)
      GO TO (120,121),MODD
  120 CONTINUE
      K = IWBB+MHM1-1
      W(K) = W(K)-W(I-1)
      W(IWBC-1) = W(IWBC-1)+W(IWBB-1)
      GO TO 122
  121 W(IWBB-1) = W(K+1)
  122 CONTINUE
      GO TO 107
C
C     REVERSE COLUMNS WHEN NPEROD = 4.
C
  123 IREV = 1
      NBY2 = N/2
  124 DO 126 J=1,NBY2
         MSKIP = N+1-J
         DO 125 I=1,M
            A1 = Y(I,J)
            Y(I,J) = Y(I,MSKIP)
            Y(I,MSKIP) = A1
  125    CONTINUE
  126 CONTINUE
      GO TO (110,113),IREV
  127 CONTINUE
      DO 132 J=1,N
         DO 128 I=1,MHM1
            MHMI = MH-I
            MHPI = MH+I
            W(MHMI) = .5D0*(Y(MHPI,J)+Y(I,J))
            W(MHPI) = .5D0*(Y(MHPI,J)-Y(I,J))
  128    CONTINUE
         W(MH) = .5D0*Y(MH,J)
         GO TO (130,129),MODD
  129    W(M) = .5D0*Y(M,J)
  130    CONTINUE
         DO 131 I=1,M
            Y(I,J) = W(I)
  131    CONTINUE
  132 CONTINUE
  133 CONTINUE
C
C     RETURN STORAGE REQUIREMENTS FOR W ARRAY.
C
      W(1) = IPSTOR+IWP-1
      RETURN
      END
      SUBROUTINE HWSSS1 (TS,TF,M,MBDCND,BDTS,BDTF,PS,PF,N,NBDCND,BDPS,
     1                   BDPF,ELMBDA,F,IDIMF,PERTRB,AM,BM,CM,SN,SS,
     2                   SINT,D)
C     =================================================================
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION       F(IDIMF,*) ,BDTS(*)    ,BDTF(*)    ,BDPS(*)    ,
     1                BDPF(*)    ,AM(*)      ,BM(*)      ,CM(*)      ,
     2                SS(*)      ,SN(*)      ,D(*)       ,SINT(*)
C
      PI = PIMACH(DUM)
      TPI = PI+PI
      HPI = PI/2.D0
      MP1 = M+1
      NP1 = N+1
      FN = N
      FM = M
      DTH = (TF-TS)/FM
      HDTH = DTH/2.D0
      TDT = DTH+DTH
      DPHI = (PF-PS)/FN
      TDP = DPHI+DPHI
      DPHI2 = DPHI*DPHI
      EDP2 = ELMBDA*DPHI2
      DTH2 = DTH*DTH
      CP = 4.D0/(FN*DTH2)
      WP = FN*SIN(HDTH)/4.D0
      DO 102 I=1,MP1
         FIM1 = I-1
         THETA = FIM1*DTH+TS
         SINT(I) = SIN(THETA)
         IF (SINT(I)) 101,102,101
  101    T1 = 1.D0/(DTH2*SINT(I))
         AM(I) = T1*SIN(THETA-HDTH)
         CM(I) = T1*SIN(THETA+HDTH)
         BM(I) = -AM(I)-CM(I)+ELMBDA
  102 CONTINUE
      INP = 0
      ISP = 0
C
C BOUNDARY CONDITION AT THETA=TS
C
      MBR = MBDCND+1
      GO TO (103,104,104,105,105,106,106,104,105,106),MBR
  103 ITS = 1
      GO TO 107
  104 AT = AM(2)
      ITS = 2
      GO TO 107
  105 AT = AM(1)
      ITS = 1
      CM(1) = AM(1)+CM(1)
      GO TO 107
  106 AT = AM(2)
      INP = 1
      ITS = 2
C
C BOUNDARY CONDITION THETA=TF
C
  107 GO TO (108,109,110,110,109,109,110,111,111,111),MBR
  108 ITF = M
      GO TO 112
  109 CT = CM(M)
      ITF = M
      GO TO 112
  110 CT = CM(M+1)
      AM(M+1) = AM(M+1)+CM(M+1)
      ITF = M+1
      GO TO 112
  111 ITF = M
      ISP = 1
      CT = CM(M)
C
C COMPUTE HOMOGENEOUS SOLUTION WITH SOLUTION AT POLE EQUAL TO ONE
C
  112 ITSP = ITS+1
      ITFM = ITF-1
      WTS = SINT(ITS+1)*AM(ITS+1)/CM(ITS)
      WTF = SINT(ITF-1)*CM(ITF-1)/AM(ITF)
      MUNK = ITF-ITS+1
      IF (ISP) 116,116,113
  113 D(ITS) = CM(ITS)/BM(ITS)
      DO 114 I=ITSP,M
         D(I) = CM(I)/(BM(I)-AM(I)*D(I-1))
  114 CONTINUE
      SS(M) = -D(M)
      IID = M-ITS
      DO 115 II=1,IID
         I = M-II
         SS(I) = -D(I)*SS(I+1)
  115 CONTINUE
      SS(M+1) = 1.D0
  116 IF (INP) 120,120,117
  117 SN(1) = 1.D0
      D(ITF) = AM(ITF)/BM(ITF)
      IID = ITF-2
      DO 118 II=1,IID
         I = ITF-II
         D(I) = AM(I)/(BM(I)-CM(I)*D(I+1))
  118 CONTINUE
      SN(2) = -D(2)
      DO 119 I=3,ITF
         SN(I) = -D(I)*SN(I-1)
  119 CONTINUE
C
C BOUNDARY CONDITIONS AT PHI=PS
C
  120 NBR = NBDCND+1
      WPS = 1.D0
      WPF = 1.D0
      GO TO (121,122,122,123,123),NBR
  121 JPS = 1
      GO TO 124
  122 JPS = 2
      GO TO 124
  123 JPS = 1
      WPS = .5D0
C
C BOUNDARY CONDITION AT PHI=PF
C
  124 GO TO (125,126,127,127,126),NBR
  125 JPF = N
      GO TO 128
  126 JPF = N
      GO TO 128
  127 WPF = .5D0
      JPF = N+1
  128 JPSP = JPS+1
      JPFM = JPF-1
      NUNK = JPF-JPS+1
      FJJ = JPFM-JPSP+1
C
C SCALE COEFFICIENTS FOR SUBROUTINE GENBUN
C
      DO 129 I=ITS,ITF
         CF = DPHI2*SINT(I)*SINT(I)
         AM(I) = CF*AM(I)
         BM(I) = CF*BM(I)
         CM(I) = CF*CM(I)
  129 CONTINUE
      AM(ITS) = 0.D0
      CM(ITF) = 0.D0
      ISING = 0
      GO TO (130,138,138,130,138,138,130,138,130,130),MBR
  130 GO TO (131,138,138,131,138),NBR
  131 IF (ELMBDA) 138,132,132
  132 ISING = 1
      SUM = WTS*WPS+WTS*WPF+WTF*WPS+WTF*WPF
      IF (INP) 134,134,133
  133 SUM = SUM+WP
  134 IF (ISP) 136,136,135
  135 SUM = SUM+WP
  136 SUM1 = 0.D0
      DO 137 I=ITSP,ITFM
         SUM1 = SUM1+SINT(I)
  137 CONTINUE
      SUM = SUM+FJJ*(SUM1+WTS+WTF)
      SUM = SUM+(WPS+WPF)*SUM1
      HNE = SUM
  138 GO TO (146,142,142,144,144,139,139,142,144,139),MBR
  139 IF (NBDCND-3) 146,140,146
  140 YHLD = F(1,JPS)-4.D0/(FN*DPHI*DTH2)*(BDPF(2)-BDPS(2))
      DO 141 J=1,NP1
         F(1,J) = YHLD
  141 CONTINUE
      GO TO 146
  142 DO 143 J=JPS,JPF
         F(2,J) = F(2,J)-AT*F(1,J)
  143 CONTINUE
      GO TO 146
  144 DO 145 J=JPS,JPF
         F(1,J) = F(1,J)+TDT*BDTS(J)*AT
  145 CONTINUE
  146 GO TO (154,150,152,152,150,150,152,147,147,147),MBR
  147 IF (NBDCND-3) 154,148,154
  148 YHLD = F(M+1,JPS)-4.D0/(FN*DPHI*DTH2)*(BDPF(M)-BDPS(M))
      DO 149 J=1,NP1
         F(M+1,J) = YHLD
  149 CONTINUE
      GO TO 154
  150 DO 151 J=JPS,JPF
         F(M,J) = F(M,J)-CT*F(M+1,J)
  151 CONTINUE
      GO TO 154
  152 DO 153 J=JPS,JPF
         F(M+1,J) = F(M+1,J)-TDT*BDTF(J)*CT
  153 CONTINUE
  154 GO TO (159,155,155,157,157),NBR
  155 DO 156 I=ITS,ITF
         F(I,2) = F(I,2)-F(I,1)/(DPHI2*SINT(I)*SINT(I))
  156 CONTINUE
      GO TO 159
  157 DO 158 I=ITS,ITF
         F(I,1) = F(I,1)+TDP*BDPS(I)/(DPHI2*SINT(I)*SINT(I))
  158 CONTINUE
  159 GO TO (164,160,162,162,160),NBR
  160 DO 161 I=ITS,ITF
         F(I,N) = F(I,N)-F(I,N+1)/(DPHI2*SINT(I)*SINT(I))
  161 CONTINUE
      GO TO 164
  162 DO 163 I=ITS,ITF
         F(I,N+1) = F(I,N+1)-TDP*BDPF(I)/(DPHI2*SINT(I)*SINT(I))
  163 CONTINUE
  164 CONTINUE
      PERTRB = 0.D0
      IF (ISING) 165,176,165
  165 SUM = WTS*WPS*F(ITS,JPS)+WTS*WPF*F(ITS,JPF)+WTF*WPS*F(ITF,JPS)+
     1      WTF*WPF*F(ITF,JPF)
      IF (INP) 167,167,166
  166 SUM = SUM+WP*F(1,JPS)
  167 IF (ISP) 169,169,168
  168 SUM = SUM+WP*F(M+1,JPS)
  169 DO 171 I=ITSP,ITFM
         SUM1 = 0.D0
         DO 170 J=JPSP,JPFM
            SUM1 = SUM1+F(I,J)
  170    CONTINUE
         SUM = SUM+SINT(I)*SUM1
  171 CONTINUE
      SUM1 = 0.D0
      SUM2 = 0.D0
      DO 172 J=JPSP,JPFM
         SUM1 = SUM1+F(ITS,J)
         SUM2 = SUM2+F(ITF,J)
  172 CONTINUE
      SUM = SUM+WTS*SUM1+WTF*SUM2
      SUM1 = 0.D0
      SUM2 = 0.D0
      DO 173 I=ITSP,ITFM
         SUM1 = SUM1+SINT(I)*F(I,JPS)
         SUM2 = SUM2+SINT(I)*F(I,JPF)
  173 CONTINUE
      SUM = SUM+WPS*SUM1+WPF*SUM2
      PERTRB = SUM/HNE
      DO 175 J=1,NP1
         DO 174 I=1,MP1
            F(I,J) = F(I,J)-PERTRB
  174    CONTINUE
  175 CONTINUE
C
C SCALE RIGHT SIDE FOR SUBROUTINE GENBUN
C
  176 DO 178 I=ITS,ITF
         CF = DPHI2*SINT(I)*SINT(I)
         DO 177 J=JPS,JPF
            F(I,J) = CF*F(I,J)
  177    CONTINUE
  178 CONTINUE
      CALL GENBUN (NBDCND,NUNK,1,MUNK,AM(ITS),BM(ITS),CM(ITS),IDIMF,
     1             F(ITS,JPS),IERROR,D)
      IF (ISING) 186,186,179
  179 IF (INP) 183,183,180
  180 IF (ISP) 181,181,186
  181 DO 182 J=1,NP1
         F(1,J) = 0.D0
  182 CONTINUE
      GO TO 209
  183 IF (ISP) 186,186,184
  184 DO 185 J=1,NP1
         F(M+1,J) = 0.D0
  185 CONTINUE
      GO TO 209
  186 IF (INP) 193,193,187
  187 SUM = WPS*F(ITS,JPS)+WPF*F(ITS,JPF)
      DO 188 J=JPSP,JPFM
         SUM = SUM+F(ITS,J)
  188 CONTINUE
      DFN = CP*SUM
      DNN = CP*((WPS+WPF+FJJ)*(SN(2)-1.D0))+ELMBDA
      DSN = CP*(WPS+WPF+FJJ)*SN(M)
      IF (ISP) 189,189,194
  189 CNP = (F(1,1)-DFN)/DNN
      DO 191 I=ITS,ITF
         HLD = CNP*SN(I)
         DO 190 J=JPS,JPF
            F(I,J) = F(I,J)+HLD
  190    CONTINUE
  191 CONTINUE
      DO 192 J=1,NP1
         F(1,J) = CNP
  192 CONTINUE
      GO TO 209
  193 IF (ISP) 209,209,194
  194 SUM = WPS*F(ITF,JPS)+WPF*F(ITF,JPF)
      DO 195 J=JPSP,JPFM
         SUM = SUM+F(ITF,J)
  195 CONTINUE
      DFS = CP*SUM
      DSS = CP*((WPS+WPF+FJJ)*(SS(M)-1.D0))+ELMBDA
      DNS = CP*(WPS+WPF+FJJ)*SS(2)
      IF (INP) 196,196,200
  196 CSP = (F(M+1,1)-DFS)/DSS
      DO 198 I=ITS,ITF
         HLD = CSP*SS(I)
         DO 197 J=JPS,JPF
            F(I,J) = F(I,J)+HLD
  197    CONTINUE
  198 CONTINUE
      DO 199 J=1,NP1
         F(M+1,J) = CSP
  199 CONTINUE
      GO TO 209
  200 RTN = F(1,1)-DFN
      RTS = F(M+1,1)-DFS
      IF (ISING) 202,202,201
  201 CSP = 0.D0
      CNP = RTN/DNN
      GO TO 205
  202 IF (ABS(DNN)-ABS(DSN)) 204,204,203
  203 DEN = DSS-DNS*DSN/DNN
      RTS = RTS-RTN*DSN/DNN
      CSP = RTS/DEN
      CNP = (RTN-CSP*DNS)/DNN
      GO TO 205
  204 DEN = DNS-DSS*DNN/DSN
      RTN = RTN-RTS*DNN/DSN
      CSP = RTN/DEN
      CNP = (RTS-DSS*CSP)/DSN
  205 DO 207 I=ITS,ITF
         HLD = CNP*SN(I)+CSP*SS(I)
         DO 206 J=JPS,JPF
            F(I,J) = F(I,J)+HLD
  206    CONTINUE
  207 CONTINUE
      DO 208 J=1,NP1
         F(1,J) = CNP
         F(M+1,J) = CSP
  208 CONTINUE
  209 IF (NBDCND) 212,210,212
  210 DO 211 I=1,MP1
         F(I,JPF+1) = F(I,JPS)
  211 CONTINUE
  212 RETURN
      END
      SUBROUTINE HWSSSP (TS,TF,M,MBDCND,BDTS,BDTF,PS,PF,N,NBDCND,BDPS,
     1                   BDPF,ELMBDA,F,IDIMF,PERTRB,IERROR,W)
C     =================================================================
C     SUBROUTINE HWSSSP SOLVES A FINITE DIFFERENCE APPROXIMATION TO THE
C     HELMHOLTZ EQUATION IN SPHERICAL COORDINATES AND ON THE SURFACE OF
C     THE UNIT SPHERE (RADIUS OF 1):
C
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION       F(IDIMF,*) ,BDTS(*)    ,BDTF(*)    ,BDPS(*)    ,
     1                BDPF(*)    ,W(*)
      NBR = NBDCND+1
      PI = PIMACH(DUM)
      TPI = 2.D0*PI
      IERROR = 0
      IF (TS.LT.0.D0 .OR. TF.GT.PI) IERROR = 1
      IF (TS .GE. TF) IERROR = 2
      IF (MBDCND.LT.1 .OR. MBDCND.GT.9) IERROR = 3
      IF (PS.LT.0.D0 .OR. PF.GT.TPI) IERROR = 4
      IF (PS .GE. PF) IERROR = 5
      IF (N .LT. 5) IERROR = 6
      IF (M .LT. 5) IERROR = 7
      IF (NBDCND.LT.0 .OR. NBDCND.GT.4) IERROR = 8
      IF (ELMBDA .GT. 0.D0) IERROR = 9
      IF (IDIMF .LT. M+1) IERROR = 10
      IF ((NBDCND.EQ.1 .OR. NBDCND.EQ.2 .OR. NBDCND.EQ.4) .AND.
     1    MBDCND.GE.5) IERROR = 11
      IF (TS.EQ.0. .AND.
     1    (MBDCND.EQ.3 .OR. MBDCND.EQ.4 .OR. MBDCND.EQ.8)) IERROR = 12
      IF (TF.EQ.PI .AND.
     1    (MBDCND.EQ.2 .OR. MBDCND.EQ.3 .OR. MBDCND.EQ.6)) IERROR = 13
      IF ((MBDCND.EQ.5 .OR. MBDCND.EQ.6 .OR. MBDCND.EQ.9) .AND.
     1    TS.NE.0.) IERROR = 14
      IF (MBDCND.GE.7 .AND. TF.NE.PI) IERROR = 15
      IF (IERROR.NE.0 .AND. IERROR.NE.9) RETURN
      CALL HWSSS1 (TS,TF,M,MBDCND,BDTS,BDTF,PS,PF,N,NBDCND,BDPS,BDPF,
     1             ELMBDA,F,IDIMF,PERTRB,W,W(M+2),W(2*M+3),W(3*M+4),
     2             W(4*M+5),W(5*M+6),W(6*M+7))
      W(1) = W(6*M+7)+DFLOAT(6*(M+1))
      RETURN
      END
C=================================================================
      SUBROUTINE INVLAP (NLO,NLA,F,COSLAT,PSI,VOR)
C=================================================================
C
C CALCULATE  PSI = INVERSE LAPLACIAN OF VOR
C USING SUBROUTINE HWSSSP OF NCAR LIBRARY
C " >>NCARL DOUBLE " REQUIRED
C
C=================================================================
      PARAMETER (NW=7000)
      REAL*4     VOR(NLO+1,NLA),PSI(NLO+1,NLA)
      REAL*8     F(NLA,NLO+1),W(NW),COSLAT(NLA)
      COMMON     W
      REAL*8 TS,TF,BDTS,BDTF,PS,PF,BDPS,BDPF,ELAMDA,PTB,PI,DUM
      REAL*8 COLAT,DLAT,GLM,WEIT
      REAL*8 PIMACH
C=================================================================
      PI = PIMACH(DUM)
      TS = 0.D0
      TF = PI
      PS = 0.D0
      PF = 2.D0*PI
      ELAMDA = 0.D0
      MBDCND = 9
      NBDCND = 0
C
      DLAT = PI/DFLOAT(NLA-1)
      DO 50 LA = 1, NLA
        COLAT = (LA-1)*DLAT
        COSLAT(LA) = SIN(COLAT)
   50 CONTINUE
C
      DO 200 LO = 1, NLO+1
        DO 100 LA = 1, NLA
          F(LA,LO) = VOR(LO,LA)
  100   CONTINUE
  200 CONTINUE
C
      CALL HWSSSP (TS,TF,NLA-1,MBDCND,BDTS,BDTF,
     &             PS,PF,NLO,  NBDCND,BDPS,BDPF,
     &             ELAMDA, F, NLA, PTB, IER, W)
C
      GLM   = 0.D0
      WEIT  = 0.D0
      DO 240 LA = 1, NLA
        DO 220 LO = 1, NLO
          GLM  = GLM + F(LA,LO)*COSLAT(LA)
          WEIT = WEIT + COSLAT(LA)
  220   CONTINUE
  240 CONTINUE
      GLM = GLM/WEIT
      WRITE (6,*) '*INVLAP* HWSSSP RETURNED, IERROR=',IER,'PERTRB=',PTB,
     &                     'GLOBAL MEAN RESULT = ',GLM
C
      DO 400 LA = 1, NLA
        DO 300 LO = 1, NLO+1
          PSI(LO,LA) = (F(LA,LO)-GLM)
  300   CONTINUE
  400 CONTINUE
      RETURN
      END
      SUBROUTINE MERGE (TCOS,I1,M1,I2,M2,I3)
C     =================================================================
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION       TCOS(*)
C
      J1 = 1
      J2 = 1
      J = I3
      IF (M1 .EQ. 0) GO TO 107
      IF (M2 .EQ. 0) GO TO 104
  101 J = J+1
      L = J1+I1
      X = TCOS(L)
      L = J2+I2
      Y = TCOS(L)
      IF (X-Y) 102,102,103
  102 TCOS(J) = X
      J1 = J1+1
      IF (J1 .GT. M1) GO TO 106
      GO TO 101
  103 TCOS(J) = Y
      J2 = J2+1
      IF (J2 .LE. M2) GO TO 101
      IF (J1 .GT. M1) GO TO 109
  104 K = J-J1+1
      DO 105 J=J1,M1
         M = K+J
         L = J+I1
         TCOS(M) = TCOS(L)
  105 CONTINUE
      GO TO 109
  106 CONTINUE
      IF (J2 .GT. M2) GO TO 109
  107 K = J-J2+1
      DO 108 J=J2,M2
         M = K+J
         L = J+I2
         TCOS(M) = TCOS(L)
  108 CONTINUE
  109 CONTINUE
      RETURN
      END
      FUNCTION PIMACH (DUM)
C     =================================================================
      IMPLICIT REAL*8(A-H,O-Z)
C
C     THIS SUBPROGRAM SUPPLIES THE VALUE OF THE CONSTANT PI CORRECT TO
C     MACHINE PRECISION WHERE
C
C     PI=3.1415926535897932384626433832795028841971693993751058209749446
C
      PIMACH = 3.14159265358979323846D0
      RETURN
      END
      SUBROUTINE POISD2 (MR,NR,ISTAG,BA,BB,BC,Q,IDIMQ,B,W,D,TCOS,P)
C     =================================================================
C     SUBROUTINE TO SOLVE POISSON"S EQUATION FOR DIRICHLET BOUNDARY
C     CONDITIONS.
C
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION       Q(IDIMQ,*) ,BA(*)      ,BB(*)      ,BC(*)      ,
     1                TCOS(*)    ,B(*)       ,D(*)       ,W(*)       ,
     2                P(*)
      M = MR
      N = NR
      JSH = 0
      FI = 1.D0/DFLOAT(ISTAG)
      IP = -M
      IPSTOR = 0
      GO TO (101,102),ISTAG
  101 KR = 0
      IRREG = 1
      IF (N .GT. 1) GO TO 106
      TCOS(1) = 0.D0
      GO TO 103
  102 KR = 1
      JSTSAV = 1
      IRREG = 2
      IF (N .GT. 1) GO TO 106
      TCOS(1) = -1.D0
  103 DO 104 I=1,M
         B(I) = Q(I,1)
  104 CONTINUE
      CALL TRIX (1,0,M,BA,BB,BC,B,TCOS,D,W)
      DO 105 I=1,M
         Q(I,1) = B(I)
  105 CONTINUE
      GO TO 183
  106 LR = 0
      DO 107 I=1,M
         P(I) = 0.D0
  107 CONTINUE
      NUN = N
      JST = 1
      JSP = N
C
C     IRREG = 1 WHEN NO IRREGULARITIES HAVE OCCURRED, OTHERWISE IT IS 2.
C
  108 L = 2*JST
      NODD = 2-2*((NUN+1)/2)+NUN
C
C     NODD = 1 WHEN NUN IS ODD, OTHERWISE IT IS 2.
C
      GO TO (110,109),NODD
  109 JSP = JSP-L
      GO TO 111
  110 JSP = JSP-JST
      IF (IRREG .NE. 1) JSP = JSP-L
  111 CONTINUE
C
C     REGULAR REDUCTION
C
      CALL COSGEN (JST,1,0.5D0,0.0D0,TCOS)
      IF (L .GT. JSP) GO TO 118
      DO 117 J=L,JSP,L
         JM1 = J-JSH
         JP1 = J+JSH
         JM2 = J-JST
         JP2 = J+JST
         JM3 = JM2-JSH
         JP3 = JP2+JSH
         IF (JST .NE. 1) GO TO 113
         DO 112 I=1,M
            B(I) = 2.D0*Q(I,J)
            Q(I,J) = Q(I,JM2)+Q(I,JP2)
  112    CONTINUE
         GO TO 115
  113    DO 114 I=1,M
            T = Q(I,J)-Q(I,JM1)-Q(I,JP1)+Q(I,JM2)+Q(I,JP2)
            B(I) = T+Q(I,J)-Q(I,JM3)-Q(I,JP3)
            Q(I,J) = T
  114    CONTINUE
  115    CONTINUE
         CALL TRIX (JST,0,M,BA,BB,BC,B,TCOS,D,W)
         DO 116 I=1,M
            Q(I,J) = Q(I,J)+B(I)
  116    CONTINUE
  117 CONTINUE
C
C     REDUCTION FOR LAST UNKNOWN
C
  118 GO TO (119,136),NODD
  119 GO TO (152,120),IRREG
C
C     ODD NUMBER OF UNKNOWNS
C
  120 JSP = JSP+L
      J = JSP
      JM1 = J-JSH
      JP1 = J+JSH
      JM2 = J-JST
      JP2 = J+JST
      JM3 = JM2-JSH
      GO TO (123,121),ISTAG
  121 CONTINUE
      IF (JST .NE. 1) GO TO 123
      DO 122 I=1,M
         B(I) = Q(I,J)
         Q(I,J) = 0.D0
  122 CONTINUE
      GO TO 130
  123 GO TO (124,126),NODDPR
  124 DO 125 I=1,M
         IP1 = IP+I
         B(I) = .5D0*(Q(I,JM2)-Q(I,JM1)-Q(I,JM3))+P(IP1)+Q(I,J)
  125 CONTINUE
      GO TO 128
  126 DO 127 I=1,M
      B(I) = .5D0*(Q(I,JM2)-Q(I,JM1)-Q(I,JM3))+Q(I,JP2)-Q(I,JP1)+Q(I,J)
  127 CONTINUE
  128 DO 129 I=1,M
         Q(I,J) = .5D0*(Q(I,J)-Q(I,JM1)-Q(I,JP1))
  129 CONTINUE
  130 CALL TRIX (JST,0,M,BA,BB,BC,B,TCOS,D,W)
      IP = IP+M
      IPSTOR = MAX0(IPSTOR,IP+M)
      DO 131 I=1,M
         IP1 = IP+I
         P(IP1) = Q(I,J)+B(I)
         B(I) = Q(I,JP2)+P(IP1)
  131 CONTINUE
      IF (LR .NE. 0) GO TO 133
      DO 132 I=1,JST
         KRPI = KR+I
         TCOS(KRPI) = TCOS(I)
  132 CONTINUE
      GO TO 134
  133 CONTINUE
      CALL COSGEN (LR,JSTSAV,0.D0,FI,TCOS(JST+1))
      CALL MERGE (TCOS,0,JST,JST,LR,KR)
  134 CONTINUE
      CALL COSGEN (KR,JSTSAV,0.0D0,FI,TCOS)
      CALL TRIX (KR,KR,M,BA,BB,BC,B,TCOS,D,W)
      DO 135 I=1,M
         IP1 = IP+I
         Q(I,J) = Q(I,JM2)+B(I)+P(IP1)
  135 CONTINUE
      LR = KR
      KR = KR+L
      GO TO 152
C
C     EVEN NUMBER OF UNKNOWNS
C
  136 JSP = JSP+L
      J = JSP
      JM1 = J-JSH
      JP1 = J+JSH
      JM2 = J-JST
      JP2 = J+JST
      JM3 = JM2-JSH
      GO TO (137,138),IRREG
  137 CONTINUE
      JSTSAV = JST
      IDEG = JST
      KR = L
      GO TO 139
  138 CALL COSGEN (KR,JSTSAV,0.0D0,FI,TCOS)
      CALL COSGEN (LR,JSTSAV,0.0D0,FI,TCOS(KR+1))
      IDEG = KR
      KR = KR+JST
  139 IF (JST .NE. 1) GO TO 141
      IRREG = 2
      DO 140 I=1,M
         B(I) = Q(I,J)
         Q(I,J) = Q(I,JM2)
  140 CONTINUE
      GO TO 150
  141 DO 142 I=1,M
         B(I) = Q(I,J)+.5D0*(Q(I,JM2)-Q(I,JM1)-Q(I,JM3))
  142 CONTINUE
      GO TO (143,145),IRREG
  143 DO 144 I=1,M
         Q(I,J) = Q(I,JM2)+.5D0*(Q(I,J)-Q(I,JM1)-Q(I,JP1))
  144 CONTINUE
      IRREG = 2
      GO TO 150
  145 CONTINUE
      GO TO (146,148),NODDPR
  146 DO 147 I=1,M
         IP1 = IP+I
         Q(I,J) = Q(I,JM2)+P(IP1)
  147 CONTINUE
      IP = IP-M
      GO TO 150
  148 DO 149 I=1,M
         Q(I,J) = Q(I,JM2)+Q(I,J)-Q(I,JM1)
  149 CONTINUE
  150 CALL TRIX (IDEG,LR,M,BA,BB,BC,B,TCOS,D,W)
      DO 151 I=1,M
         Q(I,J) = Q(I,J)+B(I)
  151 CONTINUE
  152 NUN = NUN/2
      NODDPR = NODD
      JSH = JST
      JST = 2*JST
      IF (NUN .GE. 2) GO TO 108
C
C     START SOLUTION.
C
      J = JSP
      DO 153 I=1,M
         B(I) = Q(I,J)
  153 CONTINUE
      GO TO (154,155),IRREG
  154 CONTINUE
      CALL COSGEN (JST,1,0.5D0,0.0D0,TCOS)
      IDEG = JST
      GO TO 156
  155 KR = LR+JST
      CALL COSGEN (KR,JSTSAV,0.0D0,FI,TCOS)
      CALL COSGEN (LR,JSTSAV,0.0D0,FI,TCOS(KR+1))
      IDEG = KR
  156 CONTINUE
      CALL TRIX (IDEG,LR,M,BA,BB,BC,B,TCOS,D,W)
      JM1 = J-JSH
      JP1 = J+JSH
      GO TO (157,159),IRREG
  157 DO 158 I=1,M
         Q(I,J) = .5D0*(Q(I,J)-Q(I,JM1)-Q(I,JP1))+B(I)
  158 CONTINUE
      GO TO 164
  159 GO TO (160,162),NODDPR
  160 DO 161 I=1,M
         IP1 = IP+I
         Q(I,J) = P(IP1)+B(I)
  161 CONTINUE
      IP = IP-M
      GO TO 164
  162 DO 163 I=1,M
         Q(I,J) = Q(I,J)-Q(I,JM1)+B(I)
  163 CONTINUE
  164 CONTINUE
C
C     START BACK SUBSTITUTION.
C
      JST = JST/2
      JSH = JST/2
      NUN = 2*NUN
      IF (NUN .GT. N) GO TO 183
      DO 182 J=JST,N,L
         JM1 = J-JSH
         JP1 = J+JSH
         JM2 = J-JST
         JP2 = J+JST
         IF (J .GT. JST) GO TO 166
         DO 165 I=1,M
            B(I) = Q(I,J)+Q(I,JP2)
  165    CONTINUE
         GO TO 170
  166    IF (JP2 .LE. N) GO TO 168
         DO 167 I=1,M
            B(I) = Q(I,J)+Q(I,JM2)
  167    CONTINUE
         IF (JST .LT. JSTSAV) IRREG = 1
         GO TO (170,171),IRREG
  168    DO 169 I=1,M
            B(I) = Q(I,J)+Q(I,JM2)+Q(I,JP2)
  169    CONTINUE
  170    CONTINUE
         CALL COSGEN (JST,1,0.5D0,0.0D0,TCOS)
         IDEG = JST
         JDEG = 0
         GO TO 172
  171    IF (J+L .GT. N) LR = LR-JST
         KR = JST+LR
         CALL COSGEN (KR,JSTSAV,0.0D0,FI,TCOS)
         CALL COSGEN (LR,JSTSAV,0.0D0,FI,TCOS(KR+1))
         IDEG = KR
         JDEG = LR
  172    CONTINUE
         CALL TRIX (IDEG,JDEG,M,BA,BB,BC,B,TCOS,D,W)
         IF (JST .GT. 1) GO TO 174
         DO 173 I=1,M
            Q(I,J) = B(I)
  173    CONTINUE
         GO TO 182
  174    IF (JP2 .GT. N) GO TO 177
  175    DO 176 I=1,M
            Q(I,J) = .5D0*(Q(I,J)-Q(I,JM1)-Q(I,JP1))+B(I)
  176    CONTINUE
         GO TO 182
  177    GO TO (175,178),IRREG
  178    IF (J+JSH .GT. N) GO TO 180
         DO 179 I=1,M
            IP1 = IP+I
            Q(I,J) = B(I)+P(IP1)
  179    CONTINUE
         IP = IP-M
         GO TO 182
  180    DO 181 I=1,M
            Q(I,J) = B(I)+Q(I,J)-Q(I,JM1)
  181    CONTINUE
  182 CONTINUE
      L = L/2
      GO TO 164
  183 CONTINUE
C
C     RETURN STORAGE REQUIREMENTS FOR P VECTORS.
C
      W(1) = IPSTOR
      RETURN
      END
      SUBROUTINE POISN2 (M,N,ISTAG,MIXBND,A,BB,C,Q,IDIMQ,B,B2,B3,W,W2,
     1                   W3,D,TCOS,P)
C     =================================================================
C     SUBROUTINE TO SOLVE POISSON"S EQUATION WITH NEUMANN BOUNDARY
C     CONDITIONS.
C
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION       A(*)       ,BB(*)      ,C(*)       ,Q(IDIMQ,*) ,
     1                B(*)       ,B2(*)      ,B3(*)      ,W(*)       ,
     2                W2(*)      ,W3(*)      ,D(*)       ,TCOS(*)    ,
     3                K(4)       ,P(*)
      EQUIVALENCE     (K(1),K1)  ,(K(2),K2)  ,(K(3),K3)  ,(K(4),K4)
      FISTAG = 3-ISTAG
      FNUM = 1.D0/DFLOAT(ISTAG)
      FDEN = 0.5D0*DFLOAT(ISTAG-1)
      MR = M
      IP = -MR
      IPSTOR = 0
      I2R = 1
      JR = 2
      NR = N
      NLAST = N
      KR = 1
      LR = 0
      GO TO (101,103),ISTAG
  101 CONTINUE
      DO 102 I=1,MR
         Q(I,N) = .5D0*Q(I,N)
  102 CONTINUE
      GO TO (103,104),MIXBND
  103 IF (N .LE. 3) GO TO 155
  104 CONTINUE
      JR = 2*I2R
      NROD = 1
      IF ((NR/2)*2 .EQ. NR) NROD = 0
      GO TO (105,106),MIXBND
  105 JSTART = 1
      GO TO 107
  106 JSTART = JR
      NROD = 1-NROD
  107 CONTINUE
      JSTOP = NLAST-JR
      IF (NROD .EQ. 0) JSTOP = JSTOP-I2R
      CALL COSGEN (I2R,1,0.5D0,0.0D0,TCOS)
      I2RBY2 = I2R/2
      IF (JSTOP .GE. JSTART) GO TO 108
      J = JR
      GO TO 116
  108 CONTINUE
C
C     REGULAR REDUCTION.
C
      DO 115 J=JSTART,JSTOP,JR
         JP1 = J+I2RBY2
         JP2 = J+I2R
         JP3 = JP2+I2RBY2
         JM1 = J-I2RBY2
         JM2 = J-I2R
         JM3 = JM2-I2RBY2
         IF (J .NE. 1) GO TO 109
         JM1 = JP1
         JM2 = JP2
         JM3 = JP3
  109    CONTINUE
         IF (I2R .NE. 1) GO TO 111
         IF (J .EQ. 1) JM2 = JP2
         DO 110 I=1,MR
            B(I) = 2.D0*Q(I,J)
            Q(I,J) = Q(I,JM2)+Q(I,JP2)
  110    CONTINUE
         GO TO 113
  111    CONTINUE
         DO 112 I=1,MR
            FI = Q(I,J)
            Q(I,J) = Q(I,J)-Q(I,JM1)-Q(I,JP1)+Q(I,JM2)+Q(I,JP2)
            B(I) = FI+Q(I,J)-Q(I,JM3)-Q(I,JP3)
  112    CONTINUE
  113    CONTINUE
         CALL TRIX (I2R,0,MR,A,BB,C,B,TCOS,D,W)
         DO 114 I=1,MR
            Q(I,J) = Q(I,J)+B(I)
  114    CONTINUE
C
C     END OF REDUCTION FOR REGULAR UNKNOWNS.
C
  115 CONTINUE
C
C     BEGIN SPECIAL REDUCTION FOR LAST UNKNOWN.
C
      J = JSTOP+JR
  116 NLAST = J
      JM1 = J-I2RBY2
      JM2 = J-I2R
      JM3 = JM2-I2RBY2
      IF (NROD .EQ. 0) GO TO 128
C
C     ODD NUMBER OF UNKNOWNS
C
      IF (I2R .NE. 1) GO TO 118
      DO 117 I=1,MR
         B(I) = FISTAG*Q(I,J)
         Q(I,J) = Q(I,JM2)
  117 CONTINUE
      GO TO 126
  118 DO 119 I=1,MR
         B(I) = Q(I,J)+.5D0*(Q(I,JM2)-Q(I,JM1)-Q(I,JM3))
  119 CONTINUE
      IF (NRODPR .NE. 0) GO TO 121
      DO 120 I=1,MR
         II = IP+I
         Q(I,J) = Q(I,JM2)+P(II)
  120 CONTINUE
      IP = IP-MR
      GO TO 123
  121 CONTINUE
      DO 122 I=1,MR
         Q(I,J) = Q(I,J)-Q(I,JM1)+Q(I,JM2)
  122 CONTINUE
  123 IF (LR .EQ. 0) GO TO 124
      CALL COSGEN (LR,1,0.5D0,FDEN,TCOS(KR+1))
      GO TO 126
  124 CONTINUE
      DO 125 I=1,MR
         B(I) = FISTAG*B(I)
  125 CONTINUE
  126 CONTINUE
      CALL COSGEN (KR,1,0.5D0,FDEN,TCOS)
      CALL TRIX (KR,LR,MR,A,BB,C,B,TCOS,D,W)
      DO 127 I=1,MR
         Q(I,J) = Q(I,J)+B(I)
  127 CONTINUE
      KR = KR+I2R
      GO TO 151
  128 CONTINUE
C
C     EVEN NUMBER OF UNKNOWNS
C
      JP1 = J+I2RBY2
      JP2 = J+I2R
      IF (I2R .NE. 1) GO TO 135
      DO 129 I=1,MR
         B(I) = Q(I,J)
  129 CONTINUE
      CALL TRIX (1,0,MR,A,BB,C,B,TCOS,D,W)
      IP = 0
      IPSTOR = MR
      GO TO (133,130),ISTAG
  130 DO 131 I=1,MR
         P(I) = B(I)
         B(I) = B(I)+Q(I,N)
  131 CONTINUE
      TCOS(1) = 1.D0
      TCOS(2) = 0.D0
      CALL TRIX (1,1,MR,A,BB,C,B,TCOS,D,W)
      DO 132 I=1,MR
         Q(I,J) = Q(I,JM2)+P(I)+B(I)
  132 CONTINUE
      GO TO 150
  133 CONTINUE
      DO 134 I=1,MR
         P(I) = B(I)
         Q(I,J) = Q(I,JM2)+2.D0*Q(I,JP2)+3.D0*B(I)
  134 CONTINUE
      GO TO 150
  135 CONTINUE
      DO 136 I=1,MR
         B(I) = Q(I,J)+.5D0*(Q(I,JM2)-Q(I,JM1)-Q(I,JM3))
  136 CONTINUE
      IF (NRODPR .NE. 0) GO TO 138
      DO 137 I=1,MR
         II = IP+I
         B(I) = B(I)+P(II)
  137 CONTINUE
      GO TO 140
  138 CONTINUE
      DO 139 I=1,MR
         B(I) = B(I)+Q(I,JP2)-Q(I,JP1)
  139 CONTINUE
  140 CONTINUE
      CALL TRIX (I2R,0,MR,A,BB,C,B,TCOS,D,W)
      IP = IP+MR
      IPSTOR = MAX0(IPSTOR,IP+MR)
      DO 141 I=1,MR
         II = IP+I
         P(II) = B(I)+.5D0*(Q(I,J)-Q(I,JM1)-Q(I,JP1))
         B(I) = P(II)+Q(I,JP2)
  141 CONTINUE
      IF (LR .EQ. 0) GO TO 142
      CALL COSGEN (LR,1,0.5D0,FDEN,TCOS(I2R+1))
      CALL MERGE (TCOS,0,I2R,I2R,LR,KR)
      GO TO 144
  142 DO 143 I=1,I2R
         II = KR+I
         TCOS(II) = TCOS(I)
  143 CONTINUE
  144 CALL COSGEN (KR,1,0.5D0,FDEN,TCOS)
      IF (LR .NE. 0) GO TO 145
      GO TO (146,145),ISTAG
  145 CONTINUE
      CALL TRIX (KR,KR,MR,A,BB,C,B,TCOS,D,W)
      GO TO 148
  146 CONTINUE
      DO 147 I=1,MR
         B(I) = FISTAG*B(I)
  147 CONTINUE
  148 CONTINUE
      DO 149 I=1,MR
         II = IP+I
         Q(I,J) = Q(I,JM2)+P(II)+B(I)
  149 CONTINUE
  150 CONTINUE
      LR = KR
      KR = KR+JR
  151 CONTINUE
      GO TO (152,153),MIXBND
  152 NR = (NLAST-1)/JR+1
      IF (NR .LE. 3) GO TO 155
      GO TO 154
  153 NR = NLAST/JR
      IF (NR .LE. 1) GO TO 192
  154 I2R = JR
      NRODPR = NROD
      GO TO 104
  155 CONTINUE
C
C      BEGIN SOLUTION
C
      J = 1+JR
      JM1 = J-I2R
      JP1 = J+I2R
      JM2 = NLAST-I2R
      IF (NR .EQ. 2) GO TO 184
      IF (LR .NE. 0) GO TO 170
      IF (N .NE. 3) GO TO 161
C
C     CASE N = 3.
C
      GO TO (156,168),ISTAG
  156 CONTINUE
      DO 157 I=1,MR
         B(I) = Q(I,2)
  157 CONTINUE
      TCOS(1) = 0.D0
      CALL TRIX (1,0,MR,A,BB,C,B,TCOS,D,W)
      DO 158 I=1,MR
         Q(I,2) = B(I)
         B(I) = 4.D0*B(I)+Q(I,1)+2.D0*Q(I,3)
  158 CONTINUE
      TCOS(1) = -2.D0
      TCOS(2) = 2.D0
      I1 = 2
      I2 = 0
      CALL TRIX (I1,I2,MR,A,BB,C,B,TCOS,D,W)
      DO 159 I=1,MR
         Q(I,2) = Q(I,2)+B(I)
         B(I) = Q(I,1)+2.D0*Q(I,2)
  159 CONTINUE
      TCOS(1) = 0.D0
      CALL TRIX (1,0,MR,A,BB,C,B,TCOS,D,W)
      DO 160 I=1,MR
         Q(I,1) = B(I)
  160 CONTINUE
      JR = 1
      I2R = 0
      GO TO 194
C
C     CASE N = 2**P+1
C
  161 CONTINUE
      GO TO (162,170),ISTAG
  162 CONTINUE
      DO 163 I=1,MR
         B(I) = Q(I,J)+.5D0*Q(I,1)-Q(I,JM1)+Q(I,NLAST)-Q(I,JM2)
  163 CONTINUE
      CALL COSGEN (JR,1,0.5D0,0.0D0,TCOS)
      CALL TRIX (JR,0,MR,A,BB,C,B,TCOS,D,W)
      DO 164 I=1,MR
         Q(I,J) = .5D0*(Q(I,J)-Q(I,JM1)-Q(I,JP1))+B(I)
         B(I) = Q(I,1)+2.D0*Q(I,NLAST)+4.D0*Q(I,J)
  164 CONTINUE
      JR2 = 2*JR
      CALL COSGEN (JR,1,0.0D0,0.0D0,TCOS)
      DO 165 I=1,JR
         I1 = JR+I
         I2 = JR+1-I
         TCOS(I1) = -TCOS(I2)
  165 CONTINUE
      CALL TRIX (JR2,0,MR,A,BB,C,B,TCOS,D,W)
      DO 166 I=1,MR
         Q(I,J) = Q(I,J)+B(I)
         B(I) = Q(I,1)+2.D0*Q(I,J)
  166 CONTINUE
      CALL COSGEN (JR,1,0.5D0,0.0D0,TCOS)
      CALL TRIX (JR,0,MR,A,BB,C,B,TCOS,D,W)
      DO 167 I=1,MR
         Q(I,1) = .5D0*Q(I,1)-Q(I,JM1)+B(I)
  167 CONTINUE
      GO TO 194
C
C     CASE OF GENERAL N WITH NR = 3 .
C
  168 DO 169 I=1,MR
         B(I) = Q(I,2)
         Q(I,2) = 0.D0
         B2(I) = Q(I,3)
         B3(I) = Q(I,1)
  169 CONTINUE
      JR = 1
      I2R = 0
      J = 2
      GO TO 177
  170 CONTINUE
      DO 171 I=1,MR
         B(I) = .5D0*Q(I,1)-Q(I,JM1)+Q(I,J)
  171 CONTINUE
      IF (NROD .NE. 0) GO TO 173
      DO 172 I=1,MR
         II = IP+I
         B(I) = B(I)+P(II)
  172 CONTINUE
      GO TO 175
  173 DO 174 I=1,MR
         B(I) = B(I)+Q(I,NLAST)-Q(I,JM2)
  174 CONTINUE
  175 CONTINUE
      DO 176 I=1,MR
         T = .5D0*(Q(I,J)-Q(I,JM1)-Q(I,JP1))
         Q(I,J) = T
         B2(I) = Q(I,NLAST)+T
         B3(I) = Q(I,1)+2.D0*T
  176 CONTINUE
  177 CONTINUE
      K1 = KR+2*JR-1
      K2 = KR+JR
      TCOS(K1+1) = -2.D0
      K4 = K1+3-ISTAG
      CALL COSGEN (K2+ISTAG-2,1,0.0D0,FNUM,TCOS(K4))
      K4 = K1+K2+1
      CALL COSGEN (JR-1,1,0.0D0,1.0D0,TCOS(K4))
      CALL MERGE (TCOS,K1,K2,K1+K2,JR-1,0)
      K3 = K1+K2+LR
      CALL COSGEN (JR,1,0.5D0,0.0D0,TCOS(K3+1))
      K4 = K3+JR+1
      CALL COSGEN (KR,1,0.5D0,FDEN,TCOS(K4))
      CALL MERGE (TCOS,K3,JR,K3+JR,KR,K1)
      IF (LR .EQ. 0) GO TO 178
      CALL COSGEN (LR,1,0.5D0,FDEN,TCOS(K4))
      CALL MERGE (TCOS,K3,JR,K3+JR,LR,K3-LR)
      CALL COSGEN (KR,1,0.5D0,FDEN,TCOS(K4))
  178 K3 = KR
      K4 = KR
      CALL TRI3 (MR,A,BB,C,K,B,B2,B3,TCOS,D,W,W2,W3)
      DO 179 I=1,MR
         B(I) = B(I)+B2(I)+B3(I)
  179 CONTINUE
      TCOS(1) = 2.D0
      CALL TRIX (1,0,MR,A,BB,C,B,TCOS,D,W)
      DO 180 I=1,MR
         Q(I,J) = Q(I,J)+B(I)
         B(I) = Q(I,1)+2.D0*Q(I,J)
  180 CONTINUE
      CALL COSGEN (JR,1,0.5D0,0.0D0,TCOS)
      CALL TRIX (JR,0,MR,A,BB,C,B,TCOS,D,W)
      IF (JR .NE. 1) GO TO 182
      DO 181 I=1,MR
         Q(I,1) = B(I)
  181 CONTINUE
      GO TO 194
  182 CONTINUE
      DO 183 I=1,MR
         Q(I,1) = .5D0*Q(I,1)-Q(I,JM1)+B(I)
  183 CONTINUE
      GO TO 194
  184 CONTINUE
      IF (N .NE. 2) GO TO 188
C
C     CASE  N = 2
C
      DO 185 I=1,MR
         B(I) = Q(I,1)
  185 CONTINUE
      TCOS(1) = 0.D0
      CALL TRIX (1,0,MR,A,BB,C,B,TCOS,D,W)
      DO 186 I=1,MR
         Q(I,1) = B(I)
         B(I) = 2.D0*(Q(I,2)+B(I))*FISTAG
  186 CONTINUE
      TCOS(1) = -FISTAG
      TCOS(2) = 2.D0
      CALL TRIX (2,0,MR,A,BB,C,B,TCOS,D,W)
      DO 187 I=1,MR
         Q(I,1) = Q(I,1)+B(I)
  187 CONTINUE
      JR = 1
      I2R = 0
      GO TO 194
  188 CONTINUE
C
C     CASE OF GENERAL N AND NR = 2 .
C
      DO 189 I=1,MR
         II = IP+I
         B3(I) = 0.D0
         B(I) = Q(I,1)+2.D0*P(II)
         Q(I,1) = .5D0*Q(I,1)-Q(I,JM1)
         B2(I) = 2.D0*(Q(I,1)+Q(I,NLAST))
  189 CONTINUE
      K1 = KR+JR-1
      TCOS(K1+1) = -2.D0
      K4 = K1+3-ISTAG
      CALL COSGEN (KR+ISTAG-2,1,0.0D0,FNUM,TCOS(K4))
      K4 = K1+KR+1
      CALL COSGEN (JR-1,1,0.0D0,1.0D0,TCOS(K4))
      CALL MERGE (TCOS,K1,KR,K1+KR,JR-1,0)
      CALL COSGEN (KR,1,0.5D0,FDEN,TCOS(K1+1))
      K2 = KR
      K4 = K1+K2+1
      CALL COSGEN (LR,1,0.5D0,FDEN,TCOS(K4))
      K3 = LR
      K4 = 0
      CALL TRI3 (MR,A,BB,C,K,B,B2,B3,TCOS,D,W,W2,W3)
      DO 190 I=1,MR
         B(I) = B(I)+B2(I)
  190 CONTINUE
      TCOS(1) = 2.D0
      CALL TRIX (1,0,MR,A,BB,C,B,TCOS,D,W)
      DO 191 I=1,MR
         Q(I,1) = Q(I,1)+B(I)
  191 CONTINUE
      GO TO 194
  192 DO 193 I=1,MR
         B(I) = Q(I,NLAST)
  193 CONTINUE
      GO TO 196
  194 CONTINUE
C
C     START BACK SUBSTITUTION.
C
      J = NLAST-JR
      DO 195 I=1,MR
         B(I) = Q(I,NLAST)+Q(I,J)
  195 CONTINUE
  196 JM2 = NLAST-I2R
      IF (JR .NE. 1) GO TO 198
      DO 197 I=1,MR
         Q(I,NLAST) = 0.D0
  197 CONTINUE
      GO TO 202
  198 CONTINUE
      IF (NROD .NE. 0) GO TO 200
      DO 199 I=1,MR
         II = IP+I
         Q(I,NLAST) = P(II)
  199 CONTINUE
      IP = IP-MR
      GO TO 202
  200 DO 201 I=1,MR
         Q(I,NLAST) = Q(I,NLAST)-Q(I,JM2)
  201 CONTINUE
  202 CONTINUE
      CALL COSGEN (KR,1,0.5D0,FDEN,TCOS)
      CALL COSGEN (LR,1,0.5D0,FDEN,TCOS(KR+1))
      IF (LR .NE. 0) GO TO 204
      DO 203 I=1,MR
         B(I) = FISTAG*B(I)
  203 CONTINUE
  204 CONTINUE
      CALL TRIX (KR,LR,MR,A,BB,C,B,TCOS,D,W)
      DO 205 I=1,MR
         Q(I,NLAST) = Q(I,NLAST)+B(I)
  205 CONTINUE
      NLASTP = NLAST
  206 CONTINUE
      JSTEP = JR
      JR = I2R
      I2R = I2R/2
      IF (JR .EQ. 0) GO TO 222
      GO TO (207,208),MIXBND
  207 JSTART = 1+JR
      GO TO 209
  208 JSTART = JR
  209 CONTINUE
      KR = KR-JR
      IF (NLAST+JR .GT. N) GO TO 210
      KR = KR-JR
      NLAST = NLAST+JR
      JSTOP = NLAST-JSTEP
      GO TO 211
  210 CONTINUE
      JSTOP = NLAST-JR
  211 CONTINUE
      LR = KR-JR
      CALL COSGEN (JR,1,0.5D0,0.0D0,TCOS)
      DO 221 J=JSTART,JSTOP,JSTEP
         JM2 = J-JR
         JP2 = J+JR
         IF (J .NE. JR) GO TO 213
         DO 212 I=1,MR
            B(I) = Q(I,J)+Q(I,JP2)
  212    CONTINUE
         GO TO 215
  213    CONTINUE
         DO 214 I=1,MR
            B(I) = Q(I,J)+Q(I,JM2)+Q(I,JP2)
  214    CONTINUE
  215    CONTINUE
         IF (JR .NE. 1) GO TO 217
         DO 216 I=1,MR
            Q(I,J) = 0.D0
  216    CONTINUE
         GO TO 219
  217    CONTINUE
         JM1 = J-I2R
         JP1 = J+I2R
         DO 218 I=1,MR
            Q(I,J) = .5D0*(Q(I,J)-Q(I,JM1)-Q(I,JP1))
  218    CONTINUE
  219    CONTINUE
         CALL TRIX (JR,0,MR,A,BB,C,B,TCOS,D,W)
         DO 220 I=1,MR
            Q(I,J) = Q(I,J)+B(I)
  220    CONTINUE
  221 CONTINUE
      NROD = 1
      IF (NLAST+I2R .LE. N) NROD = 0
      IF (NLASTP .NE. NLAST) GO TO 194
      GO TO 206
  222 CONTINUE
C
C     RETURN STORAGE REQUIREMENTS FOR P VECTORS.
C
      W(1) = IPSTOR
      RETURN
      END
      SUBROUTINE POISP2 (M,N,A,BB,C,Q,IDIMQ,B,B2,B3,W,W2,W3,D,TCOS,P)
C     =================================================================
C     SUBROUTINE TO SOLVE POISSON EQUATION WITH PERIODIC BOUNDARY
C     CONDITIONS.
C
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION       A(*)       ,BB(*)      ,C(*)       ,Q(IDIMQ,*) ,
     1                B(*)       ,B2(*)      ,B3(*)      ,W(*)       ,
     2                W2(*)      ,W3(*)      ,D(*)       ,TCOS(*)    ,
     3                P(*)
      MR = M
      NR = (N+1)/2
      NRM1 = NR-1
      IF (2*NR .NE. N) GO TO 107
C
C     EVEN NUMBER OF UNKNOWNS
C
      DO 102 J=1,NRM1
         NRMJ = NR-J
         NRPJ = NR+J
         DO 101 I=1,MR
            S = Q(I,NRMJ)-Q(I,NRPJ)
            T = Q(I,NRMJ)+Q(I,NRPJ)
            Q(I,NRMJ) = S
            Q(I,NRPJ) = T
  101    CONTINUE
  102 CONTINUE
      DO 103 I=1,MR
         Q(I,NR) = 2.D0*Q(I,NR)
         Q(I,N) = 2.D0*Q(I,N)
  103 CONTINUE
      CALL POISD2 (MR,NRM1,1,A,BB,C,Q,IDIMQ,B,W,D,TCOS,P)
      IPSTOR = W(1)
      CALL POISN2 (MR,NR+1,1,1,A,BB,C,Q(1,NR),IDIMQ,B,B2,B3,W,W2,W3,D,
     1             TCOS,P)
      IPSTOR = MAX0(IPSTOR,INT(W(1)))
      DO 105 J=1,NRM1
         NRMJ = NR-J
         NRPJ = NR+J
         DO 104 I=1,MR
            S = .5D0*(Q(I,NRPJ)+Q(I,NRMJ))
            T = .5D0*(Q(I,NRPJ)-Q(I,NRMJ))
            Q(I,NRMJ) = S
            Q(I,NRPJ) = T
  104    CONTINUE
  105 CONTINUE
      DO 106 I=1,MR
         Q(I,NR) = .5D0*Q(I,NR)
         Q(I,N) = .5D0*Q(I,N)
  106 CONTINUE
      GO TO 118
  107 CONTINUE
C
C     ODD  NUMBER OF UNKNOWNS
C
      DO 109 J=1,NRM1
         NRPJ = N+1-J
         DO 108 I=1,MR
            S = Q(I,J)-Q(I,NRPJ)
            T = Q(I,J)+Q(I,NRPJ)
            Q(I,J) = S
            Q(I,NRPJ) = T
  108    CONTINUE
  109 CONTINUE
      DO 110 I=1,MR
         Q(I,NR) = 2.D0*Q(I,NR)
  110 CONTINUE
      LH = NRM1/2
      DO 112 J=1,LH
         NRMJ = NR-J
         DO 111 I=1,MR
            S = Q(I,J)
            Q(I,J) = Q(I,NRMJ)
            Q(I,NRMJ) = S
  111    CONTINUE
  112 CONTINUE
      CALL POISD2 (MR,NRM1,2,A,BB,C,Q,IDIMQ,B,W,D,TCOS,P)
      IPSTOR = W(1)
      CALL POISN2 (MR,NR,2,1,A,BB,C,Q(1,NR),IDIMQ,B,B2,B3,W,W2,W3,D,
     1             TCOS,P)
      IPSTOR = MAX0(IPSTOR,INT(W(1)))
      DO 114 J=1,NRM1
         NRPJ = NR+J
         DO 113 I=1,MR
            S = .5D0*(Q(I,NRPJ)+Q(I,J))
            T = .5D0*(Q(I,NRPJ)-Q(I,J))
            Q(I,NRPJ) = T
            Q(I,J) = S
  113    CONTINUE
  114 CONTINUE
      DO 115 I=1,MR
         Q(I,NR) = .5D0*Q(I,NR)
  115 CONTINUE
      DO 117 J=1,LH
         NRMJ = NR-J
         DO 116 I=1,MR
            S = Q(I,J)
            Q(I,J) = Q(I,NRMJ)
            Q(I,NRMJ) = S
  116    CONTINUE
  117 CONTINUE
  118 CONTINUE
C
C     RETURN STORAGE REQUIREMENTS FOR P VECTORS.
C
      W(1) = IPSTOR
      RETURN
      END
      SUBROUTINE STORE (X)
C     =================================================================
      IMPLICIT REAL*8(A-H,O-Z)
      COMMON /VALUE/  V
      V = X
      RETURN
      END
      SUBROUTINE TRI3 (M,A,B,C,K,Y1,Y2,Y3,TCOS,D,W1,W2,W3)
C     =================================================================
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION       A(*)       ,B(*)       ,C(*)       ,K(4)       ,
     1                TCOS(*)    ,Y1(*)      ,Y2(*)      ,Y3(*)      ,
     2                D(*)       ,W1(*)      ,W2(*)      ,W3(*)
C
C     SUBROUTINE TO SOLVE THREE LINEAR SYSTEMS WHOSE COMMON COEFFICIENT
C     MATRIX IS A RATIONAL FUNCTION IN THE MATRIX GIVEN BY
C
C                  TRIDIAGONAL (...,A(I),B(I),C(I),...)
C
      MM1 = M-1
      K1 = K(1)
      K2 = K(2)
      K3 = K(3)
      K4 = K(4)
      F1 = K1+1
      F2 = K2+1
      F3 = K3+1
      F4 = K4+1
      K2K3K4 = K2+K3+K4
      IF (K2K3K4 .EQ. 0) GO TO 101
      L1 = F1/F2
      L2 = F1/F3
      L3 = F1/F4
      LINT1 = 1
      LINT2 = 1
      LINT3 = 1
      KINT1 = K1
      KINT2 = KINT1+K2
      KINT3 = KINT2+K3
  101 CONTINUE
      DO 115 N=1,K1
         X = TCOS(N)
         IF (K2K3K4 .EQ. 0) GO TO 107
         IF (N .NE. L1) GO TO 103
         DO 102 I=1,M
            W1(I) = Y1(I)
  102    CONTINUE
  103    IF (N .NE. L2) GO TO 105
         DO 104 I=1,M
            W2(I) = Y2(I)
  104    CONTINUE
  105    IF (N .NE. L3) GO TO 107
         DO 106 I=1,M
            W3(I) = Y3(I)
  106    CONTINUE
  107    CONTINUE
         Z = 1.D0/(B(1)-X)
         D(1) = C(1)*Z
         Y1(1) = Y1(1)*Z
         Y2(1) = Y2(1)*Z
         Y3(1) = Y3(1)*Z
         DO 108 I=2,M
            Z = 1.D0/(B(I)-X-A(I)*D(I-1))
            D(I) = C(I)*Z
            Y1(I) = (Y1(I)-A(I)*Y1(I-1))*Z
            Y2(I) = (Y2(I)-A(I)*Y2(I-1))*Z
            Y3(I) = (Y3(I)-A(I)*Y3(I-1))*Z
  108    CONTINUE
         DO 109 IP=1,MM1
            I = M-IP
            Y1(I) = Y1(I)-D(I)*Y1(I+1)
            Y2(I) = Y2(I)-D(I)*Y2(I+1)
            Y3(I) = Y3(I)-D(I)*Y3(I+1)
  109    CONTINUE
         IF (K2K3K4 .EQ. 0) GO TO 115
         IF (N .NE. L1) GO TO 111
         I = LINT1+KINT1
         XX = X-TCOS(I)
         DO 110 I=1,M
            Y1(I) = XX*Y1(I)+W1(I)
  110    CONTINUE
         LINT1 = LINT1+1
         L1 = (DFLOAT(LINT1)*F1)/F2
  111    IF (N .NE. L2) GO TO 113
         I = LINT2+KINT2
         XX = X-TCOS(I)
         DO 112 I=1,M
            Y2(I) = XX*Y2(I)+W2(I)
  112    CONTINUE
         LINT2 = LINT2+1
         L2 = (DFLOAT(LINT2)*F1)/F3
  113    IF (N .NE. L3) GO TO 115
         I = LINT3+KINT3
         XX = X-TCOS(I)
         DO 114 I=1,M
            Y3(I) = XX*Y3(I)+W3(I)
  114    CONTINUE
         LINT3 = LINT3+1
         L3 = (DFLOAT(LINT3)*F1)/F4
  115 CONTINUE
      RETURN
      END
      SUBROUTINE TRIX (IDEGBR,IDEGCR,M,A,B,C,Y,TCOS,D,W)
C     =================================================================
C     SUBROUTINE TO SOLVE A SYSTEM OF LINEAR EQUATIONS WHERE THE
C     COEFFICIENT MATRIX IS A RATIONAL FUNCTION IN THE MATRIX GIVEN BY
C     TRIDIAGONAL  ( . . . , A(I), B(I), C(I), . . . ).
C
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION       A(*)       ,B(*)       ,C(*)       ,Y(*)       ,
     1                TCOS(*)    ,D(*)       ,W(*)
      MM1 = M-1
      FB = IDEGBR+1
      FC = IDEGCR+1
      L = FB/FC
      LINT = 1
      DO 108 K=1,IDEGBR
         X = TCOS(K)
         IF (K .NE. L) GO TO 102
         I = IDEGBR+LINT
         XX = X-TCOS(I)
         DO 101 I=1,M
            W(I) = Y(I)
            Y(I) = XX*Y(I)
  101    CONTINUE
  102    CONTINUE
         Z = 1.D0/(B(1)-X)
         D(1) = C(1)*Z
         Y(1) = Y(1)*Z
         DO 103 I=2,MM1
            Z = 1.D0/(B(I)-X-A(I)*D(I-1))
            D(I) = C(I)*Z
            Y(I) = (Y(I)-A(I)*Y(I-1))*Z
  103    CONTINUE
         Z = B(M)-X-A(M)*D(MM1)
         IF (Z .NE. 0.D0) GO TO 104
         Y(M) = 0.D0
         GO TO 105
  104    Y(M) = (Y(M)-A(M)*Y(MM1))/Z
  105    CONTINUE
         DO 106 IP=1,MM1
            I = M-IP
            Y(I) = Y(I)-D(I)*Y(I+1)
  106    CONTINUE
         IF (K .NE. L) GO TO 108
         DO 107 I=1,M
            Y(I) = Y(I)+W(I)
  107    CONTINUE
         LINT = LINT+1
         L = (DFLOAT(LINT)*FB)/FC
  108 CONTINUE
      RETURN
      END
c==================================================================
      subroutine uv2div 
     o		    (div   ,
     i		     u     , v     , nlo   , nla   )
c==================================================================
c     purpose
c          calculate divergence from (u,v)
c     *************************************************************
c     arguments
c     nlo:number of grid points in longitudinal direction(input)
c     nla:number of grid points in latitudinal direction(input)
c     u(nlo+1,nla):wind in longitudinal direction(input) 
c     v(nlo+1,nla):wind in latitudinal direction(input)
c     div(nlo+1,nla):divergence(output)
c     phi(nla):latitudes of the grid points(work array)
c     coa(nla):cosin of phi(j)(work array)
c     delx(nla):work array
c     dely(nla):work array
c==================================================================
      real       u   ( nlo+1, nla ), v   ( nlo+1, nla ),
     &           div ( nlo+1, nla )
      real       phi ( nla ), cosphi(nla), sinphi(nla),
     &           delx( nla ), dely( nla )
c===================================================================
c
c     calculate parameters
c     note !!!: from north to south
c
      call const 
     i	        (nlo   , nla   ,
     o           a     , pai   ,
     o	         phi   , cosphi, sinphi,
     o           dlamda, dphi  , delx  , dely , delyc)
c     
      do 3000 j=2,nla-1
      do 2000 i=2,nlo
         div(i,j)=(u(i+1,j)-u(i-1,j))/delx(j)-
     &            (v(i,j+1)*cosphi(j+1)-v(i,j-1)*cosphi(j-1))/dely(j)
 2000 continue
 3000 continue
c
c      < divergence in the first and end points >
c  
      do 3001 j=2,nla-1
         div(1,j)=(u(2,j)-u(nlo,j))/delx(j)-
     &            (v(1,j+1)*cosphi(j+1)-v(1,j-1)*cosphi(j-1))/dely(j)
         div(nlo+1,j) = div(1,j)
 3001 continue

c
c       < divergence in the north and south poles >
c


      do 2001 i=1,nlo+1
c         div(i,1)  =-(v(i,2)*cosphi(2)-v(i,1)*cosphi(1))/(dely(2)/2.0)
c         div(i,nla)=-(v(i,nla)*cosphi(nla)-v(i,nla-1)*cosphi(nla-1))/
c     &               (dely(nla-1)/2.0)
          div(i,1)=div(i,2)
          div(i,nla)=div(i,nla-1)
 2001 continue

      return
      end
c================================================================
      subroutine uv2psi
     o                (psi,vor,
     i                 u,v,nlo,nla,
     w                 f  ,coslat)
c================================================================
c     purpose
c          calculate streamfunction from (u,v)
c     *************************************************************
c     arguments
c     nlo:number of grid points in longitudinal direction(input)
c     nla:number of grid points in latitudinal direction(input)
c     u(nlo+1,nla):wind in longitudinal direction(input) 
c     v(nlo+1,nla):wind in latitudinal direction(input)
c     vor(nlo+1,nla):vorticity(output)
c     psi(nlo+1,nla):streamfunction(output)
c     phi(nla):latitudes of the grid points(work array)
c     f(nla,nlo+1):work array
c     coslat(nla):work array
c================================================================
      real       vor   ( nlo+1, nla ), psi   ( nlo+1, nla ),
     &           u     ( nlo+1, nla ), v     ( nlo+1, nla )
      real*8     f   ( nla, nlo+1 ), coslat ( nla )
c===================================================================
      er=6.37122e+6
c 
c          < compute vorticity >
c
      call uv2vor 
     o           (vor   ,
     i		  u     , v     ,  nlo   , nla   )
c
c         < compute streamfunction > 
c 
      call invlap 
     &           (nlo,nla,f,coslat,psi,vor)
c
      do 2000 j=1,nla
         do 1000 i=1,nlo+1
           psi(i,j) = +er*er*psi(i,j)
 1000    continue
 2000 continue
c
      return
      end
c===================================================================
      subroutine uv2vor 
     o		    (vor   ,
     i		     u     , v     ,    nlo   , nla   )
c===================================================================
c     purpose
c          calculate vorticity from (u,v)
c     *************************************************************
c     arguments
c     nlo:number of grid points in longitudinal direction(input)
c     nla:number of grid points in latitudinal direction(input)
c     u(nlo+1,nla):wind in longitudinal direction(input) 
c     v(nlo+1,nla):wind in latitudinal direction(input)
c     vor(nlo+1,nla):vorticity(output)
c     phi(nla):latitudes of the grid points(work array)
c     coa(nla):cosin of phi(j)(work array)
c     delx(nla):work array
c     dely(nla):work array
c===================================================================
      real       u   ( nlo+1, nla ), v   ( nlo+1, nla ),
     &           vor ( nlo+1, nla )
      real       phi ( nla ), cosphi(nla), sinphi(nla),
     &           delx( nla ), dely( nla )
c===================================================================
c
c     calculate parameters
c     note !!!: from north to south
c
      call const 
     i	        (nlo   , nla   ,
     o           a     , pai   ,
     o	         phi   , cosphi, sinphi,
     o           dlamda, dphi  , delx  , dely , delyc)
c     
      do 3000 j=2,nla-1
      do 2000 i=2,nlo
          vor(i,j)=(v(i+1,j)-v(i-1,j))/delx(j)+
     &             (u(i,j+1)*cosphi(j+1)-u(i,j-1)*cosphi(j-1))/dely(j)
 2000 continue
 3000 continue
c  
c     < vorticity in the first and end points >
c 
      do 3001 j=2,nla-1
         vor(1,j)=(v(2,j)-v(nlo,j))/delx(j)+
     &            (u(1,j+1)*cosphi(j+1)-u(1,j-1)*cosphi(j-1))/dely(j)
        vor(nlo+1,j) = vor(1,j)
 3001 continue
c
c      < vorticity in the north and south poles >
c
      do 2001 i=1,nlo+1
c         vor(i,1)  = (u(i,2)*cosphi(2)-u(i,1)*cosphi(1))/(dely(2)/2.0)
c         vor(i,nla)= (u(i,nla)*cosphi(nla)-u(i,nla-1)*cosphi(nla-1))/
c     &               (dely(nla-1)/2.0)
          vor(i,1)  =vor(i,2)
          vor(i,nla)=vor(i,nla-1)
 2001 continue
      return
      end

c================================================================
      subroutine uv2xai
     o                (xai,div,
     i                 u,v,nlo,nla,
     w                 f,coslat)
c================================================================
c     purpose
c          calculate velocity potential from (u,v)
c     *************************************************************
c     arguments
c     nlo:number of grid points in longitudinal direction(input)
c     nla:number of grid points in latitudinal direction(input)
c     u(nlo+1,nla):wind in longitudinal direction(input) 
c     v(nlo+1,nla):wind in latitudinal direction(input)
c     div(nlo+1,nla):vorticity(output)
c     xai(nlo+1,nla):Velocity potential(output)
c     f(nla,nlo+1):work array
c     coslat(nla):work array
c================================================================
      real    div   ( nlo+1, nla ), xai   ( nlo+1, nla ),
     &        u     ( nlo+1, nla ), v     ( nlo+1, nla )
      real*8  f     ( nla, nlo+1 ), coslat ( nla )
c===================================================================
      er=6.37122e+6
c
c          < compute vorticity >
c
      call uv2div 
     o           (div   ,
     i		  u     , v     ,  nlo   , nla   )
c
c         < compute velocity potential > 

c 
      call invlap 
     &           (nlo,nla,f,coslat,xai,div)
c

      do 2000 j=1,nla
         do 1000 i=1,nlo+1
           xai(i,j) = +er*er*xai(i,j)
 1000    continue
 2000 continue
c
      return
      end
