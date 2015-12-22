C gfortran -fpic -g -c agree2.f -o agree2.o; gcc -shared -o agree2.so agree2.o -lgfortran -lm -lquadmath -L/usr/lib/R/lib -lR

C  PROGRAM TO CALCULATE AGREEMENT STATISTICS AND  VARIANCES
C  ON INDIVIDUAL SUBJECTS AND OVER A GROUP OF SUBJECTS FOR
C  M(=>2) OBSERVERS
C  USING NOMINAL DATA AND 2 MEASURES OF DISAGREEMENT FOR
C  ORDINAL DATA
C  CASE 1 IS ASSUMING THAT OBSERVERS ARE DISTINCT
C  CASE 2 IS ASSUMING THAT THE OBSERVERS ARE HOMOGENEOUS
C  VARIABLES:
C    NRATER           NO. OF OBSERVERS PER DATA SET
C    NSCORE           NO. OF CATEGORIES - 1,2,....,NSCORE
C    NSUBJ            NO. OF SUBJECTS (ITEMS)
C    P1(J,K)          PROBABILITY OF OBSERVER J GIVING SCORE K
C                     WHEN OBSERVERS ARE DISTINCT
C    P2(K)            PROBABILITY OF SCORE K WHEN OBSERVERS ARE
C                     HOMOGENEOUS
C    W1(J,K)          WEIGHTED AVERAGE OF DS FOR OBSERVER J,
C                     SCORE K
C    W2(K)            WEIGHTED AVREAGE OF DS FOR SCORE K WHEN
C                     OBSERVERS ARE HOMOGENEOUS
C    D(J)             AMOUNT OF DISAGREEMENT FOR SUBJECT J
C    S1(J)            CHANCE-CORRECTED AGREEMENT STATISTIC FOR
C                     SUBJECT J WHEN OBSERVERS ARE DISTINCT
C    S2(J)            CHANCE-CORRECTED AGREEMENT STATISTIC FOR
C                     SUBJECT J WHEN OBSERVERS ARE HOMOGENEOUS
C                     S(J) = 1 - D(J)/EXPDEL
C    DELTA(J,K)       J<K AMOUNT OF DISAGREEMENT EXPECTED BY CHANCE
C                     FOR OBSERVERS J AND K
C                     J>K AMOUNT OF DISAGREEMENT EXPECTED BY CHANCE
C                     FOR OBSERVERS J AND K WHEN OBSERVERS ARE
C                     HOMOGENEOUS
C    EXPD1            AMOUNT OF DISAGREEMENT EXPECTED BY CHANCE
C                     IN NULL CASE WHEN OBSERVERS ARE DISTINCT
C    EXPD2            AMOUNT OF DISAGREEMENT EXPECTED BY CHANCE
C                     WHEN OBSERVERS ARE HOMOGENEOUS
C    DBAR             AVERAGE VALUE OF D OVER ALL SUBJECTS
C                     DBAR = (D(1)+D(2)+.....+D(N))/N
C    SAV1             CHANCE-CORRECTED AGREEMENT STATISTIC OVER
C                     ALL SUBJECTS WHEN OBSERVERS ARE DISTINCT
C    SAV2             CHANCE-CORRECTED AGREEMENT STATISTIC OVER
C                     ALL SUBJECTS WHEN OBSERVERS ARE
C                     HOMOGENEOUS
C                     SAV = 1 - DBAR/EXPDEL
C    VAR0S1           NULL VARIANCE OF S WHEN OBSERVERS ARE
C                     DISTINCT
C    VAR0S2           NULL VARIANCE OF S WHEN OBSERVERS ARE
C                     HOMOGENEOUS
C    VARS1            UNCONSTRAINED VARIANCE OF S WHEN
C                     OBSERVERS ARE DISTINCT
C    VARS2            UNCONSTRAINED VARIANCE OF S WHEN
C                     OBSERVERS ARE HOMOGENEOUS
C    V0SAV1           NULL VARIANCE OF SAV WHEN OBSERVERS ARE
C                     DISTINCT
C    V0SAV2           NULL VARIANCE OF SAV WHEN OBSERVERS ARE
C                     HOMOGENEOUS
C    VSAV1            UNCONSTRAINED VARIANCE OF SAV WHEN
C                     OBSERVERS ARE DISTINCT
C    VSAV2            UNCONSTRAINED VARIANCE OF SAV WHEN
C                     OBSERVERS ARE HOMOGENEOUS
C    RESP(I,J)        RESPONSE FOR OBSERVER I ON SUBJECT J
C    R(I)             SCORE GIVEN BY I'TH OBSERVER
C    SCORE(I)         SCORE ASSOCIATED WITH ITH CATEGORY

C  DIMENSIONS:
C    P1(NRATER,NSCORE)
C    P2(NSCORE)
C    W1(NRATER,NSCORE)
C    W2(NRATER,NSCORE)
C    D(NSUBJ)
C    S1(NSUBJ)
C    S2(NSUBJ)
C    RESP(NRATER,NSUBJ)
C    R(NRATER)
C    SCORE(NSCORE)
C    DELTA(NRATER,NRATER)

C  MAIN PROGRAM
      subroutine oconnell(nrater,nsubj,nscore,i,resp,s1,s2,p1,p2,w1,w2
     $     ,score,delta,d,expd1,expd2,dbar,vars1,var0s1,vsav1,v0sav1
     $     ,vars2,var0s2,vsav2,v0sav2,sav1,sav2)
      integer nrater,nsubj,nscore,i
      INTEGER R(nrater),RESP(nrater,nsubj)
      DOUBLE PRECISION S1(nsubj),S2(nsubj),P1(nrater,nscore),P2(nscore)
     $     ,W1(nrater,nscore)
      DOUBLE PRECISION W2(nrater,nscore),SCORE(nscore),delta(nrater
     $     ,nrater),D(nsubj)
      double precision EXPD1,EXPD2,DBAR
      double precision VARS1,VAR0S1,VSAV1,V0SAV1,VARS2,VAR0S2
      double precision VSAV2,V0SAV2,sav1,sav2
C  ZERO OBSERVERS' MARGINAL FREQUENCIES
      DO 60 K=1,NSCORE
      DO 70 J=1,NRATER
      P1(J,K) = 0.
70    CONTINUE
      P2(K) = 0.
60    CONTINUE
      DO 130 J=1,NSUBJ
      DO 140 K=1,NRATER
      R(K) = RESP(K,J)
      P1(K,R(K)) = P1(K,R(K)) + 1.
      P2(R(K)) = P2(R(K)) + 1.
140   CONTINUE
130   CONTINUE
      X = FLOAT(NSUBJ)
      Y = X*FLOAT(NRATER)
      DO 150 K=1,NSCORE
      DO 160 L=1,NRATER
      P1(L,K) = P1(L,K)/X
160   CONTINUE
      P2(K) = P2(K)/Y
 150  CONTINUE
      DO 540 J=1,NSUBJ
      DO 550 K=1,NRATER
      R(K) = RESP(K,J)
550   CONTINUE
c$$$      CALL CALD(D(J),R,nscore,nrater,i,score)
      D(J) = 0.
      K = 1
 551  L = K+1
 552  DIS = R(K)-R(L)
      A = SCORE(R(K)) - SCORE(R(L))
      IF(I.EQ.1) GO TO 553
      IF(I.EQ.2) GO TO 554
      d(j) = d(j) + A*A
      GO TO 555
 553  IF(A.NE.0) d(j)=d(j)+1.
      GO TO 555
 554  d(j) = d(j) + ABS(A)
 555  L = L+1
      IF(L.LE.NRATER) GO TO 552
      K = K+1
      IF(K.LT.NRATER) GO TO 551
540   CONTINUE
C  CALCULATE EXPECTED VALUE OF D AND THE FIRST TERMS OF VARS
C  AND VARSAV
c$$$      CALL EXPECT(nscore,nsubj,nrater,i,score,p1,p2,w1,w2,expd1,expd2
c$$$     $     ,delta)
      DO 1310 K=1,NSCORE
      DO 1320 J=1,NRATER
      W1(J,K) = 0.
      W2(J,K) = 0.
 1320 CONTINUE
 1310 CONTINUE
      EXPD1 = 0.
      EXPD2 = 0.
      DO 1330 J=1,NRATER
      DO 1335 K=1,NRATER
      DELTA(J,K) = 0.
 1335 CONTINUE
 1330 CONTINUE
      J = 1
 1340 JDASH = J+1
 1345 DO 1350 K=1,NSCORE
      Y = 0.
      Z = 0.
      DO 1360 KDASH=1,NSCORE
C  CALCULATE D FOR THIS PAIR OF SCORES
      A = SCORE(K) - SCORE(KDASH)
      DIS = K - KDASH
      IF(I.EQ.1) GO TO 1370
      IF(I.EQ.2) GO TO 1380
      B = A*A
      GO TO 1390
 1370 B = 0.
      IF(A.NE.0) B = 1
      GO TO 1390
 1380 B = ABS(A)
 1390 Y = Y + B*P1(JDASH,KDASH)
      Z = Z + B*P2(KDASH)
 1360 CONTINUE
      DELTA(J,JDASH) = DELTA(J,JDASH) + Y*P1(J,K)
      DELTA(JDASH,J) = DELTA(JDASH,J) + Z*P2(K)
 1350 CONTINUE
      EXPD1 = EXPD1 + DELTA(J,JDASH)
      EXPD2 = EXPD2 + DELTA(JDASH,J)
      IF(JDASH.EQ.NRATER) GO TO 1400
      JDASH = JDASH + 1
      GO TO 1345
 1400 IF(J.EQ.NRATER-1) GO TO 1410
      J = J+1
      GO TO 1340
 1410 DO 1420 J=1,NRATER
      DO 1430 K=1,NSCORE
      DO 1440 KDASH=1,NSCORE
      A = SCORE(K) - SCORE(KDASH)
      DIS = K - KDASH
      IF(I.EQ.1) GO TO 1450
      IF(I.EQ.2) GO TO 1460
      B = A*A
      GO TO 1470
 1450 B = 0.
      IF(A.NE.0) B=1
      GO TO 1470
 1460 B = ABS(A)
 1470 W1(J,K) = W1(J,K) + B*P1(J,KDASH)
      W2(J,K) = W2(J,K) + B*P2(KDASH)
 1440 CONTINUE
 1430 CONTINUE
 1420 CONTINUE
C  CALCULATE DBAR, S1(J) AND S2(J) FOR EACH OF THE SUBJECTS
C  (ITEMS)
      DBAR = 0.
      DO 560 J=1,NSUBJ
      DBAR = DBAR + D(J)
      S1(J) = 1. - D(J)/EXPD1
      S2(J) = 1. - D(J)/EXPD2
560   CONTINUE
      DBAR = DBAR/X
      SAV1 = 1. - DBAR/EXPD1
      SAV2 = 1. - DBAR/EXPD2
C  CALCULATE THE VARIANCES OF S(I) AND SAV
      SSAV1 = 0.
      SSAV2 = 0.
      SSI1 = 0.
      SSI2 = 0.
      S0SAV1 = 0.
      S0SI1 = 0.
      S0SAV2 = 0.
      S0SI2 = 0.
      S0SI21 = 0.
      S0SI22 = 0.
      DEL1 = 0.
      DEL2 = 0.
c     Is the next statement now necessary?
      X = FLOAT(NSUBJ)
      DO 1110 J=1,NSUBJ
      SDJI1 = 0.
      SDJI2 = 0.
      DO 1115 K=1,NRATER
      DO 1120 KDASH=1,NRATER
      IF(KDASH.EQ.K) GO TO 1120
      SDJI1 = SDJI1 + W1(K,RESP(KDASH,J))
      SDJI2 = SDJI2 + W2(K,RESP(KDASH,J))
 1120 CONTINUE
 1115 CONTINUE
      TERM = DBAR*SDJI1 - EXPD1*D(J)
      SSAV1 = SSAV1 + TERM*TERM
      TERM = DBAR*SDJI2 - EXPD2*D(J)
      SSAV2 = SSAV2 + TERM*TERM
      TERM = DBAR*SDJI1/X
      TERM2 = TERM - EXPD1*D(J)
      SSI1 = SSI1 + TERM2*TERM2 + (X-1)*TERM*TERM
      TERM = DBAR*SDJI2/X
      TERM2 = TERM - EXPD2*D(J)
      SSI2 = SSI2 + TERM2*TERM2 + (X-1)*TERM*TERM
 1110 CONTINUE
      DD = EXPD1*DBAR
      D2 = EXPD1*EXPD1
      D4 = D2*D2
      VSAV1 = SSAV1/X - DD*DD
      VSAV1 = VSAV1/(X*D4)
      VARS1 = SSI1/X - DD*DD
      VARS1 = VARS1/D4
      DD = EXPD2*DBAR
      D2 = EXPD2*EXPD2
      D4 = D2*D2
      VSAV2 = SSAV2/X - DD*DD
      VSAV2 = VSAV2/(X*D4)
      VARS2 = SSI2/X - DD*DD
      VARS2 = VARS2/D4
      J = 1
 1122 JDASH = J+1
 1124 DO 1130 K=1,NSCORE
      DO 1140 KDASH=1,NSCORE
C  CALCULATE D FOR THIS PAIR OF SCORES
      A = SCORE(K) - SCORE(KDASH)
      DIS = K - KDASH
      IF(I.EQ.1) GO TO 1150
      IF(I.EQ.2) GO TO 1160
      B = A*A
      GO TO 1170
 1150 B = 0.
      IF(A.NE.0) B=1
      GO TO 1170
 1160 B = ABS(A)
 1170 TERM = W1(J,KDASH) + W1(JDASH,K)
      TERM2 = TERM/X
      TEMP = TERM - B
      TEMP = TEMP*TEMP
      S0SAV1 = S0SAV1 + P1(J,K)*P1(JDASH,KDASH)*TEMP
      TEMP = TERM2 - B
      TEMP = TEMP*TEMP
      TEMP2 = (X-1)*TERM2*TERM2
      S0SI1 = S0SI1 + P1(J,K)*P1(JDASH,KDASH)*(TEMP+TEMP2)
      TERM = W2(J,KDASH) + W2(JDASH,K)
      TERM2 = TERM/X
      TEMP = TERM - B
      TEMP = TEMP*TEMP
      S0SAV2 = S0SAV2 + P2(K)*P2(KDASH)*TEMP
      TEMP = TERM2 - B
      TEMP = TEMP*TEMP
      TEMP2 = (X-1)*TERM2*TERM2
      S0SI2 = S0SI2 + P2(K)*P2(KDASH)*(TEMP+TEMP2)
 1140 CONTINUE
 1130 CONTINUE
      DO 1180 IA=1,NRATER
      IF(IA.EQ.J) GO TO 1180
      IF(IA.EQ.JDASH) GO TO 1180
      SUM1 = 0.
      SUM2 = 0.
      DO 1190 K=1,NSCORE
      TERM = P1(J,K)*W1(JDASH,K)
      TERM2 = P1(JDASH,K)*W1(J,K)
      SUM1 = SUM1 + W1(IA,K)*(TERM+TERM2)
      TERM = P2(K)*W2(JDASH,K)
      TERM2 = P2(K)*W2(J,K)
      SUM2 = SUM2 + W2(IA,K)*(TERM+TERM2)
 1190 CONTINUE
      DELAJ1 = DELTA(IA,J)
      DELAJD1 = DELTA(IA,JDASH)
      DELAJ2 = DELTA(J,IA)
      DELAJD2 = DELTA(JDASH,IA)
      IF(IA.GT.J) THEN
         DELAJ1 = DELTA(J,IA)
         DELAJ2 = DELTA(IA,J)
      END IF
      IF(IA.GT.JDASH) THEN
         DELAJD1 = DELTA(JDASH,IA)
         DELAJD2 = DELTA(IA,JDASH)
      END IF
      TERM = DELTA(J,JDASH)*(DELAJ1+DELAJD1)
      S0SI21 = S0SI21 + SUM1 - TERM
      TERM = DELTA(JDASH,J)*(DELAJ2+DELAJD2)
      S0SI22 = S0SI22 + SUM2 - TERM
 1180 CONTINUE
      DEL1 = DEL1 + DELTA(J,JDASH)*DELTA(J,JDASH)
      DEL2 = DEL2 + DELTA(JDASH,J)*DELTA(JDASH,J)
      IF(JDASH.EQ.NRATER) GO TO 1200
      JDASH = JDASH+1
      GO TO 1124
 1200 IF(J.EQ.NRATER-1) GO TO 1210
      J = J+1
      GO TO 1122
 1210 D2 = EXPD1*EXPD1
      V0SAV1 = S0SAV1 - DEL1
      V0SAV1 = V0SAV1/(X*D2)
      VAR0S1 = S0SI1 + (X-1)*S0SI21/X
      VAR0S1 = VAR0S1 - DEL1
      VAR0S1 = VAR0S1/D2
      D2 = EXPD2*EXPD2
      V0SAV2 = S0SAV2 - DEL2
      V0SAV2 = V0SAV2/(X*D2)
      VAR0S2 = S0SI2 + (X-1)*S0SI22/X
      VAR0S2 = VAR0S2 - DEL2
      VAR0S2 = VAR0S2/D2
      return
      END


