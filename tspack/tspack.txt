C
C                        TSPT1
C                       10/05/98
C
C   This program provides a quick easily-verified test of
C all TSPACK modules except the high-level interface
C routines TSPxx and TSVALx.  The primary test uses data
C values taken from the quadratic function f(x) = x**2 for
C which knot-derivative estimation and interpolation with no
C tension is exact.  Since this test function is positive,
C strictly increasing, and convex, zero tension is suffici-
C ent to preserve these properties in the data.  The
C abscissae are taken to be a set of N points uniformly
C distributed in [0,1]
C
C   The maximum absolute error in the derivative estimates
C is printed for each of the following.
C
C      YPC1 (local approximation of knot-derivatives)
C      YPC2 with end conditions from YPC1 estimates
C      YPC2 with specified endpoint first derivatives
C      YPC2 with specified endpoint second derivatives
C      YPC2 with end conditions computed by ENDSLP
C
C   The maximum error in a tension factor is printed for
C each of the following.  These should all be zero since no
C tension is required in any of the cases.
C
C      SIGBI (minimum tension factor required to satisfy
C             constraints defined by array B)
C      SIGS (minimum tension factor required to preserve
C            monotonicity and convexity on each interval)
C      SIG0 with a lower bound of Y(I) on (X(I),X(I+1))
C      SIG1 with a lower bound of YP(I) on (X(I),X(I+1))
C      SIG2 (minimum tension factor required to preserve
C            convexity on each interval)
C
C   The maximum absolute error on a set of 3*(N-1)+1 points
C uniformly distributed in [0,1] is printed for each of the
C following evaluation routines.
C
C      HVAL (function values)
C      HPVAL (first derivative values)
C      HPPVAL (second derivative values)
C      TSINTL (integral from 0 to T for each evaluation
C              point T)
C
C   The following four routines are used to construct para-
C metric tension spline fits to N points uniformly
C distributed on the unit circle.  In the case of SMCRV,
C periodic end conditions are used to fit cos(A) and natural
C end conditions are used for sin(A), where A is the parame-
C ter value (angle).  Thus, both B2TRI and B2TRIP are exer-
C cized.  The weights W are based on a standard deviation of
C EPS, for machine precision EPS, so that the smoothing
C curves are close to the interpolatory curves.  The maximum
C distance from the circle to a set of 3*(N-1)+1 evaluation
C points (angles uniformly distributed in [0,2*PI]) is
C printed for each routine.
C
C      YPC1P (local approximations to knot-derivatives with
C             periodic end conditions)
C      YPC2P (global approximations to knot-derivatives with
C             periodic end conditions)
C      SIGBP (minimum tension factor required to satisfy
C             bounds on distance between line segments and
C             planar curve segments)
C      SMCRV (global approximations to both function values
C             and first derivatives at the knots)
C
C   The final test consists of computing a set of parameter
C values associated with the N points on the unit circle.
C These are computed by both ARCL2D and ARCL3D with constant
C Z values, and the maximum difference is printed.
C
C      ARCL2D (cumulative arc lengths for a planar curve)
C      ARCL3D (cumulative arc lengths for a space curve)
C
      INTEGER I, IC, IER, IFL, IM1, IP1, ISL, K, LWK, N,
     .        NM1, NPTS
      PARAMETER (N=33, NM1=N-1, LWK=10*N)
      INTEGER ICFLG(NM1)
      LOGICAL PERIOD
      DOUBLE PRECISION A(N), X(N), XP(N), Y(N), YP(N),
     .                 XS(N), YS(N), SIGMA(NM1), W(N),
     .                 WK(LWK), YPTRUE(N), B(5,N),
     .                 BL(NM1), BU(NM1)
      DOUBLE PRECISION AI, BMAX, BND, BV1, BVN, DA, DSMAX,
     .                 DT, EPS, EPSP1, ERR, ERR0, ERR1,
     .                 ERR2, ERRI, ERRT, H, HI, HIT, HP,
     .                 HPP, HPPT, HPT, HT, SIG, SM,
     .                 SMTOL, T, TOL, TWOPI, XI, XT, YT
      DOUBLE PRECISION HPPVAL, HPVAL, HVAL, SIG0, SIG1,
     .                 SIG2, STORE, TSINTL
C
C Print a heading.
C
      WRITE (*,100)
  100 FORMAT (///,17X,'TSPT1 (TSPACK TEST PROGRAM)'///
     .        2X,'THE NAME OF EACH TESTED MODULE IS ',
     .        'FOLLOWED BY THE MAXIMUM'/
     .        2X,'ABSOLUTE ERROR.  THIS SHOULD BE A ',
     .        'SMALL MULTIPLE OF THE'/
     .        2X,'MACHINE PRECISION EXCEPT IN THE ',
     .        'CASE OF YPC1P, YPC2P AND SMCRV.'///)
C
C Compute the machine precision EPS.
C
      EPS = 1.D0
    1 EPS = EPS/2.D0
        EPSP1 = STORE(EPS+1.D0)
        IF (EPSP1 .GT. 1.D0) GO TO 1
      EPS = 2.D0*EPS
C
C Compute abscissae X uniformly distributed in [0,1],
C   ordinates Y from Y(I) = f(X(I)), for f(x) = x**2,
C   true derivative values YPTRUE from f'(x) = 2*x,
C   zero tension factors SIGMA, bounds B for SIGBI, and
C   uniform weights W for SMCRV.
C
      H = 1.D0/DBLE(NM1)
      DO 2 I = 1,N
        XI = DBLE(I-1)*H
        X(I) = XI
        Y(I) = XI*XI
        YPTRUE(I) = 2.D0*XI
        IF (I .LT. N) THEN
          SIGMA(I) = 0.D0
          B(1,I) = 1.D0
          B(2,I) = Y(I)
          B(3,I) = 2.D0
          B(4,I) = YPTRUE(I)
          B(5,I) = 1.D0
        ENDIF
        W(I) = 1.D0/(EPS*EPS)
    2   CONTINUE
C
C Test YPC1 (ISL=-1) and YPC2 with all four types of end
C   conditions.  The knot-derivative estimates should agree
C   with YPTRUE in all four cases.
C
      BV1 = 0.D0
      BVN = 2.D0
      DO 4 ISL = -1,3
        IF (ISL .LT. 0) THEN
C
C * YPC1 test:
C
          CALL YPC1 (N,X,Y, YP,IER)
        ELSE
          IF (ISL .EQ. 2) BV1 = 2.D0
C
C * YPC2 test:
C
          CALL YPC2 (N,X,Y,SIGMA,ISL,ISL,BV1,BVN,WK, YP,IER)
        ENDIF
        IF (IER .NE. 0) GO TO 31
        ERR = 0.D0
        DO 3 I = 1,N
          ERR = MAX( ERR, ABS(YP(I)-YPTRUE(I)) )
    3     CONTINUE
        IF (ISL .LT. 0) THEN
          WRITE (*,110) ERR
  110     FORMAT (15X,'YPC1',18X,D8.2)
        ELSE
          WRITE (*,120) ISL, ERR
  120     FORMAT (15X,'YPC2, ISL1=ISL2=',I1,5X,D8.2)
        ENDIF
    4   CONTINUE
C
C Test SIGBI, SIGS, SIG0, SIG1, and SIG2 with bounds
C   which are satisfied by the cubic (zero tension).
C
      TOL = 0.D0
C
C * SIGBI test:
C
      I = N
      BMAX = 3.D0
      CALL SIGBI (N,X,Y,YPTRUE,TOL,B,BMAX, SIGMA, ICFLG,
     .            DSMAX,IER)
      IF (IER .LT. 0) GO TO 32
      DO 20 I = 1,NM1
        IF (ICFLG(I) .NE. 0) GO TO 32
   20   CONTINUE
      WRITE (*,130) DSMAX
  130 FORMAT (15X,'SIGBI',17X,D8.2)
C
C * SIGS test:
C
      I = 0
      CALL SIGS (N,X,Y,YPTRUE,TOL, SIGMA, DSMAX,IER)
      IF (IER .LT. 0) GO TO 32
      WRITE (*,135) DSMAX
  135 FORMAT (15X,'SIGS',18X,D8.2)
C
C * SIG0 test:
C
      IFL = -1
      ERR = 0.D0
      DO 5 I = 1,NM1
        IP1 = I+1
        BND = B(2,I)
        SIG = SIG0 (X(I),X(IP1),Y(I),Y(IP1),YPTRUE(I),
     .              YPTRUE(IP1),IFL,BND,TOL, IER)
        IF (IER .NE. 0) GO TO 33
        ERR = MAX(ERR,SIG)
    5   CONTINUE
      WRITE (*,140) ERR
  140 FORMAT (15X,'SIG0',18X,D8.2)
C
C * SIG1 test:
C
      ERR = 0.D0
      DO 6 I = 1,NM1
        IP1 = I+1
        BND = B(4,I)
        SIG = SIG1 (X(I),X(IP1),Y(I),Y(IP1),YPTRUE(I),
     .              YPTRUE(IP1),IFL,BND,TOL, IER)
        IF (IER .NE. 0) GO TO 34
        ERR = MAX(ERR,SIG)
    6   CONTINUE
      WRITE (*,150) ERR
  150 FORMAT (15X,'SIG1',18X,D8.2)
C
C * SIG2 test:
C
      ERR = 0.D0
      IFL = 1
      DO 7 I = 1,NM1
        IP1 = I+1
        SIG = SIG2 (X(I),X(IP1),Y(I),Y(IP1),YPTRUE(I),
     .              YPTRUE(IP1),IFL,TOL, IER)
        IF (IER .NE. 0) GO TO 35
        ERR = MAX(ERR,SIG)
    7   CONTINUE
      WRITE (*,160) ERR
  160 FORMAT (15X,'SIG2',18X,D8.2)
C
C Test the evaluation routines HVAL, HPVAL, HPPVAL, and
C   TSINTL with SIGMA = 0.
C
      DO 8 I = 1,NM1
        SIGMA(I) = 0.D0
    8   CONTINUE
C
C   The number of evaluation points NPTS is taken to be
C     three per subinterval uniformly distributed in [0,1].
C
      ERR0 = 0.D0
      ERR1 = 0.D0
      ERR2 = 0.D0
      ERRI = 0.D0
      NPTS = 3*NM1 + 1
      DT = 1.D0/DBLE(NPTS-1)
      DO 9 K = 1,NPTS
        T = DBLE(K-1)*DT
C
C * HVAL test:
C
        H = HVAL (T,N,X,Y,YPTRUE,SIGMA, IER)
        IF (IER .LT. 0) GO TO 36
        HT = T*T
        ERR0 = MAX(ERR0,ABS(HT-H))
C
C * HPVAL test:
C
        HP = HPVAL (T,N,X,Y,YPTRUE,SIGMA, IER)
        IF (IER .LT. 0) GO TO 37
        HPT = 2.D0*T
        ERR1 = MAX(ERR1,ABS(HPT-HP))
C
C * HPPVAL test:
C
        HPP = HPPVAL (T,N,X,Y,YPTRUE,SIGMA, IER)
        IF (IER .LT. 0) GO TO 38
        HPPT = 2.D0
        ERR2 = MAX(ERR2,ABS(HPPT-HPP))
C
C * TSINTL test:
C
        HI = TSINTL (0.D0,T,N,X,Y,YPTRUE,SIGMA, IER)
        IF (IER .LT. 0) GO TO 39
        HIT = T**3/3.D0
        ERRI = MAX(ERRI,ABS(HIT-HI))
C
    9   CONTINUE
      WRITE (*,170) ERR0
  170 FORMAT (15X,'HVAL',18X,D8.2)
      WRITE (*,180) ERR1
  180 FORMAT (15X,'HPVAL',17X,D8.2)
      WRITE (*,190) ERR2
  190 FORMAT (15X,'HPPVAL',16X,D8.2)
      WRITE (*,200) ERRI
  200 FORMAT (15X,'TSINTL',16X,D8.2)
C
C Test YPC1P, YPC2P, SIGBP, and SMCRV (with both periodic
C   and natural end conditions) by approximating a circle
C   with parametric tension splines.
C
C   Compute N points on the unit circle along with upper and
C     lower bounds B for SIGBP.  Parameter values are angles
C     A uniformly distributed in [0,2*PI].  In each interval
C     the upper bound (distance toward the center) is EPS,
C     and the lower bound is the orthogonal distance from
C     the midpoint of the line segment to the circle.
C
      TWOPI = 8.D0*ATAN(1.D0)
      DA = TWOPI/DBLE(NM1)
      DO 10 I = 1,N
        AI = DBLE(I-1)*DA
        A(I) = AI
        X(I) = COS(AI)
        Y(I) = SIN(AI)
        IF (I .EQ. 1) GO TO 10
C
        IM1 = I - 1
        T = (A(IM1)+AI)/2.D0
        BU(IM1) = EPS
        BL(IM1) = -SQRT( (COS(T)-(X(IM1)+X(I))/2.D0)**2 +
     .                   (SIN(T)-(Y(IM1)+Y(I))/2.D0)**2 )
   10   CONTINUE
C
C * YPC1P, YPC2P, SIGBP test:
C
      DO 12 IC = 1,2
        IF (IC .EQ. 1) THEN
          CALL YPC1P (N,A,X, XP,IER)
          IF (IER .NE. 0) GO TO 40
          CALL YPC1P (N,A,Y, YP,IER)
          IF (IER .NE. 0) GO TO 40
        ELSE
          CALL YPC2P (N,A,X,SIGMA,WK, XP,IER)
          IF (IER .NE. 0) GO TO 40
          CALL YPC2P (N,A,Y,SIGMA,WK, YP,IER)
          IF (IER .NE. 0) GO TO 40
        ENDIF
C
C   Compute bounds-constrained tension factors SIGMA.
C     Note that, in the case of YPC2P, XP and YP are
C     not recomputed following the change in SIGMA,
C     and the splines are therefore not C-2.
C
        BMAX = 1.D0
        CALL SIGBP (N,X,Y,XP,YP,TOL,BL,BU,
     .              BMAX, SIGMA, DSMAX,IER)
        I = -1
        IF (IER .LT. 0) GO TO 32
C
        ERR = 0
        NPTS = 3*NM1 + 1
        DT = TWOPI/DBLE(NPTS-1)
        DO 11 K = 1,NPTS
          T = DBLE(K-1)*DT
          XT = HVAL (T,N,A,X,XP,SIGMA, IER)
          IF (IER .LT. 0) GO TO 36
          YT = HVAL (T,N,A,Y,YP,SIGMA, IER)
          IF (IER .LT. 0) GO TO 36
          ERRT = SQRT( (COS(T)-XT)**2 + (SIN(T)-YT)**2 )
          ERR = MAX(ERR,ERRT)
   11     CONTINUE
        IF (IC .EQ. 1) THEN
          WRITE (*,210) ERR
  210     FORMAT (15X,'YPC1P',17X,D8.2)
        ELSE
          WRITE (*,220) ERR
  220     FORMAT (15X,'YPC2P',17X,D8.2)
        ENDIF
   12   CONTINUE
C
C * SMCRV test:
C
      SM = DBLE(N)
      SMTOL = SQRT(2.D0/DBLE(N))
      PERIOD = .TRUE.
      CALL SMCRV (N,A,X,SIGMA,PERIOD,W,SM,SMTOL,WK, XS,
     .            XP,IER)
      IF (IER .NE. 0) GO TO 41
      PERIOD = .FALSE.
      CALL SMCRV (N,A,Y,SIGMA,PERIOD,W,SM,SMTOL,WK, YS,
     .            YP,IER)
      IF (IER .NE. 0) GO TO 41
      ERR = 0
      DO 13 K = 1,NPTS
        T = DBLE(K-1)*DT
        XT = HVAL (T,N,A,XS,XP,SIGMA, IER)
        IF (IER .LT. 0) GO TO 36
        YT = HVAL (T,N,A,YS,YP,SIGMA, IER)
        IF (IER .LT. 0) GO TO 36
        ERRT = SQRT( (COS(T)-XT)**2 + (SIN(T)-YT)**2 )
        ERR = MAX(ERR,ERRT)
   13   CONTINUE
      WRITE (*,230) ERR
  230 FORMAT (15X,'SMCRV',17X,D8.2)
C
C Test ARCL2D and ARCL3D by computing the differences
C   between cumulative arc lengths for the data sets
C   (X,Y) and (X,Y,W), where W(I) = 1 for all I.
C
      DO 14 I = 1,N
        W(I) = 1.D0
   14   CONTINUE
C
C * ARCL2D/ARCL3D test:
C
      CALL ARCL2D (N,X,Y, A,IER)
      IF (IER .NE. 0) GO TO 42
      CALL ARCL3D (N,X,Y,W, WK,IER)
      IF (IER .NE. 0) GO TO 43
      ERR = 0.D0
      DO 15 I = 1,N
        ERR = MAX( ERR, ABS(A(I)-WK(I)) )
   15   CONTINUE
      WRITE (*,250) ERR
  250 FORMAT (15X,'ARCL2D/ARCL3D',9X,D8.2)
      STOP
C
C Error in YPC1 or YPC2.
C
   31 IF (ISL .LT. 0) THEN
        WRITE (*,310) IER
  310   FORMAT (///10X,'*** ERROR IN YPC1.  IER = ',I2)
      ELSE
        WRITE (*,315) ISL, IER
  315   FORMAT (///10X,'*** ERROR IN YPC2 WITH ISL1 = ',
     .          'ISLN = ',I1,'.  IER = ',I2)
      ENDIF
      STOP
C
C Error in SIGBP (I=-1), SIGBI (I>0), or SIGS (I=0).
C
   32 IF (I .EQ. N) THEN
        WRITE (*,320) IER
  320   FORMAT (///10X,'*** ERROR IN SIGBI.  IER = ',I3)
      ELSEIF (I .GT. 0) THEN
        WRITE (*,322) ICFLG(I), I
  322   FORMAT (///10X,'*** ERROR IN SIGBI.  CONSTRAINT ',I1,
     .          ' IN INTERVAL ',I2,' IS INVALID.')
      ELSEIF (I .EQ. 0) THEN
        WRITE (*,324) IER
  324   FORMAT (///10X,'*** ERROR IN SIGS.  IER = ',I3)
      ELSE
        WRITE (*,326) IER
  326   FORMAT (///10X,'*** ERROR IN SIGBP.  IER = ',I3)
      ENDIF
      STOP
C
C Error in SIG0.
C
   33 WRITE (*,330) I, IER
  330 FORMAT (///10X,'*** ERROR IN SIG0, INTERVAL ',I2,
     .        '.  IER = ',I2)
      STOP
C
C Error in SIG1.
C
   34 WRITE (*,340) I, IER
  340 FORMAT (///10X,'*** ERROR IN SIG1, INTERVAL ',I2,
     .        '.  IER = ',I2)
      STOP
C
C Error in SIG2.
C
   35 WRITE (*,350) I, IER
  350 FORMAT (///10X,'*** ERROR IN SIG2, INTERVAL ',I2,
     .        '.  IER = ',I2)
      STOP
C
C Error in HVAL.
C
   36 WRITE (*,360) T, IER
  360 FORMAT (///10X,'*** ERROR IN HVAL AT T = ',D9.3,
     .        '.  IER = ',I2)
      STOP
C
C Error in HPVAL.
C
   37 WRITE (*,370) T, IER
  370 FORMAT (///10X,'*** ERROR IN HPVAL AT T = ',D9.3,
     .        '.  IER = ',I2)
      STOP
C
C Error in HPPVAL.
C
   38 WRITE (*,380) T, IER
  380 FORMAT (///10X,'*** ERROR IN HPPVAL AT T = ',D9.3,
     .        '.  IER = ',I2)
      STOP
C
C Error in TSINTL.
C
   39 WRITE (*,390) T, IER
  390 FORMAT (///10X,'*** ERROR IN TSINTL (0,T), T = ',D9.3,
     .        '.  IER = ',I2)
      STOP
C
C Error in YPC1P or YPC2P.
C
   40 IF (IC .EQ. 1) THEN
        WRITE (*,400) IER
  400   FORMAT (///10X,'*** ERROR IN YPC1P.  IER = ',I2)
      ELSE
        WRITE (*,405) IER
  405   FORMAT (///10X,'*** ERROR IN YPC2P.  IER = ',I2)
      ENDIF
      STOP
C
C Error in SMCRV.
C
   41 WRITE (*,410) IER
  410 FORMAT (///10X,'*** ERROR IN SMCRV.  IER = ',I3)
      STOP
C
C Error in ARCL2D.
C
   42 WRITE (*,420) IER
  420 FORMAT (///10X,'*** ERROR IN ARCL2D.  IER = ',I2)
      STOP
C
C Error in ARCL3D.
C
   43 WRITE (*,430) IER
  430 FORMAT (///10X,'*** ERROR IN ARCL3D.  IER = ',I2)
      STOP
      END
      SUBROUTINE ARCL2D (N,X,Y, T,IER)
      INTEGER N, IER
      DOUBLE PRECISION X(N), Y(N), T(N)
C
C***********************************************************
C
C                                                From TSPACK
C                                            Robert J. Renka
C                                  Dept. of Computer Science
C                                       Univ. of North Texas
C                                           renka@cs.unt.edu
C                                                   06/10/92
C
C   Given an ordered sequence of points (X,Y) defining a
C polygonal curve in the plane, this subroutine computes the
C sequence T of cumulative arc lengths along the curve:
C T(1) = 0 and, for 2 .LE. K .LE. N, T(K) is the sum of
C Euclidean distances between (X(I-1),Y(I-1)) and (X(I),Y(I))
C for I = 2,...,K.  A closed curve corresponds to X(1) =
C X(N) and Y(1) = Y(N), and more generally, duplicate points
C are permitted but must not be adjacent.  Thus, T contains
C a strictly increasing sequence of values which may be used
C as parameters for fitting a smooth curve to the sequence
C of points.
C
C On input:
C
C       N = Number of points defining the curve.  N .GE. 2.
C
C       X,Y = Arrays of length N containing the coordinates
C             of the points.
C
C The above parameters are not altered by this routine.
C
C       T = Array of length at least N.
C
C On output:
C
C       T = Array containing cumulative arc lengths defined
C           above unless IER > 0.
C
C       IER = Error indicator:
C             IER = 0 if no errors were encountered.
C             IER = 1 if N < 2.
C             IER = I if X(I-1) = X(I) and Y(I-1) = Y(I) for
C                     some I in the range 2,...,N.
C
C Modules required by ARCL2D:  None
C
C Intrinsic function called by ARCL2D:  SQRT
C
C***********************************************************
C
      INTEGER I, NN
      DOUBLE PRECISION DS
C
      NN = N
      IF (NN .LT. 2) GO TO 2
      T(1) = 0.D0
      DO 1 I = 2,NN
        DS = (X(I)-X(I-1))**2 + (Y(I)-Y(I-1))**2
        IF (DS .EQ. 0.D0) GO TO 3
        T(I) = T(I-1) + SQRT(DS)
    1   CONTINUE
      IER = 0
      RETURN
C
C N is outside its valid range.
C
    2 IER = 1
      RETURN
C
C Points I-1 and I coincide.
C
    3 IER = I
      RETURN
      END
      SUBROUTINE ARCL3D (N,X,Y,Z, T,IER)
      INTEGER N, IER
      DOUBLE PRECISION X(N), Y(N), Z(N), T(N)
C
C***********************************************************
C
C                                                From TSPACK
C                                            Robert J. Renka
C                                  Dept. of Computer Science
C                                       Univ. of North Texas
C                                           renka@cs.unt.edu
C                                                   06/10/92
C
C   Given an ordered sequence of points (X,Y,Z) defining a
C polygonal curve in 3-space, this subroutine computes the
C sequence T of cumulative arc lengths along the curve:
C T(1) = 0 and, for 2 .LE. K .LE. N, T(K) is the sum of
C Euclidean distances between (X(I-1),Y(I-1),Z(I-1)) and
C (X(I),Y(I),Z(I)) for I = 2,...,K.  A closed curve corre-
C sponds to X(1) = X(N), Y(1) = Y(N), and Z(1) = Z(N).  More
C generally, duplicate points are permitted but must not be
C adjacent.  Thus, T contains a strictly increasing sequence
C of values which may be used as parameters for fitting a
C smooth curve to the sequence of points.
C
C On input:
C
C       N = Number of points defining the curve.  N .GE. 2.
C
C       X,Y,Z = Arrays of length N containing the coordi-
C               nates of the points.
C
C The above parameters are not altered by this routine.
C
C       T = Array of length at least N.
C
C On output:
C
C       T = Array containing cumulative arc lengths defined
C           above unless IER > 0.
C
C       IER = Error indicator:
C             IER = 0 if no errors were encountered.
C             IER = 1 if N < 2.
C             IER = I if X(I-1) = X(I), Y(I-1) = Y(I), and
C                     Z(I-1) = Z(I) for some I in the range
C                     2,...,N.
C
C Modules required by ARCL3D:  None
C
C Intrinsic function called by ARCL3D:  SQRT
C
C***********************************************************
C
      INTEGER I, NN
      DOUBLE PRECISION DS
C
      NN = N
      IF (NN .LT. 2) GO TO 2
      T(1) = 0.D0
      DO 1 I = 2,NN
        DS = (X(I)-X(I-1))**2 + (Y(I)-Y(I-1))**2 +
     .       (Z(I)-Z(I-1))**2
        IF (DS .EQ. 0.D0) GO TO 3
        T(I) = T(I-1) + SQRT(DS)
    1 CONTINUE
      IER = 0
      RETURN
C
C N is outside its valid range.
C
    2 IER = 1
      RETURN
C
C Points I-1 and I coincide.
C
    3 IER = I
      RETURN
      END
      SUBROUTINE B2TRI (N,X,Y,W,P,D,SD,T11,T12,T21,T22, YS,
     .                  YP,IER)
      INTEGER N, IER
      DOUBLE PRECISION X(N), Y(N), W(N), P, D(N), SD(N),
     .                 T11(N), T12(N), T21(N), T22(N),
     .                 YS(N), YP(N)
C
C***********************************************************
C
C                                                From TSPACK
C                                            Robert J. Renka
C                                  Dept. of Computer Science
C                                       Univ. of North Texas
C                                           renka@cs.unt.edu
C                                                   06/10/92
C
C   This subroutine solves the order 2N symmetric positive-
C definite block tridiagonal linear system associated with
C minimizing the quadratic functional Q(YS,YP) described in
C Subroutine SMCRV.
C
C On input:
C
C       N = Number of data points.  N .GE. 2.
C
C       X,Y,W = Arrays of length N containing abscissae,
C               data values, and positive weights, respect-
C               ively.  The abscissae must be strictly in-
C               creasing.
C
C       P = Positive smoothing parameter defining Q.
C
C       D,SD = Arrays of length N-1 containing positive ma-
C              trix entries.  Letting DX and SIG denote the
C              width and tension factor associated with the
C              interval (X(I),X(I+1)), D(I) = SIG*(SIG*
C              COSHM(SIG) - SINHM(SIG))/(DX*E) and SD(I) =
C              SIG*SINHM(SIG)/(DX*E) where E = SIG*SINH(SIG)
C              - 2*COSHM(SIG).
C
C The above parameters are not altered by this routine.
C
C       T11,T12,T21,T22 = Arrays of length N-1 used as
C                         temporary work space.
C
C On output:
C
C       YS,YP = Arrays of length N containing solution com-
C               ponents:  function and derivative values,
C               respectively, at the abscissae.
C
C       IER = Error indicator:
C             IER = 0 if no errors were encountered.
C             IER = 1 if N or P is outside its valid range
C                     on input.
C             Note that no test is made for a nonpositive
C             value of X(I+1)-X(I), W(I), D(I), or SD(I).
C
C Modules required by B2TRI:  None
C
C***********************************************************
C
      INTEGER I, IM1, NM1, NN
      DOUBLE PRECISION D11I, D12I, D22I, DEN, DI, DIM1, DX,
     .                 PP, R1, R2, S11I, S11IM1, S12I,
     .                 S12IM1, S22I, S22IM1
C
      NN = N
      NM1 = NN - 1
      PP = P
      IER = 1
      IF (NN .LT. 2  .OR.  PP .LE. 0.D0) RETURN
C
C The forward elimination step consists of scaling a row by
C   the inverse of its diagonal block and eliminating the
C   subdiagonal block.  The superdiagonal is stored in T and
C   the right hand side in YS,YP.  For J = 11, 12, and 22,
C   SJI and SJIM1 denote the elements in position J of the
C   superdiagonal block in rows I and I-1, respectively.
C   Similarly, DJI denotes an element in the diagonal block
C   of row I.
C
C Initialize for I = 2.
C
      DX = X(2) - X(1)
      DIM1 = D(1)
      S22IM1 = SD(1)
      S12IM1 = (DIM1 + S22IM1)/DX
      S11IM1 = -2.D0*S12IM1/DX
      R1 = PP*W(1)
      D11I = R1 - S11IM1
      D12I = S12IM1
      D22I = DIM1
      DEN = D11I*D22I - D12I*D12I
      T11(1) = (D22I*S11IM1 + D12I*S12IM1)/DEN
      T12(1) = (D22I*S12IM1 - D12I*S22IM1)/DEN
      T21(1) = -(D12I*S11IM1 + D11I*S12IM1)/DEN
      T22(1) = (D11I*S22IM1 - D12I*S12IM1)/DEN
      R1 = R1*Y(1)/DEN
      YS(1) = D22I*R1
      YP(1) = -D12I*R1
C
C I = 2,...,N-1:
C
      DO 1 I = 2,NM1
        IM1 = I - 1
        DX = X(I+1) - X(I)
        DI = D(I)
        S22I = SD(I)
        S12I = (DI + S22I)/DX
        S11I = -2.D0*S12I/DX
        R1 = PP*W(I)
        D11I = R1 - S11IM1 - S11I - (S11IM1*T11(IM1) -
     .         S12IM1*T21(IM1))
        D12I = S12I - S12IM1 - (S11IM1*T12(IM1) - S12IM1*
     .         T22(IM1))
        D22I = DIM1 + DI - (S12IM1*T12(IM1)+S22IM1*T22(IM1))
        DEN = D11I*D22I - D12I*D12I
        T11(I) = (D22I*S11I + D12I*S12I)/DEN
        T12(I) = (D22I*S12I - D12I*S22I)/DEN
        T21(I) = -(D12I*S11I + D11I*S12I)/DEN
        T22(I) = (D11I*S22I - D12I*S12I)/DEN
        R1 = R1*Y(I) - S11IM1*YS(IM1) + S12IM1*YP(IM1)
        R2 = -S12IM1*YS(IM1) - S22IM1*YP(IM1)
        YS(I) = (D22I*R1 - D12I*R2)/DEN
        YP(I) = (D11I*R2 - D12I*R1)/DEN
        DIM1 = DI
        S22IM1 = S22I
        S12IM1 = S12I
        S11IM1 = S11I
    1   CONTINUE
C
C I = N:
C
      R1 = PP*W(NN)
      D11I = R1 - S11IM1 - (S11IM1*T11(NM1)-S12IM1*T21(NM1))
      D12I = -S12IM1 - (S11IM1*T12(NM1) - S12IM1*T22(NM1))
      D22I = DIM1 - (S12IM1*T12(NM1) + S22IM1*T22(NM1))
      DEN = D11I*D22I - D12I*D12I
      R1 = R1*Y(NN) - S11IM1*YS(NM1) + S12IM1*YP(NM1)
      R2 = -S12IM1*YS(NM1) - S22IM1*YP(NM1)
      YS(NN) = (D22I*R1 - D12I*R2)/DEN
      YP(NN) = (D11I*R2 - D12I*R1)/DEN
C
C Back solve the system.
C
      DO 2 I = NM1,1,-1
        YS(I) = YS(I) - (T11(I)*YS(I+1) + T12(I)*YP(I+1))
        YP(I) = YP(I) - (T21(I)*YS(I+1) + T22(I)*YP(I+1))
    2   CONTINUE
      IER = 0
      RETURN
      END
      SUBROUTINE B2TRIP (N,X,Y,W,P,D,SD,T11,T12,T21,T22,U11,
     .                   U12,U21,U22, YS,YP,IER)
      INTEGER N, IER
      DOUBLE PRECISION X(N), Y(N), W(N), P, D(N), SD(N),
     .                 T11(N), T12(N), T21(N), T22(N),
     .                 U11(N), U12(N), U21(N), U22(N),
     .                 YS(N), YP(N)
C
C***********************************************************
C
C                                                From TSPACK
C                                            Robert J. Renka
C                                  Dept. of Computer Science
C                                       Univ. of North Texas
C                                           renka@cs.unt.edu
C                                                   06/10/92
C
C   This subroutine solves the order 2(N-1) symmetric posi-
C tive-definite linear system associated with minimizing the
C quadratic functional Q(YS,YP) (described in Subroutine
C SMCRV) with periodic end conditions.  The matrix is block
C tridiagonal except for nonzero blocks in the upper right
C and lower left corners.
C
C On input:
C
C       N = Number of data points.  N .GE. 3.
C
C       X = Array of length N containing a strictly in-
C           creasing sequence of abscissae.
C
C       Y,W = Arrays of length N-1 containing data values
C             and positive weights, respectively, associated
C             with the first N-1 abscissae.
C
C       P = Positive smoothing parameter defining Q.
C
C       D,SD = Arrays of length N-1 containing positive ma-
C              trix elements.  Letting DX and SIG denote the
C              width and tension factor associated with the
C              interval (X(I),X(I+1)), D(I) = SIG*(SIG*
C              COSHM(SIG) - SINHM(SIG))/(DX*E) and SD(I) =
C              SIG*SINHM(SIG)/(DX*E) where E = SIG*SINH(SIG)
C              - 2*COSHM(SIG).
C
C The above parameters are not altered by this routine.
C
C       T11,T12,T21,T22,U11,U12,U21,U22 = Arrays of length
C                                         N-2 used as temp-
C                                         orary work space.
C
C On output:
C
C       YS,YP = Arrays of length N containing solution com-
C               ponents:  function and derivative values,
C               respectively, at the abscissae.  YS(N) =
C               YS(1) and YP(N) = YP(1).
C
C       IER = Error indicator:
C             IER = 0 if no errors were encountered.
C             IER = 1 if N or P is outside its valid range
C                     on input.
C             Note that no test is made for a nonpositive
C             value of X(I+1)-X(I), W(I), D(I), or SD(I).
C
C Modules required by B2TRIP:  None
C
C***********************************************************
C
      INTEGER I, IM1, IP1, NM1, NM2, NM3, NN
      DOUBLE PRECISION D11I, D12I, D22I, DEN, DI, DIM1,
     .                 DNM1, DX, PP, R1, R2, S11I, S11IM1,
     .                 S11NM1, S12I, S12IM1, S12NM1, S22I,
     .                 S22IM1, S22NM1, SU11, SU12, SU21,
     .                 SU22, YPNM1, YSNM1
C
      NN = N
      NM1 = NN - 1
      NM2 = NN - 2
      NM3 = NN - 3
      PP = P
      IER = 1
      IF (NN .LT. 3  .OR.  PP .LE. 0.D0) RETURN
C
C The forward elimination step consists of scaling a row by
C   the inverse of its diagonal block and eliminating the
C   subdiagonal block for the first N-2 rows.  The super-
C   diagonal is stored in T, the negative of the last column
C   in U, and the right hand side in YS,YP.  For J = 11, 12,
C   and 22, SJI and SJIM1 denote the elements in position J
C   of the superdiagonal block in rows I and I-1, respect-
C   ively.  Similarly, DJI denotes an element in the diago-
C   nal block of row I.
C
C I = 1:
C
      DX = X(NN) - X(NM1)
      DNM1 = D(NM1)
      S22NM1 = SD(NM1)
      S12NM1 = -(DNM1 + S22NM1)/DX
      S11NM1 = 2.D0*S12NM1/DX
      DX = X(2) - X(1)
      DI = D(1)
      S22I = SD(1)
      S12I = (DI + S22I)/DX
      S11I = -2.D0*S12I/DX
      R1 = PP*W(1)
      D11I = R1 - S11NM1 - S11I
      D12I = S12I + S12NM1
      D22I = DNM1 + DI
      DEN = D11I*D22I - D12I*D12I
      T11(1) = (D22I*S11I + D12I*S12I)/DEN
      T12(1) = (D22I*S12I - D12I*S22I)/DEN
      T21(1) = -(D12I*S11I + D11I*S12I)/DEN
      T22(1) = (D11I*S22I - D12I*S12I)/DEN
      U11(1) = -(D22I*S11NM1 + D12I*S12NM1)/DEN
      U12(1) = (D12I*S22NM1 - D22I*S12NM1)/DEN
      U21(1) = (D12I*S11NM1 + D11I*S12NM1)/DEN
      U22(1) = (D12I*S12NM1 - D11I*S22NM1)/DEN
      R1 = R1*Y(1)/DEN
      YS(1) = D22I*R1
      YP(1) = -D12I*R1
      IF (NN .EQ. 3) GO TO 2
C
C I = 2,...,N-2:
C
      DO 1 I = 2,NM2
        IM1 = I - 1
        DIM1 = DI
        S22IM1 = S22I
        S12IM1 = S12I
        S11IM1 = S11I
        DX = X(I+1) - X(I)
        DI = D(I)
        S22I = SD(I)
        S12I = (DI + S22I)/DX
        S11I = -2.D0*S12I/DX
        R1 = PP*W(I)
        D11I = R1 - S11IM1 - S11I - (S11IM1*T11(IM1) -
     .         S12IM1*T21(IM1))
        D12I = S12I - S12IM1 - (S11IM1*T12(IM1) - S12IM1*
     .         T22(IM1))
        D22I = DIM1 + DI - (S12IM1*T12(IM1)+S22IM1*T22(IM1))
        DEN = D11I*D22I - D12I*D12I
        T11(I) = (D22I*S11I + D12I*S12I)/DEN
        T12(I) = (D22I*S12I - D12I*S22I)/DEN
        T21(I) = -(D12I*S11I + D11I*S12I)/DEN
        T22(I) = (D11I*S22I - D12I*S12I)/DEN
        SU11 = S11IM1*U11(IM1) - S12IM1*U21(IM1)
        SU12 = S11IM1*U12(IM1) - S12IM1*U22(IM1)
        SU21 = S12IM1*U11(IM1) + S22IM1*U21(IM1)
        SU22 = S12IM1*U12(IM1) + S22IM1*U22(IM1)
        U11(I) = (D12I*SU21 - D22I*SU11)/DEN
        U12(I) = (D12I*SU22 - D22I*SU12)/DEN
        U21(I) = (D12I*SU11 - D11I*SU21)/DEN
        U22(I) = (D12I*SU12 - D11I*SU22)/DEN
        R1 = R1*Y(I) - S11IM1*YS(IM1) + S12IM1*YP(IM1)
        R2 = -S12IM1*YS(IM1) - S22IM1*YP(IM1)
        YS(I) = (D22I*R1 - D12I*R2)/DEN
        YP(I) = (D11I*R2 - D12I*R1)/DEN
    1   CONTINUE
C
C The backward elimination step zeros the first N-3 blocks
C   of the superdiagonal.  For I = N-2,N-3,...,1, T(I) and
C   (YS(I),YP(I)) are overwritten with the negative of the
C   last column and the new right hand side, respectively.
C
    2 T11(NM2) = U11(NM2) - T11(NM2)
      T12(NM2) = U12(NM2) - T12(NM2)
      T21(NM2) = U21(NM2) - T21(NM2)
      T22(NM2) = U22(NM2) - T22(NM2)
      DO 3 I = NM3,1,-1
        IP1 = I + 1
        YS(I) = YS(I) - T11(I)*YS(IP1) - T12(I)*YP(IP1)
        YP(I) = YP(I) - T21(I)*YS(IP1) - T22(I)*YP(IP1)
        T11(I) = U11(I) - T11(I)*T11(IP1) - T12(I)*T21(IP1)
        T12(I) = U12(I) - T11(I)*T12(IP1) - T12(I)*T22(IP1)
        T21(I) = U21(I) - T21(I)*T11(IP1) - T22(I)*T21(IP1)
        T22(I) = U22(I) - T21(I)*T12(IP1) - T22(I)*T22(IP1)
    3   CONTINUE
C
C Solve the last equation for YS(N-1),YP(N-1).  SJI = SJNM2
C   and DJI = DJNM1.
C
      R1 = PP*W(NM1)
      D11I = R1 - S11I - S11NM1 + S11NM1*T11(1) -
     .       S12NM1*T21(1) + S11I*T11(NM2) - S12I*T21(NM2)
      D12I = -S12NM1 - S12I + S11NM1*T12(1) - S12NM1*T22(1)
     .       + S11I*T12(NM2) - S12I*T22(NM2)
      D22I = DI + DNM1 + S12NM1*T12(1) + S22NM1*T22(1) +
     .       S12I*T12(NM2) + S22I*T22(NM2)
      DEN = D11I*D22I - D12I*D12I
      R1 = R1*Y(NM1) - S11NM1*YS(1) + S12NM1*YP(1) -
     .     S11I*YS(NM2) + S12I*YP(NM2)
      R2 = -S12NM1*YS(1) - S22NM1*YP(1) - S12I*YS(NM2) -
     .     S22I*YP(NM2)
      YSNM1 = (D22I*R1 - D12I*R2)/DEN
      YPNM1 = (D11I*R2 - D12I*R1)/DEN
      YS(NM1) = YSNM1
      YP(NM1) = YPNM1
C
C Back substitute for the remainder of the solution
C   components.
C
      DO 4 I = 1,NM2
        YS(I) = YS(I) + T11(I)*YSNM1 + T12(I)*YPNM1
        YP(I) = YP(I) + T21(I)*YSNM1 + T22(I)*YPNM1
    4   CONTINUE
C
C YS(N) = YS(1) and YP(N) = YP(1).
C
      YS(NN) = YS(1)
      YP(NN) = YP(1)
      IER = 0
      RETURN
      END
      DOUBLE PRECISION FUNCTION ENDSLP (X1,X2,X3,Y1,Y2,Y3,
     .                                  SIGMA)
      DOUBLE PRECISION X1, X2, X3, Y1, Y2, Y3, SIGMA
C
C***********************************************************
C
C                                                From TSPACK
C                                            Robert J. Renka
C                                  Dept. of Computer Science
C                                       Univ. of North Texas
C                                           renka@cs.unt.edu
C                                                   11/17/96
C
C   Given data values associated with a strictly increasing
C or decreasing sequence of three abscissae X1, X2, and X3,
C this function returns a derivative estimate at X1 based on
C the tension spline H(x) which interpolates the data points
C and has third derivative equal to zero at X1.  Letting S1
C denote the slope defined by the first two points, the est-
C mate is obtained by constraining the derivative of H at X1
C so that it has the same sign as S1 and its magnitude is
C at most 3*abs(S1).  If SIGMA = 0, H(x) is quadratic and
C the derivative estimate is identical to the value computed
C by Subroutine YPC1 at the first point (or the last point
C if the abscissae are decreasing).
C
C On input:
C
C       X1,X2,X3 = Abscissae satisfying either X1 < X2 < X3
C                  or X1 > X2 > X3.
C
C       Y1,Y2,Y3 = Data values associated with the abscis-
C                  sae.  H(X1) = Y1, H(X2) = Y2, and H(X3)
C                  = Y3.
C
C       SIGMA = Tension factor associated with H in inter-
C               val (X1,X2) or (X2,X1).
C
C Input parameters are not altered by this function.
C
C On output:
C
C       ENDSLP = (Constrained) derivative of H at X1, or
C                zero if the abscissae are not strictly
C                monotonic.
C
C Module required by ENDSLP:  SNHCSH
C
C Intrinsic functions called by ENDSLP:  ABS, EXP, MAX, MIN
C
C***********************************************************
C
      DOUBLE PRECISION COSHM1, COSHMS, DUMMY, DX1, DXS, S1,
     .                 SIG1, SIGS, T
C
      DX1 = X2 - X1
      DXS = X3 - X1
      IF (DX1*(DXS-DX1) .LE. 0.D0) GO TO 2
      SIG1 = ABS(SIGMA)
      IF (SIG1 .LT. 1.D-9) THEN
C
C SIGMA = 0:  H is the quadratic interpolant.
C
        T = (DX1/DXS)**2
        GO TO 1
      ENDIF
      SIGS = SIG1*DXS/DX1
      IF (SIGS .LE. .5D0) THEN
C
C 0 < SIG1 < SIGS .LE. .5:  compute approximations to
C   COSHM1 = COSH(SIG1)-1 and COSHMS = COSH(SIGS)-1.
C
        CALL SNHCSH (SIG1, DUMMY,COSHM1,DUMMY)
        CALL SNHCSH (SIGS, DUMMY,COSHMS,DUMMY)
        T = COSHM1/COSHMS
      ELSE
C
C SIGS > .5:  compute T = COSHM1/COSHMS.
C
        T = EXP(SIG1-SIGS)*((1.D0-EXP(-SIG1))/
     .                      (1.D0-EXP(-SIGS)))**2
      ENDIF
C
C The derivative of H at X1 is
C   T = ((Y3-Y1)*COSHM1-(Y2-Y1)*COSHMS)/
C       (DXS*COSHM1-DX1*COSHMS).
C
C ENDSLP = T unless T*S1 < 0 or abs(T) > 3*abs(S1).
C
    1 T = ((Y3-Y1)*T-Y2+Y1)/(DXS*T-DX1)
      S1 = (Y2-Y1)/DX1
      IF (S1 .GE. 0.D0) THEN
        ENDSLP = MIN(MAX(0.D0,T), 3.D0*S1)
      ELSE
        ENDSLP = MAX(MIN(0.D0,T), 3.D0*S1)
      ENDIF
      RETURN
C
C Error in the abscissae.
C
    2 ENDSLP = 0.D0
      RETURN
      END
      DOUBLE PRECISION FUNCTION HPPVAL (T,N,X,Y,YP,
     .                                  SIGMA, IER)
      INTEGER N, IER
      DOUBLE PRECISION T, X(N), Y(N), YP(N), SIGMA(N)
C
C***********************************************************
C
C                                                From TSPACK
C                                            Robert J. Renka
C                                  Dept. of Computer Science
C                                       Univ. of North Texas
C                                           renka@cs.unt.edu
C                                                   11/17/96
C
C   This function evaluates the second derivative HPP of a
C Hermite interpolatory tension spline H at a point T.
C
C On input:
C
C       T = Point at which HPP is to be evaluated.  Extrap-
C           olation is performed if T < X(1) or T > X(N).
C
C       N = Number of data points.  N .GE. 2.
C
C       X = Array of length N containing the abscissae.
C           These must be in strictly increasing order:
C           X(I) < X(I+1) for I = 1,...,N-1.
C
C       Y = Array of length N containing data values.
C           H(X(I)) = Y(I) for I = 1,...,N.
C
C       YP = Array of length N containing first deriva-
C            tives.  HP(X(I)) = YP(I) for I = 1,...,N, where
C            HP denotes the derivative of H.
C
C       SIGMA = Array of length N-1 containing tension fac-
C               tors whose absolute values determine the
C               balance between cubic and linear in each
C               interval.  SIGMA(I) is associated with int-
C               erval (I,I+1) for I = 1,...,N-1.
C
C Input parameters are not altered by this function.
C
C On output:
C
C       IER = Error indicator:
C             IER = 0 if no errors were encountered and
C                     X(1) .LE. T .LE. X(N).
C             IER = 1 if no errors were encountered and
C                     extrapolation was necessary.
C             IER = -1 if N < 2.
C             IER = -2 if the abscissae are not in strictly
C                      increasing order.  (This error will
C                      not necessarily be detected.)
C
C       HPPVAL = Second derivative value HPP(T), or zero if
C                IER < 0.
C
C Modules required by HPPVAL:  INTRVL, SNHCSH
C
C Intrinsic functions called by HPPVAL:  ABS, EXP
C
C***********************************************************
C
      INTEGER I, IP1
      DOUBLE PRECISION B1, B2, CM, CM2, CMM, COSH2, D1, D2,
     .                 DUMMY, DX, E, E1, E2, EMS, S, SB1,
     .                 SB2, SBIG, SIG, SINH2, SM, SM2, TM
      INTEGER INTRVL
C
      DATA SBIG/85.D0/
      IF (N .LT. 2) GO TO 1
C
C Find the index of the left end of an interval containing
C   T.  If T < X(1) or T > X(N), extrapolation is performed
C   using the leftmost or rightmost interval.
C
      IF (T .LT. X(1)) THEN
        I = 1
        IER = 1
      ELSEIF (T .GT. X(N)) THEN
        I = N-1
        IER = 1
      ELSE
        I = INTRVL (T,N,X)
        IER = 0
      ENDIF
      IP1 = I + 1
C
C Compute interval width DX, local coordinates B1 and B2,
C   and second differences D1 and D2.
C
      DX = X(IP1) - X(I)
      IF (DX .LE. 0.D0) GO TO 2
      B1 = (X(IP1) - T)/DX
      B2 = 1.D0 - B1
      S = (Y(IP1)-Y(I))/DX
      D1 = S - YP(I)
      D2 = YP(IP1) - S
      SIG = ABS(SIGMA(I))
      IF (SIG .LT. 1.D-9) THEN
C
C SIG = 0:  H is the Hermite cubic interpolant.
C
        HPPVAL = (D1 + D2 + 3.D0*(B2-B1)*(D2-D1))/DX
      ELSEIF (SIG .LE. .5D0) THEN
C
C 0 .LT. SIG .LE. .5:  use approximations designed to avoid
C   cancellation error in the hyperbolic functions.
C
        SB2 = SIG*B2
        CALL SNHCSH (SIG, SM,CM,CMM)
        CALL SNHCSH (SB2, SM2,CM2,DUMMY)
        SINH2 = SM2 + SB2
        COSH2 = CM2 + 1.D0
        E = SIG*SM - CMM - CMM
        HPPVAL = SIG*((CM*SINH2-SM*COSH2)*(D1+D2) +
     .              SIG*(CM*COSH2-(SM+SIG)*SINH2)*D1)/(DX*E)
      ELSE
C
C SIG > .5:  use negative exponentials in order to avoid
C   overflow.  Note that EMS = EXP(-SIG).  In the case of
C   extrapolation (negative B1 or B2), H is approximated by
C   a linear function if -SIG*B1 or -SIG*B2 is large.
C
        SB1 = SIG*B1
        SB2 = SIG - SB1
        IF (-SB1 .GT. SBIG  .OR.  -SB2 .GT. SBIG) THEN
          HPPVAL = 0.D0
        ELSE
          E1 = EXP(-SB1)
          E2 = EXP(-SB2)
          EMS = E1*E2
          TM = 1.D0 - EMS
          E = TM*(SIG*(1.D0+EMS) - TM - TM)
          HPPVAL = SIG*(SIG*((E1*EMS+E2)*D1+(E1+E2*EMS)*D2)-
     .                       TM*(E1+E2)*(D1+D2))/(DX*E)
        ENDIF
      ENDIF
      RETURN
C
C N is outside its valid range.
C
    1 HPPVAL = 0.D0
      IER = -1
      RETURN
C
C X(I) .GE. X(I+1).
C
    2 HPPVAL = 0.D0
      IER = -2
      RETURN
      END
      DOUBLE PRECISION FUNCTION HPVAL (T,N,X,Y,YP,
     .                                 SIGMA, IER)
      INTEGER N, IER
      DOUBLE PRECISION T, X(N), Y(N), YP(N), SIGMA(N)
C
C***********************************************************
C
C                                                From TSPACK
C                                            Robert J. Renka
C                                  Dept. of Computer Science
C                                       Univ. of North Texas
C                                           renka@cs.unt.edu
C                                                   11/17/96
C
C   This function evaluates the first derivative HP of a
C Hermite interpolatory tension spline H at a point T.
C
C On input:
C
C       T = Point at which HP is to be evaluated.  Extrapo-
C           lation is performed if T < X(1) or T > X(N).
C
C       N = Number of data points.  N .GE. 2.
C
C       X = Array of length N containing the abscissae.
C           These must be in strictly increasing order:
C           X(I) < X(I+1) for I = 1,...,N-1.
C
C       Y = Array of length N containing data values.
C           H(X(I)) = Y(I) for I = 1,...,N.
C
C       YP = Array of length N containing first deriva-
C            tives.  HP(X(I)) = YP(I) for I = 1,...,N.
C
C       SIGMA = Array of length N-1 containing tension fac-
C               tors whose absolute values determine the
C               balance between cubic and linear in each
C               interval.  SIGMA(I) is associated with int-
C               erval (I,I+1) for I = 1,...,N-1.
C
C Input parameters are not altered by this function.
C
C On output:
C
C       IER = Error indicator:
C             IER = 0  if no errors were encountered and
C                      X(1) .LE. T .LE. X(N).
C             IER = 1  if no errors were encountered and
C                      extrapolation was necessary.
C             IER = -1 if N < 2.
C             IER = -2 if the abscissae are not in strictly
C                      increasing order.  (This error will
C                      not necessarily be detected.)
C
C       HPVAL = Derivative value HP(T), or zero if IER < 0.
C
C Modules required by HPVAL:  INTRVL, SNHCSH
C
C Intrinsic functions called by HPVAL:  ABS, EXP
C
C***********************************************************
C
      INTEGER I, IP1
      DOUBLE PRECISION B1, B2, CM, CM2, CMM, D1, D2, DUMMY,
     .                 DX, E, E1, E2, EMS, S, S1, SB1, SB2,
     .                 SBIG, SIG, SINH2, SM, SM2, TM
      INTEGER INTRVL
C
      DATA SBIG/85.D0/
      IF (N .LT. 2) GO TO 1
C
C Find the index of the left end of an interval containing
C   T.  If T < X(1) or T > X(N), extrapolation is performed
C   using the leftmost or rightmost interval.
C
      IF (T .LT. X(1)) THEN
        I = 1
        IER = 1
      ELSEIF (T .GT. X(N)) THEN
        I = N-1
        IER = 1
      ELSE
        I = INTRVL (T,N,X)
        IER = 0
      ENDIF
      IP1 = I + 1
C
C Compute interval width DX, local coordinates B1 and B2,
C   and second differences D1 and D2.
C
      DX = X(IP1) - X(I)
      IF (DX .LE. 0.D0) GO TO 2
      B1 = (X(IP1) - T)/DX
      B2 = 1.D0 - B1
      S1 = YP(I)
      S = (Y(IP1)-Y(I))/DX
      D1 = S - S1
      D2 = YP(IP1) - S
      SIG = ABS(SIGMA(I))
      IF (SIG .LT. 1.D-9) THEN
C
C SIG = 0:  H is the Hermite cubic interpolant.
C
        HPVAL = S1 + B2*(D1 + D2 - 3.D0*B1*(D2-D1))
      ELSEIF (SIG .LE. .5D0) THEN
C
C 0 .LT. SIG .LE. .5:  use approximations designed to avoid
C   cancellation error in the hyperbolic functions.
C
        SB2 = SIG*B2
        CALL SNHCSH (SIG, SM,CM,CMM)
        CALL SNHCSH (SB2, SM2,CM2,DUMMY)
        SINH2 = SM2 + SB2
        E = SIG*SM - CMM - CMM
        HPVAL = S1 + ((CM*CM2-SM*SINH2)*(D1+D2) +
     .                SIG*(CM*SINH2-(SM+SIG)*CM2)*D1)/E
      ELSE
C
C SIG > .5:  use negative exponentials in order to avoid
C   overflow.  Note that EMS = EXP(-SIG).  In the case of
C   extrapolation (negative B1 or B2), H is approximated by
C   a linear function if -SIG*B1 or -SIG*B2 is large.
C
        SB1 = SIG*B1
        SB2 = SIG - SB1
        IF (-SB1 .GT. SBIG  .OR.  -SB2 .GT. SBIG) THEN
          HPVAL = S
        ELSE
          E1 = EXP(-SB1)
          E2 = EXP(-SB2)
          EMS = E1*E2
          TM = 1.D0 - EMS
          E = TM*(SIG*(1.D0+EMS) - TM - TM)
          HPVAL = S + (TM*((E2-E1)*(D1+D2) + TM*(D1-D2)) +
     .            SIG*((E1*EMS-E2)*D1 + (E1-E2*EMS)*D2))/E
        ENDIF
      ENDIF
      RETURN
C
C N is outside its valid range.
C
    1 HPVAL = 0.D0
      IER = -1
      RETURN
C
C X(I) .GE. X(I+1).
C
    2 HPVAL = 0.D0
      IER = -2
      RETURN
      END
      DOUBLE PRECISION FUNCTION HVAL (T,N,X,Y,YP,SIGMA, IER)
      INTEGER N, IER
      DOUBLE PRECISION T, X(N), Y(N), YP(N), SIGMA(N)
C
C***********************************************************
C
C                                                From TSPACK
C                                            Robert J. Renka
C                                  Dept. of Computer Science
C                                       Univ. of North Texas
C                                           renka@cs.unt.edu
C                                                   11/17/96
C
C   This function evaluates a Hermite interpolatory tension
C spline H at a point T.  Note that a large value of SIGMA
C may cause underflow.  The result is assumed to be zero.
C
C   Given arrays X, Y, YP, and SIGMA of length NN, if T is
C known to lie in the interval (X(I),X(J)) for some I < J,
C a gain in efficiency can be achieved by calling this
C function with N = J+1-I (rather than NN) and the I-th
C components of the arrays (rather than the first) as par-
C ameters.
C
C On input:
C
C       T = Point at which H is to be evaluated.  Extrapo-
C           lation is performed if T < X(1) or T > X(N).
C
C       N = Number of data points.  N .GE. 2.
C
C       X = Array of length N containing the abscissae.
C           These must be in strictly increasing order:
C           X(I) < X(I+1) for I = 1,...,N-1.
C
C       Y = Array of length N containing data values.
C           H(X(I)) = Y(I) for I = 1,...,N.
C
C       YP = Array of length N containing first deriva-
C            tives.  HP(X(I)) = YP(I) for I = 1,...,N, where
C            HP denotes the derivative of H.
C
C       SIGMA = Array of length N-1 containing tension fac-
C               tors whose absolute values determine the
C               balance between cubic and linear in each
C               interval.  SIGMA(I) is associated with int-
C               erval (I,I+1) for I = 1,...,N-1.
C
C Input parameters are not altered by this function.
C
C On output:
C
C       IER = Error indicator:
C             IER = 0  if no errors were encountered and
C                      X(1) .LE. T .LE. X(N).
C             IER = 1  if no errors were encountered and
C                      extrapolation was necessary.
C             IER = -1 if N < 2.
C             IER = -2 if the abscissae are not in strictly
C                      increasing order.  (This error will
C                      not necessarily be detected.)
C
C       HVAL = Function value H(T), or zero if IER < 0.
C
C Modules required by HVAL:  INTRVL, SNHCSH
C
C Intrinsic functions called by HVAL:  ABS, EXP
C
C***********************************************************
C
      INTEGER I, IP1
      DOUBLE PRECISION B1, B2, CM, CM2, CMM, D1, D2, DUMMY,
     .                 DX, E, E1, E2, EMS, S, S1, SB1, SB2,
     .                 SBIG, SIG, SM, SM2, TM, TP, TS, U, Y1
      INTEGER INTRVL
C
      DATA SBIG/85.D0/
      IF (N .LT. 2) GO TO 1
C
C Find the index of the left end of an interval containing
C   T.  If T < X(1) or T > X(N), extrapolation is performed
C   using the leftmost or rightmost interval.
C
      IF (T .LT. X(1)) THEN
        I = 1
        IER = 1
      ELSEIF (T .GT. X(N)) THEN
        I = N-1
        IER = 1
      ELSE
        I = INTRVL (T,N,X)
        IER = 0
      ENDIF
      IP1 = I + 1
C
C Compute interval width DX, local coordinates B1 and B2,
C   and second differences D1 and D2.
C
      DX = X(IP1) - X(I)
      IF (DX .LE. 0.D0) GO TO 2
      U = T - X(I)
      B2 = U/DX
      B1 = 1.D0 - B2
      Y1 = Y(I)
      S1 = YP(I)
      S = (Y(IP1)-Y1)/DX
      D1 = S - S1
      D2 = YP(IP1) - S
      SIG = ABS(SIGMA(I))
      IF (SIG .LT. 1.D-9) THEN
C
C SIG = 0:  H is the Hermite cubic interpolant.
C
        HVAL = Y1 + U*(S1 + B2*(D1 + B1*(D1-D2)))
      ELSEIF (SIG .LE. .5D0) THEN
C
C 0 .LT. SIG .LE. .5:  use approximations designed to avoid
C   cancellation error in the hyperbolic functions.
C
        SB2 = SIG*B2
        CALL SNHCSH (SIG, SM,CM,CMM)
        CALL SNHCSH (SB2, SM2,CM2,DUMMY)
        E = SIG*SM - CMM - CMM
        HVAL = Y1 + S1*U + DX*((CM*SM2-SM*CM2)*(D1+D2) +
     .                         SIG*(CM*CM2-(SM+SIG)*SM2)*D1)
     .                         /(SIG*E)
      ELSE
C
C SIG > .5:  use negative exponentials in order to avoid
C   overflow.  Note that EMS = EXP(-SIG).  In the case of
C   extrapolation (negative B1 or B2), H is approximated by
C   a linear function if -SIG*B1 or -SIG*B2 is large.
C
        SB1 = SIG*B1
        SB2 = SIG - SB1
        IF (-SB1 .GT. SBIG  .OR.  -SB2 .GT. SBIG) THEN
          HVAL = Y1 + S*U
        ELSE
          E1 = EXP(-SB1)
          E2 = EXP(-SB2)
          EMS = E1*E2
          TM = 1.D0 - EMS
          TS = TM*TM
          TP = 1.D0 + EMS
          E = TM*(SIG*TP - TM - TM)
          HVAL = Y1 + S*U + DX*(TM*(TP-E1-E2)*(D1+D2) + SIG*
     .                         ((E2+EMS*(E1-2.D0)-B1*TS)*D1+
     .                        (E1+EMS*(E2-2.D0)-B2*TS)*D2))/
     .                        (SIG*E)
        ENDIF
      ENDIF
      RETURN
C
C N is outside its valid range.
C
    1 HVAL = 0.D0
      IER = -1
      RETURN
C
C X(I) .GE. X(I+1).
C
    2 HVAL = 0.D0
      IER = -2
      RETURN
      END
      INTEGER FUNCTION INTRVL (T,N,X)
      INTEGER N
      DOUBLE PRECISION T, X(N)
C
C***********************************************************
C
C                                                From TSPACK
C                                            Robert J. Renka
C                                  Dept. of Computer Science
C                                       Univ. of North Texas
C                                           renka@cs.unt.edu
C                                                   06/10/92
C
C   This function returns the index of the left end of an
C interval (defined by an increasing sequence X) which
C contains the value T.  The method consists of first test-
C ing the interval returned by a previous call, if any, and
C then using a binary search if necessary.
C
C On input:
C
C       T = Point to be located.
C
C       N = Length of X.  N .GE. 2.
C
C       X = Array of length N assumed (without a test) to
C           contain a strictly increasing sequence of
C           values.
C
C Input parameters are not altered by this function.
C
C On output:
C
C       INTRVL = Index I defined as follows:
C
C                  I = 1    if  T .LT. X(2) or N .LE. 2,
C                  I = N-1  if  T .GE. X(N-1), and
C                  X(I) .LE. T .LT. X(I+1) otherwise.
C
C Modules required by INTRVL:  None
C
C***********************************************************
C
      INTEGER IH, IL, K
      DOUBLE PRECISION TT
C
      SAVE IL
      DATA IL/1/
      TT = T
      IF (IL .GE. 1  .AND.  IL .LT. N) THEN
        IF (X(IL) .LE. TT  .AND.  TT .LT. X(IL+1)) GO TO 2
      ENDIF
C
C Initialize low and high indexes.
C
      IL = 1
      IH = N
C
C Binary search:
C
    1 IF (IH .LE. IL+1) GO TO 2
        K = (IL+IH)/2
        IF (TT .LT. X(K)) THEN
          IH = K
        ELSE
          IL = K
        ENDIF
        GO TO 1
C
C X(IL) .LE. T .LT. X(IL+1)  or  (T .LT. X(1) and IL=1)
C                            or  (T .GE. X(N) and IL=N-1)
C
    2 INTRVL = IL
      RETURN
      END
      DOUBLE PRECISION FUNCTION SIG0 (X1,X2,Y1,Y2,Y1P,Y2P,
     .                                IFL,HBND,TOL, IER)
      INTEGER IFL, IER
      DOUBLE PRECISION X1, X2, Y1, Y2, Y1P, Y2P, HBND, TOL
C
C***********************************************************
C
C                                                From TSPACK
C                                            Robert J. Renka
C                                  Dept. of Computer Science
C                                       Univ. of North Texas
C                                           renka@cs.unt.edu
C                                                   11/18/96
C
C   Given a pair of abscissae with associated ordinates and
C slopes, this function determines the smallest (nonnega-
C tive) tension factor SIGMA such that the Hermite interpo-
C latory tension spline H(x), defined by SIGMA and the data,
C is bounded (either above or below) by HBND for all x in
C (X1,X2).
C
C On input:
C
C       X1,X2 = Abscissae.  X1 < X2.
C
C       Y1,Y2 = Values of H at X1 and X2.
C
C       Y1P,Y2P = Derivative values of H at X1 and X2.
C
C       IFL = Option indicator:
C             IFL = -1 if HBND is a lower bound on H.
C             IFL = 1 if HBND is an upper bound on H.
C
C       HBND = Bound on H.  If IFL = -1, HBND .LE. min(Y1,
C              Y2).  If IFL = 1, HBND .GE. max(Y1,Y2).
C
C       TOL = Tolerance whose magnitude determines how close
C             SIGMA is to its optimal value when nonzero
C             finite tension is necessary and sufficient to
C             satisfy the constraint.  For a lower bound,
C             SIGMA is chosen so that HBND .LE. HMIN .LE.
C             HBND + abs(TOL), where HMIN is the minimum
C             value of H in the interval, and for an upper
C             bound, the maximum of H satisfies HBND -
C             abs(TOL) .LE. HMAX .LE. HBND.  Thus, the con-
C             straint is satisfied but possibly with more
C             tension than necessary.
C
C Input parameters are not altered by this function.
C
C On output:
C
C       IER = Error indicator:
C             IER = 0 if no errors were encountered and the
C                     constraint can be satisfied with fin-
C                     ite tension.
C             IER = 1 if no errors were encountered but in-
C                     finite tension is required to satisfy
C                     the constraint (e.g., IFL = -1, HBND
C                     = Y1, and Y1P < 0.).
C             IER = -1 if X2 .LE. X1 or abs(IFL) .NE. 1.
C             IER = -2 if HBND is outside its valid range
C                      on input.
C
C       SIG0 = Minimum tension factor defined above unless
C              IER < 0, in which case SIG0 = -1.  If IER =
C              1, SIG0 = 85, resulting in an approximation
C              to the linear interpolant of the endpoint
C              values.  Note, however, that SIG0 may be
C              larger than 85 if IER = 0.
C
C Modules required by SIG0:  SNHCSH, STORE
C
C Intrinsic functions called by SIG0:  ABS, DBLE, EXP, LOG,
C                                        MAX, MIN, SIGN,
C                                        SQRT
C
C***********************************************************
C
      INTEGER LUN, NIT
      DOUBLE PRECISION A, A0, AA, B, B0, BND, C, C1, C2,
     .                 COSHM, COSHMM, D, D0, D1PD2, D2,
     .                 DMAX, DSIG, DX, E, EMS, F, F0, FMAX,
     .                 FNEG, FTOL, R, RF, RSIG, RTOL, S, S1,
     .                 S2, SBIG, SCM, SIG, SINHM, SNEG,
     .                 SSINH, SSM, STOL, T, T0, T1, T2, TM,
     .                 Y1L, Y2L
      DOUBLE PRECISION STORE
C
      DATA SBIG/85.D0/,  LUN/-1/
C
C Store local parameters and test for errors.
C
      RF = DBLE(IFL)
      DX = X2 - X1
      IF (ABS(RF) .NE. 1.D0  .OR.  DX .LE. 0.D0) GO TO 8
      Y1L = Y1
      Y2L = Y2
      BND = HBND
C
C Test for a valid constraint.
C
      IF ((RF .LT. 0.D0  .AND.  MIN(Y1L,Y2L) .LT. BND)
     .    .OR.  (RF .GT. 0.D0  .AND.
     .           BND .LT. MAX(Y1L,Y2L))) GO TO 9
C
C Test for infinite tension required.
C
      S1 = Y1P
      S2 = Y2P
      IF ((Y1L .EQ. BND  .AND.  RF*S1 .GT. 0.D0)  .OR.
     .    (Y2L .EQ. BND  .AND.  RF*S2 .LT. 0.D0)) GO TO 7
C
C Test for SIG = 0 sufficient.
C
      SIG = 0.D0
      IF (RF*S1 .LE. 0.D0  .AND.  RF*S2 .GE. 0.D0) GO TO 6
C
C   Compute coefficients A0 and B0 of the Hermite cubic in-
C     terpolant H0(x) = Y2 - DX*(S2*R + B0*R**2 + A0*R**3/3)
C     where R = (X2-x)/DX.
C
      S = (Y2L-Y1L)/DX
      T0 = 3.D0*S - S1 - S2
      A0 = 3.D0*(S-T0)
      B0 = T0 - S2
      D0 = T0*T0 - S1*S2
C
C   H0 has local extrema in (X1,X2) iff S1*S2 < 0 or
C     (T0*(S1+S2) < 0 and D0 .GE. 0).
C
      IF (S1*S2 .GE. 0.D0  .AND.  (T0*(S1+S2) .GE. 0.D0
     .    .OR.  D0 .LT. 0.D0)) GO TO 6
      IF (A0 .EQ. 0.D0) THEN
C
C   H0 is quadratic and has an extremum at R = -S2/(2*B0).
C     H0(R) = Y2 + DX*S2**2/(4*B0).  Note that A0 = 0 im-
C     plies 2*B0 = S1-S2, and S1*S2 < 0 implies B0 .NE. 0.
C     Also, the extremum is a min iff HBND is a lower bound.
C
        F0 = (BND - Y2L - DX*S2*S2/(4.D0*B0))*RF
      ELSE
C
C   A0 .NE. 0 and H0 has extrema at R = (-B0 +/- SQRT(D0))/
C     A0 = S2/(-B0 -/+ SQRT(D0)), where the negative root
C     corresponds to a min.  The expression for R is chosen
C     to avoid cancellation error.  H0(R) = Y2 + DX*(S2*B0 +
C     2*D0*R)/(3*A0).
C
        T = -B0 - SIGN(SQRT(D0),B0)
        R = T/A0
        IF (RF*B0 .GT. 0.D0) R = S2/T
        F0 = (BND - Y2L - DX*(S2*B0+2.D0*D0*R)/(3.D0*A0))*RF
      ENDIF
C
C   F0 .GE. 0 iff SIG = 0 is sufficient to satisfy the
C     constraint.
C
      IF (F0 .GE. 0.D0) GO TO 6
C
C Find a zero of F(SIG) = (BND-H(R))*RF where the derivative
C   of H, HP, vanishes at R.  F is a nondecreasing function,
C   F(0) < 0, and F = FMAX for SIG sufficiently large.
C
C Initialize parameters for the secant method.  The method
C   uses three points:  (SG0,F0), (SIG,F), and (SNEG,FNEG),
C   where SG0 and SNEG are defined implicitly by DSIG = SIG
C   - SG0 and DMAX = SIG - SNEG.  SG0 is initially zero and
C   SNEG is initialized to a sufficiently large value that
C   FNEG > 0.  This value is used only if the initial value
C   of F is negative.
C
      FMAX = MAX(1.D-3,MIN(ABS(Y1L-BND),ABS(Y2L-BND)))
      T = MAX(ABS(Y1L-BND),ABS(Y2L-BND))
      SIG = DX*MAX(ABS(S1),ABS(S2))/T
      DMAX = SIG*(1.D0-T/FMAX)
      SNEG = SIG - DMAX
      IF (LUN .GE. 0  .AND.  RF .LT. 0.D0)
     .   WRITE (LUN,100) F0, FMAX, SNEG
      IF (LUN .GE. 0  .AND.  RF .GT. 0.D0)
     .   WRITE (LUN,110) F0, FMAX, SNEG
  100 FORMAT (//1X,'SIG0 (LOWER BOUND) -- F(0) = ',D15.8,
     .        ', FMAX = ',D15.8/1X,46X,'SNEG = ',D15.8/)
  110 FORMAT (//1X,'SIG0 (UPPER BOUND) -- F(0) = ',D15.8,
     .        ', FMAX = ',D15.8/1X,46X,'SNEG = ',D15.8/)
      DSIG = SIG
      FNEG = FMAX
      D2 = S2 - S
      D1PD2 = S2 - S1
      NIT = 0
C
C Compute an absolute tolerance FTOL = abs(TOL), and a
C   relative tolerance RTOL = 100*MACHEPS.
C
      FTOL = ABS(TOL)
      RTOL = 1.D0
    1 RTOL = RTOL/2.D0
        IF (STORE(RTOL+1.D0) .GT. 1.D0) GO TO 1
      RTOL = RTOL*200.D0
C
C Top of loop:  compute F.
C
    2 EMS = EXP(-SIG)
      IF (SIG .LE. .5D0) THEN
C
C   SIG .LE. .5:  use approximations designed to avoid can-
C                 cellation error (associated with small
C                 SIG) in the modified hyperbolic functions.
C
        CALL SNHCSH (SIG, SINHM,COSHM,COSHMM)
        C1 = SIG*COSHM*D2 - SINHM*D1PD2
        C2 = SIG*(SINHM+SIG)*D2 - COSHM*D1PD2
        A = C2 - C1
        AA = A/EMS
        E = SIG*SINHM - COSHMM - COSHMM
      ELSE
C
C   SIG > .5:  scale SINHM and COSHM by 2*EXP(-SIG) in order
C              to avoid overflow.
C
        TM = 1.D0 - EMS
        SSINH = TM*(1.D0+EMS)
        SSM = SSINH - 2.D0*SIG*EMS
        SCM = TM*TM
        C1 = SIG*SCM*D2 - SSM*D1PD2
        C2 = SIG*SSINH*D2 - SCM*D1PD2
        AA = 2.D0*(SIG*TM*D2 + (TM-SIG)*D1PD2)
        A = EMS*AA
        E = SIG*SSINH - SCM - SCM
      ENDIF
C
C   HP(R) = S2 - (C1*SINH(SIG*R) - C2*COSHM(SIG*R))/E = 0
C     for ESR = (-B +/- SQRT(D))/A = C/(-B -/+ SQRT(D)),
C     where ESR = EXP(SIG*R), A = C2-C1, D = B**2 - A*C, and
C     B and C are defined below.
C
      B = E*S2 - C2
      C = C2 + C1
      D = B*B - A*C
      F = 0.D0
      IF (AA*C .EQ. 0.D0  .AND.  B .EQ. 0.D0) GO TO 3
      F = FMAX
      IF (D .LT. 0.D0) GO TO 3
      T1 = SQRT(D)
      T = -B - SIGN(T1,B)
      RSIG = 0.D0
      IF (RF*B .LT. 0.D0  .AND.  AA .NE. 0.) THEN
        IF (T/AA .GT. 0.D0) RSIG = SIG + LOG(T/AA)
      ENDIF
      IF ((RF*B .GT. 0.D0  .OR.  AA .EQ. 0.D0)  .AND.
     .    C/T .GT. 0.D0) RSIG = LOG(C/T)
      IF ((RSIG .LE. 0.D0  .OR.  RSIG .GE. SIG)  .AND.
     .    B .NE. 0.D0) GO TO 3
C
C   H(R) = Y2 - DX*(B*SIG*R + C1 + RF*SQRT(D))/(SIG*E).
C
      F = (BND - Y2L + DX*(B*RSIG+C1+RF*T1)/(SIG*E))*RF
C
C   Update the number of iterations NIT.
C
    3 NIT = NIT + 1
      IF (LUN .GE. 0) WRITE (LUN,120) NIT, SIG, F
  120 FORMAT (1X,3X,I2,' -- SIG = ',D15.8,', F = ',
     .        D15.8)
      IF (F0*F .LT. 0.D0) THEN
C
C   F0*F < 0.  Update (SNEG,FNEG) to (SG0,F0) so that F
C     and FNEG always have opposite signs.  If SIG is
C     closer to SNEG than SG0, then swap (SNEG,FNEG) with
C     (SG0,F0).
C
        T1 = DMAX
        T2 = FNEG
        DMAX = DSIG
        FNEG = F0
        IF (ABS(DSIG) .GT. ABS(T1)) THEN
          DSIG = T1
          F0 = T2
        ENDIF
      ENDIF
C
C   Test for convergence.
C
      STOL = RTOL*SIG
      IF (ABS(DMAX) .LE. STOL  .OR.  (F .GE. 0.D0  .AND.
     .    F .LE. FTOL)  .OR.  ABS(F) .LE. RTOL) GO TO 6
C
C   Test for F0 = F = FMAX or F < 0 on the first iteration.
C
      IF (F0 .NE. F  .AND.  (NIT .GT. 1  .OR.  F .GT. 0.D0))
     .   GO TO 5
C
C   F*F0 > 0 and either the new estimate would be outside of
C     the bracketing interval of length abs(DMAX) or F < 0
C     on the first iteration.  Reset (SG0,F0) to (SNEG,FNEG).
C
    4 DSIG = DMAX
      F0 = FNEG
C
C   Compute the change in SIG by linear interpolation be-
C     tween (SG0,F0) and (SIG,F).
C
    5 DSIG = -F*DSIG/(F-F0)
      IF (LUN .GE. 0) WRITE (LUN,130) DSIG
  130 FORMAT (1X,8X,'DSIG = ',D15.8)
      IF ( ABS(DSIG) .GT. ABS(DMAX)  .OR.
     .     DSIG*DMAX .GT. 0. ) GO TO 4
C
C   Restrict the step-size such that abs(DSIG) .GE. STOL/2.
C     Note that DSIG and DMAX have opposite signs.
C
      IF (ABS(DSIG) .LT. STOL/2.D0)
     .  DSIG = -SIGN(STOL/2.D0,DMAX)
C
C   Bottom of loop:  update SIG, DMAX, and F0.
C
      SIG = SIG + DSIG
      DMAX = DMAX + DSIG
      F0 = F
      GO TO 2
C
C No errors encountered and SIGMA finite.
C
    6 IER = 0
      SIG0 = SIG
      RETURN
C
C Infinite tension required.
C
    7 IER = 1
      SIG0 = SBIG
      RETURN
C
C Error in an input parameter.
C
    8 IER = -1
      SIG0 = -1.D0
      RETURN
C
C Invalid constraint.
C
    9 IER = -2
      SIG0 = -1.D0
      RETURN
      END
      DOUBLE PRECISION FUNCTION SIG1 (X1,X2,Y1,Y2,Y1P,Y2P,
     .                                IFL,HPBND,TOL, IER)
      INTEGER IFL, IER
      DOUBLE PRECISION X1, X2, Y1, Y2, Y1P, Y2P, HPBND, TOL
C
C***********************************************************
C
C                                                From TSPACK
C                                            Robert J. Renka
C                                  Dept. of Computer Science
C                                       Univ. of North Texas
C                                           renka@cs.unt.edu
C                                                   11/17/96
C
C   Given a pair of abscissae with associated ordinates and
C slopes, this function determines the smallest (nonnega-
C tive) tension factor SIGMA such that the derivative HP(x)
C of the Hermite interpolatory tension spline H(x), defined
C by SIGMA and the data, is bounded (either above or below)
C by HPBND for all x in (X1,X2).
C
C On input:
C
C       X1,X2 = Abscissae.  X1 < X2.
C
C       Y1,Y2 = Values of H at X1 and X2.
C
C       Y1P,Y2P = Values of HP at X1 and X2.
C
C       IFL = Option indicator:
C             IFL = -1 if HPBND is a lower bound on HP.
C             IFL = 1 if HPBND is an upper bound on HP.
C
C       HPBND = Bound on HP.  If IFL = -1, HPBND .LE.
C               min(Y1P,Y2P,S) for S = (Y2-Y1)/(X2-X1).  If
C               IFL = 1, HPBND .GE. max(Y1P,Y2P,S).
C
C       TOL = Tolerance whose magnitude determines how close
C             SIGMA is to its optimal value when nonzero
C             finite tension is necessary and sufficient to
C             satisfy the constraint.  For a lower bound,
C             SIGMA is chosen so that HPBND .LE. HPMIN .LE.
C             HPBND + abs(TOL), where HPMIN is the minimum
C             value of HP in the interval, and for an upper
C             bound, the maximum of HP satisfies HPBND -
C             abs(TOL) .LE. HPMAX .LE. HPBND.  Thus, the
C             constraint is satisfied but possibly with more
C             tension than necessary.
C
C Input parameters are not altered by this function.
C
C On output:
C
C       IER = Error indicator:
C             IER = 0 if no errors were encountered and the
C                     constraint can be satisfied with fin-
C                     ite tension.
C             IER = 1 if no errors were encountered but in-
C                     finite tension is required to satisfy
C                     the constraint (e.g., IFL = -1, HPBND
C                     = S, and Y1P > S).
C             IER = -1 if X2 .LE. X1 or abs(IFL) .NE. 1.
C             IER = -2 if HPBND is outside its valid range
C                      on input.
C
C       SIG1 = Minimum tension factor defined above unless
C              IER < 0, in which case SIG1 = -1.  If IER =
C              1, SIG1 = 85, resulting in an approximation
C              to the linear interpolant of the endpoint
C              values.  Note, however, that SIG1 may be
C              larger than 85 if IER = 0.
C
C Modules required by SIG1:  SNHCSH, STORE
C
C Intrinsic functions called by SIG1:  ABS, DBLE, EXP, MAX,
C                                        MIN, SIGN, SQRT
C
C***********************************************************
C
      INTEGER LUN, NIT
      DOUBLE PRECISION A, A0, B0, BND, C0, C1, C2, COSHM,
     .                 COSHMM, D0, D1, D1PD2, D2, DMAX,
     .                 DSIG, DX, E, EMS, EMS2, F, F0, FMAX,
     .                 FNEG, FTOL, RF, RTOL, S, S1, S2,
     .                 SBIG, SIG, SINH, SINHM, STOL, T0, T1,
     .                 T2, TM
      DOUBLE PRECISION STORE
C
      DATA SBIG/85.D0/,  LUN/-1/
C
C Store local parameters and test for errors.
C
      RF = DBLE(IFL)
      DX = X2 - X1
      IF (ABS(RF) .NE. 1.D0  .OR.  DX .LE. 0.D0) GO TO 7
      S1 = Y1P
      S2 = Y2P
      S = (Y2-Y1)/DX
      BND = HPBND
C
C Test for a valid constraint.
C
      IF ((RF .LT. 0.D0  .AND.  MIN(S1,S2,S) .LT. BND)
     .    .OR.  (RF .GT. 0.D0  .AND.
     .           BND .LT. MAX(S1,S2,S))) GO TO 8
C
C Test for infinite tension required.
C
      IF (S .EQ. BND  .AND.  (S1 .NE. S  .OR.  S2 .NE. S))
     .   GO TO 6
C
C Test for SIG = 0 sufficient.  The Hermite cubic interpo-
C   land H0 has derivative HP0(x) = S2 + 2*B0*R + A0*R**2,
C   where R = (X2-x)/DX.
C
      SIG = 0.D0
      T0 = 3.D0*S - S1 - S2
      B0 = T0 - S2
      C0 = T0 - S1
      A0 = -B0 - C0
C
C   HP0(R) has an extremum (at R = -B0/A0) in (0,1) iff
C     B0*C0 > 0 and the third derivative of H0 has the
C     sign of A0.
C
      IF (B0*C0 .LE. 0.D0  .OR.  A0*RF .GT. 0.D0) GO TO 5
C
C   A0*RF < 0 and HP0(R) = -D0/A0 at R = -B0/A0.
C
      D0 = T0*T0 - S1*S2
      F0 = (BND + D0/A0)*RF
      IF (F0 .GE. 0.D0) GO TO 5
C
C Find a zero of F(SIG) = (BND-HP(R))*RF, where HP has an
C   extremum at R.  F has a unique zero, F(0) = F0 < 0, and
C   F = (BND-S)*RF > 0 for SIG sufficiently large.
C
C Initialize parameters for the secant method.  The method
C   uses three points:  (SG0,F0), (SIG,F), and (SNEG,FNEG),
C   where SG0 and SNEG are defined implicitly by DSIG = SIG
C   - SG0 and DMAX = SIG - SNEG.  SG0 is initially zero and
C   SIG is initialized to the zero of (BND - (SIG*S-S1-S2)/
C   (SIG-2.))*RF -- a value for which F(SIG) .GE. 0 and
C   F(SIG) = 0 for SIG sufficiently large that 2*SIG is in-
C   significant relative to EXP(SIG).
C
      FMAX = (BND-S)*RF
      IF (LUN .GE. 0  .AND.  RF .LT. 0.D0)
     .  WRITE (LUN,100) F0, FMAX
      IF (LUN .GE. 0  .AND.  RF .GT. 0.D0)
     .  WRITE (LUN,110) F0, FMAX
  100 FORMAT (//1X,'SIG1 (LOWER BOUND) -- F(0) = ',D15.8,
     .        ', FMAX = ',D15.8/)
  110 FORMAT (//1X,'SIG1 (UPPER BOUND) -- F(0) = ',D15.8,
     .        ', FMAX = ',D15.8/)
      SIG = 2.D0 - A0/(3.D0*(BND-S))
      IF (STORE(SIG*EXP(-SIG)+.5D0) .EQ. .5D0) GO TO 5
      DSIG = SIG
      DMAX = -2.D0*SIG
      FNEG = FMAX
      D1 = S - S1
      D2 = S2 - S
      D1PD2 = D1 + D2
      NIT = 0
C
C Compute an absolute tolerance FTOL = abs(TOL), and a
C   relative tolerance RTOL = 100*MACHEPS.
C
      FTOL = ABS(TOL)
      RTOL = 1.D0
    1 RTOL = RTOL/2.D0
        IF (STORE(RTOL+1.D0) .GT. 1.D0) GO TO 1
      RTOL = RTOL*200.D0
C
C Top of loop:  compute F.
C
    2 IF (SIG .LE. .5D0) THEN
C
C   Use approximations designed to avoid cancellation error
C     (associated with small SIG) in the modified hyperbolic
C     functions.
C
        CALL SNHCSH (SIG, SINHM,COSHM,COSHMM)
        C1 = SIG*COSHM*D2 - SINHM*D1PD2
        C2 = SIG*(SINHM+SIG)*D2 - COSHM*D1PD2
        A = C2 - C1
        E = SIG*SINHM - COSHMM - COSHMM
      ELSE
C
C   Scale SINHM and COSHM by 2*EXP(-SIG) in order to avoid
C     overflow.
C
        EMS = EXP(-SIG)
        EMS2 = EMS + EMS
        TM = 1.D0 - EMS
        SINH = TM*(1.D0+EMS)
        SINHM = SINH - SIG*EMS2
        COSHM = TM*TM
        C1 = SIG*COSHM*D2 - SINHM*D1PD2
        C2 = SIG*SINH*D2 - COSHM*D1PD2
        A = EMS2*(SIG*TM*D2 + (TM-SIG)*D1PD2)
        E = SIG*SINH - COSHM - COSHM
      ENDIF
C
C   The second derivative of H(R) has a zero at EXP(SIG*R) =
C     SQRT((C2+C1)/A) and R is in (0,1) and well-defined
C     iff HPP(X1)*HPP(X2) < 0.
C
      F = FMAX
      T1 = A*(C2+C1)
      IF (T1 .GE. 0.D0) THEN
        IF (C1*(SIG*COSHM*D1 - SINHM*D1PD2) .LT. 0.D0) THEN
C
C   HP(R) = (B+SIGN(A)*SQRT(A*C))/E at the critical value
C     of R, where A = C2-C1, B = E*S2-C2, and C = C2+C1.
C     NOTE THAT RF*A < 0.
C
          F = (BND - (E*S2-C2 - RF*SQRT(T1))/E)*RF
        ENDIF
      ENDIF
C
C   Update the number of iterations NIT.
C
      NIT = NIT + 1
      IF (LUN .GE. 0) WRITE (LUN,120) NIT, SIG, F
  120 FORMAT (1X,3X,I2,' -- SIG = ',D15.8,', F = ',
     .        D15.8)
      IF (F0*F .LT. 0.D0) THEN
C
C   F0*F < 0.  Update (SNEG,FNEG) to (SG0,F0) so that F
C     and FNEG always have opposite signs.  If SIG is closer
C     to SNEG than SG0 and abs(F) < abs(FNEG), then swap
C     (SNEG,FNEG) with (SG0,F0).
C
        T1 = DMAX
        T2 = FNEG
        DMAX = DSIG
        FNEG = F0
        IF ( ABS(DSIG) .GT. ABS(T1)  .AND.
     .       ABS(F) .LT. ABS(T2)          ) THEN
C
          DSIG = T1
          F0 = T2
        ENDIF
      ENDIF
C
C   Test for convergence.
C
      STOL = RTOL*SIG
      IF (ABS(DMAX) .LE. STOL  .OR.  (F .GE. 0.D0  .AND.
     .    F .LE. FTOL)  .OR.  ABS(F) .LE. RTOL) GO TO 5
      IF (F0*F .LT. 0  .OR.  ABS(F) .LT. ABS(F0)) GO TO 4
C
C   F*F0 > 0 and the new estimate would be outside of the
C     bracketing interval of length abs(DMAX).  Reset
C     (SG0,F0) to (SNEG,FNEG).
C
    3 DSIG = DMAX
      F0 = FNEG
C
C   Compute the change in SIG by linear interpolation be-
C     tween (SG0,F0) and (SIG,F).
C
    4 DSIG = -F*DSIG/(F-F0)
      IF (LUN .GE. 0) WRITE (LUN,130) DSIG
  130 FORMAT (1X,8X,'DSIG = ',D15.8)
      IF ( ABS(DSIG) .GT. ABS(DMAX)  .OR.
     .     DSIG*DMAX .GT. 0. ) GO TO 3
C
C   Restrict the step-size such that abs(DSIG) .GE. STOL/2.
C     Note that DSIG and DMAX have opposite signs.
C
      IF (ABS(DSIG) .LT. STOL/2.D0)
     .  DSIG = -SIGN(STOL/2.D0,DMAX)
C
C   Bottom of loop:  update SIG, DMAX, and F0.
C
      SIG = SIG + DSIG
      DMAX = DMAX + DSIG
      F0 = F
      GO TO 2
C
C No errors encountered and SIGMA finite.
C
    5 IER = 0
      SIG1 = SIG
      RETURN
C
C Infinite tension required.
C
    6 IER = 1
      SIG1 = SBIG
      RETURN
C
C Error in an input parameter.
C
    7 IER = -1
      SIG1 = -1.D0
      RETURN
C
C Invalid constraint.
C
    8 IER = -2
      SIG1 = -1.D0
      RETURN
      END
      DOUBLE PRECISION FUNCTION SIG2 (X1,X2,Y1,Y2,Y1P,Y2P,
     .                                IFL,TOL, IER)
      INTEGER IFL, IER
      DOUBLE PRECISION X1, X2, Y1, Y2, Y1P, Y2P, TOL
C
C***********************************************************
C
C                                                From TSPACK
C                                            Robert J. Renka
C                                  Dept. of Computer Science
C                                       Univ. of North Texas
C                                           renka@cs.unt.edu
C                                                   07/08/92
C
C   Given a pair of abscissae with associated ordinates and
C slopes, this function determines the smallest (nonnega-
C tive) tension factor SIGMA such that the Hermite interpo-
C latory tension spline H(x) preserves convexity (or con-
C cavity) of the data;  i.e.,
C
C   Y1P .LE. S .LE. Y2P implies HPP(x) .GE. 0  or
C   Y1P .GE. S .GE. Y2P implies HPP(x) .LE. 0
C
C for all x in the open interval (X1,X2), where S = (Y2-Y1)/
C (X2-X1) and HPP denotes the second derivative of H.  Note,
C however, that infinite tension is required if Y1P = S or
C Y2P = S (unless Y1P = Y2P = S).
C
C On input:
C
C       X1,X2 = Abscissae.  X1 < X2.
C
C       Y1,Y2 = Values of H at X1 and X2.
C
C       Y1P,Y2P = Derivative values of H at X1 and X2.
C
C       IFL = Option indicator (sign of HPP):
C             IFL = -1 if HPP is to be bounded above by 0.
C             IFL = 1 if HPP is to be bounded below by 0
C                     (preserve convexity of the data).
C
C       TOL = Tolerance whose magnitude determines how close
C             SIGMA is to its optimal value when nonzero
C             finite tension is necessary and sufficient to
C             satisfy convexity or concavity.  In the case
C             of convexity, SIGMA is chosen so that 0 .LE.
C             HPPMIN .LE. abs(TOL), where HPPMIN is the min-
C             imum value of HPP in the interval.  In the
C             case of concavity, the maximum value of HPP
C             satisfies -abs(TOL) .LE. HPPMAX .LE. 0.  Thus,
C             the constraint is satisfied but possibly with
C             more tension than necessary.
C
C Input parameters are not altered by this function.
C
C On output:
C
C       IER = Error indicator:
C             IER = 0 if no errors were encountered and fin-
C                     ite tension is sufficient to satisfy
C                     the constraint.
C             IER = 1 if no errors were encountered but in-
C                     finite tension is required to satisfy
C                     the constraint.
C             IER = -1 if X2 .LE. X1 or abs(IFL) .NE. 1.
C             IER = -2 if the constraint cannot be satis-
C                      fied:  the sign of S-Y1P or Y2P-S
C                      does not agree with IFL.
C
C       SIG2 = Tension factor defined above unless IER < 0,
C              in which case SIG2 = -1.  If IER = 1, SIG2
C              is set to 85, resulting in an approximation
C              to the linear interpolant of the endpoint
C              values.  Note, however, that SIG2 may be
C              larger than 85 if IER = 0.
C
C Modules required by SIG2:  SNHCSH, STORE
C
C Intrinsic functions called by SIG2:  ABS, EXP, MAX, MIN,
C                                        SQRT
C
C***********************************************************
C
      INTEGER LUN, NIT
      DOUBLE PRECISION COSHM, D1, D2, DSIG, DUMMY, DX, EMS,
     .                 F, FP, FTOL, RTOL, S, SBIG, SIG,
     .                 SINHM, SSM, T, T1, TP1
      DOUBLE PRECISION STORE
C
      DATA SBIG/85.D0/, LUN/-1/
C
C Test for an errors in the input parameters.
C
      DX = X2 - X1
      IF (ABS(IFL) .NE. 1.D0  .OR.  DX .LE. 0.D0) GO TO 5
C
C Compute the slope and second differences, and test for
C   an invalid constraint.
C
      S = (Y2-Y1)/DX
      D1 = S - Y1P
      D2 = Y2P - S
      IF ((IFL .GT. 0.D0  .AND.  MIN(D1,D2) .LT. 0.D0)
     .    .OR.  (IFL .LT. 0.D0  .AND.
     .           MAX(D1,D2) .GT. 0.D0)) GO TO 6
C
C Test for infinite tension required.
C
      IF (D1*D2 .EQ. 0.D0  .AND.  D1 .NE. D2) GO TO 4
C
C Test for SIG = 0 sufficient.
C
      SIG = 0.D0
      IF (D1*D2 .EQ. 0.D0) GO TO 3
      T = MAX(D1/D2,D2/D1)
      IF (T .LE. 2.D0) GO TO 3
C
C Find a zero of F(SIG) = SIG*COSHM(SIG)/SINHM(SIG) - (T+1).
C   Since the derivative of F vanishes at the origin, a
C   quadratic approximation is used to obtain an initial
C   estimate for the Newton method.
C
      TP1 = T + 1.D0
      SIG = SQRT(10.D0*T-20.D0)
      NIT = 0
C
C   Compute an absolute tolerance FTOL = abs(TOL) and a
C     relative tolerance RTOL = 100*MACHEPS.
C
      FTOL = ABS(TOL)
      RTOL = 1.D0
    1 RTOL = RTOL/2.D0
        IF (STORE(RTOL+1.D0) .GT. 1.D0) GO TO 1
      RTOL = RTOL*200.D0
C
C Evaluate F and its derivative FP.
C
    2 IF (SIG .LE. .5D0) THEN
C
C   Use approximations designed to avoid cancellation error
C     in the hyperbolic functions.
C
        CALL SNHCSH (SIG, SINHM,COSHM,DUMMY)
        T1 = COSHM/SINHM
        FP = T1 + SIG*(SIG/SINHM - T1*T1 + 1.D0)
      ELSE
C
C   Scale SINHM and COSHM by 2*exp(-SIG) in order to avoid
C     overflow.
C
        EMS = EXP(-SIG)
        SSM = 1.D0 - EMS*(EMS+SIG+SIG)
        T1 = (1.D0-EMS)*(1.D0-EMS)/SSM
        FP = T1 + SIG*(2.D0*SIG*EMS/SSM - T1*T1 + 1.D0)
      ENDIF
C
      F = SIG*T1 - TP1
      IF (LUN .GE. 0) WRITE (LUN,100) SIG, F, FP
  100 FORMAT (1X,'SIG2 -- SIG = ',D15.8,', F(SIG) = ',
     .        D15.8/1X,29X,'FP(SIG) = ',D15.8)
      NIT = NIT + 1
C
C   Test for convergence.
C
      IF (FP .LE. 0.D0) GO TO 3
      DSIG = -F/FP
      IF (ABS(DSIG) .LE. RTOL*SIG  .OR.  (F .GE. 0.D0 .AND.
     .    F .LE. FTOL)  .OR.  ABS(F) .LE. RTOL) GO TO 3
C
C   Update SIG.
C
      SIG = SIG + DSIG
      GO TO 2
C
C No errors encountered, and SIGMA is finite.
C
    3 IER = 0
      SIG2 = SIG
      RETURN
C
C Infinite tension required.
C
    4 IER = 1
      SIG2 = SBIG
      RETURN
C
C X2 .LE. X1 or abs(IFL) .NE. 1.
C
    5 IER = -1
      SIG2 = -1.D0
      RETURN
C
C The constraint cannot be satisfied.
C
    6 IER = -2
      SIG2 = -1.D0
      RETURN
      END
      SUBROUTINE SIGBI (N,X,Y,YP,TOL,B,BMAX, SIGMA, ICFLG,
     .                  DSMAX,IER)
      INTEGER N, ICFLG(N), IER
      DOUBLE PRECISION X(N), Y(N), YP(N), TOL, B(5,N), BMAX,
     .                 SIGMA(N), DSMAX
C
C***********************************************************
C
C                                                From TSPACK
C                                            Robert J. Renka
C                                  Dept. of Computer Science
C                                       Univ. of North Texas
C                                           renka@cs.unt.edu
C                                                   06/10/92
C
C   Given a set of abscissae X with associated data values Y
C and derivatives YP, this subroutine determines the small-
C est (nonnegative) tension factors SIGMA such that the Her-
C mite interpolatory tension spline H(x) satisfies a set of
C user-specified constraints.
C
C   SIGBI may be used in conjunction with Subroutine YPC2
C (or YPC2P) in order to produce a C-2 interpolant which
C satisfies the constraints.  This is achieved by calling
C YPC2 with SIGMA initialized to the zero vector, and then
C alternating calls to SIGBI with calls to YPC2 until the
C change in SIGMA is small (refer to the parameter descrip-
C tions for SIGMA, DSMAX and IER), or the maximum relative
C change in YP is bounded by a tolerance (a reasonable value
C is .01).  A similar procedure may be used to produce a C-2
C shape-preserving smoothing curve (Subroutine SMCRV).
C
C   Refer to Subroutine SIGS for a means of selecting mini-
C mum tension factors to preserve shape properties of the
C data.
C
C On input:
C
C       N = Number of data points.  N .GE. 2.
C
C       X = Array of length N containing a strictly in-
C           creasing sequence of abscissae:  X(I) < X(I+1)
C           for I = 1,...,N-1.
C
C       Y = Array of length N containing data values (or
C           function values computed by SMCRV) associated
C           with the abscissae.  H(X(I)) = Y(I) for I =
C           1,...,N.
C
C       YP = Array of length N containing first derivatives
C            of H at the abscissae.  Refer to Subroutines
C            YPC1, YPC1P, YPC2, YPC2P, and SMCRV.
C
C       TOL = Tolerance whose magnitude determines how close
C             each tension factor is to its optimal value
C             when nonzero finite tension is necessary and
C             sufficient to satisfy a constraint.  Refer to
C             functions SIG0, SIG1, and SIG2.  TOL should be
C             set to 0 for optimal tension.
C
C       B = Array dimensioned 5 by N-1 containing bounds or
C           flags which define the constraints.  For I = 1
C           to N-1, column I defines the constraints associ-
C           ated with interval I (X(I),X(I+1)) as follows:
C
C             B(1,I) is an upper bound on H
C             B(2,I) is a lower bound on H
C             B(3,I) is an upper bound on HP
C             B(4,I) is a lower bound on HP
C             B(5,I) specifies the required sign of HPP
C
C           where HP and HPP denote the first and second
C           derivatives of H, respectively.  A null con-
C           straint is specified by abs(B(K,I)) .GE. BMAX
C           for K < 5, or B(5,I) = 0:   B(1,I) .GE. BMAX,
C           B(2,I) .LE. -BMAX, B(3,I) .GE. BMAX, B(4,I) .LE.
C           -BMAX, or B(5,I) = 0.  Any positive value of
C           B(5,I) specifies that H should be convex, a
C           negative values specifies that H should be con-
C           cave, and 0 specifies that no restriction be
C           placed on HPP.  Refer to Functions SIG0, SIG1,
C           and SIG2 for definitions of valid constraints.
C
C       BMAX = User-defined value of infinity which, when
C              used as an upper bound in B (or when its
C              negative is used as a lower bound), specifies
C              that no constraint is to be enforced.
C
C The above parameters are not altered by this routine.
C
C       SIGMA = Array of length N-1 containing minimum val-
C               ues of the tension factors.  SIGMA(I) is as-
C               sociated with interval (I,I+1) and SIGMA(I)
C               .GE. 0 for I = 1,...,N-1.  SIGMA should be
C               set to the zero vector if minimal tension
C               is desired, and should be unchanged from a
C               previous call in order to ensure convergence
C               of the C-2 iterative procedure.
C
C       ICFLG = Array of length .GE. N-1.
C
C On output:
C
C       SIGMA = Array containing tension factors for which
C               H(x) satisfies the constraints defined by B,
C               with the restriction that SIGMA(I) .LE. 85
C               for all I (unless the input value is larger).
C               The factors are as small as possible (within
C               the tolerance), but not less than their
C               input values.  If infinite tension is re-
C               quired in interval (X(I),X(I+1)), then
C               SIGMA(I) = 85 (and H is an approximation to
C               the linear interpolant on the interval),
C               and if no constraint is specified in the
C               interval, then SIGMA(I) = 0 (unless the
C               input value is positive), and thus H is
C               cubic.  Invalid constraints are treated as
C               null constraints.
C
C       ICFLG = Array of invalid constraint flags associated
C               with intervals.  For I = 1 to N-1, ICFLG(I)
C               is a 5-bit value b5b4b3b2b1, where bK = 1 if
C               and only if constraint K cannot be satis-
C               fied.  Thus, all constraints in interval I
C               are satisfied if and only if ICFLG(I) = 0
C               (and IER .GE. 0).
C
C       DSMAX = Maximum increase in a component of SIGMA
C               from its input value.  The increase is a
C               relative change if the input value is
C               positive, and an absolute change otherwise.
C
C       IER = Error indicator and information flag:
C             IER = I if no errors (other than invalid con-
C                     straints) were encountered and I
C                     components of SIGMA were altered from
C                     their input values for 0 .LE. I .LE.
C                     N-1.
C             IER = -1 if N < 2.  SIGMA and ICFLG are not
C                      altered in this case.
C             IER = -I if X(I) .LE. X(I-1) for some I in the
C                      range 2,...,N.  SIGMA(J) and ICFLG(J)
C                      are unaltered for J .GE. I-1 in this
C                      case.
C
C Modules required by SIGBI:  SIG0, SIG1, SIG2, SNHCSH,
C                               STORE
C
C Intrinsic functions called by SIGBI:  ABS, MAX, MIN
C
C***********************************************************
C
      INTEGER I, ICFK, ICNT, IERR, IFL, K, NM1
      DOUBLE PRECISION BMX, BND, DSIG, DSM, S, SBIG, SIG,
     .                 SIGIN
      DOUBLE PRECISION SIG0, SIG1, SIG2
C
      DATA SBIG/85.D0/
      NM1 = N - 1
      IF (NM1 .LT. 1) GO TO 4
      BMX = BMAX
C
C Initialize change counter ICNT and maximum change DSM for
C   loop on intervals.
C
      ICNT = 0
      DSM = 0.D0
      DO 3 I = 1,NM1
        IF (X(I) .GE. X(I+1)) GO TO 5
        ICFLG(I) = 0
C
C Loop on constraints for interval I.  SIG is set to the
C   largest tension factor required to satisfy all five
C   constraints.  ICFK = 2**(K-1) is the increment for
C   ICFLG(I) when constraint K is invalid.
C
        SIG = 0.D0
        ICFK = 1
        DO 2 K = 1,5
          BND = B(K,I)
          IF (K .LT. 5  .AND.  ABS(BND) .GE. BMX) GO TO 1
          IF (K .LE. 2) THEN
            IFL = 3 - 2*K
            S = SIG0 (X(I),X(I+1),Y(I),Y(I+1),YP(I),YP(I+1),
     .                IFL,BND,TOL, IERR)
          ELSEIF (K .LE. 4) THEN
            IFL = 7 - 2*K
            S = SIG1 (X(I),X(I+1),Y(I),Y(I+1),YP(I),YP(I+1),
     .                IFL,BND,TOL, IERR)
          ELSE
            IF (BND .EQ. 0.D0) GO TO 1
            IFL = -1
            IF (BND .GT. 0.D0) IFL = 1
            S = SIG2 (X(I),X(I+1),Y(I),Y(I+1),YP(I),YP(I+1),
     .                IFL,TOL, IERR)
          ENDIF
          IF (IERR .EQ. -2) THEN
C
C   An invalid constraint was encountered.  Increment
C     ICFLG(I).
C
            ICFLG(I) = ICFLG(I) + ICFK
          ELSE
C
C   Update SIG.
C
            SIG = MAX(SIG,S)
          ENDIF
C
C   Bottom of loop on constraints K:  update ICFK.
C
    1     ICFK = 2*ICFK
    2     CONTINUE
C
C Bottom of loop on intervals:  update SIGMA(I), ICNT, and
C   DSM if necessary.
C
        SIG = MIN(SIG,SBIG)
        SIGIN = SIGMA(I)
        IF (SIG .GT. SIGIN) THEN
          SIGMA(I) = SIG
          ICNT = ICNT + 1
          DSIG = SIG-SIGIN
          IF (SIGIN .GT. 0.D0) DSIG = DSIG/SIGIN
          DSM = MAX(DSM,DSIG)
        ENDIF
    3   CONTINUE
C
C No errors (other than invalid constraints) encountered.
C
      DSMAX = DSM
      IER = ICNT
      RETURN
C
C N < 2.
C
    4 DSMAX = 0.D0
      IER = -1
      RETURN
C
C X(I) .GE. X(I+1).
C
    5 DSMAX = DSM
      IER = -(I+1)
      RETURN
      END
      SUBROUTINE SIGBP (N,X,Y,XP,YP,TOL,BL,BU,
     .                  BMAX, SIGMA, DSMAX,IER)
      INTEGER N, IER
      DOUBLE PRECISION X(N), Y(N), XP(N), YP(N), TOL, BL(N),
     .                 BU(N), BMAX, SIGMA(N), DSMAX
C
C***********************************************************
C
C                                                From TSPACK
C                                            Robert J. Renka
C                                  Dept. of Computer Science
C                                       Univ. of North Texas
C                                           renka@cs.unt.edu
C                                                   11/18/96
C
C   Given an ordered sequence of points C(I) = (X(I),Y(I))
C with associated derivative vectors CP(I) = (XP(I),YP(I)),
C this subroutine determines the smallest (nonnegative) ten-
C sion factors SIGMA such that a parametric planar curve
C C(t) satisfies a set of user-specified constraints.  The
C components x(t) and y(t) of C(t) are the Hermite interpo-
C latory tension splines defined by the data and tension
C factors:  C(t(I)) = C(I) and C'(t(I)) = CP(I) for para-
C meter values t(1), t(2), ..., t(N).  In each subinterval
C [t1,t2], the signed perpendicular distance from the
C corresponding line segment C1-C2 to the curve C(t) is
C given by the vector cross product
C
C     d(t) = (C2-C1)/DC X (C(t)-C1)
C
C where DC = abs(C2-C1) is the length of the line segment.
C The associated tension factor SIGMA is chosen to satisfy
C an upper bound on the maximum of d(t) and a lower bound on
C the minimum of d(t) over t in [t1,t2].  Thus, the upper
C bound is associated with distance to the left of the line
C segment as viewed from C1 toward C2.  Note that the curve
C is assumed to be parameterized by arc length (Subroutine
C ARCL2D) so that t2-t1 = DC.  If this is not the case, the
C required bounds should be scaled by DC/(t2-t1) to obtain
C the input parameters BL and BU.
C
C   SIGBP may be used in conjunction with Subroutine YPC2
C (or YPC2P) in order to produce a C-2 interpolant which
C satisfies the constraints.  This is achieved by calling
C YPC2 with SIGMA initialized to the zero vector, and then
C alternating calls to SIGBP with calls to YPC2 until the
C change in SIGMA is small (refer to the parameter descrip-
C tions for SIGMA, DSMAX and IER), or the maximum relative
C change in YP is bounded by a tolerance (a reasonable value
C is .01).  A similar procedure may be used to produce a C-2
C shape-preserving smoothing curve (Subroutine SMCRV).
C
C On input:
C
C       N = Number of data points.  N .GE. 2.
C
C       X,Y = Arrays of length N containing the Cartesian
C             coordinates of the points C(I), I = 1 to N.
C
C       XP,YP = Arrays of length N containing the components
C               of the derivative (velocity) vectors CP(I).
C               Refer to Subroutines YPC1, YPC1P, YPC2,
C               YPC2P, and SMCRV.
C
C       TOL = Tolerance whose magnitude determines how close
C             each tension factor SIGMA is to its optimal
C             value when nonzero finite tension is necessary
C             and sufficient to satisfy a constraint.
C             SIGMA(I) is chosen so that BL(I) .LE. dmin
C             .LE. BL(I) + abs(TOL) and BU(I) - abs(TOL)
C             .LE. dmax .LE. BU(I), where dmin and dmax are
C             the minimum and maximum values of d(t) in the
C             interval [t(I),t(I+1)].  Thus, a large toler-
C             ance might increase execution efficiency but
C             may result in more tension than is necessary.
C             TOL may be set to 0 for optimal tension.
C
C       BL,BU = Arrays of length N-1 containing lower and
C               upper bounds, respectively, which define
C               the constraints as described above.  BL(I)
C               < 0 and BU(I) > 0 for I = 1 to N-1.  A null
C               straint is specified by BL(I) .LE. -BMAX or
C               BU(I) .GE. BMAX.
C
C       BMAX = User-defined value of infinity which, when
C              used as an upper bound in BU (or when its
C              negative is used as a lower bound in BL),
C              specifies that no constraint is to be en-
C              forced.
C
C The above parameters are not altered by this routine.
C
C       SIGMA = Array of length N-1 containing minimum val-
C               ues of the tension factors.  SIGMA(I) is as-
C               sociated with interval (I,I+1) and SIGMA(I)
C               .GE. 0 for I = 1,...,N-1.  SIGMA should be
C               set to the zero vector if minimal tension
C               is desired, and should be unchanged from a
C               previous call in order to ensure convergence
C               of the C-2 iterative procedure.
C
C On output:
C
C       SIGMA = Array containing tension factors for which
C               d(t) satisfies the constraints defined by
C               BL and BU, with the restriction that
C               SIGMA(I) .LE. 85 for all I (unless the input
C               value is larger).  The factors are as small
C               as possible (within the tolerance), but not
C               less than their input values.  If no con-
C               straint is specified in interval I, then
C               SIGMA(I) = 0 (unless the input value is
C               positive), and thus x(t) and y(t) are cubic
C               polynomials.
C
C       DSMAX = Maximum increase in a component of SIGMA
C               from its input value.  The increase is a
C               relative change if the input value is
C               positive, and an absolute change otherwise.
C
C       IER = Error indicator and information flag:
C             IER = I if no errors were encountered and I
C                     components of SIGMA were altered from
C                     their input values for 0 .LE. I .LE.
C                     N-1.
C             IER = -1 if N < 2.  SIGMA is not altered in
C                      this case.
C             IER = -I if BL(I-1) .GE. 0 or BU(I-1) .LE. 0
C                      for some I in the range 2 to N.
C                      SIGMA(J) is unaltered for J .GE. I-1
C                      in this case.
C
C Modules required by SIGBP:  SNHCSH, STORE
C
C Intrinsic functions called by SIGBP:  ABS, EXP, LOG, MAX,
C                                         MIN, SIGN, SQRT
C
C***********************************************************
C
      INTEGER I, ICNT, IP1, LUN, NIT, NM1
      DOUBLE PRECISION A, A1, A2, AA, B, B0, BHI, BLO, BMX,
     .                 C, COSHM, COSHMM, D, D0, DM, DMAX,
     .                 DP, DSIG, DSM, DX, DY, E, EB, EMS, F,
     .                 F0, FMAX, FNEG, FTOL, RM, RP, RSM,
     .                 RSP, RTOL, S, SBIG, SIG, SIGIN, SINH,
     .                 SINHM, SNEG, STOL, T, T1, T2, TM, V1,
     .                 V2, V2M1
      DOUBLE PRECISION STORE
C
      DATA SBIG/85.D0/,  LUN/-1/
      NM1 = N - 1
      IF (NM1 .LT. 1) GO TO 8
      BMX = BMAX
C
C Compute an absolute tolerance FTOL = abs(TOL), and a
C   relative tolerance RTOL = 100*MACHEPS.
C
      FTOL = ABS(TOL)
      RTOL = 1.D0
    1 RTOL = RTOL/2.D0
        IF (STORE(RTOL+1.D0) .GT. 1.D0) GO TO 1
      RTOL = RTOL*200.D0
C
C Initialize change counter ICNT and maximum change DSM for
C   loop on intervals.
C
      ICNT = 0
      DSM = 0.D0
      DO 7 I = 1,NM1
        IP1 = I + 1
        BLO = BL(I)
        BHI = BU(I)
        SIGIN = SIGMA(I)
        IF (LUN .GE. 0) WRITE (LUN,100) I, BLO, BHI, SIGIN
  100   FORMAT (//1X,'SIGBP -- INTERVAL',I4,', BL = ',D10.3,
     .          ', BU = ',D10.3,', SIGIN = ',D15.8)
        IF (BLO .GE. 0.D0  .OR.  BHI .LE. 0.D0) GO TO 9
        IF (SIGIN .GE. SBIG) GO TO 7
C
C Initialize SIG to 0 and test for a null constraint.
C
        SIG = 0.D0
        IF (BLO .LE. -BMX  .AND.  BHI .GE. BMX) GO TO 6
C
C Test for SIG = 0 sufficient.
C
C   The signed orthogonal distance is d(b) = b*(1-b)*
C     (b*V1 - (1-b)*V2), where b = (t2-t)/(t2-t1),
C     V1 = (C2-C1) X CP(1), and V2 = (C2-C1) X CP(2).
C
        DX = X(IP1) - X(I)
        DY = Y(IP1) - Y(I)
        V1 = DX*YP(I) - DY*XP(I)
        V2 = DX*YP(IP1) - DY*XP(IP1)
C
C   Set DP and DM to the maximum and minimum values of d(b)
C     for b in [0,1].  Note that DP .GE. 0 and DM .LE. 0.
C
        S = V1 + V2
        IF (S .EQ. 0.D0) THEN
C
C   The derivative d'(b) is zero at the midpoint b = .5.
C
          IF (V1 .GE. 0.D0) THEN
            DP = V1/4.D0
            DM = 0.D0
          ELSE
            DP = 0.D0
            DM = V1/4.D0
          ENDIF
        ELSE
C
C   Set RP/RM to the roots of the quadratic equation d'(b) =
C     (B0 +/- SQRT(D0))/(3*S) = V2/(B0 -/+ SQRT(D0)) = 0,
C     where B0 = V1 + 2*V2 and D0 = V1**2 + V1*V2 + V2**2.
C     The expression is chosen to avoid cancellation error.
C
          B0 = S + V2
          D0 = S*S - V1*V2
          T = B0 + SIGN(SQRT(D0),B0)
          IF (B0 .GE. 0.D0) THEN
            RP = T/(3.D0*S)
            RM = V2/T
          ELSE
            RP = V2/T
            RM = T/(3.D0*S)
          ENDIF
          IF (V1 .LE. 0.D0  .AND.  V2 .GE. 0.D0) THEN
C
C   The maximum is DP = 0 at the endpoints.
C
            DP = 0.D0
          ELSE
            DP = RP*(1.D0-RP)*(RP*S - V2)
          ENDIF
          IF (V1 .GE. 0.D0  .AND.  V2 .LE. 0.D0) THEN
C
C   The minimum is DM = 0 at the endpoints.
C
            DM = 0.D0
          ELSE
            DM = RM*(1.D0-RM)*(RM*S - V2)
          ENDIF
        ENDIF
C
C   SIG = 0 is sufficient to satisfy the constraints iff
C     DP .LE. BHI and DM .GE. BLO iff F0 .GE. 0.
C
        F0 = MIN(BHI-DP,DM-BLO)
        IF (F0 .GE. 0.D0) GO TO 6
C
C Find a zero of F(SIG) = min(BHI-DP,DM-BLO), where DP and
C   DM are the maximum and minimum values of d(b).  F is an
C   increasing function, F(0) = F0 < 0, and F = FMAX =
C   min(BHI,-BLO) for SIG sufficiently large.  Note that F
C   has a discontinuity in its first derivative if the
C   curves BHI-DP and DM-BLO (as functions of SIG) inter-
C   sect, and the rate of convergence of the zero finder is
C   reduced to linear if such an intersection occurs near
C   the zero of F.
C
C Initialize parameters for the secant method.  The method
C   uses three points:  (SG0,F0), (SIG,F), and (SNEG,FNEG),
C   where SG0 and SNEG are defined implicitly by DSIG = SIG
C   - SG0 and DMAX = SIG - SNEG.  SG0 is initially zero and
C   SNEG is initialized to a sufficiently large value that
C   FNEG > 0.  This value is used only if the initial value
C   of F is negative.
C
        T = MIN(BHI,-BLO)
        FMAX = MAX(1.D-3,T)
        SIG = MAX(ABS(V1),ABS(V2))/T
        DMAX = SIG*(1.D0-T/FMAX)
        SNEG = SIG - DMAX
        IF (LUN .GE. 0) WRITE (LUN,110) F0, FMAX, SNEG
  110   FORMAT (//1X,'F(0) = ',D15.8,', FMAX = ',D15.8,
     .          ', SNEG = ',D15.8/)
        DSIG = SIG
        FNEG = FMAX
        V2M1 = V2 - V1
        NIT = 0
C
C Top of loop:  compute F.
C
    2   EMS = EXP(-SIG)
        IF (SIG .LE. .5D0) THEN
C
C   SIG .LE. .5:  use approximations designed to avoid can-
C                 cellation error (associated with small
C                 SIG) in the modified hyperbolic functions.
C
          CALL SNHCSH (SIG, SINHM,COSHM,COSHMM)
          SINH = SINHM + SIG
          A1 = SIG*COSHM*V2 - SINHM*V2M1
          A2 = SIG*SINH*V2 - COSHM*V2M1
          A = A2 - A1
          AA = A/EMS
          E = SIG*SINHM - COSHMM - COSHMM
        ELSE
C
C   SIG > .5:  scale SINHM and COSHM by 2*EXP(-SIG) in order
C              to avoid overflow.
C
          TM = 1.D0 - EMS
          SINH = TM*(1.D0+EMS)
          SINHM = SINH - 2.D0*SIG*EMS
          COSHM = TM*TM
          A1 = SIG*COSHM*V2 - SINHM*V2M1
          A2 = SIG*SINH*V2 - COSHM*V2M1
          AA = 2.D0*(SIG*TM*V2 + (TM-SIG)*V2M1)
          A = EMS*AA
          E = SIG*SINH - COSHM - COSHM
        ENDIF
        IF (S .EQ. 0.D0) THEN
C
C   The derivative d'(b) is zero at the midpoint b = .5.
C
          EB = SIG*COSHM - SINHM - SINHM
          IF (V1 .GE. 0.D0) THEN
            DP = E*V1/(SIG*(SQRT(EB*EB-E*E)+EB))
            DM = 0.D0
          ELSE
            DP = 0.D0
            DM = E*V1/(SIG*(SQRT(EB*EB-E*E)+EB))
          ENDIF
        ELSE
C
C   d'(b)*DC = V2 - (A1*sinh(SIG*b) - A2*coshm(SIG*b))/E = 0
C     for ESB = (-B +/- sqrt(D))/A = C/(-B -/+ sqrt(D)),
C     where ESB = exp(SIG*b), A = A2-A1, D = B**2 - A*C, and
C     B and C are defined below.
C
          B = -COSHM*S
          C = A2 + A1
          D = B*B - A*C
          F = FMAX
          IF (D .LT. 0.D0) GO TO 3
          T1 = SQRT(D)
          T = -B - SIGN(T1,B)
C
          RSP = 0.D0
          IF (B .LT. 0.D0  .AND.  AA .NE. 0.D0) THEN
            IF (T/AA .GT. 0.D0) RSP = SIG + LOG(T/AA)
          ENDIF
          IF ((B .GT. 0.D0  .OR.  AA .EQ. 0.D0)  .AND.
     .        C/T .GT. 0.D0) RSP = LOG(C/T)
          IF ((RSP .LE. 0.D0  .OR.  RSP .GE. SIG)  .AND.
     .        B .NE. 0.D0) THEN
C
C   The maximum is DP = 0 at the endpoints.
C
            DP = 0.D0
          ELSE
            DP = -(B*RSP+A1+T1)/(SIG*E)
          ENDIF
C
          RSM = 0.D0
          IF (B .GT. 0.D0  .AND.  AA .NE. 0.D0) THEN
            IF (T/AA .GT. 0.D0) RSM = SIG + LOG(T/AA)
          ENDIF
          IF ((B .LT. 0.D0  .OR.  AA .EQ. 0.D0)  .AND.
     .        C/T .GT. 0.D0) RSM = LOG(C/T)
          IF ((RSM .LE. 0.D0  .OR.  RSM .GE. SIG)  .AND.
     .        B .NE. 0.D0) THEN
C
C   The minimum is DM = 0 at the endpoints.
C
            DM = 0.D0
          ELSE
            DM = -(B*RSM+A1-T1)/(SIG*E)
          ENDIF
        ENDIF
C
        F = MIN(BHI-DP,DM-BLO)
C
C   Update the number of iterations NIT.
C
    3   NIT = NIT + 1
        IF (LUN .GE. 0) WRITE (LUN,120) NIT, SIG, F
  120   FORMAT (1X,3X,I2,' -- SIG = ',D15.8,', F = ',
     .          D15.8)
        IF (F0*F .LT. 0.D0) THEN
C
C   F0*F < 0.  Update (SNEG,FNEG) to (SG0,F0) so that F
C     and FNEG always have opposite signs.  If SIG is
C     closer to SNEG than SG0, then swap (SNEG,FNEG) with
C     (SG0,F0).
C
          T1 = DMAX
          T2 = FNEG
          DMAX = DSIG
          FNEG = F0
          IF (ABS(DSIG) .GT. ABS(T1)) THEN
            DSIG = T1
            F0 = T2
          ENDIF
        ENDIF
C
C   Test for convergence.
C
        STOL = RTOL*SIG
        IF (ABS(DMAX) .LE. STOL  .OR.  (F .GE. 0.D0  .AND.
     .      F .LE. FTOL)  .OR.  ABS(F) .LE. RTOL) GO TO 6
C
C   Test for F0 = F = FMAX or F < 0 on the first iteration.
C
        IF (F0 .NE. F  .AND.  (NIT .GT. 1  .OR.
     .                         F .GT. 0.D0))    GO TO 5
C
C   F*F0 > 0 and either the new estimate would be outside of
C     the bracketing interval of length abs(DMAX) or F < 0
C     on the first iteration.  Reset (SG0,F0) to (SNEG,FNEG).
C
    4   DSIG = DMAX
        F0 = FNEG
C
C   Compute the change in SIG by linear interpolation be-
C     tween (SG0,F0) and (SIG,F).
C
    5   DSIG = -F*DSIG/(F-F0)
        IF (LUN .GE. 0) WRITE (LUN,130) DSIG
  130   FORMAT (1X,8X,'DSIG = ',D15.8)
        IF ( ABS(DSIG) .GT. ABS(DMAX)  .OR.
     .       DSIG*DMAX .GT. 0. ) GO TO 4
C
C   Restrict the step-size such that abs(DSIG) .GE. STOL/2.
C     Note that DSIG and DMAX have opposite signs.
C
        IF (ABS(DSIG) .LT. STOL/2.D0)
     .    DSIG = -SIGN(STOL/2.D0,DMAX)
C
C   Bottom of loop:  update SIG, DMAX, and F0.
C
        SIG = SIG + DSIG
        DMAX = DMAX + DSIG
        F0 = F
        GO TO 2
C
C Bottom of loop on intervals:  update SIGMA(I), ICNT, and
C   DSM if necessary.
C
    6   SIG = MIN(SIG,SBIG)
        IF (SIG .GT. SIGIN) THEN
          SIGMA(I) = SIG
          ICNT = ICNT + 1
          DSIG = SIG-SIGIN
          IF (SIGIN .GT. 0.D0) DSIG = DSIG/SIGIN
          DSM = MAX(DSM,DSIG)
        ENDIF
    7   CONTINUE
C
C No errors encountered.
C
      DSMAX = DSM
      IER = ICNT
      RETURN
C
C N < 2.
C
    8 DSMAX = 0.D0
      IER = -1
      RETURN
C
C BL(I) .GE. 0 or BU(I) .LE. 0.
C
    9 DSMAX = DSM
      IER = -(IP1)
      RETURN
      END
      SUBROUTINE SIGS (N,X,Y,YP,TOL, SIGMA, DSMAX,IER)
      INTEGER N, IER
      DOUBLE PRECISION X(N), Y(N), YP(N), TOL, SIGMA(N),
     .                 DSMAX
C
C***********************************************************
C
C                                                From TSPACK
C                                            Robert J. Renka
C                                  Dept. of Computer Science
C                                       Univ. of North Texas
C                                           renka@cs.unt.edu
C                                                   11/17/96
C
C   Given a set of abscissae X with associated data values Y
C and derivatives YP, this subroutine determines the small-
C est (nonnegative) tension factors SIGMA such that the Her-
C mite interpolatory tension spline H(x) preserves local
C shape properties of the data.  In an interval (X1,X2) with
C data values Y1,Y2 and derivatives YP1,YP2, the properties
C of the data are
C
C       Monotonicity:  S, YP1, and YP2 are nonnegative or
C                        nonpositive,
C  and
C       Convexity:     YP1 .LE. S .LE. YP2  or  YP1 .GE. S
C                        .GE. YP2,
C
C where S = (Y2-Y1)/(X2-X1).  The corresponding properties
C of H are constant sign of the first and second deriva-
C tives, respectively.  Note that, unless YP1 = S = YP2, in-
C finite tension is required (and H is linear on the inter-
C val) if S = 0 in the case of monotonicity, or if YP1 = S
C or YP2 = S in the case of convexity.
C
C   SIGS may be used in conjunction with Subroutine YPC2
C (or YPC2P) in order to produce a C-2 interpolant which
C preserves the shape properties of the data.  This is
C achieved by calling YPC2 with SIGMA initialized to the
C zero vector, and then alternating calls to SIGS with
C calls to YPC2 until the change in SIGMA is small (refer to
C the parameter descriptions for SIGMA, DSMAX and IER), or
C the maximum relative change in YP is bounded by a toler-
C ance (a reasonable value is .01).  A similar procedure may
C be used to produce a C-2 shape-preserving smoothing curve
C (Subroutine SMCRV).
C
C   Refer to Subroutine SIGBI for a means of selecting mini-
C mum tension factors to satisfy more general constraints.
C
C On input:
C
C       N = Number of data points.  N .GE. 2.
C
C       X = Array of length N containing a strictly in-
C           creasing sequence of abscissae:  X(I) < X(I+1)
C           for I = 1,...,N-1.
C
C       Y = Array of length N containing data values (or
C           function values computed by SMCRV) associated
C           with the abscissae.  H(X(I)) = Y(I) for I =
C           1,...,N.
C
C       YP = Array of length N containing first derivatives
C            of H at the abscissae.  Refer to Subroutines
C            YPC1, YPC1P, YPC2, YPC2P, and SMCRV.
C
C       TOL = Tolerance whose magnitude determines how close
C             each tension factor is to its optimal value
C             when nonzero finite tension is necessary and
C             sufficient to satisfy the constraint:
C             abs(TOL) is an upper bound on the magnitude
C             of the smallest (nonnegative) or largest (non-
C             positive) value of the first or second deriva-
C             tive of H in the interval.  Thus, the con-
C             straint is satisfied, but possibly with more
C             tension than necessary.  TOL should be set to
C             0 for optimal tension.
C
C The above parameters are not altered by this routine.
C
C       SIGMA = Array of length N-1 containing minimum val-
C               ues of the tension factors.  SIGMA(I) is as-
C               sociated with interval (I,I+1) and SIGMA(I)
C               .GE. 0 for I = 1,...,N-1.  SIGMA should be
C               set to the zero vector if minimal tension
C               is desired, and should be unchanged from a
C               previous call in order to ensure convergence
C               of the C-2 iterative procedure.
C
C On output:
C
C       SIGMA = Array containing tension factors for which
C               H(x) preserves the properties of the data,
C               with the restriction that SIGMA(I) .LE. 85
C               for all I (unless the input value is larger).
C               The factors are as small as possible (within
C               the tolerance), but not less than their
C               input values.  If infinite tension is re-
C               quired in interval (X(I),X(I+1)), then
C               SIGMA(I) = 85 (and H is an approximation to
C               the linear interpolant on the interval),
C               and if neither property is satisfied by the
C               data, then SIGMA(I) = 0 (unless the input
C               value is positive), and thus H is cubic in
C               the interval.
C
C       DSMAX = Maximum increase in a component of SIGMA
C               from its input value.  The increase is a
C               relative change if the input value is
C               nonzero, and an absolute change otherwise.
C
C       IER = Error indicator and information flag:
C             IER = I if no errors were encountered and I
C                     components of SIGMA were altered from
C                     their input values for 0 .LE. I .LE.
C                     N-1.
C             IER = -1 if N < 2.  SIGMA is not altered in
C                      this case.
C             IER = -I if X(I) .LE. X(I-1) for some I in the
C                      range 2,...,N.  SIGMA(J-1) is unal-
C                      tered for J = I,...,N in this case.
C
C Modules required by SIGS:  SNHCSH, STORE
C
C Intrinsic functions called by SIGS:  ABS, EXP, MAX, MIN,
C                                        SIGN, SQRT
C
C***********************************************************
C
      INTEGER I, ICNT, IP1, LUN, NIT, NM1
      DOUBLE PRECISION A, C1, C2, COSHM, COSHMM, D0, D1,
     .                 D1D2, D1PD2, D2, DMAX, DSIG, DSM, DX,
     .                 E, EMS, EMS2, F, F0, FMAX, FNEG, FP,
     .                 FTOL, RTOL, S, S1, S2, SBIG, SCM,
     .                 SGN, SIG, SIGIN, SINHM, SSINH, SSM,
     .                 STOL, T, T0, T1, T2, TM, TP1
      DOUBLE PRECISION STORE
C
      DATA SBIG/85.D0/,  LUN/-1/
      NM1 = N - 1
      IF (NM1 .LT. 1) GO TO 9
C
C Compute an absolute tolerance FTOL = abs(TOL) and a
C   relative tolerance RTOL = 100*MACHEPS.
C
      FTOL = ABS(TOL)
      RTOL = 1.D0
    1 RTOL = RTOL/2.D0
        IF (STORE(RTOL+1.D0) .GT. 1.D0) GO TO 1
      RTOL = RTOL*200.D0
C
C Initialize change counter ICNT and maximum change DSM for
C   loop on intervals.
C
      ICNT = 0
      DSM = 0.D0
      DO 8 I = 1,NM1
        IF (LUN .GE. 0) WRITE (LUN,100) I
  100   FORMAT (//1X,'SIGS -- INTERVAL',I4)
        IP1 = I + 1
        DX = X(IP1) - X(I)
        IF (DX .LE. 0.D0) GO TO 10
        SIGIN = SIGMA(I)
        IF (SIGIN .GE. SBIG) GO TO 8
C
C Compute first and second differences.
C
        S1 = YP(I)
        S2 = YP(IP1)
        S = (Y(IP1)-Y(I))/DX
        D1 = S - S1
        D2 = S2 - S
        D1D2 = D1*D2
C
C Test for infinite tension required to satisfy either
C   property.
C
        SIG = SBIG
        IF ((D1D2 .EQ. 0.D0  .AND.  S1 .NE. S2)  .OR.
     .      (S .EQ. 0.D0  .AND.  S1*S2 .GT. 0.D0)) GO TO 7
C
C Test for SIGMA = 0 sufficient.  The data satisfies convex-
C   ity iff D1D2 .GE. 0, and D1D2 = 0 implies S1 = S = S2.
C
        SIG = 0.D0
        IF (D1D2 .LT. 0.D0) GO TO 3
        IF (D1D2 .EQ. 0.D0) GO TO 7
        T = MAX(D1/D2,D2/D1)
        IF (T .LE. 2.D0) GO TO 7
        TP1 = T + 1.D0
C
C Convexity:  Find a zero of F(SIG) = SIG*COSHM(SIG)/
C   SINHM(SIG) - TP1.
C
C   F(0) = 2-T < 0, F(TP1) .GE. 0, the derivative of F
C     vanishes at SIG = 0, and the second derivative of F is
C     .2 at SIG = 0.  A quadratic approximation is used to
C     obtain a starting point for the Newton method.
C
        SIG = SQRT(10.D0*T-20.D0)
        NIT = 0
C
C   Top of loop:
C
    2   IF (SIG .LE. .5D0) THEN
          CALL SNHCSH (SIG, SINHM,COSHM,COSHMM)
          T1 = COSHM/SINHM
          FP = T1 + SIG*(SIG/SINHM - T1*T1 + 1.D0)
        ELSE
C
C   Scale SINHM and COSHM by 2*EXP(-SIG) in order to avoid
C     overflow with large SIG.
C
          EMS = EXP(-SIG)
          SSM = 1.D0 - EMS*(EMS+SIG+SIG)
          T1 = (1.D0-EMS)*(1.D0-EMS)/SSM
          FP = T1 + SIG*(2.D0*SIG*EMS/SSM - T1*T1 + 1.D0)
        ENDIF
C
        F = SIG*T1 - TP1
        IF (LUN .GE. 0) WRITE (LUN,110) SIG, F, FP
  110   FORMAT (5X,'CONVEXITY -- SIG = ',D15.8,
     .          ', F(SIG) = ',D15.8/1X,35X,'FP(SIG) = ',
     .          D15.8)
        NIT = NIT + 1
C
C   Test for convergence.
C
        IF (FP .LE. 0.D0) GO TO 7
        DSIG = -F/FP
        IF (ABS(DSIG) .LE. RTOL*SIG  .OR.  (F .GE. 0.D0
     .      .AND.  F .LE. FTOL)  .OR.  ABS(F) .LE. RTOL)
     .    GO TO 7
C
C   Update SIG.
C
        SIG = SIG + DSIG
        GO TO 2
C
C Convexity cannot be satisfied.  Monotonicity can be satis-
C   fied iff S1*S .GE. 0 and S2*S .GE. 0 since S .NE. 0.
C
    3   IF (S1*S .LT. 0.D0  .OR.  S2*S .LT. 0.D0) GO TO 7
        T0 = 3.D0*S - S1 - S2
        D0 = T0*T0 - S1*S2
C
C SIGMA = 0 is sufficient for monotonicity iff S*T0 .GE. 0
C   or D0 .LE. 0.
C
        IF (D0 .LE. 0.D0  .OR.  S*T0 .GE. 0.D0) GO TO 7
C
C Monotonicity:  find a zero of F(SIG) = SIGN(S)*HP(R),
C   where HPP(R) = 0 and HP, HPP denote derivatives of H.
C   F has a unique zero, F(0) < 0, and F approaches abs(S)
C   as SIG increases.
C
C   Initialize parameters for the secant method.  The method
C     uses three points:  (SG0,F0), (SIG,F), and
C     (SNEG,FNEG), where SG0 and SNEG are defined implicitly
C     by DSIG = SIG - SG0 and DMAX = SIG - SNEG.
C
        SGN = SIGN(1.D0,S)
        SIG = SBIG
        FMAX = SGN*(SIG*S-S1-S2)/(SIG-2.D0)
        IF (FMAX .LE. 0.D0) GO TO 7
        STOL = RTOL*SIG
        F = FMAX
        F0 = SGN*D0/(3.D0*(D1-D2))
        FNEG = F0
        DSIG = SIG
        DMAX = SIG
        D1PD2 = D1 + D2
        NIT = 0
C
C   Top of loop:  compute the change in SIG by linear
C     interpolation.
C
    4   DSIG = -F*DSIG/(F-F0)
        IF (LUN .GE. 0) WRITE (LUN,120) DSIG
  120   FORMAT (5X,'MONOTONICITY -- DSIG = ',D15.8)
        IF ( ABS(DSIG) .GT. ABS(DMAX)  .OR.
     .       DSIG*DMAX .GT. 0. ) GO TO 6
C
C   Restrict the step-size such that abs(DSIG) .GE. STOL/2.
C     Note that DSIG and DMAX have opposite signs.
C
        IF (ABS(DSIG) .LT. STOL/2.D0)
     .    DSIG = -SIGN(STOL/2.D0,DMAX)
C
C   Update SIG, F0, and F.
C
        SIG = SIG + DSIG
        F0 = F
        IF (SIG .LE. .5D0) THEN
C
C   Use approximations to the hyperbolic functions designed
C     to avoid cancellation error with small SIG.
C
          CALL SNHCSH (SIG, SINHM,COSHM,COSHMM)
          C1 = SIG*COSHM*D2 - SINHM*D1PD2
          C2 = SIG*(SINHM+SIG)*D2 - COSHM*D1PD2
          A = C2 - C1
          E = SIG*SINHM - COSHMM - COSHMM
        ELSE
C
C   Scale SINHM and COSHM by 2*EXP(-SIG) in order to avoid
C     overflow with large SIG.
C
          EMS = EXP(-SIG)
          EMS2 = EMS + EMS
          TM = 1.D0 - EMS
          SSINH = TM*(1.D0+EMS)
          SSM = SSINH - SIG*EMS2
          SCM = TM*TM
          C1 = SIG*SCM*D2 - SSM*D1PD2
          C2 = SIG*SSINH*D2 - SCM*D1PD2
C
C   R is in (0,1) and well-defined iff HPP(X1)*HPP(X2) < 0.
C
          F = FMAX
          IF (C1*(SIG*SCM*D1 - SSM*D1PD2) .GE. 0.D0) GO TO 5
          A = EMS2*(SIG*TM*D2 + (TM-SIG)*D1PD2)
          IF (A*(C2+C1) .LT. 0.D0) GO TO 5
          E = SIG*SSINH - SCM - SCM
        ENDIF
C
        F = (SGN*(E*S2-C2) + SQRT(A*(C2+C1)))/E
C
C   Update number of iterations NIT.
C
    5   NIT = NIT + 1
        IF (LUN .GE. 0) WRITE (LUN,130) NIT, SIG, F
  130   FORMAT (1X,10X,I2,' -- SIG = ',D15.8,', F = ',
     .          D15.8)
C
C   Test for convergence.
C
        STOL = RTOL*SIG
        IF ( ABS(DMAX) .LE. STOL  .OR.  (F .GE. 0.D0  .AND.
     .      F .LE. FTOL)  .OR.  ABS(F) .LE. RTOL ) GO TO 7
        DMAX = DMAX + DSIG
        IF ( F0*F .GT. 0.D0  .AND.  ABS(F) .GE. ABS(F0) )
     .     GO TO 6
        IF (F0*F .LE. 0.D0) THEN
C
C   F and F0 have opposite signs.  Update (SNEG,FNEG) to
C     (SG0,F0) so that F and FNEG always have opposite
C     signs.  If SIG is closer to SNEG than SG0 and abs(F) <
C     abs(FNEG), then swap (SNEG,FNEG) with (SG0,F0).
C
          T1 = DMAX
          T2 = FNEG
          DMAX = DSIG
          FNEG = F0
          IF ( ABS(DSIG) .GT. ABS(T1)  .AND.
     .         ABS(F) .LT. ABS(T2)          ) THEN
C
            DSIG = T1
            F0 = T2
          ENDIF
        ENDIF
        GO TO 4
C
C   Bottom of loop:  F0*F > 0 and the new estimate would
C     be outside of the bracketing interval of length
C     abs(DMAX).  Reset (SG0,F0) to (SNEG,FNEG).
C
    6   DSIG = DMAX
        F0 = FNEG
        GO TO 4
C
C  Update SIGMA(I), ICNT, and DSM if necessary.
C
    7   SIG = MIN(SIG,SBIG)
        IF (SIG .GT. SIGIN) THEN
          SIGMA(I) = SIG
          ICNT = ICNT + 1
          DSIG = SIG-SIGIN
          IF (SIGIN .GT. 0.D0) DSIG = DSIG/SIGIN
          DSM = MAX(DSM,DSIG)
        ENDIF
    8   CONTINUE
C
C No errors encountered.
C
      DSMAX = DSM
      IER = ICNT
      RETURN
C
C N < 2.
C
    9 DSMAX = 0.D0
      IER = -1
      RETURN
C
C X(I+1) .LE. X(I).
C
   10 DSMAX = DSM
      IER = -IP1
      RETURN
      END
      SUBROUTINE SMCRV (N,X,Y,SIGMA,PERIOD,W,SM,SMTOL,
     .                  WK, YS,YP,IER)
      LOGICAL PERIOD
      INTEGER N, IER
      DOUBLE PRECISION X(N), Y(N), SIGMA(N), W(N), SM,
     .                 SMTOL, WK(N,10), YS(N), YP(N)
C
C***********************************************************
C
C                                                From TSPACK
C                                            Robert J. Renka
C                                  Dept. of Computer Science
C                                       Univ. of North Texas
C                                           renka@cs.unt.edu
C                                                   10/05/98
C
C   Given a sequence of abscissae X with associated data
C values Y and tension factors SIGMA, this routine deter-
C mines a set of function values YS and first derivatives YP
C associated with a Hermite interpolatory tension spline
C H(x) which smoothes the data.  H(x) has two continuous
C derivatives for all x and satisfies either natural or per-
C iodic end conditions.  The values and derivatives are
C chosen to minimize a quadratic functional Q1(YS,YP) sub-
C ject to the constraint Q2(YS) .LE. SM for Q2(YS) =
C (Y-YS)**T*W*(Y-YS), where **T denotes transpose and W is a
C diagonal matrix of positive weights.
C
C   Functions HVAL, HPVAL, HPPVAL, and TSINTL may be called
C to compute values, derivatives, and integrals of H.  The
C function values YS must be used as data values in those
C subprograms.
C
C   The smoothing procedure is an extension of the method
C for cubic spline smoothing due to C. Reinsch:  Numer.
C Math., 10 (1967) and 16 (1971).  Q1 is defined as the sum
C of integrals over the intervals (X(I),X(I+1)) of HPP**2 +
C (SIGMA(I)/DX)**2*(HP-S)**2, where DX = X(I+1)-X(I), HP and
C HPP denote first and second derivatives of H, and S =
C (YS(I+1)-YS(I))/DX.  Introducing a smoothing parameter P,
C and assuming the constraint is active, the problem is
C equivalent to minimizing Q(P,YS,YP) = Q1(YS,YP) +
C P*(Q2(YS)-SM).  The secant method is used to find a zero
C of G(P) = 1/SQRT(Q2) - 1/SQRT(SM), where YS and YP satisfy
C the order 2N symmetric positive-definite linear system
C obtained by setting the gradient of Q (treated as a func-
C tion of YS and YP) to zero.
C
C   Note that the interpolation problem corresponding to
C YS = Y, SM = 0, and P infinite is solved by Subroutine
C YPC2 or YPC2P.
C
C On input:
C
C       N = Number of data points.  N .GE. 2 if PERIOD =
C           FALSE, and N .GE. 3 if PERIOD = TRUE.
C
C       X = Array of length N containing a strictly in-
C           creasing sequence of abscissae:  X(I) < X(I+1)
C           for I = 1,...,N-1.
C
C       Y = Array of length N containing data values assoc-
C           iated with the abscissae.  If PERIOD = TRUE, it
C           is assumed that Y(N) = Y(1).
C
C       SIGMA = Array of length N-1 containing tension
C               factors.  SIGMA(I) is associated with inter-
C               val (X(I),X(I+1)) for I = 1,...,N-1.  If
C               SIGMA(I) = 0, H is cubic, and as SIGMA in-
C               creases, H approaches linear in the inter-
C               val.
C
C       PERIOD = Periodic end condition flag:
C                PERIOD = .F. if H is to satisfy natural end
C                             conditions:  zero second der-
C                             ivatives at X(1) and X(N).
C                PERIOD = .T. if H is to satisfy periodic
C                             end conditions:  the values
C                             and first two derivatives at
C                             X(1) agree with those at X(N),
C                             and a period thus has length
C                             X(N)-X(1).
C
C       W = Array of length N containing positive weights
C           associated with the data values.  The recommend-
C           ed value of W(I) is 1/DY**2, where DY is the
C           standard deviation associated with Y(I).  If
C           nothing is known about the errors in Y, a con-
C           stant (estimated value) should be used for DY.
C           If PERIOD = TRUE, it is assumed that W(N) =
C           W(1).
C
C       SM = Positive parameter specifying an upper bound on
C            Q2(YS).  H(x) is linear (and Q2 is minimized)
C            if SM is sufficiently large that the constraint
C            is not active.  It is recommended that SM sat-
C            isfy N-SQRT(2N) .LE. SM .LE. N+SQRT(2N) and
C            SM = N is reasonable if W(I) = 1/DY**2.
C
C       SMTOL = Parameter in the range (0,1) specifying the
C               relative error allowed in satisfying the
C               constraint:  the constraint is assumed to
C               be satisfied if SM*(1-SMTOL) .LE. Q2 .LE.
C               SM*(1+SMTOL).  A reasonable value for SMTOL
C               is SQRT(2/N).
C
C The above parameters are not altered by this routine.
C
C       WK = Work space of length at least 6N if PERIOD =
C            FALSE, and 10N if PERIOD = TRUE.
C
C On output:
C
C       YS = Array of length N containing values of H at the
C            abscissae unless IER < 0.  YS(N) = YS(1) if
C            PERIOD = TRUE.
C
C       YP = Array of length N containing first derivative
C            values of H at the abscissae unless IER < 0.
C            YP(N) = YP(1) if PERIOD = TRUE.
C
C       IER = Error indicator:
C             IER = 0 if no errors were encountered and the
C                     constraint is active:  Q2(YS) is ap-
C                     proximately equal to SM.
C             IER = 1 if no errors were encountered but the
C                     constraint is not active:  YS and YP
C                     are the values and derivatives of the
C                     linear function (constant function if
C                     PERIOD = TRUE) which minimizes Q2, and
C                     Q1 = 0.
C             IER = -1 if N, W, SM, or SMTOL is outside its
C                      valid range.  YS and YP are unaltered
C                      in this case.
C             IER = -I if X(I) .LE. X(I-1) for some I in the
C                      range 2,...,N.  YS and YP are unal-
C                      tered in this case.
C
C Modules required by SMCRV:  B2TRI or B2TRIP, SNHCSH,
C                               YPCOEF
C
C Intrinsic functions called by SMCRV:  ABS, SQRT
C
C***********************************************************
C
      INTEGER I, IERR, ITER, LUN, NM1, NN
      LOGICAL PER
      DOUBLE PRECISION C11, C12, C22, D, DMAX, DP, DX, G,
     .                 G0, GNEG, H0, HP, P, P0, Q2, Q2MAX,
     .                 Q2MIN, R1, R2, S, SD, SIG, WI, WIXI,
     .                 XI, YI
C
      DATA    LUN/-1/
      NN = N
      PER = PERIOD
C
C Test for errors, and compute the components of the system
C   (normal equations) for the weighted least squares linear
C   fit.
C
      IER = -1
      IF (NN .LT. 2  .OR.  (NN .LT. 3  .AND.  PER)  .OR.
     .    SM .LE. 0.D0  .OR.  SMTOL .LE. 0.D0  .OR.
     .    SMTOL .GE. 1.D0) RETURN
      C11 = 0.D0
      C12 = 0.D0
      C22 = 0.D0
      R1 = 0.D0
      R2 = 0.D0
      XI = X(1) - 1.D0
      DO 1 I = 1,NN
        WI = W(I)
        IF (WI .LE. 0.D0) RETURN
        IF (X(I) .LE. XI) THEN
          IER = -I
          RETURN
        ENDIF
        XI = X(I)
        YI = Y(I)
        C22 = C22 + WI
        R2 = R2 + WI*YI
        IF (.NOT. PER) THEN
          WIXI = WI*XI
          C11 = C11 + WIXI*XI
          C12 = C12 + WIXI
          R1 = R1 + WIXI*YI
        ENDIF
    1   CONTINUE
C
C Solve the system for (HP,H0), where HP is the derivative
C   (constant) and H0 = H(0).
C
      IF (PER) THEN
        H0 = R2/C22
        HP = 0.D0
      ELSE
        H0 = (C11*R2-C12*R1)/(C11*C22-C12*C12)
        HP = (R1 - C12*H0)/C11
      ENDIF
C
C Store function values and derivatives, and accumulate
C   Q2 = (Y-YS)**T*W*(Y-YS).
C
      Q2 = 0.D0
      DO 2 I = 1,NN
        YS(I) = HP*X(I) + H0
        YP(I) = HP
        Q2 = Q2 + W(I)*(Y(I)-YS(I))**2
    2   CONTINUE
C
C Compute bounds on Q2 defined by SMTOL, and test for the
C   constraint satisfied by the linear fit.
C
      Q2MIN = SM*(1.D0 - SMTOL)
      Q2MAX = SM*(1.D0 + SMTOL)
      IF (Q2 .LE. Q2MAX) THEN
C
C   The constraint is satisfied by a linear function.
C
        IER = 1
        IF (LUN .GE. 0) WRITE (LUN,100)
  100   FORMAT (///1X,'SMCRV -- THE CONSTRAINT IS NOT ',
     .          'ACTIVE AND THE FIT IS LINEAR.'/)
        RETURN
      ENDIF
C
C Compute the matrix components for the linear system.
C
      IER = 0
      NM1 = NN - 1
      DO 3 I = 1,NM1
        DX = X(I+1) - X(I)
        SIG = ABS(SIGMA(I))
        CALL YPCOEF (SIG,DX, D,SD)
        WK(I,1) = D
        WK(I,2) = SD
    3   CONTINUE
C
C Compute G0 = G(0), and print a heading.
C
      S = 1.D0/SQRT(SM)
      G0 = 1.D0/SQRT(Q2) - S
      IF (LUN .GE. 0) WRITE (LUN,110) SM, SMTOL, G0
  110 FORMAT (///1X,3X,'SMCRV -- SM = ',D10.4,', SMTOL = ',
     .        D14.8,', G(0) = ',D15.8///)
C
C G(P) is strictly increasing and concave, and G(0) < 0.
C
C Initialize parameters for the secant method.  The method
C   uses three points:  (P0,G0), (P,G), and (PNEG,GNEG),
C   where P0 and PNEG are defined implicitly by DP = P - P0
C   and DMAX = P - PNEG.
C
      P = 10.D0*SM
      DP = P
      DMAX = 0.D0
      ITER = 0
C
C Top of loop:  compute G and print a message.  For each
C               secant iteration, the following values are
C               printed:  P, G(P), and DP, where DP is the
C               change in P computed by linear interpolation
C               between the current point (P,G) and a previ-
C               ous point.
C
C
    4 IF (.NOT. PER) THEN
        CALL B2TRI (NN,X,Y,W,P,WK,WK(1,2),WK(1,3),WK(1,4),
     .              WK(1,5),WK(1,6), YS,YP,IERR)
      ELSE
        CALL B2TRIP (NN,X,Y,W,P,WK,WK(1,2),WK(1,3),WK(1,4),
     .               WK(1,5),WK(1,6),WK(1,7),WK(1,8),
     .               WK(1,9),WK(1,10), YS,YP,IERR)
      ENDIF
      Q2 = 0.D0
      DO 5 I = 1,NN
        Q2 = Q2 + W(I)*(Y(I)-YS(I))**2
    5   CONTINUE
      G = 1.D0/SQRT(Q2) - S
      ITER = ITER + 1
      IF (LUN .GE. 0) THEN
        P0 = P - DP
        IF (LUN .GE. 0) WRITE (LUN,120) ITER, P, G, P0, G0
  120   FORMAT (/1X,I2,' -- P = ',D15.8,',  G = ',D15.8/
     .          6X,'P0 = ',D15.8,', G0 = ',D15.8)
      ENDIF
C
C   Test for convergence.
C
      IF ( G .EQ. G0  .OR.  (Q2MIN .LE. Q2  .AND.
     .                       Q2 .LE. Q2MAX) )      RETURN
      IF (DMAX .NE. 0.D0  .OR.  G .GT. 0.D0) GO TO 6
C
C   Increase P until G(P) > 0.
C
      P = 10.D0*P
      DP = P
      GO TO 4
C
C   A bracketing interval [P0,P] has been found.
C
    6 IF (G0*G .LE. 0.D0) THEN
C
C   G0*G < 0.  Update (PNEG,GNEG) to (P0,G0) so that G
C     and GNEG always have opposite signs.
C
        DMAX = DP
        GNEG = G0
      ENDIF
C
C   Compute the change in P by linear interpolation between
C     (P0,G0) and (P,G).
C
    7 DP = -G*DP/(G-G0)
      IF (LUN .GE. 0) WRITE (LUN,130) DP
  130 FORMAT (1X ,5X,'DP = ',D15.8)
      IF (ABS(DP) .GT. ABS(DMAX)) THEN
C
C   G0*G > 0, and the new estimate would be outside of the
C     bracketing interval of length abs(DMAX).  Reset
C     (P0,G0) to (PNEG,GNEG).
C
        DP = DMAX
        G0 = GNEG
        GO TO 7
      ENDIF
C
C   Bottom of loop:  update P, DMAX, and G0.
C
      P = P + DP
      DMAX = DMAX + DP
      G0 = G
      GO TO 4
      END
      SUBROUTINE SNHCSH (X, SINHM,COSHM,COSHMM)
      DOUBLE PRECISION X, SINHM, COSHM, COSHMM
C
C***********************************************************
C
C                                                From TSPACK
C                                            Robert J. Renka
C                                  Dept. of Computer Science
C                                       Univ. of North Texas
C                                           renka@cs.unt.edu
C                                                   11/20/96
C
C   This subroutine computes approximations to the modified
C hyperbolic functions defined below with relative error
C bounded by 3.4E-20 for a floating point number system with
C sufficient precision.
C
C   Note that the 21-digit constants in the data statements
C below may not be acceptable to all compilers.
C
C On input:
C
C       X = Point at which the functions are to be
C           evaluated.
C
C X is not altered by this routine.
C
C On output:
C
C       SINHM = sinh(X) - X.
C
C       COSHM = cosh(X) - 1.
C
C       COSHMM = cosh(X) - 1 - X*X/2.
C
C Modules required by SNHCSH:  None
C
C Intrinsic functions called by SNHCSH:  ABS, EXP
C
C***********************************************************
C
      DOUBLE PRECISION AX, EXPX, F, P, P1, P2, P3, P4, Q,
     .                 Q1, Q2, Q3, Q4, XC, XS, XSD2, XSD4
C
      DATA P1/-3.51754964808151394800D5/,
     .     P2/-1.15614435765005216044D4/,
     .     P3/-1.63725857525983828727D2/,
     .     P4/-7.89474443963537015605D-1/
      DATA Q1/-2.11052978884890840399D6/,
     .     Q2/3.61578279834431989373D4/,
     .     Q3/-2.77711081420602794433D2/,
     .     Q4/1.D0/
      AX = ABS(X)
      XS = AX*AX
      IF (AX .LE. .5D0) THEN
C
C Approximations for small X:
C
        XC = X*XS
        P = ((P4*XS+P3)*XS+P2)*XS+P1
        Q = ((Q4*XS+Q3)*XS+Q2)*XS+Q1
        SINHM = XC*(P/Q)
        XSD4 = .25D0*XS
        XSD2 = XSD4 + XSD4
        P = ((P4*XSD4+P3)*XSD4+P2)*XSD4+P1
        Q = ((Q4*XSD4+Q3)*XSD4+Q2)*XSD4+Q1
        F = XSD4*(P/Q)
        COSHMM = XSD2*F*(F+2.D0)
        COSHM = COSHMM + XSD2
      ELSE
C
C Approximations for large X:
C
        EXPX = EXP(AX)
        SINHM = -(((1.D0/EXPX+AX)+AX)-EXPX)/2.D0
        IF (X .LT. 0.D0) SINHM = -SINHM
        COSHM = ((1.D0/EXPX-2.D0)+EXPX)/2.D0
        COSHMM = COSHM - XS/2.D0
      ENDIF
      RETURN
      END
      DOUBLE PRECISION FUNCTION STORE (X)
      DOUBLE PRECISION X
C
C***********************************************************
C
C                                                From TSPACK
C                                            Robert J. Renka
C                                  Dept. of Computer Science
C                                       Univ. of North Texas
C                                           renka@cs.unt.edu
C                                                   06/10/92
C
C   This function forces its argument X to be stored in a
C memory location, thus providing a means of determining
C floating point number characteristics (such as the machine
C precision) when it is necessary to avoid computation in
C high precision registers.
C
C On input:
C
C       X = Value to be stored.
C
C X is not altered by this function.
C
C On output:
C
C       STORE = Value of X after it has been stored and
C               possibly truncated or rounded to the single
C               precision word length.
C
C Modules required by STORE:  None
C
C***********************************************************
C
      DOUBLE PRECISION Y
      COMMON/STCOM/Y
      Y = X
      STORE = Y
      RETURN
      END
      DOUBLE PRECISION FUNCTION TSINTL (A,B,N,X,Y,YP,
     .                                  SIGMA, IER)
      INTEGER N, IER
      DOUBLE PRECISION A, B, X(N), Y(N), YP(N), SIGMA(N)
C
C***********************************************************
C
C                                                From TSPACK
C                                            Robert J. Renka
C                                  Dept. of Computer Science
C                                       Univ. of North Texas
C                                           renka@cs.unt.edu
C                                                   11/17/96
C
C   This function computes the integral from A to B of a
C Hermite interpolatory tension spline H.
C
C On input:
C
C       A,B = Lower and upper limits of integration, re-
C             spectively.  Note that -TSINTL(B,A,...) =
C             TSINTL(A,B,...).
C
C       N = Number of data points.  N .GE. 2.
C
C       X = Array of length N containing the abscissae.
C           These must be in strictly increasing order:
C           X(I) < X(I+1) for I = 1,...,N-1.
C
C       Y = Array of length N containing data values.
C           H(X(I)) = Y(I) for I = 1,...,N.
C
C       YP = Array of length N containing first deriva-
C            tives.  HP(X(I)) = YP(I) for I = 1,...,N, where
C            HP denotes the derivative of H.
C
C       SIGMA = Array of length N-1 containing tension fac-
C               tors whose absolute values determine the
C               balance between cubic and linear in each
C               interval.  SIGMA(I) is associated with int-
C               erval (I,I+1) for I = 1,...,N-1.
C
C Input parameters are not altered by this function.
C
C On output:
C
C       IER = Error indicator:
C             IER = 0  if no errors were encountered and
C                      X(1) .LE. T .LE. X(N) for T = A and
C                      T = B, or A = B.
C             IER = 1  if no errors were encountered but
C                      extrapolation was necessary:  A or B
C                      not in the interval (X(1),X(N)).
C             IER = -1 IF N < 2.
C             IER = -2 if the abscissae are not in strictly
C                      increasing order.  Only those in or
C                      adjacent to the interval of integra-
C                      tion are tested.
C
C       TSINTL = Integral of H from A to B, or zero if
C                IER < 0.
C
C Modules required by TSINTL:  INTRVL, SNHCSH
C
C Intrinsic functions called by TSINTL:  ABS, EXP, MAX, MIN
C
C***********************************************************
C
      INTEGER I, IL, ILP1, IMAX, IMIN, IP1, IU, IUP1
      DOUBLE PRECISION B1, B2, CM, CM1, CM2, CMM, CMM1,
     .                 CMM2, D1, D2, DX, E, E1, E2, EMS, S,
     .                 S1, S2, SB1, SB2, SBIG, SIG, SM, SM1,
     .                 SM2, SUM, T, TM, TP, U, XL, XU, Y1,
     .                 Y2
      INTEGER INTRVL
C
      DATA SBIG/85.D0/
      IF (N .LT. 2) GO TO 7
C
C Accumulate the integral from XL to XU in SUM.
C
      XL = MIN(A,B)
      XU = MAX(A,B)
      SUM = 0.D0
      IER = 0
      IF (XL .EQ. XU) GO TO 6
C
C Find left-end indexes of intervals containing XL and XU.
C   If XL < X(1) or XU > X(N), extrapolation is performed
C   using the leftmost or rightmost interval.
C
      IL = INTRVL (XL,N,X)
      IU = INTRVL (XU,N,X)
      IF (XL .LT. X(1)  .OR.  XU .GT. X(N)) IER = 1
      ILP1 = IL + 1
      IMIN = IL
      IF (XL .EQ. X(IL)) GO TO 2
C
C Compute the integral from XL to X(IL+1).
C
      DX = X(ILP1) - X(IL)
      IF (DX .LE. 0.D0) GO TO 8
      U = X(ILP1) - XL
      IF (U .EQ. 0.D0) GO TO 1
      B1 = U/DX
      Y2 = Y(ILP1)
      S = (Y2-Y(IL))/DX
      S2 = YP(ILP1)
      D1 = S - YP(IL)
      D2 = S2 - S
      SIG = ABS(SIGMA(IL))
      IF (SIG .LT. 1.D-9) THEN
C
C   SIG = 0.
C
        SUM = SUM + U*(Y2 - U*(6.D0*S2 - B1*(4.D0*D2 +
     .                 (3.D0*B1-4.D0)*(D1-D2)))/12.D0)
      ELSEIF (SIG .LE. .5D0) THEN
C
C   0 .LT. SIG .LE. .5.
C
        SB1 = SIG*B1
        CALL SNHCSH (SIG, SM,CM,CMM)
        CALL SNHCSH (SB1, SM1,CM1,CMM1)
        E = SIG*SM - CMM - CMM
        SUM = SUM + U*(Y2 - S2*U/2.D0) + ((CM*CMM1-SM*SM1)*
     .             (D1+D2) + SIG*(CM*SM1-(SM+SIG)*CMM1)*D2)/
     .             ((SIG/DX)**2*E)
      ELSE
C
C   SIG > .5.
C
        SB1 = SIG*B1
        SB2 = SIG - SB1
        IF (-SB1 .GT. SBIG  .OR.  -SB2 .GT. SBIG) THEN
          SUM = SUM + U*(Y2 - S*U/2.D0)
        ELSE
          E1 = EXP(-SB1)
          E2 = EXP(-SB2)
          EMS = E1*E2
          TM = 1.D0 - EMS
          TP = 1.D0 + EMS
          T = SB1*SB1/2.D0 + 1.D0
          E = TM*(SIG*TP - TM - TM)
          SUM = SUM +U*(Y2 - S2*U/2.D0)+(SIG*TM*(TP*T-E1-E2-
     .               TM*SB1)*D2 - (TM*(TM*T-E1+E2-TP*SB1) +
     .               SIG*(E1*EMS-E2+2.D0*SB1*EMS))*(D1+D2))/
     .               ((SIG/DX)**2*E)
        ENDIF
      ENDIF
C
C Test for XL and XU in the same interval.
C
    1 IMIN = ILP1
      IF (IL .EQ. IU) GO TO 5
C
C Add in the integral from X(IMIN) to X(J) for J =
C   Max(IL+1,IU).
C
    2 IMAX = MAX(IL,IU-1)
      DO 3 I = IMIN,IMAX
        IP1 = I + 1
        DX = X(IP1) - X(I)
        IF (DX .LE. 0.D0) GO TO 8
        SIG = ABS(SIGMA(I))
        IF (SIG .LT. 1.D-9) THEN
C
C   SIG = 0.
C
          SUM = SUM + DX*((Y(I)+Y(IP1))/2.D0 -
     .                    DX*(YP(IP1)-YP(I))/12.D0)
        ELSEIF (SIG .LE. .5D0) THEN
C
C   0 .LT. SIG .LE. .5.
C
          CALL SNHCSH (SIG, SM,CM,CMM)
          E = SIG*SM - CMM - CMM
          SUM = SUM + DX*(Y(I)+Y(IP1) - DX*E*(YP(IP1)-YP(I))
     .                /(SIG*SIG*CM))/2.D0
        ELSE
C
C   SIG > .5.
C
          EMS = EXP(-SIG)
          SUM = SUM + DX*(Y(I)+Y(IP1) - DX*(SIG*(1.D0+EMS)/
     .                (1.D0-EMS)-2.D0)*(YP(IP1)-YP(I))/
     .                (SIG*SIG))/2.D0
        ENDIF
    3   CONTINUE
C
C Add in the integral from X(IU) to XU if IU > IL.
C
      IF (IL .EQ. IU) GO TO 4
      IUP1 = IU + 1
      DX = X(IUP1) - X(IU)
      IF (DX .LE. 0.D0) GO TO 8
      U = XU - X(IU)
      IF (U .EQ. 0.D0) GO TO 6
      B2 = U/DX
      Y1 = Y(IU)
      S = (Y(IUP1)-Y1)/DX
      S1 = YP(IU)
      D1 = S - S1
      D2 = YP(IUP1) - S
      SIG = ABS(SIGMA(IU))
      IF (SIG .LT. 1.D-9) THEN
C
C   SIG = 0.
C
        SUM = SUM + U*(Y1 + U*(6.D0*S1 + B2*(4.D0*D1 +
     .                (4.D0-3.D0*B2)*(D1-D2)))/12.D0)
      ELSEIF (SIG .LE. .5D0) THEN
C
C   0 .LT. SIG .LE. .5.
C
        SB2 = SIG*B2
        CALL SNHCSH (SIG, SM,CM,CMM)
        CALL SNHCSH (SB2, SM2,CM2,CMM2)
        E = SIG*SM - CMM - CMM
        SUM = SUM + U*(Y1 + S1*U/2.D0) + ((CM*CMM2-SM*SM2)*
     .              (D1+D2) + SIG*(CM*SM2-(SM+SIG)*CMM2)*D1)
     .              /((SIG/DX)**2*E)
      ELSE
C
C   SIG > .5.
C
        SB2 = SIG*B2
        SB1 = SIG - SB2
        IF (-SB1 .GT. SBIG  .OR.  -SB2 .GT. SBIG) THEN
          SUM = SUM + U*(Y1 + S*U/2.D0)
        ELSE
          E1 = EXP(-SB1)
          E2 = EXP(-SB2)
          EMS = E1*E2
          TM = 1.D0 - EMS
          TP = 1.D0 + EMS
          T = SB2*SB2/2.D0 + 1.D0
          E = TM*(SIG*TP - TM - TM)
          SUM = SUM +U*(Y1 + S1*U/2.D0)+(SIG*TM*(TP*T-E1-E2-
     .               TM*SB2)*D1 - (TM*(TM*T-E2+E1-TP*SB2) +
     .               SIG*(E2*EMS-E1+2.D0*SB2*EMS))*(D1+D2))/
     .               ((SIG/DX)**2*E)
        ENDIF
      ENDIF
      GO TO 6
C
C IL = IU and SUM contains the integral from XL to X(IL+1).
C   Subtract off the integral from XU to X(IL+1).  DX and
C   SIG were computed above.
C
    4 Y2 = Y(ILP1)
      S = (Y2-Y(IL))/DX
      S2 = YP(ILP1)
      D1 = S - YP(IL)
      D2 = S2 - S
C
    5 U = X(ILP1) - XU
      IF (U .EQ. 0.D0) GO TO 6
      B1 = U/DX
      IF (SIG .LT. 1.D-9) THEN
C
C   SIG = 0.
C
        SUM = SUM - U*(Y2 - U*(6.D0*S2 - B1*(4.D0*D2 +
     .              (3.D0*B1-4.D0)*(D1-D2)))/12.D0)
      ELSEIF (SIG .LE. .5D0) THEN
C
C   0 .LT. SIG .LE. .5.
C
        SB1 = SIG*B1
        CALL SNHCSH (SIG, SM,CM,CMM)
        CALL SNHCSH (SB1, SM1,CM1,CMM1)
        E = SIG*SM - CMM - CMM
        SUM = SUM - U*(Y2 - S2*U/2.D0) - ((CM*CMM1-SM*SM1)*
     .              (D1+D2) + SIG*(CM*SM1-(SM+SIG)*CMM1)*D2)
     .              /((SIG/DX)**2*E)
      ELSE
C
C   SIG > .5.
C
        SB1 = SIG*B1
        SB2 = SIG - SB1
        IF (-SB1 .GT. SBIG  .OR.  -SB2 .GT. SBIG) THEN
          SUM = SUM - U*(Y2 - S*U/2.D0)
        ELSE
          E1 = EXP(-SB1)
          E2 = EXP(-SB2)
          EMS = E1*E2
          TM = 1.D0 - EMS
          TP = 1.D0 + EMS
          T = SB1*SB1/2.D0 + 1.D0
          E = TM*(SIG*TP - TM - TM)
          SUM = SUM -U*(Y2 - S2*U/2.D0)-(SIG*TM*(TP*T-E1-E2-
     .               TM*SB1)*D2 - (TM*(TM*T-E1+E2-TP*SB1) +
     .               SIG*(E1*EMS-E2+2.D0*SB1*EMS))*(D1+D2))/
     .               ((SIG/DX)**2*E)
        ENDIF
      ENDIF
C
C No errors were encountered.  Adjust the sign of SUM.
C
    6 IF (XL .EQ. B) SUM = -SUM
      TSINTL = SUM
      RETURN
C
C N < 2.
C
    7 IER = -1
      TSINTL = 0.D0
      RETURN
C
C Abscissae not strictly increasing.
C
    8 IER = -2
      TSINTL = 0.D0
      RETURN
      END
      SUBROUTINE TSPBI (N,X,Y,NCD,IENDC,PER,B,BMAX,LWK, WK,
     .                  YP, SIGMA,ICFLG,IER)
      INTEGER N, NCD, IENDC, LWK, ICFLG(N), IER
      LOGICAL PER
      DOUBLE PRECISION X(N), Y(N), B(5,N), BMAX, WK(LWK),
     .                 YP(N), SIGMA(N)
C
C***********************************************************
C
C                                                From TSPACK
C                                            Robert J. Renka
C                                  Dept. of Computer Science
C                                       Univ. of North Texas
C                                           renka@cs.unt.edu
C                                                   07/08/92
C
C   This subroutine computes a set of parameter values which
C define a Hermite interpolatory tension spline H(x).  The
C parameters consist of knot derivative values YP computed
C by Subroutine YPC1, YPC1P, YPC2, or YPC2P, and tension
C factors SIGMA chosen to satisfy user-specified constraints
C (by Subroutine SIGBI).  Refer to Subroutine TSPSI for an
C alternative method of computing tension factors.
C
C   Refer to Subroutine TSPSS for a means of computing
C parameters which define a smoothing curve rather than an
C interpolatory curve.
C
C   The tension spline may be evaluated by Subroutine TSVAL1
C or Functions HVAL (values), HPVAL (first derivatives),
C HPPVAL (second derivatives), and TSINTL (integrals).
C
C On input:
C
C       N = Number of data points.  N .GE. 2 and N .GE. 3 if
C           PER = TRUE.
C
C       X = Array of length N containing a strictly in-
C           creasing sequence of abscissae:  X(I) < X(I+1)
C           for I = 1,...,N-1.
C
C       Y = Array of length N containing data values asso-
C           ciated with the abscissae.  H(X(I)) = Y(I) for
C           I = 1,...,N.  If NCD = 1 and PER = TRUE, Y(N)
C           is set to Y(1).
C
C       NCD = Number of continuous derivatives at the knots.
C             NCD = 1 or NCD = 2.  If NCD = 1, the YP values
C             are computed by local monotonicity-constrained
C             quadratic fits.  Otherwise, a linear system is
C             solved for the derivative values which result
C             in second derivative continuity.  This re-
C             quires iterating on calls to YPC2 or YPC2P and
C             calls to SIGBI, and generally results in more
C             nonzero tension factors (hence more expensive
C             evaluation).
C
C       IENDC = End condition indicator for NCD = 2 and PER
C               = FALSE (or dummy parameter otherwise):
C               IENDC = 0 if YP(1) and YP(N) are to be com-
C                         puted by monotonicity-constrained
C                         parabolic fits to the first three
C                         and last three points, respective-
C                         ly.  This is identical to the
C                         values computed by YPC1.
C               IENDC = 1 if the first derivatives of H at
C                         X(1) and X(N) are user-specified
C                         in YP(1) and YP(N), respectively.
C               IENDC = 2 if the second derivatives of H at
C                         X(1) and X(N) are user-specified
C                         in YP(1) and YP(N), respectively.
C               IENDC = 3 if the end conditions are to be
C                         computed by Subroutine ENDSLP and
C                         vary with SIGMA(1) and SIGMA(N-1).
C
C       PER = Logical variable with value TRUE if and only
C             H(x) is to be a periodic function with period
C             X(N)-X(1).  It is assumed without a test that
C             Y(N) = Y(1) in this case.  On output, YP(N) =
C             YP(1).  If H(x) is one of the components of a
C             parametric curve, this option may be used to
C             obtained a closed curve.
C
C       B = Array dimensioned 5 by N-1 containing bounds or
C           flags which define the constraints.  For I = 1
C           to N-1, column I defines the constraints associ-
C           ated with interval (X(I),X(I+1)) as follows:
C
C             B(1,I) is an upper bound on H
C             B(2,I) is a lower bound on H
C             B(3,I) is an upper bound on HP
C             B(4,I) is a lower bound on HP
C             B(5,I) specifies the required sign of HPP
C
C           where HP and HPP denote the first and second
C           derivatives of H, respectively.  A null con-
C           straint is specified by abs(B(K,I)) .GE. BMAX
C           for K < 5, or B(5,I) = 0:   B(1,I) .GE. BMAX,
C           B(2,I) .LE. -BMAX, B(3,I) .GE. BMAX, B(4,I) .LE.
C           -BMAX, or B(5,I) = 0.  Any positive value of
C           B(5,I) specifies that H should be convex, a
C           negative values specifies that H should be con-
C           cave, and 0 specifies that no restriction be
C           placed on HPP.  Refer to Functions SIG0, SIG1,
C           and SIG2 for definitions of valid constraints.
C
C       BMAX = User-defined value of infinity which, when
C              used as an upper bound in B (or when when
C              its negative is used as a lower bound),
C              specifies that no constraint is to be en-
C              forced.
C
C       LWK = Length of work space WK:
C             LWK GE 2N-2 if NCD = 2 and PER = FALSE
C             LWK GE 3N-3 if NCD = 2 and PER = TRUE
C
C   The above parameters, except possibly Y(N), are not
C altered by this routine.
C
C       WK = Array of length at least LWK to be used as
C            temporary work space.
C
C       YP = Array of length .GE. N containing end condition
C            values in positions 1 and N if NCD = 2 and
C            IENDC = 1 or IENDC = 2.
C
C       SIGMA = Array of length .GE. N-1.
C
C       ICFLG = Array of length .GE. N-1.
C
C On output:
C
C       WK = Array containing convergence parameters in the
C            first two locations if IER > 0 (NCD = 2 and
C            no error other than invalid constraints was
C            encountered):
C            WK(1) = Maximum relative change in a component
C                    of YP on the last iteration.
C            WK(2) = Maximum relative change in a component
C                    of SIGMA on the last iteration.
C
C       YP = Array containing derivatives of H at the
C            abscissae.  YP is not altered if -3 < IER < 0,
C            and YP is only partially defined if IER = -4.
C
C       SIGMA = Array containing tension factors for which
C               H(x) satisfies the constraints defined by B.
C               SIGMA(I) is associated with interval (X(I),
C               X(I+1)) for I = 1,...,N-1.  If infinite ten-
C               sion is required in interval I, then
C               SIGMA(I) = 85 (and H is an approximation to
C               the linear interpolant on the interval),
C               and if no constraint is specified in the
C               interval, then SIGMA(I) = 0, and thus H is
C               cubic.  Invalid constraints are treated as
C               null constraints.  SIGMA is not altered if
C               -3 < IER < 0 (unless IENDC is invalid), and
C               SIGMA is the zero vector if IER = -4 or
C               IENDC (if used) is invalid.
C
C       ICFLG = Array of invalid constraint flags associated
C               with intervals.  For I = 1 to N-1, ICFLG(I)
C               is a 5-bit value b5b4b3b2b1, where bK = 1 if
C               and only if constraint K cannot be satis-
C               fied.  Thus, all constraints in interval I
C               are satisfied if and only if ICFLG(I) = 0
C               (and IER .GE. 0).  ICFLG is not altered if
C               IER < 0.
C
C       IER = Error indicator or iteration count:
C             IER = IC .GE. 0 if no errors were encountered
C                      (other than invalid constraints) and
C                      IC calls to SIGBI and IC+1 calls to
C                      YPC1, YPC1P, YPC2 or YPC2P were
C                      employed.  (IC = 0 if NCD = 1).
C             IER = -1 if N, NCD, or IENDC is outside its
C                      valid range.
C             IER = -2 if LWK is too small.
C             IER = -4 if the abscissae X are not strictly
C                      increasing.
C
C Modules required by TSPBI:  ENDSLP, SIG0, SIG1, SIG2,
C                               SIGBI, SNHCSH, STORE,
C                               YPCOEF, YPC1, YPC1P, YPC2,
C                               YPC2P
C
C Intrinsic functions called by TSPBI:  ABS, MAX
C
C***********************************************************
C
      INTEGER I, ICNT, IERR, ITER, MAXIT, NM1, NN
      LOGICAL LOOP2
      DOUBLE PRECISION DSMAX, DYP, DYPTOL, E, STOL, YP1, YPN
C
      DATA STOL/0.D0/,  MAXIT/49/,  DYPTOL/.01D0/
C
C Convergence parameters:
C
C   STOL = Absolute tolerance for SIGBI.
C   MAXIT = Maximum number of YPC2/SIGBI iterations for
C             each loop if NCD = 2.
C   DYPTOL = Bound on the maximum relative change in a
C              component of YP defining convergence of
C              the YPC2/SIGBI iteration when NCD = 2.
C
      NN = N
      NM1 = NN - 1
C
C Test for invalid input parameters N, NCD, or LWK.
C
      IF (NN .LT. 2  .OR.  (PER  .AND.  NN .LT. 3)  .OR.
     .    NCD .LT. 1  .OR.  NCD .GT. 2) GO TO 11
      IF ( NCD .EQ. 2  .AND.  (LWK .LT. 2*NM1  .OR.
     .     (PER  .AND.  LWK .LT. 3*NM1)) ) GO TO 12
C
C Initialize iteration count ITER, and initialize SIGMA to
C   zeros.
C
      ITER = 0
      DO 1 I = 1,NM1
        SIGMA(I) = 0.D0
    1   CONTINUE
      IF (NCD .EQ. 1) THEN
C
C NCD = 1.
C
        IF (.NOT. PER) THEN
          CALL YPC1 (NN,X,Y, YP,IERR)
        ELSE
          CALL YPC1P (NN,X,Y, YP,IERR)
        ENDIF
        IF (IERR .NE. 0) GO TO 14
C
C   Compute tension factors.
C
        CALL SIGBI (NN,X,Y,YP,STOL,B,BMAX, SIGMA, ICFLG,
     .              DSMAX,IERR)
        GO TO 10
      ENDIF
C
C NCD = 2.
C
      IF (.NOT. PER) THEN
C
C   Nonperiodic case:  call YPC2 and test for IENDC or X
C     invalid.
C
        YP1 = YP(1)
        YPN = YP(NN)
        CALL YPC2 (NN,X,Y,SIGMA,IENDC,IENDC,YP1,YPN,
     .             WK, YP,IERR)
        IF (IERR .EQ. 1) GO TO 11
        IF (IERR .GT. 1) GO TO 14
      ELSE
C
C   Periodic fit:  call YPC2P.
C
        CALL YPC2P (NN,X,Y,SIGMA,WK, YP,IERR)
        IF (IERR .GT. 1) GO TO 14
      ENDIF
      LOOP2 = .FALSE.
C
C   Iterate on calls to SIGBI and YPC2 (or YPC2P).  The
C     first N-1 WK locations are used to store the deriva-
C     tive estimates YP from the previous iteration.
C
C   LOOP2 is TRUE iff tension factors are not allowed to
C         decrease between iterations (loop 1 failed to
C         converge with MAXIT iterations).
C   DYP is the maximum relative change in a component of YP.
C   ICNT is the number of tension factors which were altered
C        by SIGBI.
C   DSMAX is the maximum relative change in a component of
C         SIGMA.
C
    2 DO 6 ITER = 1,MAXIT
        DYP = 0.D0
        DO 3 I = 2,NM1
          WK(I) = YP(I)
    3     CONTINUE
        CALL SIGBI (NN,X,Y,YP,STOL,B,BMAX, SIGMA, ICFLG,
     .              DSMAX,ICNT)
        IF (.NOT. PER) THEN
          CALL YPC2 (NN,X,Y,SIGMA,IENDC,IENDC,YP1,YPN,
     .               WK(NN), YP,IERR)
        ELSE
          CALL YPC2P (NN,X,Y,SIGMA,WK(NN), YP,IERR)
        ENDIF
        DO 4 I = 2,NM1
          E = ABS(YP(I)-WK(I))
          IF (WK(I) .NE. 0.D0) E = E/ABS(WK(I))
          DYP = MAX(DYP,E)
    4     CONTINUE
        IF (ICNT .EQ. 0  .OR.  DYP .LE. DYPTOL) GO TO 7
        IF (.NOT. LOOP2) THEN
C
C   Loop 1:  reinitialize SIGMA to zeros.
C
          DO 5 I = 1,NM1
            SIGMA(I) = 0.D0
    5       CONTINUE
        ENDIF
    6   CONTINUE
C
C   The loop failed to converge within MAXIT iterations.
C
      ITER = MAXIT
      IF (.NOT. LOOP2) THEN
        LOOP2 = .TRUE.
        GO TO 2
      ENDIF
C
C Store convergence parameters.
C
    7 WK(1) = DYP
      WK(2) = DSMAX
      IF (LOOP2) ITER = ITER + MAXIT
C
C No error encountered.
C
   10 IER = ITER
      RETURN
C
C Invalid input parameter N, NCD, or IENDC.
C
   11 IER = -1
      RETURN
C
C LWK too small.
C
   12 IER = -2
      RETURN
C
C Abscissae are not strictly increasing.
C
   14 IER = -4
      RETURN
      END
      SUBROUTINE TSPBP (N,X,Y,NCD,IENDC,PER,BL,BU,BMAX,
     .                  LWK, WK, T,XP,YP,SIGMA,IER)
      INTEGER N, NCD, IENDC, LWK, IER
      LOGICAL PER
      DOUBLE PRECISION X(N), Y(N), BL(N), BU(N), BMAX,
     .                 WK(LWK), T(N), XP(N), YP(N), SIGMA(N)
C
C***********************************************************
C
C                                                From TSPACK
C                                            Robert J. Renka
C                                  Dept. of Computer Science
C                                       Univ. of North Texas
C                                           renka@cs.unt.edu
C                                                   07/08/92
C
C   This subroutine computes a set of values which define a
C parametric planar curve C(t) = (H1(t),H2(t)) whose compo-
C nents are Hermite interpolatory tension splines.  The
C output values consist of parameters (knots) T computed by
C ARCL2D, knot derivative values XP and YP computed by Sub-
C routine YPC1, YPC1P, YPC2, or YPC2P, and tension factors
C SIGMA chosen (by Subroutine SIGBP) to satisfy user-
C specified bounds on the distance between C(t) and the
C polygonal curve associated with the data points (refer to
C BL and BU below).
C
C   Refer to Subroutine TSPSP for an alternative method of
C computing tension factors.
C
C   The tension splines may be evaluated by Subroutine
C TSVAL2 or Functions HVAL (values), HPVAL (first deriva-
C tives), HPPVAL (second derivatives), and TSINTL
C (integrals).
C
C On input:
C
C       N = Number of knots and data points.  N .GE. 2 and
C           N .GE. 3 if PER = TRUE.
C
C       X,Y = Arrays of length N containing the Cartesian
C             coordinates of an ordered sequence of data
C             points C(I), I = 1 to N, such that C(I) .NE.
C             C(I+1).  C(t) is constrained to pass through
C             these points.  In the case of a closed curve
C             (PER = TRUE), the first and last points should
C             coincide.  (X(N) and Y(N) are set to X(1) and
C             Y(1) if NCD = 1, but not altered if NCD = 2,
C             in this case.)
C
C       NCD = Number of continuous derivatives at the knots.
C             NCD = 1 or NCD = 2.  If NCD = 1, XP and YP are
C             computed by local monotonicity-constrained
C             quadratic fits.  Otherwise, a linear system is
C             solved for the derivative values which result
C             in second derivative continuity.  This re-
C             quires iterating on calls to YPC2 or YPC2P and
C             calls to SIGBP, and generally results in more
C             nonzero tension factors (hence more expensive
C             evaluation).
C
C       IENDC = End condition indicator for NCD = 2 and PER
C               = FALSE (or dummy parameter otherwise):
C               IENDC = 0 if XP(1), XP(N), YP(1), and YP(N)
C                         are to be computed by monotonicity-
C                         constrained parabolic fits (YPC1).
C               IENDC = 1 if the first derivatives of H1 at
C                         the left and right endpoints are
C                         user-specified in XP(1) and XP(N),
C                         respectively, and the first deriv-
C                         atives of H2 at the ends are
C                         specified in YP(1) and YP(N).
C               IENDC = 2 if the second derivatives of H1
C                         and H2 at the endpoints are user-
C                         specified in XP(1), XP(N), YP(1),
C                         and YP(N).
C               IENDC = 3 if the end conditions are to be
C                         computed by Subroutine ENDSLP and
C                         vary with SIGMA(1) and SIGMA(N-1).
C
C       PER = Logical variable with value TRUE if and only
C             a closed curve is to be constructed -- H1(t)
C             and H2(t) are to be periodic functions with
C             period T(N)-T(1), where T(1) and T(N) are the
C             parameter values associated with the first and
C             last data points.  It is assumed that X(N) =
C             X(1) and Y(N) = Y(1) in this case, and, on
C             output, XP(N) = XP(1) and YP(N) = YP(1).
C
C       BL,BU = Arrays of length N-1 containing (for each
C               knot subinterval) lower and upper bounds,
C               respectively, on the signed perpendicular
C               distance d(t) = (C2-C1)/DC X (C(t)-C1),
C               where C1 and C2 are the ordered data points
C               associated with the interval, and DC is the
C               interval length (and length of the line seg-
C               ment C1-C2).  Note that d(t) > 0 iff C(t)
C               lies strictly to the left of the line seg-
C               ment as viewed from C1 toward C2.  For I =
C               1 to N-1, SIGMA(I) is chosen to be as small
C               as possible within the constraint that
C               BL(I) .LE. d(t) .LE. BU(I) for all t in the
C               interval.  BL(I) < 0 and BU(I) > 0 for I = 1
C               to N-1.  A null constraint is specified by
C               BL(I) .LE. -BMAX or BU(I) .GE. BMAX.
C
C       BMAX = User-defined value of infinity which, when
C              used as an upper bound in BU (or when its
C              negative is used as a lower bound in BL),
C              specifies that no constraint is to be en-
C              forced.
C
C       LWK = Length of work space WK:
C             LWK GE 3N-3 if NCD = 2 and PER = FALSE
C             LWK GE 4N-4 if NCD = 2 and PER = TRUE
C
C   The above parameters, except possibly X(N) and Y(N), are
C not altered by this routine.
C
C       WK = Array of length .GE. LWK to be used as tempor-
C            ary work space.
C
C       T = Array of length .GE. N.
C
C       XP = Array of length .GE. N containing end condition
C            values in positions 1 and N if NCD = 2 and
C            IENDC = 1 or IENDC = 2.
C
C       YP = Array of length .GE. N containing end condition
C            values in positions 1 and N if NCD = 2 and
C            IENDC = 1 or IENDC = 2.
C
C       SIGMA = Array of length .GE. N-1.
C
C On output:
C
C       WK = Array containing convergence parameters in the
C            first two locations if IER > 0 (NCD = 2 and
C            no error was encountered):
C            WK(1) = Maximum relative change in a component
C                    of XP or YP on the last iteration.
C            WK(2) = Maximum relative change in a component
C                    of SIGMA on the last iteration.
C
C       T = Array containing parameter values computed by
C           Subroutine ARCL2D unless IER = -1 or IER = -2.
C           T is only partially defined if IER = -4.
C
C       XP = Array containing derivatives of H1 at the
C            knots.  XP is not altered if -5 < IER < 0,
C            and XP is only partially defined if IER = -6.
C
C       YP = Array containing derivatives of H2 at the
C            knots.  YP is not altered if -5 < IER < 0,
C            and YP is only partially defined if IER = -6.
C
C       SIGMA = Array containing tension factors for which
C               C(t) satisfies the constraints defined by
C               BL and BU.  SIGMA(I) is associated with
C               interval (T(I),T(I+1)) for I = 1,...,N-1.
C               SIGMA(I) is limited to 85 (in which case
C               C(t) is indistinguishable from the line
C               segment associated with the interval), and
C               if no constraint is specified in the
C               interval, then SIGMA(I) = 0, and thus H1 and
C               H2 are cubic functions of t.  SIGMA is not
C               altered if -5 < IER < 0 (unless IENDC in
C               invalid), and SIGMA is the zero vector if
C               IER = -6 or IENDC (if used) is invalid.
C
C       IER = Error indicator or iteration count:
C             IER = IC .GE. 0 if no errors were encountered
C                      and IC calls to SIGBP and IC+1 calls
C                      to YPC1, YPC1P, YPC2 or YPC2P were
C                      employed.  (IC = 0 if NCD = 1).
C             IER = -1 if N, NCD, or IENDC is outside its
C                      valid range.
C             IER = -2 if LWK is too small.
C             IER = -4 if a pair of adjacent data points
C                      coincide:  X(I) = X(I+1) and Y(I) =
C                      Y(I+1) for some I in the range 1 to
C                      N-1.
C             IER = -5 if BL(I) .GE. 0 or BU(I) .LE. 0 for
C                      some I in the range 1 to N-1.
C                      SIGMA(J) = 0 for J .GE. I in this
C                      case.
C             IER = -6 if invalid knots T were returned by
C                      ARCL2D.  This should not occur.
C
C Modules required by TSPBP:  ARCL2D, ENDSLP, SIGBP, SNHCSH,
C                               STORE, YPCOEF, YPC1, YPC1P,
C                               YPC2, YPC2P
C
C Intrinsic functions called by TSPBP:  ABS, MAX
C
C***********************************************************
C
      INTEGER I, ICNT, IERR, ITER, MAXIT, N2M1, NM1, NN
      LOGICAL LOOP2
      DOUBLE PRECISION DSMAX, DYP, DYPTOL, EX, EY, STOL,
     .                 XP1, XPN, YP1, YPN
C
      DATA STOL/0.D0/,  MAXIT/49/,  DYPTOL/.01D0/
C
C Convergence parameters:
C
C   STOL = Absolute tolerance for SIGBP.
C   MAXIT = Maximum number of YPC2/SIGBP iterations for each
C             loop in NCD = 2.
C   DYPTOL = Bound on the maximum relative change in a
C              component of XP or YP defining convergence
C              of the YPC2/SIGBP iteration when NCD = 2.
C
      NN = N
      NM1 = NN - 1
C
C Test for invalid input parameters N, NCD, or LWK.
C
      N2M1 = 2*NN - 1
      IF (NN .LT. 2  .OR.  (PER  .AND.  NN .LT. 3)  .OR.
     .    NCD .LT. 1  .OR.  NCD .GT. 2) GO TO 11
      IF ( NCD .EQ. 2  .AND.  (LWK .LT. 3*NM1  .OR.
     .     (PER  .AND.  LWK .LT. 4*NM1)) ) GO TO 12
C
C Compute the sequence of parameter values T.
C
      CALL ARCL2D (NN,X,Y, T,IERR)
      IF (IERR .GT. 0) GO TO 14
C
C Initialize iteration count ITER, and initialize SIGMA to
C   zeros.
C
      ITER = 0
      DO 1 I = 1,NM1
        SIGMA(I) = 0.D0
    1   CONTINUE
      IF (NCD .EQ. 1) THEN
C
C NCD = 1.
C
        IF (.NOT. PER) THEN
          CALL YPC1 (NN,T,X, XP,IERR)
          CALL YPC1 (NN,T,Y, YP,IERR)
        ELSE
          CALL YPC1P (NN,T,X, XP,IERR)
          CALL YPC1P (NN,T,Y, YP,IERR)
        ENDIF
        IF (IERR .NE. 0) GO TO 16
C
C   Compute tension factors.
C
        CALL SIGBP (NN,X,Y,XP,YP,STOL,BL,BU,
     .              BMAX, SIGMA, DSMAX,IERR)
        IF (IERR .LT. 0) GO TO 15
        GO TO 10
      ENDIF
C
C NCD = 2.
C
      IF (.NOT. PER) THEN
C
C   Nonperiodic case:  call YPC2 and test for IENDC invalid.
C
        XP1 = XP(1)
        XPN = XP(NN)
        CALL YPC2 (NN,T,X,SIGMA,IENDC,IENDC,XP1,XPN,
     .             WK, XP,IERR)
        YP1 = YP(1)
        YPN = YP(NN)
        CALL YPC2 (NN,T,Y,SIGMA,IENDC,IENDC,YP1,YPN,
     .             WK, YP,IERR)
        IF (IERR .EQ. 1) GO TO 11
        IF (IERR .GT. 1) GO TO 16
      ELSE
C
C   Periodic fit:  call YPC2P.
C
        CALL YPC2P (NN,T,X,SIGMA,WK, XP,IERR)
        CALL YPC2P (NN,T,Y,SIGMA,WK, YP,IERR)
        IF (IERR .NE. 0) GO TO 16
      ENDIF
      LOOP2 = .FALSE.
C
C   Iterate on calls to SIGBP and YPC2 (or YPC2P).  The
C     first 2N-2 WK locations are used to store the deriva-
C     tive estimates XP and YP from the previous iteration.
C
C   LOOP2 is TRUE iff tension factors are not allowed to
C         decrease between iterations (loop 1 failed to
C         converge with MAXIT iterations).
C   DYP is the maximum relative change in a component of XP
C       or YP.
C   ICNT is the number of tension factors which were altered
C        by SIGBP.
C   DSMAX is the maximum relative change in a component of
C         SIGMA.
C
    2 DO 6 ITER = 1,MAXIT
        DYP = 0.D0
        DO 3 I = 2,NM1
          WK(I) = XP(I)
          WK(NM1+I) = YP(I)
    3     CONTINUE
        CALL SIGBP (NN,X,Y,XP,YP,STOL,BL,BU,
     .              BMAX, SIGMA, DSMAX,ICNT)
        IF (ICNT .LT. 0) GO TO 15
        IF (.NOT. PER) THEN
          CALL YPC2 (NN,T,X,SIGMA,IENDC,IENDC,XP1,XPN,
     .               WK(N2M1), XP,IERR)
          CALL YPC2 (NN,T,Y,SIGMA,IENDC,IENDC,YP1,YPN,
     .               WK(N2M1), YP,IERR)
        ELSE
          CALL YPC2P (NN,T,X,SIGMA,WK(N2M1), XP,IERR)
          CALL YPC2P (NN,T,Y,SIGMA,WK(N2M1), YP,IERR)
        ENDIF
        DO 4 I = 2,NM1
          EX = ABS(XP(I)-WK(I))
          IF (WK(I) .NE. 0.D0) EX = EX/ABS(WK(I))
          EY = ABS(YP(I)-WK(NM1+I))
          IF (WK(NM1+I) .NE. 0.D0) EY = EY/ABS(WK(NM1+I))
          DYP = MAX(DYP,EX,EY)
    4     CONTINUE
        IF (ICNT .EQ. 0  .OR.  DYP .LE. DYPTOL) GO TO 7
        IF (.NOT. LOOP2) THEN
C
C   Loop 1:  reinitialize SIGMA to zeros.
C
          DO 5 I = 1,NM1
            SIGMA(I) = 0.D0
    5       CONTINUE
        ENDIF
    6   CONTINUE
C
C   The loop failed to converge within MAXIT iterations.
C
      ITER = MAXIT
      IF (.NOT. LOOP2) THEN
        LOOP2 = .FALSE.
        GO TO 2
      ENDIF
C
C Store convergence parameters.
C
    7 WK(1) = DYP
      WK(2) = DSMAX
      IF (LOOP2) ITER = ITER + MAXIT
C
C No error encountered.
C
   10 IER = ITER
      RETURN
C
C Invalid input parameter N, NCD, or IENDC.
C
   11 IER = -1
      RETURN
C
C LWK too small.
C
   12 IER = -2
      RETURN
C
C Adjacent duplicate data points encountered.
C
   14 IER = -4
      RETURN
C
C Invalid constraint encountered by SIGBP.
C
   15 IER = -5
      RETURN
C
C Error flag returned by YPC1, YPC1P, YPC2, or YPC2P:
C   T is not strictly increasing.
C
   16 IER = -6
      RETURN
      END
      SUBROUTINE TSPSI (N,X,Y,NCD,IENDC,PER,UNIFRM,LWK, WK,
     .                  YP,SIGMA, IER)
      INTEGER N, NCD, IENDC, LWK, IER
      LOGICAL PER, UNIFRM
      DOUBLE PRECISION X(N), Y(N), WK(LWK), YP(N), SIGMA(N)
C
C***********************************************************
C
C                                                From TSPACK
C                                            Robert J. Renka
C                                  Dept. of Computer Science
C                                       Univ. of North Texas
C                                           renka@cs.unt.edu
C                                                   07/08/92
C
C   This subroutine computes a set of parameter values which
C define a Hermite interpolatory tension spline H(x).  The
C parameters consist of knot derivative values YP computed
C by Subroutine YPC1, YPC1P, YPC2, or YPC2P, and tension
C factors SIGMA computed by Subroutine SIGS (unless UNIFRM =
C TRUE).  Alternative methods for computing SIGMA are pro-
C vided by Subroutine TSPBI and Functions SIG0, SIG1, and
C SIG2.
C
C   Refer to Subroutine TSPSS for a means of computing
C parameters which define a smoothing curve rather than an
C interpolatory curve.
C
C   The tension spline may be evaluated by Subroutine TSVAL1
C or Functions HVAL (values), HPVAL (first derivatives),
C HPPVAL (second derivatives), and TSINTL (integrals).
C
C On input:
C
C       N = Number of data points.  N .GE. 2 and N .GE. 3 if
C           PER = TRUE.
C
C       X = Array of length N containing a strictly in-
C           creasing sequence of abscissae:  X(I) < X(I+1)
C           for I = 1,...,N-1.
C
C       Y = Array of length N containing data values asso-
C           ciated with the abscissae.  H(X(I)) = Y(I) for
C           I = 1,...,N.  If NCD = 1 and PER = TRUE, Y(N)
C           is set to Y(1).
C
C       NCD = Number of continuous derivatives at the knots.
C             NCD = 1 or NCD = 2.  If NCD = 1, the YP values
C             are computed by local monotonicity-constrained
C             quadratic fits.  Otherwise, a linear system is
C             solved for the derivative values which result
C             in second derivative continuity.  Unless
C             UNIFRM = TRUE, this requires iterating on
C             calls to YPC2 or YPC2P and calls to SIGS, and
C             generally results in more nonzero tension
C             factors (hence more expensive evaluation).
C
C       IENDC = End condition indicator for NCD = 2 and PER
C               = FALSE (or dummy parameter otherwise):
C               IENDC = 0 if YP(1) and YP(N) are to be com-
C                         puted by monotonicity-constrained
C                         parabolic fits to the first three
C                         and last three points, respective-
C                         ly.  This is identical to the
C                         values computed by YPC1.
C               IENDC = 1 if the first derivatives of H at
C                         X(1) and X(N) are user-specified
C                         in YP(1) and YP(N), respectively.
C               IENDC = 2 if the second derivatives of H at
C                         X(1) and X(N) are user-specified
C                         in YP(1) and YP(N), respectively.
C               IENDC = 3 if the end conditions are to be
C                         computed by Subroutine ENDSLP and
C                         vary with SIGMA(1) and SIGMA(N-1).
C
C       PER = Logical variable with value TRUE if and only
C             H(x) is to be a periodic function with period
C             X(N)-X(1).  It is assumed without a test that
C             Y(N) = Y(1) in this case.  On output, YP(N) =
C             YP(1).  If H(x) is one of the components of a
C             parametric curve, this option may be used to
C             obtained a closed curve.
C
C       UNIFRM = Logical variable with value TRUE if and
C                only if constant (uniform) tension is to be
C                used.  The tension factor must be input in
C                SIGMA(1) in this case and must be in the
C                range 0 to 85.  If SIGMA(1) = 0, H(x) is
C                piecewise cubic (a cubic spline if NCD =
C                2), and as SIGMA increases, H(x) approaches
C                the piecewise linear interpolant.  If
C                UNIFRM = FALSE, tension factors are chosen
C                (by SIGS) to preserve local monotonicity
C                and convexity of the data.  This often
C                improves the appearance of the curve over
C                the piecewise cubic fit.
C
C       LWK = Length of work space WK:  no work space is
C             needed if NCD = 1; at least N-1 locations
C             are required if NCD = 2; another N-1 locations
C             are required if PER = TRUE; and an additional
C             N-1 locations are required for the convergence
C             test if SIGS is called (UNIFRM = FALSE):
C
C             LWK GE 0    if NCD=1
C             LWK GE N-1  if NCD=2, PER=FALSE, UNIFRM=TRUE
C             LWK GE 2N-2 if NCD=2, PER=TRUE,  UNIFRM=TRUE
C             LWK GE 2N-2 if NCD=2, PER=FALSE, UNIFRM=FALSE
C             LWK GE 3N-3 if NCD=2, PER=TRUE,  UNIFRM=FALSE
C
C   The above parameters, except possibly Y(N), are not
C altered by this routine.
C
C       WK = Array of length at least LWK to be used as
C            temporary work space.
C
C       YP = Array of length .GE. N containing end condition
C            values in positions 1 and N if NCD = 2 and
C            IENDC = 1 or IENDC = 2.
C
C       SIGMA = Array of length .GE. N-1 containing a ten-
C               sion factor (0 to 85) in the first position
C               if UNIFRM = TRUE.
C
C On output:
C
C       WK = Array containing convergence parameters in the
C            first two locations if IER > 0 (NCD = 2 and
C            UNIFRM = FALSE):
C            WK(1) = Maximum relative change in a component
C                    of YP on the last iteration.
C            WK(2) = Maximum relative change in a component
C                    of SIGMA on the last iteration.
C
C       YP = Array containing derivatives of H at the
C            abscissae.  YP is not altered if -4 < IER < 0,
C            and YP is only partially defined if IER = -4.
C
C       SIGMA = Array containing tension factors.  SIGMA(I)
C               is associated with interval (X(I),X(I+1))
C               for I = 1,...,N-1.  SIGMA is not altered if
C               -4 < IER < 0 (unless IENDC is invalid), and
C               SIGMA is constant (not optimal) if IER = -4
C               or IENDC (if used) is invalid.
C
C       IER = Error indicator or iteration count:
C             IER = IC .GE. 0 if no errors were encountered
C                      and IC calls to SIGS and IC+1 calls
C                      to YPC1, YPC1P, YPC2 or YPC2P were
C                      employed.  (IC = 0 if NCD = 1).
C             IER = -1 if N, NCD, or IENDC is outside its
C                      valid range.
C             IER = -2 if LWK is too small.
C             IER = -3 if UNIFRM = TRUE and SIGMA(1) is out-
C                      side its valid range.
C             IER = -4 if the abscissae X are not strictly
C                      increasing.
C
C Modules required by TSPSI:  ENDSLP, SIGS, SNHCSH, STORE,
C                               YPCOEF, YPC1, YPC1P, YPC2,
C                               YPC2P
C
C Intrinsic functions called by TSPSI:  ABS, MAX
C
C***********************************************************
C
      INTEGER I, ICNT, IERR, ITER, MAXIT, NM1, NN
      DOUBLE PRECISION DSMAX, DYP, DYPTOL, E, SBIG, SIG,
     .                 STOL, YP1, YPN
C
      DATA SBIG/85.D0/
C
C Convergence parameters:
C
C   STOL = Absolute tolerance for SIGS
C   MAXIT = Maximum number of YPC2/SIGS iterations
C   DYPTOL = Bound on the maximum relative change in a
C              component of YP defining convergence of
C              the YPC2/SIGS iteration when NCD = 2 and
C              UNIFRM = FALSE
C
      DATA STOL/0.D0/,  MAXIT/99/,  DYPTOL/.01D0/
C
C Test for invalid input parameters (other than X and
C   IENDC).
C
      NN = N
      NM1 = NN - 1
      IF (NN .LT. 2  .OR.  (PER  .AND.  NN .LT. 3)  .OR.
     .    NCD .LT. 1  .OR.  NCD .GT. 2) GO TO 11
      IF (UNIFRM) THEN
        IF ( NCD .EQ. 2  .AND.  (LWK .LT. NM1  .OR.
     .       (PER  .AND.  LWK .LT. 2*NM1)) ) GO TO 12
        SIG = SIGMA(1)
        IF (SIG .LT. 0.D0  .OR.  SIG .GT. SBIG) GO TO 13
      ELSE
        IF ( NCD .EQ. 2  .AND.  (LWK .LT. 2*NM1  .OR.
     .       (PER  .AND.  LWK .LT. 3*NM1)) ) GO TO 12
        SIG = 0.D0
      ENDIF
C
C Initialize iteration count ITER, and store uniform
C   tension factors, or initialize SIGMA to zeros.
C
      ITER = 0
      DO 1 I = 1,NM1
        SIGMA(I) = SIG
    1   CONTINUE
      IF (NCD .EQ. 1) THEN
C
C NCD = 1.
C
        IF (.NOT. PER) THEN
          CALL YPC1 (NN,X,Y, YP,IERR)
        ELSE
          CALL YPC1P (NN,X,Y, YP,IERR)
        ENDIF
        IF (IERR .NE. 0) GO TO 14
        IF (.NOT. UNIFRM) THEN
C
C   Call SIGS for UNIFRM = FALSE.
C
          CALL SIGS (NN,X,Y,YP,STOL, SIGMA, DSMAX,IERR)
        ENDIF
        GO TO 10
      ENDIF
C
C NCD = 2.
C
      IF (.NOT. PER) THEN
C
C   Nonperiodic case:  call YPC2 and test for IENDC or X
C     invalid.
C
        YP1 = YP(1)
        YPN = YP(NN)
        CALL YPC2 (NN,X,Y,SIGMA,IENDC,IENDC,YP1,YPN,
     .             WK, YP,IERR)
        IF (IERR .EQ. 1) GO TO 11
        IF (IERR .GT. 1) GO TO 14
      ELSE
C
C   Periodic fit:  call YPC2P.
C
        CALL YPC2P (NN,X,Y,SIGMA,WK, YP,IERR)
        IF (IERR .GT. 1) GO TO 14
      ENDIF
      IF (UNIFRM) GO TO 10
C
C   Iterate on calls to SIGS and YPC2 (or YPC2P).  The first
C     N-1 WK locations are used to store the derivative
C     estimates YP from the previous iteration.
C
C   DYP is the maximum relative change in a component of YP.
C   ICNT is the number of tension factors which were
C        increased by SIGS.
C   DSMAX is the maximum relative change in a component of
C         SIGMA.
C
      DO 4 ITER = 1,MAXIT
        DYP = 0.D0
        DO 2 I = 2,NM1
          WK(I) = YP(I)
    2     CONTINUE
        CALL SIGS (NN,X,Y,YP,STOL, SIGMA, DSMAX,ICNT)
        IF (.NOT. PER) THEN
          CALL YPC2 (NN,X,Y,SIGMA,IENDC,IENDC,YP1,YPN,
     .               WK(NN), YP,IERR)
        ELSE
          CALL YPC2P (NN,X,Y,SIGMA,WK(NN), YP,IERR)
        ENDIF
        DO 3 I = 2,NM1
          E = ABS(YP(I)-WK(I))
          IF (WK(I) .NE. 0.D0) E = E/ABS(WK(I))
          DYP = MAX(DYP,E)
    3     CONTINUE
        IF (ICNT .EQ. 0  .OR.  DYP .LE. DYPTOL) GO TO 5
    4   CONTINUE
      ITER = MAXIT
C
C Store convergence parameters in WK.
C
    5 WK(1) = DYP
      WK(2) = DSMAX
C
C No error encountered.
C
   10 IER = ITER
      RETURN
C
C Invalid input parameter N, NCD, or IENDC.
C
   11 IER = -1
      RETURN
C
C LWK too small.
C
   12 IER = -2
      RETURN
C
C UNIFRM = TRUE and SIGMA(1) outside its valid range.
C
   13 IER = -3
      RETURN
C
C Abscissae are not strictly increasing.
C
   14 IER = -4
      RETURN
      END
      SUBROUTINE TSPSP (N,ND,X,Y,Z,NCD,IENDC,PER,UNIFRM,
     .                  LWK, WK, T,XP,YP,ZP,SIGMA,IER)
      INTEGER N, ND, NCD, IENDC, LWK, IER
      LOGICAL PER, UNIFRM
      DOUBLE PRECISION X(N), Y(N), Z(N), WK(LWK), T(N),
     .                 XP(N), YP(N), ZP(N), SIGMA(N)
C
C***********************************************************
C
C                                                From TSPACK
C                                            Robert J. Renka
C                                  Dept. of Computer Science
C                                       Univ. of North Texas
C                                           renka@cs.unt.edu
C                                                   07/08/92
C
C   This subroutine computes a set of values which define a
C parametric planar curve C(t) = (H1(t),H2(t)) or space
C curve C(t) = (H1(t),H2(t),H3(t)) whose components are Her-
C mite interpolatory tension splines.  The output values
C consist of parameters (knots) T computed by ARCL2D or
C ARCL3D, knot derivative values XP, YP, (and ZP) computed
C by Subroutine YPC1, YPC1P, YPC2, or YPC2P, and tension
C factors SIGMA computed by Subroutine SIGS (unless UNIFRM =
C TRUE).
C
C   Refer to Subroutine TSPSP for an alternative method of
C computing tension factors in the case of a planar curve.
C
C   The tension splines may be evaluated by Subroutine
C TSVAL2 (or TSVAL3) or Functions HVAL (values), HPVAL
C (first derivatives), HPPVAL (second derivatives), and
C TSINTL (integrals).
C
C On input:
C
C       N = Number of knots and data points.  N .GE. 2 and
C           N .GE. 3 if PER = TRUE.
C
C       ND = Number of dimensions:
C            ND = 2 if a planar curve is to be constructed.
C            ND = 3 if a space curve is to be constructed.
C
C       X,Y,Z = Arrays of length N containing the Cartesian
C               coordinates of an ordered sequence of data
C               points C(I), I = 1 to N, such that C(I) .NE.
C               C(I+1).  C(t) is constrained to pass through
C               these points.  Z is an unused dummy parame-
C               ter if ND = 2.  In the case of a closed curve
C               (PER = TRUE), the first and last points
C               should coincide.  In this case, X(N), Y(N),
C               (and Z(N)) are set to X(1), Y(1), (and Z(1))
C               if NCD = 1, but not altered if NCD = 2.
C
C       NCD = Number of continuous derivatives at the knots.
C             NCD = 1 or NCD = 2.  If NCD = 1, XP, YP, (and
C             ZP) are computed by local monotonicity-
C             constrained quadratic fits.  Otherwise, a
C             linear system is solved for the derivative
C             values which result in second derivative con-
C             tinuity.  Unless UNIFRM = FALSE, this requires
C             iterating on calls to YPC2 or YPC2P and calls
C             to SIGS, and generally results in more nonzero
C             tension factors (hence more expensive evalua-
C             tion).
C
C       IENDC = End condition indicator for NCD = 2 and PER
C               = FALSE (or dummy parameter otherwise):
C               IENDC = 0 if XP(1), XP(N), YP(1), YP(N) (and
C                         ZP(1) and ZP(N)) are to be com-
C                         puted by monotonicity-constrained
C                         parabolic fits (YPC1).
C               IENDC = 1 if the first derivatives of H1 at
C                         the left and right endpoints are
C                         user-specified in XP(1) and XP(N),
C                         respectively, the first deriva-
C                         tives of H2 at the ends are
C                         specified in YP(1) and YP(N), and,
C                         if ND = 3, the first derivatives
C                         of H3 are specified in ZP(1) and
C                         ZP(N).
C               IENDC = 2 if the second derivatives of H1,
C                         H2, (and H3) at the endpoints are
C                         user-specified in XP(1), XP(N),
C                         YP(1), YP(N), (ZP(1), and ZP(N)).
C               IENDC = 3 if the end conditions are to be
C                         computed by Subroutine ENDSLP and
C                         vary with SIGMA(1) and SIGMA(N-1).
C
C       PER = Logical variable with value TRUE if and only
C             a closed curve is to be constructed -- H1(t),
C             H2(t), (and H3(t)) are to be periodic func-
C             tions with period T(N)-T(1), where T(1) and
C             T(N) are the parameter values associated with
C             the first and last data points.  It is assumed
C             in this case that X(N) = X(1), Y(N) = Y(1)
C             and, if ND = 3, Z(N) = Z(1), and, on output,
C             XP(N) = XP(1), YP(N) = YP(1), (and ZP(N) =
C             ZP(1) if ND = 3).
C
C       UNIFRM = Logical variable with value TRUE if and
C                only if constant (uniform) tension is to be
C                used.  The tension factor must be input in
C                SIGMA(1) in this case and must be in the
C                range 0 to 85.  If SIGMA(1) = 0, H(t) is
C                piecewise cubic (a cubic spline if NCD =
C                2), and as SIGMA increases, H(t) approaches
C                the piecewise linear interpolant, where H
C                is H1, H2, or H3.  If UNIFRM = FALSE,
C                tension factors are chosen (by SIGS) to
C                preserve local monotonicity and convexity
C                of the data.  This often improves the
C                appearance of the curve over the piecewise
C                cubic fitting functions.
C
C       LWK = Length of work space WK:  no work space is
C             needed if NCD = 1; at least N-1 locations
C             are required if NCD = 2; another N-1 locations
C             are required if PER = TRUE; and an additional
C             ND*(N-1) locations are required for the con-
C             vergence test if SIGS is called (UNIFRM =
C             FALSE):
C               If NCD=1 then LWK = 0 (not tested).
C               If NCD=2 then
C
C             LWK GE N-1          if PER=FALSE, UNIFRM=TRUE
C             LWK GE 2N-2         if PER=TRUE,  UNIFRM=TRUE
C             LWK GE (ND+1)*(N-1) if PER=FALSE, UNIFRM=FALSE
C             LWK GE (ND+2)*(N-1) if PER=TRUE,  UNIFRM=FALSE
C
C   The above parameters, except possibly X(N), Y(N), and
C Z(N), are not altered by this routine.
C
C       WK = Array of length .GE. LWK to be used as tempor-
C            ary work space.
C
C       T = Array of length .GE. N.
C
C       XP = Array of length .GE. N containing end condition
C            values in positions 1 and N if NCD = 2 and
C            IENDC = 1 or IENDC = 2.
C
C       YP = Array of length .GE. N containing end condition
C            values in positions 1 and N if NCD = 2 and
C            IENDC = 1 or IENDC = 2.
C
C       ZP = Dummy argument if ND=2, or, if ND=3, array of
C            length .GE. N containing end condition values
C            in positions 1 and N if NCD = 2 and IENDC = 1
C            or IENDC = 2.
C
C       SIGMA = Array of length .GE. N-1 containing a ten-
C               sion factor (0 to 85) in the first position
C               if UNIFRM = TRUE.
C
C On output:
C
C       WK = Array containing convergence parameters in the
C            first two locations if IER > 0 (NCD = 2 and
C            UNIFRM = FALSE):
C            WK(1) = Maximum relative change in a component
C                    of XP, YP, or ZP on the last iteration.
C            WK(2) = Maximum relative change in a component
C                    of SIGMA on the last iteration.
C
C       T = Array containing parameter values computed by
C           Subroutine ARCL2D or ARCL3D unless -4 < IER < 0.
C           T is only partially defined if IER = -4.
C
C       XP = Array containing derivatives of H1 at the
C            knots.  XP is not altered if -5 < IER < 0,
C            and XP is only partially defined if IER = -6.
C
C       YP = Array containing derivatives of H2 at the
C            knots.  YP is not altered if -5 < IER < 0,
C            and YP is only partially defined if IER = -6.
C
C       ZP = Array containing derivatives of H3 at the knots
C            if ND=3.  ZP is not altered if -5 < IER < 0,
C            and ZP is only partially defined if IER = -6.
C
C       SIGMA = Array containing tension factors.  SIGMA(I)
C               is associated with interval (T(I),T(I+1))
C               for I = 1,...,N-1.  SIGMA is not altered if
C               -5 < IER < 0 (unless IENDC is invalid), and
C               SIGMA is constant (not optimal) if IER = -6
C               or IENDC (if used) is invalid.
C
C       IER = Error indicator or iteration count:
C             IER = IC .GE. 0 if no errors were encountered
C                      and IC calls to SIGS and IC+1 calls
C                      to YPC1, YPC1P, YPC2 or YPC2P were
C                      employed.  (IC = 0 if NCD = 1).
C             IER = -1 if N, NCD, or IENDC is outside its
C                      valid range.
C             IER = -2 if LWK is too small.
C             IER = -3 if UNIFRM = TRUE and SIGMA(1) is out-
C                      side its valid range.
C             IER = -4 if a pair of adjacent data points
C                      coincide:  X(I) = X(I+1), Y(I) =
C                      Y(I+1), (and Z(I) = Z(I+1)) for some
C                      I in the range 1 to N-1.
C             IER = -6 if invalid knots T were returned by
C                      ARCL2D or ARCL3D.  This should not
C                      occur.
C
C Modules required by TSPSP:  ARCL2D, ARCL3D, ENDSLP, SIGS,
C                               SNHCSH, STORE, YPCOEF, YPC1,
C                               YPC1P, YPC2, YPC2P
C
C Intrinsic functions called by TSPSP:  ABS, MAX
C
C***********************************************************
C
      INTEGER I, ICNT, IERR, ITER, IW1, MAXIT, N2M2, NM1, NN
      LOGICAL SCURV
      DOUBLE PRECISION DSMAX, DYP, DYPTOL, EX, EY, EZ, SBIG,
     .                 SIG, STOL, XP1, XPN, YP1, YPN, ZP1,
     .                 ZPN
C
      DATA SBIG/85.D0/
C
C Convergence parameters:
C
C   STOL = Absolute tolerance for SIGS
C   MAXIT = Maximum number of YPC2/SIGS iterations
C   DYPTOL = Bound on the maximum relative change in a com-
C              ponent of XP, YP, or ZP defining convergence
C              of the YPC2/SIGS iteration when NCD = 2 and
C              UNIFRM = FALSE
C
      DATA STOL/0.D0/,  MAXIT/99/,  DYPTOL/.01D0/
C
C Test for invalid input parameters N, NCD, or LWK.
C
      NN = N
      NM1 = NN - 1
      N2M2 = 2*NM1
      IF (NN .LT. 2  .OR.  (PER  .AND.  NN .LT. 3)  .OR.
     .    NCD .LT. 1  .OR.  NCD .GT. 2) GO TO 11
      IF (UNIFRM) THEN
        IF ( NCD .EQ. 2  .AND.  (LWK .LT. NM1  .OR.
     .       (PER  .AND.  LWK .LT. N2M2)) ) GO TO 12
        SIG = SIGMA(1)
        IF (SIG .LT. 0.D0  .OR.  SIG .GT. SBIG) GO TO 13
      ELSE
        IF ( NCD .EQ. 2  .AND.  (LWK .LT. (ND+1)*NM1  .OR.
     .       (PER  .AND.  LWK .LT. (ND+2)*NM1)) ) GO TO 12
        SIG = 0.D0
      ENDIF
C
C Compute the sequence of parameter values T.
C
      SCURV = ND .EQ. 3
      IF (.NOT. SCURV) THEN
        CALL ARCL2D (NN,X,Y, T,IERR)
      ELSE
        CALL ARCL3D (NN,X,Y,Z, T,IERR)
      ENDIF
      IF (IERR .GT. 0) GO TO 14
C
C Initialize iteration count ITER, and store uniform
C   tension factors, or initialize SIGMA to zeros.
C
      ITER = 0
      DO 1 I = 1,NM1
        SIGMA(I) = SIG
    1   CONTINUE
      IF (NCD .EQ. 1) THEN
C
C NCD = 1.
C
        IF (.NOT. PER) THEN
          CALL YPC1 (NN,T,X, XP,IERR)
          CALL YPC1 (NN,T,Y, YP,IERR)
          IF (SCURV) CALL YPC1 (NN,T,Z, ZP,IERR)
        ELSE
          CALL YPC1P (NN,T,X, XP,IERR)
          CALL YPC1P (NN,T,Y, YP,IERR)
          IF (SCURV) CALL YPC1P (NN,T,Z, ZP,IERR)
        ENDIF
        IF (IERR .NE. 0) GO TO 16
        IF (.NOT. UNIFRM) THEN
C
C   Call SIGS for UNIFRM = FALSE.
C
          CALL SIGS (NN,T,X,XP,STOL, SIGMA, DSMAX,IERR)
          CALL SIGS (NN,T,Y,YP,STOL, SIGMA, DSMAX,IERR)
          IF (SCURV) CALL SIGS (NN,T,Z,ZP,STOL, SIGMA,
     .                           DSMAX,IERR)
        ENDIF
        GO TO 10
      ENDIF
C
C NCD = 2.
C
      IF (.NOT. PER) THEN
C
C   Nonperiodic case:  call YPC2 and test for IENDC invalid.
C
        XP1 = XP(1)
        XPN = XP(NN)
        CALL YPC2 (NN,T,X,SIGMA,IENDC,IENDC,XP1,XPN,
     .             WK, XP,IERR)
        YP1 = YP(1)
        YPN = YP(NN)
        CALL YPC2 (NN,T,Y,SIGMA,IENDC,IENDC,YP1,YPN,
     .             WK, YP,IERR)
        IF (SCURV) THEN
          ZP1 = ZP(1)
          ZPN = ZP(NN)
          CALL YPC2 (NN,T,Z,SIGMA,IENDC,IENDC,ZP1,ZPN,
     .               WK, ZP,IERR)
        ENDIF
        IF (IERR .EQ. 1) GO TO 11
        IF (IERR .GT. 1) GO TO 16
      ELSE
C
C   Periodic fit:  call YPC2P.
C
        CALL YPC2P (NN,T,X,SIGMA,WK, XP,IERR)
        CALL YPC2P (NN,T,Y,SIGMA,WK, YP,IERR)
        IF (SCURV) CALL YPC2P (NN,T,Z,SIGMA,WK, ZP,IERR)
        IF (IERR .NE. 0) GO TO 16
      ENDIF
      IF (UNIFRM) GO TO 10
C
C   Iterate on calls to SIGS and YPC2 (or YPC2P).  The
C     first ND*(N-1) WK locations are used to store the
C     derivative estimates XP, YP, (and ZP) from the
C     previous iteration.  IW1 is the first free WK location
C     following the stored derivatives.
C
C   DYP is the maximum relative change in a component of XP,
C       YP, or ZP.
C   ICNT is the number of tension factors which were
C        increased by SIGS.
C   DSMAX is the maximum relative change in a component of
C         SIGMA.
C
      IW1 = ND*NM1 + 1
      DO 5 ITER = 1,MAXIT
        DYP = 0.D0
        DO 2 I = 2,NM1
          WK(I) = XP(I)
          WK(NM1+I) = YP(I)
    2     CONTINUE
        IF (SCURV) THEN
          DO 3 I = 2,NM1
            WK(N2M2+I) = ZP(I)
    3       CONTINUE
        ENDIF
        CALL SIGS (NN,T,X,XP,STOL, SIGMA, DSMAX,ICNT)
        CALL SIGS (NN,T,Y,YP,STOL, SIGMA, DSMAX,ICNT)
        IF (SCURV) CALL SIGS (NN,T,Z,ZP,STOL, SIGMA, DSMAX,
     .                        ICNT)
        IF (.NOT. PER) THEN
          CALL YPC2 (NN,T,X,SIGMA,IENDC,IENDC,XP1,XPN,
     .               WK(IW1), XP,IERR)
          CALL YPC2 (NN,T,Y,SIGMA,IENDC,IENDC,YP1,YPN,
     .               WK(IW1), YP,IERR)
          IF (SCURV) CALL YPC2 (NN,T,Z,SIGMA,IENDC,IENDC,
     .                          ZP1,ZPN,WK(IW1), ZP,IERR)
        ELSE
          CALL YPC2P (NN,T,X,SIGMA,WK(IW1), XP,IERR)
          CALL YPC2P (NN,T,Y,SIGMA,WK(IW1), YP,IERR)
          IF (SCURV) CALL YPC2P (NN,T,Z,SIGMA,WK(IW1), ZP,
     .                           IERR)
        ENDIF
        EZ = 0.D0
        DO 4 I = 2,NM1
          EX = ABS(XP(I)-WK(I))
          IF (WK(I) .NE. 0.D0) EX = EX/ABS(WK(I))
          EY = ABS(YP(I)-WK(NM1+I))
          IF (WK(NM1+I) .NE. 0.D0) EY = EY/ABS(WK(NM1+I))
          IF (SCURV) THEN
            EZ = ABS(ZP(I)-WK(N2M2+I))
            IF (WK(N2M2+I) .NE. 0.D0)
     .        EZ = EZ/ABS(WK(N2M2+I))
          ENDIF
          DYP = MAX(DYP,EX,EY,EZ)
    4     CONTINUE
        IF (ICNT .EQ. 0  .OR.  DYP .LE. DYPTOL) GO TO 6
    5   CONTINUE
      ITER = MAXIT
C
C Store convergence parameters in WK.
C
    6 WK(1) = DYP
      WK(2) = DSMAX
C
C No error encountered.
C
   10 IER = ITER
      RETURN
C
C Invalid input parameter N, NCD, or IENDC.
C
   11 IER = -1
      RETURN
C
C LWK too small.
C
   12 IER = -2
      RETURN
C
C UNIFRM = TRUE and SIGMA(1) outside its valid range.
C
   13 IER = -3
      RETURN
C
C Adjacent duplicate data points encountered.
C
   14 IER = -4
      RETURN
C
C Error flag returned by YPC1, YPC1P, YPC2, or YPC2P:
C   T is not strictly increasing.
C
   16 IER = -6
      RETURN
      END
      SUBROUTINE TSPSS (N,X,Y,PER,UNIFRM,W,SM,SMTOL,LWK, WK,
     .                  SIGMA,YS,YP, NIT,IER)
      INTEGER N, LWK, NIT, IER
      LOGICAL PER, UNIFRM
      DOUBLE PRECISION X(N), Y(N), W(N), SM, SMTOL, WK(LWK),
     .                 SIGMA(N), YS(N), YP(N)
C
C***********************************************************
C
C                                                From TSPACK
C                                            Robert J. Renka
C                                  Dept. of Computer Science
C                                       Univ. of North Texas
C                                           renka@cs.unt.edu
C                                                   06/10/92
C
C   This subroutine computes a set of parameter values which
C define a smoothing tension spline H(x).  The parameters
C consist of knot values YS and derivatives YP computed
C by Subroutine SMCRV, and tension factors SIGMA computed by
C Subroutine SIGS (unless UNIFRM = TRUE).  The Hermite
C interpolatory tension spline H(x) defined by the knot
C values and derivatives has two continuous derivatives and
C satisfies either natural or periodic end conditions.
C
C   The tension spline may be evaluated by Subroutine TSVAL1
C or Functions HVAL (values), HPVAL (first derivatives),
C HPPVAL (second derivatives), and TSINTL (integrals).
C
C On input:
C
C       N = Number of data points.  N .GE. 2 and N .GE. 3 if
C           PER = TRUE.
C
C       X = Array of length N containing a strictly in-
C           creasing sequence of abscissae:  X(I) < X(I+1)
C           for I = 1,...,N-1.
C
C       Y = Array of length N containing data values asso-
C           ciated with the abscissae.  If PER = TRUE, it is
C           assumed that Y(N) = Y(1).
C
C       PER = Logical variable with value TRUE if and only
C             H(x) is to be a periodic function with period
C             X(N)-X(1).  It is assumed without a test that
C             Y(N) = Y(1) in this case.  On output, YP(N) =
C             YP(1) and, more generally, the values and
C             first two derivatives of H at X(1) agree with
C             those at X(N).  If H(x) is one of the compo-
C             nents of a parametric curve, this option may
C             be used to obtained a closed curve.  If PER =
C             FALSE, H satisfies natural end conditions:
C             zero second derivatives at X(1) and X(N).
C
C       UNIFRM = Logical variable with value TRUE if and
C                only if constant (uniform) tension is to be
C                used.  The tension factor must be input in
C                SIGMA(1) in this case and must be in the
C                range 0 to 85.  If SIGMA(1) = 0, H(x) is
C                a cubic spline, and as SIGMA increases,
C                H(x) approaches piecewise linear.  If
C                UNIFRM = FALSE, tension factors are chosen
C                (by SIGS) to preserve local monotonicity
C                and convexity of the data.  This may re-
C                sult in a better fit than the case of
C                uniform tension, but requires an iteration
C                on calls to SMCRV and SIGS.
C
C       W = Array of length N containing positive weights
C           associated with the data values.  The recommend-
C           ed value of W(I) is 1/DY**2, where DY is the
C           standard deviation associated with Y(I).  If
C           nothing is known about the errors in Y, a con-
C           stant (estimated value) should be used for DY.
C           If PER = TRUE, it is assumed that W(N) = W(1).
C
C       SM = Positive parameter specifying an upper bound on
C            Q2(YS), where Q2(YS) is the weighted sum of
C            squares of deviations from the data (differ-
C            ences between YS and Y).  H(x) is linear (and
C            Q2 is minimized) if SM is sufficiently large
C            that the constraint is not active.  It is
C            recommended that SM satisfy N-SQRT(2N) .LE. SM
C            .LE. N+SQRT(2N) and SM = N is reasonable if
C            W(I) = 1/DY**2.
C
C       SMTOL = Parameter in the range (0,1) specifying the
C               relative error allowed in satisfying the
C               constraint:  the constraint is assumed to
C               be satisfied if SM*(1-SMTOL) .LE. Q2 .LE.
C               SM*(1+SMTOL).  A reasonable value for SMTOL
C               is SQRT(2/N).
C
C       LWK = Length of work space WK:
C             LWK .GE. 6N   if PER=FALSE  and  UNIFRM=TRUE
C             LWK .GE. 7N   if PER=FALSE  and  UNIFRM=FALSE
C             LWK .GE. 10N  if PER=TRUE   and  UNIFRM=TRUE
C             LWK .GE. 11N  if PER=TRUE   and  UNIFRM=FALSE
C
C The above parameters are not altered by this routine.
C
C       WK = Array of length at least LWK to be used as
C            temporary work space.
C
C       SIGMA = Array of length .GE. N-1 containing a ten-
C               sion factor (0 to 85) in the first position
C               if UNIFRM = TRUE.
C
C       YS = Array of length .GE. N.
C
C       YP = Array of length .GE. N.
C
C On output:
C
C       WK = Array containing convergence parameters in the
C            first two locations if NIT > 0:
C            WK(1) = Maximum relative change in a component
C                    of YS on the last iteration.
C            WK(2) = Maximum relative change in a component
C                    of SIGMA on the last iteration.
C
C       SIGMA = Array containing tension factors.  SIGMA(I)
C               is associated with interval (X(I),X(I+1))
C               for I = 1,...,N-1.  SIGMA is not altered if
C               N is invalid or -4 < IER < -1, and SIGMA is
C               constant if IER = -1 (and N is valid) or
C               IER = -4.
C
C       YS = Array of length N containing values of H at the
C            abscissae.  YS(N) = YS(1) if PER = TRUE.  YS is
C            not altered if IER < 0.
C
C       YP = Array of length N containing first derivative
C            values of H at the abscissae.  YP(N) = YP(1)
C            if PER = TRUE.  YP is not altered if IER < 0.
C
C       NIT = Number of iterations (calls to SIGS).  NIT = 0
C             if IER < 0 or UNIFRM = TRUE.  If NIT > 0,
C             NIT+1 calls to SMCRV were employed.
C
C       IER = Error indicator:
C             IER = 0 if no errors were encountered and the
C                     constraint is active:  Q2(YS) is ap-
C                     proximately equal to SM.
C             IER = 1 if no errors were encountered but the
C                     constraint is not active:  YS and YP
C                     are the values and derivatives of the
C                     linear function (constant function if
C                     PERIOD = TRUE) which minimizes Q2, and
C                     Q1 = 0 (refer to SMCRV).
C             IER = -1 if N, W, SM, or SMTOL is outside its
C                      valid range.
C             IER = -2 if LWK is too small.
C             IER = -3 if UNIFRM = TRUE and SIGMA(1) is out-
C                      side its valid range.
C             IER = -4 if the abscissae X are not strictly
C                      increasing.
C
C Modules required by TSPSS:  B2TRI or B2TRIP, SIGS, SMCRV,
C                               SNHCSH, STORE, YPCOEF
C
C Intrinsic functions called by TSPSS:  ABS, MAX
C
C***********************************************************
C
      INTEGER I, ICNT, IERR, ITER, MAXIT, NM1, NN
      DOUBLE PRECISION DSMAX, DYS, DYSTOL, E, SBIG, SIG,
     .                 STOL
C
      DATA SBIG/85.D0/
C
C Convergence parameters:
C
C   STOL = Absolute tolerance for SIGS
C   MAXIT = Maximum number of SMCRV/SIGS iterations
C   DYSTOL = Bound on the maximum relative change in a
C              component of YS defining convergence of
C              the SMCRV/SIGS iteration when UNIFRM = FALSE
C
      DATA STOL/0.D0/,  MAXIT/99/,  DYSTOL/.01D0/
C
C Initialize NIT, and test for invalid input parameters LWK
C   and SIGMA(1).
C
      NIT = 0
      NN = N
      NM1 = NN - 1
      IF (NN .LT. 2  .OR.  (PER  .AND.  NN .LT. 3)) GO TO 11
      IF (UNIFRM) THEN
        IF ( LWK .LT. 6*NN  .OR.
     .       (PER  .AND.  LWK .LT. 10*NN) ) GO TO 12
        SIG = SIGMA(1)
        IF (SIG .LT. 0.D0  .OR.  SIG .GT. SBIG) GO TO 13
      ELSE
        IF ( LWK .LT. 7*NN  .OR.
     .       (PER  .AND.  LWK .LT. 11*NN) ) GO TO 12
        SIG = 0.D0
      ENDIF
C
C Store uniform tension factors, or initialize SIGMA to
C   zeros.
C
      DO 1 I = 1,NM1
        SIGMA(I) = SIG
    1   CONTINUE
C
C Compute smoothing curve for uniform tension.
C
      CALL SMCRV (NN,X,Y,SIGMA,PER,W,SM,SMTOL,WK, YS,YP,IER)
      IF (IER .LE. -2) IER = -4
      IF (IER .LT. 0  .OR.  UNIFRM) RETURN
C
C   Iterate on calls to SIGS and SMCRV.  The first N-1 WK
C     locations are used to store the function values YS
C     from the previous iteration.
C
C   DYS is the maximum relative change in a component of YS.
C   ICNT is the number of tension factors which were
C        increased by SIGS.
C   DSMAX is the maximum relative change in a component of
C         SIGMA.
C
      DO 4 ITER = 1,MAXIT
        DYS = 0.D0
        DO 2 I = 2,NM1
          WK(I) = YS(I)
    2     CONTINUE
        CALL SIGS (NN,X,Y,YP,STOL, SIGMA, DSMAX,ICNT)
        CALL SMCRV (NN,X,Y,SIGMA,PER,W,SM,SMTOL,WK(NN), YS,
     .              YP,IERR)
        DO 3 I = 2,NM1
          E = ABS(YS(I)-WK(I))
          IF (WK(I) .NE. 0.D0) E = E/ABS(WK(I))
          DYS = MAX(DYS,E)
    3     CONTINUE
        IF (ICNT .EQ. 0  .OR.  DYS .LE. DYSTOL) GO TO 5
    4   CONTINUE
      ITER = MAXIT
C
C No error encountered.
C
    5 WK(1) = DYS
      WK(2) = DSMAX
      NIT = ITER
      IER = IERR
      RETURN
C
C Invalid input parameter N, W, SM, or SMTOL.
C
   11 IER = -1
      RETURN
C
C LWK too small.
C
   12 IER = -2
      RETURN
C
C UNIFRM = TRUE and SIGMA(1) outside its valid range.
C
   13 IER = -3
      RETURN
      END
      SUBROUTINE TSVAL1 (N,X,Y,YP,SIGMA,IFLAG,NE,TE, V,
     .                   IER)
      INTEGER N, IFLAG, NE, IER
      DOUBLE PRECISION X(N), Y(N), YP(N), SIGMA(N), TE(NE),
     .                 V(NE)
C
C***********************************************************
C
C                                                From TSPACK
C                                            Robert J. Renka
C                                  Dept. of Computer Science
C                                       Univ. of North Texas
C                                           renka@cs.unt.edu
C                                                   06/10/92
C
C   This subroutine evaluates a Hermite interpolatory ten-
C sion spline H or its first or second derivative at a set
C of points TE.
C
C   Note that a large tension factor in SIGMA may cause
C underflow.  The result is assumed to be zero.  If not the
C default, this may be specified by either a compiler option
C or operating system option.
C
C On input:
C
C       N = Number of data points.  N .GE. 2.
C
C       X = Array of length N containing the abscissae.
C           These must be in strictly increasing order:
C           X(I) < X(I+1) for I = 1,...,N-1.
C
C       Y = Array of length N containing data values or
C           function values returned by Subroutine SMCRV.
C           Y(I) = H(X(I)) for I = 1,...,N.
C
C       YP = Array of length N containing first deriva-
C            tives.  YP(I) = HP(X(I)) for I = 1,...,N, where
C            HP denotes the derivative of H.
C
C       SIGMA = Array of length N-1 containing tension fac-
C               tors whose absolute values determine the
C               balance between cubic and linear in each
C               interval.  SIGMA(I) is associated with int-
C               erval (I,I+1) for I = 1,...,N-1.
C
C       IFLAG = Output option indicator:
C               IFLAG = 0 if values of H are to be computed.
C               IFLAG = 1 if first derivative values are to
C                         be computed.
C               IFLAG = 2 if second derivative values are to
C                         be computed.
C
C       NE = Number of evaluation points.  NE > 0.
C
C       TE = Array of length NE containing the evaluation
C            points.  The sequence should be strictly in-
C            creasing for maximum efficiency.  Extrapolation
C            is performed if a point is not in the interval
C            [X(1),X(N)].
C
C The above parameters are not altered by this routine.
C
C       V = Array of length at least NE.
C
C On output:
C
C       V = Array of function, first derivative, or second
C           derivative values at the evaluation points un-
C           less IER < 0.  If IER = -1, V is not altered.
C           If IER = -2, V may be only partially defined.
C
C       IER = Error indicator:
C             IER = 0  if no errors were encountered and
C                      no extrapolation occurred.
C             IER > 0  if no errors were encountered but
C                      extrapolation was required at IER
C                      points.
C             IER = -1 if N < 2, IFLAG < 0, IFLAG > 2, or
C                      NE < 1.
C             IER = -2 if the abscissae are not in strictly
C                      increasing order.  (This error will
C                      not necessarily be detected.)
C
C Modules required by TSVAL1:  HPPVAL, HPVAL, HVAL, INTRVL,
C                                SNHCSH
C
C***********************************************************
C
      INTEGER I, IERR, IFLG, NVAL, NX
      DOUBLE PRECISION HPPVAL, HPVAL, HVAL
C
      IFLG = IFLAG
      NVAL = NE
C
C Test for invalid input.
C
      IF (N .LT. 2  .OR.  IFLG .LT. 0  .OR.  IFLG .GT. 2
     .    .OR.  NVAL .LT. 1) GO TO 2
C
C Initialize the number of extrapolation points NX and
C   loop on evaluation points.
C
      NX = 0
      DO 1 I = 1,NVAL
        IF (IFLG .EQ. 0) THEN
          V(I) = HVAL (TE(I),N,X,Y,YP,SIGMA, IERR)
        ELSEIF (IFLG .EQ. 1) THEN
          V(I) = HPVAL (TE(I),N,X,Y,YP,SIGMA, IERR)
        ELSE
          V(I) = HPPVAL (TE(I),N,X,Y,YP,SIGMA, IERR)
        ENDIF
        IF (IERR .LT. 0) GO TO 3
        NX = NX + IERR
    1   CONTINUE
C
C No errors encountered.
C
      IER = NX
      RETURN
C
C N, IFLAG, or NE is outside its valid range.
C
    2 IER = -1
      RETURN
C
C X is not strictly increasing.
C
    3 IER = -2
      RETURN
      END
      SUBROUTINE TSVAL2 (N,T,X,Y,XP,YP,SIGMA,IFLAG,NE,
     .                   TE, VX,VY,IER)
      INTEGER N, IFLAG, NE, IER
      DOUBLE PRECISION T(N), X(N), Y(N), XP(N), YP(N),
     .                 SIGMA(N), TE(NE), VX(NE), VY(NE)
C
C***********************************************************
C
C                                                From TSPACK
C                                            Robert J. Renka
C                                  Dept. of Computer Science
C                                       Univ. of North Texas
C                                           renka@cs.unt.edu
C                                                   06/10/92
C
C   This subroutine returns values or derivatives of a pair
C of Hermite interpolatory tension splines H1 and H2 which
C form the components of a parametric planar curve C(t) =
C (H1(t),H2(t)).  Refer to Subroutines TSPBP and TSPSP.
C
C   Note that a large tension factor in SIGMA may cause
C underflow.  The result is assumed to be zero.  If not the
C default, this may be specified by either a compiler option
C or operating system option.
C
C On input:
C
C       N = Number of data points.  N .GE. 2.
C
C       T = Array of length N containing a strictly in-
C           creasing sequence of abscissae (parameter
C           values).  Refer to Subroutine ARCL2D.
C
C       X = Array of length N containing data values or
C           function values returned by Subroutine SMCRV.
C           X(I) = H1(T(I)) for I = 1,...,N.
C
C       Y = Array of length N containing data values or
C           function values returned by Subroutine SMCRV.
C           Y(I) = H2(T(I)) for I = 1,...,N.
C
C       XP = Array of length N containing first deriva-
C            tives.  XP(I) = H1P(T(I)) for I = 1,...,N,
C            where H1P denotes the derivative of H1.
C
C       YP = Array of length N containing first deriva-
C            tives.  YP(I) = H2P(T(I)) for I = 1,...,N,
C            where H2P denotes the derivative of H2.
C
C   Note that C(T(I)) = (X(I),Y(I)) and CP(T(I)) = (XP(I),
C YP(I)), I = 1,...,N, are data (control) points and deriva-
C tive (velocity) vectors, respectively.
C
C       SIGMA = Array of length N-1 containing tension fac-
C               tors whose absolute values determine the
C               balance between cubic and linear in each
C               interval.  SIGMA(I) is associated with int-
C               erval (I,I+1) for I = 1,...,N-1.
C
C       IFLAG = Output option indicator:
C               IFLAG = 0 if values of H1 and H2 (points on
C                         the curve) are to be computed.
C               IFLAG = 1 if first derivative vectors are to
C                         be computed.  Unit tangent vectors
C                         can be obtained by normalizing
C                         these to unit vectors.
C               IFLAG = 2 if second derivative (accelera-
C                         tion) vectors are to be computed.
C                         Given a unit tangent vector U and
C                         a second derivative vector V, the
C                         corresponding curvature vector
C                         can be computed as the cross
C                         product U X V X U.
C
C       NE = Number of evaluation points.  NE > 0.
C
C       TE = Array of length NE containing the evaluation
C            points.  The sequence should be strictly in-
C            creasing for maximum efficiency.  Extrapolation
C            is performed if a point is not in the interval
C            [T(1),T(N)].
C
C The above parameters are not altered by this routine.
C
C       VX,VY = Arrays of length at least NE.
C
C On output:
C
C       VX,VY = Arrays containing values, first derivatives,
C               or second derivatives of H1 and H2, respec-
C               tively, at the evaluation points (unless
C               IER < 0).  If IER = -1, VX and VY are not
C               altered.  If IER = -2, VX and VY may be only
C               partially defined.
C
C       IER = Error indicator:
C             IER = 0  if no errors were encountered and
C                      no extrapolation occurred.
C             IER > 0  if no errors were encountered but
C                      extrapolation was required at IER
C                      points.
C             IER = -1 if N < 2, IFLAG < 0, IFLAG > 2, or
C                      NE < 1.
C             IER = -2 if the abscissae are not in strictly
C                      increasing order.  (This error will
C                      not necessarily be detected.)
C
C Modules required by TSVAL2:  HPPVAL, HPVAL, HVAL, INTRVL,
C                                SNHCSH
C
C***********************************************************
C
      INTEGER I, IERR, IFLG, NVAL, NX
      DOUBLE PRECISION HPPVAL, HPVAL, HVAL
C
      IFLG = IFLAG
      NVAL = NE
C
C Test for invalid input.
C
      IF (N .LT. 2  .OR.  IFLG .LT. 0  .OR.  IFLG .GT. 2
     .    .OR.  NVAL .LT. 1) GO TO 2
C
C Initialize the number of extrapolation points NX and
C   loop on evaluation points.
C
      NX = 0
      DO 1 I = 1,NVAL
        IF (IFLG .EQ. 0) THEN
          VX(I) = HVAL (TE(I),N,T,X,XP,SIGMA, IERR)
          VY(I) = HVAL (TE(I),N,T,Y,YP,SIGMA, IERR)
        ELSEIF (IFLG .EQ. 1) THEN
          VX(I) = HPVAL (TE(I),N,T,X,XP,SIGMA, IERR)
          VY(I) = HPVAL (TE(I),N,T,Y,YP,SIGMA, IERR)
        ELSE
          VX(I) = HPPVAL (TE(I),N,T,X,XP,SIGMA, IERR)
          VY(I) = HPPVAL (TE(I),N,T,Y,YP,SIGMA, IERR)
        ENDIF
        IF (IERR .LT. 0) GO TO 3
        NX = NX + IERR
    1   CONTINUE
C
C No errors encountered.
C
      IER = NX
      RETURN
C
C N, IFLAG, or NE is outside its valid range.
C
    2 IER = -1
      RETURN
C
C T is not strictly increasing.
C
    3 IER = -2
      RETURN
      END
      SUBROUTINE TSVAL3 (N,T,X,Y,Z,XP,YP,ZP,SIGMA,IFLAG,NE,
     .                   TE, VX,VY,VZ,IER)
      INTEGER N, IFLAG, NE, IER
      DOUBLE PRECISION T(N), X(N), Y(N), Z(N), XP(N), YP(N),
     .                 ZP(N), SIGMA(N), TE(NE), VX(NE),
     .                 VY(NE), VZ(NE)
C
C***********************************************************
C
C                                                From TSPACK
C                                            Robert J. Renka
C                                  Dept. of Computer Science
C                                       Univ. of North Texas
C                                           renka@cs.unt.edu
C                                                   06/10/92
C
C   This subroutine returns values or derivatives of three
C Hermite interpolatory tension splines H1, H2, and H3 which
C form the components of a parametric space curve C(t) =
C (H1(t),H2(t),H3(t)).  Refer to Subroutines TSPBP and
C TSPSP.
C
C   Note that a large tension factor in SIGMA may cause
C underflow.  The result is assumed to be zero.  If not the
C default, this may be specified by either a compiler option
C or operating system option.
C
C On input:
C
C       N = Number of data points.  N .GE. 2.
C
C       T = Array of length N containing a strictly in-
C           creasing sequence of abscissae (parameter
C           values).  Refer to Subroutine ARCL3D.
C
C       X = Array of length N containing data values or
C           function values returned by Subroutine SMCRV.
C           X(I) = H1(T(I)) for I = 1,...,N.
C
C       Y = Array of length N containing data values or
C           function values returned by Subroutine SMCRV.
C           Y(I) = H2(T(I)) for I = 1,...,N.
C
C       Z = Array of length N containing data values or
C           function values returned by Subroutine SMCRV.
C           Z(I) = H3(T(I)) for I = 1,...,N.
C
C       XP = Array of length N containing first deriva-
C            tives.  XP(I) = H1P(T(I)) for I = 1,...,N,
C            where H1P denotes the derivative of H1.
C
C       YP = Array of length N containing first deriva-
C            tives.  YP(I) = H2P(T(I)) for I = 1,...,N,
C            where H2P denotes the derivative of H2.
C
C       ZP = Array of length N containing first deriva-
C            tives.  ZP(I) = H3P(T(I)) for I = 1,...,N,
C            where H3P denotes the derivative of H3.
C
C   Note that C(T(I)) = (X(I),Y(I),Z(I)) and CP(T(I)) =
C (XP(I),YP(I),ZP(I)), I = 1,...,N, are data (control)
C points and derivative (velocity) vectors, respectively.
C
C       SIGMA = Array of length N-1 containing tension fac-
C               tors whose absolute values determine the
C               balance between cubic and linear in each
C               interval.  SIGMA(I) is associated with int-
C               erval (I,I+1) for I = 1,...,N-1.
C
C       IFLAG = Output option indicator:
C               IFLAG = 0 if values of H1, H2, and H3
C                         (points on the curve) are to be
C                         computed.
C               IFLAG = 1 if first derivative vectors are to
C                         be computed.  Unit tangent vectors
C                         can be obtained by normalizing
C                         these to unit vectors.
C               IFLAG = 2 if second derivative (accelera-
C                         tion) vectors are to be computed.
C                         Given a unit tangent vector U and
C                         a second derivative vector V, the
C                         corresponding curvature vector
C                         can be computed as the cross
C                         product U X V X U.
C
C       NE = Number of evaluation points.  NE > 0.
C
C       TE = Array of length NE containing the evaluation
C            points.  The sequence should be strictly in-
C            creasing for maximum efficiency.  Extrapolation
C            is performed if a point is not in the interval
C            [T(1),T(N)].
C
C The above parameters are not altered by this routine.
C
C       VX,VY,VZ = Arrays of length at least NE.
C
C On output:
C
C       VX,VY,VZ = Arrays containing values, first deriva-
C                  tives, or second derivatives of H1, H2,
C                  and H3, respectively, at the evaluation
C                  points (unless IER < 0).  If IER = -1,
C                  VX, VY, and VZ are not altered.  If IER
C                  = -2, VX, VY, and VZ may be only partial-
C                  ly defined.
C
C       IER = Error indicator:
C             IER = 0  if no errors were encountered and
C                      no extrapolation occurred.
C             IER > 0  if no errors were encountered but
C                      extrapolation was required at IER
C                      points.
C             IER = -1 if N < 2, IFLAG < 0, IFLAG > 2, or
C                      NE < 1.
C             IER = -2 if the abscissae are not in strictly
C                      increasing order.  (This error will
C                      not necessarily be detected.)
C
C Modules required by TSVAL3:  HPPVAL, HPVAL, HVAL, INTRVL,
C                                SNHCSH
C
C***********************************************************
C
      INTEGER I, IERR, IFLG, NVAL, NX
      DOUBLE PRECISION HPPVAL, HPVAL, HVAL
C
      IFLG = IFLAG
      NVAL = NE
C
C Test for invalid input.
C
      IF (N .LT. 2  .OR.  IFLG .LT. 0  .OR.  IFLG .GT. 2
     .    .OR.  NVAL .LT. 1) GO TO 2
C
C Initialize the number of extrapolation points NX and
C   loop on evaluation points.
C
      NX = 0
      DO 1 I = 1,NVAL
        IF (IFLG .EQ. 0) THEN
          VX(I) = HVAL (TE(I),N,T,X,XP,SIGMA, IERR)
          VY(I) = HVAL (TE(I),N,T,Y,YP,SIGMA, IERR)
          VZ(I) = HVAL (TE(I),N,T,Z,ZP,SIGMA, IERR)
        ELSEIF (IFLG .EQ. 1) THEN
          VX(I) = HPVAL (TE(I),N,T,X,XP,SIGMA, IERR)
          VY(I) = HPVAL (TE(I),N,T,Y,YP,SIGMA, IERR)
          VZ(I) = HPVAL (TE(I),N,T,Z,ZP,SIGMA, IERR)
        ELSE
          VX(I) = HPPVAL (TE(I),N,T,X,XP,SIGMA, IERR)
          VY(I) = HPPVAL (TE(I),N,T,Y,YP,SIGMA, IERR)
          VZ(I) = HPPVAL (TE(I),N,T,Z,ZP,SIGMA, IERR)
        ENDIF
        IF (IERR .LT. 0) GO TO 3
        NX = NX + IERR
    1   CONTINUE
C
C No errors encountered.
C
      IER = NX
      RETURN
C
C N, IFLAG, or NE is outside its valid range.
C
    2 IER = -1
      RETURN
C
C T is not strictly increasing.
C
    3 IER = -2
      RETURN
      END
      SUBROUTINE YPC1 (N,X,Y, YP,IER)
      INTEGER N, IER
      DOUBLE PRECISION X(N), Y(N), YP(N)
C
C***********************************************************
C
C                                                From TSPACK
C                                            Robert J. Renka
C                                  Dept. of Computer Science
C                                       Univ. of North Texas
C                                           renka@cs.unt.edu
C                                                   06/10/92
C
C   This subroutine employs a three-point quadratic interpo-
C lation method to compute local derivative estimates YP
C associated with a set of data points.  The interpolation
C formula is the monotonicity-constrained parabolic method
C described in the reference cited below.  A Hermite int-
C erpolant of the data values and derivative estimates pre-
C serves monotonicity of the data.  Linear interpolation is
C used if N = 2.  The method is invariant under a linear
C scaling of the coordinates but is not additive.
C
C On input:
C
C       N = Number of data points.  N .GE. 2.
C
C       X = Array of length N containing a strictly in-
C           creasing sequence of abscissae:  X(I) < X(I+1)
C           for I = 1,...,N-1.
C
C       Y = Array of length N containing data values asso-
C           ciated with the abscissae.
C
C Input parameters are not altered by this routine.
C
C On output:
C
C       YP = Array of length N containing estimated deriv-
C            atives at the abscissae unless IER .NE. 0.
C            YP is not altered if IER = 1, and is only par-
C            tially defined if IER > 1.
C
C       IER = Error indicator:
C             IER = 0 if no errors were encountered.
C             IER = 1 if N < 2.
C             IER = I if X(I) .LE. X(I-1) for some I in the
C                     range 2,...,N.
C
C Reference:  J. M. Hyman, "Accurate Monotonicity-preserving
C               Cubic Interpolation",  LA-8796-MS, Los
C               Alamos National Lab, Feb. 1982.
C
C Modules required by YPC1:  None
C
C Intrinsic functions called by YPC1:  ABS, MAX, MIN, SIGN
C
C***********************************************************
C
      INTEGER I, NM1
      DOUBLE PRECISION ASI, ASIM1, DX2, DXI, DXIM1, S2, SGN,
     .                 SI, SIM1, T
C
      NM1 = N - 1
      IF (NM1 .LT. 1) GO TO 2
      I = 1
      DXI = X(2) - X(1)
      IF (DXI .LE. 0.D0) GO TO 3
      SI = (Y(2)-Y(1))/DXI
      IF (NM1 .EQ. 1) THEN
C
C Use linear interpolation for N = 2.
C
        YP(1) = SI
        YP(2) = SI
        IER = 0
        RETURN
      ENDIF
C
C N .GE. 3.  YP(1) = S1 + DX1*(S1-S2)/(DX1+DX2) unless this
C   results in YP(1)*S1 .LE. 0 or abs(YP(1)) > 3*abs(S1).
C
      I = 2
      DX2 = X(3) - X(2)
      IF (DX2 .LE. 0.D0) GO TO 3
      S2 = (Y(3)-Y(2))/DX2
      T = SI + DXI*(SI-S2)/(DXI+DX2)
      IF (SI .GE. 0.D0) THEN
        YP(1) = MIN(MAX(0.D0,T), 3.D0*SI)
      ELSE
        YP(1) = MAX(MIN(0.D0,T), 3.D0*SI)
      ENDIF
C
C YP(I) = (DXIM1*SI+DXI*SIM1)/(DXIM1+DXI) subject to the
C   constraint that YP(I) has the sign of either SIM1 or
C   SI, whichever has larger magnitude, and abs(YP(I)) .LE.
C   3*min(abs(SIM1),abs(SI)).
C
      DO 1 I = 2,NM1
        DXIM1 = DXI
        DXI = X(I+1) - X(I)
        IF (DXI .LE. 0.D0) GO TO 3
        SIM1 = SI
        SI = (Y(I+1)-Y(I))/DXI
        T = (DXIM1*SI+DXI*SIM1)/(DXIM1+DXI)
        ASIM1 = ABS(SIM1)
        ASI = ABS(SI)
        SGN = SIGN(1.D0,SI)
        IF (ASIM1 .GT. ASI) SGN = SIGN(1.D0,SIM1)
        IF (SGN .GT. 0.D0) THEN
          YP(I) = MIN(MAX(0.D0,T),3.D0*MIN(ASIM1,ASI))
        ELSE
          YP(I) = MAX(MIN(0.D0,T),-3.D0*MIN(ASIM1,ASI))
        ENDIF
    1   CONTINUE
C
C YP(N) = SNM1 + DXNM1*(SNM1-SNM2)/(DXNM2+DXNM1) subject to
C   the constraint that YP(N) has the sign of SNM1 and
C   abs(YP(N)) .LE. 3*abs(SNM1).  Note that DXI = DXNM1 and
C   SI = SNM1.
C
      T = SI + DXI*(SI-SIM1)/(DXIM1+DXI)
      IF (SI .GE. 0.D0) THEN
        YP(N) = MIN(MAX(0.D0,T), 3.D0*SI)
      ELSE
        YP(N) = MAX(MIN(0.D0,T), 3.D0*SI)
      ENDIF
      IER = 0
      RETURN
C
C N is outside its valid range.
C
    2 IER = 1
      RETURN
C
C X(I+1) .LE. X(I).
C
    3 IER = I + 1
      RETURN
      END
      SUBROUTINE YPC1P (N,X,Y, YP,IER)
      INTEGER N, IER
      DOUBLE PRECISION X(N), Y(N), YP(N)
C
C***********************************************************
C
C                                                From TSPACK
C                                            Robert J. Renka
C                                  Dept. of Computer Science
C                                       Univ. of North Texas
C                                           renka@cs.unt.edu
C                                                   06/10/92
C
C   This subroutine employs a three-point quadratic interpo-
C lation method to compute local derivative estimates YP
C associated with a set of N data points (X(I),Y(I)).  It
C is assumed that Y(N) = Y(1), and YP(N) = YP(1) on output.
C Thus, a Hermite interpolant H(x) defined by the data
C points and derivative estimates is periodic with period
C X(N)-X(1).  The derivative-estimation formula is the
C monotonicity-constrained parabolic fit described in the
C reference cited below:  H(x) is monotonic in intervals in
C which the data is monotonic.  The method is invariant
C under a linear scaling of the coordinates but is not
C additive.
C
C On input:
C
C       N = Number of data points.  N .GE. 3.
C
C       X = Array of length N containing a strictly in-
C           creasing sequence of abscissae:  X(I) < X(I+1)
C           for I = 1,...,N-1.
C
C       Y = Array of length N containing data values asso-
C           ciated with the abscissae.  Y(N) is set to Y(1)
C           on output unless IER = 1.
C
C   Input parameters, other than Y(N), are not altered by
C this routine.
C
C On output:
C
C       YP = Array of length N containing estimated deriv-
C            atives at the abscissae unless IER .NE. 0.
C            YP is not altered if IER = 1, and is only par-
C            tially defined if IER > 1.
C
C       IER = Error indicator:
C             IER = 0 if no errors were encountered.
C             IER = 1 if N < 3.
C             IER = I if X(I) .LE. X(I-1) for some I in the
C                     range 2,...,N.
C
C Reference:  J. M. Hyman, "Accurate Monotonicity-preserving
C               Cubic Interpolation",  LA-8796-MS, Los
C               Alamos National Lab, Feb. 1982.
C
C Modules required by YPC1P:  None
C
C Intrinsic functions called by YPC1P:  ABS, MAX, MIN, SIGN
C
C***********************************************************
C
      INTEGER I, NM1
      DOUBLE PRECISION ASI, ASIM1, DXI, DXIM1, SGN, SI,
     .                 SIM1, T
C
      NM1 = N - 1
      IF (NM1 .LT. 2) GO TO 2
      Y(N) = Y(1)
C
C Initialize for loop on interior points.
C
      I = 1
      DXI = X(2) - X(1)
      IF (DXI .LE. 0.D0) GO TO 3
      SI = (Y(2)-Y(1))/DXI
C
C YP(I) = (DXIM1*SI+DXI*SIM1)/(DXIM1+DXI) subject to the
C   constraint that YP(I) has the sign of either SIM1 or
C   SI, whichever has larger magnitude, and abs(YP(I)) .LE.
C   3*min(abs(SIM1),abs(SI)).
C
      DO 1 I = 2,NM1
        DXIM1 = DXI
        DXI = X(I+1) - X(I)
        IF (DXI .LE. 0.D0) GO TO 3
        SIM1 = SI
        SI = (Y(I+1)-Y(I))/DXI
        T = (DXIM1*SI+DXI*SIM1)/(DXIM1+DXI)
        ASIM1 = ABS(SIM1)
        ASI = ABS(SI)
        SGN = SIGN(1.D0,SI)
        IF (ASIM1 .GT. ASI) SGN = SIGN(1.D0,SIM1)
        IF (SGN .GT. 0.D0) THEN
          YP(I) = MIN(MAX(0.D0,T),3.D0*MIN(ASIM1,ASI))
        ELSE
          YP(I) = MAX(MIN(0.D0,T),-3.D0*MIN(ASIM1,ASI))
        ENDIF
    1   CONTINUE
C
C YP(N) = YP(1), I = 1, and IM1 = N-1.
C
      DXIM1 = DXI
      DXI = X(2) - X(1)
      SIM1 = SI
      SI = (Y(2) - Y(1))/DXI
      T = (DXIM1*SI + DXI*SIM1)/(DXIM1+DXI)
      ASIM1 = ABS(SIM1)
      ASI = ABS(SI)
      SGN = SIGN(1.D0,SI)
      IF (ASIM1 .GT. ASI) SGN = SIGN(1.D0,SIM1)
      IF (SGN .GT. 0.D0) THEN
        YP(1) = MIN(MAX(0.D0,T), 3.D0*MIN(ASIM1,ASI))
      ELSE
        YP(1) = MAX(MIN(0.D0,T), -3.D0*MIN(ASIM1,ASI))
      ENDIF
      YP(N) = YP(1)
C
C No error encountered.
C
      IER = 0
      RETURN
C
C N is outside its valid range.
C
    2 IER = 1
      RETURN
C
C X(I+1) .LE. X(I).
C
    3 IER = I + 1
      RETURN
      END
      SUBROUTINE YPC2 (N,X,Y,SIGMA,ISL1,ISLN,BV1,BVN,
     .                 WK, YP,IER)
      INTEGER N, ISL1, ISLN, IER
      DOUBLE PRECISION X(N), Y(N), SIGMA(N), BV1, BVN,
     .                 WK(N), YP(N)
C
C***********************************************************
C
C                                                From TSPACK
C                                            Robert J. Renka
C                                  Dept. of Computer Science
C                                       Univ. of North Texas
C                                           renka@cs.unt.edu
C                                                   06/10/92
C
C   This subroutine solves a linear system for a set of
C first derivatives YP associated with a Hermite interpola-
C tory tension spline H(x).  The derivatives are chosen so
C that H(x) has two continuous derivatives for all x and H
C satisfies user-specified end conditions.
C
C On input:
C
C       N = Number of data points.  N .GE. 2.
C
C       X = Array of length N containing a strictly in-
C           creasing sequence of abscissae:  X(I) < X(I+1)
C           for I = 1,...,N-1.
C
C       Y = Array of length N containing data values asso-
C           ciated with the abscissae.  H(X(I)) = Y(I) for
C           I = 1,...,N.
C
C       SIGMA = Array of length N-1 containing tension
C               factors.  SIGMA(I) is associated with inter-
C               val (X(I),X(I+1)) for I = 1,...,N-1.  If
C               SIGMA(I) = 0, H is the Hermite cubic interp-
C               olant of the data values and computed deriv-
C               atives at X(I) and X(I+1), and if all
C               tension factors are zero, H is the C-2 cubic
C               spline interpolant which satisfies the end
C               conditions.
C
C       ISL1 = Option indicator for the condition at X(1):
C              ISL1 = 0 if YP(1) is to be estimated inter-
C                       nally by a constrained parabolic
C                       fit to the first three points.
C                       This is identical to the method used
C                       by Subroutine YPC1.  BV1 is not used
C                       in this case.
C              ISL1 = 1 if the first derivative of H at X(1)
C                       is specified by BV1.
C              ISL1 = 2 if the second derivative of H at
C                       X(1) is specified by BV1.
C              ISL1 = 3 if YP(1) is to be estimated inter-
C                       nally from the derivative of the
C                       tension spline (using SIGMA(1))
C                       which interpolates the first three
C                       data points and has third derivative
C                       equal to zero at X(1).  Refer to
C                       ENDSLP.  BV1 is not used in this
C                       case.
C
C       ISLN = Option indicator for the condition at X(N):
C              ISLN = 0 if YP(N) is to be estimated inter-
C                       nally by a constrained parabolic
C                       fit to the last three data points.
C                       This is identical to the method used
C                       by Subroutine YPC1.  BVN is not used
C                       in this case.
C              ISLN = 1 if the first derivative of H at X(N)
C                       is specified by BVN.
C              ISLN = 2 if the second derivative of H at
C                       X(N) is specified by BVN.
C              ISLN = 3 if YP(N) is to be estimated inter-
C                       nally from the derivative of the
C                       tension spline (using SIGMA(N-1))
C                       which interpolates the last three
C                       data points and has third derivative
C                       equal to zero at X(N).  Refer to
C                       ENDSLP.  BVN is not used in this
C                       case.
C
C       BV1,BVN = Boundary values or dummy parameters as
C                 defined by ISL1 and ISLN.
C
C The above parameters are not altered by this routine.
C
C       WK = Array of length at least N-1 to be used as
C            temporary work space.
C
C       YP = Array of length .GE. N.
C
C On output:
C
C       YP = Array containing derivatives of H at the
C            abscissae.  YP is not defined if IER .NE. 0.
C
C       IER = Error indicator:
C             IER = 0 if no errors were encountered.
C             IER = 1 if N, ISL1, or ISLN is outside its
C                     valid range.
C             IER = I if X(I) .LE. X(I-1) for some I in the
C                     range 2,...,N.
C
C Modules required by YPC2:  ENDSLP, SNHCSH, YPCOEF
C
C Intrinsic function called by YPC2:  ABS
C
C***********************************************************
C
      INTEGER I, NM1, NN
      DOUBLE PRECISION D, D1, D2, DX, R1, R2, S, SD1, SD2,
     .                 SIG, YP1, YPN
      DOUBLE PRECISION ENDSLP
C
      NN = N
      IF (NN .LT. 2  .OR.  ISL1 .LT. 0  .OR.  ISL1 .GT. 3
     .    .OR.  ISLN .LT. 0  .OR.  ISLN .GT. 3) GO TO 3
      NM1 = NN - 1
C
C Set YP1 and YPN to the endpoint values.
C
      IF (ISL1 .EQ. 0) THEN
        IF (NN .GT. 2) YP1 = ENDSLP (X(1),X(2),X(3),Y(1),
     .                               Y(2),Y(3),0.D0)
      ELSEIF (ISL1 .NE. 3) THEN
        YP1 = BV1
      ELSE
        IF (NN .GT. 2) YP1 = ENDSLP (X(1),X(2),X(3),Y(1),
     .                               Y(2),Y(3),SIGMA(1))
      ENDIF
      IF (ISLN .EQ. 0) THEN
        IF (NN .GT. 2) YPN = ENDSLP (X(NN),X(NM1),X(NN-2),
     .                            Y(NN),Y(NM1),Y(NN-2),0.D0)
      ELSEIF (ISLN .NE. 3) THEN
        YPN = BVN
      ELSE
        IF (NN .GT. 2) YPN = ENDSLP (X(NN),X(NM1),X(NN-2),
     .                      Y(NN),Y(NM1),Y(NN-2),SIGMA(NM1))
      ENDIF
C
C Solve the symmetric positive-definite tridiagonal linear
C   system.  The forward elimination step consists of div-
C   iding each row by its diagonal entry, then introducing a
C   zero below the diagonal.  This requires saving only the
C   superdiagonal (in WK) and the right hand side (in YP).
C
      I = 1
      DX = X(2) - X(1)
      IF (DX .LE. 0.D0) GO TO 4
      S = (Y(2)-Y(1))/DX
      IF (NN .EQ. 2) THEN
        IF (ISL1 .EQ. 0  .OR.  ISL1 .EQ. 3) YP1 = S
        IF (ISLN .EQ. 0  .OR.  ISLN .EQ. 3) YPN = S
      ENDIF
C
C Begin forward elimination.
C
      SIG = ABS(SIGMA(1))
      CALL YPCOEF (SIG,DX, D1,SD1)
      R1 = (SD1+D1)*S
      WK(1) = 0.D0
      YP(1) = YP1
      IF (ISL1 .EQ. 2) THEN
        WK(1) = SD1/D1
        YP(1) = (R1-YP1)/D1
      ENDIF
      DO 1 I = 2,NM1
        DX = X(I+1) - X(I)
        IF (DX .LE. 0.D0) GO TO 4
        S = (Y(I+1)-Y(I))/DX
        SIG = ABS(SIGMA(I))
        CALL YPCOEF (SIG,DX, D2,SD2)
        R2 = (SD2+D2)*S
        D = D1 + D2 - SD1*WK(I-1)
        WK(I) = SD2/D
        YP(I) = (R1 + R2 - SD1*YP(I-1))/D
        D1 = D2
        SD1 = SD2
        R1 = R2
    1   CONTINUE
      D = D1 - SD1*WK(NM1)
      YP(NN) = YPN
      IF (ISLN .EQ. 2) YP(NN) = (R1 + YPN - SD1*YP(NM1))/D
C
C Back substitution:
C
      DO 2 I = NM1,1,-1
        YP(I) = YP(I) - WK(I)*YP(I+1)
    2   CONTINUE
      IER = 0
      RETURN
C
C Invalid integer input parameter.
C
    3 IER = 1
      RETURN
C
C Abscissae out of order or duplicate points.
C
    4 IER = I + 1
      RETURN
      END
      SUBROUTINE YPC2P (N,X,Y,SIGMA,WK, YP,IER)
      INTEGER N, IER
      DOUBLE PRECISION X(N), Y(N), SIGMA(N), WK(*), YP(N)
C
C***********************************************************
C
C                                                From TSPACK
C                                            Robert J. Renka
C                                  Dept. of Computer Science
C                                       Univ. of North Texas
C                                           renka@cs.unt.edu
C                                                   06/10/92
C
C   This subroutine solves a linear system for a set of
C first derivatives YP associated with a Hermite interpola-
C tory tension spline H(x).  The derivatives are chosen so
C that H(x) has two continuous derivatives for all x, and H
C satisfies periodic end conditions:  first and second der-
C ivatives of H at X(1) agree with those at X(N), and thus
C the length of a period is X(N) - X(1).  It is assumed that
C Y(N) = Y(1), and Y(N) is not referenced.
C
C On input:
C
C       N = Number of data points.  N .GE. 3.
C
C       X = Array of length N containing a strictly in-
C           creasing sequence of abscissae:  X(I) < X(I+1)
C           for I = 1,...,N-1.
C
C       Y = Array of length N containing data values asso-
C           ciated with the abscissae.  H(X(I)) = Y(I) for
C           I = 1,...,N.
C
C       SIGMA = Array of length N-1 containing tension
C               factors.  SIGMA(I) is associated with inter-
C               val (X(I),X(I+1)) for I = 1,...,N-1.  If
C               SIGMA(I) = 0, H is the Hermite cubic interp-
C               olant of the data values and computed deriv-
C               atives at X(I) and X(I+1), and if all
C               tension factors are zero, H is the C-2 cubic
C               spline interpolant which satisfies the end
C               conditions.
C
C The above parameters are not altered by this routine.
C
C       WK = Array of length at least 2N-2 to be used as
C            temporary work space.
C
C       YP = Array of length .GE. N.
C
C On output:
C
C       YP = Array containing derivatives of H at the
C            abscissae.  YP is not defined if IER .NE. 0.
C
C       IER = Error indicator:
C             IER = 0 if no errors were encountered.
C             IER = 1 if N is outside its valid range.
C             IER = I if X(I) .LE. X(I-1) for some I in the
C                     range 2,...,N.
C
C Modules required by YPC2P:  SNHCSH, YPCOEF
C
C Intrinsic function called by YPC2P:  ABS
C
C***********************************************************
C
      INTEGER I, NM1, NM2, NM3, NN, NP1, NPI
      DOUBLE PRECISION D, D1, D2, DIN, DNM1, DX, R1, R2,
     .                 RNM1, S, SD1, SD2, SDNM1, SIG, YPNM1
C
      NN = N
      IF (NN .LT. 3) GO TO 4
      NM1 = NN - 1
      NM2 = NN - 2
      NM3 = NN - 3
      NP1 = NN + 1
C
C The system is order N-1, symmetric, positive-definite, and
C   tridiagonal except for nonzero elements in the upper
C   right and lower left corners.  The forward elimination
C   step zeros the subdiagonal and divides each row by its
C   diagonal entry for the first N-2 rows.  The superdiago-
C   nal is stored in WK(I), the negative of the last column
C   (fill-in) in WK(N+I), and the right hand side in YP(I)
C   for I = 1,...,N-2.
C
      I = NM1
      DX = X(NN) - X(NM1)
      IF (DX .LE. 0.D0) GO TO 5
      S = (Y(1)-Y(NM1))/DX
      SIG = ABS(SIGMA(NM1))
      CALL YPCOEF (SIG,DX, DNM1,SDNM1)
      RNM1 = (SDNM1+DNM1)*S
      I = 1
      DX = X(2) - X(1)
      IF (DX .LE. 0.D0) GO TO 5
      S = (Y(2)-Y(1))/DX
      SIG = ABS(SIGMA(1))
      CALL YPCOEF (SIG,DX, D1,SD1)
      R1 = (SD1+D1)*S
      D = DNM1 + D1
      WK(1) = SD1/D
      WK(NP1) = -SDNM1/D
      YP(1) = (RNM1+R1)/D
      DO 1 I = 2,NM2
        DX = X(I+1) - X(I)
        IF (DX .LE. 0.D0) GO TO 5
        S = (Y(I+1)-Y(I))/DX
        SIG = ABS(SIGMA(I))
        CALL YPCOEF (SIG,DX, D2,SD2)
        R2 = (SD2+D2)*S
        D = D1 + D2 - SD1*WK(I-1)
        DIN = 1.D0/D
        WK(I) = SD2*DIN
        NPI = NN + I
        WK(NPI) = -SD1*WK(NPI-1)*DIN
        YP(I) = (R1 + R2 - SD1*YP(I-1))*DIN
        SD1 = SD2
        D1 = D2
        R1 = R2
    1   CONTINUE
C
C The backward elimination step zeros the superdiagonal
C   (first N-3 elements).  WK(I) and YP(I) are overwritten
C   with the negative of the last column and the new right
C   hand side, respectively, for I = N-2, N-3, ..., 1.
C
      NPI = NN + NM2
      WK(NM2) = WK(NPI) - WK(NM2)
      DO 2 I = NM3,1,-1
        YP(I) = YP(I) - WK(I)*YP(I+1)
        NPI = NN + I
        WK(I) = WK(NPI) - WK(I)*WK(I+1)
    2   CONTINUE
C
C Solve the last equation for YP(N-1).
C
      YPNM1 = (R1 + RNM1 - SDNM1*YP(1) - SD1*YP(NM2))/
     .        (D1 + DNM1 + SDNM1*WK(1) + SD1*WK(NM2))
C
C Back substitute for the remainder of the solution
C   components.
C
      YP(NM1) = YPNM1
      DO 3 I = 1,NM2
        YP(I) = YP(I) + WK(I)*YPNM1
    3   CONTINUE
C
C YP(N) = YP(1).
C
      YP(N) = YP(1)
      IER = 0
      RETURN
C
C N is outside its valid range.
C
    4 IER = 1
      RETURN
C
C Abscissae out of order or duplicate points.
C
    5 IER = I + 1
      RETURN
      END
      SUBROUTINE YPCOEF (SIGMA,DX, D,SD)
      DOUBLE PRECISION SIGMA, DX, D, SD
C
C***********************************************************
C
C                                                From TSPACK
C                                            Robert J. Renka
C                                  Dept. of Computer Science
C                                       Univ. of North Texas
C                                           renka@cs.unt.edu
C                                                   11/17/96
C
C   This subroutine computes the coefficients of the deriva-
C tives in the symmetric diagonally dominant tridiagonal
C system associated with the C-2 derivative estimation pro-
C cedure for a Hermite interpolatory tension spline.
C
C On input:
C
C       SIGMA = Nonnegative tension factor associated with
C               an interval.
C
C       DX = Positive interval width.
C
C Input parameters are not altered by this routine.
C
C On output:
C
C       D = Component of the diagonal term associated with
C           the interval.  D = SIG*(SIG*COSHM(SIG) -
C           SINHM(SIG))/(DX*E), where SIG = SIGMA and E =
C           SIG*SINH(SIG) - 2*COSHM(SIG).
C
C       SD = Subdiagonal (superdiagonal) term.  SD = SIG*
C            SINHM(SIG)/E.
C
C Module required by YPCOEF:  SNHCSH
C
C Intrinsic function called by YPCOEF:  EXP
C
C***********************************************************
C
      DOUBLE PRECISION COSHM, COSHMM, E, EMS, SCM, SIG,
     .                 SINHM, SSINH, SSM
C
      SIG = SIGMA
      IF (SIG .LT. 1.D-9) THEN
C
C SIG = 0:  cubic interpolant.
C
        D = 4.D0/DX
        SD = 2.D0/DX
      ELSEIF (SIG .LE. .5D0) THEN
C
C 0 .LT. SIG .LE. .5:  use approximations designed to avoid
C                      cancellation error in the hyperbolic
C                      functions when SIGMA is small.
C
        CALL SNHCSH (SIG, SINHM,COSHM,COSHMM)
        E = (SIG*SINHM - COSHMM - COSHMM)*DX
        D = SIG*(SIG*COSHM-SINHM)/E
        SD = SIG*SINHM/E
      ELSE
C
C SIG > .5:  scale SINHM and COSHM by 2*EXP(-SIG) in order
C            to avoid overflow when SIGMA is large.
C
        EMS = EXP(-SIG)
        SSINH = 1.D0 - EMS*EMS
        SSM = SSINH - 2.D0*SIG*EMS
        SCM = (1.D0-EMS)*(1.D0-EMS)
        E = (SIG*SSINH - SCM - SCM)*DX
        D = SIG*(SIG*SCM-SSM)/E
        SD = SIG*SSM/E
      ENDIF
      RETURN
      END

