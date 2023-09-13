!
!=======================================================================
!
      SUBROUTINE SAXPY(N,SA,SX,INCX,SY,INCY)

!
!     CONSTANT TIMES A VECTOR PLUS A VECTOR.
!     USES UNROLLED LOOP FOR INCREMENTS EQUAL TO ONE.
!     JACK DONGARRA, LINPACK, 3/11/78.
!
      REAL SX(1),SY(1),SA
      INTEGER I,INCX,INCY,IX,IY,M,MP1,N
!
      IF(N.LE.0)RETURN
      IF (SA .EQ. 0.0) RETURN
      IF(INCX.EQ.1.AND.INCY.EQ.1)GO TO 20
!
!        CODE FOR UNEQUAL INCREMENTS OR EQUAL INCREMENTS
!          NOT EQUAL TO 1
!
      IX = 1
      IY = 1
      IF(INCX.LT.0)IX = (-N+1)*INCX + 1
      IF(INCY.LT.0)IY = (-N+1)*INCY + 1
      DO 10 I = 1,N
        SY(IY) = SY(IY) + SA*SX(IX)
        IX = IX + INCX
        IY = IY + INCY
   10 CONTINUE
      RETURN
!
!        CODE FOR BOTH INCREMENTS EQUAL TO 1
!
!
!        CLEAN-UP LOOP
!
   20 M = MOD(N,4)
      IF( M .EQ. 0 ) GO TO 40
      DO 30 I = 1,M
        SY(I) = SY(I) + SA*SX(I)
   30 CONTINUE
      IF( N .LT. 4 ) RETURN
   40 MP1 = M + 1
      DO 50 I = MP1,N,4
        SY(I) = SY(I) + SA*SX(I)
        SY(I + 1) = SY(I + 1) + SA*SX(I + 1)
        SY(I + 2) = SY(I + 2) + SA*SX(I + 2)
        SY(I + 3) = SY(I + 3) + SA*SX(I + 3)
   50 CONTINUE
      RETURN
      END SUBROUTINE
!
!=======================================================================
!
      SUBROUTINE CFTFAX(N,IFAX,TRIGS)

      DIMENSION IFAX(13),TRIGS(1)
!
!     THIS ROUTINE WAS MODIFIED FROM TEMPERTON'S ORIGINAL
!     BY DAVE FULKER.  IT NO LONGER PRODUCES FACTORS IN ASCENDING
!     ORDER, AND THERE ARE NONE OF THE ORIGINAL 'MODE' OPTIONS.
!
! ON INPUT     N
!               THE LENGTH OF EACH COMPLEX TRANSFORM TO BE PERFORMED
!
!               N MUST BE GREATER THAN 1 AND CONTAIN NO PRIME
!               FACTORS GREATER THAN 5.
!
! ON OUTPUT    IFAX
!               IFAX(1)
!                 THE NUMBER OF FACTORS CHOSEN OR -99 IN CASE OF ERROR
!               IFAX(2) THRU IFAX( IFAX(1)+1 )
!                 THE FACTORS OF N IN THE FOLLOWIN ORDER:  APPEARING
!                 FIRST ARE AS MANY FACTORS OF 4 AS CAN BE OBTAINED.
!                 SUBSEQUENT FACTORS ARE PRIMES, AND APPEAR IN
!                 ASCENDING ORDER, EXCEPT FOR MULTIPLE FACTORS.
!
!              TRIGS
!               2N SIN AND COS VALUES FOR USE BY THE TRANSFORM ROUTINE
!
      CALL FACT(N,IFAX)
      K = IFAX(1)
      IF (K .LT. 1 .OR. IFAX(K+1) .GT. 5) IFAX(1) = -99
      IF (IFAX(1) .LE. 0 )THEN
        WRITE(*,1900)N
1900  FORMAT(' FFTFAX - Invalid N=',I20)
        RETURN
        ENDIF
      CALL CFTRIG (N, TRIGS)
      RETURN
      END
!
!=======================================================================
!
      SUBROUTINE CFTRIG(N,TRIGS)

      DIMENSION TRIGS(1)
      PI=2.0*ASIN(1.0)
      DEL=(PI+PI)/FLOAT(N)
      L=N+N
      DO 10 I=1,L,2
      ANGLE=0.5*FLOAT(I-1)*DEL
      TRIGS(I)=COS(ANGLE)
      TRIGS(I+1)=SIN(ANGLE)
   10 CONTINUE
      RETURN
      END
!
!=======================================================================
!
      SUBROUTINE RFFTMLT(A,WORK,TRIGS,IFAX,INC,JUMP,N,LOT,ISGN)

!
!     "RFFTMLT" - MULTIPLE REAL/HALF-COMPLEX PERIODIC
!     FAST FOURIER TRANSFORM
!
!     PROCEDURE USED TO CONVERT TO HALF-LENGTH COMPLEX TRANSFORM
!     IS GIVEN BY COOLEY, LEWIS AND WELCH (J. SOUND VIB., VOL. 12
!     (1970), 315-337)
!
!     A IS THE ARRAY CONTAINING INPUT AND OUTPUT DATA
!     WORK IS AN AREA OF SIZE (N+1)*LOT
!     TRIGS IS A PREVIOUSLY PREPARED LIST OF TRIG FUNCTION VALUES
!     IFAX IS A PREVIOUSLY PREPARED LIST OF FACTORS OF N/2
!     INC IS THE INCREMENT WITHIN EACH DATA 'VECTOR'
!         (E.G. INC=1 FOR CONSECUTIVELY STORED DATA)
!     JUMP IS THE INCREMENT BETWEEN THE START OF EACH DATA VECTOR
!     N IS THE LENGTH OF THE DATA VECTORS
!     LOT IS THE NUMBER OF DATA VECTORS
!     ISGN = +1 FOR TRANSFORM FROM SPECTRAL TO GRIDPOINT
!           = -1 FOR TRANSFORM FROM GRIDPOINT TO SPECTRAL
!
!     ORDERING OF COEFFICIENTS:
!         A(0),B(0),A(1),B(1),A(2),B(2),...,A(N/2),B(N/2)
!         WHERE B(0)=B(N/2)=0; (N+2) LOCATIONS REQUIRED
!
!     ORDERING OF DATA:
!         X(0),X(1),X(2),...,X(N-1)
!
!     VECTORIZATION IS ACHIEVED ON CRAY BY DOING THE TRANSFORMS IN
!     PARALLEL
!
!     *** N.B. N IS ASSUMED TO BE AN EVEN NUMBER
!
!     DEFINITION OF TRANSFORMS:
!     -------------------------
!
!     ISGN=+1: X(J)=SUM(K=0,...,N-1)(C(K)*EXP(2*I*J*K*PI/N))
!         WHERE C(K)=A(K)+I*B(K) AND C(N-K)=A(K)-I*B(K)
!
!     ISGN=-1: A(K)=(1/N)*SUM(J=0,...,N-1)(X(J)*COS(2*J*K*PI/N))
!               B(K)=-(1/N)*SUM(J=0,...,N-1)(X(J)*SIN(2*J*K*PI/N))
!
      DIMENSION A(N),WORK(N),TRIGS(N),IFAX(1)
      NFAX=IFAX(1)
      NX=N+1
      NH=N/2
      INK=INC+INC
      IF (ISGN.EQ.+1) GO TO 30
!
!     IF NECESSARY, TRANSFER DATA TO WORK AREA
      IGO=50
      IF (MOD(NFAX,2).EQ.1) GOTO 40
      IBASE=1
      JBASE=1
      DO 20 L=1,LOT
      I=IBASE
      J=JBASE
! 
      DO 10 M=1,N
      WORK(J)=A(I)
      I=I+INC
      J=J+1
   10 CONTINUE
      IBASE=IBASE+JUMP
      JBASE=JBASE+NX
   20 CONTINUE
!
      IGO=60
      GO TO 40
!
!     PREPROCESSING (ISGN=+1)
!     ------------------------
!
   30 CONTINUE
      CALL FFT99A(A,WORK,TRIGS,INC,JUMP,N,LOT)
      IGO=60
!
!     COMPLEX TRANSFORM
!     -----------------
!
   40 CONTINUE
      IA=1
      LA=1
      DO 80 K=1,NFAX
      IF (IGO.EQ.60) GO TO 60
   50 CONTINUE
      CALL VPASSM(A(IA),A(IA+INC),WORK(1),WORK(2),TRIGS,INK,2,JUMP,NX,LOT,NH,IFAX(K+1),LA)
      IGO=60
      GO TO 70
   60 CONTINUE
      CALL VPASSM(WORK(1),WORK(2),A(IA),A(IA+INC),TRIGS,2,INK,NX,JUMP,LOT,NH,IFAX(K+1),LA)
      IGO=50
   70 CONTINUE
      LA=LA*IFAX(K+1)
   80 CONTINUE
!
      IF (ISGN.EQ.-1) GO TO 130
!
!     IF NECESSARY, TRANSFER DATA FROM WORK AREA
      IF (MOD(NFAX,2).EQ.1) GO TO 110
      IBASE=1
      JBASE=1
      DO 100 L=1,LOT
      I=IBASE
      J=JBASE
! 
      DO 90 M=1,N
      A(J)=WORK(I)
      I=I+1
      J=J+INC
   90 CONTINUE
      IBASE=IBASE+NX
      JBASE=JBASE+JUMP
  100 CONTINUE
!
!     FILL IN ZEROS AT END
  110 CONTINUE
      IB=N*INC+1
! 
      DO 120 L=1,LOT
      A(IB)=0.0
      A(IB+INC)=0.0
      IB=IB+JUMP
  120 CONTINUE
      GO TO 140
!
!     POSTPROCESSING (ISGN=-1):
!     --------------------------
!
  130 CONTINUE
      CALL FFT99B(WORK,A,TRIGS,INC,JUMP,N,LOT)
!
  140 CONTINUE
      RETURN
      END
!
!=======================================================================
!
      SUBROUTINE FFTFAX(N,IFAX,TRIGS)

      DIMENSION IFAX(13),TRIGS(1)
!
! MODE 3 IS USED FOR REAL/HALF-COMPLEX TRANSFORMS.  IT IS POSSIBLE
! TO DO COMPLEX/COMPLEX TRANSFORMS WITH OTHER VALUES OF MODE, BUT
! DOCUMENTATION OF THE DETAILS WERE NOT AVAILABLE WHEN THIS ROUTINE
! WAS WRITTEN.
!
      DATA MODE /3/
      CALL FAX (IFAX, N, MODE)
      I = IFAX(1)
      IF (IFAX(I+1) .GT. 5 .OR. N .LE. 4) IFAX(1) = -99
      IF (IFAX(1) .LE. 0 )THEN
        WRITE(*,1900)N
1900    FORMAT(' FFTFAX - Invalid N=',I20)
        RETURN
        ENDIF
      CALL FFTRIG (TRIGS, N, MODE)
      RETURN
      END
!
!=======================================================================
!
      SUBROUTINE FFTRIG(TRIGS,N,MODE)

      DIMENSION TRIGS(1)
      PI=2.0*ASIN(1.0)
      IMODE=IABS(MODE)
      NN=N
      IF (IMODE.GT.1.AND.IMODE.LT.6) NN=N/2
      DEL=(PI+PI)/FLOAT(NN)
      L=NN+NN
      DO 10 I=1,L,2
      ANGLE=0.5*FLOAT(I-1)*DEL
      TRIGS(I)=COS(ANGLE)
      TRIGS(I+1)=SIN(ANGLE)
   10 CONTINUE
      IF (IMODE.EQ.1) RETURN
      IF (IMODE.EQ.8) RETURN
      DEL=0.5*DEL
      NH=(NN+1)/2
      L=NH+NH
      LA=NN+NN
      DO 20 I=1,L,2
      ANGLE=0.5*FLOAT(I-1)*DEL
      TRIGS(LA+I)=COS(ANGLE)
      TRIGS(LA+I+1)=SIN(ANGLE)
   20 CONTINUE
      IF (IMODE.LE.3) RETURN
      DEL=0.5*DEL
      LA=LA+NN
      IF (MODE.EQ.5) GO TO 40
      DO 30 I=2,NN
      ANGLE=FLOAT(I-1)*DEL
      TRIGS(LA+I)=2.0*SIN(ANGLE)
   30 CONTINUE
      RETURN
   40 CONTINUE
      DEL=0.5*DEL
      DO 50 I=2,N
      ANGLE=FLOAT(I-1)*DEL
      TRIGS(LA+I)=SIN(ANGLE)
   50 CONTINUE
      RETURN
      END
!
!=======================================================================
!
      SUBROUTINE CFFTMLT(AR,AI,WORK,TRIGS,IFAX,INC,JUMP,N,LOT,ISGN)

!
! PURPOSE      PERFORMS MULTIPLE FAST FOURIER TRANSFORMS.  THIS PACKAGE
!              WILL PERFORM A NUMBER OF SIMULTANEOUS COMPLEX PERIODIC
!              FOURIER TRANSFORMS OR CORRESPONDING INVERSE TRANSFORMS.
!              THAT IS, GIVEN A SET OF COMPLEX GRIDPOINT VECTORS, THE
!              PACKAGE RETURNS A SET OF COMPLEX FOURIER
!              COEFFICIENT VECTORS, OR VICE VERSA.  THE LENGTH OF THE
!              TRANSFORMS MUST BE A NUMBER GREATER THAN 1 THAT HAS
!              NO PRIME FACTORS OTHER THAN 2, 3, AND 5.
!
!              THE PACKAGE CFFT99 CONTAINS SEVERAL USER-LEVEL ROUTINES:
!
!            SUBROUTINE CFTFAX
!                AN INITIALIZATION ROUTINE THAT MUST BE CALLED ONCE
!                BEFORE A SEQUENCE OF CALLS TO CFFT99
!                (PROVIDED THAT N IS NOT CHANGED).
!
!            SUBROUTINE CFFT99
!                THE ACTUAL TRANSFORM ROUTINE ROUTINE, CABABLE OF
!                PERFORMING BOTH THE TRANSFORM AND ITS INVERSE.
!                HOWEVER, AS THE TRANSFORMS ARE NOT NORMALIZED,
!                THE APPLICATION OF A TRANSFORM FOLLOWED BY ITS
!                INVERSE WILL YIELD THE ORIGINAL VALUES MULTIPLIED
!                BY N.
!
!
! ACCESS       *FORTRAN,P=XLIB,SN=CFFT99
!
!
! USAGE        LET N BE OF THE FORM 2**P * 3**Q * 5**R, WHERE P .GE. 0,
!              Q .GE. 0, AND R .GE. 0.  THEN A TYPICAL SEQUENCE OF
!              CALLS TO TRANSFORM A GIVEN SET OF COMPLEX VECTORS OF
!              LENGTH N TO A SET OF (UNSCALED) COMPLEX FOURIER
!              COEFFICIENT VECTORS OF LENGTH N IS
!
!                   DIMENSION IFAX(13),TRIGS(2*N)
!                   COMPLEX A(...), WORK(...)
!
!                   CALL CFTFAX (N, IFAX, TRIGS)
!                   CALL CFFT99 (A,WORK,TRIGS,IFAX,INC,JUMP,N,LOT,ISGN)
!
!              THE OUTPUT VECTORS OVERWRITE THE INPUT VECTORS, AND
!              THESE ARE STORED IN A.  WITH APPROPRIATE CHOICES FOR
!              THE OTHER ARGUMENTS, THESE VECTORS MAY BE CONSIDERED
!              EITHER THE ROWS OR THE COLUMNS OF THE ARRAY A.
!              SEE THE INDIVIDUAL WRITE-UPS FOR CFTFAX AND
!              CFFT99 BELOW, FOR A DETAILED DESCRIPTION OF THE
!              ARGUMENTS.
!
! HISTORY      THE PACKAGE WAS WRITTEN BY CLIVE TEMPERTON AT ECMWF IN
!              NOVEMBER, 1978.  IT WAS MODIFIED, DOCUMENTED, AND TESTED
!              FOR NCAR BY RUSS REW IN SEPTEMBER, 1980.  IT WAS
!              FURTHER MODIFIED FOR THE FULLY COMPLEX CASE BY DAVE
!              FULKER IN NOVEMBER, 1980.
!
!-----------------------------------------------------------------------
!
! SUBROUTINE CFTFAX (N,IFAX,TRIGS)
!
! PURPOSE      A SET-UP ROUTINE FOR CFFT99.  IT NEED ONLY BE
!              CALLED ONCE BEFORE A SEQUENCE OF CALLS TO CFFT99,
!              PROVIDED THAT N IS NOT CHANGED.
!
! ARGUMENT     IFAX(13),TRIGS(2*N)
! DIMENSIONS
!
! ARGUMENTS
!
! ON INPUT     N
!               AN EVEN NUMBER GREATER THAN 1 THAT HAS NO PRIME FACTOR
!               GREATER THAN 5.  N IS THE LENGTH OF THE TRANSFORMS (SEE
!               THE DOCUMENTATION FOR CFFT99 FOR THE DEFINITION OF
!               THE TRANSFORMS).
!
!              IFAX
!               AN INTEGER ARRAY.  THE NUMBER OF ELEMENTS ACTUALLY USED
!               WILL DEPEND ON THE FACTORIZATION OF N.  DIMENSIONING
!               IFAX FOR 13 SUFFICES FOR ALL N LESS THAN 1 MILLION.
!
!              TRIGS
!               A REAL ARRAY OF DIMENSION 2*N
!
! ON OUTPUT    IFAX
!               CONTAINS THE FACTORIZATION OF N.  IFAX(1) IS THE
!               NUMBER OF FACTORS, AND THE FACTORS THEMSELVES ARE STORED
!               IN IFAX(2),IFAX(3),...  IF N HAS ANY PRIME FACTORS
!               GREATER THAN 5, IFAX(1) IS SET TO -99.
!
!              TRIGS
!               AN ARRAY OF TRIGONOMETRI! FUNCTION VALUES SUBSEQUENTLY
!               USED BY THE CFT ROUTINES.
!
!-----------------------------------------------------------------------
!
! SUBROUTINE CFFT99 (A,WORK,TRIGS,IFAX,INC,JUMP,N,LOT,ISGN)
!
! PURPOSE      PERFORM A NUMBER OF SIMULTANEOUS (UNNORMALIZED) COMPLEX
!              PERIODIC FOURIER TRANSFORMS OR CORRESPONDING INVERSE
!              TRANSFORMS.  GIVEN A SET OF COMPLEX GRIDPOINT
!              VECTORS, THE PACKAGE RETURNS A SET OF
!              COMPLEX FOURIER COEFFICIENT VECTORS, OR VICE
!              VERSA.  THE LENGTH OF THE TRANSFORMS MUST BE A
!              NUMBER HAVING NO PRIME FACTORS OTHER THAN
!              2, 3, AND 5.  THIS ROUTINE IS
!              OPTIMIZED FOR USE ON THE CRAY-1.
!
! ARGUMENT     COMPLEX A(N*INC+(LOT-1)*JUMP), WORK(N*LOT)
! DIMENSIONS   REAL TRIGS(2*N), INTEGER IFAX(13)
!
! ARGUMENTS
!
! ON INPUT     A
!               A COMPLEX ARRAY OF LENGTH N*INC+(LOT-1)*JUMP CONTAINING
!               THE INPUT GRIDPOINT OR COEFFICIENT VECTORS.  THIS ARRAY
!               OVERWRITTEN BY THE RESULTS.
!
!               N.B. ALTHOUGH THE ARRAY A IS USUALLY CONSIDERED TO BE OF
!               TYPE COMPLEX IN THE CALLING PROGRAM, IT IS TREATED AS
!               REAL WITHIN THE TRANSFORM PACKAGE.  THIS REQUIRES THAT
!               SUCH TYPE CONFLICTS ARE PERMITTED IN THE USER"S
!               ENVIRONMENT, AND THAT THE STORAGE OF COMPLEX NUMBERS
!               MATCHES THE ASSUMPTIONS OF THIS ROUTINE.  THIS ROUTINE
!               ASSUMES THAT THE REAL AND IMAGINARY PORTIONS OF A
!               COMPLEX NUMBER OCCUPY ADJACENT ELEMENTS OF MEMORY.  IF
!               THESE CONDITIONS ARE NOT MET, THE USER MUST TREAT THE
!               ARRAY A AS REAL (AND OF TWICE THE ABOVE LENGTH), AND
!               WRITE THE CALLING PROGRAM TO TREAT THE REAL AND
!               IMAGINARY PORTIONS EXPLICITLY.
!
!              WORK
!               A COMPLEX WORK ARRAY OF LENGTH N*LOT OR A REAL ARRAY
!               OF LENGTH 2*N*LOT.  SEE N.B. ABOVE.
!
!              TRIGS
!               AN ARRAY SET UP BY CFTFAX, WHICH MUST BE CALLED FIRST.
!
!              IFAX
!               AN ARRAY SET UP BY CFTFAX, WHICH MUST BE CALLED FIRST.
!
!
!               N.B. IN THE FOLLOWING ARGUMENTS, INCREMENTS ARE MEASURED
!               IN WORD PAIRS, BECAUSE EACH COMPLEX ELEMENT IS ASSUMED
!               TO OCCUPY AN ADJACENT PAIR OF WORDS IN MEMORY.
!
!              INC
!               THE INCREMENT (IN WORD PAIRS) BETWEEN SUCCESSIVE ELEMENT
!               OF EACH (COMPLEX) GRIDPOINT OR COEFFICIENT VECTOR
!               (E.G.  INC=1 FOR CONSECUTIVELY STORED DATA).
!
!              JUMP
!               THE INCREMENT (IN WORD PAIRS) BETWEEN THE FIRST ELEMENTS
!               OF SUCCESSIVE DATA OR COEFFICIENT VECTORS.  ON THE CRAY-
!               TRY TO ARRANGE DATA SO THAT JUMP IS NOT A MULTIPLE OF 8
!               (TO AVOID MEMORY BANK CONFLICTS).  FOR CLARIFICATION OF
!               INC AND JUMP, SEE THE EXAMPLES BELOW.
!
!              N
!               THE LENGTH OF EACH TRANSFORM (SEE DEFINITION OF
!               TRANSFORMS, BELOW).
!
!              LOT
!               THE NUMBER OF TRANSFORMS TO BE DONE SIMULTANEOUSLY.
!
!              ISGN
!               = -1 FOR A TRANSFORM FROM GRIDPOINT VALUES TO FOURIER
!                    COEFFICIENTS.
!               = +1 FOR A TRANSFORM FROM FOURIER COEFFICIENTS TO
!                    GRIDPOINT VALUES.
!
! ON OUTPUT    A
!               IF ISGN = -1, AND LOT GRIDPOINT VECTORS ARE SUPPLIED,
!               EACH CONTAINING THE COMPLEX SEQUENCE:
!
!               G(0),G(1), ... ,G(N-1)  (N COMPLEX VALUES)
!
!               THEN THE RESULT CONSISTS OF LOT COMPLEX VECTORS EACH
!               CONTAINING THE CORRESPONDING N COEFFICIENT VALUES:
!
!               C(0),C(1), ... ,C(N-1)  (N COMPLEX VALUES)
!
!               DEFINED BY:
!                 C(K) = SUM(J=0,...,N-1)( G(J)*EXP(-2*I*J*K*PI/N) )
!                 WHERE I = SQRT(-1)
!
!
!               IF ISGN = +1, AND LOT COEFFICIENT VECTORS ARE SUPPLIED,
!               EACH CONTAINING THE COMPLEX SEQUENCE:
!
!               C(0),C(1), ... ,C(N-1)  (N COMPLEX VALUES)
!
!               THEN THE RESULT CONSISTS OF LOT COMPLEX VECTORS EACH
!               CONTAINING THE CORRESPONDING N GRIDPOINT VALUES:
!
!               G(0),G(1), ... ,G(N-1)  (N COMPLEX VALUES)
!
!               DEFINED BY:
!                 G(J) = SUM(K=0,...,N-1)( G(K)*EXP(+2*I*J*K*PI/N) )
!                 WHERE I = SQRT(-1)
!
!
!               A CALL WITH ISGN=-1 FOLLOWED BY A CALL WITH ISGN=+1
!               (OR VICE VERSA) RETURNS THE ORIGINAL DATA, MULTIPLIED
!               BY THE FACTOR N.
!
!
! EXAMPLE       GIVEN A 64 BY 9 GRID OF COMPLEX VALUES, STORED IN
!               A 66 BY 9 COMPLEX ARRAY, A, COMPUTE THE TWO DIMENSIONAL
!               FOURIER TRANSFORM OF THE GRID.  FROM TRANSFORM THEORY,
!               IT IS KNOWN THAT A TWO DIMENSIONAL TRANSFORM CAN BE
!               OBTAINED BY FIRST TRANSFORMING THE GRID ALONG ONE
!               DIRECTION, THEN TRANSFORMING THESE RESULTS ALONG THE
!               ORTHOGONAL DIRECTION.
!
!               COMPLEX A(66,9), WORK(64,9)
!               REAL TRIGS1(128), TRIGS2(18)
!               INTEGER IFAX1(13), IFAX2(13)
!
!               SET UP THE IFAX AND TRIGS ARRAYS FOR EACH DIRECTION:
!
!               CALL CFTFAX(64, IFAX1, TRIGS1)
!               CALL CFTFAX( 9, IFAX2, TRIGS2)
!
!               IN THIS CASE, THE COMPLEX VALUES OF THE GRID ARE
!               STORED IN MEMORY AS FOLLOWS (USING U AND V TO
!               DENOTE THE REAL AND IMAGINARY COMPONENTS, AND
!               ASSUMING CONVENTIONAL FORTRAN STORAGE):
!
!   U(1,1), V(1,1), U(2,1), V(2,1),  ...  U(64,1), V(64,1), 4 NULLS,
!
!   U(1,2), V(1,2), U(2,2), V(2,2),  ...  U(64,2), V(64,2), 4 NULLS,
!
!   .       .       .       .         .   .        .        .
!   .       .       .       .         .   .        .        .
!   .       .       .       .         .   .        .        .
!
!   U(1,9), V(1,9), U(2,9), V(2,9),  ...  U(64,9), V(64,9), 4 NULLS.
!
!               WE CHOOSE (ARBITRARILY) TO TRANSORM FIRST ALONG THE
!               DIRECTION OF THE FIRST SUBSCRIPT.  THUS WE DEFINE
!               THE LENGTH OF THE TRANSFORMS, N, TO BE 64, THE
!               NUMBER OF TRANSFORMS, LOT, TO BE 9, THE INCREMENT
!               BETWEEN ELEMENTS OF EACH TRANSFORM, INC, TO BE 1,
!               AND THE INCREMENT BETWEEN THE STARTING POINTS
!               FOR EACH TRANSFORM, JUMP, TO BE 66 (THE FIRST
!               DIMENSION OF A).
!
!               CALL CFFT99( A, WORK, TRIGS1, IFAX1, 1, 66, 64, 9, -1)
!
!               TO TRANSFORM ALONG THE DIRECTION OF THE SECOND SUBSCRIPT
!               THE ROLES OF THE INCREMENTS ARE REVERSED.  THUS WE DEFIN
!               THE LENGTH OF THE TRANSFORMS, N, TO BE 9, THE
!               NUMBER OF TRANSFORMS, LOT, TO BE 64, THE INCREMENT
!               BETWEEN ELEMENTS OF EACH TRANSFORM, INC, TO BE 66,
!               AND THE INCREMENT BETWEEN THE STARTING POINTS
!               FOR EACH TRANSFORM, JUMP, TO BE 1
!
!               CALL CFFT99( A, WORK, TRIGS2, IFAX2, 66, 1, 9, 64, -1)
!
!               THESE TWO SEQUENTIAL STEPS RESULTS IN THE TWO-DIMENSIONA
!               FOURIER COEFFICIENT ARRAY OVERWRITING THE INPUT
!               GRIDPOINT ARRAY, A.  THE SAME TWO STEPS APPLIED AGAIN
!               WITH ISGN = +1 WOULD RESULT IN THE RECONSTRUCTION OF
!               THE GRIDPOINT ARRAY (MULTIPLIED BY A FACTOR OF 64*9).
!
!
!-----------------------------------------------------------------------
      DIMENSION AR(N),AI(N),WORK(N),TRIGS(N),IFAX(N)
!
!     SUBROUTINE "CFFT99" - MULTIPLE FAST COMPLEX FOURIER TRANSFORM
!
!     A IS THE ARRAY CONTAINING INPUT AND OUTPUT DATA
!     WORK IS AN AREA OF SIZE N*LOT
!     TRIGS IS A PREVIOUSLY PREPARED LIST OF TRIG FUNCTION VALUES
!     IFAX IS A PREVIOUSLY PREPARED LIST OF FACTORS OF N
!     INC IS THE INCREMENT WITHIN EACH DATA 'VECTOR'
!         (E.G. INC=1 FOR CONSECUTIVELY STORED DATA)
!     JUMP IS THE INCREMENT BETWEEN THE START OF EACH DATA VECTOR
!     N IS THE LENGTH OF THE DATA VECTORS
!     LOT IS THE NUMBER OF DATA VECTORS
!     ISGN = +1 FOR TRANSFORM FROM SPECTRAL TO GRIDPOINT
!           = -1 FOR TRANSFORM FROM GRIDPOINT TO SPECTRAL
!
!
!     VECTORIZATION IS ACHIEVED ON CRAY BY DOING THE TRANSFORMS IN
!     PARALLEL.
!
      NN = N+N
      INK=INC+INC
      JUM = JUMP+JUMP
      NFAX=IFAX(1)
      JNK = 2
      JST = 2
      IF (ISGN.GE.0) GO TO 30
!
!     THE INNERMOST TEMPERTON ROUTINES HAVE NO FACILITY FOR THE
!     FORWARD (ISGN = -1) TRANSFORM.  THEREFORE, THE INPUT MUST BE
!     REARRANGED AS FOLLOWS:
!
!     THE ORDER OF EACH INPUT VECTOR,
!
!     G(0), G(1), G(2), ... , G(N-2), G(N-1)
!
!     IS REVERSED (EXCLUDING G(0)) TO YIELD
!
!     G(0), G(N-1), G(N-2), ... , G(2), G(1).
!
!     WITHIN THE TRANSFORM, THE CORRESPONDING EXPONENTIAL MULTIPLIER
!     IS THEN PRECISELY THE CONJUGATE OF THAT FOR THE NORMAL
!     ORDERING.  THUS THE FORWARD (ISGN = -1) TRANSFORM IS
!     ACCOMPLISHED
!
!     FOR NFAX ODD, THE INPUT MUST BE TRANSFERRED TO THE WORK ARRAY,
!     AND THE REARRANGEMENT CAN BE DONE DURING THE MOVE.
!
      JNK = -2
      JST = NN-2
      IF (MOD(NFAX,2).EQ.1) GOTO 40
!
!     FOR NFAX EVEN, THE REARRANGEMENT MUST BE APPLIED DIRECTLY TO
!     THE INPUT ARRAY.  THIS CAN BE DONE BY SWAPPING ELEMENTS.
!
      IBASE = 1
      ILAST = (N-1)*INC
      NH = N/2
      DO 20 L=1,LOT
      I1 = IBASE+INC
      I2 = IBASE+ILAST
! 
      DO 10 M=1,NH
!     SWAP REAL AND IMAGINARY PORTIONS
      HREAL = AR(I1)
      HIMAG = AI(I1)
      AR(I1) = AR(I2)
      AI(I1) = AI(I2)
      AR(I2) = HREAL
      AI(I2) = HIMAG
      I1 = I1+INC
      I2 = I2-INC
   10 CONTINUE
      IBASE = IBASE+JUMP
   20 CONTINUE
      GOTO 100
!
   30 CONTINUE
      IF (MOD(NFAX,2).EQ.0) GOTO 100
!
   40 CONTINUE
!
!     DURING THE TRANSFORM PROCESS, NFAX STEPS ARE TAKEN, AND THE
!     RESULTS ARE STORED ALTERNATELY IN WORK AND IN A.  IF NFAX IS
!     ODD, THE INPUT DATA ARE FIRST MOVED TO WORK SO THAT THE FINAL
!     RESULT (AFTER NFAX STEPS) IS STORED IN ARRAY A.
!
      IBASE=1
      JBASE=1
      DO 60 L=1,LOT
!     MOVE REAL AND IMAGINARY PORTIONS OF ELEMENT ZERO
      WORK(JBASE) = AR(IBASE)
      WORK(JBASE+1) = AI(IBASE)
      I=IBASE+INC
      J=JBASE+JST
! 
      DO 50 M=2,N
!     MOVE REAL AND IMAGINARY PORTIONS OF OTHER ELEMENTS (POSSIBLY IN
!     REVERSE ORDER, DEPENDING ON JST AND JNK)
      WORK(J) = AR(I)
      WORK(J+1) = AI(I)
      I=I+INC
      J=J+JNK
   50 CONTINUE
      IBASE=IBASE+JUMP
      JBASE=JBASE+NN
   60 CONTINUE
!
  100 CONTINUE
!
!     PERFORM THE TRANSFORM PASSES, ONE PASS FOR EACH FACTOR.  DURING
!     EACH PASS THE DATA ARE MOVED FROM A TO WORK OR FROM WORK TO A.
!
!     FOR NFAX EVEN, THE FIRST PASS MOVES FROM A TO WORK
      IGO = 110
!     FOR NFAX ODD, THE FIRST PASS MOVES FROM WORK TO A
      IF (MOD(NFAX,2).EQ.1) IGO = 120
      LA=1
      DO 140 K=1,NFAX
      IF (IGO.EQ.120) GO TO 120
  110 CONTINUE
      CALL VPASSM(AR,AI,WORK(1),WORK(2),TRIGS,INC,2,JUMP,NN,LOT,N,IFAX(K+1),LA)
      IGO=120
      GO TO 130
  120 CONTINUE
      CALL VPASSM(WORK(1),WORK(2),AR,AI,TRIGS,2,INC,NN,JUMP,LOT,N,IFAX(K+1),LA)
      IGO=110
  130 CONTINUE
      LA=LA*IFAX(K+1)
  140 CONTINUE
!
!     AT THIS POINT THE FINAL TRANSFORM RESULT IS STORED IN A.
!
      RETURN
      END
!
!=======================================================================
!
      SUBROUTINE  SCOPY(N,SX,INCX,SY,INCY)

!
!     COPIES A VECTOR, X, TO A VECTOR, Y.
!     USES UNROLLED LOOPS FOR INCREMENTS EQUAL TO 1.
!     JACK DONGARRA, LINPACK, 3/11/78.
!
      REAL SX(1),SY(1)
      INTEGER I,INCX,INCY,IX,IY,M,MP1,N
!
      IF(N.LE.0)RETURN
      IF(INCX.EQ.1.AND.INCY.EQ.1)GO TO 20
!
!        CODE FOR UNEQUAL INCREMENTS OR EQUAL INCREMENTS
!          NOT EQUAL TO 1
!
      IX = 1
      IY = 1
      IF(INCX.LT.0)IX = (-N+1)*INCX + 1
      IF(INCY.LT.0)IY = (-N+1)*INCY + 1
      DO 10 I = 1,N
        SY(IY) = SX(IX)
        IX = IX + INCX
        IY = IY + INCY
   10 CONTINUE
      RETURN
!
!        CODE FOR BOTH INCREMENTS EQUAL TO 1
!
!
!        CLEAN-UP LOOP
!
   20 M = MOD(N,7)
      IF( M .EQ. 0 ) GO TO 40
      DO 30 I = 1,M
        SY(I) = SX(I)
   30 CONTINUE
      IF( N .LT. 7 ) RETURN
   40 MP1 = M + 1
      DO 50 I = MP1,N,7
        SY(I) = SX(I)
        SY(I + 1) = SX(I + 1)
        SY(I + 2) = SX(I + 2)
        SY(I + 3) = SX(I + 3)
        SY(I + 4) = SX(I + 4)
        SY(I + 5) = SX(I + 5)
        SY(I + 6) = SX(I + 6)
   50 CONTINUE
      RETURN
      END
!
!=======================================================================
!
      REAL FUNCTION SSUM(N,SX,INCX)
!
!        TAKES THE SUM OF THE VALUES OF A VECTOR.
!        USES UNROLLED LOOPS FOR INCREMENT EQUAL TO ONE.
!        JACK DONGARRA, LINPACK, 3/11/78.
!
      REAL SX(1),STEMP
      INTEGER I,INCX,M,MP1,N,IX
!
      SSUM = 0.0E0
      STEMP = 0.0E0
      IF (N .LT. 0) STOP
      IF (N .EQ. 0) RETURN
      IF (INCX .EQ. 1) GO TO 20
!
!        CODE FOR INCREMENT NOT EQUAL TO 1
!
      IX = 1
      IF (INCX .LE. 0) IX = (-N+1)*INCX + 1
      DO 10 I = 1,N
        STEMP = STEMP + SX(IX)
        IX = IX + INCX
   10 CONTINUE
      SSUM = STEMP
      RETURN
!
!        CODE FOR INCREMENT EQUAL TO 1
!
!        CLEAN-UP LOOP
!
   20 M = MOD(N,6)
      IF (M .EQ. 0) GO TO 40
      DO 30 I = 1,M
        STEMP = STEMP + SX(I)
   30 CONTINUE
      IF (N .LT. 6) GO TO 60
   40 MP1 = M + 1
      DO 50 I = MP1,N,6
        STEMP = STEMP + SX(I) + SX(I+1) + SX(I+2) + SX(I+3) + SX(I+4) + SX(I+5)
   50 CONTINUE
   60 SSUM = STEMP
      RETURN
      END
!
!=======================================================================
!
      SUBROUTINE CCOPY(N,CX,INCX,CY,INCY)

!
!     COPIES A VECTOR, X, TO A VECTOR, Y.
!     JACK DONGARRA, LINPACK, 3/11/78.
!
      COMPLEX CX(1),CY(1)
      INTEGER I,INCX,INCY,IX,IY,N
!
      IF(N.LE.0)RETURN
      IF(INCX.EQ.1.AND.INCY.EQ.1)GO TO 20
!
!        CODE FOR UNEQUAL INCREMENTS OR EQUAL INCREMENTS
!          NOT EQUAL TO 1
!
      IX = 1
      IY = 1
      IF(INCX.LT.0)IX = (-N+1)*INCX + 1
      IF(INCY.LT.0)IY = (-N+1)*INCY + 1
      DO 10 I = 1,N
        CY(IY) = CX(IX)
        IX = IX + INCX
        IY = IY + INCY
   10 CONTINUE
      RETURN
!
!        CODE FOR BOTH INCREMENTS EQUAL TO 1
!
   20 DO 30 I = 1,N
        CY(I) = CX(I)
   30 CONTINUE
      RETURN
      END
!
!=======================================================================
!
      REAL FUNCTION SNRM2 ( N, SX, INCX)
      INTEGER          NEXT
      REAL   SX(1),  CUTLO, CUTHI, HITEST, SUM, XMAX, ZERO, ONE
      DATA   ZERO, ONE /0.0E0, 1.0E0/
!
!     EUCLIDEAN NORM OF THE N-VECTOR STORED IN SX() WITH STORAGE
!     INCREMENT INCX .
!     IF    N .LE. 0 RETURN WITH RESULT = 0.
!     IF N .GE. 1 THEN INCX MUST BE .GE. 1
!
!           C.L.LAWSON, 1978 JAN 08
!
!     FOUR PHASE METHOD     USING TWO BUILT-IN CONSTANTS THAT ARE
!     HOPEFULLY APPLICABLE TO ALL MACHINES.
!         CUTLO = MAXIMUM OF  SQRT(U/EPS)  OVER ALL KNOWN MACHINES.
!         CUTHI = MINIMUM OF  SQRT(V)      OVER ALL KNOWN MACHINES.
!     WHERE
!         EPS = SMALLEST NO. SUCH THAT EPS + 1. .GT. 1.
!         U   = SMALLEST POSITIVE NO.   (UNDERFLOW LIMIT)
!         V   = LARGEST  NO.            (OVERFLOW  LIMIT)
!
!     BRIEF OUTLINE OF ALGORITHM..
!
!     PHASE 1    SCANS ZERO COMPONENTS.
!     MOVE TO PHASE 2 WHEN A COMPONENT IS NONZERO AND .LE. CUTLO
!     MOVE TO PHASE 3 WHEN A COMPONENT IS .GT. CUTLO
!     MOVE TO PHASE 4 WHEN A COMPONENT IS .GE. CUTHI/M
!     WHERE M = N FOR X() REAL AND M = 2*N FOR COMPLEX.
!
!     VALUES FOR CUTLO AND CUTHI..
!     FROM THE ENVIRONMENTAL PARAMETERS LISTED IN THE IMSL CONVERTER
!     DOCUMENT THE LIMITING VALUES ARE AS FOLLOWS..
!     CUTLO, S.P.   U/EPS = 2**(-102) FOR  HONEYWELL.  CLOSE SECONDS ARE
!                   UNIVAC AND DEC AT 2**(-103)
!                   THUS CUTLO = 2**(-51) = 4.44089E-16
!     CUTHI, S.P.   V = 2**127 FOR UNIVAC, HONEYWELL, AND DEC.
!                   THUS CUTHI = 2**(63.5) = 1.30438E19
!     CUTLO, D.P.   U/EPS = 2**(-67) FOR HONEYWELL AND DEC.
!                   THUS CUTLO = 2**(-33.5) = 8.23181D-11
!     CUTHI, D.P.   SAME AS S.P.  CUTHI = 1.30438D19
!     DATA CUTLO, CUTHI / 8.232D-11,  1.304D19 /
!     DATA CUTLO, CUTHI / 4.441E-16,  1.304E19 /
      DATA CUTLO, CUTHI / 4.441E-16,  1.304E19 /
!
      IF(N .GT. 0) GO TO 10
         SNRM2  = ZERO
         GO TO 300
!
   10 ASSIGN 30 TO NEXT
      SUM = ZERO
      NN = N * INCX
!                                                 BEGIN MAIN LOOP
      I = 1
   20    GO TO NEXT,(30, 50, 70, 110)
   30 IF( ABS(SX(I)) .GT. CUTLO) GO TO 85
      ASSIGN 50 TO NEXT
      XMAX = ZERO
!
!                        PHASE 1.  SUM IS ZERO
!
   50 IF( SX(I) .EQ. ZERO) GO TO 200
      IF( ABS(SX(I)) .GT. CUTLO) GO TO 85
!
!                                PREPARE FOR PHASE 2.
      ASSIGN 70 TO NEXT
      GO TO 105
!
!                                PREPARE FOR PHASE 4.
!
  100 I = J
      ASSIGN 110 TO NEXT
      SUM = (SUM / SX(I)) / SX(I)
  105 XMAX = ABS(SX(I))
      GO TO 115
!
!                   PHASE 2.  SUM IS SMALL.
!                             SCALE TO AVOID DESTRUCTIVE UNDERFLOW.
!
   70 IF( ABS(SX(I)) .GT. CUTLO ) GO TO 75
!
!                     COMMON CODE FOR PHASES 2 AND 4.
!                     IN PHASE 4 SUM IS LARGE.  SCALE TO AVOID OVERFLOW.
!
  110 IF( ABS(SX(I)) .LE. XMAX ) GO TO 115
         SUM = ONE + SUM * (XMAX / SX(I))**2
         XMAX = ABS(SX(I))
         GO TO 200
!
  115 SUM = SUM + (SX(I)/XMAX)**2
      GO TO 200
!
!
!                  PREPARE FOR PHASE 3.
!
   75 SUM = (SUM * XMAX) * XMAX
!
!
!     FOR REAL OR D.P. SET HITEST = CUTHI/N
!     FOR COMPLEX      SET HITEST = CUTHI/(2*N)
!
   85 HITEST = CUTHI/FLOAT( N )
!
!                   PHASE 3.  SUM IS MID-RANGE.  NO SCALING.
!
      DO 95 J =I,NN,INCX
      IF(ABS(SX(J)) .GE. HITEST) GO TO 100
   95    SUM = SUM + SX(J)**2
      SNRM2 = SQRT( SUM )
      GO TO 300
!
  200 CONTINUE
      I = I + INCX
      IF ( I .LE. NN ) GO TO 20
!
!              END OF MAIN LOOP.
!
!              COMPUTE SQUARE ROOT AND ADJUST FOR SCALING.
!
      SNRM2 = XMAX * SQRT(SUM)
  300 CONTINUE
      RETURN
      END
!
!=======================================================================
!
      SUBROUTINE FACT(N,IFAX)

!     FACTORIZATION ROUTINE THAT FIRST EXTRACTS ALL FACTORS OF 4
      DIMENSION IFAX(13)
      IF (N.GT.1) GO TO 10
      IFAX(1) = 0
      IF (N.LT.1) IFAX(1) = -99
      RETURN
   10 NN=N
      K=1
!     TEST FOR FACTORS OF 4
   20 IF (MOD(NN,4).NE.0) GO TO 30
      K=K+1
      IFAX(K)=4
      NN=NN/4
      IF (NN.EQ.1) GO TO 80
      GO TO 20
!     TEST FOR EXTRA FACTOR OF 2
   30 IF (MOD(NN,2).NE.0) GO TO 40
      K=K+1
      IFAX(K)=2
      NN=NN/2
      IF (NN.EQ.1) GO TO 80
!     TEST FOR FACTORS OF 3
   40 IF (MOD(NN,3).NE.0) GO TO 50
      K=K+1
      IFAX(K)=3
      NN=NN/3
      IF (NN.EQ.1) GO TO 80
      GO TO 40
!     NOW FIND REMAINING FACTORS
   50 L=5
      MAX = SQRT(FLOAT(NN))
      INC=2
!     INC ALTERNATELY TAKES ON VALUES 2 AND 4
   60 IF (MOD(NN,L).NE.0) GO TO 70
      K=K+1
      IFAX(K)=L
      NN=NN/L
      IF (NN.EQ.1) GO TO 80
      GO TO 60
   70 IF (L.GT.MAX) GO TO 75
      L=L+INC
      INC=6-INC
      GO TO 60
   75 K = K+1
      IFAX(K) = NN
   80 IFAX(1)=K-1
!     IFAX(1) NOW CONTAINS NUMBER OF FACTORS
      RETURN
      END
!
!=======================================================================
!
      SUBROUTINE FAX(IFAX,N,MODE)

      DIMENSION IFAX(10)
      NN=N
      IF (IABS(MODE).EQ.1) GO TO 10
      IF (IABS(MODE).EQ.8) GO TO 10
      NN=N/2
      IF ((NN+NN).EQ.N) GO TO 10
      IFAX(1)=-99
      RETURN
   10 K=1
!     TEST FOR FACTORS OF 4
   20 IF (MOD(NN,4).NE.0) GO TO 30
      K=K+1
      IFAX(K)=4
      NN=NN/4
      IF (NN.EQ.1) GO TO 80
      GO TO 20
!     TEST FOR EXTRA FACTOR OF 2
   30 IF (MOD(NN,2).NE.0) GO TO 40
      K=K+1
      IFAX(K)=2
      NN=NN/2
      IF (NN.EQ.1) GO TO 80
!     TEST FOR FACTORS OF 3
   40 IF (MOD(NN,3).NE.0) GO TO 50
      K=K+1
      IFAX(K)=3
      NN=NN/3
      IF (NN.EQ.1) GO TO 80
      GO TO 40
!     NOW FIND REMAINING FACTORS
   50 L=5
      INC=2
!     INC ALTERNATELY TAKES ON VALUES 2 AND 4
   60 IF (MOD(NN,L).NE.0) GO TO 70
      K=K+1
      IFAX(K)=L
      NN=NN/L
      IF (NN.EQ.1) GO TO 80
      GO TO 60
   70 L=L+INC
      INC=6-INC
      GO TO 60
   80 IFAX(1)=K-1
!     IFAX(1) CONTAINS NUMBER OF FACTORS
      NFAX=IFAX(1)
!     SORT FACTORS INTO ASCENDING ORDER
      IF (NFAX.EQ.1) GO TO 110
      DO 100 II=2,NFAX
      ISTOP=NFAX+2-II
      DO 90 I=2,ISTOP
      IF (IFAX(I+1).GE.IFAX(I)) GO TO 90
      ITEM=IFAX(I)
      IFAX(I)=IFAX(I+1)
      IFAX(I+1)=ITEM
   90 CONTINUE
  100 CONTINUE
  110 CONTINUE
      RETURN
      END
!
!=======================================================================
!
      SUBROUTINE FFT99A(A,WORK,TRIGS,INC,JUMP,N,LOT)

      DIMENSION A(N),WORK(N),TRIGS(N)
!
!     FFT99A - PREPROCESSING STEP, ISGN=+1
!     (SPECTRAL TO GRIDPOINT TRANSFORM)
!
      NH=N/2
      NX=N+1
      INK=INC+INC
!
!     A(0) AND A(N/2)
      IA=1
      IB=N*INC+1
      JA=1
      JB=2
! 
      DO 10 L=1,LOT
      WORK(JA)=A(IA)+A(IB)
      WORK(JB)=A(IA)-A(IB)
      IA=IA+JUMP
      IB=IB+JUMP
      JA=JA+NX
      JB=JB+NX
   10 CONTINUE
!
!     REMAINING WAVENUMBERS
      IABASE=2*INC+1
      IBBASE=(N-2)*INC+1
      JABASE=3
      JBBASE=N-1
!
      DO 30 K=3,NH,2
      IA=IABASE
      IB=IBBASE
      JA=JABASE
      JB=JBBASE
      C=TRIGS(N+K)
      S=TRIGS(N+K+1)
! 
      DO 20 L=1,LOT
      WORK(JA)=(A(IA)+A(IB))-(S*(A(IA)-A(IB))+C*(A(IA+INC)+A(IB+INC)))
      WORK(JB)=(A(IA)+A(IB))+(S*(A(IA)-A(IB))+C*(A(IA+INC)+A(IB+INC)))
      WORK(JA+1)=(C*(A(IA)-A(IB))-S*(A(IA+INC)+A(IB+INC)))+(A(IA+INC)-A(IB+INC))
      WORK(JB+1)=(C*(A(IA)-A(IB))-S*(A(IA+INC)+A(IB+INC)))-(A(IA+INC)-A(IB+INC))
      IA=IA+JUMP
      IB=IB+JUMP
      JA=JA+NX
      JB=JB+NX
   20 CONTINUE
      IABASE=IABASE+INK
      IBBASE=IBBASE-INK
      JABASE=JABASE+2
      JBBASE=JBBASE-2
   30 CONTINUE
!
      IF (IABASE.NE.IBBASE) GO TO 50
!     WAVENUMBER N/4 (IF IT EXISTS)
      IA=IABASE
      JA=JABASE
! 
      DO 40 L=1,LOT
      WORK(JA)=2.0*A(IA)
      WORK(JA+1)=-2.0*A(IA+INC)
      IA=IA+JUMP
      JA=JA+NX
   40 CONTINUE
!
   50 CONTINUE
      RETURN
      END
!
!=======================================================================
!
      SUBROUTINE FFT99B(WORK,A,TRIGS,INC,JUMP,N,LOT)

      DIMENSION WORK(N),A(N),TRIGS(N)
!
!     FFT99B - POSTPROCESSING STEP, ISGN=-1
!     (GRIDPOINT TO SPECTRAL TRANSFORM)
!
      NH=N/2
      NX=N+1
      INK=INC+INC
!
!     A(0) AND A(N/2)
      SCALE=1.0/FLOAT(N)
      IA=1
      IB=2
      JA=1
      JB=N*INC+1
! 
      DO 10 L=1,LOT
      A(JA)=SCALE*(WORK(IA)+WORK(IB))
      A(JB)=SCALE*(WORK(IA)-WORK(IB))
      A(JA+INC)=0.0
      A(JB+INC)=0.0
      IA=IA+NX
      IB=IB+NX
      JA=JA+JUMP
      JB=JB+JUMP
   10 CONTINUE
!
!     REMAINING WAVENUMBERS
      SCALE=0.5*SCALE
      IABASE=3
      IBBASE=N-1
      JABASE=2*INC+1
      JBBASE=(N-2)*INC+1
!
      DO 30 K=3,NH,2
      IA=IABASE
      IB=IBBASE
      JA=JABASE
      JB=JBBASE
      C=TRIGS(N+K)
      S=TRIGS(N+K+1)
! 
      DO 20 L=1,LOT
      A(JA)=SCALE*((WORK(IA)+WORK(IB))+(C*(WORK(IA+1)+WORK(IB+1))+S*(WORK(IA)-WORK(IB))))
      A(JB)=SCALE*((WORK(IA)+WORK(IB))-(C*(WORK(IA+1)+WORK(IB+1))+S*(WORK(IA)-WORK(IB))))
      A(JA+INC)=SCALE*((C*(WORK(IA)-WORK(IB))-S*(WORK(IA+1)+WORK(IB+1)))+(WORK(IB+1)-WORK(IA+1)))
      A(JB+INC)=SCALE*((C*(WORK(IA)-WORK(IB))-S*(WORK(IA+1)+WORK(IB+1)))-(WORK(IB+1)-WORK(IA+1)))
      IA=IA+NX
      IB=IB+NX
      JA=JA+JUMP
      JB=JB+JUMP
   20 CONTINUE
      IABASE=IABASE+2
      IBBASE=IBBASE-2
      JABASE=JABASE+INK
      JBBASE=JBBASE-INK
   30 CONTINUE
!
      IF (IABASE.NE.IBBASE) GO TO 50
!     WAVENUMBER N/4 (IF IT EXISTS)
      IA=IABASE
      JA=JABASE
      SCALE=2.0*SCALE
! 
      DO 40 L=1,LOT
      A(JA)=SCALE*WORK(IA)
      A(JA+INC)=-SCALE*WORK(IA+1)
      IA=IA+NX
      JA=JA+JUMP
   40 CONTINUE
!
   50 CONTINUE
      RETURN
      END
!
!=======================================================================
!
      SUBROUTINE VPASSM(A,B,C,D,TRIGS,INC1,INC2,INC3,INC4,LOT,N,IFAC,LA)

      DIMENSION A(N),B(N),C(N),D(N),TRIGS(N)
!
!     "VPASSM" - MULTIPLE VERSION OF "VPASSA"
!     PERFORMS ONE PASS THROUGH DATA
!     AS PART OF MULTIPLE COMPLEX FFT ROUTINE
!     A IS FIRST REAL INPUT VECTOR
!     B IS FIRST IMAGINARY INPUT VECTOR
!     C IS FIRST REAL OUTPUT VECTOR
!     D IS FIRST IMAGINARY OUTPUT VECTOR
!     TRIGS IS PRECALCULATED TABLE OF SINES " COSINES
!     INC1 IS ADDRESSING INCREMENT FOR A AND B
!     INC2 IS ADDRESSING INCREMENT FOR C AND D
!     INC3 IS ADDRESSING INCREMENT BETWEEN A"S & B"S
!     INC4 IS ADDRESSING INCREMENT BETWEEN C"S & D"S
!     LOT IS THE NUMBER OF VECTORS
!     N IS LENGTH OF VECTORS
!     IFAC IS CURRENT FACTOR OF N
!     LA IS PRODUCT OF PREVIOUS FACTORS
!
      DATA SIN36/0.587785252292473/,COS36/0.809016994374947/,SIN72/0.951056516295154/,COS72/0.309016994374947/,   &
           SIN60/0.866025403784437/
!
      M=N/IFAC
      IINK=M*INC1
      JINK=LA*INC2
      JUMP=(IFAC-1)*JINK
      IBASE=0
      JBASE=0
      IGO=IFAC-1
      IF (IGO.GT.4) RETURN
      GO TO (10,50,90,130),IGO
!
!     CODING FOR FACTOR 2
!
   10 IA=1
      JA=1
      IB=IA+IINK
      JB=JA+JINK
      DO 20 L=1,LA
      I=IBASE
      J=JBASE
! 
      DO 15 IJK=1,LOT
      C(JA+J)=A(IA+I)+A(IB+I)
      D(JA+J)=B(IA+I)+B(IB+I)
      C(JB+J)=A(IA+I)-A(IB+I)
      D(JB+J)=B(IA+I)-B(IB+I)
      I=I+INC3
      J=J+INC4
   15 CONTINUE
      IBASE=IBASE+INC1
      JBASE=JBASE+INC2
   20 CONTINUE
      IF (LA.EQ.M) RETURN
      LA1=LA+1
      JBASE=JBASE+JUMP
      DO 40 K=LA1,M,LA
      KB=K+K-2
      C1=TRIGS(KB+1)
      S1=TRIGS(KB+2)
      DO 30 L=1,LA
      I=IBASE
      J=JBASE
! 
      DO 25 IJK=1,LOT
      C(JA+J)=A(IA+I)+A(IB+I)
      D(JA+J)=B(IA+I)+B(IB+I)
      C(JB+J)=C1*(A(IA+I)-A(IB+I))-S1*(B(IA+I)-B(IB+I))
      D(JB+J)=S1*(A(IA+I)-A(IB+I))+C1*(B(IA+I)-B(IB+I))
      I=I+INC3
      J=J+INC4
   25 CONTINUE
      IBASE=IBASE+INC1
      JBASE=JBASE+INC2
   30 CONTINUE
      JBASE=JBASE+JUMP
   40 CONTINUE
      RETURN
!
!     CODING FOR FACTOR 3
!
   50 IA=1
      JA=1
      IB=IA+IINK
      JB=JA+JINK
      IC=IB+IINK
      JC=JB+JINK
      DO 60 L=1,LA
      I=IBASE
      J=JBASE
! 
      DO 55 IJK=1,LOT
      C(JA+J)=A(IA+I)+(A(IB+I)+A(IC+I))
      D(JA+J)=B(IA+I)+(B(IB+I)+B(IC+I))
      C(JB+J)=(A(IA+I)-0.5*(A(IB+I)+A(IC+I)))-(SIN60*(B(IB+I)-B(IC+I)))
      C(JC+J)=(A(IA+I)-0.5*(A(IB+I)+A(IC+I)))+(SIN60*(B(IB+I)-B(IC+I)))
      D(JB+J)=(B(IA+I)-0.5*(B(IB+I)+B(IC+I)))+(SIN60*(A(IB+I)-A(IC+I)))
      D(JC+J)=(B(IA+I)-0.5*(B(IB+I)+B(IC+I)))-(SIN60*(A(IB+I)-A(IC+I)))
      I=I+INC3
      J=J+INC4
   55 CONTINUE
      IBASE=IBASE+INC1
      JBASE=JBASE+INC2
   60 CONTINUE
      IF (LA.EQ.M) RETURN
      LA1=LA+1
      JBASE=JBASE+JUMP
      DO 80 K=LA1,M,LA
      KB=K+K-2
      KC=KB+KB
      C1=TRIGS(KB+1)
      S1=TRIGS(KB+2)
      C2=TRIGS(KC+1)
      S2=TRIGS(KC+2)
      DO 70 L=1,LA
      I=IBASE
      J=JBASE
! 
      DO 65 IJK=1,LOT
      C(JA+J)=A(IA+I)+(A(IB+I)+A(IC+I))
      D(JA+J)=B(IA+I)+(B(IB+I)+B(IC+I))
      C(JB+J)=C1*((A(IA+I)-0.5*(A(IB+I)+A(IC+I)))-(SIN60*(B(IB+I)-B(IC+I))))   &
         -S1*((B(IA+I)-0.5*(B(IB+I)+B(IC+I)))+(SIN60*(A(IB+I)-A(IC+I))))
      D(JB+J)=S1*((A(IA+I)-0.5*(A(IB+I)+A(IC+I)))-(SIN60*(B(IB+I)-B(IC+I))))   &
         +C1*((B(IA+I)-0.5*(B(IB+I)+B(IC+I)))+(SIN60*(A(IB+I)-A(IC+I))))
      C(JC+J)=C2*((A(IA+I)-0.5*(A(IB+I)+A(IC+I)))+(SIN60*(B(IB+I)-B(IC+I))))   &
         -S2*((B(IA+I)-0.5*(B(IB+I)+B(IC+I)))-(SIN60*(A(IB+I)-A(IC+I))))
      D(JC+J)=S2*((A(IA+I)-0.5*(A(IB+I)+A(IC+I)))+(SIN60*(B(IB+I)-B(IC+I))))   &
         +C2*((B(IA+I)-0.5*(B(IB+I)+B(IC+I)))-(SIN60*(A(IB+I)-A(IC+I))))
      I=I+INC3
      J=J+INC4
   65 CONTINUE
      IBASE=IBASE+INC1
      JBASE=JBASE+INC2
   70 CONTINUE
      JBASE=JBASE+JUMP
   80 CONTINUE
      RETURN
!
!     CODING FOR FACTOR 4
!
   90 IA=1
      JA=1
      IB=IA+IINK
      JB=JA+JINK
      IC=IB+IINK
      JC=JB+JINK
      ID=IC+IINK
      JD=JC+JINK
      DO 100 L=1,LA
      I=IBASE
      J=JBASE
! 
      DO 95 IJK=1,LOT
      C(JA+J)=(A(IA+I)+A(IC+I))+(A(IB+I)+A(ID+I))
      C(JC+J)=(A(IA+I)+A(IC+I))-(A(IB+I)+A(ID+I))
      D(JA+J)=(B(IA+I)+B(IC+I))+(B(IB+I)+B(ID+I))
      D(JC+J)=(B(IA+I)+B(IC+I))-(B(IB+I)+B(ID+I))
      C(JB+J)=(A(IA+I)-A(IC+I))-(B(IB+I)-B(ID+I))
      C(JD+J)=(A(IA+I)-A(IC+I))+(B(IB+I)-B(ID+I))
      D(JB+J)=(B(IA+I)-B(IC+I))+(A(IB+I)-A(ID+I))
      D(JD+J)=(B(IA+I)-B(IC+I))-(A(IB+I)-A(ID+I))
      I=I+INC3
      J=J+INC4
   95 CONTINUE
      IBASE=IBASE+INC1
      JBASE=JBASE+INC2
  100 CONTINUE
      IF (LA.EQ.M) RETURN
      LA1=LA+1
      JBASE=JBASE+JUMP
      DO 120 K=LA1,M,LA
      KB=K+K-2
      KC=KB+KB
      KD=KC+KB
      C1=TRIGS(KB+1)
      S1=TRIGS(KB+2)
      C2=TRIGS(KC+1)
      S2=TRIGS(KC+2)
      C3=TRIGS(KD+1)
      S3=TRIGS(KD+2)
      DO 110 L=1,LA
      I=IBASE
      J=JBASE
! 
      DO 105 IJK=1,LOT
      C(JA+J)=(A(IA+I)+A(IC+I))+(A(IB+I)+A(ID+I))
      D(JA+J)=(B(IA+I)+B(IC+I))+(B(IB+I)+B(ID+I))
      C(JC+J)=C2*((A(IA+I)+A(IC+I))-(A(IB+I)+A(ID+I)))   &
         -S2*((B(IA+I)+B(IC+I))-(B(IB+I)+B(ID+I)))
      D(JC+J)=S2*((A(IA+I)+A(IC+I))-(A(IB+I)+A(ID+I)))   &
         +C2*((B(IA+I)+B(IC+I))-(B(IB+I)+B(ID+I)))
      C(JB+J)=C1*((A(IA+I)-A(IC+I))-(B(IB+I)-B(ID+I)))   &
         -S1*((B(IA+I)-B(IC+I))+(A(IB+I)-A(ID+I)))
      D(JB+J)=S1*((A(IA+I)-A(IC+I))-(B(IB+I)-B(ID+I)))   &
         +C1*((B(IA+I)-B(IC+I))+(A(IB+I)-A(ID+I)))
      C(JD+J)=C3*((A(IA+I)-A(IC+I))+(B(IB+I)-B(ID+I)))   &
         -S3*((B(IA+I)-B(IC+I))-(A(IB+I)-A(ID+I)))
      D(JD+J)=S3*((A(IA+I)-A(IC+I))+(B(IB+I)-B(ID+I)))   &
         +C3*((B(IA+I)-B(IC+I))-(A(IB+I)-A(ID+I)))
      I=I+INC3
      J=J+INC4
  105 CONTINUE
      IBASE=IBASE+INC1
      JBASE=JBASE+INC2
  110 CONTINUE
      JBASE=JBASE+JUMP
  120 CONTINUE
      RETURN
!
!     CODING FOR FACTOR 5
!
  130 IA=1
      JA=1
      IB=IA+IINK
      JB=JA+JINK
      IC=IB+IINK
      JC=JB+JINK
      ID=IC+IINK
      JD=JC+JINK
      IE=ID+IINK
      JE=JD+JINK
      DO 140 L=1,LA
      I=IBASE
      J=JBASE
! 
      DO 135 IJK=1,LOT
      C(JA+J)=A(IA+I)+(A(IB+I)+A(IE+I))+(A(IC+I)+A(ID+I))
      D(JA+J)=B(IA+I)+(B(IB+I)+B(IE+I))+(B(IC+I)+B(ID+I))
      C(JB+J)=(A(IA+I)+COS72*(A(IB+I)+A(IE+I))-COS36*(A(IC+I)+A(ID+I)))   &
        -(SIN72*(B(IB+I)-B(IE+I))+SIN36*(B(IC+I)-B(ID+I)))
      C(JE+J)=(A(IA+I)+COS72*(A(IB+I)+A(IE+I))-COS36*(A(IC+I)+A(ID+I)))   &
        +(SIN72*(B(IB+I)-B(IE+I))+SIN36*(B(IC+I)-B(ID+I)))
      D(JB+J)=(B(IA+I)+COS72*(B(IB+I)+B(IE+I))-COS36*(B(IC+I)+B(ID+I)))   &
        +(SIN72*(A(IB+I)-A(IE+I))+SIN36*(A(IC+I)-A(ID+I)))
      D(JE+J)=(B(IA+I)+COS72*(B(IB+I)+B(IE+I))-COS36*(B(IC+I)+B(ID+I)))   &
        -(SIN72*(A(IB+I)-A(IE+I))+SIN36*(A(IC+I)-A(ID+I)))
      C(JC+J)=(A(IA+I)-COS36*(A(IB+I)+A(IE+I))+COS72*(A(IC+I)+A(ID+I)))   &
        -(SIN36*(B(IB+I)-B(IE+I))-SIN72*(B(IC+I)-B(ID+I)))
      C(JD+J)=(A(IA+I)-COS36*(A(IB+I)+A(IE+I))+COS72*(A(IC+I)+A(ID+I)))   &
        +(SIN36*(B(IB+I)-B(IE+I))-SIN72*(B(IC+I)-B(ID+I)))
      D(JC+J)=(B(IA+I)-COS36*(B(IB+I)+B(IE+I))+COS72*(B(IC+I)+B(ID+I)))   &
        +(SIN36*(A(IB+I)-A(IE+I))-SIN72*(A(IC+I)-A(ID+I)))
      D(JD+J)=(B(IA+I)-COS36*(B(IB+I)+B(IE+I))+COS72*(B(IC+I)+B(ID+I)))   &
        -(SIN36*(A(IB+I)-A(IE+I))-SIN72*(A(IC+I)-A(ID+I)))
      I=I+INC3
      J=J+INC4
  135 CONTINUE
      IBASE=IBASE+INC1
      JBASE=JBASE+INC2
  140 CONTINUE
      IF (LA.EQ.M) RETURN
      LA1=LA+1
      JBASE=JBASE+JUMP
      DO 160 K=LA1,M,LA
      KB=K+K-2
      KC=KB+KB
      KD=KC+KB
      KE=KD+KB
      C1=TRIGS(KB+1)
      S1=TRIGS(KB+2)
      C2=TRIGS(KC+1)
      S2=TRIGS(KC+2)
      C3=TRIGS(KD+1)
      S3=TRIGS(KD+2)
      C4=TRIGS(KE+1)
      S4=TRIGS(KE+2)
      DO 150 L=1,LA
      I=IBASE
      J=JBASE
! 
      DO 145 IJK=1,LOT
      C(JA+J)=A(IA+I)+(A(IB+I)+A(IE+I))+(A(IC+I)+A(ID+I))
      D(JA+J)=B(IA+I)+(B(IB+I)+B(IE+I))+(B(IC+I)+B(ID+I))
      C(JB+J)=C1*((A(IA+I)+COS72*(A(IB+I)+A(IE+I))-COS36*(A(IC+I)+A(ID+I)))   &
            -(SIN72*(B(IB+I)-B(IE+I))+SIN36*(B(IC+I)-B(ID+I))))               &
         -S1*((B(IA+I)+COS72*(B(IB+I)+B(IE+I))-COS36*(B(IC+I)+B(ID+I)))       &
            +(SIN72*(A(IB+I)-A(IE+I))+SIN36*(A(IC+I)-A(ID+I))))
      D(JB+J)=S1*((A(IA+I)+COS72*(A(IB+I)+A(IE+I))-COS36*(A(IC+I)+A(ID+I)))   &
            -(SIN72*(B(IB+I)-B(IE+I))+SIN36*(B(IC+I)-B(ID+I))))               & 
         +C1*((B(IA+I)+COS72*(B(IB+I)+B(IE+I))-COS36*(B(IC+I)+B(ID+I)))       & 
            +(SIN72*(A(IB+I)-A(IE+I))+SIN36*(A(IC+I)-A(ID+I))))
      C(JE+J)=C4*((A(IA+I)+COS72*(A(IB+I)+A(IE+I))-COS36*(A(IC+I)+A(ID+I)))   &
            +(SIN72*(B(IB+I)-B(IE+I))+SIN36*(B(IC+I)-B(ID+I))))               &
         -S4*((B(IA+I)+COS72*(B(IB+I)+B(IE+I))-COS36*(B(IC+I)+B(ID+I)))       &
            -(SIN72*(A(IB+I)-A(IE+I))+SIN36*(A(IC+I)-A(ID+I))))
      D(JE+J)=S4*((A(IA+I)+COS72*(A(IB+I)+A(IE+I))-COS36*(A(IC+I)+A(ID+I)))   &
            +(SIN72*(B(IB+I)-B(IE+I))+SIN36*(B(IC+I)-B(ID+I))))               &
         +C4*((B(IA+I)+COS72*(B(IB+I)+B(IE+I))-COS36*(B(IC+I)+B(ID+I)))       &
            -(SIN72*(A(IB+I)-A(IE+I))+SIN36*(A(IC+I)-A(ID+I))))
      C(JC+J)=C2*((A(IA+I)-COS36*(A(IB+I)+A(IE+I))+COS72*(A(IC+I)+A(ID+I)))   &
            -(SIN36*(B(IB+I)-B(IE+I))-SIN72*(B(IC+I)-B(ID+I))))               &
         -S2*((B(IA+I)-COS36*(B(IB+I)+B(IE+I))+COS72*(B(IC+I)+B(ID+I)))       &
            +(SIN36*(A(IB+I)-A(IE+I))-SIN72*(A(IC+I)-A(ID+I))))
      D(JC+J)=S2*((A(IA+I)-COS36*(A(IB+I)+A(IE+I))+COS72*(A(IC+I)+A(ID+I)))   &
            -(SIN36*(B(IB+I)-B(IE+I))-SIN72*(B(IC+I)-B(ID+I))))               &
         +C2*((B(IA+I)-COS36*(B(IB+I)+B(IE+I))+COS72*(B(IC+I)+B(ID+I)))       &
            +(SIN36*(A(IB+I)-A(IE+I))-SIN72*(A(IC+I)-A(ID+I))))
      C(JD+J)=C3*((A(IA+I)-COS36*(A(IB+I)+A(IE+I))+COS72*(A(IC+I)+A(ID+I)))   &
            +(SIN36*(B(IB+I)-B(IE+I))-SIN72*(B(IC+I)-B(ID+I))))               &  
         -S3*((B(IA+I)-COS36*(B(IB+I)+B(IE+I))+COS72*(B(IC+I)+B(ID+I)))       &
            -(SIN36*(A(IB+I)-A(IE+I))-SIN72*(A(IC+I)-A(ID+I))))
      D(JD+J)=S3*((A(IA+I)-COS36*(A(IB+I)+A(IE+I))+COS72*(A(IC+I)+A(ID+I)))   &
            +(SIN36*(B(IB+I)-B(IE+I))-SIN72*(B(IC+I)-B(ID+I))))               &
         +C3*((B(IA+I)-COS36*(B(IB+I)+B(IE+I))+COS72*(B(IC+I)+B(ID+I)))       &
            -(SIN36*(A(IB+I)-A(IE+I))-SIN72*(A(IC+I)-A(ID+I))))
      I=I+INC3
      J=J+INC4
  145 CONTINUE
      IBASE=IBASE+INC1
      JBASE=JBASE+INC2
  150 CONTINUE
      JBASE=JBASE+JUMP
  160 CONTINUE
      RETURN
      END

      subroutine rcsexp(lfft,ifax,trig)
      dimension trig(1)
      integer ifax(1)
      call fftfax(lfft,ifax,trig)
      return
      end subroutine

      subroutine ccexp(lfft,ifax,trig)
      dimension trig(1)
      integer ifax(1)
      call cftfax(lfft,ifax,trig)
      return
      end subroutine
      subroutine rcosexp(lfft,sine,cosine,ifax,trig)
      dimension sine(1),cosine(1)
      dimension trig(1)
      integer ifax(1)
      
      tpi=8.*atan(1.)
      xnum = tpi/float(2*lfft)

      call fftfax(lfft,ifax,trig)
      do i=1,lfft
        sine(i)   = sin(float(i-1)*xnum)
        cosine(i) = cos(float(i-1)*xnum)
      end do 
      return
      end subroutine
      subroutine ccosexp(lfft,sine,cosine,ifax,trig)
      dimension sine(1),cosine(1)
      dimension trig(1)
      integer ifax(1) 

      tpi=8.*atan(1.)
      xnum = tpi/float(2*lfft)

      call cftfax(lfft,ifax,trig)
      do i=1,lfft
        sine(i)   = sin(float(i-1)*xnum)
        cosine(i) = cos(float(i-1)*xnum)
      end do 
      return
      end subroutine
      subroutine rsinexp(lfft,sine,cosine,ifax,trig)
      dimension sine(1),cosine(1)
      dimension trig(1)
      integer ifax(1) 

      tpi=8.*atan(1.)
      xnum = tpi/float(2*lfft)

      call fftfax(lfft,ifax,trig)
      do i=1,lfft
        sine(i)   = sin(float(i-1)*xnum)
        cosine(i) = cos(float(i-1)*xnum)
      end do 
       
      end subroutine
      subroutine csinexp(lfft,sine,cosine,ifax,trig)
      dimension sine(1),cosine(1)
      dimension trig(1)
      integer ifax(1)

      tpi=8.*atan(1.)
      xnum = tpi/float(2*lfft)

      call cftfax(lfft,ifax,trig)
      do i=1,lfft
        sine(i)   = sin(float(i-1)*xnum)
        cosine(i) = cos(float(i-1)*xnum)
      end do 
       
      end subroutine


      subroutine csr(a,b,work,lfft,n,ifax,trig)
      dimension a(1),b(1),work(1)
      integer ifax(1)
      dimension trig(1)
!cccccccccc   symmetric complex to real  fft  
!     Note: lfft must be even. n is the number of ffts to be 
!     performed. The input , a , is a packed array and the 
!     output is returned in packed form to a . 
!       The inverse of csr is rcs. IS = +1 always for csr. 
!ccccccc    first , unpack the data  cccccccccccccccccccccccccccc

      do k=1,n 

         loca = (k-1)*lfft       ! starting location or a
         locb = (k-1)*(lfft+2)   ! starting location for b
         loc1 = locb  + 2        ! imag part of first mode 
         loc2 = locb  + (lfft+1) ! real part of last mode 

         do i=1,lfft

            b(i+locb) = a(i+loca)

         end do


         b(loc2) = b(loc1)  !real part of last mode equals
         !imaginary part of first mode

         b(loc1) = 0.0   !imag part of first mode is zero
         b(loc2+1) = 0.0 !imag part of last mode is zero

      end do
!ccccccc    now perform the complex to real transform on b  cccccc
            is =   +1
            inc =   1
            jump =   lfft +2
            lot =   n 

         call rfftmlt(b,work,trig,ifax,inc,jump,lfft,lot,is)

!ccccccc   now pack the result  ccccccccccccccccccccccccccccccccccc
         do k =1,n 

            loca = (k-1)*lfft       ! starting location for a  
            locb = (k-1)*(lfft+2)   ! starting location for b 

            do i=1,lfft

               a(i+loca)  =  b(i+locb) 

            end do
         end do
           
      end subroutine 
      subroutine ccfft(a,br,bi,work,lfft,n,is,ifax,trig)
      dimension a(1),br(1),bi(1),work(1)
      integer ifax(1)
      dimension trig(1)
!cccccccccc  complex to complex fft  
!    Note: lfft must be even and designates the fft length in (complex)
!    word pairs. N is the number of ffts to be performed.
!    performed. 
!    IS = -1 or +1 for forward or inverse . 
!ccccccc    Must avoid jump being a multiple of 2 ccccccccccccccccccc
              lfft2 = 2*lfft
              
              
              if (is.eq.-1) then 
                 fac = 1.0
              else
                 fac = (1./float(lfft))
              endif


              do k=1,n 
                 loca = (k-1)*lfft2     ! starting location or a
                 locb = (k-1)*(lfft+1)  ! starting location for br and bi

                 do i=1,lfft
                    locare  = loca + (2*i-1)
                    locaim   = locare+1

                    br(locb +i) = a(locare)
                    bi(locb +i) = a(locaim)
                 end do
              end do
!ccccccc    now perform complex to complex transform on b  cccccc
            inc =   1
            jump =   lfft +1       
            lot = n

         call cfftmlt(br,bi,work,trig,ifax,inc,jump,lfft,lot,is)
                           
!ccccccc   now pack the result  ccccccccccccccccccccccccccccccccccc
         do k=1,n 
            loca = (k-1)*lfft2     ! starting location or a
            locb = (k-1)*(lfft+1)  ! starting location for br and bi 

            do i=1,lfft
               locare  = loca + (2*i-1)
               locaim   = locare+1

               a(locare)=br(locb +i)*fac
               a(locaim)=bi(locb +i)*fac

            end do
         end do
        
      end subroutine 
      subroutine ccheb(a,b,br,bi,work,sumre,sumim,lfft,n,is,ifax,trig,sine,cosine)
      dimension a(1),b(1),br(1),bi(1),work(1)
      dimension trig(1),sine(1),cosine(1)
      dimension sumre(1),sumim(1)
      integer ifax(1) 
!           is = -1 real to Chebyshev
!           is = +1 Chebyshev to real                    
           lfftd2 = lfft/2
           lfft2  = 2*lfft
           fac1 = 4.0
           fac2 = 1./(8.*(float(lfft)))
!ccccccccc     preprocess the complex Chebyshev coefficients  cc
             if (is.eq.+1) then 

                
                do k=1,n

                   loca = (k-1)*(lfft2+2)

                   a(loca+1) = 2.0*a(loca+1)
                   a(loca+2) = 2.0*a(loca+2)
                   a(loca+lfft2+1) = 2.0*a(loca+lfft2+1)
                   a(loca+lfft2+2) = 2.0*a(loca+lfft2+2)

                   do  i=1, lfft+1

                      locr =  loca + 2*i-1 
                      locim = locr+1


                      a(locr)= a(locr)*fac2
                      a(locim) = a(locim)*fac2

                   end do
                end do

        
             endif
!ccccccccc Compute complex cosine transform  cccccccccccccccccccccccc
                     ! initialize sums 
             
             do k=1,n 
                sumr = 0.0 
                sumi = 0.0
                loca = (k-1)*(lfft2+2)  ! locations for a 
                loc1r = loca +1          ! first real loc for a 
                loc2r = loca +(lfft2+1) ! last real location for a   
                loc1i = loca +2          ! first imag. loc for a 
                loc2i = loca +(lfft2+2) ! last imag loc for a

                do i=2,lfft

                   locar = loca + (2*i-1)  ! real locs for a 
                   locai = locar +1        ! imag locs for a 

                   sumr = sumr + a(locar)*cosine(i)
                   sumi = sumi + a(locai)*cosine(i)
                end do

                sumre(k) = sumr + a(loc1r)/2. - a(loc2r)/2. 
                sumim(k) = sumi + a(loc1i)/2. - a(loc2i)/2. 
             end do

             !  pre-processing 

             do k=1,n 

                loca = (k-1)*(lfft2+2)    ! locs for a 
                locb = (k-1)*(lfft+1)    ! locs for br and bi 

                do i=1,lfft
                   ir = 2*i-1
                   im = ir+1

                   ilowr  =  loca +ir                  ! 1:2l-1
                   iupr   =  loca +(lfft2+2) - ir     ! 2l+1:3
                   ilowim =  loca + im                 ! 2:2l
                   iupim  =  loca +(lfft2+2) - (im-2) ! 2l+2:4


                   br(locb +i) = (a(ilowr)+a(iupr))-2.0*sine(i)*(a(ilowr)-a(iupr))

                   bi(locb +i) = (a(ilowim)+a(iupim))-2.0*sine(i)*(a(ilowim)-a(iupim))

                end do
             end do

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
              ! complex to complex transform 

                       isign = -1
                       inc = 1 
                       jump = lfft +1 
                       lot = n
       call cfftmlt(br,bi,work,trig,ifax,inc,jump,lfft,lot,isign)    

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccc
       
       do k=1,n 
          loca =   (k-1)*(lfft2+2)     ! starting location or a
          locb = (k-1)*(lfft+1)  ! starting location for real part of b

          do i=1,lfft
             locare  = loca + (2*i-1)
             locaim   = locare+1

             b(locare)=br(locb +i)
             b(locaim)=bi(locb +i)
          end do
       end do

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
                        ! post-processing 

       do k=1,n 



          loca = (k-1)*(lfft2+2)
          locb = (k-1)*(lfft2+2)



          ! fill real part of a 

          a(loca+1) = 2.* b(locb+1)
          a(loca+3) = fac1*sumre(k)
          a(loca+(lfft2+1)) = 2.*b(locb+lfft+1)

          ! fill imag part of a 

          a(loca+2) = 2.* b(locb+2)
          a(loca+4) = fac1*sumim(k)
          a(loca+(lfft2+2)) = 2.*b(locb+lfft+2)
       end do

               ! fill rest of a 





       do i=2,lfftd2
          ir = 2*i-1
          im = ir+1


          ir1 = 2*ir-1    ! 5:2l-3
          im1 = ir1+2     ! 7:2l-1
          ir2 = 2*ir      ! 6:2l-2
          im2 = ir2+2     ! 8:2l



          



          do k=1,n 



             loca = (k-1)*(lfft2+2)
             locb = (k-1)*(lfft2+2)

             ilowr  = locb +ir                  ! 3:l-1
             iupr   = locb +lfft2 -(ir-2) ! 2l-1:l-1
             ilowim = locb +im                  ! 4:l
             iupim  = locb +lfft2 - (ir-3)  ! 2l:l

             a(loca+ir1) = b(ilowr)+b(iupr)
             a(loca+im1) = a(loca+im1-4) - (b(ilowim)-b(iupim))
             a(loca+ir2) = b(ilowim)+b(iupim)
             a(loca+im2) = a(loca+im2-4) + (b(ilowr)-b(iupr))

          end do
       end do

!cccccccccccccc  postprocessing on the cosine coefficients to get 
!ccccccccccccccc Chebyshev coefficients 

       if (is.eq.-1) then 


          do k=1,n
             loca = (k-1)*(lfft2+2)

             a(loca+1)=  a(loca+1)/2.0
             a(loca+2) = a(loca+2)/2.0
             a(loca+lfft2+1) = a(loca+lfft2+1)/2.0
             a(loca+lfft2+2) = a(loca+lfft2+2)/2.0

          end do

       endif
         
      end subroutine 


      subroutine rcs(a,b,work,lfft,n,ifax,trig)
      dimension a(1),b(1),work(1)
      integer ifax(1)
      dimension trig(1) 
!ccccccc      real to symmetric complex fft  
!     Note: lfft must be even. n is the number of ffts to be 
!     performed. The input , a , is a packed array and the 
!     output is returned in packed form to a . 
!       The inverse of rcs is csr. IS = -1 always for rcs.  
!ccccccc    first , unpack the data  cccccccccccccccccccccccccccc

      do k=1,n 
         loca = (k-1)*lfft     ! starting location or a
         locb = (k-1)*(lfft+2) ! starting location for b
         do i=1,lfft
            b(i+locb) = a(i+loca)
         end do
      end do
 
!ccccccc    now perform the real to complex transform on b  cccccc
            is =   -1
            inc =   1
            jump =   lfft +2
            lot =   n 


          call rfftmlt(b,work,trig,ifax,inc,jump,lfft,lot,is)

!ccccccc   now pack the result  ccccccccccccccccccccccccccccccccccc
          
      
          do k =1,n

             loca = (k-1)*lfft       !starting location for a
             locb = (k-1)*(lfft+2)   !starting location for b

             loc1 = locb + (lfft+1)   !location for real part of last mode
             loc2 = locb + 2          !location for imag part of first mode        

             b(loc2) = b(loc1)


             do i=1,lfft

                a(i+loca)  =  b(i+locb) 

             end do
          end do
        
      end subroutine 
