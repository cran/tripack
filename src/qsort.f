C
C     from an older version of TRIPACK, I need it it in
C     voronoi.f
C     A. Gebhardt <albrecht.gebhardt@uni-klu.ac.at>
C

      SUBROUTINE QSORT (N,X, IND)
      INTEGER N, IND(N)
      DOUBLE PRECISION    X(N)
C
C***********************************************************
C
C                                               From TRIPACK
C                                            Robert J. Renka
C                                  Dept. of Computer Science
C                                       Univ. of North Texas
C                                             (817) 565-2767
C                                                   09/01/88
C
C   This subroutine uses an order N*LOG(N) quick sort to
C sort the real array X into increasing order.  The algor-
C ithm is as follows.  IND is initialized to the ordered
C sequence of indexes 1,...,N, and all interchanges are
C applied to IND.  X is divided into two portions by picking
C a central element T.  The first and last elements are com-
C pared with T, and interchanges are applied as necessary so
C that the three values are in ascending order.  Inter-
C changes are then applied so that all elements greater than
C T are in the upper portion of the array and all elements
C less than T are in the lower portion.  The upper and lower
C indices of one of the portions are saved in local arrays,
C and the process is repeated recursively on the other
C portion.  When a portion is completely sorted, the process
C begins again by retrieving the indexes bounding another
C unsorted portion.
C
C
C On input:
C
C       N = Number of elements to be sorted.
C
C       X = Array of length N to be sorted.
C
C The above parameters are not altered by this routine.
C
C       IND = Array of length at least N.
C
C On output:
C
C       IND = Sequence of integers 1,...,N permuted in the
C             the same fashion that X would be.  Thus, the
C             sorted array may be stored in an array Y with
C             the assignment statements:  Y(I) = X(IND(I))
C             for I = 1 to N.  Alternatively, X may be over-
C             written with the sorted array by a call to
C             subroutine PERMUT.
C
C Modules required by QSORT:  None
C
C Intrinsic functions called by QSORT:  DBLE, INT
C
C***********************************************************
C
C NOTE -- IU and IL must be dimensioned .GE. log(N), where
C         the log has base 2.
C
C***********************************************************
C
      INTEGER IU(21), IL(21)
      INTEGER M, I, J, K, L, IJ, IT, ITT, INDX
      DOUBLE PRECISION    R, T
C
C Local parameters:
C
C IU,IL =  Temporary storage for the upper and lower
C            indexes of portions of the array X
C M =      Index for IU and IL
C I,J =    Lower and upper indexes of a portion of X
C K,L =    Indexes in the range I,...,J
C IJ =     Randomly chosen index between I and J
C IT,ITT = Temporary storage for interchanges in IND
C INDX =   Temporary index for X
C R =      Pseudo random number for generating IJ
C T =      Central element of X
C
      IF (N .LE. 0) RETURN
C
C Initialize IND, M, I, J, and R.
C
      DO 1 I = 1,N
        IND(I) = I
    1   CONTINUE
      M = 1
      I = 1
      J = N
      R = .375
C
C Top of loop --
C
    2 IF (I .GE. J) GO TO 7
      IF (R .GT. .5898437) THEN
        R = R - .21875
      ELSE
        R = R + .0390625
      ENDIF
C
C Initialize K.
C
    3 K = I
C
C Select a central element of X and save it in T.
C
      IJ = I + INT(R*DBLE(J-I))
      IT = IND(IJ)
      T = X(IT)
C
C If the first element of the array is greater than T,
C   interchange it with T.
C
      INDX = IND(I)
      IF (X(INDX) .GT. T) THEN
        IND(IJ) = INDX
        IND(I) = IT
        IT = INDX
        T = X(IT)
      ENDIF
C
C Initialize L.
C
      L = J
C
C If the last element of the array is less than T,
C   interchange it with T.
C
      INDX = IND(J)
      IF (X(INDX) .GE. T) GO TO 5
      IND(IJ) = INDX
      IND(J) = IT
      IT = INDX
      T = X(IT)
C
C If the first element of the array is greater than T,
C   interchange it with T.
C
      INDX = IND(I)
      IF (X(INDX) .LE. T) GO TO 5
      IND(IJ) = INDX
      IND(I) = IT
      IT = INDX
      T = X(IT)
      GO TO 5
C
C Interchange elements K and L.
C
    4 ITT = IND(L)
      IND(L) = IND(K)
      IND(K) = ITT
C
C Find an element in the upper part of the array which is
C   not larger than T.
C
    5 L = L - 1
        INDX = IND(L)
        IF (X(INDX) .GT. T) GO TO 5
C
C Find an element in the lower part of the array which is
C   not smaller than T.
C
    6 K = K + 1
        INDX = IND(K)
        IF (X(INDX) .LT. T) GO TO 6
C
C If K .LE. L, interchange elements K and L.
C
      IF (K .LE. L) GO TO 4
C
C Save the upper and lower subscripts of the portion of the
C   array yet to be sorted.
C
      IF (L-I .GT. J-K) THEN
        IL(M) = I
        IU(M) = L
        I = K
        M = M + 1
      ELSE
        IL(M) = K
        IU(M) = J
        J = L
        M = M + 1
      ENDIF
      GO TO 8
C
C Begin again on another unsorted portion of the array.
C
    7 M = M - 1
      IF (M .EQ. 0) RETURN
      I = IL(M)
      J = IU(M)
C
    8 IF (J-I .GE. 11) GO TO 3
      IF (I .EQ. 1) GO TO 2
      I = I - 1
C
C Sort elements I+1,...,J.  Note that 1 .LE. I .LT. J and
C   J-I .LT. 11.
C
    9 I = I + 1
      IF (I .EQ. J) GO TO 7
      INDX = IND(I+1)
      T = X(INDX)
      IT = INDX
      INDX = IND(I)
      IF (X(INDX) .LE. T) GO TO 9
      K = I
C
   10 IND(K+1) = IND(K)
        K = K - 1
        INDX = IND(K)
        IF (T .LT. X(INDX)) GO TO 10
      IND(K+1) = IT
      GO TO 9
      END
