      SUBROUTINE INHULL(XP,YP,NP,X,Y,N,LIST,LPTR,LEND,INH)
      INTEGER N, LIST(*), LPTR(*), LEND(N)
      DOUBLE PRECISION    XP(*), YP(*), X(*), Y(*)
      LOGICAL INH(*)
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
C   Modified version of BNODES from TRIPACK
C   Albrecht Gebhardt <albrecht.gebhardt@uni-klu.ac.at>
C   
C This function returns TRUE if a given point is contained
C in the convex hull of a given triangulation.
C
C On input:
C
C       PX,PY = coordinates of the point to be located
C
C       X,Y = Arrays containing the coordinates of the nodes
C             in the triangulation.
C
C       N = Number of nodes in the triangulation.  N .GE. 3.
C
C       LIST,LPTR,LEND = Data structure defining the trian-
C                        gulation.  Refer to subroutine
C                        TRMESH.
C
C The above parameters are not altered by this routine.
C
C Return value:
C
C       TRUE if (PX,PY) lies in convex hull of the triangulation.
C
C       This is determined by sequentially applying the LEFT
C       function along the boundary (counterclockwise). 
C       If all these tests return TRUE (PX,PY) lies in the convex hull.
C
C Modules required by INHULL:  LEFT
C
C***********************************************************
C     
      INTEGER K, LP, N0, NST
      INTEGER N1, NK, NKM1, IP
      LOGICAL LEFT

      DO 10, IP=1,NP
C 
C     initialize LOGICAL INH 
C
      INH(IP) = .TRUE.
C     
C     Set NST to the first boundary node encountered.
C     
      NST = 1
 1    LP = LEND(NST)
      IF (LIST(LP) .LT. 0) GO TO 2
      NST = NST + 1
      GO TO 1
C     
C     Initialization.
C     
 2    N1 = NST
      NK = NST
      K = 1
      N0 = NST
C     
C     Traverse the boundary in counterclockwise order.
C     
 3    LP = LEND(N0)
      LP = LPTR(LP)
      N0 = LIST(LP)
      IF (N0 .EQ. NST) GO TO 4
      K = K + 1
      NKM1 = NK 
      NK = N0
C     
C     Check if point (X,Y) is left of the line through NODES(K-1)->NODES(K),
C     
      INH(IP) = INH(IP) .AND. LEFT(X(NKM1),Y(NKM1),
     .                             X(NK),Y(NK),
     .                             XP(IP),YP(IP))
      GO TO 3        
C     
C     Check if point (X,Y) is left of the line through NODES(NB)->NODES(1),
C     
 4    INH(IP) = INH(IP) .AND. LEFT(X(NK),Y(NK),
     .                             X(N1),Y(N1),
     .                             XP(IP),YP(IP))

 10   CONTINUE
      RETURN
      END

      SUBROUTINE ONHULL(XP,YP,NP,X,Y,N,LIST,LPTR,LEND,ONH)
      INTEGER N, LIST(*), LPTR(*), LEND(N)
      DOUBLE PRECISION    XP(*), YP(*), X(*), Y(*)
      LOGICAL ONH(*)
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
C   Modified version of BNODES from TRIPACK
C   Albrecht Gebhardt <albrecht.gebhardt@uni-klu.ac.at>
C   
C This function returns TRUE if a given point lays
C on the convex hull of a given triangulation.
C
C On input:
C
C       PX,PY = coordinates of the point to be located
C
C       X,Y = Arrays containing the coordinates of the nodes
C             in the triangulation.
C
C       N = Number of nodes in the triangulation.  N .GE. 3.
C
C       LIST,LPTR,LEND = Data structure defining the trian-
C                        gulation.  Refer to subroutine
C                        TRMESH.
C
C The above parameters are not altered by this routine.
C
C Return value:
C
C       TRUE if (PX,PY) lies on the convex hull of the triangulation.
C
C       This is determined by sequentially applying the LEFT
C       function along the boundary (counterclockwise, to ensure 
c       (PX,PY) is IN the hull) and applying BORDER to check if it 
c       lays ON some border segment
C
C Modules required by ONHULL:  LEFT, BORDER
C
C***********************************************************
C     
      INTEGER K, LP, N0, NST
      INTEGER N1, NK, NKM1, IP
      LOGICAL LEFT, BORDER

      DO 10, IP=1,NP
C 
C     initialize LOGICAL INH 
C
      ONH(IP) = .TRUE.
C     
C     Set NST to the first boundary node encountered.
C     
      NST = 1
 1    LP = LEND(NST)
      IF (LIST(LP) .LT. 0) GO TO 2
      NST = NST + 1
      GO TO 1
C     
C     Initialization.
C     
 2    N1 = NST
      NK = NST
      K = 1
      N0 = NST
C     
C     Traverse the boundary in counterclockwise order.
C     
 3    LP = LEND(N0)
      LP = LPTR(LP)
      N0 = LIST(LP)
      IF (N0 .EQ. NST) GO TO 4
      K = K + 1
      NKM1 = NK
      NK = N0
C     
C     Check if point (X,Y) is left of the line through NODES(K-1)->NODES(K),
C     
      ONH(IP) = ONH(IP) .AND. LEFT(X(NKM1),Y(NKM1),
     .                             X(NK),Y(NK),
     .                             XP(IP),YP(IP))
      GO TO 3        
C     
C     Check if point (X,Y) is left of the line through NODES(NB)->NODES(1),
C     
 4    ONH(IP) = ONH(IP) .AND. LEFT(X(NK),Y(NK),
     .                             X(N1),Y(N1),
     .                             XP(IP),YP(IP))

      if(ONH(IP)) then
c     IP is IN the hull, now check if it is ON the hull:
         ONH(IP)=.FALSE.
         NST = 1
 11      LP = LEND(NST)
         IF (LIST(LP) .LT. 0) GO TO 12
         NST = NST + 1
         GO TO 11
C     
C     Initialization.
C     
 12      N1 = NST
         NK = NST
         K = 1
         N0 = NST
C     
C     Traverse the boundary in counterclockwise order.
C     
 13      LP = LEND(N0)
         LP = LPTR(LP)
         N0 = LIST(LP)
         IF (N0 .EQ. NST) GO TO 14
         K = K + 1
         NKM1 = NK
         NK = N0
C     
C     Check if point (X,Y) is on the line through NODES(K-1)->NODES(K),
C     
         ONH(IP) = ONH(IP) .OR. BORDER(X(NKM1),Y(NKM1),
     .                                 X(NK),Y(NK),
     .                                 XP(IP),YP(IP))
         GO TO 13        
C     
C     Check if point (X,Y) is on the line through NODES(NB)->NODES(1),
C     
 14      ONH(IP) = ONH(IP) .OR. BORDER(X(NK),Y(NK),
     .                                 X(N1),Y(N1),
     .                                 XP(IP),YP(IP))
         
      endif
      
 10   CONTINUE
      RETURN
      END


      LOGICAL FUNCTION BORDER (X1,Y1,X2,Y2,X0,Y0)
      DOUBLE PRECISION    X1, Y1, X2, Y2, X0, Y0
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
C modified version of LEFT: A. Gebhardt
C
C   This function determines whether node N0 is on the line 
C through N1-N2 
C 
C
C
C On input:
C
C       X1,Y1 = Coordinates of N1.
C
C       X2,Y2 = Coordinates of N2.
C
C       X0,Y0 = Coordinates of N0.
C
C Input parameters are not altered by this function.
C
C On output:
C
C       BORDER = .TRUE. if and only if (X0,Y0) is on the
C              directed line N1->N2.
C
C Modules required by BORDER:  None
C
C***********************************************************
C
      DOUBLE PRECISION DX1, DY1, DX2, DY2
C
C Local parameters:
C
C DX1,DY1 = X,Y components of the vector N1->N2
C DX2,DY2 = X,Y components of the vector N1->N0
C
      DX1 = X2-X1
      DY1 = Y2-Y1
      DX2 = X0-X1
      DY2 = Y0-Y1
C
C
      BORDER = DX1*DY2 .eq. DX2*DY1
      RETURN
      end
