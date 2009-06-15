      SUBROUTINE VLEFT (X1,Y1,X2,Y2,X0,Y0,N0,LEFT0)
      IMPLICIT NONE
      DOUBLE PRECISION X1, Y1, X2, Y2, X0(*), Y0(*)
      INTEGER N0
      LOGICAL LEFT0(*)
c     modified version of LEFT, accepts X0,Y0 as vektor, 
c     returns vector
C
C***********************************************************
C
C                                               From TRIPACK
C                                            Robert J. Renka
C                                  Dept. of Computer Science
C                                       Univ. of North Texas
C                                           renka@cs.unt.edu
C                                                   09/01/88
C
C   This function determines whether node N0 is to the left
C or to the right of the line through N1-N2 as viewed by an
C observer at N1 facing N2.
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
C       LEFT = .TRUE. if and only if (X0,Y0) is on or to the
C              left of the directed line N1->N2.
C
C Modules required by LEFT:  None
C
C***********************************************************
C
      INTEGER I
C
C Local parameters:
C
C DX1,DY1 = X,Y components of the vector N1->N2
C DX2,DY2 = X,Y components of the vector N1->N0
      DO 10 I=1,N0
c         DX1 = X2-X1
c         DY1 = Y2-Y1
c         DX2 = X0-X1
c         DY2 = Y0-Y1
C
C If the sign of the vector cross product of N1->N2 and
C   N1->N0 is positive, then sin(A) > 0, where A is the
C   angle between the vectors, and thus A is in the range
C   (0,180) degrees.
C
         LEFT0(I) = (X2-X1)*(Y0(I)-Y1) .GE. 
     +        (X0(I)-X1)*(Y2-Y1)
 10   CONTINUE
      RETURN
      END
