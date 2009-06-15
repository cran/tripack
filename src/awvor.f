      SUBROUTINE AWVORONOI(NCC,LCC,N,X,Y,WT,LIST,LPTR,LEND,NT,LCCC,ICCC,
     $                     LCT,LTRI,MAXIT,IER)
      IMPLICIT NONE
      INTEGER NCC, LCC(*), N, LIST(*), LPTR(*), LEND(N),
     $        NT, LCT(NCC), IER, LTRI(9,*), ICCC(8,*),MAXIT
      DOUBLE PRECISION X(*), Y(*), LCCC(14,*), WT(*)
C
C     DETERMINE AN ADDITIVE WEIGHTED  VORONOI MOSAIC. 
C       TO BE USED WITH THE FUNCTIONS FROM TRIPACK.
C
C     THIS IS DONE BY DETERMINING THE WEIGHTED CIRCUMCIRCLE CENTERS OF A 
C     DELAUNAY TRIANGULATION, USING ITS TRIANGLE NEIGBOURSHIP RELATIONS AS 
C     (initial) ADJACENCY RELATION FOR THE CIRCUMCIRCLE CENTERS.
C
C     THE TRIANGULATION PRODUCED WITH TRLIST FROM TRIPACK CAN 
C     INCLUDE EMPTY TRIANGLES WITH AREA 0 (IF COLLINEAR POINTS EXIST).
C     THESE TRIANGLES MAKE TROUBLES IF THEY OCCUR ON THE BORDER OF THE 
C     TRIANGULATION. SO WE MUST REMOVE THEM.
C     
C     A.GEBHARDT <albrecht.gebhardt@uni-klu.ac.at>
C
C     ON ENTRY:
C     ... SEE E.G. TRLIST
C
C     ON EXIT:
C     NT   - NUMBER OF TRIANGLES
C     LCCC - LIST OF CIRCUMCIRCLE CENTERS (13 by NT):
C            ROW 1,2     : X, Y COORDINATES
C            ROW 3       : TRIANGLE AREA (-1 = REMOVED TRIANGLE)
C            ROW 4       : ASPECT RATIO
C            ROW 5       : CIRCUMCIRCLE RADIUS
C            ROW 6,7,8   : HYPERBOLA PARAMETER t
C            ROW 9,10    : X, Y COORDINATES OF SECOND CCC
C            ROW 11      : SECOND CC RADUIS
C            ROW 12,13,14: HYPERBOLA PARAMETER t FOR SECOND CCC

C     ICCC   LIST OF NEIGHBOUR RELATIONS (6 by NT):
C            ROW 1-3: NEIGHBOUR TRIANGLE INDICES
C            ROW 4-6: NODE INDICES OF TRIANGLE
C            ROW 7: NO OF SOLUTIONS FOUND, 0, 1, 2
C            ROW 8: EDGE INDEX WHERE WCCC LEFT THE INTERIOR OF
C                   THIS TRIANGLE (needed for plotting on the border)
C            TODO: how to handle second ccc??
C
C     Modules required by AWVORONOI:  TRLIST, CIRCUM, QSORT, RMSHNB,
C     TRSWP

C     DETERMINE ALL CIRCUMCIRCLE CENTERS
      INTEGER NROW, IT, IP(3), NTRI, IT1, IT2, IW(3), CS,I1,I2,
     $     INT,J,IN1,IN2,IO1,IO2,IS(3),JS(3),I,LP21,ITCNT
      DOUBLE PRECISION X1, X2, X3, Y1, Y2, Y3, XC, YC, SA, AR, 
     $    XP(3), YP(3), DEPS, RC, W(3),XA,YA,XB,YB,SCAL
      LOGICAL RATIO
      
      RATIO = .TRUE.
      DEPS = 1.0D-7
      NROW = 9
      ITCNT=0

c 1    continue

      CALL TRLIST(NCC,LCC,N,LIST,LPTR,LEND,NROW,NT,LTRI,LCT,IER)
      NTRI=NT

      DO 10 IT = 1,NTRI
         CALL WCRCTR(IT,LTRI(1,IT),ICCC(1,IT),LCCC(1,IT),X,Y,WT)
 10   CONTINUE

C     REMOVE TRIANGLES WITH AREA 0 FROM THE BORDER, UPDATE NEIGBOUR RELATIONS.
C     AREA=0 MEANS COLLINEAR POINTS, 3 CASES (I'M NOT SURE IF THEY ALL WILL
C     OCCUR!):
C       ALL POINTS DIFFERENT
C       2   POINTS COINCIDE
C       ALL POINTS COINCIDE
C
 11   CONTINUE
      DO 20 IT = 1,NTRI
         IF (((LCCC(3,IT).LE.DEPS) .OR.(LCCC(4,IT).LE.DEPS)) .AND. 
     $        (ICCC(1,IT).EQ.0 .OR. 
     $         ICCC(2,IT).EQ.0 .OR. 
     $         ICCC(3,IT).EQ.0)) THEN
C     WE HAVE A 0-TRIANGLE AT THE BORDER
            XP(1)=X(ICCC(4,IT))
            XP(2)=X(ICCC(5,IT))
            XP(3)=X(ICCC(6,IT))
            YP(1)=Y(ICCC(4,IT))
            YP(2)=Y(ICCC(5,IT))
            YP(3)=Y(ICCC(6,IT))
            CALL QSORT(3,XP,IP)
C     CHECK IF ALL X AND/OR Y COORDINATES COINCIDE
            IF (XP(IP(1)) .EQ. XP(IP(3))) THEN
               CALL QSORT(3,XP,IP)
               IF (YP(IP(1)) .EQ. YP(IP(3))) THEN
C     ALL POINTS COINCIDE, ALL NEIGHBOUR TRIANGLES HAVE AREA 0
C     AND WILL BE REMOVED LATER.
C     NOTHING TO DO (?)
                  CONTINUE
               ELSE
C     X1=X2=X3, USE Y-COORDINATES 
                  IF ((YP(IP(1)).EQ.YP(IP(2))) 
     $                 .OR. (YP(IP(2)).EQ.YP(IP(3)))) THEN
C     TWO POINTS COINCIDE
                     IF (YP(IP(1)).EQ.YP(IP(2))) THEN
                        IT1=ICCC(IP(1),IT)
                        IT2=ICCC(IP(2),IT)
                        CALL RMSHNB(IT,IT1,IT2,LCCC,ICCC)
                     ELSE
                        IT1=ICCC(IP(2),IT)
                        IT2=ICCC(IP(3),IT)
                        CALL RMSHNB(IT,IT1,IT2,LCCC,ICCC)
                     ENDIF
                  ELSE
C     ALL POINTS DIFFERENT
                     IT1=ICCC(IP(1),IT)
                     IT2=ICCC(IP(2),IT)
                     CALL RMSHNB(IT,IT1,IT2,LCCC,ICCC)
                     IT1=ICCC(IP(3),IT)
                     IT2=ICCC(IP(2),IT)
                     CALL RMSHNB(IT,IT1,IT2,LCCC,ICCC)
                  ENDIF
               ENDIF
            ELSE
C     USE X COORDINATES
               IF ((XP(IP(1)).EQ.XP(IP(2))) 
     $              .OR. (XP(IP(2)).EQ.XP(IP(3)))) THEN
C     TWO POINTS COINCIDE
                  IF (XP(IP(1)).EQ.XP(IP(2))) THEN
                     IT1=ICCC(IP(1),IT)
                     IT2=ICCC(IP(2),IT)
                     CALL RMSHNB(IT,IT1,IT2,LCCC,ICCC)
                  ELSE
                     IT1=ICCC(IP(2),IT)
                     IT2=ICCC(IP(3),IT)
                     CALL RMSHNB(IT,IT1,IT2,LCCC,ICCC)
                  ENDIF
               ELSE
C     ALL POINTS DIFFERENT
                  IT1=ICCC(IP(1),IT)
                  IT2=ICCC(IP(2),IT)
                  CALL RMSHNB(IT,IT1,IT2,LCCC,ICCC)
                  IT1=ICCC(IP(2),IT)
                  IT2=ICCC(IP(3),IT)
                  CALL RMSHNB(IT,IT1,IT2,LCCC,ICCC)
               ENDIF
            ENDIF
         ENDIF         
 20   CONTINUE
C     REPEAT THE ABOVE STEPS UNTIL NO MORE AREA 0 TRIANGLES 
C     ARE LEFT ON THE BORDER.
      DO 30 IT = 1,NTRI
         IF ((LCCC(3,IT).LE.DEPS) .AND. 
     $       (ICCC(1,IT).EQ.0 .OR. 
     $        ICCC(2,IT).EQ.0 .OR. 
     $        ICCC(3,IT).EQ.0)) THEN
            GO TO 11
         ENDIF
 30   CONTINUE

C     SWAP TRIANGLES
C
C          +-----+C     +-----+C
C         /\    /      /   ,'/
C        / 1\2 /  ->  / 2-1 /
C       /    \/      /,'   /
C     A+-----+     A+-----+ 
C
C     IF WEIGTED CIRCUMCIRCLE CENTERS 1 AND 2 "BYPASS" EACH OTHER
C     (THAT IS THE ANGLE (1>2,A>B) LEAVES THE INTERVAL [-PI/2,PI/2])
C

 51   ITCNT=ITCNT+1

      DO 40 IT = 1,NTRI

         IF(LCCC(IT,3).GE.0.0D0) THEN
            DO 50 INT = 1,3
c     FILL "I1" "I2" WITH NODE INDICES OF EDGE COMMON WITH 
c     POSSIBLE NEIGHBOUR "INT":
               IF(INT.EQ.1) THEN
                  I1=2;I2=3
               ENDIF
               IF(INT.EQ.2) THEN
                  I1=3;I2=1
               ENDIF
               IF(INT.EQ.3) THEN
                  I1=1;I2=2
               ENDIF
               
C     ITERATE OVER EXISTING NEIGBOUR TRIANGLES:
               IF(ICCC(INT,IT)>0)THEN
c      write(*,*)"iccc(int,it):",ICCC(INT,IT)
C     SKIP REMOVED TRIANGLES:
                  IF(LCCC(3, ICCC(INT,IT)).GE.0.0D0) THEN
                     
C     GET CIRCUMCIRCLE CENTERS OF BOTH TRIANGLES "IT" AND 
C     NEIGHBOUR "ICCC(INT,IT)":
                     X1=LCCC(1,IT)
                     Y1=LCCC(2,IT)
                     X2=LCCC(1,ICCC(INT,IT))
                     Y2=LCCC(2,ICCC(INT,IT))
c      write(*,*)"x1,y1:",x1,y1               
c      write(*,*)"x2,y2:",x2,y2               
C     GET OPPOSITE NODES A AND C:
C     INDEX OF A = INDEX OF NEIGHBOUR "INT"
c      write(*,*)"it:",it," int:",int               
c      write(*,*)"iccc(3+int,it):",    ICCC(3+INT,IT)           
                     XA=X(ICCC(3+INT,IT))
                     YA=Y(ICCC(3+INT,IT))
c      write(*,*)"xa,ya:",xa,ya            
C     INDEX OF B:
C     SEARCH INDEX J IN  NEIGHBOUR TRIANGLE WHICH IS DIFFERENT 
C     FROM INDEX ON SHARED EDGE (IT IS THE OPPOSITE NODE B)  
                     DO 60 J=1,3
                        IF(ICCC(3+J,ICCC(INT,IT)).NE.I1 .AND.            .!!!..   
     $                       ICCC(3+J,ICCC(INT,IT)).NE.I2) THEN
                           XB=X(ICCC(3+J,ICCC(INT,IT)))
                           YB=Y(ICCC(3+J,ICCC(INT,IT)))
                        ENDIF
 60                  CONTINUE
C     SWAP IF
C     ANGLE BETWEEN VECTOR A>B and 1>2 OUTSIDE [-PI/2,PI/2]
c                     write(*,*)"scalar product:",
c     $                    (XB-XA)*(X2-X1)+(YB-YA)*(Y2-Y1)
                     SCAL=(XB-XA)*(X2-X1)+(YB-YA)*(Y2-Y1)
c                     write(*,*)"SCAL:",SCAL
                     IF(SCAL .LT. 0.0D0) THEN
                        I1=IT
                        I2=ICCC(INT,IT)
                        CALL TRSWP(I1,I2,X,Y,WT,NTRI,ICCC,LCCC,
     $                       NCC,LCC,N,LIST,LPTR,LEND,NROW,NT,LTRI,
     $                       LCT,IER )
C     RESTART (hopefully converges!, IER=1: non-convex swap skipped):
                        IF(IER.NE.1)THEN
                           IF(ITCNT.LT.MAXIT)THEN
                              GO TO 51
c                           restart:
c                           ELSE
c                              GO TO 52
                           END IF
                        END IF
                     ENDIF
                  ENDIF
               ENDIF
 50         CONTINUE
         ENDIF
 40   CONTINUE
C     SUCCESS:
      IER=0
      RETURN
C     DID not CONVERGE:
      IER=1
 52   CONTINUE



      END


      SUBROUTINE TRSWP(I1,I2,X,Y,WT,NTRI,ICCC,LCCC,
     $     NCC,LCC,N,LIST,LPTR,LEND,NROW,NT,LTRI,
     $     LCT,IER )

      IMPLICIT NONE
      INTEGER I1,I2,NCC, LCC(*), N, LIST(*), LPTR(*), LEND(N),
     $     NT, LCT(NCC), IER, LTRI(9,*), ICCC(8,*),NROW,NTRI
      DOUBLE PRECISION X(*), Y(*), LCCC(14,*), WT(*)
      
      INTEGER IN1,IN2,IO1,IO2,IS(3),JS(3),I,J,LP21,IT,ISW
      LOGICAL CONVEX
      LOGICAL LEFT
      EXTERNAL LEFT


      IER=0
C     sort indices, really necessary ????
      IF(I1.GT.I2)THEN
         ISW=I2
         I2=I1
         I1=ISW
      END IF

c      WRITE(*,*)"SWAP: ",I1," <-> ",I2

      IS(1)=ICCC(4,I1)
      IS(2)=ICCC(5,I1)
      IS(3)=ICCC(6,I1)
      JS(1)=ICCC(4,I2)
      JS(2)=ICCC(5,I2)
      JS(3)=ICCC(6,I2)
      IO1=-1
      IO2=-1
      IN1=-1
      IN2=-1

C     find common edge IO1 - IO2
      DO 45 I=1,3
         DO 46 J=1,3
            IF(IS(I).EQ.JS(J))THEN
               IF(IO1.EQ.-1)THEN
                  IO1=IS(I)
               ELSE
                  IO2=IS(I)
               ENDIF
            ENDIF
 46      CONTINUE
 45   CONTINUE

c     find opposite points IN1 and IN2

      DO 47 I=1,3
         IF((IS(I).NE.IO1).AND.(IS(I).NE.IO2))THEN
            IN1=IS(I)
         ENDIF
         IF((JS(I).NE.IO1).AND.(JS(I).NE.IO2))THEN
            IN2=JS(I)
         ENDIF

 47   CONTINUE
c      WRITE(*,*)"IO1,IO2,IN1,IN2:",IO1,IO2,IN1,IN2
      

c     sort IO1 -> IN2 -> IO2 -> IN1 counterclockwise 
c     triangle IO1-IO2-IN1 is counterclockwise ordered if:
c     X(IN1),Y(IN1) is LEFT from line X(IO1),Y(IO1) -- X(IO2),Y(IO2) 
c     triangle IO2-IO1-IN2 is counterclockwise ordered if:
c     X(IN2),Y(IN2) is LEFT from line X(IO2),Y(IO2) -- X(IO1),Y(IO1)
      IF(.NOT.LEFT(X(IO1),Y(IO1),X(IO2),Y(IO2), X(IN1),Y(IN1)) .AND.
     $ .NOT.LEFT(X(IO2),Y(IO2),X(IO1),Y(IO1), X(IN2),Y(IN2))) THEN
c     swap IO1 and IO2:
         ISW=IO1
         IO1=IO2
         IO2=ISW
      END IF

c     ensure that io1 -> io2 -> in1 -> in2 forms two counterclowise
c     enumerated triangles in the following way:
c
c           in1<-_   
c            |    -_
c     io1 ---+-------> io2
c            |
c            v
c           in2
c     swap in2 <-> in1 if in1 not left of io1->io2:
      IF(.NOT.LEFT(X(IO1),Y(IO1),X(IN2),Y(IN2), X(IO2),Y(IO2)))THEN
         ISW=IN1
         IN1=IN2
         IN2=ISW         
      END IF

C     now test for convexity of
c     quadrilateral IO1 -> IN2 -> IO2 -> IN1

      CONVEX=LEFT(X(IO1),Y(IO1),X(IN2),Y(IN2), X(IO2),Y(IO2)) .AND. 
     $       LEFT(X(IO1),Y(IO1),X(IN2),Y(IN2), X(IN1),Y(IN1)) .AND. 
     $       LEFT(X(IN2),Y(IN2),X(IO2),Y(IO2), X(IN1),Y(IN1)) .AND. 
     $       LEFT(X(IN2),Y(IN2),X(IO2),Y(IO2), X(IO1),Y(IO1)) .AND. 
     $       LEFT(X(IO2),Y(IO2),X(IN1),Y(IN1), X(IO1),Y(IO1)) .AND. 
     $       LEFT(X(IO2),Y(IO2),X(IN1),Y(IN1), X(IN2),Y(IN2)) .AND. 
     $       LEFT(X(IN1),Y(IN1),X(IO1),Y(IO1), X(IN2),Y(IN2)) .AND. 
     $       LEFT(X(IN1),Y(IN1),X(IO1),Y(IO1), X(IO2),Y(IO2))

c
      IF(CONVEX)THEN
         CALL SWAP(IN1,IN2,IO1,IO2, LIST,LPTR,
     $        LEND, LP21)
c     WRITE(*,*)"LP21:",LP21
      ELSE
c         WRITE(*,*)"SWAP CANCELLED, NOT CONVEX: ",I1," <-> ",I2
         IER=1
         RETURN
      END IF
c     goto 1


C     GET NEW TRIANGLE LISTS
      CALL TRLIST(NCC,LCC,N,LIST,LPTR,LEND,NROW,NT,LTRI,
     $     LCT,IER)

C     UPDATE CIRCUMCIRCLE CENTER DATA (+ ICCC, LCCC)
      DO 100 IT = 1,NTRI
         CALL WCRCTR(IT,LTRI(1,IT),ICCC(1,IT),LCCC(1,IT),X,Y,WT)
 100  CONTINUE

      RETURN
      END

      SUBROUTINE WCIRCUM (X,Y,W,
     $     RATIO, XC,YC,RC,TC,SA,AR,CS)

      IMPLICIT NONE
      LOGICAL RATIO
      DOUBLE PRECISION X(3), Y(3), XC(2), YC(2), RC(2),
     $           SA, AR, W(3),TC(6)
      INTEGER CS
C
C***********************************************************
C
C     modified version from CIRCUM from
C                                               From TRIPACK
C                                            Robert J. Renka
C                                  Dept. of Computer Science
C                                       Univ. of North Texas
C                                           renka@cs.unt.edu
C                                                   12/10/96
C
C     A. Gebhardt, agebhard@uni-klu.ac.at, Dec 2006
C
C     Given three vertices defining a triangle and three weights 
c     defining three circles centered in the triangle nodes, this 
C     subroutine returns 0, 1 or 2 circumcenters, circumradiii, 
C     and hyperbola parameters (as modified "mittelsenkrechte")
C
C     There are three possible situations:
C     
C     * a standard triangle and weigthing circles that result in an 
C       outer and an inner circumcircle 
C       (in the sense that the inner circumcircle touches the three
C        weight circles but does not contain them and the outer touches
C        and contains them)
C        only the "inner" circle is solution:
C
C                    .   _    
C               .      /   \ °
C                     (  +  )   
C           °          \ _ /     °
C            _      .  °  .        
C         °/   \  .         °      .
C         (  +  ).           °   _   
C        . \ _ /             . /   \. 
C                °           .(  +  ) 
C         .       .         .  \ _ /. 
C                   ° . . °        
C           °                   °
C               .           .
C                  °  .  ° 
C
C     * a triangle (not far from colinear case) and weigthing circles that 
C     result in two "inner" circumcircles 
C     (the middle weigting circle is contained in the convex hull of the 
C     two outer weigting circles):
C     both "inner" circles are solutions:
C
C            .                              .                    
C                                             _ _        
C             °                           ° ´      `                
C            _ _.                      . /           \   
C          ´      `  .          ._  °                                 
C       /           \    °  °  /   \    |      +      |   
C                             (  +  )                            
C      |      +      |         \ _ /     \           /          
C                         .  °    °    .    . _ _  ,              
C       \           /  °                   °                          
C          . _ _  ,.                          .
C                                              
C                °                              °
C              .               
C
C     * a triangle (not far from colinear case) and weigthing circles that 
C     result in two "outer" circumcircles 
C     (the convex hull of the smaller outer weigthing circles divides 
C     the large middle weigthing circle into two non connected parts)
C     none of these circles are solutions:
C
C                                       
C                                    
C                                                            
C                           .  _°_ °    °   .          
C                    .  °    ´      `           °  _.          
C                .        /           \          /   \ °     .
C     .       °_                                (  +  )   . 
C          ° /   \       |      +      |         \ _ / .    
C        .  (  +  )                              .  °       °
C           .\ _ /        \           /     .          
C      °        °    .       . _ _  , .                
C                         °    °                                 
C                                                      
C                              
C
c TODO:signed triangle area, and, optionally, the aspect ratio of the
C triangle.
C
C
C On input:
C
C       X(3), Y(3) = Cartesian coordinates of the vertices.
C       W(3) = weights (=radii)
C
C       RATIO = Logical variable with value TRUE if and only
C               if the aspect ratio is to be computed.
C
C Input parameters are not altered by this routine.
C
C On output:
C
C       XC,YC = Cartesian coordinates of the weigthed circumcenter(s)
C               (center of the circle defined by the three
C               points) unless SA = 0, in which XC and YC
C               are not altered.
C
C       RC = Circumradius (radius of the circle defined by
C            the three circles) unless SA = 0 (infinite
C            radius), in which case CR is not altered.
C            If < 0, circle does not exist.
C
C       SA = Signed triangle area with positive value if
C            and only if the vertices are specified in
C            counterclockwise order:  (X3,Y3) is strictly
C            to the left of the directed line from (X1,Y1)
C            toward (X2,Y2).
C
C       AR = Aspect ratio r/CR, where r is the radius of the
C            inscribed circle, unless RATIO = FALSE, in
C            which case AR is not altered.  AR is in the
C            range 0 to .5, with value 0 iff SA = 0 and
C            value .5 iff the vertices define an equilateral
C            triangle.
C FIXME:
C     W1,W2,W3: overwritten with T(1), T(2), T(3)
C     
C     CS = no of solutions found, 0,1,2
C
C Modules required by WCIRCUM:  None
C
C Intrinsic functions called by WCIRCUM:  ABS, SQRT
C
C***********************************************************
C
      INTEGER I, IW(6),IFLAG, I1, I2, ITRY, J, IA(6,6,4),K,
     $     LRCNT,NS(2)
      DOUBLE PRECISION DS(3), FX, FY, U(3), V(3), 
     $     XX(3), RW(25),A(3),B(3),P(3,2), Q(3,2),
     $     XA(6,6,4),YA(6,6,4),XS(2),YS(2)
      DOUBLE PRECISION FEQN,MYATAN,TRAFOX,TRAFOY,DASINH,DACOSH,
     $     ITRFOX,ITRFOY
      INTEGER DSGN, DIR(3)
      EXTERNAL FEQN,MYATAN,TRAFOX,TRAFOY,DASINH,DACOSH,
     $     ITRFOX,ITRFOY,DSGN
      
      DOUBLE PRECISION LX1,LX2,LX3,LY1,LY2,LY3,LW1,LW2,LW3,
     $     EPS,PHI(3),MX(3),MY(3),F1,F2,F3,T(3),XH,YH,XT,YT,TH1,TH2
      LOGICAL LEFT
      EXTERNAL LEFT

      DATA EPS/1.0D-8/
      COMMON /PARAMS/LX1,LX2,LX3,LY1,LY2,LY3,LW1,LW2,LW3

      LX1=X(1);LX2=X(2);LX3=X(3);LY1=Y(1);LY2=Y(2);LY3=Y(3);
      LW1=W(1);LW2=W(2);LW3=W(3)
c     W(1)=W1;W(2)=W2;W(3)=W3

C     
C     Set U(K) and V(K) to the x and y components, respectively,
C     of the directed edge opposite vertex K.
C     
      U(1) = X(3) - X(2)
      U(2) = X(1) - X(3)
      U(3) = X(2) - X(1)
      V(1) = Y(3) - Y(2)
      V(2) = Y(1) - Y(3)
      V(3) = Y(2) - Y(1)

C     the midpoints of the edge opposite vertex K.
      MX(1) = X(2) + U(1)/2.
      MX(2) = X(3) + U(2)/2.
      MX(3) = X(1) + U(3)/2.
      MY(1) = Y(2) + V(1)/2.
      MY(2) = Y(3) + V(2)/2.
      MY(3) = Y(1) + V(3)/2.

C     the slope of the edge opposite vertex K.
      PHI(1)=MYATAN(U(1),V(1))
      PHI(2)=MYATAN(U(2),V(2))
      PHI(3)=MYATAN(U(3),V(3))

c     the direction of the edge opposite vertex K
c     (from larger weight to smaller weight, seen in counterclockwise order 
c     around the triangle)
c     DIR > 0 means we use right hyperbola, otherwise left
      DIR(3)=DSGN(W(1)-W(2))
      DIR(1)=DSGN(W(2)-W(3))
      DIR(2)=DSGN(W(3)-W(1))
      
C     
C     Set SA to the signed triangle area.
C     
      SA = (U(1)*V(2) - U(2)*V(1))/2.
      IF (SA .EQ. 0.) THEN
         IF (RATIO) AR = 0.
         RETURN
      ENDIF
C     
C     Set DS(K) to the squared distance from the origin to
C     vertex K.
C     
      DS(1) = X(1)*X(1) + Y(1)*Y(1)
      DS(2) = X(2)*X(2) + Y(2)*Y(2)
      DS(3) = X(3)*X(3) + Y(3)*Y(3)


      
C     
C     Compute XC and YC as ordinary circumcircele center, take this as 
C     (possibly bad) starting point for iteration later.
C     

      FX = 0.
      FY = 0.
      DO 1 I = 1,3
         FX = FX - DS(I)*V(I)
         FY = FY + DS(I)*U(I)
 1    CONTINUE
      XC(1) = FX/(4.*SA)
      YC(1) = FY/(4.*SA)
      RC(1) = SQRT( (XC(1)-X(1))**2 + (YC(1)-Y(1))**2 )
      IF (.NOT. RATIO) RETURN

C     
C     Compute the squared edge lengths and aspect ratio.
C     
      DO 2 I = 1,3
         DS(I) = U(I)*U(I) + V(I)*V(I)
 2    CONTINUE
      AR = 2.*ABS(SA)/
     $     ( (SQRT(DS(1)) + SQRT(DS(2)) + SQRT(DS(3)))*RC(1) )


      DO 112 I=1,3
C     GET HYPERBOLA PARAMTERS A AND B
         IF(I.EQ.1) A(I)=(W(3)-W(2))/2.0D0
         IF(I.EQ.2) A(I)=(W(1)-W(3))/2.0D0
         IF(I.EQ.3) A(I)=(W(2)-W(1))/2.0D0

         B(I)=SQRT(DS(I)/4.0D0 - A(I)*A(I))

c     GET ASYMPTOTES FOR HYPERBOLA I
C     in hyperbola coordinates x',y':
c     
c     y'= (+/-) a/b *x'
c     slope: alpha'=atan(a/b)
c     slope in world coordinates:
c     alpha=alpha'+phi
c     asymptote in world coordnates:
c     y=p*x+q
c     p=tan(alpha)  
c     find point (dx,dy) on asymptote:
c     dy = p * dx +q
c     q  = dy - p*dx  !!!!
C     get PI, QI, for both asymptotes, dx=MX, dy=MY

         P(I,1)=DTAN(MYATAN(A(I),B(I)) +PHI(I))
         Q(I,1)=MY(I)-P(I,1)*(MX(I))
         P(I,2)=DTAN(MYATAN(-A(I),B(I)) +PHI(I))
         Q(I,2)=MY(I)-P(I,2)*(MX(I))
c         write(*,*)"abline(",Q(I,1),",",P(I,1),")"
c         write(*,*)"abline(",Q(I,2),",",P(I,2),")"
         
 112  CONTINUE

      DO 113 I=1,3
         
C     get intersections of hyperbola asymptotes as better starting solution,
c     

C TODO:
C         DO 114 J=I,3
         DO 114 J=I,3
            IF (J.NE.I)THEN

c     asymptote of HYP I intersects asymptote j:
c     
c     P(I)*X+Q(I)=P(J)*X+Q(J)
c     ->    X=(Q(J)-Q(I))/(P(I)-P(J))
c     ->    Y=P(I)*X+Q(I)
               IF(DABS(P(I,1)-P(J,1)).GE.EPS) THEN
                  XA(I,J,1)=(Q(J,1)-Q(I,1))/(P(I,1)-P(J,1))
                  YA(I,J,1)=P(I,1)*XA(I,J,1)+Q(I,1)
c     keep only intersections with x'>=0
c     and x' is xa in hyperbola centric coordinates:
c     and dir=1 (that is right hyperola, else dir=-1)

                  IF(DIR(I)*TRAFOX(XA(I,J,1),YA(I,J,1),
     $                 -PHI(I),-MX(I),-MY(I)).GE.0.0D0 .AND.
     $                 DIR(J)*TRAFOX(XA(I,J,1),YA(I,J,1),
     $                 -PHI(J),-MX(J),-MY(J)).GE.0.0D0) THEN
                     IA(I,J,1)=1
c           write(*,*)"points(",xa(I,J,1),", ",ya(I,J,1),")" 
                  ELSE
                     IA(I,J,1)=-1
                  END IF

               ELSE
                  IA(I,J,1)=0
               END IF

               IF(DABS(P(I,2)-P(J,1)).GE.EPS) THEN
                  XA(I,J,2)=(Q(J,1)-Q(I,2))/(P(I,2)-P(J,1))
                  YA(I,J,2)=P(I,2)*XA(I,J,2)+Q(I,2)

c     keep only intersections with x'>=0
c     and x' is xa in hyperbola centric coordinates:
                  IF(DIR(I)*TRAFOX(XA(I,J,2),YA(I,J,2),
     $                 -PHI(I),-MX(I),-MY(I)).GE.0.0D0 .AND.
     $                 DIR(J)*TRAFOX(XA(I,J,2),YA(I,J,2),
     $                 -PHI(J),-MX(J),-MY(J)).GE.0.0D0) THEN
                     IA(I,J,2)=1
c           write(*,*)"points(",xa(I,J,2),", ",ya(I,J,2),")" 
                  ELSE
                     IA(I,J,2)=-1
                  END IF
               ELSE
                  IA(I,J,2)=0
               END IF

               IF(DABS(P(I,1)-P(J,2)).GE.EPS) THEN
                  XA(I,J,3)=(Q(J,2)-Q(I,1))/(P(I,1)-P(J,2))
                  YA(I,J,3)=P(I,1)*XA(I,J,3)+Q(I,1)

c     keep only intersections with x'>=0
c     and x' is xa in hyperbola centric coordinates:
                  IF(DIR(I)*TRAFOX(XA(I,J,3),YA(I,J,3),
     $                 -PHI(I),-MX(I),-MY(I)).GE.0.0D0 .AND.
     $                 DIR(J)*TRAFOX(XA(I,J,3),YA(I,J,3),
     $                 -PHI(J),-MX(J),-MY(J)).GE.0.0D0) THEN
                     IA(I,J,3)=1
c           write(*,*)"points(",xa(I,J,3),", ",ya(I,J,3),")" 
                  ELSE
                     IA(I,J,3)=-1
                  END IF
                
               ELSE
                  IA(I,J,3)=0
               END IF

               IF(DABS(P(I,2)-P(J,2)).GE.EPS) THEN
                  XA(I,J,4)=(Q(J,2)-Q(I,2))/(P(I,2)-P(J,2))
                  YA(I,J,4)=P(I,2)*XA(I,J,4)+Q(I,2)

c     keep only intersections with x'>=0
c     and x' is xa in hyperbola centric coordinates:
                  IF(DIR(I)*TRAFOX(XA(I,J,4),YA(I,J,4),
     $                 -PHI(I),-MX(I),-MY(I)).GE.0.0D0 .AND.
     $                 DIR(J)*TRAFOX(XA(I,J,4),YA(I,J,4),
     $                 -PHI(J),-MX(J),-MY(J)).GE.0.0D0) THEN
                     IA(I,J,4)=1
c          write(*,*)"points(",xa(I,J,4),", ",ya(I,J,4),")"
                  ELSE
                     IA(I,J,4)=-1
                  END IF
               ELSE
                  IA(I,J,4)=0
               END IF
            END IF
 114     CONTINUE
 113  CONTINUE

C     Merge every three hyperbola asymptote intersections into
C     a starting solution, distinguish the three cases mentioned above
C     (0, 1 or 2 solutions)
C     
C     two cases: 
C     * circumcircle center opposite to vertex A wrt the 
C       other edge BC: intersections are 2times left 1 time right to 
C       traingle edges
C       or: left to all edges
C     * circumcircle center on same side as A wrt to edge BC:
C       asympt. intersections (and solution) 2times right 1 time left 
C       to  triangle edges
C     * FIXME: what if no solution exists??
C
      XS(1)=0.0D0;XS(2)=0.0D0;YS(1)=0.0D0;YS(2)=0.0D0;NS(1)=0;NS(2)=0
      DO 115 I=1,3
C     TODO use symmetry
         DO 116 J=1,3
            IF (J.NE.I)THEN
               DO 117 K=1,4
                  IF (IA(I,J,K).EQ.1) THEN
                     LRCNT=0
                   IF(LEFT(X(1),Y(1),X(2),Y(2),XA(I,J,K),YA(I,J,K)))THEN
                     LRCNT=LRCNT+1
                   ELSE
                     LRCNT=LRCNT-1
                   END IF

                   IF(LEFT(X(2),Y(2),X(3),Y(3),XA(I,J,K),YA(I,J,K)))THEN
                     LRCNT=LRCNT+1
                   ELSE
                     LRCNT=LRCNT-1
                   END IF

                   IF(LEFT(X(3),Y(3),X(1),Y(1),XA(I,J,K),YA(I,J,K)))THEN
                     LRCNT=LRCNT+1
                   ELSE
                     LRCNT=LRCNT-1
                   END IF

c      write(*,*)"text(",xa(I,J,k),", ",ya(I,J,k),",",LRCNT,")"
      IF(LRCNT.EQ.-1) THEN
         XS(2)=XS(2)+XA(I,J,K)
         YS(2)=YS(2)+YA(I,J,K)
         NS(2)=NS(2)+1
      END IF
      IF(LRCNT.EQ.1 .OR. LRCNT.EQ.3 ) THEN
         XS(1)=XS(1)+XA(I,J,K)
         YS(1)=YS(1)+YA(I,J,K)
         NS(1)=NS(1)+1
      END IF 
      IF(NS(1).GT.0)THEN
         XS(1)=XS(1)/NS(1)
         YS(1)=YS(1)/NS(1)
c         write(*,*)"# XS1,YS1:",XS(1),YS(1)
      END IF
      IF(NS(2).GT.0)THEN
         XS(2)=XS(2)/NS(2)
         YS(2)=YS(2)/NS(2)
c         write(*,*)"# XS2,YS2:",XS(2),YS(2)
      END IF
      
c      write(*,*)"# ",ITRFOX(XA(I,J,K),YA(I,J,K),-PHI(I),-MX(I),-MY(I)),
c     $ ", ",ITRFOY(XA(I,J,K),YA(I,J,K),-PHI(I),-MX(I),-MY(I))
                  END IF
 117           CONTINUE
            END IF
 116     CONTINUE
 115  CONTINUE


      
C     
C     Call DSOS from SLATEC to solve the 3 hyperbola intersection equation
C     
c     start with ordinary circum circle center
      XX(1)=XC(1);XX(2)=YC(1);XX(3)=RC(1); 
c     print? set to -1
      IW(1)=0
c     max iter:
      IW(2)=100
      IFLAG=-1
C     
      CALL DSOS (FEQN, 3, XX, 1.0D-9, 1.0D-16, 1.0D-16, IFLAG, RW, 25,
     $     IW, 6)
      IF (IFLAG.GT.3 .AND. IFLAG.LT.8) THEN
c         WRITE(*,*)"DSOS returned correctable error flag: ",IFLAG
c         WRITE(*,*)"... retry with better starting point"
C     TODO...
c     try the intersections of the hyperbolas with the triangle edges
         
         DO 344 ITRY=1,3
            IF(ITRY.EQ.1) THEN
               I1=2;I2=1
            ELSE IF (ITRY.EQ.1) THEN
               I1=3; I2=2
            ELSE IF (ITRY.EQ.1) THEN
               I1=1; I2=3
            END IF
            XX(1)=MX(ITRY)
            XX(2)=MY(ITRY)
c     not very usefull??:
            XX(3)=DABS((SQRT((X(I1)-X(I2))**2+(Y(I1)-Y(I2))**2)-
     $           W(I1)-W(I2))/2.0D0); 
c     print? set to -1 
            IW(1)=0
c     max iter:
            IW(2)=100
            IFLAG=-1
            CALL DSOS (FEQN, 3, XX, 1.0D-9, 1.0D-16, 1.0D-16, IFLAG, RW
     $           , 25, IW, 6)
            IF (IFLAG.GT.3 .AND. IFLAG.LT.8) THEN
c               WRITE(*,*)"#DSOS again returned correctable error flag:",
c     $              IFLAG
               CONTINUE
            ELSE IF (IFLAG.GT.7) THEN
               CS=-1
c               WRITE(*,*)"# DSOS returned error flag: ",IFLAG
               RETURN
            ELSE
C     solution:
c               WRITE(*,*)"# DSOS retry succeeded in step: ",ITRY
               GOTO 345
            END IF
            CS=-1
 344     CONTINUE
C     try with XS1,YS1 from above
         XX(1)=XS(1); XX(2)=YS(1); XX(3)=RC(1);
        CALL DSOS (FEQN, 3, XX, 1.0D-9, 1.0D-16, 1.0D-16, IFLAG, RW, 25,
     $        IW, 6)
c     solution:
        IF (IFLAG.LE.3) THEN
           GOTO 345
        END IF
C     No solution
      ELSE IF (IFLAG.GT.7) THEN
         CS=-1
c         WRITE(*,*)"# DSOS returned error flag: ",IFLAG
         RETURN
      ELSE
         GOTO 345
      END IF
      
C     solution found:
 345  CONTINUE
      CS=1; XC(1)=XX(1); YC(1)=XX(2); RC(1)=XX(3)
C     get parameter t for hyperbola intersection:

c     WRITE(*,*)"XC:",XC(1)
c     WRITE(*,*)"YC:",YC(1)

      DO 111 I=1,3
C     GET HYPERBOLA PARAMTERS A AND B
         IF(I.EQ.1) A(I)=(W(3)-W(2))/2.0D0
         IF(I.EQ.2) A(I)=(W(1)-W(3))/2.0D0
         IF(I.EQ.3) A(I)=(W(2)-W(1))/2.0D0

         B(I)=SQRT(DS(I)/4.0D0 - A(I)*A(I))


         TH1=0.0D0; TH2=0.0D0
C     TRANSFORM INTERSECTION INTO HYPERBOLA CENTRIC COORDINATES:
         XT=TRAFOX(XC(1),YC(1),-PHI(I),-MX(I),-MY(I))
         YT=TRAFOY(XC(1),YC(1),-PHI(I),-MX(I),-MY(I))
C     USE TRANSORMED Y TO GET T1:
C     y=b*sinh(t)
         IF(ABS(B(I)).GE.EPS) THEN
            TH1=DASINH(YT/B(I))
         ENDIF
c     WRITE(*,*)"T1:",TH1
C     USE TRANSFORMED X TO GET T2:
C     x=a*cosh(t) or -a*cosh(t)
         IF(ABS(A(I)).GE.EPS .AND.
     $        DABS(XT).GE.DABS(A(I))) THEN
            TH2=DACOSH(DABS(XT)/DABS(A(I)))*DSGN(YT)
         ENDIF

c     WRITE(*,*)"T2:",TH2
C     average if TH1 and TH2 have same sign:
         IF(TH1.NE.0.0D0 .AND. TH2.NE.0.0D0 
     $        .AND. TH1*TH2.GT.0.0D0) THEN 
            TC(I)=(TH1+TH2)/2.0D0
c     WRITE(*,*)"average T1,T2"
         ELSE IF(TH1.NE.0.0D0) THEN
            TC(I)=TH1
c     WRITE(*,*)"use T1"
         ELSE IF (TH2.NE.0.0D0) THEN
            TC(I)=TH2
c     WRITE(*,*)"use T2"
         ELSE
            TC(I)=0.0D0
c     WRITE(*,*)"no T?"
         END IF
c     WRITE(*,*)"set T:",TC(I)

c     T(I)=TH1
c     T(I)=TH2
         
         
c     WRITE(*,*)"MX:",MX(I)," MY:",MY(I)," Phi:",PHI(I)
c     WRITE(*,*)"XCt:",XT," YCt:",YT
c     WRITE(*,*)"XCi:",-A*DCOSH(TC(I))," YCi:",B*SINH(TC(I))

         XH=ITRFOX(-A(I)*DCOSH(TC(I)),B(I)*SINH(TC(I)),-PHI(I),
     $        -MX(I),-MY(I))
         YH=ITRFOY(-A(I)*DCOSH(TC(I)),B(I)*SINH(TC(I)),-PHI(I),
     $        -MX(I),-MY(I))
c     WRITE(*,*)"inv X:",XH
c     WRITE(*,*)"inv Y:",YH
c     WRITE(*,*)"check X:",XH-XC(1)
c     WRITE(*,*)"check Y:",YH-YC(1)

c     IF(I.EQ.1) T(I)=T(I)*DSGN(W3-W2)
c     IF(I.EQ.2) T(I)=T(I)*DSGN(W1-W3)
c     IF(I.EQ.3) T(I)=T(I)*DSGN(W2-W1)

 111  CONTINUE

C     second intersection?
      IF(NS(2).GT.0)THEN
C     start with XS1,YS1 from above
         XX(1)=XS(2); XX(2)=YS(2); XX(3)=RC(1);
c     print? set to -1
         IW(1)=0
c     max iter:
         IW(2)=100
         IFLAG=-1
C     
         CALL DSOS(FEQN, 3, XX, 1.0D-9, 1.0D-16, 1.0D-16, IFLAG, RW, 25,
     $        IW, 6)
         IF (IFLAG.GT.3 .AND. IFLAG.LT.8) THEN
c            write(*,*)"error"
         ELSE
            write(*,*)"# second intersection found!!"
            CS=2; XC(2)=XX(1); YC(2)=XX(2); RC(2)=XX(3)
            write(*,*)"circles(",XX(1),",",XX(2),",",XX(3),")"
         END IF
      END IF
      
c     W1=T(1);W2=T(2);W3=T(3)


C     IF (.NOT. RATIO) RETURN
      RETURN
      END

      DOUBLE PRECISION FUNCTION FEQN(PAR, K)

      IMPLICIT NONE
      INTEGER K
      DOUBLE PRECISION PAR(3),X1,X2,X3,Y1,Y2,Y3,W1,W2,W3

      COMMON /PARAMS/X1,X2,X3,Y1,Y2,Y3,W1,W2,W3

      IF(K.EQ.1) THEN
         FEQN=SQRT((X1-PAR(1))**2+(Y1-PAR(2))**2)-W1-PAR(3)
      ELSE IF (K.EQ.2) THEN
         FEQN=SQRT((X2-PAR(1))**2+(Y2-PAR(2))**2)-W2-PAR(3)
      ELSE IF (K.EQ.3) THEN
         FEQN=SQRT((X3-PAR(1))**2+(Y3-PAR(2))**2)-W3-PAR(3)
      END IF

      RETURN
      END

      DOUBLE PRECISION FUNCTION TRAFOX(X,Y,PHI,DX,DY)
C     rotate and shift X, Y by DX, DY, PHI,
C     return X component
      IMPLICIT NONE
      DOUBLE PRECISION X,Y,PHI,DX,DY,MYATAN
      EXTERNAL MYATAN
      TRAFOX=SQRT((X+DX)**2+(Y+DY)**2)*COS(MYATAN((X+DX),(Y+DY))+PHI)
      RETURN
      END

      DOUBLE PRECISION FUNCTION TRAFOY(X,Y,PHI,DX,DY)
C     rotate and shift X, Y by DX, DY, PHI,
C     return Y component
      IMPLICIT NONE
      DOUBLE PRECISION X,Y,PHI,DX,DY,MYATAN
      EXTERNAL MYATAN
      TRAFOY=SQRT((X+DX)**2+(Y+DY)**2)*SIN(MYATAN((X+DX),(Y+DY))+PHI)
      RETURN
      END

      DOUBLE PRECISION FUNCTION ITRFOX(X,Y,PHI,DX,DY)
C     invert rotate and shift X, Y by DX, DY, PHI,
C     return X component
      IMPLICIT NONE
      DOUBLE PRECISION X,Y,PHI,DX,DY,MYATAN
      EXTERNAL MYATAN
      ITRFOX=SQRT(X**2+Y**2)*COS(MYATAN(X,Y)-PHI)-DX
      RETURN
      END

      DOUBLE PRECISION FUNCTION ITRFOY(X,Y,PHI,DX,DY)
C     invert rotate and shift X, Y by DX, DY, PHI,
C     return Y component
      IMPLICIT NONE
      DOUBLE PRECISION X,Y,PHI,DX,DY,MYATAN
      EXTERNAL MYATAN
      ITRFOY=SQRT(X**2+Y**2)*SIN(MYATAN(X,Y)-PHI)-DY
      RETURN
      END

      DOUBLE PRECISION FUNCTION MYATAN(X,Y)
      IMPLICIT NONE
      DOUBLE PRECISION X,Y,EPS,PI
      PARAMETER (PI=3.1415926535D0)
      DATA EPS/1.0D-8/
      
      IF (ABS(X).LE.EPS) THEN
         IF(Y.GE.0) THEN 
            MYATAN=PI/2.
         ELSE
            MYATAN=-PI/2.
         END IF
      ELSE
         IF (X.GE.0) THEN
            MYATAN = ATAN(Y/X)
         ELSE
            MYATAN = ATAN(Y/X)+PI
         END IF
      END IF

      RETURN
      END

      INTEGER FUNCTION DSGN(A)
      IMPLICIT NONE
      DOUBLE PRECISION A
      IF (A.LT.0.0D0) DSGN=-1
      IF (A.EQ.0.0D0) DSGN=0
      IF (A.GT.0.0D0) DSGN=1
      RETURN
      END


      SUBROUTINE WCRCTR(IT,LTRI,ICCC,LCCC,XD,YD,WD)

C     CALL WCIRCUM FOR A GIVEN TRIANGLE IT
      
      IMPLICIT NONE
      INTEGER IT
      INTEGER LTRI(9),ICCC(8),CS
      DOUBLE PRECISION LCCC(14),XD(*),YD(*),WD(*),TC(6)
      DOUBLE PRECISION X(3),Y(3),W(3),XC(2),YC(2),RC(2),SA,AR
      LOGICAL RATIO, LEFT
      EXTERNAL LEFT

      DATA TC/0.0D0,0.0D0,0.0D0,0.0D0,0.0D0,0.0D0/
      XC(1)=0.0D0
      XC(2)=0.0D0
      YC(1)=0.0D0
      YC(2)=0.0D0
      RC(1)=0.0D0
      RC(1)=0.0D0




      W(1) = WD(LTRI(1))
      W(2) = WD(LTRI(2))
      W(3) = WD(LTRI(3))
      X(1) = XD(LTRI(1))
      Y(1) = YD(LTRI(1))
      X(2) = XD(LTRI(2))
      Y(2) = YD(LTRI(2))
      X(3) = XD(LTRI(3))
      Y(3) = YD(LTRI(3))


      RATIO=.TRUE.
      CALL WCIRCUM(X,Y,W,
     $     RATIO,XC,YC,RC,TC,SA,AR,CS)

      IF(CS.GT.0) THEN
         LCCC(1) = XC(1)
         LCCC(2) = YC(1)
         LCCC(5) = RC(1)
         LCCC(6) = TC(1)
         LCCC(7) = TC(2)
         LCCC(8) = TC(3)

         ICCC(7) = CS

         LCCC(9)  = XC(2)
         LCCC(10) = YC(2)
         LCCC(11) = RC(2)
         LCCC(12) = TC(4)
         LCCC(13) = TC(5)
         LCCC(14) = TC(6)

C     CHECK IF wccc IS OUTSIDE OF TRIANGLE, AND IF THIS IS THE CASE
C     RECORD THE INDEX OF THE EDGE WHERE IT LEFT IT
         IF(.NOT.LEFT(X(1),Y(1),X(2),Y(2), XC,YC)) ICCC(8)=LTRI(3)
         IF(.NOT.LEFT(X(2),Y(2),X(3),Y(3), XC,YC)) ICCC(8)=LTRI(1)
         IF(.NOT.LEFT(X(3),Y(3),X(1),Y(1), XC,YC)) ICCC(8)=LTRI(2)

      ELSE
         ICCC(7) = -1
      END IF
      LCCC(3) = ABS(SA)
      LCCC(4) = AR
C     ???:
      ICCC(1) = LTRI(4)
      ICCC(2) = LTRI(5)
      ICCC(3) = LTRI(6)
      ICCC(4) = LTRI(1)
      ICCC(5) = LTRI(2)
      ICCC(6) = LTRI(3)

      RETURN
      END


      DOUBLE PRECISION FUNCTION DPL(X0,Y0,A,B,C)

      IMPLICIT NONE

      DOUBLE PRECISION X0, Y0, A, B, C
C     Distance from point x0,y0 to line ax+by+c=0

      DPL=DABS(A*X0+B*Y0+C)/SQRT(A*A+B*B)
      RETURN
      END


      DOUBLE PRECISION FUNCTION DPL2(X0,Y0,X1,Y1,PHI)

      IMPLICIT NONE

      DOUBLE PRECISION X0, Y0, X1, Y1, PHI, A, B, C
C     Distance from point x0,y0 to line thru x1,y1 with slope phi

      A=DTAN(PHI)
      B=-1.0D0
      C=Y1-A*X1

      DPL2=DABS(A*X0+B*Y0+C)/SQRT(A*A+B*B)
      RETURN
      END      
