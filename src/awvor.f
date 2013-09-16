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
      DOUBLE PRECISION X1, X2, X3, Y1, Y2, Y3, XC(2), YC(2), SA, AR, 
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

