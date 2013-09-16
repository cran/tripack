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
C        and contains them all)
C        only the "inner" circle is the solution needed here:
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
C TODO:
C     Formulate criteria to distingusih the three cases!!!
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
     $     XX(3), RW(25),A(3),B(3),C(3),P(3,2), Q(3,2),
     $     XA(6,6,4),YA(6,6,4),XS(2),YS(2),RS(2)
      DOUBLE PRECISION FEQN,MYATAN,TRAFOX,TRAFOY,DASINH,DACOSH,
     $     ITRFOX,ITRFOY
      INTEGER DSGN, DIR(3), HI,HJ
      EXTERNAL FEQN,MYATAN,TRAFOX,TRAFOY,DASINH,DACOSH,
     $     ITRFOX,ITRFOY,DSGN
      
      DOUBLE PRECISION LX1,LX2,LX3,LY1,LY2,LY3,LW1,LW2,LW3,
     $     EPS,PHI(3),MX(3),MY(3),F1,F2,F3,T(3),XH,YH,XT,YT,TH1,TH2,
     $     MTX,MTY
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

C     the midpoint of the triangle
      MTX=(X(1)+X(2)+X(3))/3.0D0
      MTY=(Y(1)+Y(2)+Y(3))/3.0D0

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
C     IF (.NOT. RATIO) RETURN ??????

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
c         C(I)=DSQRT(U(I)*U(I)+V(I)*V(I)) ??? why = SQRT(DS(I))/2

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

         P(I,1)=DTAN(MYATAN(A(I),B(I))+PHI(I))
         Q(I,1)=MY(I)-P(I,1)*(MX(I))
         P(I,2)=DTAN(MYATAN(-A(I),B(I))+PHI(I))
         Q(I,2)=MY(I)-P(I,2)*(MX(I))
c        write(*,*)"abline(",Q(I,1),",",P(I,1),")"
c        write(*,*)"abline(",Q(I,2),",",P(I,2),")"
         
 112  CONTINUE

      DO 113 I=1,3
         
C     get intersections of hyperbola asymptotes as better starting solution,
c     

C TODO:
         DO 114 J=I,3
            DO 1141 K=1,4
               IF (J.NE.I) THEN
                  IF(K.EQ.1) THEN
                     HI=1;HJ=1
                  ELSE IF (K.EQ.2) THEN
                     HI=2;HJ=1
                  ELSE IF (K.EQ.2) THEN
                     HI=1;HJ=2
                  ELSE IF (K.EQ.2) THEN
                     HI=2;HJ=2
                  ENDIF
                  
c     asymptote HI of hyperbola I intersects asymptote HJ of hyperbola J:
c     
c     P(I)*X+Q(I)=P(J)*X+Q(J)
c     ->    X=(Q(J)-Q(I))/(P(I)-P(J))
c     ->    Y=P(I)*X+Q(I)
                  IF(DABS(P(I,HI)-P(J,HJ)).GE.EPS) THEN
                     XA(I,J,K)=(Q(J,HJ)-Q(I,HI))/(P(I,HI)-P(J,HJ))
                     YA(I,J,K)=P(I,HI)*XA(I,J,K)+Q(I,HI)
c     keep only intersections with x'>=0
c     and x' is xa in hyperbola centric coordinates:
c     and dir=1 (that is right hyperbola, else dir=-1)

                     IF(DIR(I)*TRAFOX(XA(I,J,K),YA(I,J,K),
     $                    -PHI(I),-MX(I),-MY(I)).GE.0.0D0 .AND.
     $                    DIR(J)*TRAFOX(XA(I,J,K),YA(I,J,K),
     $                    -PHI(J),-MX(J),-MY(J)).GE.0.0D0) THEN
                        IA(I,J,K)=1
                        
c          write(*,*)"points(",xa(I,J,K),", ",ya(I,J,K),",col='red')" 
                     ELSE
                        IA(I,J,K)=-1
                     END IF
                     
                  ELSE
                     IA(I,J,K)=0
                  END IF
               END IF
 1141       CONTINUE
 114     CONTINUE
 113  CONTINUE
      
C     ?????????????????
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
        DO 116 J=I,3
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

c     write(*,*)"text(",xa(I,J,k),", ",ya(I,J,k),",",LRCNT,")"

c     Each hyperbola intersection is generated by three hyperbolas,
c     hence it is "surrounded" by three intersections of asymptotes.
c     take the average of these three points as starting solution
c     (a single asymptotes intersection can be numerically bad if
c     asymptotes are near to parallel), NS(i) will be between 0 and 3
c     depending on the situation if an  asymptotes intersection has been 
c     recorded above or not ...
C     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!????????????
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
                END IF
 117         CONTINUE
      
c      write(*,*)"# ",ITRFOX(XA(I,J,K),YA(I,J,K),-PHI(I),-MX(I),-MY(I)),
c     $ ", ",ITRFOY(XA(I,J,K),YA(I,J,K),-PHI(I),-MX(I),-MY(I))
          END IF
 116   CONTINUE
 115  CONTINUE
      IF(NS(1).GT.0)THEN
         XS(1)=XS(1)/NS(1)
         YS(1)=YS(1)/NS(1)
         RS(1)=DSQRT((XS(2)-MTX)**2+(YS(2)-MTY)**2)
c        write(*,*)"# XS(1),YS(1),RS(1),NS(1):",XS(1),YS(1),RS(1),NS(1)
      END IF
      IF(NS(2).GT.0)THEN
         XS(2)=XS(2)/NS(2)
         YS(2)=YS(2)/NS(2)
         RS(2)=DSQRT((XS(2)-MTX)**2+(YS(2)-MTY)**2)
c        write(*,*)"# XS(2),YS(2),RS(2),NS(2):",XS(2),YS(2),RS(2),NS(2)
      END IF
      

      
C     
C     Call DSOS from SLATEC to solve the 3 hyperbola intersection equations
C     
c     start with ordinary circum circle center
      XX(1)=XC(1);XX(2)=YC(1);XX(3)=RC(1); 
c     print? set to -1
      IW(1)=-1
c     max iter:
      IW(2)=100
      IFLAG=-1
C     
c      write(*,*)"# starting values: X:",XX(1)," ,Y:",XX(2),
c     $     " ,R:",XX(3)
c      CALL DSOS (FEQN, 3, XX, EPS, EPS, EPS, IFLAG, RW, 25,
c     $     IW, 6)
c      IF (IFLAG.GT.3 .AND. IFLAG.LT.8) THEN
c         WRITE(*,*)"# DSOS returned correctable error flag: ",IFLAG
c         WRITE(*,*)"# ... retry with better starting point"
C     TODO...
c     try the intersections of the hyperbolas with the triangle edges
c         DO 344 ITRY=1,3
c            IF(ITRY.EQ.1) THEN
c               I1=2;I2=1
c            ELSE IF (ITRY.EQ.1) THEN
c               I1=3; I2=2
c            ELSE IF (ITRY.EQ.1) THEN
c               I1=1; I2=3
c            END IF
c            XX(1)=MX(ITRY)
c            XX(2)=MY(ITRY)
c     not very usefull??:
c            XX(3)=DABS((SQRT((X(I1)-X(I2))**2+(Y(I1)-Y(I2))**2)-
c     $           W(I1)-W(I2))/2.0D0); 
c     print? set to -1 
c            IW(1)=-1
c     max iter:
c            IW(2)=100
c            IFLAG=-1
c            write(*,*)"# starting values: X:",XX(1)," ,Y:",XX(2),
c     $     " ,R:",XX(3)
c
c            CALL DSOS (FEQN, 3, XX, EPS, EPS, EPS, IFLAG, RW
c     $           , 25, IW, 6)
c            IF (IFLAG.GT.3 .AND. IFLAG.LT.8) THEN
c               WRITE(*,*)"# DSOS again returned correctable error flag:",
c     $              IFLAG
c               CONTINUE
c            ELSE IF (IFLAG.GT.7) THEN
c               CS=-1
c               WRITE(*,*)"# DSOS returned FINAL error flag: ",IFLAG
c               RETURN
c            ELSE
C     solution:
c               WRITE(*,*)"# DSOS retry succeeded in step: ",ITRY
c               GOTO 345
c            END IF
c            CS=-1
c 344     CONTINUE
C     try with XS1,YS1 (intersection of asymptotes) from above
       
         XX(1)= XS(1)
         XX(2)= YS(1)
         XX(3)= RS(1)
c        write(*,*)"# starting values: X:",XX(1)," ,Y:",XX(2),
c    $     " ,R:",XX(3)

         CALL DSOS (FEQN, 3, XX,EPS,EPS,EPS, IFLAG, RW,25,
     $        IW, 6)
c     solution:
        IF (IFLAG.LE.3) THEN
           GOTO 345
        END IF
C     No solution
c      ELSE IF (IFLAG.GT.7) THEN
c         CS=-1
c         WRITE(*,*)"# DSOS returned error flag: ",IFLAG
c         RETURN
c      ELSE
c         GOTO 345
c      END IF
      
C     solution found:
 345  CONTINUE
      CS=1; XC(1)=XX(1); YC(1)=XX(2); RC(1)=XX(3)
C     get parameter t for hyperbola intersection:
c     write(*,*)"# first intersection found!!"
c     WRITE(*,*)"# solution XC:",XC(1)
c     WRITE(*,*)"# solution YC:",YC(1)
c     WRITE(*,*)"# solution RC:",RC(1)

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
C TODO: RC(1) !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
         XX(1)=XS(2) 
         XX(2)=YS(2)
         XX(3)= RS(2)

ccc   RC(1);
c     print? set to -1
         IW(1)=-1
c     max iter:
         IW(2)=100
         IFLAG=-1
c        write(*,*)"# starting values: X:",XX(1)," ,Y:",XX(2),
c    $     " ,R:",XX(3)

         CALL DSOS(FEQN, 3, XX, EPS,EPS,EPS, IFLAG, RW, 25,
     $        IW, 6)
         IF (IFLAG.GT.3 .AND. IFLAG.LT.8) THEN
c           write(*,*)"# error: second intersection not found"
c           write(*,*)"# starting values where: X:",XS(2)," ,Y:",YS(2),
c    $          " R:", RS(2)
            CONTINUE
         ELSE
c           write(*,*)"# second intersection found!!"
            CS=2; XC(2)=XX(1); YC(2)=XX(2); RC(2)=XX(3)
c            write(*,*)"circles(",XX(1),",",XX(2),",",XX(3),")"
c           WRITE(*,*)"# solution XC:",XC(2)
c           WRITE(*,*)"# solution YC:",YC(2)
c           WRITE(*,*)"# solution RC:",RC(2)
         END IF
      END IF
      
c     W1=T(1);W2=T(2);W3=T(3)


C     IF (.NOT. RATIO) RETURN
      RETURN
      END

      DOUBLE PRECISION FUNCTION FEQN(PAR, K)
C     Hyperbola equations for one triangle with weigths
C     K indicates which edge of triangle is of interest.
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
