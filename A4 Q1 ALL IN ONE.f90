PROGRAM ASSIGNMENT_1
IMPLICIT NONE
INTEGER::N
CALL SHOOTING(N)
PRINT*,' '
CALL FINITE(N)
END PROGRAM



SUBROUTINE SHOOTING (N)
IMPLICIT NONE
    REAL::K11,K12,K21,K22,K31,K32,K41,K42
    REAL::KP11,KP12,KP21,KP22,KP31,KP32,KP41,KP42
    REAL::U1(0:10),U2(0:10),V1(0:10),V2(0:10)
    REAL::H,T,Y,F,A,B,P,Q,R,X,alpha,beta
    INTEGER::I,J,N
    A=1.
    B=2.
    N=10
    alpha=1.
    beta=2.


    H=(B-A)/N
    U1(0)=alpha
    U2(0)=0
    V1(0)=0
    V2(0)=1

    DO I=0,N-1
        X=A+I*H
        K11=H*U2(I)
        K12=H*(P(X)*U2(I)+Q(X)*U1(I)+R(X))
        K21=H*(U2(I)+0.5*K12)
        K22=H*(P(X+H/2)*(U2(I)+0.5*K12)+Q(X+H/2)*(U1(I)+0.5*K11)+R(X+H/2))
        K31=H*(U2(I)+0.5*K22)
        K32=H*(P(X+H/2)*(U2(I)+0.5*K22)+Q(X+H/2)*(U1(I)+0.5*K21)+R(X+H/2))
        K41=H*(U2(I)+K32)
        K42=H*(P(X+H)*(U2(I)+K32)+Q(X+H)*(U1(I)+K31)+R(X+H))
        U1(I+1)=U1(I)+(1/6.)*(K11+2*K21+2*K31+K41)
        U2(I+1)=U2(I)+(1/6.)*(K12+2*K22+2*K32+K42)

        KP11=H*V2(I)
        KP12=H*(P(X)*V2(I)+Q(X)*V1(I))
        KP21=H*(V2(I)+0.5*KP12)
        KP22=H*(P(X+H/2)*(V2(I)+0.5*KP12)+Q(X+H/2)*(V1(I)+0.5*KP11))
        KP31=H*(V2(I)+0.5*KP22)
        KP32=H*(P(X+H/2)*(V2(I)+0.5*KP22)+Q(X+H/2)*(V1(I)+0.5*KP21))
        KP41=H*(V2(I)+KP32)
        KP42=H*(P(X+H)*(V2(I)+KP32)+Q(X+H)*(V1(I)+KP31))
        V1(I+1)=V1(I)+(1/6.)*(KP11+2*KP21+2*KP31+KP41)
        V2(I+1)=V2(I)+(1/6.)*(KP12+2*KP22+2*KP32+KP42)

    END DO
    PRINT*,'          ===================SHOOTING METHOD=============================='
    PRINT*,'          N      t                y(t)             F(X)           ERROR'


    DO I=0,N
        X=A+I*H
        Y=U1(I)+V1(I)*(BETA-U1(N))/V1(N)
        WRITE(*,*)I,X,Y,F(X),ABS(F(X)-Y)
    END DO

    PRINT*,' '
    PRINT*,'THE SOLUTION IS = ',Y

END SUBROUTINE SHOOTING

SUBROUTINE FINITE(N)
    IMPLICIT NONE
    REAL::A1,B1,ALPHA,BETA,P,Q,R,F,H,X
    INTEGER::N,I
    REAL,ALLOCATABLE::A(:),B(:),C(:),D(:),L(:),U(:),Z(:),W(:)
    A1=1.
    B1=2.
    ALPHA=1.
    BETA=2.
    N=9

    H=(B1-A1)/(N+1)

    ALLOCATE(A(0:N+1))
    ALLOCATE(B(0:N+1))
    ALLOCATE(C(0:N+1))
    ALLOCATE(D(0:N+1))
    ALLOCATE(L(0:N+1))
    ALLOCATE(U(0:N+1))
    ALLOCATE(Z(0:N+1))
    ALLOCATE(W(0:N+1))

    X=A1+H
    A(1)=2+H**2*Q(X)
    B(1)=-1+(H/2.)*P(X)
    D(1)=-H**2*R(X)+(1+(H/2)*P(X))*ALPHA


    DO I=2,N-1
        X=A1+I*H
        A(I)=2+H**2*Q(X)
        B(I)=-1+(H/2)*P(X)
        C(I)=-1-(H/2.)*P(X)
        D(I)=-H**2*R(X)
    END DO


    X=B1-H
    A(N)=2+H**2*Q(X)
    C(N)=-1-(H/2)*P(X)
    D(N)=-H**2*R(X)+(1-(H/2)*P(X))*BETA


    L(1)=A(1)
    U(1)=B(1)/A(1)
    Z(1)=D(1)/L(1)


    DO I=2,N-1
        L(I)=A(I)-C(I)*U(I-1)
        U(I)=B(I)/L(I)
        Z(I)=(D(I)-C(I)*Z(I-1))/L(I)
    END DO


    L(N)=A(N)-C(N)*U(N-1)
    Z(N)=(D(N)-C(N)*Z(N-1))/L(N)


    W(0)=ALPHA
    W(N+1)=BETA
    W(N)=Z(N)


    DO I=N-1,1,-1
        W(I)=Z(I)-U(I)*W(I+1)
    END DO
     PRINT*,'         ===============FINITE DIFFRENCE METHOD=========================='
    PRINT*,'          N     t                y(t)             F(x)            ERROR'

    DO I=0,N+1
        X=A1+I*H
        PRINT*,I,X,W(I),F(X),ABS(F(X)-W(I))
    END DO

    PRINT*,' '
    PRINT*,'THE SOLUTION IS = ',W(N+1)
END SUBROUTINE FINITE

FUNCTION P(X)
    IMPLICIT NONE
    REAL::P,X
    P=-2/X
END FUNCTION

FUNCTION Q(X)
    IMPLICIT NONE
    REAL::Q,X
    Q=2/X**2
END FUNCTION

FUNCTION R(X)
    IMPLICIT NONE
    REAL::R,X
    R=SIN(LOG(X))/X**2
END FUNCTION



FUNCTION F(X)
    IMPLICIT NONE
    REAL::F,X,C1,C2
    C1=1.1392070132
    C2=-0.03920701320
    F=C1*X+C2/X**2-(3/10.)*SIN(LOG(X))-(1/10.)*COS(LOG(X))
END FUNCTION
