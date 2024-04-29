PROGRAM FIXED_POINT2
    IMPLICIT NONE
    REAL :: G1, Y, G2, X, REGN1(2), REGN2(2), G12, G21, G11, G22
    REAL :: VAL1, VAL2, VAL3, VAL4, X0, Y0, MAXX, X1, Y1, K
    REAL :: NORM, TEMP1, TEMP2, TOL
    INTEGER :: N,I

    X0 = .25
    Y0 = .25
    G11 = 0
    G12 = 1./SQRT(5.)
    N = 2 ! NUMBER OF EQUATIONS
    TOL = 0.0000001
    VAL1 = G11
    VAL2 = G12
    VAL3 = ABS(G21(Y0))
    VAL4 = ABS(G22(Y0))
    MAXX = MAX(VAL1,VAL2,VAL3,VAL4)
    K = MAXX*N
    IF(K>=1)THEN
        WRITE(*,*) 'THE SYSTEM WILL NOT CONVERGE'
    ELSE
        DO I=1,5
            X1 = G1(Y0)
            Y1 = G2(X0, Y0)
            TEMP1 = X0
            TEMP2 = Y0
            X0 = X1
            Y0 = Y1
            NORM = SQRT(((X1 - TEMP1)**2) + ((Y1 - TEMP2)**2))
            write(*,*)X1,Y1
            IF(NORM < TOL)EXIT
        END DO
        WRITE(*,*) 'APPROXIMATIONS:'
        WRITE(*,*) 'X = ',X1
        WRITE(*,*) 'Y = ',Y1
    END IF
END PROGRAM

FUNCTION G1(Y)
    IMPLICIT NONE
    REAL :: G1,Y
    G1 = Y/SQRT(5.)
    RETURN
END FUNCTION

FUNCTION G2(X, Y)
    IMPLICIT NONE
    REAL :: G2,X, Y
    G2 = .25*(SIN(X) + COS(Y))
    RETURN
END FUNCTION

FUNCTION G21(X)
    IMPLICIT NONE
    REAL :: G21, X
    G21 = .25*COS(X)
    RETURN
END FUNCTION

FUNCTION G22(Y)
    IMPLICIT NONE
    REAL :: G22, Y
    G22 = -.25*SIN(Y)
    RETURN
END FUNCTION