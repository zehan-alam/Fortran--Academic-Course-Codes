PROGRAM FIXED_POINT
    IMPLICIT NONE
    REAL :: PH1, PH2, PH3, PH11, PH21, PH31
    REAL :: X1, X0, TEMP, TOL
    
    X0 = 2.0
    TOL = .0001

    IF(ABS(PH11(X0)) < 1.0)THEN
        WRITE(*,*) 'THE FORM (X**3 - 5)/2 CONVERGES'
        DO
            X1 = PH1(X0)
            TEMP = X0
            X0 = X1
            IF(ABS(TEMP - X1) < TOL)EXIT
        END DO
        WRITE(*,*) 'APPROXIMATE ROOT :'
        WRITE(*,*) X1
    END IF

    IF(ABS(PH21(X0)) < 1.0)THEN
        WRITE(*,*) 'THE FORM (2X + 5)**(1/3) CONVERGES'
        DO
            X1 = PH2(X0)
            TEMP = X0
            X0 = X1
            IF(ABS(TEMP - X1) < TOL)EXIT
        END DO
        WRITE(*,*) 'APPROXIMATE ROOT :'
        WRITE(*,*) X1
    END IF

    IF(ABS(PH31(X0)) < 1.0)THEN
        WRITE(*,*) 'THE FORM (2X + 5)/X**2 CONVERGES'
        DO
            X1 = PH3(X0)
            TEMP = X0
            X0 = X1
            IF(ABS(TEMP - X1) < TOL)EXIT
        END DO
        WRITE(*,*) 'APPROXIMATE ROOT :'
        WRITE(*,*) X1
    END IF
END PROGRAM

FUNCTION PH1(X)
    IMPLICIT NONE
    REAL :: PH1, X
    PH1 = ((X**3) - 5)/2
    RETURN
END FUNCTION

FUNCTION PH2(X)
    IMPLICIT NONE
    REAL :: PH2, X
    PH2 = (2*X + 5)**(1./3)
    RETURN
END FUNCTION

FUNCTION PH3(X)
    IMPLICIT NONE
    REAL :: PH3, X
    PH3 = (2*X + 5)/X**2
    RETURN
END FUNCTION

FUNCTION PH11(X)
    IMPLICIT NONE
    REAL :: PH11, X
    PH11 = 1.5*X**2
    RETURN
END FUNCTION

FUNCTION PH21(X)
    IMPLICIT NONE
    REAL :: PH21, X
    PH21 = (2./3)*(2*X + 5)**(-2./3)
    RETURN
END FUNCTION

FUNCTION PH31(X)
    IMPLICIT NONE
    REAL :: PH31, X
    PH31 = (-2./X**2) - (10./X**3)
    RETURN
END FUNCTION