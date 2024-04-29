PROGRAM NEWTON
    REAL :: F1, F2, F3, F11, F21, F31
    REAL :: A, TEMP1, TOL, X1, X0, X11, X00
    REAL :: TEMP2, X12, X02, TEMP3
    TOL = 0.000001
    X0 = 1.5
    X00 = 1.5
    A = .5
    F31 = 1

    DO
        X1 = X0 - F1(X0)/F11(X0)
        TEMP1 = X0
        X0 = X1
        IF(ABS(X1 - TEMP1) < TOL)EXIT
    END DO
    WRITE(*,*) 'APPRX. SQUARE ROOT: '
    WRITE(*,*) X1

    DO
        X11 = X00 - F2(X00)/F21(X00)
        TEMP2 = X00
        X00 = X11
        IF(ABS(X11 - TEMP2) < TOL)EXIT
    END DO
    WRITE(*,*) 'APPRX. CUBIC ROOT: '
    WRITE(*,*) X11

    DO
        X12 = X02 - F3(X02)/F31
        TEMP3 = X02
        X02 = X12
        IF(ABS(X12 - TEMP3) < TOL)EXIT
    END DO
    WRITE(*,*) 'APPRX. INVERSE: '
    WRITE(*,*) X12
END PROGRAM

FUNCTION F1(X)
    IMPLICIT NONE
    REAL :: F1, X, A
    A = 0.5
    F1 = X**2 - A
    RETURN
END FUNCTION

FUNCTION F11(X)
    IMPLICIT NONE
    REAL :: F11, X, A
    A = 0.5
    F11 = 2*X
    RETURN
END FUNCTION

FUNCTION F2(X)
    IMPLICIT NONE
    REAL :: F2, X, A
    A = 0.5
    F2 = X**3 - A
    RETURN
END FUNCTION

FUNCTION F21(X)
    IMPLICIT NONE
    REAL :: F21, X, A
    A = 0.5
    F21 = 3*X**2
    RETURN
END FUNCTION

FUNCTION F3(X)
    IMPLICIT NONE
    REAL :: F3, X, A
    A = 0.5
    F3 = X - 1/A
    RETURN
END FUNCTION