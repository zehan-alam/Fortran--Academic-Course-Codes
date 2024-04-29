PROGRAM SECANT
    IMPLICIT NONE
    REAL :: F, X1, X0, X00, TEMP, TOL
    X0 = 2.0
    X00 = 2.5
    TOL = 0.0001

    DO
        X1 = X0 - F(X0)*((X0 - X00)/(F(X0) - F(X00)))
        TEMP = X0
        X0 = X1
        X00 = TEMP
        IF(ABS(X1 - TEMP) < TOL)EXIT
    END DO
    WRITE(*,*) 'APPROXIMATE ROOT USING SECENT:'
    WRITE(*,*) X1
END PROGRAM


FUNCTION F(X)
IMPLICIT NONE
    REAL :: F, X
    F = X**3 -2*X -5
    RETURN
END FUNCTION