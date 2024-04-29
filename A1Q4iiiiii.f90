PROGRAM HERO
    IMPLICIT NONE
    REAL :: A, X1, X0, TEMP, TOL
    A = .3
    X0 = 1.5
    TOL = 0.0000000001

    
    ! (I)
    DO
        X1 = .5*(X0 + A/X0)
        TEMP = X0
        X0 = X1
        IF(ABS(X1 - TEMP) < TOL)EXIT
    END DO
    WRITE(*,*) 'APPROXIMATE SQUARE ROOT: (USING HERO''S ALGORITHM)'
    WRITE(*,*) X1


    ! (II)
    DO
        X1 = (1./3)*(2*X0 + A/X0**2)
        TEMP = X0
        X0 = X1
        IF(ABS(X1 - TEMP) < TOL)EXIT
    END DO
    WRITE(*,*) 'APPROXIMATE CUBIC ROOT: (USING SECOND ALGORITHM)'
    WRITE(*,*) X1


    ! (III)
    DO
        X1 = X0*(2 - A*X0)
        TEMP = X0
        X0 = X1
        IF(ABS(X1 - TEMP) < TOL)EXIT
    END DO
    WRITE(*,*) 'APPROXIMATION OF INVERSE OF ',A, ':'
    WRITE(*,*) X1
END PROGRAM