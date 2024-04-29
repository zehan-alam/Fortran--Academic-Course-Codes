PROGRAM SYSTEM_JACOBI
    IMPLICIT NONE
    REAL :: X1, X2, X3, X(3), NORM, TOL
    INTEGER:: COUNT1
    X = [1, 0, 0]
    X1 = X(1) 
    X2 = X(2) 
    X3 = X(3) 
    TOL = 0.0000000001
    COUNT1 = 0
    DO
        X1 = (1./12)*(1 - 3*X2 + 5*X3)
        X2 = (1./5)*(28 - X1 - 3*X3)
        X3 = (1./13)*(76 - 3*X1 - 7*X2)
        COUNT1 = COUNT1 + 1
        NORM = SQRT((X1 - X(1))**2 + (X2 - X(2))**2 + (X3 - X(3))**2)
        IF(NORM < TOL)EXIT
        X(1) = X1
        X(2) = X2
        X(3) = X3
    END DO
    
    WRITE(*,*) 'APPROXIMATE SOLUTION:'
    WRITE(*,*) 'X1 = ',X1
    WRITE(*,*) 'X2 = ',X2
    WRITE(*,*) 'X3 = ',X3
    WRITE(*,*) 'ITERATIONS: ',COUNT1
    ! GAUSS - SEIDEL'S METHOD CONVERGES FASTER THAN JACOBI'S METHOD
END PROGRAM