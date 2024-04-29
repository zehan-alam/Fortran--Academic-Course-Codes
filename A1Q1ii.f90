PROGRAM FALSE_POSITION
    IMPLICIT NONE
    REAL :: F, A, B, X1
    REAL :: TOL, TEMP

    A = 1.5
    B = 2.5
    TOL = .000000001
    TEMP = 100000000

    IF(F(A)*F(B) < 0)THEN
        DO
            X1 = (A*F(B) - B*F(A))/(F(B) - F(A))
            IF(F(X1)*F(A) < 0)THEN
                B = X1
            ELSEIF(F(X1)*F(B) < 0)THEN
                A = X1
            END IF
            IF(ABS(X1 - TEMP) < TOL)EXIT
            TEMP = X1
        END DO

        WRITE(*,*) 'APPROXIMATE ROOT USING FALSE POSITION: '
        WRITE(*,*) X1
    ELSE
        WRITE(*,*) 'CHANGE INTERVAL'
    END IF
END PROGRAM


FUNCTION F(X)
IMPLICIT NONE
    REAL :: F, X
    F = X**3 - 2*X - 5
RETURN
END FUNCTION