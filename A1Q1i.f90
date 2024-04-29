PROGRAM BISECTION
    IMPLICIT NONE
    REAL :: F, A, B, C, TEMP, TOL
    A = 1.5
    B = 2.5
    TOL = .000001
    
    IF(F(A)*F(B) < 0)THEN
        DO
            C = (A+B)/2
            IF(F(A)*F(C)<0)THEN
                B = C
                TEMP = A
            ELSE IF(F(B)*F(C) < 0)THEN
                A = C
                TEMP = B
            END IF

            IF(ABS(TEMP - C) < TOL)EXIT
        END DO

        WRITE(*,*) 'APPROXIMATE ROOT USING BISECTION: '
        WRITE(*,*) C
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