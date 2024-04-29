PROGRAM GAUSS_ELIM
    IMPLICIT NONE
    REAL :: A(3, 4), TEMP1, TEMP2, TEMP3
    INTEGER :: I, J ,K

    OPEN(UNIT=101,FILE='A2Q1INPUT.TXT')
    !OPEN(UNIT=201,FILE='A2Q1OUTPUT.TXT')

    DO I = 1, 3
        READ(101,*)(A(I, J), J = 1, 4)
    END DO

    DO I = 1, 3
        TEMP1 = A(I, I)
        DO K = I, 4
            A(I, K) = A(I, K)/TEMP1
        END DO
        DO J = I+1, 3
            TEMP2 = A(J, I)
            DO K = I, 4
                A(J, K) = A(J, K) -TEMP2*A(I, K)
            END DO
        END DO
    END DO

    DO I = 3, 1, -1
        DO J = I-1, 1, -1
            TEMP3 = A(J, I)
            DO K = 4, 1, -1
                A(J, K) =  A(J, K) -TEMP3*A(I, K)
            END DO
        END DO
    END DO

    DO I = 1, 3
        WRITE(*,*)(A(I, J), J = 1, 4)
    END DO

    WRITE(*, *) 'SOLUTION: '
    WRITE(*, *) 'X = ',A(1,4)
    WRITE(*, *) 'Y = ',A(2,4)
    WRITE(*, *) 'Z = ',A(3,4)

    CLOSE(101)
    !CLOSE(201)
END PROGRAM