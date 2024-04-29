PROGRAM POWER_METHOD
    IMPLICIT NONE
    REAL::A(3,3),X(3),Y(3),TMAX,TOL=0.001,M
    INTEGER::N=0,I,J

    OPEN(1,FILE='A2Q2input.txt')
    DO I=1,3
        READ(1,*)(A(I,J),J=1,3)
    END DO
    READ(1,*)(x(i),i=1,3)

    DO
        Y=MATMUL(A,X)
        TMAX=M
        M=MAX(Y(1),Y(2),Y(3))
        X=Y/M
        N=N+1

        WRITE(*,10)N,X(1),X(2),X(3),M
        10 FORMAT(I3,5X,'(',F7.5,',',F7.5,',',F7.5,')',5X,F10.5)
        IF(ABS(TMAX-M)<TOL) THEN
            PRINT*,' '
            PRINT*,'DOMINANT EIGENVALUE OF THE MATRIX A IS = ',M
            EXIT
        END IF
    END DO

    close(1)
END PROGRAM
