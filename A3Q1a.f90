PROGRAM EULER
    IMPLICIT NONE

    REAL :: F, T, H, A, ALPHA, W, B, N
    INTEGER :: I

    A = 1.
    B = 2.
    ALPHA = .5
    H = .2
    N = ((B - A)/H)

    W = ALPHA

    DO I = 0, INT(N - 1)
        T = A + I*H
        W = W + H*F(T, W)
    END DO

    WRITE(*,*) 'Y(2.0) = '
    WRITE(*,*) W

END PROGRAM

FUNCTION F(T, Y)
    IMPLICIT NONE

    REAL :: F, T, Y

    F = T**2 - T*Y + Y**2

END FUNCTION
