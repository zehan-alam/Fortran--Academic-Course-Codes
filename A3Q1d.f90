PROGRAM RUNGE_KUTTA
    IMPLICIT NONE

    REAL :: W, K1, K2, K3, K4, N, H, A, B, ALPHA
    REAL :: F, T, T1, W1
    INTEGER :: I

    A = 1.0
    B = 2.0
    H = 0.2
    N = (B - A)/H

    ALPHA = 0.5
    W = ALPHA

    DO I = 0, INT(N - 1)
        T = A + I*H
        K1 = H*F(T, W)
        K2 = H*F(T + H/2, W + K1/2)
        K3 = H*F(T + H/2, W + K2/2)
        K4 = H*F(T + H, W + K3)
        W = W + (1./6)*(K1 + 2*K2 + 2*K3 + K4)
    END DO

    WRITE(*,*) 'Y(2.0) ='
    WRITE(*,*) W

END PROGRAM

FUNCTION F(T1, W1)
    IMPLICIT NONE

    REAL :: F, T1, W1

    F = T1**2 - T1*W1 + W1**2

    RETURN
END FUNCTION

