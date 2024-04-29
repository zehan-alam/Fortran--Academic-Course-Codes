PROGRAM ADAMS_BASH_4
    IMPLICIT NONE

    REAL :: F, T, W, K1, K2, W2
    REAL :: K3, K4, K5, K6, A, B, H, N, ALPHA
    REAL, ALLOCATABLE :: WS(:)
    INTEGER :: I

    A = 0.0
    B = 2.0
    H = 0.25
    N = (B - A)/H
    ALLOCATE(WS(0:INT(N)))

    ALPHA = 0.5
    WS(0) = ALPHA
    W = ALPHA

    DO I = 0, 3
        T = A + I*H
        K1 = H*F(T, W)
        K3 = H*F(T + 3*H/8, W + (3./32)*K1 + (9./32)*K2)
        K4 = H*F(T + 12*H/13, W + (1932./2197)*K1 - (7200./2197)*K2 + (7296./2197)*K3)
        K5 = H*F(T + H, W + (439./216)*K1 - 8*K2 + (3680./513)*K3 - (845./4104)*K4)
        K6 = H*F(T + H/2, W - (8./27)*K1 + 2*K2 - (3544./2565)*K3 + (1859./4104)*K4 - (11./40)*K5)
        W = W + (16./135)*K1 + (6656./12825)*K3 + (28561./56430)*K4 - (9./50)*K5 + (2./55)*K6
        WS(I+1) = W
    END DO

    DO I = 4, INT(N - 1)
        T = A + I*H
        W2 = WS(I) + (H/720)*(1901*F(T, WS(I)) - 2774*F(T - H, WS(I - 1)) + 2616*F(T - 2*H, WS(I - 2)) - &
             1274*F(T - 3*H, WS(I - 3)) + 251*F(T - 4*H, WS(I - 4)))
        WS(I + 1) = W2
    END DO

    WRITE(*,*) 'Y(2.0) ='
    WRITE(*,*) WS(INT(N))

  END PROGRAM

  FUNCTION F(T1, W1)
    IMPLICIT NONE

    REAL :: F, T1, W1

    F = W1 - T1**2 + 1

    RETURN
  END FUNCTION

