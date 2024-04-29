PROGRAM ADAMS_BASH_4
    IMPLICIT NONE

    REAL :: F, T, W, T1, W1, K1, K2, W2
    REAL :: K3, K4, A, B, H, N, ALPHA
    REAL, ALLOCATABLE :: WS(:)
    INTEGER :: I

    A = 0.0
    B = 2.0
    H = 0.25
    N = (B - A)/H
    ALLOCATE(WS(0: int(N)))

    ALPHA = 0.5
    WS(0) = ALPHA
    W = ALPHA

    DO I = 0, 2
        T = A + I*H
        K1 = H*F(T, W)
        K2 = H*F(T + H/2, W + K1/2)
        K3 = H*F(T + H/2, W + K2/2)
        K4 = H*F(T + H, W + K3)
        W = W + (1./6)*(K1 + 2*K2 + 2*K3 + K4)
        WS(I+1) = W
        WRITE(*,*) WS(I + 1)
    END DO

    DO I = 3, INT(N - 1)
        T = A + I*H
        W2 = WS(I) + (H/24)*(55*F(T, WS(I)) - 59*F(T - H, WS(I - 1)) + 37*F(T - 2*H, WS(I - 2)) - &
                             9*F(T - 3*H, WS(I - 3)))
        WS(I + 1) = W2
        WRITE(*,*) WS(I + 1)
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
