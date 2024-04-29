PROGRAM ADAMS_MOULTON_3
    IMPLICIT NONE

    REAL :: F, T, W, T1, W1, K1, K2, W2
    REAL :: K3, K4, A, B, H, N, ALPHA
    REAL, ALLOCATABLE :: WS(:)
    INTEGER :: I

    A = 0.0
    B = 2.0
    H = 0.25
    N = (B - A)/H
    ALLOCATE(WS(0: N))

    ALPHA = 0.5
    WS(0) = ALPHA
    W = ALPHA

    DO I = 0, 1
        T = A + I*H
        K1 = H*F(T, W)
        K2 = H*F(T + H/2, W + K1/2)
        K3 = H*F(T + H, W - K1 + 2*K2)
        W = W + (1./6)*(K1 + 4*K2 + K3)
        WS(I+1) = W
        write(*,*)WS(I+1)
    END DO

    DO I = 2, INT(N - 1)
        T = A + I*H
        W2 = (1/(1 - (9*H/24)))*(WS(I) + 9*H/24*(1 - (T + H)**2)+(H/24)*(19*F(T, WS(I)) - &
                                 5*F(T - H, WS(I - 1)) + F(T - 2*H, WS(I - 2))))
        WS(I + 1) = W2
        write(*,*)WS(I+1)
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


