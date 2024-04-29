PROGRAM TAYLOR
    IMPLICIT NONE

    REAL :: T_N, T, W, H, N, A, B, ALPHA
    INTEGER :: I

    H = 0.2
    A = 1.0
    B = 2.0
    N = (B - A)/H

    ALPHA = .5
    W = ALPHA

    DO I = 0, INT(N - 1)
        T = A + H*I
        W = W + H*T_N(T, W)
    END DO

    WRITE(*,*) 'Y(2.0) ='
    WRITE(*,*) W

END PROGRAM

FUNCTION T_N(T, W)
    IMPLICIT NONE

    REAL :: T_N, T, W, H

    H = 0.2

    T_N = (T**2 - T*W +W**2) + (H/2)*(2*T - T**3 +3*(T**2)*W - 3*T*(W**2) - W + 2*W**3) + &
           ((H**2)/6)*(2 - 4*(T**2) + 3*(T**4) - 9*(T**3)*W + 15*(T**2)*(W**2) + 7*T*W - &
           12*T*(W**3) - 4*(W**2) + 6*(W**4))

END FUNCTION
