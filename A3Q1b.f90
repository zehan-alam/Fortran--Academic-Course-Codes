PROGRAM MODIFIED_EULER
    IMPLICIT NONE

    REAL :: M_N, T ,W, H, A, B, N, ALPHA
    INTEGER :: I

    H = 0.2
    A = 1.0
    B = 2.0
    N = (B - A)/H

    ALPHA = 0.5
    W = ALPHA

    DO I = 0, INT(N - 1)
        T = A + I*H
        W = W + M_N(T, W)
    END DO

    WRITE(*,*) 'Y(2.0) ='
    WRITE(*,*) W


END PROGRAM


FUNCTION M_N(T, W)
    IMPLICIT NONE

    REAL :: M_N, T, W, H
    H = 0.2

    M_N = 0.1*((T**2) - T*W + (W**2) + ((T + H)**2) - (T + H)*(W + 0.2*((T**2) - T*W + (W**2))) + &
                   (W + .2*((T**2) - T*W + (W**2)))**2)


    RETURN
END FUNCTION
