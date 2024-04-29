PROGRAM ADAMS_4TH_ORDER_PREDICTOR_CORRECTOR
IMPLICIT NONE
    REAL::F,A,B,H,K1,K2,K3,K4,NN
    REAL,allocatable::Y(:),T(:)
    INTEGER::N,I,J

    A=0.
    B=2.
    H=0.2
    NN=(B-A)/H
    N=INT(NN) !NUMBER OF SUBINTERVALS

    ! DEFINING THE ARRAY SIZE
    allocate(Y(0:N))
    allocate(T(0:N))

    !FINDING THE TIME ENDPOINTS OF THE SUBINTERVAL
    DO I=0,N
        T(I)=A+I*H
    END DO

    PRINT*,' N t y(t)'

    !GIVEN INITIAL APPROXIMATION
    Y(0)=0.5
    PRINT*,0,' ',T(0),Y(0)

    !USING RUNGE-KUTTA METHOD FINDING FIRST 3 APROXIMATIONS
    DO I=1,3
        K1=H*F(T(I-1),Y(I-1))
        K2=H*F(T(I-1)+H/2.,Y(I-1)+K1/2.)
        K3=H*F(T(I-1)+H/2.,Y(I-1)+K2/2.)
        K4=H*F(T(I-1)+H,Y(I-1)+K3)
        Y(I)=Y(I-1)+(K1+2*K2+2*K3+K4)/6.
        PRINT*,I,' ',T(I),Y(I)
    END DO

    ! FINDING REMAINING APPROXIMATIONS BY ADAMS 4 ORDER PREDICTOR CORRECTOR ORDER METHOD
    DO I=4,N
        !PREDICTOR (SIMILAR ADAMS BASHFORTH)
        Y(I)=Y(I-1)+(H/24.)*(55*F(T(I-1),Y(I-1))-59*F(T(I-2),Y(I-2))+37*F(T(I-3),Y(I-3))-9*F(T(I-4),Y(I-4)))

        ! CORRECTOR (SIMILAR ADAMS_MOULTONS 3 STEP)
        Y(I)=Y(I-1)+(H/24.)*(9*F(T(I),Y(I))+19*F(T(I-1),Y(I-1))-5*F(T(I-2),Y(I-2))+F(T(I-3),Y(I-3)))

        PRINT*,I,' ',T(I),Y(I)
    END DO

    print*,' '
    PRINT*,'SOLUTION IS = ',Y(N)

END PROGRAM

FUNCTION F(T,Y)
    REAL::F,T,Y
    F=Y-T**2+1
END FUNCTION