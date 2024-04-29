program assingment_3
implicit none
integer::n
call adams_bash(real(n))
cALL Adams_multon(real(n))
CALL ADAMS_PREDICT(real(n))
end program

subroutine adams_bash(n)
IMPLICIT NONE
REAL::A=0,B=2,H=0.25,F,N,R1,R2,R3,R4
REAL,ALLOCATABLE::Y(:),T(:)
INTEGER::I
N=(B-A)/H
ALLOCATE(Y(0:N))
ALLOCATE(T(0:N))
Y(0)=0.5
T(0)=A
WRITE(*,*)'FROM RK METHOD'
DO I=1,3
    T(I)=A+I*H
    R1=H*F(T(I-1),Y(I-1))
    R2=H*F(T(I-1)+H/2,Y(I-1)+R1/2)
    R3=H*F(T(I-1)+H/2,Y(I-1)+R2/2)
    R4=H*F(T(I-1)+H,Y(I-1)+R3)

    Y(I)=Y(I-1)+(R1+2*R2+2*R3+R4)/6
    WRITE(*,*)I,T(I),Y(I)
END DO
    WRITE(*,*)'ADAMS-BASHFORTH 4 STEPS===================='
    DO I=4,N
        T(I)=A+I*H
        Y(I)=Y(I-1)+(H/24)*(55*F(T(I-1),Y(I-1))-59*F(T(I-2),Y(I-2))+37*F(T(I-3),Y(I-3))-9*F(T(I-4),Y(I-4)))
        WRITE(*,*)I,T(I),Y(I)
    END DO
        PRINT*,' '
        PRINT*,'SOLUTION IS=',Y(N)
end subroutine adams_bash

 subroutine adams_multon(n)
 IMPLICIT NONE
REAL::A=0,B=2,H=0.25,F,N,R1,R2,R3,R4
REAL,ALLOCATABLE::Y(:),T(:)
INTEGER::I
N=(B-A)/H
ALLOCATE(Y(0:N))
ALLOCATE(T(0:N))
Y(0)=0.5
T(0)=A
WRITE(*,*)' '
DO I=1,3
    T(I)=A+I*H
    R1=H*F(T(I-1),Y(I-1))
    R2=H*F(T(I-1)+H/2,Y(I-1)+R1/2)
    R3=H*F(T(I-1)+H/2,Y(I-1)+R2/2)
    R4=H*F(T(I-1)+H,Y(I-1)+R3)

    Y(I)=Y(I-1)+(R1+2*R2+2*R3+R4)/6
    END DO
    WRITE(*,*)'ADAMS MOULTON 4 STEPS====================='
    DO I=4,N
        T(I)=A+I*H

        Y(I)=(Y(I-1)+(H/720.)*(251.0*(1.0-T(I)**2)+646*F(T(I-1),Y(I-1))-&
        +264*F(T(I-2),Y(I-2))+106*F(T(I-3),Y(I-3))-19*F(T(I-4),Y(I-4))))/(1.-251*H/720.)
        WRITE(*,*)I,T(I),Y(I)

        END DO
        PRINT*,' '
        PRINT*,'SOLUTION IS=',Y(N)

 end subroutine adams_multon

 SUBROUTINE ADAMS_PREDICT(N)
    IMPLICIT NONE
REAL::A=0,B=2,H=0.2,F,N,R1,R2,R3,R4
REAL,ALLOCATABLE::Y(:),T(:)
INTEGER::I
N=(B-A)/H
ALLOCATE(Y(0:N))
ALLOCATE(T(0:N))
Y(0)=0.5
T(0)=A
         WRITE(*,*)" "
            WRITE(*,*)'ADAMS-PREDICTOR-CORRECTOR======================= '


DO I=1,3
    T(I)=A+I*H
    R1=H*F(T(I-1),Y(I-1))
    R2=H*F(T(I-1)+H/2,Y(I-1)+R1/2)
    R3=H*F(T(I-1)+H/2,Y(I-1)+R2/2)
    R4=H*F(T(I-1)+H,Y(I-1)+R3)

    Y(I)=Y(I-1)+(R1+2*R2+2*R3+R4)/6
    END DO
    DO I=4,N
        T(I)=A+I*H

        Y(I)=(Y(I-1)+(H/24.)*(55*F(T(I-1),Y(I-1))-59*F(T(I-2),Y(I-2))+&
        +37*F(T(I-3),Y(I-3))-9*F(T(I-4),Y(I-4))))

        Y(I)=(Y(I-1)+(H/24.)*(9*F(T(I),Y(I))+19*F(T(I-1),Y(I-1))-&
        +5*F(T(I-2),Y(I-2))+F(T(I-3),Y(I-3))))

        WRITE(*,*)I,T(I),Y(I)

        END DO
        PRINT*,' '
        PRINT*,'SOLUTION IS=',Y(N)

    END  SUBROUTINE ADAMS_PREDICT



REAL FUNCTION F(T,Y)
IMPLICIT NONE
REAL::Y,T
F=Y-T**2+1
RETURN
END FUNCTION
