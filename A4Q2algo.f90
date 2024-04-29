PROGRAM power_method
    IMPLICIT NONE
    INTEGER,PARAMETER::n=3
    INTEGER::i,j,k,maxn,p

    REAL:: mu,tol,norm
    REAL,DIMENSION(n,n) :: a
    REAL,DIMENSION(n) :: x,xn,y

    OPEN(1,FILE='A2Q2input.txt')

    READ(1,*)((a(i,j),j=1,n),i=1,n)
    READ(1,*)(x(i),i=1,n)

    WRITE(*,7)"A=",((a(i,j),j=1,n),i=1,n)

    maxn=100
    tol=0.0001

    WRITE(*,*)'step          x1              x2           x3            mu'
    WRITE(*,8)0,(x(i),i=1,n)

    DO p=1,n
        IF(ABS(x(p))==MAXVAL(ABS(x))) EXIT
    END DO

    DO k=1,maxn
        y=MATMUL(a,x)
        mu=y(p)

        DO p=1,n
            IF(ABS(y(p))==MAXVAL(ABS(y))) EXIT
        END DO

        IF(y(p)==0) THEN
            WRITE(*,*)'eigenvector ',x
            WRITE(*,*)'A has the eigenvalue 0, select a new vector x and restart.'
            STOP
        END IF

        xn=y/y(p)
        norm=MAXVAL(ABS(x-xn))
        x=xn

        WRITE(*,8)k,(x(i),i=1,n),mu

        IF(norm<tol)THEN
            WRITE(*,9)"Dominant eigenvalue=",mu
            WRITE(*,9)"eigenvector=",(x(i),i=1,n)
            STOP
        END IF

    END DO
    WRITE(*,*)'Max number of iteration exceeded.'

    7 FORMAT(a,/,3(3(f12.7,3x),/))
    8 FORMAT(i4,4(3x,f12.7))
    9 FORMAT(a,/,3(f12.7,/))

END PROGRAM