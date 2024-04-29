program power
    implicit none

    real::A(3,3),X(3),Y(3),temp,M,tol=0.001
    integer::n=0,i,j

    OPEN(1,FILE='A2Q2input.txt')
    DO I=1,3
        READ(1,*)(A(I,J),J=1,3)
    END DO
    READ(1,*)(x(i),i=1,3)

    do
      Y=matmul(A,X)
      temp=M
      M=(max(abs(Y(1)),abs(Y(2)),abs(Y(3))))
      X=Y/M
      n=n+1

      WRITE(*,10)N,X(1),X(2),X(3),M
        10 FORMAT(I3,5X,'(',F7.5,',',F7.5,',',F7.5,')',5X,F10.5)
    
        if(abs(temp-M)<tol)then
            PRINT*,'DOMINANT EIGENVALUE OF THE MATRIX A IS = ',M
            EXIT
        end if
    end do


    close(1)
    stop
end program power