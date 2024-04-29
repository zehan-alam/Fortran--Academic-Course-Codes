program hello
  implicit none

  real::g1,g11,g12,g2,g21,g22,x0,y0,tol,n,v1,v2,v3,v4,k,x1,y1,t1,t2,norm
  integer::i
  X0 = .25
  Y0 = .25
  G11 = 0
  G12 = 1./SQRT(5.)
  N = 2 ! NUMBER OF EQUATIONS
  TOL = 0.0000001
  v1=g11
  v2=g12
  v3=abs(g21(x0))
  v4=abs(g22(y0))
  k=max(v1,v2,v3,v4)*N

  if(k>=1)then
    write(*,*)"Not converging"
  else
    do i=1,5
      x1=g1(y0)
      y1=g2(x0,y0)
      t1=x0
      t2=y0
      x0=x1
      y0=y1
      norm=sqrt((t1-x1)**2 +(t2-y1)**2)
      if(norm<tol)exit
    end do
    WRITE(*,*) 'APPROXIMATIONS:'
    WRITE(*,*) 'X = ',X1
    WRITE(*,*) 'Y = ',Y1
  end if
  stop
end program hello

FUNCTION G1(Y)
  IMPLICIT NONE
  REAL :: G1,Y
  G1 = Y/SQRT(5.)
  RETURN
END FUNCTION

FUNCTION G2(X, Y)
  IMPLICIT NONE
  REAL :: G2,X, Y
  G2 = .25*(SIN(X) + COS(Y))
  RETURN
END FUNCTION

FUNCTION G21(X)
  IMPLICIT NONE
  REAL :: G21, X
  G21 = .25*COS(X)
  RETURN
END FUNCTION

FUNCTION G22(Y)
  IMPLICIT NONE
  REAL :: G22, Y
  G22 = -.25*SIN(Y)
  RETURN
END FUNCTION