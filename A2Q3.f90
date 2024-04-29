program Nrm
    implicit none
    real::f1,f11,f12,f2,f21,f22
    real::x0(2),x1(2),y0(2),f(2),jacobi(2,2),inv_jacobi(2,2)
    integer::i,j,c=0

    x0=[2.0 , 1.0]

    do
      jacobi(1,1)=f11(x0(1),x0(2))
      jacobi(1,2)=f12(x0(1),x0(2))
      jacobi(2,1)=f21(x0(1),x0(2))
      jacobi(2,2)=f22(x0(1),x0(2))
      
      c=c+1
      if(c==1)then
        write(*,*)"Jacboian at (2,1): "
        do i=1,2
          write(*,*)(jacobi(i,j),j=1,2)
        end do
      end if

      f=[f1(x0(1),x0(2)),f2(x0(1),x0(2))]
      call inv(jacobi,inv_jacobi)
      y0=-matmul(inv_jacobi,f)

      do i=1,2
        x1(i)=x0(i)
        x0(i)=x0(i)+y0(i)
      end do
      if(abs(x0(1)-x1(1))< .0001 .and. abs(x0(2)-x1(2))< .0001)exit
    end do

    write(*,*)x0
end program

function f1(x,y)
  real::x,y,f1
  real,parameter::pi=3.1416
  f1=log(x- y**2) -sin(x*y) -sin(pi)
  return
end function

function f11(x,y)
  real::x,y,f11
  real,parameter::pi=3.1416
  f11=1/(x-y**2) - y*cos(x*y)
  return
end function

function f12(x,y)
  real::x,y,f12
  real,parameter::pi=3.1416
  f12=-(2*y)/(x-y**2) - x*cos(x*y) 
  return
end function

function f2(x,y)
  real::x,y,f2
  real,parameter::pi=3.1416
  f2=exp(x*y) + cos(x-y) -2
  return
end function

function f21(x,y)
  real::x,y,f21
  real,parameter::pi=3.1416
  f21=y*exp(x*y)- sin(x-y)
  return
end function

function f22(x,y)
  real::x,y,f22
  real,parameter::pi=3.1416
  f22=x*exp(x*y) + sin(x-y)
  return
end function

subroutine inv(jacobi,inv_jacobi)
  implicit none
  real::jacobi(2,2),inv_jacobi(2,2),det
  det=jacobi(1,1)*jacobi(2,2) - jacobi(1,2)*jacobi(2,1)
  inv_jacobi(1,1)=(1./det)*jacobi(2,2)
  inv_jacobi(1,2)=-(1./det)*jacobi(1,2)
  inv_jacobi(2,1)=-(1./det)*jacobi(2,1)
  inv_jacobi(2,2)=(1./det)*jacobi(1,1)
  return
end subroutine

subroutine MAT_MUL(inv_jacobi,f,y0)
  implicit none
  integer::i,j
  real::inv_jacobi(2,2),f(2),y0(2),sum

  do i=1,2
    sum=0
    do j=1,2
      sum=sum+ inv_jacobi(i,j)*f(j)
    end do
    y0(i)=-sum
  end do
  return
end subroutine