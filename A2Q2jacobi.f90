program jocabi
    implicit none
    real::x=1.5,y=1.5,z=-.5,A(3),tol=0.0000001,norm,c=0
    A=[1.5,1.5,-.5]

    do
        x = (1./8)*(11 +3*A(2) -2*A(3))
        y = (1./17)*(7 -2*A(1) -14*A(3))
        z = (1./10)*(-10 +2*A(1) -4*A(2))
      c=c+1
        norm=sqrt((A(1)-x)**2 +(A(2)-y)**2 + (A(3)-z)**2)
        if(norm<tol)exit
        A(1)=x
        A(2)=y
        A(3)=z
    end do
    write(*,*)x,y,z,c
end program