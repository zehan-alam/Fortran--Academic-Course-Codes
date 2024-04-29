Program karmack
    Implicit none
    Real:: s(100), i(100), r(100), b, g, dt=0.1, ds, di, dr, t
    Integer:: j, n=100
    b= 10.0
    g= 1.0
    t=0.0
    
    s(1)= 0.9
    i(1)= 0.1
    r(1)= 0.0
    
    do j=1, n+1
    ds= -1*b*s(j)*i(j)
    dr= g*i(j)
    di=-1*(ds+dr)
    print*, j-1, t, s(j), r(j), i(j)
    
    s(j+1)= s(j) + dt*ds
    r(j+1)= r(j) + dt*dr
    i(j+1)= i(j) + dt*di
    
    
    t= t+dt
    end do
    print*, j-1, t, s(j), r(j), i(j)
    
    end program