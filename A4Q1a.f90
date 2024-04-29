program linearshooting
    implicit none

    integer::i,n

    real::a,b,h,x,y,alpha,beta,p,q,r,w1,w2
    real,dimension(4,2)::k,kp
    real,dimension(2,0:1000)::u,v,w


    a=1.
    b=2.


    alpha=.5
    beta=2.
    N=10

    h=(b-a)/N
    u(1,0)=alpha
    u(2,0)=0
    v(1,0)=0
    v(2,0)=1

    do i=0,N-1
        x=a+i*h
        k(1,1) = h*u(2,i)
        k(1,2) = h*(p(x)*u(2,i) + q(x)*u(1,i) + r(x))

        k(2,1) = h*(u(2,i) + k(1,2)/2)
        k(2,2) = h*(p(x + h/2) *(u(2,i) + k(1,2)/2)+q(x + h/2))*(u(1,i) +k(1,1)/2 + r(x + h/2))

        k(3,1) = h*(u(2,i) + k(2,2)/2)
        k(3,2) = h*(p(x + h/2) *(u(2,i) + k(2,2)/2)+q(x + h/2))*(u(1,i) +k(2,1)/2 + r(x + h/2))

        k(4,1) = h*(u(2,i) + k(3,2))
        k(4,2) = h*(p(x + h) *(u(2,i) + k(3,2))+q(x + h))*(u(1,i) +k(3,1) + r(x + h))

        u(1,i+1) = u(1,i) + (k(1,1) + 2*k(2,1) + 2*k(3,1) + k(4,1))/6
        u(2,i+1) = u(2,i) + (k(1,2) + 2*k(2,2) + 2*k(3,2) + k(4,2))/6





        kp(1,1) = h*v(2,i)
        kp(1,2) = h*(p(x)*v(2,i) + q(x)*v(1,i) )

        kp(2,1) = h*(v(2,i) + kp(1,2)/2)
        kp(2,2) = h*(p(x + h/2) *(v(2,i) + kp(1,2)/2)+q(x + h/2))*(v(1,i) +kp(1,1)/2)

        kp(3,1) = h*(v(2,i) + kp(2,2)/2)
        kp(3,2) = h*(p(x + h/2) *(v(2,i) + kp(2,2)/2)+q(x + h/2))*(v(1,i) +kp(2,1)/2)

        kp(4,1) = h*(v(2,i) + kp(3,2))
        kp(4,2) = h*(p(x + h) *(v(2,i) + kp(3,2))+q(x + h))*(v(1,i) +kp(3,1))

        v(1,i+1) = v(1,i) + (kp(1,1) + 2*kp(2,1) + 2*kp(3,1) + kp(4,1))/6
        v(2,i+1) = v(2,i) + (kp(1,2) + 2*kp(2,2) + 2*kp(3,2) + kp(4,2))/6

        end do
        w(1,0)=alpha
        w(2,0)=(beta - u(1,N))/v(1,N)

        write(*,*)"     x         Exact solution  LSM Approx.     Error in LSM"
        write(*,8)a,y(a),w(1,0),abs(y(x)-w(1,0))
        do i=1,N
            W1=u(1,i)+w(2,0)*v(1,i)
            W2=u(2,i)+w(2,0)*v(2,i)
            x=a+i*h
            write(*,8) x,y(x),w1,abs(y(x)-w1)

        end do
        8 format(4(f10.5,6x))
end program

function p(x)
    real::p,x
    p=-(2./x)
end function

function q(x)
    real::q,x
    q=2./ x**2
end function

function r(x)
    real::r,x
    r=sin(log(x)) / x**2
end function

function y(x)
    real::y,x,c1,c2
    c1=-0.03920701320
    c2=1.1392070132
    y=c2*x +c1*(1./x**2) -(3./10)*sin(log(x)) -(1./10)*cos(log(x))
end function