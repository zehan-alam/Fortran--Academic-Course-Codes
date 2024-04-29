program ansto1b
    implicit none
    integer :: i, N
    real ::  h, x, alpha, delta
    real :: y, p, q, r
    real, dimension(:), allocatable :: a, b !(0 to N)
    real, dimension(:), allocatable :: c, d, l, u, z !(1 to N)
    real, dimension(:), allocatable :: w !(0 to N+1)
    
    !input
    N = 10

    allocate(a(0 : N), b(0 : N))
    allocate(c(N), d(N), l(N), u(N), z(N))
    allocate(w(0 : N + 1))

    a(0) = 1.0
    b(0) = 2.0
    alpha = 1.0
    delta = 2.0
    
    h = (b(0) - a(0)) / (N + 1)
    x = a(0) + h
    a(1) = 2.0 + h ** 2.0 * q(x)
    b(1) = - 1.0 + 0.5 * h * p(x)
    d(1) = - h ** 2.0 * r(x) + (1.0 + 0.5 * h * p(x)) * alpha

    do i = 2, N - 1
        x = a(0) + i * h
        a(i) = 2.0 + h ** 2.0 * q(x)
        b(i) = - 1.0 + 0.5 * h * p(x)
        c(i) = - 1.0 - 0.5 * h * p(x)
        d(i) = - h ** 2.0 * r(x)
    end do

    x = b(0) - h
    a(N) = 2.0 + h ** 2.0 * q(x)
    c(N) = - 1.0 - 0.5 * h * p(x)
    d(N) = - h ** 2.0 * r(x) + (1.0 - 0.5 * h * p(x)) * delta

    l(1) = a(1)
    u(1) = b(1) / a(1)
    z(1) = d(1) / l(1)

    do i = 2, N - 1
        l(i) = a(i) - c(i) * u(i - 1)
        u(i) = b(i) / l(i)
        z(i) = (d(i) - c(i) * z(i - 1)) / l(i)
    end do

    l(N) = a(N) - c(N) * u(N - 1)
    z(N) = (d(N) - c(N) * z(N - 1)) / l(N)

    w(0) = alpha
    w(N + 1) = delta
    w(N) = z(N)
    
    do i = N - 1, 1, - 1
        w(i) = z(i) - u(i) * w(i + 1)
        ! write(*,*)w(i),z(i),u(i),w(i+1)
    end do

    do i = 0, N + 1
        x = a(0) + i * h
        write(*,*) x, w(i), y(x), abs(y(x) - w(i))
    end do

end program ansto1b

real function p(x)
    real :: x
    p = -(2.0 / x)
end function p

real function q(x)
    real :: x
    q = 2.0 / x ** 2.0
end function q

real function r(x)
    real :: x
    r = sin(log(x)) / x ** 2.0
end function r

real function y(x)
    real :: x, C1, C2
    C1 = 1.1392070132
    C2 = -(0.03920701320)
    y = C1 * x + C2 * (1 / x ** 2.0) - (3.0 / 10.0) * sin(log(x)) - (1.0 / 10.0) * cos(log(x))
end function y


