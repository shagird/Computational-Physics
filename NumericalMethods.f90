function math_func (x)
    implicit none
    real*8 math_func
    real*8, intent(in) :: x 
    real*8 :: div
    div = 1.0d0 + x**2
    math_func = 1.0d0/div
end function math_func
function gauss(pi,x)
    implicit none
    real*8, intent(in) :: x, pi
    real*8 :: gauss
    gauss = (1/sqrt(2.0d0*pi)*exp(-1*(x**2)/2.0d0))
end function gauss
subroutine ez_integral(interval, ans)
    implicit none
    integer :: i, n 
    real*8 :: sum, math_func
    real*8, intent(in) :: interval
    real*8, intent(out) :: ans 
    sum = 0.0d0
    n = nint(1.0d0/interval)
    do i=1,n-1
        sum = sum + math_func(i*interval)
    end do
    ans = 4.0d0*(interval/2.0d0)*(math_func(0.d0) + math_func(1.0d0) + 2.0d0*sum)
end subroutine ez_integral
subroutine random_exp(lambda, ans)
    implicit none
    real*8, intent(in) :: lambda
    real*8, intent(out) :: ans 
    real*8 :: r 
    call random_number(r)
    ans = (-1.0d0/lambda)*log(r)
end subroutine random_exp
subroutine rand_normal(mu, sigma, z1,z2)
    implicit none
    real*8, intent(in) :: mu, sigma
    real*8, intent(out) :: z1, z2
    real*8 :: u1, u2, rad
    rad = 0.0d0 
    do while ((rad .le. 0.d0) .or. (rad .ge. 1))
        call random_number(u1)
        u1 = 2.0d0*u1 - 1.0d0
        call random_number(u2)
        u2 = 2.0d0*u2 - 1.0d0
        rad = u1**2 + u2**2
    end do
    z1 = sqrt(-2.0d0*log(rad))*(u1/sqrt(rad))
    z2 = sqrt(-2.0d0*log(rad))*(u2/sqrt(rad))
    z1 = mu + sigma*z1 
    z2 = mu + sigma*z2
end subroutine rand_normal
function func_g(x)
    implicit none
    real*8, dimension(6), intent(in) :: x
    real*8 :: func_g, xx, yy, xy
    xx = x(1)**2 + x(2)**2+ x(3)**2
    yy = x(4)**2 + x(5)**2+ x(6)**2 
    xy = (x(1) - x(4))**2 + (x(2) - x(5))**2 + (x(3) - x(6))**2
    func_g = exp(- xx - yy - 0.5d0*xy)
end function func_g
function func_imp_g(x)
    implicit none
    real*8, dimension(6), intent(in) :: x
    real*8 :: func_imp_g, xy
    xy = (x(1) - x(4))**2 + (x(2) - x(5))**2 + (x(3) - x(6))**2
    func_imp_g = exp(- 0.5d0*xy)
end function func_imp_g
!--------------------------------------------------------------------------------------------------------
program fortrine
    !4*integrate(1/(1+x^2))
    implicit none
    integer :: i,j,k,n 
    real :: mean, sq_sum, uni_rand(1:100000), sq_mean, stdev, summ1, c_k  
    real*8 :: math_func, gauss, h, summ, dx(1:6), new
    real*8, parameter :: pi = 2*asin(1.0d0)
    real*8  :: r,z1,z2,u1, u2, rad 
    real*8 :: x(6),func_g, func_imp_g, int_mc, var, siga, p(6), length, volume, ul 
    !q1a
    open(11, file='all_solutions.dat', status='new')
    h = 0.1d0
    summ = 0.0d0
    do i=1,9
        summ = summ + math_func(i*h)
    end do
    write(11,*) 'q1a: ', 4.0d0*(h/2.0d0)*(math_func(0.d0) + math_func(1.0d0) + 2.0d0*summ)
    !q1b
    open(unit=1, file ='trap_vs_interval.txt', status='new')
    dx = (/ 0.1d0, 0.01d0, 0.001d0, 0.0001d0, 0.00001d0, 0.000001d0/)
    do i=1,size(dx)
        call ez_integral(dx(i), new)
        write(1,*) nint(1.0d0/dx(i)), ' , ', abs(pi - new)
    end do
    print*, 'q1b done'
    close(1)
    !q1c
    !integrating sin(x) from0to /pi
    h = 0.000001d0!add two more zeros before 1
    summ = 0.0d0
    do i=1,nint(1.0d0/h) - 1
        summ = summ + sin(pi*i*h)
    end do
    write(11,*) 'q1c: ', (pi*h/2.0d0)*(sin(0.0d0) + sin(pi) + 2.0d0*summ)
    !q1d
    summ = 0.0d0
    do i=1,nint(1.0d0/h) - 1
        summ = summ + gauss(pi,-3.0d0 + 6.0d0*i*h)
    end do
    write(11,*) 'q1d: ', (6.0d0*h/2.0d0)*(gauss(pi,-3.0d0) + gauss(pi,3.0d0) + 2.0d0*summ)
    write(11,*) 'As expected the normalized gaussian function will give area near one for three sigma'
    !q1e
    !Error function is (2/sqrt(pi))*integral_(0 to z)(e^(-t^2))
    summ = 0.0d0
    do i=1,nint(1.0d0/h) - 1
        summ = summ + exp(-1*(i*h)**2)
    end do
    write(11,*) 'q1e: ', (2.0d0/sqrt(pi))*(h/2.0d0)*(exp(0.0d0) + exp(1.0d0) + 2.0d0*summ)
    !q2a&b
    open(2, file = 'uniform_rand_vals.txt', status='new')
    call random_number(uni_rand)
    write(2,'(f7.5)') uni_rand
    close(2)
    print*, 'q2a&b done'
    !q2d
    n = size(uni_rand)
    mean = sum(uni_rand)/real(n) 
    sq_sum = 0.0
    do i=1,n 
        sq_sum = sq_sum + uni_rand(i)**2
    end do 
    sq_mean = sq_sum/real(n)
    stdev = sqrt(sq_mean - mean**2)
    write(11,*) 'q2d mu: ', mean 
    write(11,*) 'q2d sigma: ', stdev
    !q2c
    open(3, file='ACF_vals.txt', status='new')
    do i=1,n 
        summ1 = 0.0
        k = i-1
        do j=1,n-k
            summ1 = summ1 + uni_rand(j)*uni_rand(j+k)
        end do 
        summ1 = summ1/real(n-k)
        c_k = (summ1 - mean**2)/(stdev**2)
        write(3,'(f8.6)') c_k
    end do
    close(3)
    print*, 'q2c done'
    !q3a
    call random_exp(0.01d0,r)
    write(11,*) 'q3a Exp random no.: ', r 
    open(10, file='expo_rand.txt', status='new')
    do i=1,n 
        call random_exp(0.01d0,r)
        write(10,*) r
    end do
    close(10)
    print*, 'q3a end'
    !q3b
    rad = 0.0d0 
    do while ((rad .le. 0.d0) .or. (rad .ge. 1))
        call random_number(u1)
        u1 = 2.0d0*u1 - 1.0d0
        call random_number(u2)
        u2 = 2.0d0*u2 - 1.0d0
        rad = u1**2 + u2**2
    end do
    z1 = sqrt(-2.0d0*log(rad))*(u1/sqrt(rad))
    z2 = sqrt(-2.0d0*log(rad))*(u2/sqrt(rad))
    write(11,*) 'q3b Normal random no.: ', z1, z2
    !q4a
    open(4, file='gaussian_distributed_2.txt', status='new')
    do i=1,n
        call rand_normal(0.0d0,2.0d0,z1,z2)
        write(4,'(2f10.6)') z1, z2
    end do
    close(4)
    print*, 'q4a done'
    !q4b
    open(7, file='gaussian_distributed_root2.txt', status='new')
    do i=1,n
        call rand_normal(0.0d0,1/sqrt(2.0d0),z1,z2)
        write(7,'(2f12.7)') z1, z2
    end do
    close(7)
    print*, 'q4b done'
    close(11)
    !q5
    !Brute force monte carlo integral
    !Limit our range from -5 to +5 since gaussian decays very rapidly

    n = 100 
    open(8, file='brute_force.txt', status='new')
    do while (n .le. 10**8)
        length = 5.0d0
        volume = (2.0d0*length)**6
        int_mc = 0.0d0
        var = 0.0d0
        siga = 0.0d0
        do i=1,n 
            call random_number(p)
            do j=1,6
                x(j) = -length + 2*length*p(j)
            end do
            int_mc = int_mc + func_g(x)
            siga = siga + func_g(x)**2
        end do
        int_mc = int_mc/real(n)
        siga = siga/real(n)
        var = siga - int_mc**2
        int_mc = volume*int_mc
        siga = volume*sqrt(var/real(n))
        write(8,*) n, ' : ', siga
        n = 10*n 
    end do
    close(8)
    print*, 'q5a done'
    !q5b
    n = 100
    open(9, file='imp_sampling.txt', status='new')
    do while (n .le. 10**8)
        volume = pi**3
        int_mc = 0.0d0
        var = 0.0d0
        siga = 0.0d0
        do i=1,n 
            do j=1,6
                call rand_normal(0.0d0,1.0d0/sqrt(2.0d0),x(j),ul)
            end do
            int_mc = int_mc + func_imp_g(x)
            siga = siga + func_imp_g(x)**2
        end do
        int_mc = int_mc/real(n)
        siga = siga/real(n)
        var = siga - int_mc**2
        int_mc = volume*int_mc
        siga = volume*sqrt(var/real(n))
        write(9,*) n , ' : ', siga 
        n = 10*n 
    end do
    close(9)
    print*, 'q5b done'
end program fortrine