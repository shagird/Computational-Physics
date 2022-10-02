subroutine ChangeSeed(k,l,seed_size)
    implicit none

    integer,intent(in) :: k,l,seed_size
    integer :: j,maxa, mina
    integer, dimension(1:seed_size) :: seed_new
    !integer, dimension(2) :: maxa, mina
    maxa = maxval((/l,seed_size/))
    mina = minval((/k,1/))
    do j = mina,maxa
        seed_new(j) = int(rand()*10**(10-(l-j)))
    end do
    call random_seed(put = seed_new)
end subroutine ChangeSeed

program fnew
    implicit none
    integer :: a,i,j,k,l,n,num2
    integer, allocatable:: seed(:), seed_old(:)
    real, allocatable :: for_bin(:)
    integer, dimension(1:8) :: datetime
    integer, dimension(1:3) :: arr
    integer, dimension(2) :: arr_new
    real :: r, mat(10,10), avg, summ, bw, num
    real*8 :: rd

    !q0
    k = 1
    rd = 5.2333365616519846
    r = 10.6
    print*, 'Raghav Sharma'
    print*, k,r,rd
    !q1a
    do i=1,10
        call random_number(rd)
        print*, rd
    end do
    !q1b
    open(unit = 1, file = 'test_ran.dat', status = 'new')
    do i = 1,10
        call random_number(r)
        write(1,*) r
    end do
    !q1c&d
    write(1,*) '!Changing seed and generating 10 new random numbers'
    call random_seed(size = n)
    !print*, n
    allocate(seed(n), seed_old(n))
    call random_seed(get = seed)
    seed_old = seed
    call date_and_time(values = datetime)
    seed(1:8) = datetime
    call random_seed(put = seed)
    do i = 1,10
        call random_number(r)
        write(1,*) r
    end do
    !close(1)

    !Create a new file
    open(2, file = 'test_ran_10_seeds.dat', status = 'new')
    do i=1,10
        call ChangeSeed(i,i+5,12)
        do j = 1,10
            call random_number(r)
            mat(i,j) = r
        end do
        !call random_seed(put = seed_old)
    end do

    write(2,'(10f14.10)') mat
    close(2)
    !q1e
    call random_seed(put = seed_old)
    !open(1,file='test_ran.dat', status = 'old')
    write(1,*) '!NOW calculating average of 10 random numbers'
    n = 10
    summ = 0
    do i =1,n
        call random_number(r)
        summ = summ + r
    end do
    avg = summ/n
    write(1,*) avg
    arr = (/ 100,10000,1000000 /)
    do j = 1,size(arr)
        n = arr(j)
        summ = 0
        do i =1,n
            call random_number(r)
            summ = summ + r
        end do
        avg = summ/n
        write(1,*) 'For ',n,' average: ',avg
    end do
    do j = 1,size(arr)
        n = arr(j)
        summ = 0
        do i =1,n
            call random_number(r)
            summ = summ + r
        end do
        avg = summ/n
        write(1,*) 'For ',n,' deviation: ',abs(0.5d0 - avg)
    end do
    write(1,*) 'Observation: Mean is very close to 0.5'
    write(1,*) 'And as we increase the number of observations we get closer to 0.5 which is the true mean of uniform[0, 1]'
    write(1,*) 'We also see the pattern that the deviation depends square root of number of iterations.'
    close(1)
    !q1h
    open(4, file = 'dist_sum_random.txt', status='new')
    l = 10000
    n = 10000
    allocate(for_bin(l))
    do i=1,l
        summ = 0
        do j =1,n
            call random_number(r)
            summ = summ + r
        end do
        for_bin(i) = summ
        write(4,*) summ
    end do
    close(4)

    open(41, file = 'dist_sum_random_bin_one.txt', status = 'new')!I have done binning in python this code for binning is a proof of concept.

    bw = 1
    num = minval(for_bin)
    do while (num .le. maxval(for_bin))
        num2 = 0
        do i =1,l
            if ((num < for_bin(i)) .and. (num + bw .ge. for_bin(i))) then
                num2 = num2 + 1
            end if
        end do
        write(41,*) num, num2
        num = num+bw
    end do
    close(41)

    open(7, file = 'dist_sum_random_1.txt', status='new')

    arr_new = (/10000, 100000/)
    n = 10000
    do k = 1, size(arr_new)
        do i=1,n
            summ = 0
            do j =1,arr_new(k)
                call random_number(r)
                r = 2*r - 1
                summ = summ + r
            end do
            write(7,*) summ
        end do
    end do
    close(7)
    !q1i
    l = 10000
    n = 10000
    open(8, file = 'dist_sum_random_walk.txt', status = 'new')
    do j = 1,l
        a = 0
        do i = 1,n
            call random_number(r)
            if (r < 0.5) then
                a = a - 1
            else
                a = a + 1
            end if
        end do
        write(8,*) a
    end do
    close(8)
    !q1k
    l = 100000
    n = 10000
    open(9, file = 'dist_sum_random_walk1.txt', status = 'new')
    do j = 1,l
        a = 0
        do i = 1,n
            call random_number(r)
            if (r < 0.5) then
                a = a - 1
            else
                a = a + 1
            end if
        end do
        write(9,*) a
    end do
    close(9)
    !q1l
    l = 100000
    n = 100000
    open(10, file = 'dist_sum_random_walk2.txt', status = 'new')
    do j = 1,l
        a = 0
        do i = 1,n
            call random_number(r)
            if (r < 0.5) then
                a = a - 1
            else
                a = a + 1
            end if
        end do
        write(10,*) a
    end do
    close(10)
end program fnew
