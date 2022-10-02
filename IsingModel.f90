! k_B, J_ising = 1 by default
!***************************************************************************************************************
program fortrine
    implicit none
    integer :: i,j,k,L,p,a,b,c,d,x,y,niter,time,ii,jj,kk,N,T_temp,n_eq,n_stat 
    integer, allocatable :: spin(:,:,:) 
    real*8 :: r,E,M,mag,Ei,Ef,dE,u
    real*8 :: av_m,av_e,avg_M,avg_E,av_M2,av_E2,av_abs_m,avg_abs_M,cv,chi,chi_abs,ul,av_M4  
    real*8 :: kbT,Ji=1.0d0 
    !Q1-5
    print*, 'L : '
    read*, L 
    print*, 'niter : '
    read*, niter
    allocate(spin(L,L,L))
    E = 0.0d0
    M = 0.0d0 
    N = L*L*L 
! Initializing the lattice
    p = 0.0d0 
    do k = 1,L
        do j = 1,L
            do i = 1,L
                spin(i,j,k) = 1
            end do
        end do
    end do 
! Initializing energy and magnetization
    do k =1,L 
        do j = 1,L 
            do i = 1,L 
                a = i+1; b = i-1; c = j+1; d = j-1; x = k+1; y = k-1
                if (i == L) a = 1; if (i == 1) b = L 
                if (j == L) c = 1; if (j == 1) d = L 
                if (k == L) x = 1; if (k == 1) y = L 

                E = E - Ji*(spin(i,j,k)*(spin(a,j,k)+spin(b,j,k)+spin(i,c,k)+spin(i,d,k)+spin(i,j,x)+spin(i,j,y)))
                M = M + spin(i,j,k)
            end do
        end do 
    end do 
    mag = M/(dfloat(N))
    E = E*0.5d0
    print*, 'Initial Energy: ', E 
    print*, 'Initial Energy per spin: ', E/(dfloat(N))
    print*, 'Initial magnetization: ', M 
    print*, 'Initial magnetization per spin: ', mag
! Time evolution of the system 
    open(1, file='3d_ising_val.txt', status='new')
    open(2, file='3d_ising_kB_4,9.txt', status='new')
    open(3, file='3d_ising_kB_3,9.txt', status='new')
    open(4, file='3d_ising_kB_4,05.txt', status='new')
    open(7, file='3d_ising_L8.txt', status='new')
    open(8, file='3d_ising_L9.txt', status='new')
    open(9, file='3d_ising_L10.txt', status='new')
    do time = 1,niter
        do kk = 1,L 
            do jj = 1,L 
                do ii = 1,L 
                    call random_number(r); i = int(r*dfloat(L)) + 1
                    call random_number(r); j = int(r*dfloat(L)) + 1
                    call random_number(r); k = int(r*dfloat(L)) + 1

                    a = i+1; b = i-1; c = j+1; d = j-1; x = k+1; y = k-1
                    if (i == L) a = 1; if (i == 1) b = L 
                    if (j == L) c = 1; if (j == 1) d = L 
                    if (k == L) x = 1; if (k == 1) y = L 

                    Ei = -Ji*(spin(i,j,k)*(spin(a,j,k)+spin(b,j,k)+spin(i,c,k)+spin(i,d,k)+spin(i,j,x)+spin(i,j,y)))
                    spin(i,j,k) = -spin(i,j,k)
                    Ef = -Ji*(spin(i,j,k)*(spin(a,j,k)+spin(b,j,k)+spin(i,c,k)+spin(i,d,k)+spin(i,j,x)+spin(i,j,y)))

                    dE = Ef - Ei 
                    if (dE .le. 0) then
                        E = E + dE 
                        M = M + 2.0d0*dfloat(spin(i,j,k))
                    else 
                        u = exp(-dE/(kbT))
                        call random_number(r)
                        if (r < u) then 
                            E = E + dE 
                            M = M + 2.0d0*dfloat(spin(i,j,k))
                        else 
                            spin(i,j,k) = -spin(i,j,k)!Trial flip not accepted
                        end if 
                    end if 
                end do 
            end do 
        end do 
        write(1,*) time, M/dfloat(N), E/dfloat(N)
        write(2,*) M/dfloat(N)
        write(3,*) E/dfloat(N)
        write(4,*) M/dfloat(N), E/dfloat(N)
        write(L-1,*) M/dfloat(N), E/dfloat(N)
    end do 
    close(1)
    close(2)
    close(3)
    close(4)
    close(7)
    close(8)
    close(9)
!---------------------------------------------------------------------------------------------------------------   
! Q7-10 
    print*, 'q1-5 done.'
    open(17, file='3d_ising_all_L7.txt', status='new')
    open(18, file='3d_ising_all_L8.txt', status='new')
    open(19, file='3d_ising_all_L9.txt', status='new')
    n_eq = 10000
    n_stat = 10
    do L=7,9
        E = 0.0d0
        M = 0.0d0 
        N = L*L*L
        allocate(spin(L,L,L))
        p = 0.0d0 
        do k = 1,L
            do j = 1,L
                do i = 1,L
                    spin(i,j,k) = 1
                end do
            end do
        end do 
        do k =1,L 
            do j = 1,L 
                do i = 1,L 
                    a = i+1; b = i-1; c = j+1; d = j-1; x = k+1; y = k-1
                    if (i == L) a = 1; if (i == 1) b = L 
                    if (j == L) c = 1; if (j == 1) d = L 
                    if (k == L) x = 1; if (k == 1) y = L 

                    E = E - Ji*(spin(i,j,k)*(spin(a,j,k)+spin(b,j,k)+spin(i,c,k)+spin(i,d,k)+spin(i,j,x)+spin(i,j,y)))
                    M = M + spin(i,j,k)
                end do
            end do 
        end do 
        mag = M/(dfloat(N))
        E = E*0.5d0   
        do T_temp=470,380,-2
            kbT = dfloat(T_temp)/100.0d0
            av_abs_m=0.0d0
            av_m=0.0d0;av_e=0.0d0
            avg_abs_M = 0.0d0
            avg_M = 0.0d0; avg_E = 0.0d0
            av_M2 = 0.0d0; av_E2 = 0.0d0 
            av_M4 = 0.0d0 
            chi_abs = 0.0d0 
            do time = 1,niter
                do kk = 1,L 
                    do jj = 1,L 
                        do ii = 1,L 
                            call random_number(r); i = int(r*dfloat(L)) + 1
                            call random_number(r); j = int(r*dfloat(L)) + 1
                            call random_number(r); k = int(r*dfloat(L)) + 1
        
                            a = i+1; b = i-1; c = j+1; d = j-1; x = k+1; y = k-1
                            if (i == L) a = 1; if (i == 1) b = L 
                            if (j == L) c = 1; if (j == 1) d = L 
                            if (k == L) x = 1; if (k == 1) y = L 
        
                            Ei = -Ji*(spin(i,j,k)*(spin(a,j,k)+spin(b,j,k)+spin(i,c,k)+spin(i,d,k)+spin(i,j,x)+spin(i,j,y)))
                            spin(i,j,k) = -spin(i,j,k)
                            Ef = -Ji*(spin(i,j,k)*(spin(a,j,k)+spin(b,j,k)+spin(i,c,k)+spin(i,d,k)+spin(i,j,x)+spin(i,j,y)))
        
                            dE = Ef - Ei 
                            if (dE .le. 0) then
                                E = E + dE 
                                M = M + 2.0d0*dfloat(spin(i,j,k))
                            else 
                                u = exp(-dE/(kbT))
                                call random_number(r)
                                if (r < u) then 
                                    E = E + dE 
                                    M = M + 2.0d0*dfloat(spin(i,j,k))
                                else 
                                    spin(i,j,k) = -spin(i,j,k)!Trial flip not accepted
                                end if 
                            end if 
                        end do 
                    end do 
                end do   
                if (time .ge. n_eq)then 
                    av_abs_m = av_abs_m + abs(M)/dfloat(N)
                    av_m = av_m + M/dfloat(N); av_e = av_e + E/dfloat(N)
                    avg_abs_M = avg_abs_M + abs(M)
                    avg_M = avg_M + M; avg_E = avg_E + E 
                    av_M2 = av_M2 + M**2; av_E2 = av_E2 + E**2 
                    av_M4 = av_M4 + M**4 
                end if 
            end do !ending time 
            av_abs_m = av_abs_m/dfloat(niter-n_eq)
            av_m = av_m/dfloat(niter-n_eq); av_e = av_e/dfloat(niter-n_eq)
            avg_abs_M = avg_abs_M/dfloat(niter-n_eq)
            avg_M = avg_M/dfloat(niter-n_eq); avg_E = avg_E/dfloat(niter-n_eq)
            av_M2 = av_M2/dfloat(niter-n_eq);av_E2 = av_E2/dfloat(niter-n_eq)
            av_M4 = av_M4/dfloat(niter-n_eq)
            cv = (av_E2 - avg_E**2)/(kbT**2)
            chi_abs = (av_M2 - avg_abs_M**2)/(kbT)
            chi = (av_M2 - avg_M**2)/(kbT)
            ul = 1.0d0 - av_M4/(3*av_M2**2)
            write(L+10,*) kbT,',',av_abs_m,',',av_m,',',av_e,',',chi_abs,',',chi,',',cv,',',ul 
            
        end do !ending T_temp
        close(L+10)
        deallocate(spin)
    end do !ending L 
end program fortrine 