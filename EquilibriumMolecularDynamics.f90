module params
    implicit none
    integer, parameter :: npart=3600, niter=5000, ncalc_avg = 10
    integer, parameter :: lx=20, ly=20, lz=20 
    real*8, parameter :: m=1.0d0, kbT=1.0d0, sigma=1.0d0, rc=2.50d0, rs = rc + 2.0d0, rg = 10.0d0

    real*8, parameter :: llx = dfloat(lx), lly = dfloat(ly), llz = dfloat(lz), dnpart=dfloat(npart)
    real*8 :: llxb2 = llx/2.0d0, llyb2 = lly/2.0d0, llzb2 = llz/2.0d0
    real*8, dimension(3*npart) :: pos, vel, force, old_force
    integer, dimension(npart) :: neb_no
    real*8, dimension(100) :: neb_no_g, neb_no_g_avg, gr_array_time_avg=0.0d0
    integer, dimension(npart,900) :: neb_list

    real*8, parameter :: sigma6=sigma**6, sigma12=sigma**12, eps=1.0d0, pi =3.141592653589793 
    real*8, parameter :: fc = 4.0d0*eps*(12.0d0*sigma12/(rc**13) - 6.0d0*sigma6/(rc**7))
    real*8, parameter :: uc = fc*rc + 4.0d0*eps*((sigma/rc)**12 - (sigma/rc)**6)
    real*8, parameter :: dr = 0.1d0, dt = 0.0025d0, dtb2 = 0.5d0*dt, dt2b2 = 0.5d0*dt**2.0d0, bm = 1.0d0/m
    real*8, parameter :: vel_const = dsqrt(12.0d0)*kbT, kE_act = 1.50d0*npart*kbT
    real*8, parameter :: rho = dnpart/(llx*lly*llz)

    integer :: i,j,k,c,t,ctr,ind_r 
    real*8 :: rd, avg_vx, avg_vy, avg_vz
    real*8 :: vavg_x, vavg_y, vavg_z
    real*8 :: x,y,z,x1,y1,z1,x2,y2,z2,dx,dy,dz,r,LJ,pE,kE,gy,gz,gx
    real*8 :: LJ_force 
    real*8 :: s, s1, temp_r, vol_lcl
end module params

subroutine thermostat
    use params
    implicit none
    s = kE_act/kE
    s1 = dsqrt(s)
    vel = s1*vel
    kE = 0.5d0*m*sum(vel**2)
end subroutine thermostat

program fortrine
    use params
    implicit none
    open(11, file = '14_PK.txt')
    open(12, file = '14_Moment.txt')
    open(17, file='1_PCF.txt')
    open(18, file='4a.txt')
    open(19, file='1_Velocities.txt')
! --------------------------------------------------------------------------------------------------------
    !q2
    !initializing position 
! --------------------------------------------------------------------------------------------------------

    pos = 0.0d0 
    ctr = 1
    do i=1,lx
        do j=1,ly
            do k=1,lz 
                call random_number(rd)
                if (rd<3800.0d0/8000.0d0) then 
                    if (ctr.le.npart) then 
                        pos(3*ctr-2) = i
                        pos(3*ctr-1) = j
                        pos(3*ctr) = k 
                        ctr = ctr + 1 
                    end if 
                end if 
            end do 
        end do
    end do 

    ! c = 0
    ! do i = 1,10
    !     do j = 1,10
    !         do k = 1, 12
    !             pos(c+1) = 1.5*i
    !             pos(c+2) = 1.5*j
    !             pos(c+3) = 1.5*k
    !             c = c+3
    !         end do
    !     end do
    ! end do
    ! c = 0
    
    !initializing velocity 
    vel = 0.0d0
    kE = 0.0d0 
    
    do i=1,npart
        call random_number(rd)
        vel(3*i-2) = vel_const*(rd-0.5d0)
        call random_number(rd)
        vel(3*i-1) = vel_const*(rd-0.5d0)
        call random_number(rd)
        vel(3*i) = vel_const*(rd-0.5d0)
    end do
    avg_vx = 0.0d0; avg_vy = 0.0d0; avg_vz = 0.0d0 
    do i=1,npart
        avg_vx = vel(3*i-2) + avg_vx
        avg_vy = vel(3*i-1) + avg_vy
        avg_vz = vel(3*i) + avg_vz
    end do
    avg_vx = avg_vx/dnpart
    avg_vy = avg_vy/dnpart
    avg_vz = avg_vz/dnpart

    do i=1,npart
        vel(3*i-2) = vel(3*i-2) - avg_vx
        vel(3*i-1) = vel(3*i-1) - avg_vy
        vel(3*i) =  vel(3*i) - avg_vz 
    end do 
    kE = 0.5d0*m*sum(vel**2)
    !Initial force calculation
    old_force = 0.0d0; force = 0.0d0; pE = 0.0d0 
    do i=1,npart-1 
        x1 = pos(3*i-2); y1 = pos(3*i-1); z1 = pos(3*i)
        do j=i+1,npart
            x2 = pos(3*j-2); y2 = pos(3*j-1); z2 = pos(3*j)
            
            x = x1-x2; y = y1-y2; z = z1-z2 
            if (abs(x).ge.llxb2) x = (llx - abs(x))*((-1.0d0*x)/abs(x))
            if (abs(y).ge.llyb2) y = (lly - abs(y))*((-1.0d0*y)/abs(y))
            if (abs(z).ge.llzb2) z = (llz - abs(z))*((-1.0d0*z)/abs(z)) 
            
            r = dsqrt(x**2 + y**2 + z**2)
            if (r .le. rc) then 
               LJ = 4.0d0*eps*((sigma/r)**12 - (sigma/r)**6) - uc + fc*r 
               pE = pE + LJ 
               LJ_force = 4.0d0*eps*(12.0d0*sigma12/(r**13) - 6.0d0*sigma6/(r**7)) - fc 

               force(3*i-2) = force(3*i-2) + LJ_force*(x/r)
               force(3*i-1) = force(3*i-1) + LJ_force*(y/r)
               force(3*i)   = force(3*i)   + LJ_force*(z/r)
               force(3*j-2) = force(3*j-2) - LJ_force*(x/r)
               force(3*j-1) = force(3*j-1) - LJ_force*(y/r)
               force(3*j)   = force(3*j)   - LJ_force*(z/r)
            end if 
        end do 
    end do
    
    !implementing initial neighbour list
    neb_no = 0; neb_list = 0
    do i=1,npart-1 
        x1 = pos(3*i-2); y1 = pos(3*i-1); z1 = pos(3*i)
        do j=i+1,npart
            x2 = pos(3*j-2); y2 = pos(3*j-1); z2 = pos(3*j)
            
            x = x1-x2; y = y1-y2; z = z1-z2 
            if (abs(x).ge.llxb2) x = (llx - abs(x))*((-1.0d0*x)/abs(x))
            if (abs(y).ge.llyb2) y = (lly - abs(y))*((-1.0d0*y)/abs(y))
            if (abs(z).ge.llzb2) z = (llz - abs(z))*((-1.0d0*z)/abs(z)) 
            
            r = dsqrt(x**2 + y**2 + z**2)
            if (r .le. rs) then 
                neb_no(i) = neb_no(i) + 1
                neb_list(i,neb_no(i)) = j
            end if 
        end do 
    end do 

    !updating whole dynamics    
! --------------------------------------------------------------------------------------------------------

    
    do t=1,niter
        !updating position 
        do i=1,npart
            pos(3*i-2) = pos(3*i-2) + vel(3*i-2)*dt + dt2b2*force(3*i-2)*bm
            pos(3*i-1) = pos(3*i-1) + vel(3*i-1)*dt + dt2b2*force(3*i-1)*bm
            pos(3*i) = pos(3*i) + vel(3*i)*dt + dt2b2*force(3*i)*bm
        end do
        ! pos = pos + vel*dt + force*dt2b2
        pos(3*i-2) = modulo(pos(3*i-2), llx)
        pos(3*i-1) = modulo(pos(3*i-1), lly)
        pos(3*i) = modulo(pos(3*i), llz)

        old_force = force 

        !calculate new forces using neighbour list
        pE = 0.0d0; force = 0.0d0 
        do i=1,npart 
            x1 = pos(3*i-2); y1 = pos(3*i-1); z1 = pos(3*i)
            do k=1,neb_no(i)
                j = neb_list(i,k)
                x2 = pos(3*j-2); y2 = pos(3*j-1); z2 = pos(3*j)
                x = x1-x2; y = y1-y2; z = z1-z2 

                if (abs(x).ge.llxb2) x = (llx - abs(x))*((-1.0d0*x)/abs(x))
                if (abs(y).ge.llyb2) y = (lly - abs(y))*((-1.0d0*y)/abs(y))
                if (abs(z).ge.llzb2) z = (llz - abs(z))*((-1.0d0*z)/abs(z)) 
                
                r = dsqrt(x**2 + y**2 + z**2)
                if (r .le. rc) then 
                    LJ = 4.0d0*eps*((sigma/r)**12 - (sigma/r)**6) - uc + fc*r 
                    pE = pE + LJ 
                    LJ_force = 4.0d0*eps*(12.0d0*sigma12/(r**13) - 6.0d0*sigma6/(r**7)) - fc 
    
                    force(3*i-2) = force(3*i-2) + LJ_force*(x/r)
                    force(3*i-1) = force(3*i-1) + LJ_force*(y/r)
                    force(3*i) = force(3*i) + LJ_force*(z/r)
                    force(3*j-2) = force(3*j-2) - LJ_force*(x/r)
                    force(3*j-1) = force(3*j-1) - LJ_force*(y/r)
                    force(3*j) = force(3*j) - LJ_force*(z/r)
                end if 
            end do 
        end do 

        !updating velocities 
        vavg_x = 0.0d0; vavg_y = 0.0d0; vavg_z = 0.0d0 
        
        kE = 0.0d0 
        do i=1,npart 
            vel(3*i-2) = vel(3*i-2) + dtb2*(force(3*i-2) + old_force(3*i-2))*bm
            vel(3*i-1) = vel(3*i-1) + dtb2*(force(3*i-1) + old_force(3*i-1))*bm
            vel(3*i) = vel(3*i) + dtb2*(force(3*i) + old_force(3*i))*bm
        end do 

        kE = 0.5d0*sum(vel**2)

        !applying thermostat
        if (mod(t,100).eq.0) call thermostat 

        !applying neigbour list
        if (mod(t,40).eq.0) then 
            neb_no = 0; neb_list = 0
            do i=1,npart-1 
                x1 = pos(3*i-2); y1 = pos(3*i-1); z1 = pos(3*i)
                do j=i+1,npart
                    x2 = pos(3*j-2); y2 = pos(3*j-1); z2 = pos(3*j)
                    
                    x = x1-x2; y = y1-y2; z = z1-z2 
                    if (abs(x).ge.llxb2) x = (llx - abs(x))*((-1.0d0*x)/abs(x))
                    if (abs(y).ge.llyb2) y = (lly - abs(y))*((-1.0d0*y)/abs(y))
                    if (abs(z).ge.llzb2) z = (llz - abs(z))*((-1.0d0*z)/abs(z)) 
                    
                    r = dsqrt(x**2 + y**2 + z**2)
                    if (r .le. rs) then 
                        neb_no(i) = neb_no(i) + 1
                        neb_list(i,neb_no(i)) = j
                    end if 
                end do 
            end do 
            ! !q4a
            ! if ((t.ge.20000).and.(t.le.30000)) then
            !     write(18,*) maxval(neb_no)
            ! end if
        end if
        
        !calcultating thermodynamic quantities
        if (t > 0) then 
            if (mod(t,100).eq.0) then !calculating pair correlation function
                ctr = 0
                neb_no_g_avg = 0.d0
                do i=1,npart
                x1 = pos(3*i-2); y1 = pos(3*i-1); z1 = pos(3*i)
                neb_no_g = 0.0d0
                    do j=1,npart
                        if (i.ne.j) then
                            x2 = pos(3*j-2); y2 = pos(3*j-1); z2 = pos(3*j)
                            
                            x = x1-x2; y = y1-y2; z = z1-z2 
                            if (abs(x).ge.llxb2) x = (llx - abs(x))*((-1.0d0*x)/abs(x))
                            if (abs(y).ge.llyb2) y = (lly - abs(y))*((-1.0d0*y)/abs(y))
                            if (abs(z).ge.llzb2) z = (llz - abs(z))*((-1.0d0*z)/abs(z)) 
                            
                            r = dsqrt(x**2 + y**2 + z**2)
                            if (r .le. rg) then
                                ctr = ctr + 1  
                                ind_r = ceiling(10*r)
                                neb_no_g(ind_r) = neb_no_g(ind_r) + 1.0d0
                            end if 
                        end if
                    end do 
                    neb_no_g_avg = neb_no_g_avg + neb_no_g
                end do
            end if
        !     if (mod(t,50).eq.0) then 
        !         do i=1,size(vel)                 
        !             write(19,*) vel(i)
        !         end do
        !     end if
        end if 
        gr_array_time_avg = gr_array_time_avg + neb_no_g_avg
        !writing data 
        write(11,*) pE, ',', kE
        write(12,*) m*sum(vel)
    end do
    gr_array_time_avg = gr_array_time_avg * 100 /(niter)
    do i=1,size(gr_array_time_avg)
        write(17,*) gr_array_time_avg(i)
    end do
    close(11); close(12); close(17); close(18); close(19)
    print*, 'q2 done.'
end program fortrine