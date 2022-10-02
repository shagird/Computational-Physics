program fortrine
    implicit none
    integer, parameter :: lx=100, ly=100
    real*8, dimension(0:lx-1,0:ly-1) :: A, oldA, B, oldB
    integer :: i,j,k,ii,jj
    real*8 :: bd_temp, inc_temp, pf, rhoA, rhoB
    integer :: xp,xn,yp,yn

    real*8,parameter :: dx=1.0d0,dy=1.0d0,dt=0.002d0
    real*8 :: idx2=1.0d0/dx**2,idy2=1.0d0/dy**2
    real*8 :: toler= 0.0001d0
    real*8 :: DA=1.0d0,DB=100.d0,alpha=-0.005d0,beta=10.0d0
    integer:: nsnap=1000, niter=40000, t
    real*8 :: rd

    A=0.0d0; B=0.0d0; oldA=0.0d0; oldB=0.0d0
    rhoA=0.0d0; rhoB=0.0d0
    ! do i=1,nint(0.30d0*100.0d0*dy)
    !     A(nint(1.3*i),ly/2) = 3.0d0
    ! end do
    A(:,ly/2) = 0.30d0
    ! A(lx/2,ly/2) = 0.1d0
    oldA = A
    open(11,file='A.txt'); open(12,file='B.txt')
    do t=1,niter
        do jj=0,ly-1
            xp=modulo(jj+1,ly); xn=modulo(jj-1,ly)
            do ii=0,lx-1
                yp=modulo(ii+1,lx); yn=modulo(ii-1,lx)
                A(ii,jj) = oldA(ii,jj) + (dt*DA*idx2)*(oldA(yn,jj) + oldA(yp,jj) + oldA(ii,xn) + oldA(ii,xp) - 4.0d0*oldA(ii,jj))
                A(ii,jj) = A(ii,jj) + dt*(A(ii,jj) - A(ii,jj)**3 + alpha - oldB(ii,jj))

                B(ii,jj) = oldB(ii,jj) + (dt*DB*idx2)*(oldB(yn,jj) + oldB(yp,jj) + oldB(ii,xn) + oldB(ii,xp) - 4.0d0*oldB(ii,jj))
                B(ii,jj) = B(ii,jj) + dt*beta*(oldA(ii,jj) - oldB(ii,jj))
            end do
        end do
        oldA = A 
        oldB = B
    end do
    do ii=0,lx-1
        do jj=0,ly-1
            write(11,*) ii,' ',jj,',',A(ii,jj)
            write(12,*) ii,' ',jj,',',B(ii,jj)
            rhoA = rhoA + A(ii,jj)
            rhoB = rhoB + B(ii,jj)
        end do 
    end do
    close(11); close(12)
    print*, rhoA, rhoB
end program fortrine