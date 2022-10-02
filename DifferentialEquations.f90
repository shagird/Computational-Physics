program fortrine 
implicit none
integer:: n,i,ii,nop=50,mm,cond 
real*8:: dy,dx=0.001d0,x=0.0d0,y=0.0d0,y1,dy1,y_act = 48.078d0,f0,f1,f2,f3,yt1,yt2,yt3
real*8:: v,dt,f0x,xt1,f1x,xt2,f2x,xt3,f3x,f0v,vt1,f1v,vt2,f2v,vt3,f3v
real*8, dimension(0:49):: vp,v1p,v2p,v3p,yp,y1p,y2p,y3p,f0yp,f1yp,f2yp,f3yp,f0vp,f1vp,f2vp,f3vp
real*8,parameter :: fx = 0.0d0, endx= 1.0d0,dx1=0.01d0,lt=0.0001d0
integer, parameter :: ngrid = nint((endx -fx)/dx1) + 1
real*8, parameter::d1=1.0d0/(2.0d0-10.00d0*dx1*dx1),d2=(1.d0-2.5d0*dx1),d3=(1.0d0+2.5d0*dx1),d4=-10.0d0*dx1*dx1 
real*8, dimension(ngrid) :: x1,y2,y2old
n = nint(1.550d0/dx)
open(21,file='q1.txt')
do i=1,n
    dy = (y*y + 1.0d0)*dx
    y = y+dy
    write(21,*) y 
end do 
close(21)
open(22,file='q2.txt')
print'(a21,f8.4)', 'Euler diff: ', y_act - y 
y = 0.0d0
do i=1,n
    dy = (y*y + 1.0d0)*dx/2.0d0
    y1 = y+dy
    dy = (y1*y1 + 1)*dx
    y = y+dy
    write(22,*) y
end do 
close(22)
open(23,file='q3.txt')
print'(a21,f8.4)', 'Modified Euler diff: ', y_act - y
y = 0.0d0
do i = 1,n
    dy = (y*y + 1.0d0)*dx
    y1 = y + dy
    dy1 = (y1*y1 + 1.0d0)*dx
    y = y + (dy+dy1)/2.0d0 
    write(23,*) y
end do 
close(23)
open(24,file='q4.txt')
print'(a21,f8.4)', 'Improved Euler diff: ', y_act - y 
dx = 0.01d0
y = 0.0d0
n = nint(1.550d0/dx)
do i = 1,n
    f0 = (y*y + 1.0d0)
    yt1 = y + (f0)*(dx/2.0d0)
    f1 = yt1*yt1 + 1.0d0
    yt2 = y + (f1)*(dx/2.0d0)
    f2 = yt2*yt2 + 1.0d0 
    yt3 = y + f2*dx 
    f3 = yt3*yt3 + 1.0d0 
    y = y + (dx/6.0d0)*(f0 + 2.0d0*f1 + 2.0d0*f2 + f3)
    write(24,*) y
end do 
close(24)
print'(a21,f8.4)', 'RK4 diff: ', y_act - y
print*, 'q1-4 done.'
x = 0.1d0
v = 1.90d0 
dt = 0.010d0
n = 5000
open(11,file='q5.txt')
do i =1,n
    f0x = v 
    f0v = -sin(x)
    xt1 = x + (f0x)*(dt/2.0d0)
    vt1 = v + (f0v)*(dt/2.0d0)
    f1x = vt1 
    f1v = -sin(xt1) 
    xt2 = x + (f1x)*(dt/2.0d0)
    vt2 = v + (f1v)*(dt/2.0d0)
    f2x = vt2 
    f2v = -sin(xt2)
    xt3 = x + f2x*dt  
    vt3 = v + f2v*dt
    f3x = vt3 
    f3v = -sin(xt3)
    x = x + (dt/6.0d0)*(f0x + 2.0d0*f1x + 2.0d0*f2x + f3x)    
    v = v + (dt/6.0d0)*(f0v + 2.0d0*f1v + 2.0d0*f2v + f3v)
    write(11,*) x, ',', v, ',', v*v/2.0d0 - cos(x)
end do
close(11) 
print*,'x at t=50: ', x
print*, 'q5 done.'
x = 0.0d0
v = 1.9990d0 
dt = 0.010d0
n = 5000
open(12,file='q6.txt')
do i =1,n
    f0x = v 
    f0v = -sin(x)
    xt1 = x + (f0x)*(dt/2.0d0)
    vt1 = v + (f0v)*(dt/2.0d0)
    f1x = vt1 
    f1v = -sin(xt1) 
    xt2 = x + (f1x)*(dt/2.0d0)
    vt2 = v + (f1v)*(dt/2.0d0)
    f2x = vt2 
    f2v = -sin(xt2)
    xt3 = x + f2x*dt  
    vt3 = v + f2v*dt
    f3x = vt3 
    f3v = -sin(xt3)
    x = x + (dt/6.0d0)*(f0x + 2.0d0*f1x + 2.0d0*f2x + f3x)    
    v = v + (dt/6.0d0)*(f0v + 2.0d0*f1v + 2.0d0*f2v + f3v)
    write(12,*) x, ',', v, ',', v*v/2.0d0 - cos(x)
end do 
close(12)
print*, 'x at t=50: ', x
print*, 'q6 done.'
x = 0.0d0
v = 2.9990d0 
dt = 0.010d0
n = 5000
open(13,file='q7.txt')
do i =1,n
    f0x = v 
    f0v = -sin(x)
    xt1 = x + (f0x)*(dt/2.0d0)
    vt1 = v + (f0v)*(dt/2.0d0)
    f1x = vt1 
    f1v = -sin(xt1) 
    xt2 = x + (f1x)*(dt/2.0d0)
    vt2 = v + (f1v)*(dt/2.0d0)
    f2x = vt2 
    f2v = -sin(xt2)
    xt3 = x + f2x*dt  
    vt3 = v + f2v*dt
    f3x = vt3 
    f3v = -sin(xt3)
    x = x + (dt/6.0d0)*(f0x + 2.0d0*f1x + 2.0d0*f2x + f3x)    
    v = v + (dt/6.0d0)*(f0v + 2.0d0*f1v + 2.0d0*f2v + f3v)
    write(13,*) x, ',', v, ',', v*v/2.0d0 - cos(x)
end do 
close(13)
print*, 'q7 done.'
dt = 0.02d0
vp = 0.0d0
yp = 0.0d0
yp(0) = 0.80d0
yp(25) = 0.80d0
n=2000
open(14,file='q8.txt')
do i=1,n
    do ii=0,nop-1 
        f0yp(ii) = vp(ii)
        f0vp(ii) = (yp(modulo(ii+1,nop)) + yp(modulo(ii-1,nop)) - 2.0d0*yp(ii))
    end do 
    do ii=0,nop-1
        y1p(ii) = yp(ii) + f0yp(ii)*(dt/2.0d0)
        v1p(ii) = vp(ii) + f0vp(ii)*(dt/2.0d0)
    end do 
    do ii=0,nop-1
        f1yp(ii) = v1p(ii)
        f1vp(ii) = (y1p(modulo(ii+1,nop))+y1p(modulo(ii-1,nop)) - 2.0d0*y1p(ii))
    end do 
    
    do ii=0,nop-1
        y2p(ii) = yp(ii) + f1yp(ii)*(dt/2.0d0)
        v2p(ii) = vp(ii) + f1vp(ii)*(dt/2.0d0)
    end do 
    
    do ii=0,nop-1
        f2yp(ii) = v2p(ii)
        f2vp(ii) = (y2p(modulo(ii+1,nop))+y2p(modulo(ii-1,nop)) - 2.0d0*y2p(ii))
    end do 
    do ii=0,nop-1
        y3p(ii) = yp(ii) + f2yp(ii)*dt 
        v3p(ii) = vp(ii) + f2vp(ii)*dt 
    end do 
    do ii=0,nop-1
        f3yp(ii) = v3p(ii)
        f3vp(ii) = (y3p(modulo(ii+1,nop))+y3p(modulo(ii-1,nop)) - 2.0d0*y3p(ii)) 
    end do 
    do ii=0,nop-1
        yp(ii) = yp(ii) + (dt/6.0d0)*(f0yp(ii) + 2.0d0*f1yp(ii) + 2.0d0*f2yp(ii) + f3yp(ii))
        vp(ii) = vp(ii) + (dt/6.0d0)*(f0vp(ii) + 2.0d0*f1vp(ii) + 2.0d0*f2vp(ii) + f3vp(ii))
    end do 
    write(14,*) yp(0), ',', vp(0)
end do
close(14) 
print*, 'position of 1st particle at t=40: ', yp(0)
print*, 'q8 done.'
open(15, file='q9.txt')
mm = 0
cond = 0
x1(1) = fx; x1(ngrid) = endx; y2(1) = 0.0d0; y2(ngrid) = 2.0d0
do i=2,ngrid-1
    x1(i) = x1(i-1) + dx1
    y2(i) = (y2(ngrid) - y2(1))*x1(i)/(endx - fx)
end do 
do 
    mm = mm + 1
    if(cond == 1) exit
    y2old = y2 
    do i=2,ngrid-1
        y2(i) = d1*(d2*y2old(i+1) + d3*y2(i-1) + d4*x1(i))
    end do
    cond = 1
    do i=2,ngrid-1
        if (abs(y2old(i) - y2(i)) .ge. lt) cond=0
    end do 
end do 
print*, 'y at x = 0.8: ', y2(81)
do i=1,ngrid
    write(15,*) y2(i)
end do
close(15)
end program fortrine
