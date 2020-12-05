! this a test file of 
!            hwscrt
! Standard five-point finite difference approximation to the Helmoholtz 
! Equation in cartesian coordinates using a centered finite difference grid.
! The equation has the form 
!            u_xx + u_yy + lambda u = f
! 
! Since the orginal file is only single precision, I changed all definition
! of real to real(8), then the result is double precision.
!
! author : Qinghe Wang
! data   : 2020/12/05
! email  : qinghev@gmail.com
    
 program main
    implicit none
    real(8) :: a,b,c,d
    integer :: xbd,ybd,nx, i,num
    real(8),allocatable :: error(:,:), order(:,:),enorm(:)
    integer,allocatable :: gn(:)
    
    num = 7
    allocate(error(num+1,2),order(num,2),enorm(2),gn(num))
    
    a = 0d0
    b = 1d0
    c = 0d0
    d = 1d0
    xbd = 1
    ybd = 1
    nx  = 10
    
    do i = 1,num
        gn(i) = nx*2**(i-1)
        call mytest_hwscrt(a,b,c,d,xbd,ybd,gn(i),enorm)
        error(i,:) = enorm
    enddo
    order = 0d0
    do i = 1,num-1
        order(i+1,:) = log(error(i,:)/error(i+1,:))/log(2.)
    enddo
    write(*,*) 'nx    norm2     order    norm_inf     order'
    do i = 1,num
      write(*,'(I5,4x, E9.4,4x, F9.4,4x, E9.4,4x, F9.4)') gn(i),error(i,1), &
          order(i,1),error(i,2),order(i,2) 
    end do
    deallocate(error,order,enorm,gn)
    stop
 end program main
    
    
    
    
    
 subroutine mytest_hwscrt(a,b,c,d,xbd,ybd,nx,enorm)
    implicit none
    !REAL(8)(4),parameter :: pi = 3.1415926535897932
    real(8)           :: a,b,c,d,  pi
    real(8)           :: enorm(2)
    !real(8), external :: log
    real(8)           :: h, lambda, hh, t, pertrb, norm2,norm_inf
    integer           :: nx,ny, i,j, xbd, ybd, s, ierror, nxp1,nyp1
    real(8), allocatable :: bda(:),bdb(:),bdc(:),bdd(:), work(:)
    real(8), allocatable :: f(:,:), err(:,:), u(:,:)
    
    h  = (b-a)/nx
    ny = nint((d-c)/h)
    hh = h*h
    s  = 4*(ny+1) + (13 + int(log((ny+1)*1.)/log(2.)))*(nx+1)
    lambda = 1d0
    nxp1 = nx + 1
    nyp1 = ny + 1 
    
    allocate(u(nxp1,nyp1), f(nxp1,nyp1),bda(nyp1),bdb(nyp1), &
        bdc(nxp1),bdd(nxp1),work(s),err(nxp1,nyp1))
    
    ! ---- exact solution ----
    pi = 4.*atan(1.)
    do i = 1,nxp1
        do j = 1,nyp1
            u(i,j) = exp((i-1)*h)*sin(pi*(j-1)*h)
        enddo
    enddo
    
    ! ---- right hand side ----
    f = (1d0-pi*pi+lambda)*u
    
    ! ---- boundary condition ----
    if (xbd == 1) then
        f(1,:)    = u(1,:)
        f(nxp1,:) = u(nxp1,:)
    else if (xbd == 2) then
        f(1,:)  = u(1,:)
        do i = 1,nyp1
            f(nxp1,i) = exp(b)*sin(pi*(i-1)*h)
        enddo
        bdb = f(nxp1,:)
    else if (xbd == 3) then
        do i = 1,nyp1
            f(1,i)    = exp(a)*sin(pi*(i-1)*h)
            f(nxp1,i) = exp(b)*sin(pi*(i-1)*h)
        enddo
        bda = f(1,:)
        bdb = f(nxp1,:)
    else if (xbd == 4) then
        do i = 1,nyp1
            f(1,i)    = exp(a)*sin(pi*(i-1)*h)
        enddo
        f(nxp1,:) = u(nxp1,:)
        bda = f(1,:)
    endif
    
    if (ybd == 1) then
        f(:,1)    = u(:,1)
        f(:,nyp1) = u(:,nyp1)
    else if (ybd == 2) then
        f(:,1) = u(:,1)
        do i = 1,nxp1
            f(i,nyp1) = pi*exp((i-1)*h)*cos(pi*d)
        enddo
        bdd = f(:,nyp1)
    else if (ybd == 3) then
        do i = 1,nxp1
            f(i,1)    = pi*exp((i-1)*h)*cos(pi*c)
            f(i,nyp1) = pi*exp((i-1)*h)*cos(pi*d)
        enddo
        bdc = f(:,1)
        bdd = f(:,nyp1)
    else if (ybd == 4) then
        do i = 1,nxp1
            f(i,1)    = pi*exp((i-1)*h)*cos(pi*c)
        enddo
        f(:,nyp1) = u(:,nyp1)
        bdc = f(:,1)
    endif 
    
    call hwscrt(a,b,nx,xbd,bda,bdb,c,d,ny,ybd,bdc,bdd,lambda,f,nxp1,pertrb,ierror,work)
    write(*,*) nx, ierror, work(1)
    
    err = abs(f-u)
    enorm(1) = sqrt(sum(sum(err*err,2))/(nxp1*1.)/(nyp1*1.))
    enorm(2) = maxval(maxval(err,2))
    
    deallocate(u,f,bda,bdb,bdc,bdd,work,err)
    return
 end subroutine mytest_hwscrt
    
    