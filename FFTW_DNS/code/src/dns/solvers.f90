module solvers
implicit none
! This module contains all the relevant solvers used in the main dns code
contains
    subroutine gfcn(gfns,a,wrkc,wrk1,bct,bcb,dnvb,dnvt,nb,nt,pb,pt)
    
    !---------------------------------------------------------------------!
    !                                                                     !
    !     This subroutine performs the homogeneous solution of the        !
    !     system:                                                         !
    !                                                                     !
    !      d  d(phi)   d  d(phi)   d  d(phi)      Re                      !
    !      -- ------ + -- ------ + -- ------ - -------- phi = 0           !
    !      dx   dx     dy   dy     dz   dz     theta*dt                   !
    !                                                                     !
    !      d  d(Vm)   d  d(Vm)   d  d(Vm)                                 !
    !      -- ----- + -- ----- + -- ----- =  phi                          !
    !      dx  dx     dy  dy     dz  dz                                   ! 
    !                                                                     !
    !     with BCs:                                                       !
    !                                                                     !
    !         phi(1) = pt, phi(-1) = pb, vm(1) = vm(-1) = 0               !
    !                                                                     !
    !---------------------------------------------------------------------!
    
    ! =================================================================== !
    !                           Declare Modules                           !
    ! =================================================================== !
    use grid_size
    use derivs
    !---------------------------------------------------------------------!
    
    ! =================================================================== !
    !                   Declare all variables explicitly                  !
    ! =================================================================== !
    implicit none
    
    ! Passed variables
    complex, dimension(nyp,nz,nxh) :: wrkc, gfns, wrk1
    complex, dimension(nz,nxh)     :: bct,bcb,dnvb,dnvt
    real,    dimension(nyhp,nz,3)  :: a
    real                           :: pb,pt
    integer                        :: nb,nt
    
    ! Calculation variables
    real    :: g,x
    integer :: i,j,k,ib,jp
    
    ! Solver variables
    real :: gain,ugain,theta,alpha,beta,dyde
    real :: wavx(nxh),wavz(nz) 
    
    ! Flow setup variables
    real    :: re,Uinf,R_tau,dPdx
    
    ! Simulation variables
    integer :: it
    real    :: dt
    
    !---------------------------------------------------------------------!
    
    ! =================================================================== !
    !                            Common Blocks                            !
    ! =================================================================== !
    common/solver/     gain,ugain,theta,alpha,beta,dyde
    common/waves/      wavx,wavz
    common/flow/       re,Uinf,R_tau,dPdx
    common/itime/      it
    common/dtime/      dt
    !---------------------------------------------------------------------!
    
    ! =================================================================== !
    !                         Begin Calculations                          !
    ! =================================================================== !
    
    ! Set RHS and BCs for phi equation
    g = re/(theta*dt)
    
    do k = 1,nxh
        do j = 1,nz
            bct(j,k) = cmplx(pt,pt)
            bcb(j,k) = cmplx(pb,pb)
        end do
    end do
    
    do j = nz/2+1,nz
        bct(j,1) = conjg(bct(j,1))
        bcb(j,1) = conjg(bcb(j,1))
    end do
    
    do k = 1,nxh
        do j = 1,nz
            do i = 1,nyp
                gfns(i,j,k) = (0.0,0.0)
            end do
        end do
    end do
    
    bct(1,1) = cmplx(pt,0.0)
    bcb(1,1) = cmplx(pb,0.0)
    
    ! Dirichlet BCs: ib = 0
    ib = 0
    
    ! Solve for phi using Poisson solver (subroutine solve(...))
    ! Solution is returned in variable gfns
    do k = 1,nxh ! debug
        x = wavx(k)**2
        call solve(gfns(1,1,k),x,g,dyde,bct(1,k),bcb(1,k),ib,a,wrkc)
    end do
    
    !---------------------------------------------------------------------!
    ! Set RHS and BCs for Vm equation (RHS = gfns)
    
    g = 0.0
    do k = 1,nxh
        do j = 1,nz
            bct(j,k) = (0.0,0.0)
            bcb(j,k) = (0.0,0.0)
        end do
    end do
    
    do k = 1,nxh
        x = wavx(k)**2
        call solve(gfns(1,1,k),x,g,dyde,bct(1,k),bcb(1,k),ib,a,wrkc)
    end do
    
    do i = 1,nyp
        gfns(i,1,1) = real(gfns(i,1,1))
    end do
    
    
    ! k = 1, j > 1 modes ensure conjugate symmetry
    do j = 2,nz/2
        jp = nz+2-j
        do i = 1,nyp
            gfns(i,j,1) = 0.5*(gfns(i,j,1) + conjg(gfns(i,jp,1)))
        end do
    end do
    
    do j = 2,nz/2
        jp = nz+2-j
        do i = 1,nyp
            gfns(i,jp,1) = conjg(gfns(i,j,1))
        end do
    end do
    
    !---------------------------------------------------------------------!
    ! Calculate the required wall derivatives of gfns and store in dnvt,dnvb
    
    if (nb .eq. 1) then
        call c1derbw(gfns,wrkc,1)
    else if (nb .eq. 2) then
        call cderiv(gfns,wrk1)
        call c1derbw(wrk1,wrkc,1)
    end if
    
    if (nt .eq. 1) then
        call c1dertw(gfns,wrkc,2)
    else if (nt .eq. 2) then
        call cderiv(gfns,wrk1)
        call c1dertw(wrk1,wrkc,2)
    end if
    
    do k = 1,nxh
        do j = 1,nz
            dnvt(j,k) = wrkc(2,j,k)
            dnvb(j,k) = wrkc(1,j,k)
        end do
    end do
    
    end subroutine gfcn
    
    subroutine solve(u,x,g,dyde,bctop,bcbot,ib,a,wrkc)
    
    ! =================================================================== !
    !                           Declare Modules                           !
    ! =================================================================== !
    use grid_size
    !---------------------------------------------------------------------!
    
    ! =================================================================== !
    !                   Declare all variables explicitly                  !
    ! =================================================================== !
    implicit none
    
    ! Passed variables
    complex,dimension(nyp,nz) :: u,wrkc
    complex,dimension(nz) :: bctop,bcbot
    real,   dimension(nyhp,nz,3)  :: a
    real,   dimension(nyp)        :: t
    real    :: x,g,dyde
    integer :: ib
    
    ! Calculation variables
    real :: wn(nz)
    integer :: i,j,k,ip,nrank,kzero
    
    ! Solver variables
    real :: wavx(nxh),wavz(nz) 
    real :: c(nyp)
    
    !---------------------------------------------------------------------!
    
    ! =================================================================== !
    !                            Common Blocks                            !
    ! =================================================================== !
    common/waves/      wavx,wavz,c
    !---------------------------------------------------------------------!
    
    ! =================================================================== !
    !                         Begin Calculations                          !
    ! =================================================================== !
    
    
    if (ib .eq. 0) call dircbc(t,nyp)
    
    do k = 1,nz
        wn(k) = wavz(k)**2 + x + g
    end do
    
    kzero = 0
    if(ib .ne. 0) then ! Neumann BCs
        do k = 1,nz
            if (wn(k) .le. 1.0e-20) then
                kzero = k
                wn(k) = 1.0
            end if
        end do
    end if
    
    !---------------------------------------------------------------------!
    ! Even problem
    ip = 0
    nrank = nyhp
    
    call evnc(a,wn,dyde)
    call evns(u,bctop,bcbot,ib,c,wrkc)
    
    
    if (ib .ne. 0) call neumbc(t,dyde,nrank,ip)
    
    call eqsolve(a,t,wrkc,nrank)
    
    
    ! Update even values of u
    do j = 1,nyp,2
        i = (j+1)/2
        do k = 1,nz
            u(j,k) = wrkc(i,k)
        end do
    end do
    
    !---------------------------------------------------------------------!
    ! Odd problem
    ip = 1
    nrank = nyh
    
    call oddc(a,wn,dyde)
    call odds(u,bctop,bcbot,ib,wrkc)
    
    if (ib .ne. 0) call neumbc(t,dyde,nrank,ip)
    
    call eqsolve(a,t,wrkc,nrank)
    
    ! Update odd values of u
    do j = 2,ny,2
        i = j/2
        do k = 1,nz
            u(j,k) = wrkc(i,k)
        end do
    end do
    
    if ((ib .ne. 0) .and. (kzero .ne. 0)) then
        do j = 1,nyp
            u(j,kzero) = 0.0
        end do
    end if  
    
    !---------------------------------------------------------------------!
    
    
    end subroutine solve
    
    !---------------------------------------------------------------------!
    
    subroutine dircbc(t,n)
    use grid_size
    implicit none
    
    integer :: n,j
    real,dimension(nyp) :: t
    
    do j = 1,n
        t(j)=1.0
    end do
    
    end subroutine dircbc
    
    !---------------------------------------------------------------------!
    
    subroutine neumbc(t,dyde,n,ip)
    use grid_size
    implicit none
    
    integer :: n,ip,j
    real,dimension(nyp) :: t
    real :: dyde
    
    if (ip .eq. 0) then
        do j = 1,n
            t(j) = dyde*(2*(j-1))**2
        end do
    else
        do j = 1,n
            t(j) = dyde*(2*(j-1) + 1)**2
        end do
    end if
    
    end subroutine neumbc
    
    !---------------------------------------------------------------------!
    
    subroutine evnc(a,wn,dyde)
    use grid_size
    implicit none
    
    real, dimension(nyhp,nz,3)  :: a
    real, dimension(nz) :: wn
    real :: dyde,rj1,rj2,rj3
    integer :: i,j,k
    
    do i = 1,nyh
        j = 2*i
        rj1 = 1.0/float(4*j*(j-1))
        rj2 = 1.0/float(2*(j*j-1))
        rj3 = 1.0/float(4*j*(j+1))
        do k = 1,nz
            a(i,k,1) = wn(k)*rj1
            a(i,k,2) = -1.0*dyde*dyde - wn(k)*rj2
            a(i,k,3) = wn(k)*rj3
        end do
    end do
    
    do k = 1,nz
        a(1,k,1) = 2.0*a(1,k,1)
        a(nyhm,k,3) = 0.0
        a(nyh,k,2) = -1.0*dyde*dyde
        a(nyh,k,3) = 0.0
    end do
    
    end subroutine evnc
    
    !---------------------------------------------------------------------!
    
    subroutine oddc(a,wn,dyde)
    use grid_size
    implicit none
    
    real, dimension(nyhp,nz,3)  :: a
    real, dimension(nz) :: wn
    real :: dyde,rj1,rj2,rj3
    integer :: i,j,k
    
    do i = 1,nyhm
        j = 2*i + 1
        rj1 = 1.0/float(4*j*(j-1))  
        rj2 = 1.0/float(2*(j*j-1))  
        rj3 = 1.0/float(4*j*(j+1))  
        do k = 1,nz
            a(i,k,1) = wn(k)*rj1
            a(i,k,2) = -1.0*dyde*dyde - wn(k)*rj2
            a(i,k,3) = wn(k)*rj3
        end do
    end do
    
    do k = 1,nz
        a(nyhm-1,k,3) = 0.0
        a(nyhm,k,2) = -1.0*dyde*dyde
        a(nyhm,k,3) = 0.0
    end do
    
    end subroutine oddc
    
    !---------------------------------------------------------------------!
    
    subroutine evns(s,bctop,bcbot,ib,c,wrk)
    ! Sets up RHS for even tri-diagonal problem, given the results of sbr RHS
    use grid_size
    implicit none
    
    complex, dimension(nyp,nz) :: wrk
    complex, dimension(nz)     :: bctop, bcbot
    complex :: s(0:ny,nz) ! note zero-based arrays
    real :: c(0:ny)
    
    integer :: i,j,k,jp2,jm2,ib
    real    :: rj1,rj2,rj3
    
    do i = 2,nyhm
        j = 2*i - 2
        jp2 = j+2
        jm2 = j-2
        rj1 = 1.0/float(4*j*(j-1))
        rj2 = 1.0/float(2*(j*j-1))
        rj3 = 1.0/float(4*j*(j+1))
        do k = 1,nz
            wrk(i,k) = -c(jm2)*s(jm2,k)*rj1 + s(j,k)*rj2 - s(jp2,k)*rj3
        end do
    end do
    
    j = ny-2
    jm2 = j - 2
    rj1 = 1.0/float(4*j*(j-1))
    rj2 = 1.0/float(2*(j*j-1))
    do k = 1,nz
        wrk(nyh,k) = -s(jm2,k)*rj1 + s(j,k)*rj2
    end do
    
    j = ny
    jm2 = j - 2
    rj1 = 1.0/float(4*j*(j-1))
    do k = 1,nz
        wrk(nyhp,k) = -s(jm2,k)*rj1
    end do
    
    ! BCs
    if (ib .eq. 0) then
        do k = 1,nz
            wrk(1,k) = 0.5*(bctop(k) + bcbot(k))
        end do
    else
        do k = 1,nz
            wrk(1,k) = 0.5*(bctop(k) - bcbot(k))
        end do
    end if
    
    end subroutine evns
    
    !---------------------------------------------------------------------!
    
    subroutine odds(s,bctop,bcbot,ib,wrk)
    ! Sets up RHS for even tri-diagonal problem, given the results of sbr RHS
    use grid_size
    implicit none
    
    complex, dimension(nyp,nz) :: wrk
    complex, dimension(nz)     :: bctop, bcbot
    complex :: s(0:ny,nz) ! note zero-based arrays
    
    integer :: i,j,k,jp2,jm2,ib
    real    :: rj1,rj2,rj3
    
    do i = 2,nyhm-1
        j = 2*i - 1
        jp2 = j+2
        jm2 = j-2
        rj1 = 1.0/float(4*j*(j-1))
        rj2 = 1.0/float(2*(j*j-1))
        rj3 = 1.0/float(4*j*(j+1))
        do k = 1,nz
            wrk(i,k) = -s(jm2,k)*rj1 + s(j,k)*rj2 - s(jp2,k)*rj3
        end do
    end do
    
    i = nyhm
    j = ny-3
    jm2 = j - 2
    rj1 = 1.0/float(4*j*(j-1))
    rj2 = 1.0/float(2*(j*j-1))
    do k = 1,nz
        wrk(i,k) = -s(jm2,k)*rj1 + s(j,k)*rj2
    end do
    
    i = nyh
    j = ny - 1
    jm2 = j - 2
    rj1 = 1.0/float(4*j*(j-1))
    do k = 1,nz
        wrk(i,k) = -s(jm2,k)*rj1
    end do
    
    ! BCs
    if (ib .eq. 0) then
        do k = 1,nz
            wrk(1,k) = 0.5*(bctop(k) - bcbot(k))
        end do
    else
        do k = 1,nz
            wrk(1,k) = 0.5*(bctop(k) + bcbot(k))
        end do
    end if
    
    end subroutine odds
    
    !---------------------------------------------------------------------!
    
    subroutine eqsolve(a,t,wrk,nrank)
    ! Solves odd or even problem by Gaussian elimination
    use grid_size
    
    complex, dimension(nyp,nz)     :: wrk
    real,    dimension(nyhp,nz,3)  :: a
    real,    dimension(nyp)        :: t
    integer :: nrank
    
    real, dimension(nyp,nz) :: d
    real, dimension(nz)     :: f
    integer :: i,k,m,mm1,mm2
    
    m = nrank-1
    mm1 = m-1
    mm2 = m-2
    
    do k = 1,nz
        d(m,k) = a(m,k,2)
    end do
    
    do i = mm1,1,-1
        do k = 1,nz
            d(i,k) = a(i,k,2) - a(i,k,3)*a(i+1,k,1)/d(i+1,k)
            wrk(i+1,k) = wrk(i+1,k) - a(i,k,3)*wrk(i+2,k)/d(i+1,k)
        end do
    end do
    
    ! Eliminate the top row
    
    do k = 1,nz
        f(k) = t(m) - t(nrank)*a(m,k,1)/d(m,k)
        wrk(1,k) = wrk(1,k) - t(nrank)*wrk(nrank,k)/d(m,k)
    end do
    
    do i = mm1,1,-1
        do k = 1,nz
            wrk(1,k) = wrk(1,k) - f(k)*wrk(i+1,k)/d(i,k)
            f(k) = t(i) - f(k)*a(i,k,1)/d(i,k)
        end do
    end do
    
    ! Backward substitution: put solution in wrk variable
    do k = 1,nz
        wrk(1,k) = wrk(1,k)/f(k)
    end do
    
    do i = 2,nrank
        do k = 1,nz
            wrk(i,k) = (wrk(i,k) - a(i-1,k,1)*wrk(i-1,k))/d(i-1,k)
        end do
    end do
    
    end subroutine eqsolve
    
    subroutine phirhs(v,wrkc,wrk1,fn,fnm1,bcbot,bctop,var)
    ! If theta <= 0.99, works on both odd and even parts simultaneously
    
    !---------------------------------------------------------------------!
    !                                                                     !
    ! This subroutine evaluates the RHS of the Poisson equation for the   !
    ! particular phi.                                                     !
    !                                                                     !
    ! See p. 139 jfm v.177 (1987) KMM                                     !
    !                                                                     !
    ! RHS = -(1 - theta)/theta[phi"(n) - (kx^2 +kz^2)*phi(n)]             !
    !       - re*phi(n)/(dt*theta) - (3fn - fnm1)*re/(2*theta)            !
    !                                                                     !
    !---------------------------------------------------------------------!
    
    ! =================================================================== !
    !                           Declare Modules                           !
    ! =================================================================== !
    use grid_size
    use omp_lib
    use derivs
    !---------------------------------------------------------------------!
    
    ! =================================================================== !
    !                   Declare all variables explicitly                  !
    ! =================================================================== !
    implicit none
    
    ! Passed variables
    complex, dimension(nyp,nz,nxh) :: v,wrkc,wrk1,fn,fnm1
    complex, dimension(nz,nxh)     :: bcbot,bctop
    real :: var ! Allows for different variables to call this same subroutine
    
    ! Calculation variables
    integer :: i,j,k
    real    :: w2   
    
    ! Solver variables
    real :: gain,ugain,theta,alpha,beta,dyde
    real :: wavx(nxh),wavz(nz) 
    
    ! Flow setup variables
    real    :: re,Uinf,R_tau,dPdx
    
    ! Simulation variables
    integer :: it
    real    :: dt
    !---------------------------------------------------------------------!
    
    ! =================================================================== !
    !                            Common Blocks                            !
    ! =================================================================== !
    common/solver/     gain,ugain,theta,alpha,beta,dyde
    common/waves/      wavx,wavz
    common/flow/       re,Uinf,R_tau,dPdx
    common/itime/      it
    common/dtime/      dt
    !---------------------------------------------------------------------!
    
    ! =================================================================== !
    !                         Begin Calculations                          !
    ! =================================================================== !
   
    ! First reset BCs to 0 since they may have been altered in the
    ! pressure solution
    !$omp parallel do
    do k = 1,nxh
        do j = 1,nz
            bctop(j,k) = (0.0,0.0)
            bcbot(j,k) = (0.0,0.0)
        end do
    end do
    !$omp end parallel do 
    
    !$omp parallel do 
    do k = 1,nxh
        do j = 1,nz
            do i = 1,nyp
                wrkc(i,j,k) = -re*var/(dt*theta)*v(i,j,k) - re*var/(2.0*theta)*(3.0*fn(i,j,k) - fnm1(i,j,k))
            end do
        end do
    end do
    !$omp end parallel do 
    
    if (theta .gt. 0.99) then
        do k = 1,nxh
            do j = 1,nz
                do i = 1,nyp
                    v(i,j,k) = wrkc(i,j,k)
                end do
            end do
        end do
    end if
    
    !$omp parallel do 
    do k = 1,nxh
        do j = 1,nz
            w2 = wavz(j)**2 + wavx(k)**2
            do i = 1,nyp
                wrkc(i,j,k) = wrkc(i,j,k) + (1.0 - theta)*w2*v(i,j,k)/theta
            end do
        end do
    end do
    !$omp end parallel do 
    
    ! First derivative (stored in wrk1)
    call cderiv(v,wrk1)
    
    ! Second derivative (stored in v)
    call cderiv(wrk1,v)
    
    !$omp parallel do 
    do k = 1,nxh
        do j = 1,nz
            do i = 1,nyp
                v(i,j,k) = wrkc(i,j,k) - (1.0 - theta)*v(i,j,k)/theta
            end do
        end do
    end do
    !$omp end parallel do 
   
    end subroutine phirhs
    
    !---------------------------------------------------------------------!
    
    subroutine penta(u,x,g,dyde,ib,at,bt,gtop,ab,bb,gbot,a,wrkc)
    
    ! =================================================================== !
    !                           Declare Modules                           !
    ! =================================================================== !
    use grid_size
    !---------------------------------------------------------------------!
    
    ! =================================================================== !
    !                   Declare all variables explicitly                  !
    ! =================================================================== !
    implicit none
    
    ! Passed variables
    complex, dimension(nyp,nz) :: u,wrkc
    complex, dimension(nz)     :: gtop,gbot
    real,    dimension(nyhp,nz,3)  :: a
    integer :: ib
    real    :: x,g,dyde,at,bt,ab,bb
    
    ! Calculation variables
    real, dimension(nyp) :: t
    real, dimension(nz)  :: wn
    integer :: j,k,nrank,kzero,ip
    
    ! Solver variables
    real :: wavx(nxh),wavz(nz) 
    real :: c(nyp)
    
    ! Flow setup variables
    real    :: re,Uinf,R_tau,dPdx
    
    ! Simulation variables
    integer :: it
    real    :: dt
    !---------------------------------------------------------------------!
    
    ! =================================================================== !
    !                            Common Blocks                            !
    ! =================================================================== !
    common/waves/      wavx,wavz,c
    common/flow/       re,Uinf,R_tau,dPdx
    common/itime/      it
    common/dtime/      dt
    !---------------------------------------------------------------------!
    
    ! =================================================================== !
    !                         Begin Calculations                          !
    ! =================================================================== !
   
     
    call mixedbc(at,ab,bt,bb,dyde)
    
    do k = 1,nz
        wn(k) = wavz(k)**2 + x + g
    end do
    
    kzero = 0
    if (ib .ne. 0) then ! Strictly Neumann BCs
        do k = 1,nz
            if (wn(k) .le. 1.0e-20) then
                kzero = k
                wn(k) = 1.0
            end if
        end do
    end if
    
    ! Even problem
    ip = 1
    nrank = nyhp
    call evnc(a,wn,dyde)
    call pntevns(u,c,wrkc)
    call pntsolve(a,t,wrkc,gtop,gbot,ip,nrank)
    
    ! Odd problem
    ip = 2
    nrank = nyh
    call oddc(a,wn,dyde)
    call pntodds(u,wrkc)
    call pntsolve(a,t,wrkc,gtop,gbot,ip,nrank)
    
    ! Update u
    do k = 1,nz
        do j = 1,nyp
            u(j,k) = wrkc(j,k)
        end do
    end do
    if ((ib .ne. 0) .and. (kzero .ne. 0)) then
        do j = 1,nyp
            u(j,kzero) = 0.0
        end do
    end if
    
    end subroutine penta
    
    !---------------------------------------------------------------------!
    
    subroutine mixedbc(at,ab,bt,bb,dyde)
    
    ! Use module(s)
    use grid_size
    
    ! Declare variables
    implicit none
    
    real, dimension(nyp) :: tt,tb
    real    :: at,ab,bt,bb,dyde
    integer :: j,j1
    
    ! Common block(s)
    common/robin/ tt,tb
    
    ! Begin Calculations
    do j = 1,nyp
        j1 = j-1
        tt(j) = at + dyde*bt*float(j1*j1)
        tb(j) = ab + dyde*bb*float(j1*j1)
    end do
    
    do j = 2,nyp,2
        tb(j) = -tb(j)
    end do
    
    end subroutine mixedbc
    
    !---------------------------------------------------------------------!
    
    subroutine pntevns(s,c,wrk)
    ! Sets up RHS for even tri-diagonal problem given results of sbr RHS
    
    ! Use module(s)
    use grid_size
    
    ! Delcare variables
    implicit none
    
    complex, dimension(nyp,nz) :: wrk
    complex :: s(0:ny,nz) ! note zero-based arrays
    real :: c(0:ny)
    
    integer :: i,j,k,jp2,jm2
    real    :: rj1,rj2,rj3
    
    do i = 2,nyhm
        j = 2*i - 2
        jp2 = j+2
        jm2 = j-2
        rj1 = 1.0/float(4*j*(j-1))
        rj2 = 1.0/float(2*(j*j-1))
        rj3 = 1.0/float(4*j*(j+1))
        do k = 1,nz
            wrk(i,k) = -c(jm2)*s(jm2,k)*rj1 + s(j,k)*rj2 - s(jp2,k)*rj3
        end do
    end do
    
    j = ny-2
    jm2 = j - 2
    rj1 = 1.0/float(4*j*(j-1))
    rj2 = 1.0/float(2*(j*j-1))
    do k = 1,nz
        wrk(nyh,k) = -s(jm2,k)*rj1 + s(j,k)*rj2
    end do
    
    j = ny
    jm2 = j - 2
    rj1 = 1.0/float(4*j*(j-1))
    do k = 1,nz
        wrk(nyhp,k) = -s(jm2,k)*rj1
    end do
    
    end subroutine pntevns
    
    !---------------------------------------------------------------------!
    
    subroutine pntodds(s,wrk)
    ! Sets up RHS for even tri-diagonal problem, given the results of sbr RHS
    use grid_size
    implicit none
    
    complex, dimension(nyp,nz) :: wrk
    complex :: s(0:ny,nz) ! note zero-based arrays
    
    integer :: i,j,k,jp2,jm2
    real    :: rj1,rj2,rj3
    
    do i = 2,nyhm-1
        j = 2*i - 1
        jp2 = j+2
        jm2 = j-2
        rj1 = 1.0/float(4*j*(j-1))
        rj2 = 1.0/float(2*(j*j-1))
        rj3 = 1.0/float(4*j*(j+1))
        do k = 1,nz
            wrk(i,k) = -s(jm2,k)*rj1 + s(j,k)*rj2 - s(jp2,k)*rj3
        end do
    end do
    
    i = nyhm
    j = ny-3
    jm2 = j - 2
    rj1 = 1.0/float(4*j*(j-1))
    rj2 = 1.0/float(2*(j*j-1))
    do k = 1,nz
        wrk(i,k) = -s(jm2,k)*rj1 + s(j,k)*rj2
    end do
    
    i = nyh
    j = ny - 1
    jm2 = j - 2
    rj1 = 1.0/float(4*j*(j-1))
    do k = 1,nz
        wrk(i,k) = -s(jm2,k)*rj1
    end do
    
    end subroutine pntodds
    
    !---------------------------------------------------------------------!
    
    subroutine pntsolve(a,t,wrk,gtop,gbot,ip,nrank)
    
    ! Use module(s)
    use grid_size
    
    ! Declare variables
    implicit none
    
    complex, save, dimension(nyhp,nz,2) :: gs
    complex, save, dimension(2,nz)      :: s
    real,    save, dimension(nyhp,nz,2) :: cs,ds
    real,    save, dimension(2,nz,2)    :: f
    
    complex, dimension(nyp,nz) :: wrk
    complex, dimension(nz)     :: gtop,gbot
    real,    dimension(nyhp,nz,3)  :: a
    real,    dimension(nyp)        :: t
    
    integer :: i,k,ip,nrank,ipass,js,n,m,mm1,mm2,ii   
    
    ! Begin calculations
    
    ipass = 1
    if (ip .eq. 1) then
        do k = 1,nz
            s(1,k) = gtop(k)
            s(2,k) = gbot(k)
        end do
    end if
    
    js = 2
    n = nrank
    m = n - 1
    mm1 = m - 1
    mm2 = m - 2
    
    do i = 1,m
        do k = 1,nz
            gs(i,k,ip) = wrk(i,k)
            cs(i,k,ip) = a(i,k,1)
        end do
    end do
    
    do k = 1,nz
        ds(m,k,ip) = a(m,k,2)
        gs(n,k,ip) = wrk(n,k)
    end do
    
    do i = mm1,1,-1
        do k = 1,nz
            ds(i,k,ip) = a(i,k,2) - a(i,k,3)*cs(i+1,k,ip)/ds(i+1,k,ip)
            gs(i+1,k,ip) = gs(i+1,k,ip) - a(i,k,3)*gs(i+2,k,ip)/ds(i+1,k,ip)
        end do
    end do
    
    ! Eliminate 2nd row on ipass = 1 and 1st row on ipass = 2
    
2 continue
    
    call fillt(n,ip,ipass,t)
    
    do k = 1,nz
        f(js,k,ip) = t(m) - t(n)*cs(m,k,ip)/ds(m,k,ip)
        s(js,k) = s(js,k) - t(n)*gs(n,k,ip)/ds(m,k,ip)
    end do
    
    do i = mm1,1,-1
        do k = 1,nz
            s(js,k) = s(js,k) - f(js,k,ip)*gs(i+1,k,ip)/ds(i,k,ip)
            f(js,k,ip) = t(i) - f(js,k,ip)*cs(i,k,ip)/ds(i,k,ip)
        end do
    end do
    
    if (ipass .eq. 1) then
        ipass = ipass + 1
        js = js - 1
        goto 2
    end if
    
    ! Solve the two equations:
    !   f(1,k,1)*u0 + f(1,k,2)*u1 = s1
    !   f(2,k,1)*u0 + f(2,k,2)*u1 = s2
    ! Store solution in wrk
    
    if (ip .eq. 2 .and. ipass .eq. 2) then
        
        do k = 1,nz
            wrk(1,k) = (s(1,k)*f(2,k,2) - s(2,k)*f(1,k,2))/(f(1,k,1)*f(2,k,2) - f(2,k,1)*f(1,k,2))
            wrk(2,k) = (s(2,k)*f(1,k,1) - s(1,k)*f(2,k,1))/(f(1,k,1)*f(2,k,2) - f(2,k,1)*f(1,k,2))
        end do
    
        ! Even u
        n = n + 1
        do i = 2,n
            ii = 2*i - 1
            do k = 1,nz
                wrk(ii,k) = (gs(i,k,1) - cs(i-1,k,1)*wrk(ii-2,k))/ds(i-1,k,1)
            end do
        end do
    
        ! Odd u
        n = n - 1
        do i = 2,n
            ii = 2*i
            do k = 1,nz
                wrk(ii,k) = (gs(i,k,2) - cs(i-1,k,2)*wrk(ii-2,k))/ds(i-1,k,2)
            end do
        end do
    end if
    
    end subroutine pntsolve
    
    !---------------------------------------------------------------------!
    
    subroutine fillt(nrank,ip,ipass,t)
    
    ! Use module(s)
    use grid_size
    
    ! Declare variables
    implicit none
    
    real, dimension(nyp) :: tt,tb,t
    integer :: nrank,ip,ipass
    
    integer :: j,k
    
    ! Common block(s)
    common/robin/ tt,tb
    
    ! Begin Calculations
    
    if (ip .eq. 1 .and. ipass .eq. 1) then      ! Even, first pass
        do j = 1,nrank
            k = 2*j - 1
            t(j) = tb(k)
        end do
    else if (ip .eq. 1 .and. ipass .eq. 2) then ! Even, second pass
        do j = 1,nrank
            k = 2*j - 1
            t(j) = tt(k)
        end do
    else if (ip .eq. 2 .and. ipass .eq. 1) then  ! Odd, first pass
        do j = 1,nrank
            k = 2*j
            t(j) = tb(k)
        end do
    else                                         ! Odd, second pass
        do j = 1,nrank
            k = 2*j
            t(j) = tt(k)
        end do
    end if
    
    end subroutine fillt
    
    !---------------------------------------------------------------------!
    
    subroutine uzero(u0,h1,h1l,uflag)
    
    !---------------------------------------------------------------------!
    !                                                                     !
    !  This subroutine sets up to solve the x-momentum equation for       !
    !  kx = kz = 0, the zeroth Fourier mode of streamwise velocity.       !
    !                                                                     !
    !      d                  1   d  d(U0)                                !
    !      -- (U0) = H1 + F + --  -- ------                               !
    !      dt                 re  dy   dy                                 !
    !                                                                     !
    !  Discrete form:                                                     !
    !                                                                     !
    !      U0''(n+1) - (re/theta*dt) U0(n+1) = -(Re/theta*dt) U0(n)       ! 
    !                   - re*f/theta - (re/2*theta) (3h1(n) - h1(n-1))    !
    !                   - (1-theta) U0''(n)/theta                         !
    !                                                                     !
    !---------------------------------------------------------------------!
    
    ! =================================================================== !
    !                           Declare Modules                           !
    ! =================================================================== !
    use grid_size
    use derivs
    !---------------------------------------------------------------------!
    
    ! =================================================================== !
    !                   Declare all variables explicitly                  !
    ! =================================================================== !
    implicit none
    
    ! Passed variables
    real, dimension(nyp)       :: u0,h1,h1l
!    real, dimension(nyhp,nz,3) :: a
    integer :: uflag
    
    ! Calculation variables
    real, dimension(nyp) :: t,wrk1
    real    :: g
    integer :: i
    
    ! Solver variables
    real :: gain,ugain,theta,alpha,beta,dyde
    
    ! Flow setup variables
    real    :: re,Uinf,R_tau,dPdx
    
    ! Simulation variables
    integer :: it
    real    :: dt
    
    ! BC variables
    real :: atu,btu,gtu,abu,bbu,gbu
    real :: atw,btw,gtw,abw,bbw,gbw
    !---------------------------------------------------------------------!
    
    ! =================================================================== !
    !                            Common Blocks                            !
    ! =================================================================== !
    common/solver/ gain,ugain,theta,alpha,beta,dyde
    common/flow/   re,Uinf,R_tau,dPdx
    common/itime/  it
    common/dtime/  dt
    common/u0bcs/  atu,btu,gtu,abu,bbu,gbu
    common/w0bcs/  atw,btw,gtw,abw,bbw,gbw 
    !---------------------------------------------------------------------!
    
    ! =================================================================== !
    !                         Begin Calculations                          !
    ! =================================================================== !
 
    ! Note: f = 0 so we are neglecting that term in the discrete equation
    ! I don't think it's ever needed, but if it is some day, this is where
    ! you would put it ( just add "- f(i)*re/theta" to the end )
    do i = 1,nyp
        wrk1(i) = -re*u0(i)/(dt*theta) - re*(3.0*h1(i) - h1l(i))/(2.0*theta)
    end do
    
    if (theta .gt. 0.99) then
        do i = 1,nyp
            u0(i) = wrk1(i)
        end do
    else
    
        ! Calculate 2nd derivative of u0
        call cderiv1d(u0,t) ! First  derivative
        call cderiv1d(t,u0) ! Second derivative
        
        do i = 1,nyp
            u0(i) = wrk1(i) - (1.0 - theta)*u0(i)/theta
        end do
    end if
   
    ! Call the 1D penta-diagonal solver
    ! RHS on input is u0, solution is returned in u0
    ! Solver assumes mixed BCs are of the form:
    
    ! at*u0 + bt*du0/dy = gt at y = 1
    ! ab*u0 + bb*du0/dy = gb at y = -1
    
    ! BCs have been set in subroutine setbcs
    
    g = re/(theta*dt) ! time step term
    
    ! Switch for u0 and w0
    if (uflag .eq. 1) then
        call psolve1d(u0,g,dyde,atu,btu,gtu,abu,bbu,gbu,t)
    else 
        call psolve1d(u0,g,dyde,atw,btu,gtw,abw,bbw,gbw,t)
    end if
   
    end subroutine uzero
    
    !---------------------------------------------------------------------!
    
    subroutine psolve1d(u,g,dyde,at,bt,gtop,ab,bb,gbot,t)
    
    ! =================================================================== !
    !                           Declare Modules                           !
    ! =================================================================== !
    use grid_size
    !---------------------------------------------------------------------!
    
    ! =================================================================== !
    !                   Declare all variables explicitly                  !
    ! =================================================================== !
    implicit none
    
    ! Passed variables
    real, dimension(nyp)    :: u,t
    real, dimension(nyhp,3) :: a
    real :: g,dyde,gtop,gbot,at,bt,ab,bb
    
    ! Calculation variables
    real, dimension(nyp) :: wrkc
    real    :: x
    integer :: j,ip,nrank
    
    ! Solver variables
    real :: wavx(nxh),wavz(nz) 
    real :: c(nyp)
 
    ! Flow variables
    real :: re,Uinf,R_tau,dPdx
    
    !---------------------------------------------------------------------!
    
    ! =================================================================== !
    !                            Common Blocks                            !
    ! =================================================================== !
    common/flow/   re,Uinf,R_tau,dPdx
    common/waves/  wavx,wavz,c
    !---------------------------------------------------------------------!
    
    ! =================================================================== !
    !                         Begin Calculations                          !
    ! =================================================================== !
  
    x = 0.0
    call mixedbc1d(at,ab,bt,bb,dyde)
    
    ! Even problem
    ip = 1
    nrank = nyhp
    call p1devnc(a,x,g,wavz(1),dyde)
    call p1devns(u,c,wrkc)
    call p1dsolve(a,t,wrkc,gtop,gbot,ip,nrank)
    
    ! Odd problem
    ip = 2
    nrank = nyh
    call p1doddc(a,x,g,wavz(1),dyde)
    call p1dodds(u,wrkc)
    call p1dsolve(a,t,wrkc,gtop,gbot,ip,nrank)
    
    ! Update u
    do j = 1,nyp
        u(j) = wrkc(j)
    end do
    
    end subroutine psolve1d
    
    !---------------------------------------------------------------------!
    
    subroutine mixedbc1d(at,ab,bt,bb,dyde)
    
    use grid_size
    
    implicit none
    
    real, dimension(nyp) :: tt,tb
    integer :: j,j1
    real    :: at,ab,bt,bb,dyde
    
    common/robin/ tt,tb
    
    do j = 1,nyp
        j1 = j-1
        tt(j) = at + dyde*bt*float(j1*j1)
        tb(j) = ab - dyde*bb*float(j1*j1)
    end do
    
    do j = 2,nyp,2
        tb(j) = -tb(j)
    end do
    
    end subroutine mixedbc1d
    
    !---------------------------------------------------------------------!
    
    subroutine p1devnc(a,x,g,wavz,dyde)
    use grid_size
    implicit none
    
    real,    dimension(nyhp,3)  :: a
    real :: x,g,wavz,dyde,rj1,rj2,rj3,wn
    integer :: i,j
    
    wn = wavz**2 + x + g
    do i = 1,nyh
        j = 2*i
        rj1 = 1.0/float(4*j*(j-1))
        rj2 = 1.0/float(2*(j*j-1))
        rj3 = 1.0/float(4*j*(j+1))
    
        a(i,1) = wn*rj1
        a(i,2) = -1.0*dyde*dyde - wn*rj2
        a(i,3) = wn*rj3
    end do
    
    a(1,1) = 2.0*a(1,1)
    a(nyhm,3) = 0.0
    a(nyh,2) = -1.0*dyde*dyde
    a(nyh,3) = 0.0
    
    end subroutine p1devnc
    
    !---------------------------------------------------------------------!
    
    subroutine p1doddc(a,x,g,wavz,dyde)
    use grid_size
    implicit none
    
    real, dimension(nyhp,3)  :: a
    real :: x,g,wavz,dyde,rj1,rj2,rj3,wn
    integer :: i,j
    
    wn = wavz**2 + x + g
    
    do i = 1,nyhm
        j = 2*i + 1
        rj1 = 1.0/float(4*j*(j-1))  
        rj2 = 1.0/float(2*(j*j-1))  
        rj3 = 1.0/float(4*j*(j+1))  
    
        a(i,1) = wn*rj1
        a(i,2) = -1.0*dyde*dyde - wn*rj2
        a(i,3) = wn*rj3
    end do
    
    a(nyhm-1,3) = 0.0
    a(nyhm,2) = -1.0*dyde*dyde
    a(nyhm,3) = 0.0
    
    end subroutine p1doddc
    
    !---------------------------------------------------------------------!
    
    subroutine p1devns(s,c,wrk)
    ! Sets up RHS for even tri-diagonal problem, given the results of sbr RHS
    use grid_size
    implicit none
    
    real :: s(0:ny), c(0:ny), wrk(nyp) ! note zero-based arrays
    
    integer :: i,j,jp2,jm2
    real    :: rj1,rj2,rj3
    
    do i = 2,nyhm
        j = 2*i - 2
        jp2 = j+2
        jm2 = j-2
        rj1 = 1.0/float(4*j*(j-1))
        rj2 = 1.0/float(2*(j*j-1))
        rj3 = 1.0/float(4*j*(j+1))
    
        wrk(i) = -c(jm2)*s(jm2)*rj1 + s(j)*rj2 - s(jp2)*rj3
    end do
    
    j = ny-2
    jm2 = j - 2
    rj1 = 1.0/float(4*j*(j-1))
    rj2 = 1.0/float(2*(j*j-1))
    wrk(nyh) = -s(jm2)*rj1 + s(j)*rj2
    
    j = ny
    jm2 = j - 2
    rj1 = 1.0/float(4*j*(j-1))
    wrk(nyhp) = -s(jm2)*rj1
    
    end subroutine p1devns
    
    !---------------------------------------------------------------------!
    
    subroutine p1dodds(s,wrk)
    ! Sets up RHS for even tri-diagonal problem, given the results of sbr RHS
    use grid_size
    implicit none
    
    real    :: s(0:ny),wrk(nyp) ! Note zero-based array
    
    integer :: i,j,jp2,jm2
    real    :: rj1,rj2,rj3
    
    do i = 2,nyhm-1
        j = 2*i - 1
        jp2 = j+2
        jm2 = j-2
        rj1 = 1.0/float(4*j*(j-1))
        rj2 = 1.0/float(2*(j*j-1))
        rj3 = 1.0/float(4*j*(j+1))
        wrk(i) = -s(jm2)*rj1 + s(j)*rj2 - s(jp2)*rj3
    end do
    
    i = nyhm
    j = ny-3
    jm2 = j - 2
    rj1 = 1.0/float(4*j*(j-1))
    rj2 = 1.0/float(2*(j*j-1))
    wrk(i) = -s(jm2)*rj1 + s(j)*rj2
    
    i = nyh
    j = ny - 1
    jm2 = j - 2
    rj1 = 1.0/float(4*j*(j-1))
    wrk(i) = -s(jm2)*rj1
    
    end subroutine p1dodds
    
    !---------------------------------------------------------------------!
    
    subroutine p1dsolve(a,t,wrk,gamt,gamb,ip,nrank)
    
    use grid_size
    
    implicit none
    
    real, save, dimension(nyhp,2) :: cs,ds,gs
    real, save, dimension(2,2)    :: f
    real, save, dimension(2)      :: s
    real,       dimension(nyp)    :: wrk
    real,    dimension(nyhp,3)    :: a
    real,    dimension(nyp)       :: t
    
    real    :: gamb,gamt
    integer :: ip,ipass,nrank,i,js,n,m,mm1,mm2,ii
    
    ! Begin calculations
    ipass = 1
    if (ip .eq. 1) then
        s(1) = gamt
        s(2) = gamb
    end if
    
    js = 2
    n = nrank
    m = n - 1
    mm1 = m - 1
    mm2 = m - 2
    
    do i = 1,m
        gs(i,ip) = wrk(i)
        cs(i,ip) = a(i,1)
    end do
    
    ds(m,ip) = a(m,2)
    gs(n,ip) = wrk(n)
    
    do i = mm1,1,-1
        ds(i,ip) = a(i,2) - a(i,3)*cs(i+1,ip)/ds(i+1,ip)
        gs(i+1,ip) = gs(i+1,ip) - a(i,3)*gs(i+2,ip)/ds(i+1,ip)
    end do
    
    ! Eliminate 2nd row on ipass = 1 and first row on ipass = 2
    
2 continue
    
    call fillt(n,ip,ipass,t)
    
    f(js,ip) = t(m)  - t(n)*cs(m,ip)/ds(m,ip)
    s(js)    = s(js) - t(n)*gs(n,ip)/ds(m,ip)
    
    do i = mm1,1,-1
        s(js) = s(js) - f(js,ip)*gs(i+1,ip)/ds(i,ip)
        f(js,ip) = t(i) - f(js,ip)*cs(i,ip)/ds(i,ip)
    end do
    
    if (ipass .eq. 1) then
        ipass = ipass + 1
        js = js - 1
        goto 2
    end if
    
    ! Solve the two equations:
    !   f(1,1)*u0 + f(1,2)*u1 = s1
    !   f(2,1)*u0 + f(2,2)*u1 = s2
    ! Store solution in wrk
    
    if (ip .eq. 2 .and. ipass .eq. 2) then
        
        wrk(1) = (s(1)*f(2,2) - s(2)*f(1,2))/(f(1,1)*f(2,2) - f(2,1)*f(1,2))
        wrk(2) = (s(2)*f(1,1) - s(1)*f(2,1))/(f(1,1)*f(2,2) - f(2,1)*f(1,2))
    
        ! Even u
        n = n + 1
        do i = 2,n
            ii = 2*i - 1
            wrk(ii) = (gs(i,1) - cs(i-1,1)*wrk(ii-2))/ds(i-1,1)
        end do
    
        ! Odd u
        n = n - 1
        do i = 2,n
            ii = 2*i
            wrk(ii) = (gs(i,2) - cs(i-1,2)*wrk(ii-2))/ds(i-1,2)
        end do
    end if
    
    end subroutine p1dsolve
    
    !---------------------------------------------------------------------!
    
    subroutine nonlin(fn,omz,gn,wrkc,wrk1)
    !---------------------------------------------------------------------!
    !  This subroutine calculates the nonlinear term in the 4th order     !
    !  system:                                                            !
    !                                                                     !
    !     fn = -d/dy(d/dx(H1) + d/dz(H3)) + d/dz(d/dz(H2))                !
    !                                                                     !
    !     gn = d/dz(H1) - d/dx(H3)                                        !
    !---------------------------------------------------------------------!
    ! =================================================================== !
    !                           Declare Modules                           !
    ! =================================================================== !
    use grid_size
    use derivs
    !---------------------------------------------------------------------!
    
    ! =================================================================== !
    !                   Declare all variables explicitly                  !
    ! =================================================================== !
    implicit none
    
    ! Passed variables
    complex, dimension(nyp,nz,nxh) :: fn,gn,omz,wrkc,wrk1
    
    ! Calculation variables
    complex :: im
    real    :: w2
    integer :: i,j,k
    
    ! Solver variables
    real :: wavx(nxh),wavz(nz) 
    
    !---------------------------------------------------------------------!
    
    ! =================================================================== !
    !                            Common Blocks                            !
    ! =================================================================== !
    common/waves/      wavx,wavz
    !---------------------------------------------------------------------!
    
    ! =================================================================== !
    !                         Begin Calculations                          !
    ! =================================================================== !
    
    im = (0.0,1.0)
    
    do k = 1,nxh
        do j = 1,nz
            do i = 1,nyp
                wrkc(i,j,k) = im*(wavz(j)*omz(i,j,k) + wavx(k)*gn(i,j,k))
            end do
        end do
    end do
    
    call cderiv(wrkc,wrk1)
    
    ! Calculate fn
    do k = 1,nxh
        do j = 1,nz
            w2 = -wavx(k)**2 - wavz(j)**2
            do i = 1,nyp
                fn(i,j,k) = w2*fn(i,j,k) - wrk1(i,j,k)
            end do
        end do
    end do
    
    ! Calculate gn
    do k = 1,nxh
        do j = 1,nz
            do i = 1,nyp
                gn(i,j,k) = im*(wavz(j)*gn(i,j,k) - wavx(k)*omz(i,j,k))
            end do
        end do
    end do
    
    end subroutine nonlin
    
    !---------------------------------------------------------------------!
    
    subroutine veloc(u,v,w,u0,w0,omy,wrkc)
    !---------------------------------------------------------------------!
    !  This subroutine calculates the remaining components of velocity    !
    !  from continuity and normal vorticity in spectral space.            !
    !                                                                     !
    !    U(n,kz,kx) = -((i*kx)(-dV/dy) + (i*kz)(omy))/(kx^2 + kz^2)       !
    !  for kx^2 + kz^2 > 0. Otherwise, make u0 the real part of U(i,kx=0) !
    !                                                                     !
    !    W(n,kz,kx) = -((i*kz)(-dV/dy) - (i*kx)(omy))/(kx^2 + kz^2)       !
    !  for kx^2 + kz^2 > 0. Otherwise, W(n,kz,kx) = 0.0                   !
    !                                                                     !
    !---------------------------------------------------------------------!
    ! =================================================================== !
    !                           Declare Modules                           !
    ! =================================================================== !
    use grid_size
    use derivs
    !---------------------------------------------------------------------!
    
    ! =================================================================== !
    !                   Declare all variables explicitly                  !
    ! =================================================================== !
    implicit none
    
    ! Passed variables
    complex, dimension(nyp,nz,nxh) :: u,v,w,omy,wrkc
    real,    dimension(nyp)        :: u0,w0
    
    ! Calculation variables
    complex :: im
    real    :: w2
    integer :: i,j,k
    
    ! Solver variables
    real :: wavx(nxh),wavz(nz) 
    
    !---------------------------------------------------------------------!
    
    ! =================================================================== !
    !                            Common Blocks                            !
    ! =================================================================== !
    common/waves/      wavx,wavz
    !---------------------------------------------------------------------!
    
    ! =================================================================== !
    !                         Begin Calculations                          !
    ! =================================================================== !
    
    im = (0.0,1.0)
    
    call cderiv(v,wrkc)
    
    ! Calculate for all modes except the first (k=j=1)
    do k = 2,nxh
        do j = 1,nz
            w2 = wavx(k)**2 + wavz(j)**2
            do i = 1,nyp
                u(i,j,k) = -(-im*wavx(k)*wrkc(i,j,k) + im*wavz(j)*omy(i,j,k))/w2
                w(i,j,k) = -(-im*wavz(j)*wrkc(i,j,k) - im*wavx(k)*omy(i,j,k))/w2
            end do
        end do
    end do
    
    k = 1 ! Skip only j = k = 1 mode
    do j = 2,nz
        w2 = wavz(j)**2
        do i = 1,nyp
            u(i,j,k) = -(im*wavz(j)*omy(i,j,k))/w2
            w(i,j,k) = -(-im*wavz(j)*wrkc(i,j,k))/w2
        end do
    end do
    
    ! Calculate first mode:
    ! Update the first (kz = 0) mode of streamwise velocity. The imaginary part
    ! of the first mode has been used to store the real part of the last Fourier
    ! mode (nz/2 + 1)
    
    do i = 1,nyp
        u(i,1,1) = cmplx(u0(i),0.0)
        v(i,1,1) = 0.0
        w(i,1,1) = cmplx(w0(i),0.0)
    end do
    
    end subroutine veloc
    
    !---------------------------------------------------------------------!
    
    subroutine vort(u,v,w,omx,omz)
    !---------------------------------------------------------------------!
    !  This subroutine calculates the x and z vorticity componenents:     !
    !                                                                     !
    !    omx = dw/dy - dv/dz                                              !
    !    omy = du/dz - dw/dx (calculated prior)                           !
    !    omz = dv/dx - du/dy                                              !
    !                                                                     !
    !    All quantities are spectral                                      !
    !                                                                     !
    !---------------------------------------------------------------------!
    ! =================================================================== !
    !                           Declare Modules                           !
    ! =================================================================== !
    use grid_size
    use derivs
    !---------------------------------------------------------------------!
    
    ! =================================================================== !
    !                   Declare all variables explicitly                  !
    ! =================================================================== !
    implicit none
    
    ! Passed variables
    complex, dimension(nyp,nz,nxh) :: u,v,w,omx,omz
    
    ! Calculation variables
    complex :: im
    integer :: i,j,k
    
    ! Solver variables
    real :: wavx(nxh),wavz(nz) 
    
    !---------------------------------------------------------------------!
    
    ! =================================================================== !
    !                            Common Blocks                            !
    ! =================================================================== !
    common/waves/      wavx,wavz
    !---------------------------------------------------------------------!
    
    ! =================================================================== !
    !                         Begin Calculations                          !
    ! =================================================================== !
    
    im = (0.0,1.0)
    
    call cderiv(u,omz) ! store du/dy in omz
    call cderiv(w,omx) ! store dw/dy in omx
   
    do k = 1,nxh
    do j = 1,nz
        do i = 1,nyp
            omz(i,j,k) = im*wavx(k)*v(i,j,k) - omz(i,j,k) 
            omx(i,j,k) = omx(i,j,k) - im*wavz(j)*v(i,j,k)
        end do
    end do
    end do
    
    end subroutine vort
    
    !---------------------------------------------------------------------!
    
    subroutine vcw3dp(u,v,w,omx,omy,omz,fn,gn,u11,u12,u13,u21,u22,u23,       &
                    u31,u32,u33,Lu,Lv,Lw, &
#IFDEF SCALAR
                    scalar,sclx,scly,sclz,scn,        &
#ENDIF
#IFDEF POLYMER
                    c11,c12,c13,c21,c22,c23,c31,c32,c33,                   &
                    dc111,dc112,dc113,dc211,dc212,dc213,dc311,dc312,dc313, &
                    dc121,dc122,dc123,dc221,dc222,dc223,dc321,dc322,dc323, &
                    dc131,dc132,dc133,dc231,dc232,dc233,dc331,dc332,dc333, &
                    c11n,c12n,c13n,c22n,c23n,c33n,                         &
                    str11n,str12n,str13n,str22n,str23n,str33n,             &
                    qp11,qp12,qp13,qp22,qp23,qp33,  &
#ENDIF
                    Lu_old,Lv_old,Lw_old, &
                    planZb,planXb,planY,planXf,planZf) 

    !---------------------------------------------------------------------!
    !  This subroutine calculates the nonlinear term in the rotational    !
    !  form of the N-S equations (3D calculation)                         !
    !                                                                     !
    !    h1 = (u x omega)_x                                               !
    !    h2 = (u x omega)_y                                               !
    !    h3 = (u x omega)_z                                               !
    !                                                                     !
    !---------------------------------------------------------------------!
    !  The velocity and vorticity fields coming in are assumed to be in   !
    !  spectral space in all three dimensions. We also assume that the    !
    !  appropriate modes in all fields (i.e., the first and last modes)   !
    !  are set to zero prior to this point.                               !
    !---------------------------------------------------------------------!
    ! =================================================================== !
    !                           Declare Modules                           !
    ! =================================================================== !
    use,intrinsic :: iso_c_binding
    use omp_lib
    use grid_size
    use derivs
    use helpers
    !---------------------------------------------------------------------!
    
    ! =================================================================== !
    !                   Declare all variables explicitly                  !
    ! =================================================================== !
    implicit none
  
    include 'fftw3.f03'
 
    integer, parameter :: qn = 10000
 
    ! Passed variables
    type(C_PTR) :: planZb,planXb,planY,planXf,planZf

    ! Input/Output variables
    complex, dimension(nyp,nz,nxh) :: u,v,w,omx,omy,omz,fn,gn
    complex, dimension(nyp,nz,nxh) :: u11,u12,u13,u21,u22,u23,u31,u32,u33
    complex, dimension(nyp,nz,nxh) :: Lu,Lv,Lw
#IFDEF SCALAR
    complex, dimension(nyp,nz,nxh) :: scalar,sclx,scly,sclz,scn
#ENDIF
#IFDEF POLYMER
    complex, dimension(nyp,nz,nxh) :: c11,c12,c13,c21,c22,c23,c31,c32,c33
    complex, dimension(nyp,nz,nxh) :: c11n,c12n,c13n,c22n,c23n,c33n

    complex, dimension(nyp,nz,nxh) :: dc111,dc112,dc113,dc121,dc122,dc123,dc131,dc132,dc133
    complex, dimension(nyp,nz,nxh) :: dc211,dc212,dc213,dc221,dc222,dc223,dc231,dc232,dc233
    complex, dimension(nyp,nz,nxh) :: dc311,dc312,dc313,dc321,dc322,dc323,dc331,dc332,dc333

    complex, dimension(nyp,nz,nxh) :: str11n,str12n,str13n,str22n,str23n,str33n

    complex, dimension(nyp,nz,nxh) :: qp11,qp12,qp13,qp22,qp23,qp33
    real,    dimension(nyp,mz,mx)  :: beta_poly
#ENDIF
   
    ! FFT complex arrays 
    complex(C_DOUBLE_COMPLEX), dimension(nyp,mz,mx) :: us,vs,ws,wxs,wys,wzs
    complex(C_DOUBLE_COMPLEX), dimension(nyp,mz,mx) :: u11s,u12s,u13s,u21s,u22s,u23s,u31s,u32s,u33s
    complex(C_DOUBLE_COMPLEX), dimension(nyp,mz,mx) :: Lus,Lvs,Lws

    complex(C_DOUBLE_COMPLEX), dimension(nyp,mz,mx) :: vwxs,vwys,vwzs
#IFDEF SCALAR
    complex(C_DOUBLE_COMPLEX), dimension(nyp,mz,mx) :: scs,cxs,cys,czs
    complex(C_DOUBLE_COMPLEX), dimension(nyp,mz,mx) :: vcs
#ENDIF
#IFDEF POLYMER
    complex(C_DOUBLE_COMPLEX), dimension(nyp,mz,mx) :: c11s,c12s,c13s,c21s,c22s,c23s,c31s,c32s,c33s
    complex(C_DOUBLE_COMPLEX), dimension(nyp,mz,mx) :: c11ns,c12ns,c13ns,c22ns,c23ns,c33ns

    complex(C_DOUBLE_COMPLEX), dimension(nyp,mz,mx) :: dc111s,dc112s,dc113s,dc121s,dc122s,dc123s,dc131s,dc132s,dc133s
    complex(C_DOUBLE_COMPLEX), dimension(nyp,mz,mx) :: dc211s,dc212s,dc213s,dc221s,dc222s,dc223s,dc231s,dc232s,dc233s
    complex(C_DOUBLE_COMPLEX), dimension(nyp,mz,mx) :: dc311s,dc312s,dc313s,dc321s,dc322s,dc323s,dc331s,dc332s,dc333s

    complex(C_DOUBLE_COMPLEX), dimension(nyp,mz,mx) :: str11ns,str12ns,str13ns,str22ns,str23ns,str33ns

    complex(C_DOUBLE_COMPLEX), dimension(nyp,mz,mx) :: qp11s,qp12s,qp13s,qp22s,qp23s,qp33s
#ENDIF
    
    ! Calculation variables
    ! FFTW variables - must be C-type arrays
    ! Input variables
    real, save, dimension(nyp,mz,mx) :: u_old,v_old,w_old
    real, dimension(nyp,mz,mx) :: Lu_old,Lv_old,Lw_old

    real(C_DOUBLE), dimension(nyp,mz,mx) :: up,vp,wp,wxp,wyp,wzp
    real(C_DOUBLE), dimension(nyp,mz,mx) :: u11p,u12p,u13p
    real(C_DOUBLE), dimension(nyp,mz,mx) :: u21p,u22p,u23p
    real(C_DOUBLE), dimension(nyp,mz,mx) :: u31p,u32p,u33p
    real(C_DOUBLE), dimension(nyp,mz,mx) :: Lup,Lvp,Lwp,swirl_3d
    
    real(C_DOUBLE), dimension(nyp,mz,mx) :: vwx,vwy,vwz

#IFDEF POLYMER    
    real(C_DOUBLE), dimension(nyp,mz,mx)  :: c11p,c12p,c13p
    real(C_DOUBLE), dimension(nyp,mz,mx)  :: c21p,c22p,c23p
    real(C_DOUBLE), dimension(nyp,mz,mx)  :: c31p,c32p,c33p

    real(C_DOUBLE), dimension(nyp,mz,mx)  :: trp

    real(C_DOUBLE), dimension(nyp,mz,mx)  :: dc111p, dc112p, dc113p
    real(C_DOUBLE), dimension(nyp,mz,mx)  :: dc211p, dc212p, dc213p
    real(C_DOUBLE), dimension(nyp,mz,mx)  :: dc311p, dc312p, dc313p
    
    real(C_DOUBLE), dimension(nyp,mz,mx)  :: dc121p, dc122p, dc123p
    real(C_DOUBLE), dimension(nyp,mz,mx)  :: dc221p, dc222p, dc223p
    real(C_DOUBLE), dimension(nyp,mz,mx)  :: dc321p, dc322p, dc323p
    
    real(C_DOUBLE), dimension(nyp,mz,mx)  :: dc131p, dc132p, dc133p
    real(C_DOUBLE), dimension(nyp,mz,mx)  :: dc231p, dc232p, dc233p
    real(C_DOUBLE), dimension(nyp,mz,mx)  :: dc331p, dc332p, dc333p

    real(C_DOUBLE), dimension(nyp,mz,mx)  :: c11np,c12np,c13np
    real(C_DOUBLE), dimension(nyp,mz,mx)  :: c22np,c23np,c33np
    
    real(C_DOUBLE), dimension(nyp,mz,mx)  :: str11np,str12np,str13np
    real(C_DOUBLE), dimension(nyp,mz,mx)  :: str22np,str23np,str33np
    
    real(C_DOUBLE), dimension(nyp,mz,mx)  :: qp11np,qp12np,qp13np
    real(C_DOUBLE), dimension(nyp,mz,mx)  :: qp22np,qp23np,qp33np
#ENDIF
#IFDEF SCALAR    
    real(C_DOUBLE), dimension(nyp,mz,mx) :: scp,cxp,cyp,czp,vc

    real, dimension(nyp,mz,mx) :: scsource
    real, dimension(nyp,nz,nx) :: Qcrit
    real, dimension(qn)     :: Qmin = 0.0
    integer, dimension(qn)  :: Qx,Qy,Qz 
#ENDIF
    
    integer :: i,j,k,n
    integer :: ii,jj,ipii,jpjj,idif,jdif
    real    :: delxm,delzm,pi!,vArea
    real    :: segdrag,xsegdrag,ysegdrag,zsegdrag,wdes
    real    :: cflcheck,cflmax,swirl
    
    real,dimension(nyp)   :: cfl
    real,dimension(nyp,mz,mx) :: dragx,dragy,dragz
    
    real :: xi,zj,argx,argrad,fx,fr,rsq
    real :: uxmean(mz),uzmean(nyp)
    real :: massFlux
    
    ! FFT variables
    real :: trigx(2*nx),trigz(2*nz),trigy(2*ny),sine(ny),cosine(ny)
    integer, dimension(19) :: ixfax,iyfax,izfax,ixfax32,izfax32
    real, dimension(16000) :: trigx32
    real, dimension(4000)  :: trigz32 
    
    ! Simulation control variables
    integer :: irstrt,nsteps,iprnfrq,print3d
    integer :: it
    real    :: dt
    
    ! Solver variables
    real :: gain,ugain,theta,alpha,beta,dyde
    real :: wavx(nxh),wavz(nz) 
    real :: c(nyp)
    real :: fac
    
    ! Buffer variables
    real, dimension(1200) :: bfgain,bfugain
    real, dimension(mx)   :: vdes
    integer :: bfhead,bftail,bfwidth
    real    :: slopelength
    
    ! Geometry variables
    real    :: xl,yl,zl
    integer :: kwall,kmaxsurf
    integer, dimension(nyp,mz,mx) :: imatrix
    
    ! Flow setup variables
    real, dimension(nyp,nz,nx) :: initu,initv,initw
    real    :: re,Uinf,R_tau,dPdx
    integer :: geomtype,flow_select,perturbtime
    
    ! Vortex variables
    real    :: forbeta,xcenter,ycenter,zcenter,L,rad,bdyfx
    
    ! Immersed boundary force variables
    real, dimension(nyp,mz,mx) :: fxintg,fyintg,fzintg
    real, dimension(mz2,mx2)   :: fspread

#IFDEF SCALAR
    ! Scalar variables
    real    :: sigmax,sigmay,sigmaz
    real    :: betax,betay,betaz
    real    :: xsq,ysq,zsq
    real    :: xc1,yc1,zc1
    real    :: deltaT,diff
    integer :: scl_flag,scltarg
    integer :: src_start,src_stop
    logical :: condition
#ENDIF

    ! Particle variables
    real,dimension(npart) :: xpart,ypart,zpart,swirl_part

#IFDEF POLYMER
    ! Polymer variables
    integer :: ipolyflag,itarget,ipeter
    real    :: alpha_poly,tpoly,zlmax,diffpoly,qbeta
    real    :: zbeta1
#ENDIF

    !---------------------------------------------------------------------!
    
    ! =================================================================== !
    !                            Common Blocks                            !
    ! =================================================================== !
    common/iocontrl/   irstrt,nsteps,iprnfrq,print3d
    common/solver/     gain,ugain,theta,alpha,beta,dyde
    common/ibforce/    fxintg,fyintg,fzintg,fspread
    common/waves/      wavx,wavz,c
    common/init/       initu,initv,initw
    common/buffer/     bfgain,bfugain,vdes,bfhead,bftail,bfwidth
    common/slope/      slopelength
    common/flow/       re,Uinf,R_tau,dPdx
    common/domain/     xl,yl,zl
    common/trig/       trigx,trigy,trigz,sine,cosine,trigx32,trigz32
    common/ifax/       ixfax,iyfax,izfax,ixfax32,izfax32
    common/itime/      it
    common/dtime/      dt
    common/imat/       imatrix,kwall,kmaxsurf
    common/vortexring/ forbeta,xcenter,ycenter,zcenter,L,rad,bdyfx
    common/setup/      geomtype,flow_select,perturbtime
    common/part_traj/  xpart,ypart,zpart,swirl_part
#IFDEF SCALAR
    common/scl_stuff/  sigmax,sigmay,sigmaz,deltaT,diff,scl_flag,scltarg
    common/src_time/   src_start,src_stop
#ENDIF
#IFDEF POLYMER
    common/poly_flgs/  ipolyflag,itarget,ipeter
    common/poly_var/   alpha_poly,tpoly,zlmax,diffpoly,qbeta
#ENDIF
    !---------------------------------------------------------------------!
    
    ! =================================================================== !
    !                         Begin Calculations                          !
    ! =================================================================== !
   
!    write(100+it,*) u 
!    write(200+it,*) v 
!    write(300+it,*) w 
    !---------------------------------------------------------------------!
    !                    Initialize some variables                        !
    !---------------------------------------------------------------------!
 
    pi = 2.0*acos(0.0)
    delxm = xl/float(mx-1)
    delzm = zl/float(mz-1)
   
#IFDEF SCALAR 
    !---------------------------------------------------------------------!
    !                    Calculate scalar source term                     !
    !---------------------------------------------------------------------!
    if (scl_flag .eq. 2 .and. it .ge. src_start .and. it .le. src_stop) then
        scsource = 0.0
        do n = 1,npart
            if (it .eq. irstrt) then
                open(95,file='setup/particles/particles.dat',status='old',action='read')
                do j = 1,n
                    read(95,*)
                end do
                read(95,*) xc1,yc1,zc1
                close(95)
            else if (it .gt. irstrt) then
                xc1 = xpart(n)
                yc1 = ypart(n)
                zc1 = zpart(n)
                swirl = swirl_part(n)
            end if
    
            ! Choose targeting case
            if (scltarg .eq. 1) then
                condition = swirl .gt. 2000.0 ! Target high Q
            else if (scltarg .eq. 2) then
                condition = swirl .lt. -2000.0 ! Target low Q
            else
                condition = .true. ! No targeting
            end if
    
            if (condition) then 
            do k = 1,mx
                xsq = (float(k-1)*delxm - xc1)**2
                betax = xsq/(2.0*sigmax**2)
                do j = 1,mz
                    zsq = (float(j-1)*delzm - zc1)**2
                    betaz = zsq/(2.0*sigmaz**2)
                    do i = 1,nyp
                        ysq = (ycoord(i) - yc1)**2
                        betay = ysq/(2.0*sigmay**2)
    
                        scsource(i,j,k) = scsource(i,j,k) + deltaT*exp(-(betax + betay + betaz))
                    end do
                end do
            end do
            end if
        end do
#IFDEF POLYMER
    else if (scl_flag .eq. 4 .and. it .ge. src_start .and. it .le. src_stop) then
        Qx = 0
        Qy = 0
        Qz = 0
        call findmaxQ(u11,u12,u13,u21,u22,u23,u31,u32,u33,scalar,Qmin,Qx,Qy,Qz,beta3d)
        scsource = 0.0
        n = 1
        cnt = 1
        do while (cnt .lt. 101)
            kk = Qx(n)
            ii = Qy(n)
            jj = Qz(n)
            if (beta3d(ii,jj,kk) .gt. 0.9) then 
                xc1 = delxm*(kk-1)
                yc1 = ycoord(ii)
                zc1 = delzm*(jj-1)
                do k = 1,mx
                    xsq = (float(k-1)*delxm - xc1)**2
                    betax = xsq/(2.0*sigmax**2)
                    do j = 1,mz
                        zsq = (float(j-1)*delzm - zc1)**2
                        betaz = zsq/(2.0*sigmaz**2)
                        do i = 1,nyp
                            ysq = (ycoord(i) - yc1)**2
                            betay = ysq/(2.0*sigmay**2)

                            scsource(i,j,k) = scsource(i,j,k) + deltaT*exp(-(betax + betay + betaz))
                        end do
                    end do
                end do
                cnt = cnt + 1
            end if
            n = n+1
        end do
#ENDIF
    else if (scl_flag .eq. 3 .and. it .eq. 1 .and. irstrt .eq. 1) then
        scsource = 1.0
    else
        scsource = 0.0
    end if
#ENDIF

    ! Set cfl variables to 0
    cfl = 0.0
    cflcheck = 0.0
    cflmax = 0.0
   

    
    !----------------------------------------------------------------------!
    !  Compute FFTs on all spectral data to transform into real 3/2 space  ! 
    !----------------------------------------------------------------------!
    ! We take the normal coefficients (spectral field data) and we         !
    ! interpolate them onto a 3/2 grid so we can "de-alias" the last 1/3   !
    ! spectral modes (filter them out) when we transform them back.        !
    ! We do this during the nonlinear term calculations because this is    !
    ! most likely place for a numerical instability to manifest.           !
    !                                                                      !
    ! The FFTW routines do this interpolation automatically when we        !
    ! transform using different sizes. However, it requires the smaller    !
    ! array to be "padded" with zeros in order for the algorithm to work.  !
    ! So we must copy the (nyp)x(nz)x(nxh) variables over to a larger      !
    ! array of size (nyp)x(mz)x(mx) before the FFTs can be performed.      !
    !                                                                      !
    ! It is assumed that the FFTW plans are created before the first call  !
    ! to this subroutine, so all we have to do is execute them with the    !
    ! appropriate arguments.                                               !
    !----------------------------------------------------------------------!


    ! Zero out arrays
    !$omp parallel do default(shared) private(i,j,k)
    do k = 1,mx
        do j = 1,mz
            do i = 1,nyp
                us(i,j,k) = 0.0
                vs(i,j,k) = 0.0
                ws(i,j,k) = 0.0

                wxs(i,j,k) = 0.0
                wys(i,j,k) = 0.0
                wzs(i,j,k) = 0.0

                u11s(i,j,k) = 0.0
                u12s(i,j,k) = 0.0
                u13s(i,j,k) = 0.0
                u21s(i,j,k) = 0.0
                u22s(i,j,k) = 0.0
                u23s(i,j,k) = 0.0
                u31s(i,j,k) = 0.0
                u32s(i,j,k) = 0.0
                u33s(i,j,k) = 0.0

                Lus(i,j,k) = 0.0
                Lvs(i,j,k) = 0.0
                Lws(i,j,k) = 0.0
            end do
        end do
    end do
    !$omp end parallel do

    ! Copy spectral variables into larger arrays for transforms 
    ! Also convert Chebyshev modes to cosine modes
    !$omp parallel do default(shared) private(i,j,k,jj)
    do k = 1,nxh
        do j = 1,nz 
            if (j .le. nzh) jj = j
            if (j .gt. nzh) jj = (mz-nz) + j
            do i = 1,nyp

                fac = c(i)/2.0

                ! Velocity Field
                us(i,jj,k) = u(i,j,k)*fac
                vs(i,jj,k) = v(i,j,k)*fac
                ws(i,jj,k) = w(i,j,k)*fac

                ! Vorticity Field
                wxs(i,jj,k) = omx(i,j,k)*fac
                wys(i,jj,k) = omy(i,j,k)*fac
                wzs(i,jj,k) = omz(i,j,k)*fac

                ! Velocity Gradient
                u11s(i,jj,k) = u11(i,j,k)*fac
                u12s(i,jj,k) = u12(i,j,k)*fac
                u13s(i,jj,k) = u13(i,j,k)*fac
                u21s(i,jj,k) = u21(i,j,k)*fac
                u22s(i,jj,k) = u22(i,j,k)*fac
                u23s(i,jj,k) = u23(i,j,k)*fac
                u31s(i,jj,k) = u31(i,j,k)*fac
                u32s(i,jj,k) = u32(i,j,k)*fac
                u33s(i,jj,k) = u33(i,j,k)*fac

                ! Laplacian
                Lus(i,jj,k) = Lu(i,j,k)*fac
                Lvs(i,jj,k) = Lv(i,j,k)*fac
                Lws(i,jj,k) = Lw(i,j,k)*fac
#IFDEF SCALAR
                ! Scalar and its gradient
                scs(i,jj,k) = scalar(i,j,k)*fac
                cxs(i,jj,k) = sclx(i,j,k)*fac 
                cys(i,jj,k) = scly(i,j,k)*fac
                czs(i,jj,k) = sclz(i,j,k)*fac
#ENDIF
#IFDEF POLYMER
                ! Conformation tensor
                c11s(i,jj,k) = c11(i,j,k)*fac
                c12s(i,jj,k) = c12(i,j,k)*fac
                c13s(i,jj,k) = c13(i,j,k)*fac
                c21s(i,jj,k) = c21(i,j,k)*fac
                c22s(i,jj,k) = c22(i,j,k)*fac
                c23s(i,jj,k) = c23(i,j,k)*fac
                c31s(i,jj,k) = c31(i,j,k)*fac
                c32s(i,jj,k) = c32(i,j,k)*fac
                c33s(i,jj,k) = c33(i,j,k)*fac
               
                ! Conformation tensor gradient
                dc111s(i,jj,k) = dc111(i,j,k)*fac 
                dc112s(i,jj,k) = dc112(i,j,k)*fac 
                dc113s(i,jj,k) = dc113(i,j,k)*fac 
                dc121s(i,jj,k) = dc121(i,j,k)*fac 
                dc122s(i,jj,k) = dc122(i,j,k)*fac 
                dc123s(i,jj,k) = dc123(i,j,k)*fac 
                dc131s(i,jj,k) = dc131(i,j,k)*fac 
                dc132s(i,jj,k) = dc132(i,j,k)*fac 
                dc133s(i,jj,k) = dc133(i,j,k)*fac 

                dc211s(i,jj,k) = dc211(i,j,k)*fac 
                dc212s(i,jj,k) = dc212(i,j,k)*fac 
                dc213s(i,jj,k) = dc213(i,j,k)*fac 
                dc221s(i,jj,k) = dc221(i,j,k)*fac 
                dc222s(i,jj,k) = dc222(i,j,k)*fac 
                dc223s(i,jj,k) = dc223(i,j,k)*fac 
                dc231s(i,jj,k) = dc231(i,j,k)*fac 
                dc232s(i,jj,k) = dc232(i,j,k)*fac 
                dc233s(i,jj,k) = dc233(i,j,k)*fac 

                dc311s(i,jj,k) = dc311(i,j,k)*fac 
                dc312s(i,jj,k) = dc312(i,j,k)*fac 
                dc313s(i,jj,k) = dc313(i,j,k)*fac 
                dc321s(i,jj,k) = dc321(i,j,k)*fac 
                dc322s(i,jj,k) = dc322(i,j,k)*fac 
                dc323s(i,jj,k) = dc323(i,j,k)*fac 
                dc331s(i,jj,k) = dc331(i,j,k)*fac 
                dc332s(i,jj,k) = dc332(i,j,k)*fac 
                dc333s(i,jj,k) = dc333(i,j,k)*fac 
#ENDIF
            end do
        end do
    end do
    !$omp end parallel do

!    write(110+it,*) us
!    write(210+it,*) vs
!    write(310+it,*) ws
    ! NOTE:
    !   The following routines notably use a 1D FFT in each direction,
    !   looped over the indices of the non-transformed variable. FFTW
    !   does have routines to do multiple FFTs like the legacy DNS code,
    !   however, it requires the data to be "shuffled" before and after
    !   each transform, resulting in a total of 3x the memory allocation
    !   and at best a 2x slow-down in the overall transform procedure.
    !   Perhaps there is a way to do this more efficiently, but I'm not
    !   aware of it. - REK 5/15/24

    ! Complex --> Complex Transform (z-direction)
    !$omp parallel do default(shared) private(i,k) 
    do k = 1,nxh
        do i = 1,nyp
            ! Velocity Field
            call fftw_execute_dft(planZb,us(i,:,k),us(i,:,k))
            call fftw_execute_dft(planZb,vs(i,:,k),vs(i,:,k))
            call fftw_execute_dft(planZb,ws(i,:,k),ws(i,:,k))

            ! Vorticity Field
            call fftw_execute_dft(planZb,wxs(i,:,k),wxs(i,:,k))
            call fftw_execute_dft(planZb,wys(i,:,k),wys(i,:,k))
            call fftw_execute_dft(planZb,wzs(i,:,k),wzs(i,:,k))
            
            ! Velocity Gradient
            call fftw_execute_dft(planZb,u11s(i,:,k),u11s(i,:,k))
            call fftw_execute_dft(planZb,u12s(i,:,k),u12s(i,:,k))
            call fftw_execute_dft(planZb,u13s(i,:,k),u13s(i,:,k))
            call fftw_execute_dft(planZb,u21s(i,:,k),u21s(i,:,k))
            call fftw_execute_dft(planZb,u22s(i,:,k),u22s(i,:,k))
            call fftw_execute_dft(planZb,u23s(i,:,k),u23s(i,:,k))
            call fftw_execute_dft(planZb,u31s(i,:,k),u31s(i,:,k))
            call fftw_execute_dft(planZb,u32s(i,:,k),u32s(i,:,k))
            call fftw_execute_dft(planZb,u33s(i,:,k),u33s(i,:,k))

            if (npart .gt. 0) then ! only used in particle tracking
            ! Laplacian
            call fftw_execute_dft(planZb,Lus(i,:,k),Lus(i,:,k))
            call fftw_execute_dft(planZb,Lvs(i,:,k),Lvs(i,:,k))
            call fftw_execute_dft(planZb,Lws(i,:,k),Lws(i,:,k))
            end if
#IFDEF SCALAR
            ! Scalar field and its gradient
            call fftw_execute_dft(planZb,scs(i,:,k),scs(i,:,k))
            call fftw_execute_dft(planZb,cxs(i,:,k),cxs(i,:,k))
            call fftw_execute_dft(planZb,cys(i,:,k),cys(i,:,k))
            call fftw_execute_dft(planZb,czs(i,:,k),czs(i,:,k))
            
#ENDIF
#IFDEF POLYMER
            if (it .ge. (src_start-1)) then
            ! Conformation tensor
            call fftw_execute_dft(planZb,c11s(i,:,k),c11s(i,:,k))
            call fftw_execute_dft(planZb,c12s(i,:,k),c12s(i,:,k))
            call fftw_execute_dft(planZb,c13s(i,:,k),c13s(i,:,k))
            call fftw_execute_dft(planZb,c21s(i,:,k),c21s(i,:,k))
            call fftw_execute_dft(planZb,c22s(i,:,k),c22s(i,:,k))
            call fftw_execute_dft(planZb,c23s(i,:,k),c23s(i,:,k))
            call fftw_execute_dft(planZb,c31s(i,:,k),c31s(i,:,k))
            call fftw_execute_dft(planZb,c32s(i,:,k),c32s(i,:,k))
            call fftw_execute_dft(planZb,c33s(i,:,k),c33s(i,:,k))

            ! Conformation tensor gradient
            call fftw_execute_dft(planZb,dc111s(i,:,k),dc111s(i,:,k))
            call fftw_execute_dft(planZb,dc112s(i,:,k),dc112s(i,:,k))
            call fftw_execute_dft(planZb,dc113s(i,:,k),dc113s(i,:,k))
            call fftw_execute_dft(planZb,dc121s(i,:,k),dc121s(i,:,k))
            call fftw_execute_dft(planZb,dc122s(i,:,k),dc122s(i,:,k))
            call fftw_execute_dft(planZb,dc123s(i,:,k),dc123s(i,:,k))
            call fftw_execute_dft(planZb,dc131s(i,:,k),dc131s(i,:,k))
            call fftw_execute_dft(planZb,dc132s(i,:,k),dc132s(i,:,k))
            call fftw_execute_dft(planZb,dc133s(i,:,k),dc133s(i,:,k))

            call fftw_execute_dft(planZb,dc211s(i,:,k),dc211s(i,:,k))
            call fftw_execute_dft(planZb,dc212s(i,:,k),dc212s(i,:,k))
            call fftw_execute_dft(planZb,dc213s(i,:,k),dc213s(i,:,k))
            call fftw_execute_dft(planZb,dc221s(i,:,k),dc221s(i,:,k))
            call fftw_execute_dft(planZb,dc222s(i,:,k),dc222s(i,:,k))
            call fftw_execute_dft(planZb,dc223s(i,:,k),dc223s(i,:,k))
            call fftw_execute_dft(planZb,dc231s(i,:,k),dc231s(i,:,k))
            call fftw_execute_dft(planZb,dc232s(i,:,k),dc232s(i,:,k))
            call fftw_execute_dft(planZb,dc233s(i,:,k),dc233s(i,:,k))

            call fftw_execute_dft(planZb,dc311s(i,:,k),dc311s(i,:,k))
            call fftw_execute_dft(planZb,dc312s(i,:,k),dc312s(i,:,k))
            call fftw_execute_dft(planZb,dc313s(i,:,k),dc313s(i,:,k))
            call fftw_execute_dft(planZb,dc321s(i,:,k),dc321s(i,:,k))
            call fftw_execute_dft(planZb,dc322s(i,:,k),dc322s(i,:,k))
            call fftw_execute_dft(planZb,dc323s(i,:,k),dc323s(i,:,k))
            call fftw_execute_dft(planZb,dc331s(i,:,k),dc331s(i,:,k))
            call fftw_execute_dft(planZb,dc332s(i,:,k),dc332s(i,:,k))
            call fftw_execute_dft(planZb,dc333s(i,:,k),dc333s(i,:,k))
            end if
#ENDIF
        end do
    end do
    !$omp end parallel do

!    write(120+it,*) us
!    write(220+it,*) vs
!    write(320+it,*) ws

    ! Complex --> Real Transform (x-direction)
    !$omp parallel do default(shared) private(i,j) 
    do j = 1,mz
        do i = 1,nyp
            ! Velocity Field
            call fftw_execute_dft_c2r(planXb,us(i,j,:),up(i,j,:))
            call fftw_execute_dft_c2r(planXb,vs(i,j,:),vp(i,j,:))
            call fftw_execute_dft_c2r(planXb,ws(i,j,:),wp(i,j,:))

            ! Vorticity Field
            call fftw_execute_dft_c2r(planXb,wxs(i,j,:),wxp(i,j,:))
            call fftw_execute_dft_c2r(planXb,wys(i,j,:),wyp(i,j,:))
            call fftw_execute_dft_c2r(planXb,wzs(i,j,:),wzp(i,j,:))
            
            ! Velocity Gradient
            call fftw_execute_dft_c2r(planXb,u11s(i,j,:),u11p(i,j,:))
            call fftw_execute_dft_c2r(planXb,u12s(i,j,:),u12p(i,j,:))
            call fftw_execute_dft_c2r(planXb,u13s(i,j,:),u13p(i,j,:))
            call fftw_execute_dft_c2r(planXb,u21s(i,j,:),u21p(i,j,:))
            call fftw_execute_dft_c2r(planXb,u22s(i,j,:),u22p(i,j,:))
            call fftw_execute_dft_c2r(planXb,u23s(i,j,:),u23p(i,j,:))
            call fftw_execute_dft_c2r(planXb,u31s(i,j,:),u31p(i,j,:))
            call fftw_execute_dft_c2r(planXb,u32s(i,j,:),u32p(i,j,:))
            call fftw_execute_dft_c2r(planXb,u33s(i,j,:),u33p(i,j,:))

            if (npart .gt. 0) then ! only used in particle tracking
            ! Laplacian
            call fftw_execute_dft_c2r(planXb,Lus(i,j,:),Lup(i,j,:))
            call fftw_execute_dft_c2r(planXb,Lvs(i,j,:),Lvp(i,j,:))
            call fftw_execute_dft_c2r(planXb,Lws(i,j,:),Lwp(i,j,:))
            end if
#IFDEF SCALAR
            ! Scalar field and its gradient
            call fftw_execute_dft_c2r(planXb,scs(i,j,:),scp(i,j,:))
            call fftw_execute_dft_c2r(planXb,cxs(i,j,:),cxp(i,j,:))
            call fftw_execute_dft_c2r(planXb,cys(i,j,:),cyp(i,j,:))
            call fftw_execute_dft_c2r(planXb,czs(i,j,:),czp(i,j,:))
            
#ENDIF
#IFDEF POLYMER
            if (it .ge. (src_start-1)) then
            ! Conformation tensor
            call fftw_execute_dft_c2r(planXb,c11s(i,j,:),c11p(i,j,:))
            call fftw_execute_dft_c2r(planXb,c12s(i,j,:),c12p(i,j,:))
            call fftw_execute_dft_c2r(planXb,c13s(i,j,:),c13p(i,j,:))
            call fftw_execute_dft_c2r(planXb,c21s(i,j,:),c21p(i,j,:))
            call fftw_execute_dft_c2r(planXb,c22s(i,j,:),c22p(i,j,:))
            call fftw_execute_dft_c2r(planXb,c23s(i,j,:),c23p(i,j,:))
            call fftw_execute_dft_c2r(planXb,c31s(i,j,:),c31p(i,j,:))
            call fftw_execute_dft_c2r(planXb,c32s(i,j,:),c32p(i,j,:))
            call fftw_execute_dft_c2r(planXb,c33s(i,j,:),c33p(i,j,:))

            ! Conformation tensor gradient
            call fftw_execute_dft_c2r(planXb,dc111s(i,j,:),dc111p(i,j,:))
            call fftw_execute_dft_c2r(planXb,dc112s(i,j,:),dc112p(i,j,:))
            call fftw_execute_dft_c2r(planXb,dc113s(i,j,:),dc113p(i,j,:))
            call fftw_execute_dft_c2r(planXb,dc121s(i,j,:),dc121p(i,j,:))
            call fftw_execute_dft_c2r(planXb,dc122s(i,j,:),dc122p(i,j,:))
            call fftw_execute_dft_c2r(planXb,dc123s(i,j,:),dc123p(i,j,:))
            call fftw_execute_dft_c2r(planXb,dc131s(i,j,:),dc131p(i,j,:))
            call fftw_execute_dft_c2r(planXb,dc132s(i,j,:),dc132p(i,j,:))
            call fftw_execute_dft_c2r(planXb,dc133s(i,j,:),dc133p(i,j,:))

            call fftw_execute_dft_c2r(planXb,dc211s(i,j,:),dc211p(i,j,:))
            call fftw_execute_dft_c2r(planXb,dc212s(i,j,:),dc212p(i,j,:))
            call fftw_execute_dft_c2r(planXb,dc213s(i,j,:),dc213p(i,j,:))
            call fftw_execute_dft_c2r(planXb,dc221s(i,j,:),dc221p(i,j,:))
            call fftw_execute_dft_c2r(planXb,dc222s(i,j,:),dc222p(i,j,:))
            call fftw_execute_dft_c2r(planXb,dc223s(i,j,:),dc223p(i,j,:))
            call fftw_execute_dft_c2r(planXb,dc231s(i,j,:),dc231p(i,j,:))
            call fftw_execute_dft_c2r(planXb,dc232s(i,j,:),dc232p(i,j,:))
            call fftw_execute_dft_c2r(planXb,dc233s(i,j,:),dc233p(i,j,:))

            call fftw_execute_dft_c2r(planXb,dc311s(i,j,:),dc311p(i,j,:))
            call fftw_execute_dft_c2r(planXb,dc312s(i,j,:),dc312p(i,j,:))
            call fftw_execute_dft_c2r(planXb,dc313s(i,j,:),dc313p(i,j,:))
            call fftw_execute_dft_c2r(planXb,dc321s(i,j,:),dc321p(i,j,:))
            call fftw_execute_dft_c2r(planXb,dc322s(i,j,:),dc322p(i,j,:))
            call fftw_execute_dft_c2r(planXb,dc323s(i,j,:),dc323p(i,j,:))
            call fftw_execute_dft_c2r(planXb,dc331s(i,j,:),dc331p(i,j,:))
            call fftw_execute_dft_c2r(planXb,dc332s(i,j,:),dc332p(i,j,:))
            call fftw_execute_dft_c2r(planXb,dc333s(i,j,:),dc333p(i,j,:))
            end if
#ENDIF
        end do
    end do
    !$omp end parallel do

!    write(130+it,*) up
!    write(230+it,*) vp
!    write(330+it,*) wp
    ! Initialize the force field array for IBF
    dragx = 0.0
    dragy = 0.0
    dragz = 0.0

    ! Real --> Real Transform (y-transform) Uses DCT-I to transform to/from cosine grid

    ! These loops also compute the nonlinear terms and then begin the process of
    ! transforming the nonlinear output arrays to spectral space
    !$omp parallel do default(shared) private(i,j,k,cflcheck,wdes,zj,rsq,argrad,fr,xi,argx,fx,ii,ipii,idif,jdif,segdrag,jpjj,jj) schedule(dynamic)
    do k = 1,mx
        do j = 1,mz
            ! Velocity Field
            call fftw_execute_r2r(planY,up(:,j,k),up(:,j,k))
            call fftw_execute_r2r(planY,vp(:,j,k),vp(:,j,k))
            call fftw_execute_r2r(planY,wp(:,j,k),wp(:,j,k))

            ! Vorticity Field
            call fftw_execute_r2r(planY,wxp(:,j,k),wxp(:,j,k))
            call fftw_execute_r2r(planY,wyp(:,j,k),wyp(:,j,k))
            call fftw_execute_r2r(planY,wzp(:,j,k),wzp(:,j,k))
            
            ! Velocity Gradient
            call fftw_execute_r2r(planY,u11p(:,j,k),u11p(:,j,k))
            call fftw_execute_r2r(planY,u12p(:,j,k),u12p(:,j,k))
            call fftw_execute_r2r(planY,u13p(:,j,k),u13p(:,j,k))
            call fftw_execute_r2r(planY,u21p(:,j,k),u21p(:,j,k))
            call fftw_execute_r2r(planY,u22p(:,j,k),u22p(:,j,k))
            call fftw_execute_r2r(planY,u23p(:,j,k),u23p(:,j,k))
            call fftw_execute_r2r(planY,u31p(:,j,k),u31p(:,j,k))
            call fftw_execute_r2r(planY,u32p(:,j,k),u32p(:,j,k))
            call fftw_execute_r2r(planY,u33p(:,j,k),u33p(:,j,k))

            if (npart .gt. 0) then ! only used in particle tracking
            ! Laplacian
            call fftw_execute_r2r(planY,Lup(:,j,k),Lup(:,j,k))
            call fftw_execute_r2r(planY,Lvp(:,j,k),Lvp(:,j,k))
            call fftw_execute_r2r(planY,Lwp(:,j,k),Lwp(:,j,k))
            end if
#IFDEF SCALAR
            ! Scalar field and its gradient
            call fftw_execute_r2r(planY,scp(:,j,k),scp(:,j,k))
            call fftw_execute_r2r(planY,cxp(:,j,k),cxp(:,j,k))
            call fftw_execute_r2r(planY,cyp(:,j,k),cyp(:,j,k))
            call fftw_execute_r2r(planY,czp(:,j,k),czp(:,j,k))
#ENDIF
#IFDEF POLYMER
            if (it .ge. (src_start-1)) then
            ! Conformation tensor
            call fftw_execute_r2r(planY,c11p(:,j,k),c11p(:,j,k))
            call fftw_execute_r2r(planY,c12p(:,j,k),c12p(:,j,k))
            call fftw_execute_r2r(planY,c13p(:,j,k),c13p(:,j,k))
            call fftw_execute_r2r(planY,c21p(:,j,k),c21p(:,j,k))
            call fftw_execute_r2r(planY,c22p(:,j,k),c22p(:,j,k))
            call fftw_execute_r2r(planY,c23p(:,j,k),c23p(:,j,k))
            call fftw_execute_r2r(planY,c31p(:,j,k),c31p(:,j,k))
            call fftw_execute_r2r(planY,c32p(:,j,k),c32p(:,j,k))
            call fftw_execute_r2r(planY,c33p(:,j,k),c33p(:,j,k))

            ! Conformation tensor gradYnt
            call fftw_execute_r2r(planY,dc111p(:,j,k),dc111p(:,j,k))
            call fftw_execute_r2r(planY,dc112p(:,j,k),dc112p(:,j,k))
            call fftw_execute_r2r(planY,dc113p(:,j,k),dc113p(:,j,k))
            call fftw_execute_r2r(planY,dc121p(:,j,k),dc121p(:,j,k))
            call fftw_execute_r2r(planY,dc122p(:,j,k),dc122p(:,j,k))
            call fftw_execute_r2r(planY,dc123p(:,j,k),dc123p(:,j,k))
            call fftw_execute_r2r(planY,dc131p(:,j,k),dc131p(:,j,k))
            call fftw_execute_r2r(planY,dc132p(:,j,k),dc132p(:,j,k))
            call fftw_execute_r2r(planY,dc133p(:,j,k),dc133p(:,j,k))

            call fftw_execute_r2r(planY,dc211p(:,j,k),dc211p(:,j,k))
            call fftw_execute_r2r(planY,dc212p(:,j,k),dc212p(:,j,k))
            call fftw_execute_r2r(planY,dc213p(:,j,k),dc213p(:,j,k))
            call fftw_execute_r2r(planY,dc221p(:,j,k),dc221p(:,j,k))
            call fftw_execute_r2r(planY,dc222p(:,j,k),dc222p(:,j,k))
            call fftw_execute_r2r(planY,dc223p(:,j,k),dc223p(:,j,k))
            call fftw_execute_r2r(planY,dc231p(:,j,k),dc231p(:,j,k))
            call fftw_execute_r2r(planY,dc232p(:,j,k),dc232p(:,j,k))
            call fftw_execute_r2r(planY,dc233p(:,j,k),dc233p(:,j,k))

            call fftw_execute_r2r(planY,dc311p(:,j,k),dc311p(:,j,k))
            call fftw_execute_r2r(planY,dc312p(:,j,k),dc312p(:,j,k))
            call fftw_execute_r2r(planY,dc313p(:,j,k),dc313p(:,j,k))
            call fftw_execute_r2r(planY,dc321p(:,j,k),dc321p(:,j,k))
            call fftw_execute_r2r(planY,dc322p(:,j,k),dc322p(:,j,k))
            call fftw_execute_r2r(planY,dc323p(:,j,k),dc323p(:,j,k))
            call fftw_execute_r2r(planY,dc331p(:,j,k),dc331p(:,j,k))
            call fftw_execute_r2r(planY,dc332p(:,j,k),dc332p(:,j,k))
            call fftw_execute_r2r(planY,dc333p(:,j,k),dc333p(:,j,k))
            end if
#ENDIF

            ! At this point, each y-column of the variables is appropriately
            ! transformed to 3/2 physical space and can be used to compute
            ! the nonlinear terms

            !---------------------------------------------------------------------!
            !        Compute the cross product of velocity and vorticity          !
            !---------------------------------------------------------------------!
   
            do i = 1,nyp 
                vwx(i,j,k) =  vp(i,j,k)*wzp(i,j,k) - wp(i,j,k)*wyp(i,j,k)
                vwy(i,j,k) = -up(i,j,k)*wzp(i,j,k) + wp(i,j,k)*wxp(i,j,k)
                vwz(i,j,k) =  up(i,j,k)*wyp(i,j,k) - vp(i,j,k)*wxp(i,j,k)
#IFDEF SCALAR
                 vc(i,j,k) = -(up(i,j,k)*cxp(i,j,k) + vp(i,j,k)*cyp(i,j,k) + wp(i,j,k)*czp(i,j,k))
#ENDIF
#IFDEF POLYMER    
                if (it .ge. src_start-1) then
                c11np(i,j,k)=c11p(i,j,k)*u11p(i,j,k)+c12p(i,j,k)*u12p(i,j,k)+c13p(i,j,k)*u13p(i,j,k)+  &
                             c11p(i,j,k)*u11p(i,j,k)+c21p(i,j,k)*u12p(i,j,k)+c31p(i,j,k)*u13p(i,j,k)-  &
                            (up(i,j,k)*dc111p(i,j,k)+vp(i,j,k)*dc112p(i,j,k)+wp(i,j,k)*dc113p(i,j,k))            
                  
                c12np(i,j,k)=c11p(i,j,k)*u21p(i,j,k)+c12p(i,j,k)*u22p(i,j,k)+c13p(i,j,k)*u23p(i,j,k)+   &
                             c12p(i,j,k)*u11p(i,j,k)+c22p(i,j,k)*u12p(i,j,k)+c32p(i,j,k)*u13p(i,j,k)-   &
                            (up(i,j,k)*dc121p(i,j,k)+vp(i,j,k)*dc122p(i,j,k)+wp(i,j,k)*dc123p(i,j,k))
                  
                c13np(i,j,k)=c11p(i,j,k)*u31p(i,j,k)+c12p(i,j,k)*u32p(i,j,k)+c13p(i,j,k)*u33p(i,j,k)+   &
                             c13p(i,j,k)*u11p(i,j,k)+c23p(i,j,k)*u12p(i,j,k)+c33p(i,j,k)*u13p(i,j,k)-   &
                             (up(i,j,k)*dc131p(i,j,k)+vp(i,j,k)*dc132p(i,j,k)+wp(i,j,k)*dc133p(i,j,k))
                  
                c22np(i,j,k)=c21p(i,j,k)*u21p(i,j,k)+c22p(i,j,k)*u22p(i,j,k)+c23p(i,j,k)*u23p(i,j,k)+   &
                             c12p(i,j,k)*u21p(i,j,k)+c22p(i,j,k)*u22p(i,j,k)+c32p(i,j,k)*u23p(i,j,k)-   &
                             (up(i,j,k)*dc221p(i,j,k)+vp(i,j,k)*dc222p(i,j,k)+wp(i,j,k)*dc223p(i,j,k))
                  
                c23np(i,j,k)=c21p(i,j,k)*u31p(i,j,k)+c22p(i,j,k)*u32p(i,j,k)+c23p(i,j,k)*u33p(i,j,k)+   &
                             c13p(i,j,k)*u21p(i,j,k)+c23p(i,j,k)*u22p(i,j,k)+c33p(i,j,k)*u23p(i,j,k)-   &
                             (up(i,j,k)*dc231p(i,j,k)+vp(i,j,k)*dc232p(i,j,k)+wp(i,j,k)*dc233p(i,j,k))

                c33np(i,j,k)=c31p(i,j,k)*u31p(i,j,k)+c32p(i,j,k)*u32p(i,j,k)+c33p(i,j,k)*u33p(i,j,k)+   &         
                             c13p(i,j,k)*u31p(i,j,k)+c23p(i,j,k)*u32p(i,j,k)+c33p(i,j,k)*u33p(i,j,k)-   &
                            (up(i,j,k)*dc331p(i,j,k)+vp(i,j,k)*dc332p(i,j,k)+wp(i,j,k)*dc333p(i,j,k))
     
                trp(i,j,k) = c11p(i,j,k) + c22p(i,j,k) + c33p(i,j,k)
     
                if (ipeter .eq. 0) then  ! peterlin function is 1.0
    
                    str11np(i,j,k)= c11p(i,j,k)
                    str12np(i,j,k)= c12p(i,j,k)
                    str13np(i,j,k)= c13p(i,j,k)
                    str22np(i,j,k)= c22p(i,j,k)
                    str23np(i,j,k)= c23p(i,j,k)
                    str33np(i,j,k)= c33p(i,j,k)
              
                else
               
                    str11np(i,j,k)=((zlmax**2 - 3.0)/(zlmax**2 - trp(i,j,k)))*c11p(i,j,k)
                    str12np(i,j,k)=((zlmax**2 - 3.0)/(zlmax**2 - trp(i,j,k)))*c12p(i,j,k)
                    str13np(i,j,k)=((zlmax**2 - 3.0)/(zlmax**2 - trp(i,j,k)))*c13p(i,j,k)
                    str22np(i,j,k)=((zlmax**2 - 3.0)/(zlmax**2 - trp(i,j,k)))*c22p(i,j,k)
                    str23np(i,j,k)=((zlmax**2 - 3.0)/(zlmax**2 - trp(i,j,k)))*c23p(i,j,k)
                    str33np(i,j,k)=((zlmax**2 - 3.0)/(zlmax**2 - trp(i,j,k)))*c33p(i,j,k)
               
                end if
               
                ! Add brownian motion terms
                str11np(i,j,k)=str11np(i,j,k)-1.0
                str22np(i,j,k)=str22np(i,j,k)-1.0
                str33np(i,j,k)=str33np(i,j,k)-1.0

                ! Polymer model
                if (itarget .eq. 0) then ! polymer ocean
                    beta_poly(i,j,k) = qbeta
                else if (itarget .eq. 1) then
                    ! Nonlinear model:
                    ! beta = exp(-alpha*gamma) --> beta_poly
                    ! alpha = 2.6e-03 PPM --> alpha_poly
                    ! gamma = scalar concentration (PPM) --> scp

                    beta_poly(i,j,k) = exp(-alpha_poly*abs(scp(i,j,k)))

                else if (itarget .eq. 2) then ! Linear polymer model
                    ! Linear model:
                    ! beta = (alpha*|gamma|)
                    beta_poly(i,j,k) = 1.0/(alpha_poly*abs(scp(i,j,k)) + 1.0)
                end if

                zbeta1 = (1.0 - beta_poly(i,j,k))/(re*beta_poly(i,j,k)*tpoly) ! = (nu_0 - nu_s)

                qp11np(i,j,k) = zbeta1*str11np(i,j,k)
                qp12np(i,j,k) = zbeta1*str12np(i,j,k)
                qp13np(i,j,k) = zbeta1*str13np(i,j,k)
                qp22np(i,j,k) = zbeta1*str22np(i,j,k)
                qp23np(i,j,k) = zbeta1*str23np(i,j,k)
                qp33np(i,j,k) = zbeta1*str33np(i,j,k)

                else
                    beta_poly(i,j,k) = 1.0 ! for printing
                end if ! it >= src_start-1
#ENDIF
                cflcheck = (abs(up(i,j,k))/delxm + abs(vp(i,j,k))/seght(i) + abs(wp(i,j,k))/delzm)*dt
                if(cflcheck .gt. cfl(i)) cfl(i) = cflcheck

                !---------------------------------------------------------------------!
                ! Compute immersed boundary force terms and add to nonlinear term     !
                ! The type of forcing is based on the value of imatrix at the grid    !
                ! point. Several of these cases probably require more coding and      !
                ! specification in the geometry file, but they have been unused since !
                ! I've used the code                                                  !
                !---------------------------------------------------------------------!
    
                if (imatrix(i,j,k) .eq. 1 .or. (imatrix(i,j,k) .eq. 4 .and. it .le. perturbtime)) then ! Generate solid surfaces
                    fxintg(i,j,k) = fxintg(i,j,k) + up(i,j,k)
                    fyintg(i,j,k) = fyintg(i,j,k) + vp(i,j,k)
                    fzintg(i,j,k) = fzintg(i,j,k) + wp(i,j,k)
                    dragx(i,j,k) = (-ugain*up(i,j,k)) - gain*fxintg(i,j,k)
                    dragy(i,j,k) = (-ugain*vp(i,j,k)) - gain*fyintg(i,j,k)
                    dragz(i,j,k) = (-ugain*wp(i,j,k)) - gain*fzintg(i,j,k)
        
                else if (imatrix(i,j,k) .eq. 3) then ! Spanwise-damping textures
                    fzintg(i,j,k) = fzintg(i,j,k) + wp(i,j,k)
                    dragz(i,j,k) = -ugain*wp(i,j,k) - gain*fzintg(i,j,k)
        
                else if (imatrix(i,j,k) .eq. 8) then ! Spanwise moving wall
                    if (k .le. (bfwidth*3 + 1)) then
                        wdes = (0.3*Uinf)/slopelength*(k - (bfwidth*3 + 1 - slopelength))
                    else
                        wdes = 0.3*Uinf
                    end if
    
                    fxintg(i,j,k) = fxintg(i,j,k) + up(i,j,k)
                    fyintg(i,j,k) = fyintg(i,j,k) + vp(i,j,k)
                    fzintg(i,j,k) = fzintg(i,j,k) + wp(i,j,k)-wdes
                    dragx(i,j,k)    = -ugain*up(i,j,k) - gain*fxintg(i,j,k)
                    dragy(i,j,k)    = -ugain*vp(i,j,k) - gain*fyintg(i,j,k)
                    dragz(i,j,k)    = -ugain*(wp(i,j,k)-wdes) - gain*fzintg(i,j,k)
                    
                else if (imatrix(i,j,k) .eq. 6) then ! Buffer zone
                    ! NOTE: initu and initw should technically be transformed first, but since
                    ! the buffer region has only been used for x- and z-constant flows, using
                    ! them as they are should be fine
                    fxintg(i,j,k) = fxintg(i,j,k) + up(i,j,k) - initu(i,j,k)
                    fyintg(i,j,k) = fyintg(i,j,k) + vp(i,j,k) - initv(i,j,k)
                    fzintg(i,j,k) = fzintg(i,j,k) + wp(i,j,k) - initw(i,j,k)
                    vwx(i,j,k) = vwx(i,j,k) + (-bfugain(k-bfhead+1)*(up(i,j,k) - initu(i,j,k)) - (bfgain(k-bfhead+1) * fxintg(i,j,k)))
                    vwy(i,j,k) = vwy(i,j,k) + (-bfugain(k-bfhead+1)*(vp(i,j,k) - initv(i,j,k)) - (bfgain(k-bfhead+1) * fyintg(i,j,k)))
                    vwz(i,j,k) = vwz(i,j,k) + (-bfugain(k-bfhead+1)*(wp(i,j,k) - initw(i,j,k)) - (bfgain(k-bfhead+1) * fzintg(i,j,k)))
    
                else if (imatrix(i,j,k) .eq. 7) then ! Suction region
                    fyintg(i,j,k) = fyintg(i,j,k) + (vp(i,j,k)-vdes(i))
                    vwy(i,j,k) = vwy(i,j,k) + (-gain*fyintg(i,j,k)) + (-ugain*(vp(i,j,k)-vdes(i)))
                    fxintg(i,j,k) = fxintg(i,j,k) + (up(i,j,k)-uinf)
                    vwx(i,j,k) = vwx(i,j,k) + (-gain*fxintg(i,j,k)) + (-ugain*(up(i,j,k)-uinf))
                    fzintg(i,j,k) = fzintg(i,j,k) + (wp(i,j,k)-initw(i,j,k))
                    vwz(i,j,k) = vwz(i,j,k) + (-gain*fzintg(i,j,k)) + (-ugain*(wp(i,j,k)-initw(i,j,k)))
     
                else if (imatrix(i,j,k) .eq. 11) then ! Constant x-force between wall and suction layer
                    dragx(i,j,k) = 1.0
    
                else if (imatrix(i,j,k) .eq. 0 .and. flow_select .eq. 1) then ! Constant pressure gradient
                    vwx(i,j,k) = vwx(i,j,k) - dPdx 
                
                else if (flow_select .eq. 5 .and. it .le. perturbtime) then ! Vortex ring
                    zj = float(j-1)*delzm
                    rsq = (ycoord(i) - ycenter)**2 + (zj - zcenter)**2
                    argrad = forbeta*(rad - sqrt(rsq))
                    fr = 0.5*(1.0 + tanh(argrad))
    
                    xi = float(k-1)*delxm
                    argx = forbeta*(L - abs(xi - xcenter))
                    fx = 0.5*(1.0 + tanh(argx))
  
                    vwx(i,j,k) = vwx(i,j,k) + bdyfx*fx*fr
    
                end if
#IFDEF SCALAR
                if (scl_flag .ge. 2) then
                    vc(i,j,k) = vc(i,j,k) + scsource(i,j,k)/dt
                end if
#ENDIF

        !---------------------------------------------------------------------!
        !            Apply the force field to the solid surface               !
        !---------------------------------------------------------------------!
                if (i .le. kmaxsurf) then
                do ii = -2, 2
                    ipii = k + ii
                    idif = 1 + iabs(ii)        
                    if(ipii .lt. 1 ) ipii = ipii + mx         
                    if(ipii .gt. mx) ipii = ipii - mx
                    do jj = -2, 2
                        jdif = 1 + iabs(jj)
                        segdrag = fspread(jdif,idif)
        
                        if (j .ge. 3 .or. j .le. mz-2) then 
                            jpjj = j + jj
                            xsegdrag = segdrag*dragx(i,j,k)
                            ysegdrag = segdrag*dragy(i,j,k)
                            zsegdrag = segdrag*dragz(i,j,k)
                            vwx(i,jpjj,ipii) = vwx(i,jpjj,ipii) + xsegdrag
                            vwy(i,jpjj,ipii) = vwy(i,jpjj,ipii) + ysegdrag
                            vwz(i,jpjj,ipii) = vwz(i,jpjj,ipii) + zsegdrag
                        else if (j .eq. 1 .or. j .eq. 2) then 
                            jpjj = j + jj
                            if(jpjj .lt. 1 ) jpjj = jpjj + mz
                            xsegdrag = segdrag*dragx(i,j,k)
                            ysegdrag = segdrag*dragy(i,j,k)
                            zsegdrag = segdrag*dragz(i,j,k)
                            vwx(i,jpjj,ipii) = vwx(i,jpjj,ipii) + xsegdrag
                            vwy(i,jpjj,ipii) = vwy(i,jpjj,ipii) + ysegdrag
                            vwz(i,jpjj,ipii) = vwz(i,jpjj,ipii) + zsegdrag
                        else 
                            jpjj = j + jj
                            if(jpjj .gt. mz) jpjj = jpjj - mz
                            xsegdrag = segdrag*dragx(i,j,k)
                            ysegdrag = segdrag*dragy(i,j,k)
                            zsegdrag = segdrag*dragz(i,j,k)
                            vwx(i,jpjj,ipii) = vwx(i,jpjj,ipii) + xsegdrag
                            vwy(i,jpjj,ipii) = vwy(i,jpjj,ipii) + ysegdrag
                            vwz(i,jpjj,ipii) = vwz(i,jpjj,ipii) + zsegdrag
                        end if
                    end do ! jj
                end do ! ii
                end if

            end do ! y-planes


        !---------------------------------------------------------------------!
        !              Now transform vxw back to spectral space               !
        !---------------------------------------------------------------------!
        ! We continue in the same loop to do the DCT-I again (it is its own inverse)    
        ! But this time we need to normalize by 2*ny

            ! velocity x vorticity
            call fftw_execute_r2r(planY,vwx(:,j,k),vwx(:,j,k))
            call fftw_execute_r2r(planY,vwy(:,j,k),vwy(:,j,k))
            call fftw_execute_r2r(planY,vwz(:,j,k),vwz(:,j,k))
#IFDEF SCALAR
            ! velocity x scalar
            call fftw_execute_r2r(planY,vc(:,j,k),vc(:,j,k))
#ENDIF
#IFDEF POLYMER
            if (it .ge. (src_start-1)) then
            ! Conformation tensor
            call fftw_execute_r2r(planY,c11np(:,j,k),c11np(:,j,k))
            call fftw_execute_r2r(planY,c12np(:,j,k),c12np(:,j,k))
            call fftw_execute_r2r(planY,c13np(:,j,k),c13np(:,j,k))
            call fftw_execute_r2r(planY,c22np(:,j,k),c22np(:,j,k))
            call fftw_execute_r2r(planY,c23np(:,j,k),c23np(:,j,k))
            call fftw_execute_r2r(planY,c33np(:,j,k),c33np(:,j,k))

            ! Polymer stress
            call fftw_execute_r2r(planY,str11np(:,j,k),str11np(:,j,k))
            call fftw_execute_r2r(planY,str12np(:,j,k),str12np(:,j,k))
            call fftw_execute_r2r(planY,str13np(:,j,k),str13np(:,j,k))
            call fftw_execute_r2r(planY,str22np(:,j,k),str22np(:,j,k))
            call fftw_execute_r2r(planY,str23np(:,j,k),str23np(:,j,k))
            call fftw_execute_r2r(planY,str33np(:,j,k),str33np(:,j,k))

            call fftw_execute_r2r(planY,qp11np(:,j,k),qp11np(:,j,k))
            call fftw_execute_r2r(planY,qp12np(:,j,k),qp12np(:,j,k))
            call fftw_execute_r2r(planY,qp13np(:,j,k),qp13np(:,j,k))
            call fftw_execute_r2r(planY,qp22np(:,j,k),qp22np(:,j,k))
            call fftw_execute_r2r(planY,qp23np(:,j,k),qp23np(:,j,k))
            call fftw_execute_r2r(planY,qp33np(:,j,k),qp33np(:,j,k))
            end if
#ENDIF
        end do ! j
    end do ! k
    !$omp end parallel do

!    write(140+it,*) up
!    write(240+it,*) vp
!    write(340+it,*) wp

    ! Normalize then convert to Chebyshev modes
    vwx = vwx/float(ny)
    vwy = vwy/float(ny)
    vwz = vwz/float(ny)

    vwx(1,:,:) = vwx(1,:,:)/2.0
    vwy(1,:,:) = vwy(1,:,:)/2.0
    vwz(1,:,:) = vwz(1,:,:)/2.0
    vwx(nyp,:,:) = vwx(nyp,:,:)/2.0
    vwy(nyp,:,:) = vwy(nyp,:,:)/2.0
    vwz(nyp,:,:) = vwz(nyp,:,:)/2.0
#IFDEF SCALAR
    vc  =  vc/float(ny)

    vc(1,:,:) = vc(1,:,:)/2.0
#ENDIF
#IFDEF POLYMER  
    if (it .ge. (src_start-1)) then
    c11np = c11np/float(ny)
    c12np = c12np/float(ny)
    c13np = c13np/float(ny)
    c21np = c21np/float(ny)
    c22np = c22np/float(ny)
    c33np = c33np/float(ny)

    str11np = str11np/float(ny)
    str12np = str12np/float(ny)
    str13np = str13np/float(ny)
    str21np = str21np/float(ny)
    str22np = str22np/float(ny)
    str33np = str33np/float(ny)

    qp11np = qp11np/float(ny)
    qp12np = qp12np/float(ny)
    qp13np = qp13np/float(ny)
    qp21np = qp21np/float(ny)
    qp22np = qp22np/float(ny)
    qp33np = qp33np/float(ny)

    c11np(1,:,:) = c11np(1,:,:)/2.0
    c12np(1,:,:) = c12np(1,:,:)/2.0
    c13np(1,:,:) = c13np(1,:,:)/2.0
    c22np(1,:,:) = c22np(1,:,:)/2.0
    c23np(1,:,:) = c23np(1,:,:)/2.0
    c33np(1,:,:) = c33np(1,:,:)/2.0
    c11np(nyp,:,:) = c11np(nyp,:,:)/2.0
    c12np(nyp,:,:) = c12np(nyp,:,:)/2.0
    c13np(nyp,:,:) = c13np(nyp,:,:)/2.0
    c22np(nyp,:,:) = c22np(nyp,:,:)/2.0
    c23np(nyp,:,:) = c23np(nyp,:,:)/2.0
    c33np(nyp,:,:) = c33np(nyp,:,:)/2.0

    str11np(1,:,:)   = str11np(1,:,:)/2.0
    str12np(1,:,:)   = str12np(1,:,:)/2.0
    str13np(1,:,:)   = str13np(1,:,:)/2.0
    str22np(1,:,:)   = str22np(1,:,:)/2.0
    str23np(1,:,:)   = str23np(1,:,:)/2.0
    str33np(1,:,:)   = str33np(1,:,:)/2.0
    str11np(nyp,:,:) = str11np(nyp,:,:)/2.0
    str12np(nyp,:,:) = str12np(nyp,:,:)/2.0
    str13np(nyp,:,:) = str13np(nyp,:,:)/2.0
    str22np(nyp,:,:) = str22np(nyp,:,:)/2.0
    str23np(nyp,:,:) = str23np(nyp,:,:)/2.0
    str33np(nyp,:,:) = str33np(nyp,:,:)/2.0

    qp11np(1,:,:)   = qp11np(1,:,:)/2.0
    qp12np(1,:,:)   = qp12np(1,:,:)/2.0
    qp13np(1,:,:)   = qp13np(1,:,:)/2.0
    qp22np(1,:,:)   = qp22np(1,:,:)/2.0
    qp23np(1,:,:)   = qp23np(1,:,:)/2.0
    qp33np(1,:,:)   = qp33np(1,:,:)/2.0
    qp11np(nyp,:,:) = qp11np(nyp,:,:)/2.0
    qp12np(nyp,:,:) = qp12np(nyp,:,:)/2.0
    qp13np(nyp,:,:) = qp13np(nyp,:,:)/2.0
    qp22np(nyp,:,:) = qp22np(nyp,:,:)/2.0
    qp23np(nyp,:,:) = qp23np(nyp,:,:)/2.0
    qp33np(nyp,:,:) = qp33np(nyp,:,:)/2.0
    end if
#ENDIF
    
    ! Real --> Complex Transform (x-direction)
    !$omp parallel do default(shared) private(i,j)
    do j = 1,mz
        do i = 1,nyp
            call fftw_execute_dft_r2c(planXf,vwx(i,j,:),vwxs(i,j,:))
            call fftw_execute_dft_r2c(planXf,vwy(i,j,:),vwys(i,j,:))
            call fftw_execute_dft_r2c(planXf,vwz(i,j,:),vwzs(i,j,:))
#IFDEF SCALAR
            call fftw_execute_dft_r2c(planXf,vc(i,j,:),vcs(i,j,:))
#ENDIF
#IFDEF POLYMER
            if (it .ge. (src_start-1)) then
            call fftw_execute_dft_r2c(planXf,c11np(i,j,:),c11ns(i,j,:))
            call fftw_execute_dft_r2c(planXf,c12np(i,j,:),c12ns(i,j,:))
            call fftw_execute_dft_r2c(planXf,c13np(i,j,:),c13ns(i,j,:))
            call fftw_execute_dft_r2c(planXf,c22np(i,j,:),c22ns(i,j,:))
            call fftw_execute_dft_r2c(planXf,c23np(i,j,:),c23ns(i,j,:))
            call fftw_execute_dft_r2c(planXf,c33np(i,j,:),c33ns(i,j,:))

            call fftw_execute_dft_r2c(planXf,str11np(i,j,:),str11ns(i,j,:))
            call fftw_execute_dft_r2c(planXf,str12np(i,j,:),str12ns(i,j,:))
            call fftw_execute_dft_r2c(planXf,str13np(i,j,:),str13ns(i,j,:))
            call fftw_execute_dft_r2c(planXf,str22np(i,j,:),str22ns(i,j,:))
            call fftw_execute_dft_r2c(planXf,str23np(i,j,:),str23ns(i,j,:))
            call fftw_execute_dft_r2c(planXf,str33np(i,j,:),str33ns(i,j,:))

            call fftw_execute_dft_r2c(planXf,qp11np(i,j,:),qp11s(i,j,:))
            call fftw_execute_dft_r2c(planXf,qp12np(i,j,:),qp12s(i,j,:))
            call fftw_execute_dft_r2c(planXf,qp13np(i,j,:),qp13s(i,j,:))
            call fftw_execute_dft_r2c(planXf,qp22np(i,j,:),qp22s(i,j,:))
            call fftw_execute_dft_r2c(planXf,qp23np(i,j,:),qp23s(i,j,:))
            call fftw_execute_dft_r2c(planXf,qp33np(i,j,:),qp33s(i,j,:))
            end if
#ENDIF
        end do
    end do
    !$omp end parallel do

    ! Normalize values by mx
    vwxs = vwxs/float(mx)
    vwys = vwys/float(mx)
    vwzs = vwzs/float(mx)
#IFDEF SCALAR
    vcs  =  vcs/float(mx)
#ENDIF
#IFDEF POLYMER  
    if (it .ge. (src_start-1)) then
    c11ns = c11ns/float(mx)
    c12ns = c12ns/float(mx)
    c13ns = c13ns/float(mx)
    c21ns = c21ns/float(mx)
    c22ns = c22ns/float(mx)
    c33ns = c33ns/float(mx)

    str11ns = str11ns/float(mx)
    str12ns = str12ns/float(mx)
    str13ns = str13ns/float(mx)
    str21ns = str21ns/float(mx)
    str22ns = str22ns/float(mx)
    str33ns = str33ns/float(mx)

    qp11s = qp11s/float(mx)
    qp12s = qp12s/float(mx)
    qp13s = qp13s/float(mx)
    qp21s = qp21s/float(mx)
    qp22s = qp22s/float(mx)
    qp33s = qp33s/float(mx)
    end if
#ENDIF


    ! Complex --> Complex Transform (z-direction)
    !$omp parallel do default(shared) private(i,k)
    do k = 1,nxh
        do i = 1,nyp
            call fftw_execute_dft(planZf,vwxs(i,:,k),vwxs(i,:,k))
            call fftw_execute_dft(planZf,vwys(i,:,k),vwys(i,:,k))
            call fftw_execute_dft(planZf,vwzs(i,:,k),vwzs(i,:,k))
#IFDEF SCALAR
            call fftw_execute_dft(planZf,vcs(i,:,k),vcs(i,:,k))
#ENDIF
#IFDEF POLYMER
            if (it .ge. (src_start-1)) then
            call fftw_execute_dft(planZf,c11ns(i,:,k),c11ns(i,:,k))
            call fftw_execute_dft(planZf,c12ns(i,:,k),c12ns(i,:,k))
            call fftw_execute_dft(planZf,c13ns(i,:,k),c13ns(i,:,k))
            call fftw_execute_dft(planZf,c22ns(i,:,k),c22ns(i,:,k))
            call fftw_execute_dft(planZf,c23ns(i,:,k),c23ns(i,:,k))
            call fftw_execute_dft(planZf,c33ns(i,:,k),c33ns(i,:,k))

            call fftw_execute_dft(planZf,str11ns(i,:,k),str11ns(i,:,k))
            call fftw_execute_dft(planZf,str12ns(i,:,k),str12ns(i,:,k))
            call fftw_execute_dft(planZf,str13ns(i,:,k),str13ns(i,:,k))
            call fftw_execute_dft(planZf,str22ns(i,:,k),str22ns(i,:,k))
            call fftw_execute_dft(planZf,str23ns(i,:,k),str23ns(i,:,k))
            call fftw_execute_dft(planZf,str33ns(i,:,k),str33ns(i,:,k))

            call fftw_execute_dft(planZf,qp11s(i,:,k),qp11s(i,:,k))
            call fftw_execute_dft(planZf,qp12s(i,:,k),qp12s(i,:,k))
            call fftw_execute_dft(planZf,qp13s(i,:,k),qp13s(i,:,k))
            call fftw_execute_dft(planZf,qp22s(i,:,k),qp22s(i,:,k))
            call fftw_execute_dft(planZf,qp23s(i,:,k),qp23s(i,:,k))
            call fftw_execute_dft(planZf,qp33s(i,:,k),qp33s(i,:,k))
            end if
#ENDIF
        end do
    end do
    !$omp end parallel do

    ! Normalize values by mz
    vwxs = vwxs/float(mz)
    vwys = vwys/float(mz)
    vwzs = vwzs/float(mz)
#IFDEF SCALAR
    vcs  =  vcs/float(mz)
#ENDIF
#IFDEF POLYMER  
    if (it .ge. (src_start-1)) then
    c11ns = c11ns/float(mz)
    c12ns = c12ns/float(mz)
    c13ns = c13ns/float(mz)
    c21ns = c21ns/float(mz)
    c22ns = c22ns/float(mz)
    c33ns = c33ns/float(mz)

    str11ns = str11ns/float(mz)
    str12ns = str12ns/float(mz)
    str13ns = str13ns/float(mz)
    str21ns = str21ns/float(mz)
    str22ns = str22ns/float(mz)
    str33ns = str33ns/float(mz)

    qp11s = qp11s/float(mz)
    qp12s = qp12s/float(mz)
    qp13s = qp13s/float(mz)
    qp21s = qp21s/float(mz)
    qp22s = qp22s/float(mz)
    qp33s = qp33s/float(mz)
    end if
#ENDIF

    ! Fill in regular spectral variables | We cut out highest 1/3 of modes
    ! for de-aliasing
    !$omp parallel do default(shared) private(i,j,k,jj)
    do k = 1,nxh
        do j = 1,nz
            if (j .le. nzh) jj = j
            if (j .gt. nzh) jj = (mz-nz) + j
            do i = 1,nyp
                gn(i,j,k) = vwxs(i,jj,k)
                fn(i,j,k) = vwys(i,jj,k)
               omz(i,j,k) = vwzs(i,jj,k)
#IFDEF SCALAR
               scn(i,j,k) = vcs(i,jj,k)
#ENDIF
#IFDEF POLYMER
                if (it .ge. src_start-1) then
                c11n(i,j,k) = c11ns(i,jj,k)
                c12n(i,j,k) = c12ns(i,jj,k)
                c13n(i,j,k) = c13ns(i,jj,k)
                c22n(i,j,k) = c22ns(i,jj,k)
                c23n(i,j,k) = c23ns(i,jj,k)
                c33n(i,j,k) = c33ns(i,jj,k)

                str11n(i,j,k) = str11ns(i,jj,k)
                str12n(i,j,k) = str12ns(i,jj,k)
                str13n(i,j,k) = str13ns(i,jj,k)
                str22n(i,j,k) = str22ns(i,jj,k)
                str23n(i,j,k) = str23ns(i,jj,k)
                str33n(i,j,k) = str33ns(i,jj,k)

                qp11(i,j,k) = qp11s(i,jj,k)
                qp12(i,j,k) = qp12s(i,jj,k)
                qp13(i,j,k) = qp13s(i,jj,k)
                qp22(i,j,k) = qp22s(i,jj,k)
                qp23(i,j,k) = qp23s(i,jj,k)
                qp33(i,j,k) = qp33s(i,jj,k)
                end if
#ENDIF
            end do
        end do
    end do
    !$omp end parallel do


    !---------------------------------------------------------------------!
    !     Calculate swirl criterion and write data for visualization      !
    !---------------------------------------------------------------------!

    !$omp parallel do default(shared) private(i,j,k,swirl)
    do k = 1,mx
        do j = 1,mz
            do i = 1,nyp
!                call calcswirl(u11p(i,j,k),u21p(i,j,k),u31p(i,j,k),u12p(i,j,k),u22p(i,j,k), &
!                               u32p(i,j,k),u13p(i,j,k),u23p(i,j,k),u33p(i,j,k),swirl)
    
                swirl_3d(i,j,k) = swirl
            end do
        end do
    end do
    !$omp end parallel do

    ! Process 3D variables (write outputs in physical space)
    if (print3d .ne. 0) then
    
        if ((mod(it,iprnfrq) .eq. 0 .and. it .ne. 0) .or. it .eq. 1) then
            
            ! Write output files
            if (print3d .eq. 1) then ! Write output in ASCII format
#IFDEF POLYMER
                call write_flowfield_ascii(up,vp,wp,wxp,wyp,wzp,swirl_3d,beta_poly,float(imatrix))
#ELIF DEFINED SCALAR
                call write_flowfield_ascii(up,vp,wp,wxp,wyp,wzp,swirl_3d,scp,float(imatrix))
#ELSE
!                call write_flowfield_ascii(up,vp,wp,wxp,wyp,wzp,swirl_3d,float(imatrix))
#ENDIF
#IFDEF OUTPUTFORM
            else if (print3d .eq. 3) then ! Write output in Tecplot binary (.szplt)
#IFDEF POLYMER
                call write_flowfield_plt(up,vp,wp,wxp,wyp,wzp,swirl_3d,beta_poly)
#ELIF DEFINED SCALAR                                         
                call write_flowfield_plt(up,vp,wp,wxp,wyp,wzp,swirl_3d,scp)
#ELSE                                                        
                call write_flowfield_plt(up,vp,wp,wxp,wyp,wzp,swirl_3d)
#ENDIF
#ENDIF
            else if (print3d .eq. 2) then ! Write outputs specifically for FTLE
!                call write_FTLE_output(up,vp,wp)
            else
                write(*,*) 'Warning: Unknown print type. No output data will be written.'
            end if
        
            write(*,*) '    Done!'
        end if
    end if
   
    !---------------------------------------------------------------------!
    !                   Particle Tracking Implementation                  !
    !---------------------------------------------------------------------!
    
    if (npart .ne. 0) then
        call part_track(up,vp,wp,wxp,wyp,wzp,u_old,v_old,w_old,  &
                        u11p,u12p,u13p,u21p,u22p,u23p,u31p, &
                        u32p,u33p,Lup,Lvp,Lwp,Lu_old,Lv_old,Lw_old)!,vArea)
    
        ! Save velocity and Laplacian for time derivatives inside particle integration
        u_old = up
        v_old = vp
        w_old = wp
    
        Lu_old = Lup
        Lv_old = Lvp
        Lw_old = Lwp
    end if
   
    ! Calculate Q-criterion at particle locations
    do n = 1,npart
        call fluid_interp1(xpart(n),ypart(n),zpart(n),swirl_3d,swirl_part(n))
    end do
 
    !---------------------------------------------------------------------!
    ! Write Mean U Data, calculate mass flux, & KE/enstrophy if relevant  !
    !---------------------------------------------------------------------!
  
    ! Make these into subroutines later...  
!    if (flow_select .eq. 1 .or. flow_select .eq. 4) then ! Only relevant for wall-bounded turbulence
!        write(*,*) 'Writing mean U data...'
!    
!        ! Mean velocity
!!        if (irstrt .eq. it) then
!!            open(71,file = 'outputs/mean_u_data.dat')
!!        else
!!            open(71,file = 'outputs/mean_u_data.dat', position = 'append')
!!        end if
!    
!        do i = 1,nyp
!            do j = 1,mz
!                uxmean(j) = sum(up(i,j,:))/mx
!            end do
!            uzmean(i) = sum(uxmean)/mz
!!            write(71,"(*(e14.6,1x))") uzmean
!        end do
!    
!!        close(71)
!        write(*,*) '    Done!'
!    end if
    
    ! Compute mass flux: Integrate <U> over y
    if (irstrt .eq. it) then
        open(72,file='outputs/mass_flux')
    else
        open(72,file='outputs/mass_flux',position='append')
    end if
    
    massFlux = 0.0
    do i = 1,ny
        massFlux = massFlux + 1.0/yl*0.5*(uzmean(i+1) + uzmean(i))*abs(ycoord(i+1) - ycoord(i)) ! Bulk velocity
    end do
    
    write(72,*) massFlux
    close(72)

#IFDEF SCALAR
    ! Calculate global enstrophy, TKE, and S-gamma correlation
        call writeoutputs(up,vp,wp,wxp,wyp,wzp, &
#IFDEF POLYMER
                          beta_poly, &
#ELSE
                          scp, &
#ENDIF
                          u11p,u22p,u33p,uzmean)
#ENDIF
#IFDEF POLYMER
    if (scl_flag .eq. 2 .and. it .le. src_stop .and. it .ge. src_start) then
        call calc_total_beta(it,delxm,delzm,scp,beta_poly)
    end if

    ! Calculate correlation between beta and vortex structure (swirl/Q)
    call correlate_vars(beta_poly,swirl_3d,it)
#ENDIF
 
    !---------------------------------------------------------------------!
    !               Check for numerical instabilities                     !
    !---------------------------------------------------------------------!
    
    do i = 1,nyp
        if (cfl(i) > cflmax) cflmax = cfl(i)
    end do
    
    write(*,*) 'max CFL = ',cflmax
    
    ! Check CFL condition
    if (cflmax .gt. 1.0) then
        write(*,*) 'CFL failure at time step ',it
        stop
    end if
    
    end subroutine vcw3dp
end module solvers


