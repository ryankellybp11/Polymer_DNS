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
    common/solver/     gain,ugain,theta,alpha,beta,dyde
    common/waves/      wavx,wavz,c
    common/flow/       re,Uinf,R_tau,dPdx
    common/itime/      it,dt
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
    call odds(u,bctop,bcbot,ib,c,wrkc)
    
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
    
    subroutine odds(s,bctop,bcbot,ib,c,wrk)
    ! Sets up RHS for even tri-diagonal problem, given the results of sbr RHS
    use grid_size
    implicit none
    
    complex, dimension(nyp,nz) :: wrk
    complex, dimension(nz)     :: bctop, bcbot
    complex :: s(0:ny,nz) ! note zero-based arrays
    real    :: c(nyp)
    
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
    integer :: i,j,k,m,mm1,mm2
    
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
    common/solver/     gain,ugain,theta,alpha,beta,dyde
    common/waves/      wavx,wavz,c
    common/flow/       re,Uinf,R_tau,dPdx
    common/itime/      it,dt
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
    integer :: i,j,k,nrank,kzero,ip
    
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
    common/itime/      it,dt
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
    call pntodds(u,c,wrkc)
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
    
    subroutine pntodds(s,c,wrk)
    ! Sets up RHS for even tri-diagonal problem, given the results of sbr RHS
    use grid_size
    implicit none
    
    complex, dimension(nyp,nz) :: wrk
    complex :: s(0:ny,nz) ! note zero-based arrays
    real    :: c(0:ny)
    
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
    
    integer :: i,j,k,ip,nrank,ipass,js,n,m,mm1,mm2,ii   
    
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
    
    integer :: i,j,k
    
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
    integer :: i,j,k
    
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
    common/itime/  it,dt
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
    real,    dimension(nyp)     :: t
    real :: x,g,wavz,dyde,rj1,rj2,rj3,wn
    integer :: i,j,k
    
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
    real, dimension(nyp)        :: t
    real :: x,g,wavz,dyde,rj1,rj2,rj3,wn
    integer :: i,j,k
    
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
    
    subroutine vort(u,v,w,omx,omz,wrkc)
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
    complex, dimension(nyp,nz,nxh) :: u,v,w,omx,omz,wrkc
    
    ! Calculation variables
    complex :: im
    integer :: i,j,k
    
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
    
    subroutine vcw3d(u,v,w,omx,omy,omz,fn,gn,u11,u12,u13,u21,u22,u23,u31, &
                     u32,u33,Lu,Lv,Lw,scalar,sclx,scly,sclz,scn,          &
                     Lu_old,Lv_old,Lw_old)
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
    !  real space in y and Fourier in x and z. Before this subroutine is  !
    !  called, we must to a y-transformation on all arrays containing     !
    !  velocity, vorticity, or velocity derivatives. Also, we assume that !
    !  the appropriate modes in all fields (i.e., the first and last      !
    !  modes) are set to zero prior to this point.                        !
    !---------------------------------------------------------------------!
    ! =================================================================== !
    !                           Declare Modules                           !
    ! =================================================================== !
    use omp_lib
    use grid_size
    use derivs
    use helpers
    !---------------------------------------------------------------------!
    
    ! =================================================================== !
    !                   Declare all variables explicitly                  !
    ! =================================================================== !
    implicit none
    
    ! Passed variables
    complex, dimension(nyp,nz,nxh) :: u,v,w,omx,omy,omz,fn,gn
    complex, dimension(nyp,nz,nxh) :: u11,u12,u13,u21,u22,u23,u31,u32,u33
    complex, dimension(nyp,nz,nxh) :: Lu,Lv,Lw
    complex, dimension(nyp,nz,nxh) :: scalar,sclx,scly,sclz,scn
    
    ! Calculation variables
    real, save, dimension(nyp,mz,mx) :: u_old,v_old,w_old
    real, dimension(nyp,mz,mx) :: Lu_old,Lv_old,Lw_old
    real, dimension(nyp,mz,mx) :: up3d,vp3d,wp3d,wx3d,wy3d,wz3d
    real, dimension(nyp,mz,mx) :: u11p3d,u12p3d,u13p3d
    real, dimension(nyp,mz,mx) :: u21p3d,u22p3d,u23p3d
    real, dimension(nyp,mz,mx) :: u31p3d,u32p3d,u33p3d
    real, dimension(nyp,mz,mx) :: Lup3d,Lvp3d,Lwp3d,swirl_3d
    real, dimension(nyp,mz,mx) :: scp3d,scsource
    
    real, dimension(mzp,mxp2)  :: up,vp,wp,wx,wy,wz,vwx,vwy,vwz
    real, dimension(mzp,mxp2)  :: u11p,u12p,u13p
    real, dimension(mzp,mxp2)  :: u21p,u22p,u23p
    real, dimension(mzp,mxp2)  :: u31p,u32p,u33p
    real, dimension(mzp,mxp2)  :: Lup,Lvp,Lwp
    real, dimension(mzp,mxp2)  :: scp,cx,cy,cz,vc
    
    integer :: i,j,k,i1,i2,inc,isgn,jump,lot,n
    integer :: ii,jj,ipii,jpjj,idif,jdif
    real    :: delx,delz,delxm,delzm,pi
    real    :: segdrag,xsegdrag,ysegdrag,zsegdrag,wdes
    real    :: cflcheck,cflmax,swirl
    
    real,dimension(nwrk)  :: wrk
    real,dimension(nyp)   :: cfl
    real,dimension(mz,mx) :: dragx,dragy,dragz
    
    real :: xi,zj,argx,argrad,fx,fr,rsq
    real :: uxmean(mz),uzmean(nyp)
    real :: Lx,Ly,Lz,massFlux,sumens,KE,scl_total
    
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
    real    :: vArea
    
    ! Immersed boundary force variables
    real, dimension(nyp,mz,mx) :: fxintg,fyintg,fzintg
    real, dimension(mz2,mx2)   :: fspread

    ! Scalar variables
    real    :: sigmax,sigmay,sigmaz
    real    :: betax,betay,betaz
    real    :: xsq,ysq,zsq
    real    :: xc1,yc1,zc1
    real    :: deltaT,diff
    integer :: scl_flag
    integer :: src_start,src_stop

    ! Particle variables
    real,save,dimension(npart) :: xpart,ypart,zpart

    !---------------------------------------------------------------------!
    
    ! =================================================================== !
    !                            Common Blocks                            !
    ! =================================================================== !
    common/iocontrl/   irstrt,nsteps,iprnfrq,print3d
    common/solver/     gain,ugain,theta,alpha,beta,dyde
    common/ibforce/    fxintg,fyintg,fzintg,fspread
    common/waves/      wavx,wavz,c
    common/init/       initu,initv,initw
    common/buffer/     bfgain,bfugain,vdes,bfhead,bftail,bfwidth,slopelength
    common/flow/       re,Uinf,R_tau,dPdx
    common/domain/     xl,yl,zl
    common/trig/       trigx,trigy,trigz,sine,cosine,trigx32,trigz32
    common/ifax/       ixfax,iyfax,izfax,ixfax32,izfax32
    common/itime/      it,dt
    common/imat/       imatrix,kwall,kmaxsurf
    common/vortexring/ forbeta,xcenter,ycenter,zcenter,L,rad,bdyfx
    common/setup/      geomtype,flow_select,perturbtime
    common/scl_stuff/  sigmax,sigmay,sigmaz,deltaT,diff,scl_flag
    common/part_traj/  xpart,ypart,zpart
    common/src_time/   src_start,src_stop
    !---------------------------------------------------------------------!
    
    ! =================================================================== !
    !                         Begin Calculations                          !
    ! =================================================================== !
    
    !---------------------------------------------------------------------!
    !                    Initialize some variables                        !
    !---------------------------------------------------------------------!
   
    pi = 2.0*acos(0.0)
    delx = xl/float(nx-1)
    delz = zl/float(nz-1)
    delxm = xl/float(mx-1)
    delzm = zl/float(mz-1)
    
    cfl = 0.0
    cflmax = 0.0
    
    !---------------------------------------------------------------------!
    !                    Calculate scalar source term                     !
    !---------------------------------------------------------------------!
    if (scl_flag .eq. 2 .and. it .ge. src_start .and. it .le. src_stop) then
        scsource = 0.0
        do n = 1,npart
        if (it .eq. 0) then
            open(95,file='setup/particles/particles.dat',status='old',action='read')
            do j = 1,n
                read(95,*)
            end do
            read(95,*) xc1,yc1,zc1
            close(95)
        else
            xc1 = xpart(n)
            yc1 = ypart(n)
            zc1 = zpart(n)
        end if
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
        end do
    else if (scl_flag .eq. 3 .and. it .eq. 1 .and. irstrt .eq. 1) then
        scsource = 1.0
    else
        scsource = 0.0
    end if

    !---------------------------------------------------------------------!
    !         Calculate (u x omega), iterating over y-planes              !
    !---------------------------------------------------------------------!
    
    !$omp parallel do default(shared) private(i,j,k,ii,i1,i2,ipii,idif,jj, &
    !$omp    jpjj,jdif,segdrag,xsegdrag,ysegdrag,zsegdrag,cflcheck,        &
    !$omp    up,vp,wp,wx,wy,wz,vwx,vwy,vwz,dragx,dragy,dragz,wrk,inc,isgn, &
    !$omp    u11p,u12p,u13p,u21p,u22p,u23p,u31p,u32p,u33p,Lup,Lvp,Lwp,     &
    !$omp    jump,lot,fx,fr,xi,zj,argrad,rsq,argx,scp,cx,cy,cz,vc) schedule(dynamic)
    
    do k = 1,nyp
        ! Initialize the velocities and vorticities on the kth x-z plane 
        ! with zeros. We need to pad the extra 1/2 modes with zeros by 3/2 rule
    
        do j = 1,mzp
            do i = 1,mxp2
    
                up(j,i) = 0.0
                vp(j,i) = 0.0
                wp(j,i) = 0.0
    
                wx(j,i) = 0.0
                wy(j,i) = 0.0
                wz(j,i) = 0.0
    
                vwx(j,i) = 0.0
                vwy(j,i) = 0.0
                vwz(j,i) = 0.0
    
                ! Used in particle tracking
                u11p(j,i) = 0.0
                u12p(j,i) = 0.0
                u13p(j,i) = 0.0
                u21p(j,i) = 0.0
                u22p(j,i) = 0.0
                u23p(j,i) = 0.0
                u31p(j,i) = 0.0
                u32p(j,i) = 0.0
                u33p(j,i) = 0.0
          
                Lup(j,i) = 0.0
                Lvp(j,i) = 0.0
                Lwp(j,i) = 0.0  
   
                scp(j,i) = 0.0
                 cx(j,i) = 0.0 
                 cy(j,i) = 0.0 
                 cz(j,i) = 0.0 
                 vc(j,i) = 0.0 
            end do
        end do
    
        do j = 1,nz
            if (j .le. nzh) jj = j
            if (j .gt. nzh) jj = (mz-nz) + j
            
            do i = 1,nxh
                i1 = 2*(i-1) + 1
                i2 = 2*i
    
                up(jj,i1) =  real(u(k,j,i))
                up(jj,i2) = aimag(u(k,j,i))
                vp(jj,i1) =  real(v(k,j,i))
                vp(jj,i2) = aimag(v(k,j,i))
                wp(jj,i1) =  real(w(k,j,i))
                wp(jj,i2) = aimag(w(k,j,i))
    
                wx(jj,i1) =  real(omx(k,j,i))
                wx(jj,i2) = aimag(omx(k,j,i))
                wy(jj,i1) =  real(omy(k,j,i))
                wy(jj,i2) = aimag(omy(k,j,i))
                wz(jj,i1) =  real(omz(k,j,i))
                wz(jj,i2) = aimag(omz(k,j,i))
    
                u11p(jj,i1) =  real(u11(k,j,i))
                u11p(jj,i2) = aimag(u11(k,j,i))
                u12p(jj,i1) =  real(u12(k,j,i))
                u12p(jj,i2) = aimag(u12(k,j,i))
                u13p(jj,i1) =  real(u13(k,j,i))
                u13p(jj,i2) = aimag(u13(k,j,i))
                u21p(jj,i1) =  real(u21(k,j,i))
                u21p(jj,i2) = aimag(u21(k,j,i))
                u22p(jj,i1) =  real(u22(k,j,i))
                u22p(jj,i2) = aimag(u22(k,j,i))
                u23p(jj,i1) =  real(u23(k,j,i))
                u23p(jj,i2) = aimag(u23(k,j,i))
                u31p(jj,i1) =  real(u31(k,j,i))
                u31p(jj,i2) = aimag(u31(k,j,i))
                u32p(jj,i1) =  real(u32(k,j,i))
                u32p(jj,i2) = aimag(u32(k,j,i))
                u33p(jj,i1) =  real(u33(k,j,i))
                u33p(jj,i2) = aimag(u33(k,j,i))
                
                Lup(jj,i1) =  real(Lu(k,j,i))
                Lup(jj,i2) = aimag(Lu(k,j,i))
                Lvp(jj,i1) =  real(Lv(k,j,i))
                Lvp(jj,i2) = aimag(Lv(k,j,i))
                Lwp(jj,i1) =  real(Lw(k,j,i))
                Lwp(jj,i2) = aimag(Lw(k,j,i))
    
                scp(jj,i1) =  real(scalar(k,j,i))
                scp(jj,i2) = aimag(scalar(k,j,i))
                cx(jj,i1)  =  real(sclx(k,j,i))
                cx(jj,i2)  = aimag(sclx(k,j,i))
                cy(jj,i1)  =  real(scly(k,j,i))
                cy(jj,i2)  = aimag(scly(k,j,i))
                cz(jj,i1)  =  real(sclz(k,j,i))
                cz(jj,i2)  = aimag(sclz(k,j,i))
            end do
        end do
   
    !---------------------------------------------------------------------!
    !        Transform (interpolate) into 3/2 grid physical space         !
    !---------------------------------------------------------------------!
    
        inc  = 1
        isgn = 1
        jump = 2*mzp
        lot  = nx/2
    
        call cfftmlt(up(1,1),up(1,2),wrk,trigz32,izfax32,inc,jump,mz,lot,isgn)
        call cfftmlt(vp(1,1),vp(1,2),wrk,trigz32,izfax32,inc,jump,mz,lot,isgn)
        call cfftmlt(wp(1,1),wp(1,2),wrk,trigz32,izfax32,inc,jump,mz,lot,isgn)
        call cfftmlt(wx(1,1),wx(1,2),wrk,trigz32,izfax32,inc,jump,mz,lot,isgn)
        call cfftmlt(wy(1,1),wy(1,2),wrk,trigz32,izfax32,inc,jump,mz,lot,isgn)
        call cfftmlt(wz(1,1),wz(1,2),wrk,trigz32,izfax32,inc,jump,mz,lot,isgn)
    
        call cfftmlt(u11p(1,1),u11p(1,2),wrk,trigz32,izfax32,inc,jump,mz,lot,isgn)
        call cfftmlt(u12p(1,1),u12p(1,2),wrk,trigz32,izfax32,inc,jump,mz,lot,isgn)
        call cfftmlt(u13p(1,1),u13p(1,2),wrk,trigz32,izfax32,inc,jump,mz,lot,isgn)
        call cfftmlt(u21p(1,1),u21p(1,2),wrk,trigz32,izfax32,inc,jump,mz,lot,isgn)
        call cfftmlt(u22p(1,1),u22p(1,2),wrk,trigz32,izfax32,inc,jump,mz,lot,isgn)
        call cfftmlt(u23p(1,1),u23p(1,2),wrk,trigz32,izfax32,inc,jump,mz,lot,isgn)
        call cfftmlt(u31p(1,1),u31p(1,2),wrk,trigz32,izfax32,inc,jump,mz,lot,isgn)
        call cfftmlt(u32p(1,1),u32p(1,2),wrk,trigz32,izfax32,inc,jump,mz,lot,isgn)
        call cfftmlt(u33p(1,1),u33p(1,2),wrk,trigz32,izfax32,inc,jump,mz,lot,isgn)
     
        call cfftmlt(Lup(1,1),Lup(1,2),wrk,trigz32,izfax32,inc,jump,mz,lot,isgn)
        call cfftmlt(Lvp(1,1),Lvp(1,2),wrk,trigz32,izfax32,inc,jump,mz,lot,isgn)
        call cfftmlt(Lwp(1,1),Lwp(1,2),wrk,trigz32,izfax32,inc,jump,mz,lot,isgn)
    
        call cfftmlt(scp(1,1),scp(1,2),wrk,trigz32,izfax32,inc,jump,mz,lot,isgn)
        call cfftmlt(cx(1,1),cx(1,2),wrk,trigz32,izfax32,inc,jump,mz,lot,isgn)
        call cfftmlt(cy(1,1),cy(1,2),wrk,trigz32,izfax32,inc,jump,mz,lot,isgn)
        call cfftmlt(cz(1,1),cz(1,2),wrk,trigz32,izfax32,inc,jump,mz,lot,isgn)

        do j = 1,mz
            up(j,nxp) = up(j,2)
            vp(j,nxp) = vp(j,2)
            wp(j,nxp) = wp(j,2)
            wx(j,nxp) = wx(j,2)
            wy(j,nxp) = wy(j,2)
            wz(j,nxp) = wz(j,2)
            
            u11p(j,nxp) = u11p(j,2)
            u12p(j,nxp) = u12p(j,2)
            u13p(j,nxp) = u13p(j,2)
            u21p(j,nxp) = u21p(j,2)
            u22p(j,nxp) = u22p(j,2)
            u23p(j,nxp) = u23p(j,2)
            u31p(j,nxp) = u31p(j,2)
            u32p(j,nxp) = u32p(j,2)
            u33p(j,nxp) = u33p(j,2)
            
            Lup(j,nxp) = Lup(j,2) 
            Lvp(j,nxp) = Lvp(j,2) 
            Lwp(j,nxp) = Lwp(j,2) 
    
            scp(j,nxp) = scp(j,2)
             cx(j,nxp) =  cx(j,2)
             cy(j,nxp) =  cy(j,2)
             cz(j,nxp) =  cz(j,2)
            
            up(j,2) = 0.0
            vp(j,2) = 0.0
            wp(j,2) = 0.0
            wx(j,2) = 0.0
            wy(j,2) = 0.0
            wz(j,2) = 0.0
            
            u11p(j,2) = 0.0
            u12p(j,2) = 0.0
            u13p(j,2) = 0.0
            u21p(j,2) = 0.0
            u22p(j,2) = 0.0
            u23p(j,2) = 0.0
            u31p(j,2) = 0.0
            u32p(j,2) = 0.0
            u33p(j,2) = 0.0
            
            Lup(j,2) = 0.0
            Lvp(j,2) = 0.0
            Lwp(j,2) = 0.0

            scp(j,2) = 0.0
             cx(j,2) = 0.0
             cy(j,2) = 0.0
             cz(j,2) = 0.0
            
            up(j,nxp2) = up(j,2)
            vp(j,nxp2) = vp(j,2)
            wp(j,nxp2) = wp(j,2)
            wx(j,nxp2) = wx(j,2)
            wy(j,nxp2) = wy(j,2)
            wz(j,nxp2) = wz(j,2)
            
            u11p(j,nxp2) = u11p(j,2)
            u12p(j,nxp2) = u12p(j,2)
            u13p(j,nxp2) = u13p(j,2)
            u21p(j,nxp2) = u21p(j,2)
            u22p(j,nxp2) = u22p(j,2)
            u23p(j,nxp2) = u23p(j,2)
            u31p(j,nxp2) = u31p(j,2)
            u32p(j,nxp2) = u32p(j,2)
            u33p(j,nxp2) = u33p(j,2)
            
            Lup(j,nxp2) = Lup(j,2) 
            Lvp(j,nxp2) = Lvp(j,2) 
            Lwp(j,nxp2) = Lwp(j,2) 

            scp(j,nxp) = scp(j,2)
             cx(j,nxp) =  cx(j,2)
             cy(j,nxp) =  cy(j,2)
             cz(j,nxp) =  cz(j,2)
        end do
            
    
        isgn = 1
        inc  = mzp
        jump = 1
        lot  = mz
    
        call rfftmlt(up,wrk,trigx32,ixfax32,inc,jump,mx,lot,isgn)
        call rfftmlt(vp,wrk,trigx32,ixfax32,inc,jump,mx,lot,isgn)
        call rfftmlt(wp,wrk,trigx32,ixfax32,inc,jump,mx,lot,isgn)
        call rfftmlt(wx,wrk,trigx32,ixfax32,inc,jump,mx,lot,isgn)
        call rfftmlt(wy,wrk,trigx32,ixfax32,inc,jump,mx,lot,isgn)
        call rfftmlt(wz,wrk,trigx32,ixfax32,inc,jump,mx,lot,isgn)
    
        call rfftmlt(u11p,wrk,trigx32,ixfax32,inc,jump,mx,lot,isgn)
        call rfftmlt(u12p,wrk,trigx32,ixfax32,inc,jump,mx,lot,isgn)
        call rfftmlt(u13p,wrk,trigx32,ixfax32,inc,jump,mx,lot,isgn)
        call rfftmlt(u21p,wrk,trigx32,ixfax32,inc,jump,mx,lot,isgn)
        call rfftmlt(u22p,wrk,trigx32,ixfax32,inc,jump,mx,lot,isgn)
        call rfftmlt(u23p,wrk,trigx32,ixfax32,inc,jump,mx,lot,isgn)
        call rfftmlt(u31p,wrk,trigx32,ixfax32,inc,jump,mx,lot,isgn)
        call rfftmlt(u32p,wrk,trigx32,ixfax32,inc,jump,mx,lot,isgn)
        call rfftmlt(u33p,wrk,trigx32,ixfax32,inc,jump,mx,lot,isgn)
     
        call rfftmlt(Lup,wrk,trigx32,ixfax32,inc,jump,mx,lot,isgn)
        call rfftmlt(Lvp,wrk,trigx32,ixfax32,inc,jump,mx,lot,isgn)
        call rfftmlt(Lwp,wrk,trigx32,ixfax32,inc,jump,mx,lot,isgn)
    
        call rfftmlt(scp,wrk,trigx32,ixfax32,inc,jump,mx,lot,isgn)
        call rfftmlt(cx ,wrk,trigx32,ixfax32,inc,jump,mx,lot,isgn)
        call rfftmlt(cy ,wrk,trigx32,ixfax32,inc,jump,mx,lot,isgn)
        call rfftmlt(cz ,wrk,trigx32,ixfax32,inc,jump,mx,lot,isgn)
    !---------------------------------------------------------------------!
    !        Compute the cross product of velocity and vorticity          !
    !---------------------------------------------------------------------!
    
        do j = 1,mz
            do i = 1,mx
                vwx(j,i) =  vp(j,i)*wz(j,i) - wp(j,i)*wy(j,i)
                vwy(j,i) = -up(j,i)*wz(j,i) + wp(j,i)*wx(j,i)
                vwz(j,i) =  up(j,i)*wy(j,i) - vp(j,i)*wx(j,i)

                 vc(j,i) = -(up(j,i)*cx(j,i) + vp(j,i)*cy(j,i) + wp(j,i)*cz(j,i))
    
                cflcheck = (abs(up(j,i))/delxm + abs(vp(j,i))/seght(k) + abs(wp(j,i))/delzm)*dt
                if(cflcheck .gt. cfl(k)) then
                     cfl(k) = cflcheck
                end if
            end do
        end do
    
    !---------------------------------------------------------------------!
    ! Compute immersed boundary force terms and add to nonlinear term     !
    ! The type of forcing is based on the value of imatrix at the grid    !
    ! point. Several of these cases probably require more coding and      !
    ! specification in the geometry file, but they have been unused since !
    ! I've used the code                                                  !
    !---------------------------------------------------------------------!
    
        ! Initialize the force field array
        dragx = 0.0
        dragy = 0.0
        dragz = 0.0
   
        do j = 1,mz
            do i = 1,mx
                if (imatrix(k,j,i) .eq. 1 .or. (imatrix(k,j,i) .eq. 4 .and. it .le. perturbtime)) then ! Generate solid surfaces
                    fxintg(k,j,i) = fxintg(k,j,i) + up(j,i)
                    fyintg(k,j,i) = fyintg(k,j,i) + vp(j,i)
                    fzintg(k,j,i) = fzintg(k,j,i) + wp(j,i)
                    dragx(j,i) = (-ugain*up(j,i)) - gain*fxintg(k,j,i)
                    dragy(j,i) = (-ugain*vp(j,i)) - gain*fyintg(k,j,i)
                    dragz(j,i) = (-ugain*wp(j,i)) - gain*fzintg(k,j,i)
        
                else if (imatrix(k,j,i) .eq. 3) then ! Spanwise-damping textures
                    fzintg(k,j,i) = fzintg(k,j,i) + wp(j,i)
                    dragz(j,i) = -ugain*wp(j,i) - gain*fzintg(k,j,i)
        
                else if (imatrix(k,j,i) .eq. 8) then ! Spanwise moving wall
                    if (i .le. (bfwidth*3 + 1)) then
                        wdes = (0.3*Uinf)/slopelength*(i - (bfwidth*3 + 1 - slopelength))
                    else
                        wdes = 0.3*Uinf
                    end if
    
                    fxintg(k,j,i) = fxintg(k,j,i) + up(j,i)
                    fyintg(k,j,i) = fyintg(k,j,i) + vp(j,i)
                    fzintg(k,j,i) = fzintg(k,j,i) + wp(j,i)-wdes
                    dragx(j,i)    = -ugain*up(j,i) - gain*fxintg(k,j,i)
                    dragy(j,i)    = -ugain*vp(j,i) - gain*fyintg(k,j,i)
                    dragz(j,i)    = -ugain*(wp(j,i)-wdes) - gain*fzintg(k,j,i)
                    
                else if (imatrix(k,j,i) .eq. 6) then ! Buffer zone
                    ! NOTE: initu and initw should technically be transformed first, but since
                    ! the buffer region has only been used for x- and z-constant flows, using
                    ! them as they are should be fine
                    fxintg(k,j,i) = fxintg(k,j,i) + up(j,i) - initu(k,j,i)
                    fyintg(k,j,i) = fyintg(k,j,i) + vp(j,i) - initv(k,j,i)
                    fzintg(k,j,i) = fzintg(k,j,i) + wp(j,i) - initw(k,j,i)
                    vwx(j,i) = vwx(j,i) + (-bfugain(i-bfhead+1)*(up(j,i) - initu(k,j,i)) - (bfgain(i-bfhead+1) * fxintg(k,j,i)))
                    vwy(j,i) = vwy(j,i) + (-bfugain(i-bfhead+1)*(vp(j,i) - initv(k,j,i)) - (bfgain(i-bfhead+1) * fyintg(k,j,i)))
                    vwz(j,i) = vwz(j,i) + (-bfugain(i-bfhead+1)*(wp(j,i) - initw(k,j,i)) - (bfgain(i-bfhead+1) * fzintg(k,j,i)))
    
                else if (imatrix(k,j,i) .eq. 7) then ! Suction region
                    fyintg(k,j,i) = fyintg(k,j,i) + (vp(j,i)-vdes(i))
                    vwy(j,i) = vwy(j,i) + (-gain*fyintg(k,j,i)) + (-ugain*(vp(j,i)-vdes(i)))
                    fxintg(k,j,i) = fxintg(k,j,i) + (up(j,i)-uinf)
                    vwx(j,i) = vwx(j,i) + (-gain*fxintg(k,j,i)) + (-ugain*(up(j,i)-uinf))
                    fzintg(k,j,i) = fzintg(k,j,i) + (wp(j,i)-initw(k,j,i))
                    vwz(j,i) = vwz(j,i) + (-gain*fzintg(k,j,i)) + (-ugain*(wp(j,i)-initw(k,j,i)))
     
                else if (imatrix(k,j,i) .eq. 11) then ! Constant x-force between wall and suction layer
                    dragx(j,i) = 1.0
    
                else if (imatrix(k,j,i) .eq. 0 .and. flow_select .eq. 1) then ! Constant pressure gradient
                    vwx(j,i) = vwx(j,i) - dPdx 
                
                else if (flow_select .eq. 5 .and. it .le. perturbtime) then ! Vortex ring
                    zj = float(j-1)*delzm
                    rsq = (ycoord(k) - ycenter)**2 + (zj - zcenter)**2
                    argrad = forbeta*(rad - sqrt(rsq))
                    fr = 0.5*(1.0 + tanh(argrad))
    
                    xi = float(i-1)*delxm
                    argx = forbeta*(L - abs(xi - xcenter))
                    fx = 0.5*(1.0 + tanh(argx))
  
                    vwx(j,i) = vwx(j,i) + bdyfx*fx*fr
    
                end if

                if (scl_flag .ge. 2) then
                    vc(j,i) = vc(j,i) + scsource(k,j,i)/dt
                end if
            end do ! i
        end do ! j
    
        !---------------------------------------------------------------------!
        !            Apply the force field to the solid surface               !
        !---------------------------------------------------------------------!
        if(k .le. kmaxsurf) then
            do i = 1,mx
                do ii = -2, 2
                    ipii = i + ii
                    idif = 1 + iabs(ii)        
                    if(ipii .lt. 1 ) ipii = ipii + mx         
                    if(ipii .gt. mx) ipii = ipii - mx
                    do jj = -2, 2
                        jdif = 1 + iabs(jj)
                        segdrag = fspread(jdif,idif)
    
                        do j = 3, mzm2
                            jpjj = j + jj
                            xsegdrag = segdrag*dragx(j,i)
                            ysegdrag = segdrag*dragy(j,i)
                            zsegdrag = segdrag*dragz(j,i)
                            vwx(jpjj,ipii) = vwx(jpjj,ipii) + xsegdrag
                            vwy(jpjj,ipii) = vwy(jpjj,ipii) + ysegdrag
                            vwz(jpjj,ipii) = vwz(jpjj,ipii) + zsegdrag
                        end do
    
                        do j = 1, 2
                            jpjj = j + jj
                            if(jpjj .lt. 1 ) jpjj = jpjj + mz
                            xsegdrag = segdrag*dragx(j,i)
                            ysegdrag = segdrag*dragy(j,i)
                            zsegdrag = segdrag*dragz(j,i)
                            vwx(jpjj,ipii) = vwx(jpjj,ipii) + xsegdrag
                            vwy(jpjj,ipii) = vwy(jpjj,ipii) + ysegdrag
                            vwz(jpjj,ipii) = vwz(jpjj,ipii) + zsegdrag
                        end do
    
                        do j = mzm1, mz
                            jpjj = j + jj
                            if(jpjj .gt. mz) jpjj = jpjj - mz
                            xsegdrag = segdrag*dragx(j,i)
                            ysegdrag = segdrag*dragy(j,i)
                            zsegdrag = segdrag*dragz(j,i)
                            vwx(jpjj,ipii) = vwx(jpjj,ipii) + xsegdrag
                            vwy(jpjj,ipii) = vwy(jpjj,ipii) + ysegdrag
                            vwz(jpjj,ipii) = vwz(jpjj,ipii) + zsegdrag
                        end do
    
                    end do ! jj
                end do ! ii
            end do ! i
        end if
    
        !---------------------------------------------------------------------!
        ! Save variables to print (or use in particle integration)            !
        ! Note: There's a weird artifact of how the variables are transformed !
        ! which flips the y-direction. In order to print out correctly, all   !
        ! terms which include v or an odd number of derivatives in y must be  !
        ! made negative. I know, I don't really get it either but it works.   !
        !---------------------------------------------------------------------!
    
        up3d(k,1:mz,1:mx) = up(1:mz,1:mx)
        vp3d(k,1:mz,1:mx) = -vp(1:mz,1:mx)
        wp3d(k,1:mz,1:mx) = wp(1:mz,1:mx)
        wx3d(k,1:mz,1:mx) = -wx(1:mz,1:mx)
        wy3d(k,1:mz,1:mx) = wy(1:mz,1:mx)
        wz3d(k,1:mz,1:mx) = -wz(1:mz,1:mx)
        u11p3d(k,1:mz,1:mx) = u11p(1:mz,1:mx) 
        u12p3d(k,1:mz,1:mx) = -u12p(1:mz,1:mx) 
        u13p3d(k,1:mz,1:mx) = u13p(1:mz,1:mx) 
        u21p3d(k,1:mz,1:mx) = -u21p(1:mz,1:mx) 
        u22p3d(k,1:mz,1:mx) = u22p(1:mz,1:mx) 
        u23p3d(k,1:mz,1:mx) = -u23p(1:mz,1:mx) 
        u31p3d(k,1:mz,1:mx) = u31p(1:mz,1:mx) 
        u32p3d(k,1:mz,1:mx) = -u32p(1:mz,1:mx) 
        u33p3d(k,1:mz,1:mx) = u33p(1:mz,1:mx) 
    
        Lup3d(k,1:mz,1:mx) = Lup(1:mz,1:mx) 
        Lvp3d(k,1:mz,1:mx) = Lvp(1:mz,1:mx) 
        Lwp3d(k,1:mz,1:mx) = Lwp(1:mz,1:mx) 
   
        scp3d(k,1:mz,1:mx) = scp(1:mz,1:mx)
 
        !---------------------------------------------------------------------!
        !              Now transform vxw back to spectral space               !
        !---------------------------------------------------------------------!
    
        isgn = -1
        inc  = mzp
        jump = 1
        lot  = mz
    
        call rfftmlt(vwx,wrk,trigx32,ixfax32,inc,jump,mx,lot,isgn)
        call rfftmlt(vwy,wrk,trigx32,ixfax32,inc,jump,mx,lot,isgn)
        call rfftmlt(vwz,wrk,trigx32,ixfax32,inc,jump,mx,lot,isgn)
        call rfftmlt(vc ,wrk,trigx32,ixfax32,inc,jump,mx,lot,isgn) 
        
        ! Wipe out unwanted last mode
        do j = 1,mz
            vwx(j,2) = 0.0
            vwy(j,2) = 0.0
            vwz(j,2) = 0.0
             vc(j,2) = 0.0
        end do
        
        isgn = -1
        inc  = 1
        jump = 2*mzp
        lot  = nx/2
    
        call cfftmlt(vwx(1,1),vwx(1,2),wrk,trigz32,izfax32,inc,jump,mz,lot,isgn)
        call cfftmlt(vwy(1,1),vwy(1,2),wrk,trigz32,izfax32,inc,jump,mz,lot,isgn)
        call cfftmlt(vwz(1,1),vwz(1,2),wrk,trigz32,izfax32,inc,jump,mz,lot,isgn)
        call cfftmlt(vc(1,1) ,vc(1,2) ,wrk,trigz32,izfax32,inc,jump,mz,lot,isgn)
    
        do j = 1,nz
            if (j .le. nzh) then
                jj = j
            else if (j .gt. nzh) then
                jj = (mz-nz) + j
            end if
            do i = 1,nxh
                i1 = 2*(i-1) + 1
                i2 = 2*i
                gn(k,j,i)  = cmplx(vwx(jj,i1)*rmz,vwx(jj,i2)*rmz)
                fn(k,j,i)  = cmplx(vwy(jj,i1)*rmz,vwy(jj,i2)*rmz)
                omz(k,j,i) = cmplx(vwz(jj,i1)*rmz,vwz(jj,i2)*rmz)
                scn(k,j,i) = cmplx( vc(jj,i1)*rmz, vc(jj,i2)*rmz)
            end do
        end do
    
    end do ! k
    !$omp end parallel do
    !                -----  end of loop over normal(y) planes  -----
   
    !---------------------------------------------------------------------!
    !     Calculate swirl criterion and write data for visualization      !
    !---------------------------------------------------------------------!
   
    ! Calculate swirl each time step to calculate area of vortex core 
    !$omp parallel do private(i,j,k,swirl)
    do i = 1,nyp
        do j = 1,mz
            do k = 1,mx

                call calcswirl(u11p3d(i,j,k),u21p3d(i,j,k),u31p3d(i,j,k),u12p3d(i,j,k),u22p3d(i,j,k), &
                               u32p3d(i,j,k),u13p3d(i,j,k),u23p3d(i,j,k),u33p3d(i,j,k),swirl)
    
                swirl_3d(i,j,k) = swirl
            end do
        end do
    end do
    !$omp end parallel do

    call vortArea(swirl_3d,vArea) ! Calculate area of vortex based on swirl criterion

    if (print3d .ne. 0) then
    
        if ((mod(it,iprnfrq) .eq. 0 .and. it .ne. 0) .or. it .eq. 1) then
            ! Calculate swirl
!            write(*,*) 'Calculating swirl and printing output data...'

 
            ! Write output files
            if (print3d .eq. 1) then
                call write_flowfield_ascii(up3d,vp3d,wp3d,wx3d,wy3d,wz3d,swirl_3d,real(imatrix))
            else if (print3d .eq. 3) then
                if (scl_flag .ne. 0) then
                    call write_flowfield_plt(up3d,vp3d,wp3d,wx3d,wy3d,wz3d,swirl_3d,real(imatrix),scp3d)
                else
                    call write_flowfield_plt(up3d,vp3d,wp3d,wx3d,wy3d,wz3d,swirl_3d,real(imatrix))
                end if
            else
                write(*,*) 'Warning: Unrecognized print type. No output data will be written.'
            end if
        
!            write(*,*) '    Done!'
        end if
    end if
    
    !---------------------------------------------------------------------!
    !                   Particle Tracking Implementation                  !
    !---------------------------------------------------------------------!
    
    if (npart .ne. 0) then
        call part_track(up3d,vp3d,wp3d,wx3d,wy3d,wz3d,u_old,v_old,w_old,  &
                        u11p3d,u12p3d,u13p3d,u21p3d,u22p3d,u23p3d,u31p3d, &
                        u32p3d,u33p3d,Lup3d,Lvp3d,Lwp3d,Lu_old,Lv_old,Lw_old,vArea)
    
        ! Save velocity and Laplacian for time derivatives inside particle integration
        u_old = up3d
        v_old = vp3d
        w_old = wp3d
    
        Lu_old = Lup3d
        Lv_old = Lvp3d
        Lw_old = Lwp3d
    end if
    
    !---------------------------------------------------------------------!
    ! Write Mean U Data, calculate mass flux, & KE/enstrophy if relevant  !
    !---------------------------------------------------------------------!
    
!    if (flow_select .eq. 1 .or. flow_select .eq. 4) then ! Only relevant for wall-bounded turbulence
!        write(*,*) 'Writing mean U data...'
!    
!        ! Mean velocity
!        if (irstrt .eq. it) then
!            open(71,file = 'outputs/mean_u_data.dat')
!        else
!            open(71,file = 'outputs/mean_u_data.dat', position = 'append')
!        end if
!    
!        do i = 1,nyp
!            do j = 1,mz
!                uxmean(j) = sum(up3d(i,j,:))/mx
!            end do
!            uzmean(i) = sum(uxmean)/mz
!            write(71,"(*(e14.6,1x))") uzmean
!        end do
!    
!        close(71)
!        write(*,*) '    Done!'
!    end if
    
    ! Compute mass flux: Integrate <U> over y
!    if (irstrt .eq. it) then
!        open(72,file='outputs/mass_flux')
!    else
!        open(72,file='outputs/mass_flux',position='append')
!    end if
!    
!    massFlux = 0.0
!    do i = 1,ny
!        massFlux = massFlux + 1.0/yl*0.5*(uzmean(i+1) + uzmean(i))*(ycoord(i+1) - ycoord(i)) ! Bulk velocity
!    end do
!    
!    write(72,*) massFlux
!    close(72)
    
    ! Compute average KE and Enstrophy for Vortex Ring
    if (mod(flow_select,5) .eq. 0) then
        sumens = 0.0
        KE = 0.0
        scl_total = 0.0
    
        do i = 1,nyp
            ! Calc y length
            if (i .eq. 1) then
                Ly = (ycoord(2) - ycoord(1))/2.0
            else if (i .eq. nyp) then
                Ly = (ycoord(nyp) - ycoord(ny))/2.0
            else
                Ly = (ycoord(i+1) - ycoord(i))/2.0 + (ycoord(i) - ycoord(i-1))/2.0
            end if
    
            do j = 1,mz
                ! Calc z length
                if (j .eq. 1 .or. j .eq. mz) then
                    Lz = delzm/2.0
                else
                    Lz = delzm
                end if
    
                do k = 1,mx
                    ! Calc x length
                    if (k .eq. 1 .or. k .eq. mx) then
                        Lx = delxm/2.0
                    else
                        Lx = delxm
                    end if
    
                    sumens = sumens + (wx3d(i,j,k)**2 + wy3d(i,j,k)**2 + wz3d(i,j,k)**2)*Lx*Ly*Lz
                    KE = KE + 0.5*(up3d(i,j,k)**2 + vp3d(i,j,k)**2 + wp3d(i,j,k)**2)*Lx*Ly*Lz
                    scl_total = scl_total + scp3d(i,j,k)*Lx*Ly*Lz
                end do
            end do
        end do
   
        if (it .eq. irstrt) then
            open(73,file='outputs/enstrophy')
            open(74,file='outputs/KE')
            open(75,file='outputs/total_scalar')
        else
            open(73,file='outputs/enstrophy',position='append')
            open(74,file='outputs/KE',position='append')
            open(75,file='outputs/total_scalar',position='append')
        end if
    
        write(73,*) sumens
        write(74,*) KE
        write(75,*) scl_total
        close(73)
        close(74)
        close(75) 
    end if
    
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
    
    end subroutine vcw3d
    
    !---------------------------------------------------------------------!
    
    subroutine vcw3dp(u,v,w,omx,omy,omz,fn,gn,u11,u12,u13,u21,u22,u23,     &
                    u31,u32,u33,Lu,Lv,Lw,scalar,sclx,scly,sclz,scn,        &
                    c11,c12,c13,c21,c22,c23,c31,c32,c33,                   &
                    dc111,dc112,dc113,dc211,dc212,dc213,dc311,dc312,dc313, &
                    dc121,dc122,dc123,dc221,dc222,dc223,dc321,dc322,dc323, &
                    dc131,dc132,dc133,dc231,dc232,dc233,dc331,dc332,dc333, &
                    c11n,c12n,c13n,c22n,c23n,c33n,                         &
                    str11n,str12n,str13n,str22n,str23n,str33n,             &
                    qp11,qp12,qp13,qp22,qp23,qp33,Lu_old,Lv_old,Lw_old,rank) 
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
    !  real space in y and Fourier in x and z. Before this subroutine is  !
    !  called, we must to a y-transformation on all arrays containing     !
    !  velocity, vorticity, or velocity derivatives. Also, we assume that !
    !  the appropriate modes in all fields (i.e., the first and last      !
    !  modes) are set to zero prior to this point.                        !
    !---------------------------------------------------------------------!
    ! =================================================================== !
    !                           Declare Modules                           !
    ! =================================================================== !
    use omp_lib
    use grid_size
    use derivs
    use helpers
    !---------------------------------------------------------------------!
    
    ! =================================================================== !
    !                   Declare all variables explicitly                  !
    ! =================================================================== !
    implicit none
    
    ! Passed variables
    complex, dimension(nyp,nz,nxh) :: u,v,w,omx,omy,omz,fn,gn
    complex, dimension(nyp,nz,nxh) :: u11,u12,u13,u21,u22,u23,u31,u32,u33
    complex, dimension(nyp,nz,nxh) :: Lu,Lv,Lw
    complex, dimension(nyp,nz,nxh) :: scalar,sclx,scly,sclz,scn
    complex, dimension(nyp,nz,nxh) :: c11,c12,c13,c21,c22,c23,c31,c32,c33
    complex, dimension(nyp,nz,nxh) :: c11n,c12n,c13n,c22n,c23n,c33n

    complex, dimension(nyp,nz,nxh) :: dc111,dc112,dc113,dc121,dc122,dc123,dc131,dc132,dc133
    complex, dimension(nyp,nz,nxh) :: dc211,dc212,dc213,dc221,dc222,dc223,dc231,dc232,dc233
    complex, dimension(nyp,nz,nxh) :: dc311,dc312,dc313,dc321,dc322,dc323,dc331,dc332,dc333

    complex, dimension(nyp,nz,nxh) :: str11n,str12n,str13n,str22n,str23n,str33n
    complex, dimension(nyp,nz,nxh) :: str11nm1,str12nm1,str13nm1,str22nm1,str23nm1,str33nm1

    complex, dimension(nyp,nz,nxh) :: qp11,qp12,qp13,qp22,qp23,qp33

    
    ! Calculation variables
    real, save, dimension(nyp,mz,mx) :: u_old,v_old,w_old
    real, dimension(nyp,mz,mx) :: Lu_old,Lv_old,Lw_old
    real, dimension(nyp,mz,mx) :: up3d,vp3d,wp3d,wx3d,wy3d,wz3d
    real, dimension(nyp,mz,mx) :: u11p3d,u12p3d,u13p3d
    real, dimension(nyp,mz,mx) :: u21p3d,u22p3d,u23p3d
    real, dimension(nyp,mz,mx) :: u31p3d,u32p3d,u33p3d
    real, dimension(nyp,mz,mx) :: Lup3d,Lvp3d,Lwp3d,swirl_3d
    real, dimension(nyp,mz,mx) :: scp3d,scsource,beta3d
    
    real, dimension(mzp,mxp2)  :: up,vp,wp,wx,wy,wz,vwx,vwy,vwz
    real, dimension(mzp,mxp2)  :: u11p,u12p,u13p
    real, dimension(mzp,mxp2)  :: u21p,u22p,u23p
    real, dimension(mzp,mxp2)  :: u31p,u32p,u33p
    real, dimension(mzp,mxp2)  :: Lup,Lvp,Lwp
    real, dimension(mzp,mxp2)  :: scp,cx,cy,cz,vc,beta_poly
    
    real, dimension(mzp,mxp2)  :: c11p,c12p,c13p
    real, dimension(mzp,mxp2)  :: c21p,c22p,c23p
    real, dimension(mzp,mxp2)  :: c31p,c32p,c33p

    real, dimension(mzp,mxp2)  :: trp

    real, dimension(mzp,mxp2)  :: dc111p, dc112p, dc113p
    real, dimension(mzp,mxp2)  :: dc211p, dc212p, dc213p
    real, dimension(mzp,mxp2)  :: dc311p, dc312p, dc313p
    
    real, dimension(mzp,mxp2)  :: dc121p, dc122p, dc123p
    real, dimension(mzp,mxp2)  :: dc221p, dc222p, dc223p
    real, dimension(mzp,mxp2)  :: dc321p, dc322p, dc323p
    
    real, dimension(mzp,mxp2)  :: dc131p, dc132p, dc133p
    real, dimension(mzp,mxp2)  :: dc231p, dc232p, dc233p
    real, dimension(mzp,mxp2)  :: dc331p, dc332p, dc333p

!    ! Enstrophy equation variables
!    real, dimension(nyp,mz,mx) :: trp3d,df1,df2,df11,df12,df22,coeff
!    real, dimension(nyp,mz,mx) :: zbeta3d,dnu1,dnu2
!    real, dimension(nyp,mz,mx) :: c11p3d,c12p3d,c22p3d
!    real, dimension(nyp,mz,mx) :: dc111p3d,dc112p3d
!    real, dimension(nyp,mz,mx) :: dc121p3d,dc122p3d
!    real, dimension(nyp,mz,mx) :: dc221p3d,dc222p3d
!    real, dimension(nyp,mz,mx) :: dc1112,dc1222,dc1211,dc2212


    real, dimension(mzp,mxp2)  :: c11np,c12np,c13np
    real, dimension(mzp,mxp2)  :: c22np,c23np,c33np
    
    real, dimension(mzp,mxp2)  :: str11np,str12np,str13np
    real, dimension(mzp,mxp2)  :: str22np,str23np,str33np
    
    real, dimension(mzp,mxp2)  :: qp11np,qp12np,qp13np
    real, dimension(mzp,mxp2)  :: qp22np,qp23np,qp33np
    
    integer :: i,j,k,i1,i2,inc,isgn,jump,lot,n,rank
    integer :: ii,jj,ipii,jpjj,idif,jdif
    real    :: delx,delz,delxm,delzm,pi,vArea
    real    :: segdrag,xsegdrag,ysegdrag,zsegdrag,wdes
    real    :: cflcheck,cflmax,swirl
    
    real,dimension(nwrk)  :: wrk
    real,dimension(nyp)   :: cfl
    real,dimension(mz,mx) :: dragx,dragy,dragz
    
    real :: xi,zj,argx,argrad,fx,fr,rsq
    real :: uxmean(mz),uzmean(nyp)
    real :: Lx,Ly,Lz,massFlux,sumens,KE,scl_total,SG
    
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

    ! Scalar variables
    real    :: sigmax,sigmay,sigmaz
    real    :: betax,betay,betaz
    real    :: xsq,ysq,zsq
    real    :: xc1,yc1,zc1
    real    :: deltaT,diff
    integer :: scl_flag
    integer :: src_start,src_stop

    ! Particle variables
    real,save,dimension(npart) :: xpart,ypart,zpart

    ! Polymer variables
    integer :: ipolyflag,itarget,ipeter
    real    :: alpha_poly,tpoly,zlmax,diffpoly,qbeta
    real    :: zbeta1

    ! 2D Vortex variables
    character(33) :: filename2
    character(10) :: dirname
    real :: r,z,z1,y1,y2,wxi,sci,argy,scj,wxj,up1
    real :: vortGamma,vortSigma,vortY,vortZ

    !---------------------------------------------------------------------!
    
    ! =================================================================== !
    !                            Common Blocks                            !
    ! =================================================================== !
    common/iocontrl/   irstrt,nsteps,iprnfrq,print3d
    common/solver/     gain,ugain,theta,alpha,beta,dyde
    common/ibforce/    fxintg,fyintg,fzintg,fspread
    common/waves/      wavx,wavz,c
    common/init/       initu,initv,initw
    common/buffer/     bfgain,bfugain,vdes,bfhead,bftail,bfwidth,slopelength
    common/flow/       re,Uinf,R_tau,dPdx
    common/domain/     xl,yl,zl
    common/trig/       trigx,trigy,trigz,sine,cosine,trigx32,trigz32
    common/ifax/       ixfax,iyfax,izfax,ixfax32,izfax32
    common/itime/      it,dt
    common/imat/       imatrix,kwall,kmaxsurf
    common/vortexring/ forbeta,xcenter,ycenter,zcenter,L,rad,bdyfx
    common/setup/      geomtype,flow_select,perturbtime
    common/scl_stuff/  sigmax,sigmay,sigmaz,deltaT,diff,scl_flag
    common/part_traj/  xpart,ypart,zpart
    common/src_time/   src_start,src_stop
    common/poly_flgs/  ipolyflag,itarget,ipeter
    common/poly_var/   alpha_poly,tpoly,zlmax,diffpoly,qbeta
    common/vort/       vortGamma,vortSigma,vortY,vortZ
    !---------------------------------------------------------------------!
    
    ! =================================================================== !
    !                         Begin Calculations                          !
    ! =================================================================== !
    
    !---------------------------------------------------------------------!
    !                    Initialize some variables                        !
    !---------------------------------------------------------------------!
   
 
    pi = 2.0*acos(0.0)
    delx = xl/float(nx)
    delz = zl/float(nz)
    delxm = xl/float(mx)
    delzm = zl/float(mz)
    
    cfl = 0.0
    cflmax = 0.0
    
    !---------------------------------------------------------------------!
    !                    Calculate scalar source term                     !
    !---------------------------------------------------------------------!
    if (scl_flag .eq. 2 .and. it .ge. src_start .and. it .le. src_stop) then
        scsource = 0.0
        do n = 1,npart
        if (it .eq. 0) then
            open(95,file='setup/particles/particles.dat',status='old',action='read')
            do j = 1,n
                read(95,*)
            end do
            read(95,*) xc1,yc1,zc1
            close(95)
        else
            xc1 = xpart(n)
            yc1 = ypart(n)
            zc1 = zpart(n)
        end if
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
        end do
    else if (scl_flag .eq. 3 .and. it .eq. 1 .and. irstrt .eq. 1) then
        scsource = 1.0
    else
        scsource = 0.0
    end if

    !---------------------------------------------------------------------!
    !         Calculate (u x omega), iterating over y-planes              !
    !---------------------------------------------------------------------!
   
    !$omp parallel do default(shared) private(i,j,k,ii,i1,i2,ipii,idif,jj,   &
    !$omp    jpjj,jdif,segdrag,xsegdrag,ysegdrag,zsegdrag,cflcheck,          &
    !$omp    up,vp,wp,wx,wy,wz,vwx,vwy,vwz,dragx,dragy,dragz,wrk,inc,isgn,   &
    !$omp    u11p,u12p,u13p,u21p,u22p,u23p,u31p,u32p,u33p,Lup,Lvp,Lwp,       &
    !$omp    jump,lot,fx,fr,xi,zj,argrad,rsq,argx,scp,cx,cy,cz,vc,beta_poly, &
    !$omp    c11p,c12p,c13p,c21p,c22p,c23p,c31p,c32p,c33p,                   &
    !$omp    dc111p,dc112p,dc113p,dc211p,dc212p,dc213p,                      &
    !$omp    dc311p,dc312p,dc313p,dc121p,dc122p,dc123p,dc221p,dc222p,dc223p, &
    !$omp    dc321p,dc322p,dc323p,dc131p,dc132p,dc133p,dc231p,dc232p,dc233p, &
    !$omp    dc331p,dc332p,dc333p,trp,c11np,c12np,c13np,c22np,c23np,c33np,   &
    !$omp    str11np,str12np,str13np,str22np,str23np,str33np,                &
    !$omp    qp11np,qp12np,qp13np,qp22np,qp23np,qp33np,zbeta1) schedule(dynamic)
    
    do k = 1,nyp
        ! Initialize the velocities and vorticities on the kth x-z plane 
        ! with zeros. We need to pad the extra 1/2 modes with zeros by 3/2 rule
    
        do j = 1,mzp
            do i = 1,mxp2
    
                up(j,i) = 0.0
                vp(j,i) = 0.0
                wp(j,i) = 0.0
    
                wx(j,i) = 0.0
                wy(j,i) = 0.0
                wz(j,i) = 0.0
    
                vwx(j,i) = 0.0
                vwy(j,i) = 0.0
                vwz(j,i) = 0.0
    
                ! Used in particle tracking
                u11p(j,i) = 0.0
                u12p(j,i) = 0.0
                u13p(j,i) = 0.0
                u21p(j,i) = 0.0
                u22p(j,i) = 0.0
                u23p(j,i) = 0.0
                u31p(j,i) = 0.0
                u32p(j,i) = 0.0
                u33p(j,i) = 0.0
          
                Lup(j,i) = 0.0
                Lvp(j,i) = 0.0
                Lwp(j,i) = 0.0  
   
                scp(j,i) = 0.0
                 cx(j,i) = 0.0 
                 cy(j,i) = 0.0 
                 cz(j,i) = 0.0 
                 vc(j,i) = 0.0 

                if (it .ge. src_start-1) then
                c11p(j,i) = 0.0
                c12p(j,i) = 0.0
                c13p(j,i) = 0.0
                c21p(j,i) = 0.0
                c22p(j,i) = 0.0
                c23p(j,i) = 0.0
                c31p(j,i) = 0.0
                c32p(j,i) = 0.0
                c33p(j,i) = 0.0
                
                dc111p(j,i) = 0.0
                dc112p(j,i) = 0.0
                dc113p(j,i) = 0.0 
                dc211p(j,i) = 0.0
                dc212p(j,i) = 0.0
                dc213p(j,i) = 0.0
                dc311p(j,i) = 0.0
                dc312p(j,i) = 0.0
                dc313p(j,i) = 0.0
            
                dc121p(j,i) = 0.0
                dc122p(j,i) = 0.0
                dc123p(j,i) = 0.0 
                dc221p(j,i) = 0.0
                dc222p(j,i) = 0.0
                dc223p(j,i) = 0.0
                dc321p(j,i) = 0.0
                dc322p(j,i) = 0.0
                dc323p(j,i) = 0.0  
            
                dc131p(j,i) = 0.0
                dc132p(j,i) = 0.0
                dc133p(j,i) = 0.0 
                dc231p(j,i) = 0.0
                dc232p(j,i) = 0.0
                dc233p(j,i) = 0.0
                dc331p(j,i) = 0.0
                dc332p(j,i) = 0.0
                dc333p(j,i) = 0.0
            
                trp(j,i) = 0.0
            
                c11np(j,i) = 0.0
                c12np(j,i) = 0.0
                c13np(j,i) = 0.0
                c22np(j,i) = 0.0
                c23np(j,i) = 0.0
                c33np(j,i) = 0.0
            
                str11np(j,i) = 0.0
                str12np(j,i) = 0.0
                str13np(j,i) = 0.0
                str22np(j,i) = 0.0
                str23np(j,i) = 0.0
                str33np(j,i) = 0.0
            
                qp11np(j,i) = 0.0
                qp12np(j,i) = 0.0
                qp13np(j,i) = 0.0
                qp22np(j,i) = 0.0
                qp23np(j,i) = 0.0
                qp33np(j,i) = 0.0
                end if
            end do
        end do
    
        do j = 1,nz
            if (j .le. nzh) jj = j
            if (j .gt. nzh) jj = (mz-nz) + j
            
            do i = 1,nxh
                i1 = 2*(i-1) + 1
                i2 = 2*i
    
                up(jj,i1) =  real(u(k,j,i))
                up(jj,i2) = aimag(u(k,j,i))
                vp(jj,i1) =  real(v(k,j,i))
                vp(jj,i2) = aimag(v(k,j,i))
                wp(jj,i1) =  real(w(k,j,i))
                wp(jj,i2) = aimag(w(k,j,i))
    
                wx(jj,i1) =  real(omx(k,j,i))
                wx(jj,i2) = aimag(omx(k,j,i))
                wy(jj,i1) =  real(omy(k,j,i))
                wy(jj,i2) = aimag(omy(k,j,i))
                wz(jj,i1) =  real(omz(k,j,i))
                wz(jj,i2) = aimag(omz(k,j,i))
    
                u11p(jj,i1) =  real(u11(k,j,i))
                u11p(jj,i2) = aimag(u11(k,j,i))
                u12p(jj,i1) =  real(u12(k,j,i))
                u12p(jj,i2) = aimag(u12(k,j,i))
                u13p(jj,i1) =  real(u13(k,j,i))
                u13p(jj,i2) = aimag(u13(k,j,i))
                u21p(jj,i1) =  real(u21(k,j,i))
                u21p(jj,i2) = aimag(u21(k,j,i))
                u22p(jj,i1) =  real(u22(k,j,i))
                u22p(jj,i2) = aimag(u22(k,j,i))
                u23p(jj,i1) =  real(u23(k,j,i))
                u23p(jj,i2) = aimag(u23(k,j,i))
                u31p(jj,i1) =  real(u31(k,j,i))
                u31p(jj,i2) = aimag(u31(k,j,i))
                u32p(jj,i1) =  real(u32(k,j,i))
                u32p(jj,i2) = aimag(u32(k,j,i))
                u33p(jj,i1) =  real(u33(k,j,i))
                u33p(jj,i2) = aimag(u33(k,j,i))
                
                Lup(jj,i1) =  real(Lu(k,j,i))
                Lup(jj,i2) = aimag(Lu(k,j,i))
                Lvp(jj,i1) =  real(Lv(k,j,i))
                Lvp(jj,i2) = aimag(Lv(k,j,i))
                Lwp(jj,i1) =  real(Lw(k,j,i))
                Lwp(jj,i2) = aimag(Lw(k,j,i))
    
                scp(jj,i1) =  real(scalar(k,j,i))
                scp(jj,i2) = aimag(scalar(k,j,i))
                cx(jj,i1)  =  real(sclx(k,j,i))
                cx(jj,i2)  = aimag(sclx(k,j,i))
                cy(jj,i1)  =  real(scly(k,j,i))
                cy(jj,i2)  = aimag(scly(k,j,i))
                cz(jj,i1)  =  real(sclz(k,j,i))
                cz(jj,i2)  = aimag(sclz(k,j,i))

                if (it .ge. src_start-1) then
                c11p(jj,i1) = real(c11(k,j,i))
                c12p(jj,i1) = real(c12(k,j,i))
                c13p(jj,i1) = real(c13(k,j,i))
                c21p(jj,i1) = real(c21(k,j,i))
                c22p(jj,i1) = real(c22(k,j,i))
                c23p(jj,i1) = real(c23(k,j,i))
                c31p(jj,i1) = real(c31(k,j,i))
                c32p(jj,i1) = real(c32(k,j,i))
                c33p(jj,i1) = real(c33(k,j,i))
            
                dc111p(jj,i1) = real(dc111(k,j,i))
                dc112p(jj,i1) = real(dc112(k,j,i))
                dc113p(jj,i1) = real(dc113(k,j,i))
                dc211p(jj,i1) = real(dc211(k,j,i))
                dc212p(jj,i1) = real(dc212(k,j,i))
                dc213p(jj,i1) = real(dc213(k,j,i))
                dc311p(jj,i1) = real(dc311(k,j,i))
                dc312p(jj,i1) = real(dc312(k,j,i))
                dc313p(jj,i1) = real(dc313(k,j,i))
            
                dc121p(jj,i1) = real(dc121(k,j,i))
                dc122p(jj,i1) = real(dc122(k,j,i))
                dc123p(jj,i1) = real(dc123(k,j,i))
                dc221p(jj,i1) = real(dc221(k,j,i))
                dc222p(jj,i1) = real(dc222(k,j,i))
                dc223p(jj,i1) = real(dc223(k,j,i))
                dc321p(jj,i1) = real(dc321(k,j,i))
                dc322p(jj,i1) = real(dc322(k,j,i))
                dc323p(jj,i1) = real(dc323(k,j,i))
            
                dc131p(jj,i1) = real(dc131(k,j,i))
                dc132p(jj,i1) = real(dc132(k,j,i))
                dc133p(jj,i1) = real(dc133(k,j,i))
                dc231p(jj,i1) = real(dc231(k,j,i))
                dc232p(jj,i1) = real(dc232(k,j,i))
                dc233p(jj,i1) = real(dc233(k,j,i))
                dc331p(jj,i1) = real(dc331(k,j,i))
                dc332p(jj,i1) = real(dc332(k,j,i))
                dc333p(jj,i1) = real(dc333(k,j,i))

                c11p(jj,i2) = aimag(c11(k,j,i))
                c12p(jj,i2) = aimag(c12(k,j,i))
                c13p(jj,i2) = aimag(c13(k,j,i))
                c21p(jj,i2) = aimag(c21(k,j,i))
                c22p(jj,i2) = aimag(c22(k,j,i))
                c23p(jj,i2) = aimag(c23(k,j,i))
                c31p(jj,i2) = aimag(c31(k,j,i))
                c32p(jj,i2) = aimag(c32(k,j,i))
                c33p(jj,i2) = aimag(c33(k,j,i))
           
                dc111p(jj,i2) = aimag(dc111(k,j,i))
                dc112p(jj,i2) = aimag(dc112(k,j,i))
                dc113p(jj,i2) = aimag(dc113(k,j,i))
                dc211p(jj,i2) = aimag(dc211(k,j,i))
                dc212p(jj,i2) = aimag(dc212(k,j,i))
                dc213p(jj,i2) = aimag(dc213(k,j,i))
                dc311p(jj,i2) = aimag(dc311(k,j,i))
                dc312p(jj,i2) = aimag(dc312(k,j,i))
                dc313p(jj,i2) = aimag(dc313(k,j,i))
           
                dc121p(jj,i2) = aimag(dc121(k,j,i))
                dc122p(jj,i2) = aimag(dc122(k,j,i))
                dc123p(jj,i2) = aimag(dc123(k,j,i))
                dc221p(jj,i2) = aimag(dc221(k,j,i))
                dc222p(jj,i2) = aimag(dc222(k,j,i))
                dc223p(jj,i2) = aimag(dc223(k,j,i))
                dc321p(jj,i2) = aimag(dc321(k,j,i))
                dc322p(jj,i2) = aimag(dc322(k,j,i))
                dc323p(jj,i2) = aimag(dc323(k,j,i))
           
                dc131p(jj,i2) = aimag(dc131(k,j,i))
                dc132p(jj,i2) = aimag(dc132(k,j,i))
                dc133p(jj,i2) = aimag(dc133(k,j,i))
                dc231p(jj,i2) = aimag(dc231(k,j,i))
                dc232p(jj,i2) = aimag(dc232(k,j,i))
                dc233p(jj,i2) = aimag(dc233(k,j,i))
                dc331p(jj,i2) = aimag(dc331(k,j,i))
                dc332p(jj,i2) = aimag(dc332(k,j,i))
                dc333p(jj,i2) = aimag(dc333(k,j,i))
                end if
            end do
        end do
   
    !---------------------------------------------------------------------!
    !        Transform (interpolate) into 3/2 grid physical space         !
    !---------------------------------------------------------------------!
    
        inc  = 1
        isgn = 1
        jump = 2*mzp
        lot  = nx/2
    
        call cfftmlt(up(1,1),up(1,2),wrk,trigz32,izfax32,inc,jump,mz,lot,isgn)
        call cfftmlt(vp(1,1),vp(1,2),wrk,trigz32,izfax32,inc,jump,mz,lot,isgn)
        call cfftmlt(wp(1,1),wp(1,2),wrk,trigz32,izfax32,inc,jump,mz,lot,isgn)
        call cfftmlt(wx(1,1),wx(1,2),wrk,trigz32,izfax32,inc,jump,mz,lot,isgn)
        call cfftmlt(wy(1,1),wy(1,2),wrk,trigz32,izfax32,inc,jump,mz,lot,isgn)
        call cfftmlt(wz(1,1),wz(1,2),wrk,trigz32,izfax32,inc,jump,mz,lot,isgn)
    
        call cfftmlt(u11p(1,1),u11p(1,2),wrk,trigz32,izfax32,inc,jump,mz,lot,isgn)
        call cfftmlt(u12p(1,1),u12p(1,2),wrk,trigz32,izfax32,inc,jump,mz,lot,isgn)
        call cfftmlt(u13p(1,1),u13p(1,2),wrk,trigz32,izfax32,inc,jump,mz,lot,isgn)
        call cfftmlt(u21p(1,1),u21p(1,2),wrk,trigz32,izfax32,inc,jump,mz,lot,isgn)
        call cfftmlt(u22p(1,1),u22p(1,2),wrk,trigz32,izfax32,inc,jump,mz,lot,isgn)
        call cfftmlt(u23p(1,1),u23p(1,2),wrk,trigz32,izfax32,inc,jump,mz,lot,isgn)
        call cfftmlt(u31p(1,1),u31p(1,2),wrk,trigz32,izfax32,inc,jump,mz,lot,isgn)
        call cfftmlt(u32p(1,1),u32p(1,2),wrk,trigz32,izfax32,inc,jump,mz,lot,isgn)
        call cfftmlt(u33p(1,1),u33p(1,2),wrk,trigz32,izfax32,inc,jump,mz,lot,isgn)
     
        call cfftmlt(Lup(1,1),Lup(1,2),wrk,trigz32,izfax32,inc,jump,mz,lot,isgn)
        call cfftmlt(Lvp(1,1),Lvp(1,2),wrk,trigz32,izfax32,inc,jump,mz,lot,isgn)
        call cfftmlt(Lwp(1,1),Lwp(1,2),wrk,trigz32,izfax32,inc,jump,mz,lot,isgn)
    
        call cfftmlt(scp(1,1),scp(1,2),wrk,trigz32,izfax32,inc,jump,mz,lot,isgn)
        call cfftmlt(cx(1,1),cx(1,2),wrk,trigz32,izfax32,inc,jump,mz,lot,isgn)
        call cfftmlt(cy(1,1),cy(1,2),wrk,trigz32,izfax32,inc,jump,mz,lot,isgn)
        call cfftmlt(cz(1,1),cz(1,2),wrk,trigz32,izfax32,inc,jump,mz,lot,isgn)

        if (it .ge. (src_start - 1)) then
        call cfftmlt(c11p(1,1),c11p(1,2),wrk,trigz32,izfax32,inc,jump,mz,lot,isgn)
        call cfftmlt(c12p(1,1),c12p(1,2),wrk,trigz32,izfax32,inc,jump,mz,lot,isgn)
        call cfftmlt(c13p(1,1),c13p(1,2),wrk,trigz32,izfax32,inc,jump,mz,lot,isgn)
        call cfftmlt(c21p(1,1),c21p(1,2),wrk,trigz32,izfax32,inc,jump,mz,lot,isgn)
        call cfftmlt(c22p(1,1),c22p(1,2),wrk,trigz32,izfax32,inc,jump,mz,lot,isgn)
        call cfftmlt(c23p(1,1),c23p(1,2),wrk,trigz32,izfax32,inc,jump,mz,lot,isgn)
        call cfftmlt(c31p(1,1),c31p(1,2),wrk,trigz32,izfax32,inc,jump,mz,lot,isgn)
        call cfftmlt(c32p(1,1),c32p(1,2),wrk,trigz32,izfax32,inc,jump,mz,lot,isgn)
        call cfftmlt(c33p(1,1),c33p(1,2),wrk,trigz32,izfax32,inc,jump,mz,lot,isgn)
     
        call cfftmlt(dc111p(1,1),dc111p(1,2),wrk,trigz32,izfax32,inc,jump,mz,lot,isgn)
        call cfftmlt(dc112p(1,1),dc112p(1,2),wrk,trigz32,izfax32,inc,jump,mz,lot,isgn)
        call cfftmlt(dc113p(1,1),dc113p(1,2),wrk,trigz32,izfax32,inc,jump,mz,lot,isgn)
        call cfftmlt(dc211p(1,1),dc211p(1,2),wrk,trigz32,izfax32,inc,jump,mz,lot,isgn)
        call cfftmlt(dc212p(1,1),dc212p(1,2),wrk,trigz32,izfax32,inc,jump,mz,lot,isgn)
        call cfftmlt(dc213p(1,1),dc213p(1,2),wrk,trigz32,izfax32,inc,jump,mz,lot,isgn)
        call cfftmlt(dc311p(1,1),dc311p(1,2),wrk,trigz32,izfax32,inc,jump,mz,lot,isgn)
        call cfftmlt(dc312p(1,1),dc312p(1,2),wrk,trigz32,izfax32,inc,jump,mz,lot,isgn)
        call cfftmlt(dc313p(1,1),dc313p(1,2),wrk,trigz32,izfax32,inc,jump,mz,lot,isgn)
     
        call cfftmlt(dc121p(1,1),dc121p(1,2),wrk,trigz32,izfax32,inc,jump,mz,lot,isgn)
        call cfftmlt(dc122p(1,1),dc122p(1,2),wrk,trigz32,izfax32,inc,jump,mz,lot,isgn)
        call cfftmlt(dc123p(1,1),dc123p(1,2),wrk,trigz32,izfax32,inc,jump,mz,lot,isgn)
        call cfftmlt(dc221p(1,1),dc221p(1,2),wrk,trigz32,izfax32,inc,jump,mz,lot,isgn)
        call cfftmlt(dc222p(1,1),dc222p(1,2),wrk,trigz32,izfax32,inc,jump,mz,lot,isgn)
        call cfftmlt(dc223p(1,1),dc223p(1,2),wrk,trigz32,izfax32,inc,jump,mz,lot,isgn)
        call cfftmlt(dc321p(1,1),dc321p(1,2),wrk,trigz32,izfax32,inc,jump,mz,lot,isgn)
        call cfftmlt(dc322p(1,1),dc322p(1,2),wrk,trigz32,izfax32,inc,jump,mz,lot,isgn)
        call cfftmlt(dc323p(1,1),dc323p(1,2),wrk,trigz32,izfax32,inc,jump,mz,lot,isgn)
     
        call cfftmlt(dc131p(1,1),dc131p(1,2),wrk,trigz32,izfax32,inc,jump,mz,lot,isgn)
        call cfftmlt(dc132p(1,1),dc132p(1,2),wrk,trigz32,izfax32,inc,jump,mz,lot,isgn)
        call cfftmlt(dc133p(1,1),dc133p(1,2),wrk,trigz32,izfax32,inc,jump,mz,lot,isgn)
        call cfftmlt(dc231p(1,1),dc231p(1,2),wrk,trigz32,izfax32,inc,jump,mz,lot,isgn)
        call cfftmlt(dc232p(1,1),dc232p(1,2),wrk,trigz32,izfax32,inc,jump,mz,lot,isgn)
        call cfftmlt(dc233p(1,1),dc233p(1,2),wrk,trigz32,izfax32,inc,jump,mz,lot,isgn)
        call cfftmlt(dc331p(1,1),dc331p(1,2),wrk,trigz32,izfax32,inc,jump,mz,lot,isgn)
        call cfftmlt(dc332p(1,1),dc332p(1,2),wrk,trigz32,izfax32,inc,jump,mz,lot,isgn)
        call cfftmlt(dc333p(1,1),dc333p(1,2),wrk,trigz32,izfax32,inc,jump,mz,lot,isgn)
        end if

        do j = 1,mz
            up(j,nxp) = up(j,2)
            vp(j,nxp) = vp(j,2)
            wp(j,nxp) = wp(j,2)
            wx(j,nxp) = wx(j,2)
            wy(j,nxp) = wy(j,2)
            wz(j,nxp) = wz(j,2)
            
            u11p(j,nxp) = u11p(j,2)
            u12p(j,nxp) = u12p(j,2)
            u13p(j,nxp) = u13p(j,2)
            u21p(j,nxp) = u21p(j,2)
            u22p(j,nxp) = u22p(j,2)
            u23p(j,nxp) = u23p(j,2)
            u31p(j,nxp) = u31p(j,2)
            u32p(j,nxp) = u32p(j,2)
            u33p(j,nxp) = u33p(j,2)
            
            Lup(j,nxp) = Lup(j,2) 
            Lvp(j,nxp) = Lvp(j,2) 
            Lwp(j,nxp) = Lwp(j,2) 
    
            scp(j,nxp) = scp(j,2)
             cx(j,nxp) =  cx(j,2)
             cy(j,nxp) =  cy(j,2)
             cz(j,nxp) =  cz(j,2)
            
            up(j,2) = 0.0
            vp(j,2) = 0.0
            wp(j,2) = 0.0
            wx(j,2) = 0.0
            wy(j,2) = 0.0
            wz(j,2) = 0.0
            
            u11p(j,2) = 0.0
            u12p(j,2) = 0.0
            u13p(j,2) = 0.0
            u21p(j,2) = 0.0
            u22p(j,2) = 0.0
            u23p(j,2) = 0.0
            u31p(j,2) = 0.0
            u32p(j,2) = 0.0
            u33p(j,2) = 0.0
            
            Lup(j,2) = 0.0
            Lvp(j,2) = 0.0
            Lwp(j,2) = 0.0

            scp(j,2) = 0.0
             cx(j,2) = 0.0
             cy(j,2) = 0.0
             cz(j,2) = 0.0
            
            up(j,nxp2) = up(j,2)
            vp(j,nxp2) = vp(j,2)
            wp(j,nxp2) = wp(j,2)
            wx(j,nxp2) = wx(j,2)
            wy(j,nxp2) = wy(j,2)
            wz(j,nxp2) = wz(j,2)
            
            u11p(j,nxp2) = u11p(j,2)
            u12p(j,nxp2) = u12p(j,2)
            u13p(j,nxp2) = u13p(j,2)
            u21p(j,nxp2) = u21p(j,2)
            u22p(j,nxp2) = u22p(j,2)
            u23p(j,nxp2) = u23p(j,2)
            u31p(j,nxp2) = u31p(j,2)
            u32p(j,nxp2) = u32p(j,2)
            u33p(j,nxp2) = u33p(j,2)
            
            Lup(j,nxp2) = Lup(j,2) 
            Lvp(j,nxp2) = Lvp(j,2) 
            Lwp(j,nxp2) = Lwp(j,2) 

            scp(j,nxp2) = scp(j,2)
             cx(j,nxp2) =  cx(j,2)
             cy(j,nxp2) =  cy(j,2)
             cz(j,nxp2) =  cz(j,2)

            if (it .ge. src_start-1) then
            c11p(j,nxp) = c11p(j,2)
            c12p(j,nxp) = c12p(j,2)
            c13p(j,nxp) = c13p(j,2)
            c21p(j,nxp) = c21p(j,2)
            c22p(j,nxp) = c22p(j,2)
            c23p(j,nxp) = c23p(j,2)
            c31p(j,nxp) = c31p(j,2)
            c32p(j,nxp) = c32p(j,2)
            c33p(j,nxp) = c33p(j,2)
                           
            dc111p(j,nxp) = dc111p(j,2)
            dc112p(j,nxp) = dc112p(j,2)
            dc113p(j,nxp) = dc113p(j,2)
            dc211p(j,nxp) = dc211p(j,2)
            dc212p(j,nxp) = dc212p(j,2)
            dc213p(j,nxp) = dc213p(j,2)
            dc311p(j,nxp) = dc311p(j,2)
            dc312p(j,nxp) = dc312p(j,2)
            dc313p(j,nxp) = dc313p(j,2)
                                     
            dc121p(j,nxp) = dc121p(j,2)
            dc122p(j,nxp) = dc122p(j,2)
            dc123p(j,nxp) = dc123p(j,2)
            dc221p(j,nxp) = dc221p(j,2)
            dc222p(j,nxp) = dc222p(j,2)
            dc223p(j,nxp) = dc223p(j,2)
            dc321p(j,nxp) = dc321p(j,2)
            dc322p(j,nxp) = dc322p(j,2)
            dc323p(j,nxp) = dc323p(j,2) 
                                    
            dc131p(j,nxp) = dc131p(j,2)
            dc132p(j,nxp) = dc132p(j,2)
            dc133p(j,nxp) = dc133p(j,2)
            dc231p(j,nxp) = dc231p(j,2)
            dc232p(j,nxp) = dc232p(j,2)
            dc233p(j,nxp) = dc233p(j,2)
            dc331p(j,nxp) = dc331p(j,2)
            dc332p(j,nxp) = dc332p(j,2)
            dc333p(j,nxp) = dc333p(j,2)

            c11p(j,2) = 0.0
            c12p(j,2) = 0.0
            c13p(j,2) = 0.0
            c21p(j,2) = 0.0
            c22p(j,2) = 0.0
            c23p(j,2) = 0.0
            c31p(j,2) = 0.0
            c32p(j,2) = 0.0
            c33p(j,2) = 0.0
             
            dc111p(j,2) = 0.0
            dc112p(j,2) = 0.0
            dc113p(j,2) = 0.0 
            dc211p(j,2) = 0.0
            dc212p(j,2) = 0.0
            dc213p(j,2) = 0.0
            dc311p(j,2) = 0.0
            dc312p(j,2) = 0.0
            dc313p(j,2) = 0.0
     
            dc121p(j,2) = 0.0
            dc122p(j,2) = 0.0
            dc123p(j,2) = 0.0 
            dc221p(j,2) = 0.0
            dc222p(j,2) = 0.0
            dc223p(j,2) = 0.0
            dc321p(j,2) = 0.0
            dc322p(j,2) = 0.0
            dc323p(j,2) = 0.0  
     
            dc131p(j,2) = 0.0
            dc132p(j,2) = 0.0
            dc133p(j,2) = 0.0 
            dc231p(j,2) = 0.0
            dc232p(j,2) = 0.0
            dc233p(j,2) = 0.0
            dc331p(j,2) = 0.0
            dc332p(j,2) = 0.0
            dc333p(j,2) = 0.0

            c11p(j,nxp2) = c11p(j,2)
            c12p(j,nxp2) = c12p(j,2)
            c13p(j,nxp2) = c13p(j,2)
            c21p(j,nxp2) = c21p(j,2)
            c22p(j,nxp2) = c22p(j,2)
            c23p(j,nxp2) = c23p(j,2)
            c31p(j,nxp2) = c31p(j,2)
            c32p(j,nxp2) = c32p(j,2)
            c33p(j,nxp2) = c33p(j,2)
                           
            dc111p(j,nxp2) = dc111p(j,2)
            dc112p(j,nxp2) = dc112p(j,2)
            dc113p(j,nxp2) = dc113p(j,2)
            dc211p(j,nxp2) = dc211p(j,2)
            dc212p(j,nxp2) = dc212p(j,2)
            dc213p(j,nxp2) = dc213p(j,2)
            dc311p(j,nxp2) = dc311p(j,2)
            dc312p(j,nxp2) = dc312p(j,2)
            dc313p(j,nxp2) = dc313p(j,2)

            dc121p(j,nxp2) = dc121p(j,2)
            dc122p(j,nxp2) = dc122p(j,2)
            dc123p(j,nxp2) = dc123p(j,2)
            dc221p(j,nxp2) = dc221p(j,2)
            dc222p(j,nxp2) = dc222p(j,2)
            dc223p(j,nxp2) = dc223p(j,2)
            dc321p(j,nxp2) = dc321p(j,2)
            dc322p(j,nxp2) = dc322p(j,2)
            dc323p(j,nxp2) = dc323p(j,2) 

            dc131p(j,nxp2) = dc131p(j,2)
            dc132p(j,nxp2) = dc132p(j,2)
            dc133p(j,nxp2) = dc133p(j,2)
            dc231p(j,nxp2) = dc231p(j,2)
            dc232p(j,nxp2) = dc232p(j,2)
            dc233p(j,nxp2) = dc233p(j,2)
            dc331p(j,nxp2) = dc331p(j,2)
            dc332p(j,nxp2) = dc332p(j,2)
            dc333p(j,nxp2) = dc333p(j,2)
            end if

        end do
            
    
        isgn = 1
        inc  = mzp
        jump = 1
        lot  = mz
    
        call rfftmlt(up,wrk,trigx32,ixfax32,inc,jump,mx,lot,isgn)
        call rfftmlt(vp,wrk,trigx32,ixfax32,inc,jump,mx,lot,isgn)
        call rfftmlt(wp,wrk,trigx32,ixfax32,inc,jump,mx,lot,isgn)
        call rfftmlt(wx,wrk,trigx32,ixfax32,inc,jump,mx,lot,isgn)
        call rfftmlt(wy,wrk,trigx32,ixfax32,inc,jump,mx,lot,isgn)
        call rfftmlt(wz,wrk,trigx32,ixfax32,inc,jump,mx,lot,isgn)
    
        call rfftmlt(u11p,wrk,trigx32,ixfax32,inc,jump,mx,lot,isgn)
        call rfftmlt(u12p,wrk,trigx32,ixfax32,inc,jump,mx,lot,isgn)
        call rfftmlt(u13p,wrk,trigx32,ixfax32,inc,jump,mx,lot,isgn)
        call rfftmlt(u21p,wrk,trigx32,ixfax32,inc,jump,mx,lot,isgn)
        call rfftmlt(u22p,wrk,trigx32,ixfax32,inc,jump,mx,lot,isgn)
        call rfftmlt(u23p,wrk,trigx32,ixfax32,inc,jump,mx,lot,isgn)
        call rfftmlt(u31p,wrk,trigx32,ixfax32,inc,jump,mx,lot,isgn)
        call rfftmlt(u32p,wrk,trigx32,ixfax32,inc,jump,mx,lot,isgn)
        call rfftmlt(u33p,wrk,trigx32,ixfax32,inc,jump,mx,lot,isgn)
     
        call rfftmlt(Lup,wrk,trigx32,ixfax32,inc,jump,mx,lot,isgn)
        call rfftmlt(Lvp,wrk,trigx32,ixfax32,inc,jump,mx,lot,isgn)
        call rfftmlt(Lwp,wrk,trigx32,ixfax32,inc,jump,mx,lot,isgn)
    
        call rfftmlt(scp,wrk,trigx32,ixfax32,inc,jump,mx,lot,isgn)
        call rfftmlt(cx ,wrk,trigx32,ixfax32,inc,jump,mx,lot,isgn)
        call rfftmlt(cy ,wrk,trigx32,ixfax32,inc,jump,mx,lot,isgn)
        call rfftmlt(cz ,wrk,trigx32,ixfax32,inc,jump,mx,lot,isgn)

        if (it .ge. src_start-1) then
        call rfftmlt(c11p,wrk,trigx32,ixfax32,inc,jump,mx,lot,isgn)
        call rfftmlt(c12p,wrk,trigx32,ixfax32,inc,jump,mx,lot,isgn)
        call rfftmlt(c13p,wrk,trigx32,ixfax32,inc,jump,mx,lot,isgn)
        call rfftmlt(c21p,wrk,trigx32,ixfax32,inc,jump,mx,lot,isgn)
        call rfftmlt(c22p,wrk,trigx32,ixfax32,inc,jump,mx,lot,isgn)
        call rfftmlt(c23p,wrk,trigx32,ixfax32,inc,jump,mx,lot,isgn)
        call rfftmlt(c31p,wrk,trigx32,ixfax32,inc,jump,mx,lot,isgn)
        call rfftmlt(c32p,wrk,trigx32,ixfax32,inc,jump,mx,lot,isgn)
        call rfftmlt(c33p,wrk,trigx32,ixfax32,inc,jump,mx,lot,isgn)
     
        call rfftmlt(dc111p,wrk,trigx32,ixfax32,inc,jump,mx,lot,isgn)
        call rfftmlt(dc112p,wrk,trigx32,ixfax32,inc,jump,mx,lot,isgn)
        call rfftmlt(dc113p,wrk,trigx32,ixfax32,inc,jump,mx,lot,isgn)
        call rfftmlt(dc211p,wrk,trigx32,ixfax32,inc,jump,mx,lot,isgn)
        call rfftmlt(dc212p,wrk,trigx32,ixfax32,inc,jump,mx,lot,isgn)
        call rfftmlt(dc213p,wrk,trigx32,ixfax32,inc,jump,mx,lot,isgn)
        call rfftmlt(dc311p,wrk,trigx32,ixfax32,inc,jump,mx,lot,isgn)
        call rfftmlt(dc312p,wrk,trigx32,ixfax32,inc,jump,mx,lot,isgn)
        call rfftmlt(dc313p,wrk,trigx32,ixfax32,inc,jump,mx,lot,isgn)
     
        call rfftmlt(dc121p,wrk,trigx32,ixfax32,inc,jump,mx,lot,isgn)
        call rfftmlt(dc122p,wrk,trigx32,ixfax32,inc,jump,mx,lot,isgn)
        call rfftmlt(dc123p,wrk,trigx32,ixfax32,inc,jump,mx,lot,isgn)
        call rfftmlt(dc221p,wrk,trigx32,ixfax32,inc,jump,mx,lot,isgn)
        call rfftmlt(dc222p,wrk,trigx32,ixfax32,inc,jump,mx,lot,isgn)
        call rfftmlt(dc223p,wrk,trigx32,ixfax32,inc,jump,mx,lot,isgn)
        call rfftmlt(dc321p,wrk,trigx32,ixfax32,inc,jump,mx,lot,isgn)
        call rfftmlt(dc322p,wrk,trigx32,ixfax32,inc,jump,mx,lot,isgn)
        call rfftmlt(dc323p,wrk,trigx32,ixfax32,inc,jump,mx,lot,isgn)
     
        call rfftmlt(dc131p,wrk,trigx32,ixfax32,inc,jump,mx,lot,isgn)
        call rfftmlt(dc132p,wrk,trigx32,ixfax32,inc,jump,mx,lot,isgn)
        call rfftmlt(dc133p,wrk,trigx32,ixfax32,inc,jump,mx,lot,isgn)
        call rfftmlt(dc231p,wrk,trigx32,ixfax32,inc,jump,mx,lot,isgn)
        call rfftmlt(dc232p,wrk,trigx32,ixfax32,inc,jump,mx,lot,isgn)
        call rfftmlt(dc233p,wrk,trigx32,ixfax32,inc,jump,mx,lot,isgn)
        call rfftmlt(dc331p,wrk,trigx32,ixfax32,inc,jump,mx,lot,isgn)
        call rfftmlt(dc332p,wrk,trigx32,ixfax32,inc,jump,mx,lot,isgn)
        call rfftmlt(dc333p,wrk,trigx32,ixfax32,inc,jump,mx,lot,isgn)
        end if  
    !---------------------------------------------------------------------!
    !        Compute the cross product of velocity and vorticity          !
    !---------------------------------------------------------------------!
 
        ! CFL Number check
        do j = 1,mz
            do i = 1,mx
                cflcheck = (abs(up(j,i))/delxm + abs(vp(j,i))/seght(k) + abs(wp(j,i))/delzm)*dt
                if(cflcheck .gt. cfl(k)) then 
                    cfl(k) = cflcheck
                end if
            end do
        end do

        do j = 1,mz
            do i = 1,mx
                vwx(j,i) =  vp(j,i)*wz(j,i) - wp(j,i)*wy(j,i)
                vwy(j,i) = -up(j,i)*wz(j,i) + wp(j,i)*wx(j,i)
                vwz(j,i) =  up(j,i)*wy(j,i) - vp(j,i)*wx(j,i)

                 vc(j,i) = -(up(j,i)*cx(j,i) + vp(j,i)*cy(j,i) + wp(j,i)*cz(j,i))
    
                if (it .ge. src_start-1) then
                c11np(j,i)=c11p(j,i)*u11p(j,i)+c12p(j,i)*u12p(j,i)+c13p(j,i)*u13p(j,i)+  &
                           c11p(j,i)*u11p(j,i)+c21p(j,i)*u12p(j,i)+c31p(j,i)*u13p(j,i)-  &
                           (up(j,i)*dc111p(j,i)+vp(j,i)*dc112p(j,i)+wp(j,i)*dc113p(j,i))            
                  
                c12np(j,i)=c11p(j,i)*u21p(j,i)+c12p(j,i)*u22p(j,i)+c13p(j,i)*u23p(j,i)+   &
                           c12p(j,i)*u11p(j,i)+c22p(j,i)*u12p(j,i)+c32p(j,i)*u13p(j,i)-   &
                           (up(j,i)*dc121p(j,i)+vp(j,i)*dc122p(j,i)+wp(j,i)*dc123p(j,i))
                  
                c13np(j,i)=c11p(j,i)*u31p(j,i)+c12p(j,i)*u32p(j,i)+c13p(j,i)*u33p(j,i)+   &
                           c13p(j,i)*u11p(j,i)+c23p(j,i)*u12p(j,i)+c33p(j,i)*u13p(j,i)-   &
                           (up(j,i)*dc131p(j,i)+vp(j,i)*dc132p(j,i)+wp(j,i)*dc133p(j,i))
                  
                c22np(j,i)=c21p(j,i)*u21p(j,i)+c22p(j,i)*u22p(j,i)+c23p(j,i)*u23p(j,i)+   &
                           c12p(j,i)*u21p(j,i)+c22p(j,i)*u22p(j,i)+c32p(j,i)*u23p(j,i)-   &
                           (up(j,i)*dc221p(j,i)+vp(j,i)*dc222p(j,i)+wp(j,i)*dc223p(j,i))
                  
                c23np(j,i)=c21p(j,i)*u31p(j,i)+c22p(j,i)*u32p(j,i)+c23p(j,i)*u33p(j,i)+   &
                           c13p(j,i)*u21p(j,i)+c23p(j,i)*u22p(j,i)+c33p(j,i)*u23p(j,i)-   &
                           (up(j,i)*dc231p(j,i)+vp(j,i)*dc232p(j,i)+wp(j,i)*dc233p(j,i))

                c33np(j,i)=c31p(j,i)*u31p(j,i)+c32p(j,i)*u32p(j,i)+c33p(j,i)*u33p(j,i)+   &         
                           c13p(j,i)*u31p(j,i)+c23p(j,i)*u32p(j,i)+c33p(j,i)*u33p(j,i)-   &
                           (up(j,i)*dc331p(j,i)+vp(j,i)*dc332p(j,i)+wp(j,i)*dc333p(j,i))
     
                trp(j,i) = c11p(j,i) + c22p(j,i) + c33p(j,i)
     
                if (ipeter .eq. 0) then  ! peterlin function is 1.0
    
                    str11np(j,i)= c11p(j,i)
                    str12np(j,i)= c12p(j,i)
                    str13np(j,i)= c13p(j,i)
                    str22np(j,i)= c22p(j,i)
                    str23np(j,i)= c23p(j,i)
                    str33np(j,i)= c33p(j,i)
              
                else
               
                    str11np(j,i)=((zlmax**2 - 3.0)/(zlmax**2 - trp(j,i)))*c11p(j,i)
                    str12np(j,i)=((zlmax**2 - 3.0)/(zlmax**2 - trp(j,i)))*c12p(j,i)
                    str13np(j,i)=((zlmax**2 - 3.0)/(zlmax**2 - trp(j,i)))*c13p(j,i)
                    str22np(j,i)=((zlmax**2 - 3.0)/(zlmax**2 - trp(j,i)))*c22p(j,i)
                    str23np(j,i)=((zlmax**2 - 3.0)/(zlmax**2 - trp(j,i)))*c23p(j,i)
                    str33np(j,i)=((zlmax**2 - 3.0)/(zlmax**2 - trp(j,i)))*c33p(j,i)
               
                end if
               
                ! Add brownian motion terms
                str11np(j,i)=str11np(j,i)-1.0
                str22np(j,i)=str22np(j,i)-1.0
                str33np(j,i)=str33np(j,i)-1.0

                ! Polymer model
                if (itarget .eq. 0) then ! polymer ocean
                    beta_poly(j,i) = qbeta
                else if (itarget .eq. 1) then
                    ! Nonlinear model:
                    ! beta = exp(-alpha*gamma), alpha = 2.6e-03 PPM --> alpha_poly
                    ! gamma - scalar concentration (PPM) --> scp

                    beta_poly(j,i) = exp(-alpha_poly*abs(scp(j,i)))
                else if (itarget .eq. 2) then ! target swirl
                    ! Use a linear relationship for beta based on swirl
                    ! do nothing for now
                    beta_poly(j,i) = 1.0
                end if

                zbeta1 = (1.0 - beta_poly(j,i))/(re*beta_poly(j,i)*tpoly) ! = (nu_0 - nu_s)/lambda

    !---------------------------------------------------------------------!
!                ! Calculate terms for enstrophy dissipation
!                coeff(k,j,i) = wx(j,i)*zbeta1
!                
!
!                ! Save as 3D variables for writing
!                zbeta3d(k,j,i) = zbeta1*tpoly
!
!                ! Cij tensor physical terms
!                 trp3d(k,j,i) = trp(j,i)
!                c11p3d(k,j,i) = c11p(j,i)
!                c12p3d(k,j,i) = c12p(j,i)
!                c22p3d(k,j,i) = c22p(j,i)
!                
!                ! d(Cij) tensor physical terms
!                dc111p3d(k,j,i) = dc111p(j,i)
!                dc112p3d(k,j,i) = dc112p(j,i)
!                dc121p3d(k,j,i) = dc121p(j,i)
!                dc122p3d(k,j,i) = dc122p(j,i)
!                dc221p3d(k,j,i) = dc221p(j,i)
!                dc222p3d(k,j,i) = dc222p(j,i)
!                
    !---------------------------------------------------------------------!

                ! Polymer stress
                qp11np(j,i) = zbeta1*str11np(j,i)
                qp12np(j,i) = zbeta1*str12np(j,i)
                qp13np(j,i) = zbeta1*str13np(j,i)
                qp22np(j,i) = zbeta1*str22np(j,i)
                qp23np(j,i) = zbeta1*str23np(j,i)
                qp33np(j,i) = zbeta1*str33np(j,i)

                else
                beta_poly(j,i) = 1.0 ! for printing
                end if ! it >= src_start-1

            end do
        end do

   
    !---------------------------------------------------------------------!
    ! Compute immersed boundary force terms and add to nonlinear term     !
    ! The type of forcing is based on the value of imatrix at the grid    !
    ! point. Several of these cases probably require more coding and      !
    ! specification in the geometry file, but they have been unused since !
    ! I've used the code                                                  !
    !---------------------------------------------------------------------!
    
        ! Initialize the force field array
        dragx = 0.0
        dragy = 0.0
        dragz = 0.0
   
        do j = 1,mz
            do i = 1,mx
                if (imatrix(k,j,i) .eq. 1 .or. (imatrix(k,j,i) .eq. 4 .and. it .le. perturbtime)) then ! Generate solid surfaces
                    fxintg(k,j,i) = fxintg(k,j,i) + up(j,i)
                    fyintg(k,j,i) = fyintg(k,j,i) + vp(j,i)
                    fzintg(k,j,i) = fzintg(k,j,i) + wp(j,i)
                    dragx(j,i) = (-ugain*up(j,i)) - gain*fxintg(k,j,i)
                    dragy(j,i) = (-ugain*vp(j,i)) - gain*fyintg(k,j,i)
                    dragz(j,i) = (-ugain*wp(j,i)) - gain*fzintg(k,j,i)
        
                else if (imatrix(k,j,i) .eq. 3) then ! Spanwise-damping textures
                    fzintg(k,j,i) = fzintg(k,j,i) + wp(j,i)
                    dragz(j,i) = -ugain*wp(j,i) - gain*fzintg(k,j,i)
        
                else if (imatrix(k,j,i) .eq. 8) then ! Spanwise moving wall
                    if (i .le. (bfwidth*3 + 1)) then
                        wdes = (0.3*Uinf)/slopelength*(i - (bfwidth*3 + 1 - slopelength))
                    else
                        wdes = 0.3*Uinf
                    end if
    
                    fxintg(k,j,i) = fxintg(k,j,i) + up(j,i)
                    fyintg(k,j,i) = fyintg(k,j,i) + vp(j,i)
                    fzintg(k,j,i) = fzintg(k,j,i) + wp(j,i)-wdes
                    dragx(j,i)    = -ugain*up(j,i) - gain*fxintg(k,j,i)
                    dragy(j,i)    = -ugain*vp(j,i) - gain*fyintg(k,j,i)
                    dragz(j,i)    = -ugain*(wp(j,i)-wdes) - gain*fzintg(k,j,i)
                    
                else if (imatrix(k,j,i) .eq. 6) then ! Buffer zone
                    ! NOTE: initu and initw should technically be transformed first, but since
                    ! the buffer region has only been used for x- and z-constant flows, using
                    ! them as they are should be fine
                    fxintg(k,j,i) = fxintg(k,j,i) + up(j,i) - initu(k,j,i)
                    fyintg(k,j,i) = fyintg(k,j,i) + vp(j,i) - initv(k,j,i)
                    fzintg(k,j,i) = fzintg(k,j,i) + wp(j,i) - initw(k,j,i)
                    vwx(j,i) = vwx(j,i) + (-bfugain(i-bfhead+1)*(up(j,i) - initu(k,j,i)) - (bfgain(i-bfhead+1) * fxintg(k,j,i)))
                    vwy(j,i) = vwy(j,i) + (-bfugain(i-bfhead+1)*(vp(j,i) - initv(k,j,i)) - (bfgain(i-bfhead+1) * fyintg(k,j,i)))
                    vwz(j,i) = vwz(j,i) + (-bfugain(i-bfhead+1)*(wp(j,i) - initw(k,j,i)) - (bfgain(i-bfhead+1) * fzintg(k,j,i)))
    
                else if (imatrix(k,j,i) .eq. 7) then ! Suction region
                    fyintg(k,j,i) = fyintg(k,j,i) + (vp(j,i)-vdes(i))
                    vwy(j,i) = vwy(j,i) + (-gain*fyintg(k,j,i)) + (-ugain*(vp(j,i)-vdes(i)))
                    fxintg(k,j,i) = fxintg(k,j,i) + (up(j,i)-uinf)
                    vwx(j,i) = vwx(j,i) + (-gain*fxintg(k,j,i)) + (-ugain*(up(j,i)-uinf))
                    fzintg(k,j,i) = fzintg(k,j,i) + (wp(j,i)-initw(k,j,i))
                    vwz(j,i) = vwz(j,i) + (-gain*fzintg(k,j,i)) + (-ugain*(wp(j,i)-initw(k,j,i)))
     
                else if (imatrix(k,j,i) .eq. 11) then ! Constant x-force between wall and suction layer
                    dragx(j,i) = 1.0
    
                else if (imatrix(k,j,i) .eq. 0 .and. flow_select .eq. 1) then ! Constant pressure gradient
                    vwx(j,i) = vwx(j,i) - dPdx 
                
                else if (flow_select .eq. 5 .and. it .le. perturbtime) then ! Vortex ring
                    zj = float(j-1)*delzm
                    rsq = (ycoord(k) - ycenter)**2 + (zj - zcenter)**2
                    argrad = forbeta*(rad - sqrt(rsq))
                    fr = 0.5*(1.0 + tanh(argrad))
    
                    xi = float(i-1)*delxm
                    argx = forbeta*(L - abs(xi - xcenter))
                    fx = 0.5*(1.0 + tanh(argx))
  
                    vwx(j,i) = vwx(j,i) + bdyfx*fx*fr
    
                end if

                if (scl_flag .ge. 2) then
                    vc(j,i) = vc(j,i) + scsource(k,j,i)/dt
                end if
            end do ! i
        end do ! j
    
        !---------------------------------------------------------------------!
        !            Apply the force field to the solid surface               !
        !---------------------------------------------------------------------!
        if(k .le. kmaxsurf) then
            do i = 1,mx
                do ii = -2, 2
                    ipii = i + ii
                    idif = 1 + iabs(ii)        
                    if(ipii .lt. 1 ) ipii = ipii + mx         
                    if(ipii .gt. mx) ipii = ipii - mx
                    do jj = -2, 2
                        jdif = 1 + iabs(jj)
                        segdrag = fspread(jdif,idif)
    
                        do j = 3, mzm2
                            jpjj = j + jj
                            xsegdrag = segdrag*dragx(j,i)
                            ysegdrag = segdrag*dragy(j,i)
                            zsegdrag = segdrag*dragz(j,i)
                            vwx(jpjj,ipii) = vwx(jpjj,ipii) + xsegdrag
                            vwy(jpjj,ipii) = vwy(jpjj,ipii) + ysegdrag
                            vwz(jpjj,ipii) = vwz(jpjj,ipii) + zsegdrag
                        end do
    
                        do j = 1, 2
                            jpjj = j + jj
                            if(jpjj .lt. 1 ) jpjj = jpjj + mz
                            xsegdrag = segdrag*dragx(j,i)
                            ysegdrag = segdrag*dragy(j,i)
                            zsegdrag = segdrag*dragz(j,i)
                            vwx(jpjj,ipii) = vwx(jpjj,ipii) + xsegdrag
                            vwy(jpjj,ipii) = vwy(jpjj,ipii) + ysegdrag
                            vwz(jpjj,ipii) = vwz(jpjj,ipii) + zsegdrag
                        end do
    
                        do j = mzm1, mz
                            jpjj = j + jj
                            if(jpjj .gt. mz) jpjj = jpjj - mz
                            xsegdrag = segdrag*dragx(j,i)
                            ysegdrag = segdrag*dragy(j,i)
                            zsegdrag = segdrag*dragz(j,i)
                            vwx(jpjj,ipii) = vwx(jpjj,ipii) + xsegdrag
                            vwy(jpjj,ipii) = vwy(jpjj,ipii) + ysegdrag
                            vwz(jpjj,ipii) = vwz(jpjj,ipii) + zsegdrag
                        end do
    
                    end do ! jj
                end do ! ii
            end do ! i
        end if
    
        !---------------------------------------------------------------------!
        ! Save variables to print (or use in particle integration)            !
        ! Note: There's a weird artifact of how the variables are transformed !
        ! which flips the y-direction. In order to print out correctly, all   !
        ! terms which include v or an odd number of derivatives in y must be  !
        ! made negative. I know, I don't really get it either but it works.   !
        !---------------------------------------------------------------------!
    
        up3d(k,1:mz,1:mx) = up(1:mz,1:mx)
        vp3d(k,1:mz,1:mx) = -vp(1:mz,1:mx)
        wp3d(k,1:mz,1:mx) = wp(1:mz,1:mx)
        wx3d(k,1:mz,1:mx) = -wx(1:mz,1:mx)
        wy3d(k,1:mz,1:mx) = wy(1:mz,1:mx)
        wz3d(k,1:mz,1:mx) = -wz(1:mz,1:mx)
        u11p3d(k,1:mz,1:mx) = u11p(1:mz,1:mx) 
        u12p3d(k,1:mz,1:mx) = -u12p(1:mz,1:mx) 
        u13p3d(k,1:mz,1:mx) = u13p(1:mz,1:mx) 
        u21p3d(k,1:mz,1:mx) = -u21p(1:mz,1:mx) 
        u22p3d(k,1:mz,1:mx) = u22p(1:mz,1:mx) 
        u23p3d(k,1:mz,1:mx) = -u23p(1:mz,1:mx) 
        u31p3d(k,1:mz,1:mx) = u31p(1:mz,1:mx) 
        u32p3d(k,1:mz,1:mx) = -u32p(1:mz,1:mx) 
        u33p3d(k,1:mz,1:mx) = u33p(1:mz,1:mx) 
    
        Lup3d(k,1:mz,1:mx) = Lup(1:mz,1:mx) 
        Lvp3d(k,1:mz,1:mx) = Lvp(1:mz,1:mx) 
        Lwp3d(k,1:mz,1:mx) = Lwp(1:mz,1:mx) 
   
        scp3d(k,1:mz,1:mx) = scp(1:mz,1:mx)

        beta3d(k,1:mz,1:mx) = beta_poly(1:mz,1:mx)
 
        !---------------------------------------------------------------------!
        !              Now transform vxw back to spectral space               !
        !---------------------------------------------------------------------!
    
        isgn = -1
        inc  = mzp
        jump = 1
        lot  = mz
    
        call rfftmlt(vwx,wrk,trigx32,ixfax32,inc,jump,mx,lot,isgn)
        call rfftmlt(vwy,wrk,trigx32,ixfax32,inc,jump,mx,lot,isgn)
        call rfftmlt(vwz,wrk,trigx32,ixfax32,inc,jump,mx,lot,isgn)
        call rfftmlt(vc ,wrk,trigx32,ixfax32,inc,jump,mx,lot,isgn) 

        if (it .ge. src_start-1) then
        call rfftmlt(c11np,wrk,trigx32,ixfax32,inc,jump,mx,lot,isgn)
        call rfftmlt(c12np,wrk,trigx32,ixfax32,inc,jump,mx,lot,isgn)
        call rfftmlt(c13np,wrk,trigx32,ixfax32,inc,jump,mx,lot,isgn)
        call rfftmlt(c22np,wrk,trigx32,ixfax32,inc,jump,mx,lot,isgn)
        call rfftmlt(c23np,wrk,trigx32,ixfax32,inc,jump,mx,lot,isgn)
        call rfftmlt(c33np,wrk,trigx32,ixfax32,inc,jump,mx,lot,isgn)
     
        call rfftmlt(str11np,wrk,trigx32,ixfax32,inc,jump,mx,lot,isgn)
        call rfftmlt(str12np,wrk,trigx32,ixfax32,inc,jump,mx,lot,isgn)
        call rfftmlt(str13np,wrk,trigx32,ixfax32,inc,jump,mx,lot,isgn)
        call rfftmlt(str22np,wrk,trigx32,ixfax32,inc,jump,mx,lot,isgn)
        call rfftmlt(str23np,wrk,trigx32,ixfax32,inc,jump,mx,lot,isgn)
        call rfftmlt(str33np,wrk,trigx32,ixfax32,inc,jump,mx,lot,isgn)
     
        call rfftmlt(qp11np,wrk,trigx32,ixfax32,inc,jump,mx,lot,isgn)
        call rfftmlt(qp12np,wrk,trigx32,ixfax32,inc,jump,mx,lot,isgn)
        call rfftmlt(qp13np,wrk,trigx32,ixfax32,inc,jump,mx,lot,isgn)
        call rfftmlt(qp22np,wrk,trigx32,ixfax32,inc,jump,mx,lot,isgn)
        call rfftmlt(qp23np,wrk,trigx32,ixfax32,inc,jump,mx,lot,isgn)
        call rfftmlt(qp33np,wrk,trigx32,ixfax32,inc,jump,mx,lot,isgn)

        do j = 1,mz
            c11np(j,2) = 0.0  
            c12np(j,2) = 0.0  
            c13np(j,2) = 0.0  
            c22np(j,2) = 0.0  
            c23np(j,2) = 0.0  
            c33np(j,2) = 0.0  
         
            str11np(j,2) = 0.0 
            str12np(j,2) = 0.0 
            str13np(j,2) = 0.0 
            str22np(j,2) = 0.0 
            str23np(j,2) = 0.0 
            str33np(j,2) = 0.0 
         
            qp11np(j,2) = 0.0 
            qp12np(j,2) = 0.0 
            qp13np(j,2) = 0.0 
            qp22np(j,2) = 0.0 
            qp23np(j,2) = 0.0 
            qp33np(j,2) = 0.0 
        end do
        end if 
        
        ! Wipe out unwanted last mode
        do j = 1,mz
            vwx(j,2) = 0.0
            vwy(j,2) = 0.0
            vwz(j,2) = 0.0
             vc(j,2) = 0.0
        end do
        
        isgn = -1
        inc  = 1
        jump = 2*mzp
        lot  = nx/2
    
        call cfftmlt(vwx(1,1),vwx(1,2),wrk,trigz32,izfax32,inc,jump,mz,lot,isgn)
        call cfftmlt(vwy(1,1),vwy(1,2),wrk,trigz32,izfax32,inc,jump,mz,lot,isgn)
        call cfftmlt(vwz(1,1),vwz(1,2),wrk,trigz32,izfax32,inc,jump,mz,lot,isgn)
        call cfftmlt(vc(1,1) ,vc(1,2) ,wrk,trigz32,izfax32,inc,jump,mz,lot,isgn)
        if (it .ge. src_start-1) then
        call cfftmlt(c11np(1,1),c11np(1,2),wrk,trigz32,izfax32,inc,jump,mz,lot,isgn)
        call cfftmlt(c12np(1,1),c12np(1,2),wrk,trigz32,izfax32,inc,jump,mz,lot,isgn)
        call cfftmlt(c13np(1,1),c13np(1,2),wrk,trigz32,izfax32,inc,jump,mz,lot,isgn)
        call cfftmlt(c22np(1,1),c22np(1,2),wrk,trigz32,izfax32,inc,jump,mz,lot,isgn)
        call cfftmlt(c23np(1,1),c23np(1,2),wrk,trigz32,izfax32,inc,jump,mz,lot,isgn)
        call cfftmlt(c33np(1,1),c33np(1,2),wrk,trigz32,izfax32,inc,jump,mz,lot,isgn)
     
        call cfftmlt(str11np(1,1),str11np(1,2),wrk,trigz32,izfax32,inc,jump,mz,lot,isgn)
        call cfftmlt(str12np(1,1),str12np(1,2),wrk,trigz32,izfax32,inc,jump,mz,lot,isgn)
        call cfftmlt(str13np(1,1),str13np(1,2),wrk,trigz32,izfax32,inc,jump,mz,lot,isgn)
        call cfftmlt(str22np(1,1),str22np(1,2),wrk,trigz32,izfax32,inc,jump,mz,lot,isgn)
        call cfftmlt(str23np(1,1),str23np(1,2),wrk,trigz32,izfax32,inc,jump,mz,lot,isgn)
        call cfftmlt(str33np(1,1),str33np(1,2),wrk,trigz32,izfax32,inc,jump,mz,lot,isgn)
     
        call cfftmlt(qp11np(1,1),qp11np(1,2),wrk,trigz32,izfax32,inc,jump,mz,lot,isgn)
        call cfftmlt(qp12np(1,1),qp12np(1,2),wrk,trigz32,izfax32,inc,jump,mz,lot,isgn)
        call cfftmlt(qp13np(1,1),qp13np(1,2),wrk,trigz32,izfax32,inc,jump,mz,lot,isgn)
        call cfftmlt(qp22np(1,1),qp22np(1,2),wrk,trigz32,izfax32,inc,jump,mz,lot,isgn)
        call cfftmlt(qp23np(1,1),qp23np(1,2),wrk,trigz32,izfax32,inc,jump,mz,lot,isgn)
        call cfftmlt(qp33np(1,1),qp33np(1,2),wrk,trigz32,izfax32,inc,jump,mz,lot,isgn)
        end if 
    
        do j = 1,nz
            if (j .le. nzh) then
                jj = j
            else if (j .gt. nzh) then
                jj = (mz-nz) + j
            end if
            do i = 1,nxh
                i1 = 2*(i-1) + 1
                i2 = 2*i
                gn(k,j,i)  = cmplx(vwx(jj,i1)*rmz,vwx(jj,i2)*rmz)
                fn(k,j,i)  = cmplx(vwy(jj,i1)*rmz,vwy(jj,i2)*rmz)
                omz(k,j,i) = cmplx(vwz(jj,i1)*rmz,vwz(jj,i2)*rmz)
                scn(k,j,i) = cmplx( vc(jj,i1)*rmz, vc(jj,i2)*rmz)

                if (ipolyflag .eq. 1 .and. it .ge. src_start-1) then
                c11n(k,j,i) = cmplx(c11np(jj,i1)*rmz,c11np(jj,i2)*rmz)
                c12n(k,j,i) = cmplx(c12np(jj,i1)*rmz,c12np(jj,i2)*rmz)
                c13n(k,j,i) = cmplx(c13np(jj,i1)*rmz,c13np(jj,i2)*rmz)
                c22n(k,j,i) = cmplx(c22np(jj,i1)*rmz,c22np(jj,i2)*rmz)
                c23n(k,j,i) = cmplx(c23np(jj,i1)*rmz,c23np(jj,i2)*rmz)
                c33n(k,j,i) = cmplx(c33np(jj,i1)*rmz,c33np(jj,i2)*rmz)
     
                str11n(k,j,i) = cmplx(str11np(jj,i1)*rmz,str11np(jj,i2)*rmz)
                str12n(k,j,i) = cmplx(str12np(jj,i1)*rmz,str12np(jj,i2)*rmz)
                str13n(k,j,i) = cmplx(str13np(jj,i1)*rmz,str13np(jj,i2)*rmz)
                str22n(k,j,i) = cmplx(str22np(jj,i1)*rmz,str22np(jj,i2)*rmz)
                str23n(k,j,i) = cmplx(str23np(jj,i1)*rmz,str23np(jj,i2)*rmz)
                str33n(k,j,i) = cmplx(str33np(jj,i1)*rmz,str33np(jj,i2)*rmz)
     
                qp11(k,j,i) = cmplx(qp11np(jj,i1)*rmz,qp11np(jj,i2)*rmz)
                qp12(k,j,i) = cmplx(qp12np(jj,i1)*rmz,qp12np(jj,i2)*rmz)
                qp13(k,j,i) = cmplx(qp13np(jj,i1)*rmz,qp13np(jj,i2)*rmz)
                qp22(k,j,i) = cmplx(qp22np(jj,i1)*rmz,qp22np(jj,i2)*rmz)
                qp23(k,j,i) = cmplx(qp23np(jj,i1)*rmz,qp23np(jj,i2)*rmz)
                qp33(k,j,i) = cmplx(qp33np(jj,i1)*rmz,qp33np(jj,i2)*rmz)
                end if
            end do
        end do
   
         
    end do ! k
    !$omp end parallel do
    !                -----  end of loop over normal(y) planes  -----
   

    !---------------------------------------------------------------------!
    !     Calculate swirl criterion and write data for visualization      !
    !---------------------------------------------------------------------!
    
    ! Calculate swirl each time step to calculate area of vortex core 
    !$omp parallel do private(i,j,k,swirl)
    do i = 1,nyp
        do j = 1,mz
            do k = 1,mx
                call calcswirl(u11p3d(i,j,k),u21p3d(i,j,k),u31p3d(i,j,k),u12p3d(i,j,k),u22p3d(i,j,k), &
                               u32p3d(i,j,k),u13p3d(i,j,k),u23p3d(i,j,k),u33p3d(i,j,k),swirl)
    
                swirl_3d(i,j,k) = swirl
            end do
        end do
    end do
    !$omp end parallel do
    
    call vortArea(swirl_3d,vArea) ! Calculate area of vortex based on swirl criterion

    if (print3d .ne. 0) then
    
        if ((mod(it,iprnfrq) .eq. 0 .and. it .ne. 0) .or. it .eq. 1) then
            
            ! Write output files
            if (print3d .eq. 1) then
                call write_flowfield_ascii(up3d,vp3d,wp3d,wx3d,wy3d,wz3d,swirl_3d,real(imatrix))
            else if (print3d .eq. 3) then
                call write_flowfield_plt(up3d,vp3d,wp3d,wx3d,wy3d,wz3d,swirl_3d,real(imatrix),beta3d,rank)
                ! In certain cases, it may be better to print both scalar and beta, but since that's
                ! not relevant at the moment, I'm just leaving beta
            else
                write(*,*) 'Warning: Unknown print type. No output data will be written.'
            end if
        
            ! Print out wx and scalar as a function of r
            write(filename2,'("outputs",i2.2,"/morestuff/fr-",i6.6,".dat")') rank,it
            open(190, file = filename2)
            do j = 1,mz/2-1
                z = delzm*(j-1)
                r = abs(z - vortZ)
                wxi = wx3d(nyh,mz/2+j,1) + wx3d(nyh,mz/2-j,1)
                do k = 1,3
                    argy = pi/4.0*float(k)
                    z1 = r*cos(argy) + vortZ
                    y1 = r*sin(argy) + vortY
                    y2 = vortY - r*sin(argy)
                   
                    if (y1 .gt. 0.0 .and. y1 .lt. yl) then 
                        call fluid_interp(0.0,y1,z1,wx3d,scp3d,up3d,wxj,scj,up1)
                        wxi = wxi + wxj
                        sci = sci + scj
                        call fluid_interp(0.0,y2,z1,wx3d,scp3d,up3d,wxj,scj,up1)
                        wxi = wxi + wxj
                        sci = sci + scj
                    end if
                end do

                wxi = wxi/8.0
                sci = sci/8.0

                write(190,"(3(e14.6,1x))") r, wxi ,sci
            end do
            close(190)

        end if
    end if
   
    !---------------------------------------------------------------------!
    !                   Particle Tracking Implementation                  !
    !---------------------------------------------------------------------!
    
    if (npart .ne. 0) then
        call part_track(up3d,vp3d,wp3d,wx3d,wy3d,wz3d,u_old,v_old,w_old,  &
                        u11p3d,u12p3d,u13p3d,u21p3d,u22p3d,u23p3d,u31p3d, &
                        u32p3d,u33p3d,Lup3d,Lvp3d,Lwp3d,Lu_old,Lv_old,Lw_old,vArea,rank)
    
        ! Save velocity and Laplacian for time derivatives inside particle integration
        u_old = up3d
        v_old = vp3d
        w_old = wp3d
    
        Lu_old = Lup3d
        Lv_old = Lvp3d
        Lw_old = Lwp3d
    end if
    
    !---------------------------------------------------------------------!
    ! Write Mean U Data, calculate mass flux, & KE/enstrophy if relevant  !
    !---------------------------------------------------------------------!
  
    ! Make these into subroutines later...  
    if (flow_select .eq. 1 .or. flow_select .eq. 4) then ! Only relevant for wall-bounded turbulence
        write(*,*) 'Writing mean U data...'
    
!        ! Mean velocity
!        if (irstrt .eq. it) then
!            open(71,file = 'outputs/mean_u_data.dat')
!        else
!            open(71,file = 'outputs/mean_u_data.dat', position = 'append')
!        end if
    
        do i = 1,nyp
            do j = 1,mz
                uxmean(j) = sum(up3d(i,j,:))/mx
            end do
            uzmean(i) = sum(uxmean)/mz
!            write(71,"(*(e14.6,1x))") uzmean
        end do
    
!        close(71)
!        write(*,*) '    Done!'
    end if
    
!    ! Compute mass flux: Integrate <U> over y
!    if (irstrt .eq. it) then
!        open(72,file='outputs/mass_flux')
!    else
!        open(72,file='outputs/mass_flux',position='append')
!    end if
    
!    massFlux = 0.0
!    do i = 1,ny
!        massFlux = massFlux + 1.0/yl*0.5*(uzmean(i+1) + uzmean(i))*(ycoord(i+1) - ycoord(i)) ! Bulk velocity
!    end do
!    
!    write(72,*) massFlux
!    close(72)
    
    ! Compute average KE and Enstrophy for Vortex Ring
    if (mod(flow_select,5) .eq. 0) then
        sumens = 0.0
        KE = 0.0
        scl_total = 0.0
   
        Lx = delxm; Lz = delzm 
        do k = 1,mx
            do j = 1,mz
                do i = 1,nyp
                    Ly = seght(i)
                    sumens = sumens + (wx3d(i,j,k)**2 + wy3d(i,j,k)**2 + wz3d(i,j,k)**2)*Lx*Ly*Lz
                    KE = KE + 0.5*(up3d(i,j,k)**2 + vp3d(i,j,k)**2 + wp3d(i,j,k)**2)*Lx*Ly*Lz
                    scl_total = scl_total + scp3d(i,j,k)*Lx*Ly*Lz
                    SG = SG + swirl_3d(i,j,k)*scp3d(i,j,k)*Lx*Ly*Lz
                end do
            end do
        end do

        ! Fix 2D bug   
        if (flow_select .ge. 10) then
            sumens = sumens/xl
            KE = KE/xl
            scl_total = scl_total/xl
            SG = SG/xl
        end if
 
        write(dirname,'("outputs",i2.2,"/")') rank
        if (it .eq. irstrt) then
            open(73,file=dirname//'enstrophy')
            open(74,file=dirname//'KE')
            open(75,file=dirname//'total_scalar')
            open(76,file=dirname//'SG')
        else
            open(73,file=dirname//'enstrophy',position='append')
            open(74,file=dirname//'KE',position='append')
            open(75,file=dirname//'total_scalar',position='append')
            open(76,file=dirname//'SG',position='append')
        end if
    
        write(73,*) sumens
        write(74,*) KE
        write(75,*) scl_total
        write(76,*) SG
        close(73)
        close(74)
        close(75) 
        close(76) 
    end if

    ! Check total amount of polymer added
!    if (scl_flag .eq. 2 .and. it .le. src_stop .and. it .ge. src_start) then
        call calc_total_beta(it,delxm,delzm,scp3d,beta3d)
!    end if

!    ! Calculate enstrophy equation terms
!    if (ipeter .eq. 1) then
!        ! First calculate necessary derivatives 
!        call calc_f_derivatives(trp3d,df1,df2,df11,df12,df22,    &
!                                zbeta3d,dnu1,dnu2,               &
!                                dc121p3d,dc1211,dc221p3d,dc2212, &
!                                dc111p3d,dc1112,dc122p3d,dc1222)
!
!        ! Calculate enstrophy terms
!!        call calc_ens_terms(coeff,df1,df2,df11,df12,df22, &
!!                            trp3d,c11p3d,c12p3d,c22p3d,   &
!!                            dc111p3d,dc112p3d,dc121p3d,   &
!!                            dc122p3d,dc221p3d,dc222p3d,   &
!!                            dc1112,dc1222,dc1211,dc2212)
!    end if

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


