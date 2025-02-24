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
    
    !omp parallel do shared(wrk,c,s) default(private) schedule(dynamic) 
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
    !omp end parallel do
    
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
    !$omp parallel do schedule(dynamic) shared(bctop,bcbot) private(j,k)
    do k = 1,nxh
        do j = 1,nz
            bctop(j,k) = (0.0,0.0)
            bcbot(j,k) = (0.0,0.0)
        end do
    end do
    !$omp end parallel do 
    
    !$omp parallel do default(shared) private(i,j,k) collapse(3)
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
    
    !$omp parallel do default(shared) private(i,j,k,w2) collapse(2) schedule(dynamic)
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
    
    !$omp parallel do default(shared) private(i,j,k) schedule(dynamic)
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
    use omp_lib
    
    ! Delcare variables
    implicit none
    
    complex, dimension(nyp,nz) :: wrk
    complex :: s(0:ny,nz) ! note zero-based arrays
    real :: c(0:ny)
    
    integer :: i,j,k,jp2,jm2
    real    :: rj1,rj2,rj3
   
    !omp parallel do shared(wrk,c,s) default(private) schedule(dynamic) 
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
    !omp end parallel do
    
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
    use omp_lib

    implicit none
    
    complex, dimension(nyp,nz) :: wrk
    complex :: s(0:ny,nz) ! note zero-based arrays
    
    integer :: i,j,k,jp2,jm2
    real    :: rj1,rj2,rj3

    !omp parallel do shared(wrk,s) default(private) schedule(dynamic)    
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
    !omp end parallel do
    
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
    use omp_lib
    
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
   
    !omp parallel do default(shared) private(i,k) schedule(dynamic) 
    do i = 1,m
        do k = 1,nz
            gs(i,k,ip) = wrk(i,k)
            cs(i,k,ip) = a(i,k,1)
        end do
    end do
    !omp end parallel do
    
    do k = 1,nz
        ds(m,k,ip) = a(m,k,2)
        gs(n,k,ip) = wrk(n,k)
    end do
   
    !omp parallel do default(shared) private(i,k) collapse(2) 
    do i = mm1,1,-1
        do k = 1,nz
            ds(i,k,ip) = a(i,k,2) - a(i,k,3)*cs(i+1,k,ip)/ds(i+1,k,ip)
            gs(i+1,k,ip) = gs(i+1,k,ip) - a(i,k,3)*gs(i+2,k,ip)/ds(i+1,k,ip)
        end do
    end do
    !omp end parallel do
    
    ! Eliminate 2nd row on ipass = 1 and 1st row on ipass = 2
    
2 continue
    
    call fillt(n,ip,ipass,t)
    
    do k = 1,nz
        f(js,k,ip) = t(m) - t(n)*cs(m,k,ip)/ds(m,k,ip)
        s(js,k) = s(js,k) - t(n)*gs(n,k,ip)/ds(m,k,ip)
    end do
   
    !omp parallel do shared(s,f,gs,cs,ds,t) private(i,k) collapse(2) 
    do i = mm1,1,-1
        do k = 1,nz
            s(js,k) = s(js,k) - f(js,k,ip)*gs(i+1,k,ip)/ds(i,k,ip)
            f(js,k,ip) = t(i) - f(js,k,ip)*cs(i,k,ip)/ds(i,k,ip)
        end do
    end do
    !omp end parallel do
    
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
    
            ! Even u
            do i = 2,n+1
            ii = 2*i - 1
                wrk(ii,k) = (gs(i,k,1) - cs(i-1,k,1)*wrk(ii-2,k))/ds(i-1,k,1)
            end do
 
            ! Odd u
            do i = 2,n
                ii = 2*i
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
                    u31,u32,u33, &
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
                    Lu,Lv,Lw)

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
 
    ! Passed variables
    type(C_PTR),save :: planZb,planXb,planXZb,planY,planXZf,planZf,planXf

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
    real,    dimension(mz,mx)      :: beta_poly
    real,    dimension(nyp,mz,mx)  :: beta3d,p12,trC
#ENDIF
   
    ! FFT complex arrays 
    complex(C_DOUBLE_COMPLEX), dimension(mz,mx) :: bcomp

    complex(C_DOUBLE_COMPLEX), dimension(mz,mx) :: us,vs,ws,wxs,wys,wzs
    complex(C_DOUBLE_COMPLEX), dimension(mz,mx) :: u11s,u12s,u13s,u21s,u22s,u23s,u31s,u32s,u33s
    complex(C_DOUBLE_COMPLEX), dimension(mz,mx) :: Lus,Lvs,Lws

    complex(C_DOUBLE_COMPLEX), dimension(mz,mx) :: vwxs,vwys,vwzs
#IFDEF SCALAR
    complex(C_DOUBLE_COMPLEX), dimension(mz,mx) :: scs,cxs,cys,czs
    complex(C_DOUBLE_COMPLEX), dimension(mz,mx) :: vcs
#ENDIF
#IFDEF POLYMER
    complex(C_DOUBLE_COMPLEX), dimension(mz,mx) :: c11s,c12s,c13s,c21s,c22s,c23s,c31s,c32s,c33s
    complex(C_DOUBLE_COMPLEX), dimension(mz,mx) :: c11ns,c12ns,c13ns,c22ns,c23ns,c33ns

    complex(C_DOUBLE_COMPLEX), dimension(mz,mx) :: dc111s,dc112s,dc113s,dc121s,dc122s,dc123s,dc131s,dc132s,dc133s
    complex(C_DOUBLE_COMPLEX), dimension(mz,mx) :: dc211s,dc212s,dc213s,dc221s,dc222s,dc223s,dc231s,dc232s,dc233s
    complex(C_DOUBLE_COMPLEX), dimension(mz,mx) :: dc311s,dc312s,dc313s,dc321s,dc322s,dc323s,dc331s,dc332s,dc333s

    complex(C_DOUBLE_COMPLEX), dimension(mz,mx) :: str11ns,str12ns,str13ns,str22ns,str23ns,str33ns

    complex(C_DOUBLE_COMPLEX), dimension(mz,mx) :: qp11s,qp12s,qp13s,qp22s,qp23s,qp33s
#ENDIF
    
    ! Calculation variables
    real, save, dimension(nyp,mz,mx) :: u_old,v_old,w_old
    real, save, dimension(nyp,mz,mx) :: Lunm1,Lvnm1,Lwnm1

    ! FFTW variables - must be C-type arrays
    ! Input variables
    real(C_DOUBLE), dimension(nyp,nz,nxh) :: areal
    real(C_DOUBLE), dimension(mz,mx) :: breal

    real(C_DOUBLE), dimension(mz,mx) :: up,vp,wp,wxp,wyp,wzp
    real(C_DOUBLE), dimension(mz,mx) :: u11p,u12p,u13p
    real(C_DOUBLE), dimension(mz,mx) :: u21p,u22p,u23p
    real(C_DOUBLE), dimension(mz,mx) :: u31p,u32p,u33p
    real(C_DOUBLE), dimension(mz,mx) :: Lup,Lvp,Lwp
    
    real(C_DOUBLE), dimension(mz,mx) :: vwx,vwy,vwz

    real(C_DOUBLE), dimension(nyp,nz,nxh) :: ur,vr,wr,wxr,wyr,wzr
    real(C_DOUBLE), dimension(nyp,nz,nxh) :: u11r,u12r,u13r,u21r,u22r,u23r,u31r,u32r,u33r
    real(C_DOUBLE), dimension(nyp,nz,nxh) :: Lur,Lvr,Lwr

    real(C_DOUBLE), dimension(nyp,nz,nxh) :: fnr,gnr,omzr

    real(C_DOUBLE), dimension(nyp,nz,nxh) :: ui,vi,wi,wxi,wyi,wzi
    real(C_DOUBLE), dimension(nyp,nz,nxh) :: u11i,u12i,u13i,u21i,u22i,u23i,u31i,u32i,u33i
    real(C_DOUBLE), dimension(nyp,nz,nxh) :: Lui,Lvi,Lwi

    real(C_DOUBLE), dimension(nyp,nz,nxh) :: fni,gni,omzi
#IFDEF POLYMER    
    real(C_DOUBLE), dimension(nyp,nz,nxh)  :: c11r,c12r,c13r
    real(C_DOUBLE), dimension(nyp,nz,nxh)  :: c21r,c22r,c23r
    real(C_DOUBLE), dimension(nyp,nz,nxh)  :: c31r,c32r,c33r

    real(C_DOUBLE), dimension(nyp,nz,nxh)  :: dc111r, dc112r, dc113r
    real(C_DOUBLE), dimension(nyp,nz,nxh)  :: dc211r, dc212r, dc213r
    real(C_DOUBLE), dimension(nyp,nz,nxh)  :: dc311r, dc312r, dc313r
    
    real(C_DOUBLE), dimension(nyp,nz,nxh)  :: dc121r, dc122r, dc123r
    real(C_DOUBLE), dimension(nyp,nz,nxh)  :: dc221r, dc222r, dc223r
    real(C_DOUBLE), dimension(nyp,nz,nxh)  :: dc321r, dc322r, dc323r
    
    real(C_DOUBLE), dimension(nyp,nz,nxh)  :: dc131r, dc132r, dc133r
    real(C_DOUBLE), dimension(nyp,nz,nxh)  :: dc231r, dc232r, dc233r
    real(C_DOUBLE), dimension(nyp,nz,nxh)  :: dc331r, dc332r, dc333r

    real(C_DOUBLE), dimension(nyp,nz,nxh)  :: c11nr,c12nr,c13nr
    real(C_DOUBLE), dimension(nyp,nz,nxh)  :: c22nr,c23nr,c33nr
    
    real(C_DOUBLE), dimension(nyp,nz,nxh)  :: str11nr,str12nr,str13nr
    real(C_DOUBLE), dimension(nyp,nz,nxh)  :: str22nr,str23nr,str33nr
    
    real(C_DOUBLE), dimension(nyp,nz,nxh)  :: qp11r,qp12r,qp13r
    real(C_DOUBLE), dimension(nyp,nz,nxh)  :: qp22r,qp23r,qp33r


    real(C_DOUBLE), dimension(nyp,nz,nxh)  :: c11i,c12i,c13i
    real(C_DOUBLE), dimension(nyp,nz,nxh)  :: c21i,c22i,c23i
    real(C_DOUBLE), dimension(nyp,nz,nxh)  :: c31i,c32i,c33i

    real(C_DOUBLE), dimension(nyp,nz,nxh)  :: dc111i, dc112i, dc113i
    real(C_DOUBLE), dimension(nyp,nz,nxh)  :: dc211i, dc212i, dc213i
    real(C_DOUBLE), dimension(nyp,nz,nxh)  :: dc311i, dc312i, dc313i
    
    real(C_DOUBLE), dimension(nyp,nz,nxh)  :: dc121i, dc122i, dc123i
    real(C_DOUBLE), dimension(nyp,nz,nxh)  :: dc221i, dc222i, dc223i
    real(C_DOUBLE), dimension(nyp,nz,nxh)  :: dc321i, dc322i, dc323i
    
    real(C_DOUBLE), dimension(nyp,nz,nxh)  :: dc131i, dc132i, dc133i
    real(C_DOUBLE), dimension(nyp,nz,nxh)  :: dc231i, dc232i, dc233i
    real(C_DOUBLE), dimension(nyp,nz,nxh)  :: dc331i, dc332i, dc333i

    real(C_DOUBLE), dimension(nyp,nz,nxh)  :: c11ni,c12ni,c13ni
    real(C_DOUBLE), dimension(nyp,nz,nxh)  :: c22ni,c23ni,c33ni
    
    real(C_DOUBLE), dimension(nyp,nz,nxh)  :: str11ni,str12ni,str13ni
    real(C_DOUBLE), dimension(nyp,nz,nxh)  :: str22ni,str23ni,str33ni
    
    real(C_DOUBLE), dimension(nyp,nz,nxh)  :: qp11i,qp12i,qp13i
    real(C_DOUBLE), dimension(nyp,nz,nxh)  :: qp22i,qp23i,qp33i


    real(C_DOUBLE), dimension(mz,mx)  :: c11p,c12p,c13p
    real(C_DOUBLE), dimension(mz,mx)  :: c21p,c22p,c23p
    real(C_DOUBLE), dimension(mz,mx)  :: c31p,c32p,c33p

    real(C_DOUBLE), dimension(mz,mx)  :: trp

    real(C_DOUBLE), dimension(mz,mx)  :: dc111p, dc112p, dc113p
    real(C_DOUBLE), dimension(mz,mx)  :: dc211p, dc212p, dc213p
    real(C_DOUBLE), dimension(mz,mx)  :: dc311p, dc312p, dc313p
    
    real(C_DOUBLE), dimension(mz,mx)  :: dc121p, dc122p, dc123p
    real(C_DOUBLE), dimension(mz,mx)  :: dc221p, dc222p, dc223p
    real(C_DOUBLE), dimension(mz,mx)  :: dc321p, dc322p, dc323p
    
    real(C_DOUBLE), dimension(mz,mx)  :: dc131p, dc132p, dc133p
    real(C_DOUBLE), dimension(mz,mx)  :: dc231p, dc232p, dc233p
    real(C_DOUBLE), dimension(mz,mx)  :: dc331p, dc332p, dc333p

    real(C_DOUBLE), dimension(mz,mx)  :: c11np,c12np,c13np
    real(C_DOUBLE), dimension(mz,mx)  :: c22np,c23np,c33np
    
    real(C_DOUBLE), dimension(mz,mx)  :: str11np,str12np,str13np
    real(C_DOUBLE), dimension(mz,mx)  :: str22np,str23np,str33np
    
    real(C_DOUBLE), dimension(mz,mx)  :: qp11np,qp12np,qp13np
    real(C_DOUBLE), dimension(mz,mx)  :: qp22np,qp23np,qp33np
#ENDIF
#IFDEF SCALAR    
    real(C_DOUBLE), dimension(nyp,nz,nxh) :: scr,cxr,cyr,czr,scnr
    real(C_DOUBLE), dimension(nyp,nz,nxh) :: sci,cxi,cyi,czi,scni
    real(C_DOUBLE), dimension(mz,mx) :: scp,cxp,cyp,czp,vc

    real, dimension(nyp,mz,mx) :: scsource
    real, dimension(nyp,nz,nx) :: Qcrit
    real, dimension(qn)     :: Qmin = 0.0
    integer, dimension(qn)  :: Qx,Qy,Qz 
#ENDIF

    ! Printing variables
    real, dimension(nyp,mz,mx) :: up3d,vp3d,wp3d,wx3d,wy3d,wz3d
    real, dimension(nyp,mz,mx) :: u11p3d,u12p3d,u13p3d
    real, dimension(nyp,mz,mx) :: u21p3d,u22p3d,u23p3d
    real, dimension(nyp,mz,mx) :: u31p3d,u32p3d,u33p3d
    real, dimension(nyp,mz,mx) :: Lup3d,Lvp3d,Lwp3d
    real, dimension(nyp,mz,mx) :: swirl_3d
#IFDEF SCALAR
    real, dimension(nyp,mz,mx) :: scp3d
#ENDIF
    
    integer :: i,j,k,n,ierr
    integer :: ii,jj,ipii,jpjj,idif,jdif
    real    :: delxm,delzm,pi!,vArea
    real    :: segdrag,xsegdrag,ysegdrag,zsegdrag,wdes
    real    :: cflcheck,cflmax,swirl
    
    real,dimension(nyp)   :: cfl
    real,dimension(mz,mx) :: dragx,dragy,dragz
    
    real :: xi,zj,argx,argrad,fx,fr,rsq
    real :: uzmean(mx),uxmean(nyp)
    real :: massFlux
    
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
    integer :: cnt,kk,tnum,kmin,kmax,jmin,jmax
    logical :: condition
#ENDIF

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
    common/itime/      it
    common/dtime/      dt
    common/imat/       imatrix,kwall,kmaxsurf
    common/vortexring/ forbeta,xcenter,ycenter,zcenter,L,rad,bdyfx
    common/setup/      geomtype,flow_select,perturbtime
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

    !---------------------------------------------------------------------!
    !                    Initialize some variables                        !
    !---------------------------------------------------------------------!

    if (it .eq. irstrt) then ! Make fftw plans
        ierr = fftw_init_threads()
        call fftw_plan_with_nthreads(1)

        planY   = fftw_plan_r2r_1d(nyp,areal,areal,FFTW_REDFT00,FFTW_PATIENT)
        planXZb = fftw_plan_dft_c2r(2,[mx,mz],bcomp,breal,FFTW_PATIENT)
        planXZf = fftw_plan_dft_r2c(2,[mx,mz],breal,bcomp,FFTW_PATIENT)

        PlanZf  =     fftw_plan_many_dft(1,[mz],nxh,bcomp,[mz],1,mz,bcomp,[mz],1,mz,FFTW_FORWARD,FFTW_PATIENT)
        PlanZb  =     fftw_plan_many_dft(1,[mz],nxh,bcomp,[mz],1,mz,bcomp,[mz],1,mz,FFTW_BACKWARD,FFTW_PATIENT)
        planXf  = fftw_plan_many_dft_r2c(1,[mx],mz, breal,[mx],mz,1,bcomp,[mx],mz,1,FFTW_PATIENT)
        planXb  = fftw_plan_many_dft_c2r(1,[mx],mz, bcomp,[mx],mz,1,breal,[mx],mz,1,FFTW_PATIENT)
    end if
 
    pi = 2.0*acos(0.0)
    delxm = xl/float(mx)
    delzm = zl/float(mz)
   
#IFDEF SCALAR 
    !---------------------------------------------------------------------!
    !                    Calculate scalar source term                     !
    !---------------------------------------------------------------------!
    if (scl_flag .eq. 2 .and. it .ge. src_start .and. it .le. src_stop) then
        scsource = 0.0 !collapse(3) schedule(dynamic)
        !$omp parallel do schedule(auto) reduction(+:scsource)
        do n = 1,npart
        xc1 = xpart(n)
        yc1 = ypart(n)
        zc1 = zpart(n)

!        ! Ignore anything where beta* > 18
!        kmin = max(1,floor(1 + (xc1-6.0*sigmax)/delxm))
!        kmax = min(mx,ceiling(1 + (xc1+6.0*sigmax)/delxm))
!
!        jmin = max(1,floor(1 + (zc1-6.0*sigmaz)/delzm))
!        jmax = min(mz,ceiling(1 + (zc1+6.0*sigmaz)/delzm))
!
!        do k = kmin,kmax
        do k = 1,mx
            xsq = (float(k-1)*delxm - xc1)**2
            betax = xsq/(2.0*sigmax**2)
            if (betax .lt. 20.0) then
!            do j = jmin,jmax
            do j = 1,mz
                zsq = (float(j-1)*delzm - zc1)**2
                betaz = zsq/(2.0*sigmaz**2)
                if (betax+betaz .lt. 20.0) then
                do i = 1,nyp

                    ysq = (ycoord(i) - yc1)**2
                    betay = ysq/(2.0*sigmay**2)
       
                    if (betax + betay + betaz .lt. 18.0) then 
                        scsource(i,j,k) = scsource(i,j,k) + deltaT*exp(-(betax + betay + betaz))
                    end if
                end do
                end if
            end do 
            end if
        end do
        end do
        !$omp end parallel do

!        do n = 1,npart
!            if (it .eq. irstrt) then
!                tnum = 100+OMP_GET_THREAD_NUM() 
!                open(tnum,file='setup/particles/particles.dat',status='old',action='read')
!                do j = 1,n
!                    read(tnum,*)
!                end do
!                read(tnum,*) xc1,yc1,zc1
!                close(tnum)
!            else if (it .gt. irstrt) then
!                xc1 = xpart(n)
!                yc1 = ypart(n)
!                zc1 = zpart(n)
!                swirl = swirl_part(n)
!            end if
!    
!            ! Choose targeting case
!            if (scltarg .eq. 1) then
!                condition = swirl .gt. 2000.0 ! Target high Q
!            else if (scltarg .eq. 2) then
!                condition = swirl .lt. -2000.0 ! Target low Q
!            else
!                condition = .true. ! No targeting
!            end if
!    
!            if (condition) then 
!            do k = 1,mx
!                xsq = (float(k-1)*delxm - xc1)**2
!                betax = xsq/(2.0*sigmax**2)
!                do j = 1,mz
!                    zsq = (float(j-1)*delzm - zc1)**2
!                    betaz = zsq/(2.0*sigmaz**2)
!                    do i = 1,nyp
!                        ysq = (ycoord(i) - yc1)**2
!                        betay = ysq/(2.0*sigmay**2)
!    
!                        scsource(i,j,k) = scsource(i,j,k) + deltaT*exp(-(betax + betay + betaz))
!                    end do
!                end do
!            end do
!            end if
!        end do
!        !$omp end parallel do
#IFDEF POLYMER
    else if (scl_flag .eq. 4 .and. it .ge. src_start .and. it .le. src_stop) then
        Qx = 0
        Qy = 0
        Qz = 0
        call findmaxQ(u11,u12,u13,u21,u22,u23,u31,u32,u33,scalar,Qmin,Qx,Qy,Qz,beta3d,planY,planZb,planXb)
        scsource = 0.0
        n = 1
        cnt = 1
        do while (cnt .lt. 1024)
            kk = Qx(n)
            ii = Qy(n)
            jj = Qz(n)
            if (beta3d(ii,jj,kk) .gt. 0.8) then 
                xc1 = delxm*(kk-1)
                yc1 = ycoord(ii)
                zc1 = delzm*(jj-1)
                do k = 1,mx
                    xsq = (float(k-1)*delxm - xc1)**2
                    betax = xsq/(2.0*sigmax**2)
                    if (betax .lt. 18.0) then
                    do j = 1,mz
                        zsq = (float(j-1)*delzm - zc1)**2
                        betaz = zsq/(2.0*sigmaz**2)
                        if (betax+betaz .lt. 18.0) then
                        do i = 1,nyp
                            ysq = (ycoord(i) - yc1)**2
                            betay = ysq/(2.0*sigmay**2)

                            if (betax + betay + betaz .lt. 18.0) then
                                scsource(i,j,k) = scsource(i,j,k) + deltaT*exp(-(betax+betay+betaz))
                            end if
                        end do
                        end if
                    end do
                    end if
                end do
                cnt = cnt + 1
            end if
            n = n+1
        end do
#ENDIF
    else if (scl_flag .eq. 5 .and. it .ge. src_start .and. it .le. src_stop) then
        call newtarget(u11,u12,u13,u21,u22,u23,u31,u32,u33,scsource,planY,planZb,planXb)
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

    ! Copy spectral variables into larger arrays for transforms 
    ! Also convert Chebyshev modes to cosine modes
    !$omp parallel do default(shared) private(i,j,k,fac) schedule(auto)
    do k = 1,nxh
        do j = 1,nz 
            do i = 1,nyp

                fac = c(i)/2.0

                ! Velocity Field
                ur(i,j,k) = real(u(i,j,k))*fac
                vr(i,j,k) = real(v(i,j,k))*fac
                wr(i,j,k) = real(w(i,j,k))*fac

                ! Vorticity Field
                wxr(i,j,k) = real(omx(i,j,k))*fac
                wyr(i,j,k) = real(omy(i,j,k))*fac
                wzr(i,j,k) = real(omz(i,j,k))*fac

                ! Velocity Gradient
                u11r(i,j,k) = real(u11(i,j,k))*fac
                u12r(i,j,k) = real(u12(i,j,k))*fac
                u13r(i,j,k) = real(u13(i,j,k))*fac
                u21r(i,j,k) = real(u21(i,j,k))*fac
                u22r(i,j,k) = real(u22(i,j,k))*fac
                u23r(i,j,k) = real(u23(i,j,k))*fac
                u31r(i,j,k) = real(u31(i,j,k))*fac
                u32r(i,j,k) = real(u32(i,j,k))*fac
                u33r(i,j,k) = real(u33(i,j,k))*fac

                ! Laplacian
                Lur(i,j,k) = real(Lu(i,j,k))*fac
                Lvr(i,j,k) = real(Lv(i,j,k))*fac
                Lwr(i,j,k) = real(Lw(i,j,k))*fac

                ! Velocity Field
                ui(i,j,k) = aimag(u(i,j,k))*fac
                vi(i,j,k) = aimag(v(i,j,k))*fac
                wi(i,j,k) = aimag(w(i,j,k))*fac

                ! Vorticity Field
                wxi(i,j,k) = aimag(omx(i,j,k))*fac
                wyi(i,j,k) = aimag(omy(i,j,k))*fac
                wzi(i,j,k) = aimag(omz(i,j,k))*fac

                ! Velocity Gradient
                u11i(i,j,k) = aimag(u11(i,j,k))*fac
                u12i(i,j,k) = aimag(u12(i,j,k))*fac
                u13i(i,j,k) = aimag(u13(i,j,k))*fac
                u21i(i,j,k) = aimag(u21(i,j,k))*fac
                u22i(i,j,k) = aimag(u22(i,j,k))*fac
                u23i(i,j,k) = aimag(u23(i,j,k))*fac
                u31i(i,j,k) = aimag(u31(i,j,k))*fac
                u32i(i,j,k) = aimag(u32(i,j,k))*fac
                u33i(i,j,k) = aimag(u33(i,j,k))*fac

                ! Laplacian
                Lui(i,j,k) = aimag(Lu(i,j,k))*fac
                Lvi(i,j,k) = aimag(Lv(i,j,k))*fac
                Lwi(i,j,k) = aimag(Lw(i,j,k))*fac
#IFDEF SCALAR
                ! Scalar and its gradient
                scr(i,j,k) = real(scalar(i,j,k))*fac
                cxr(i,j,k) = real(sclx(i,j,k))*fac 
                cyr(i,j,k) = real(scly(i,j,k))*fac
                czr(i,j,k) = real(sclz(i,j,k))*fac

                sci(i,j,k) = aimag(scalar(i,j,k))*fac
                cxi(i,j,k) = aimag(sclx(i,j,k))*fac 
                cyi(i,j,k) = aimag(scly(i,j,k))*fac
                czi(i,j,k) = aimag(sclz(i,j,k))*fac
#ENDIF
#IFDEF POLYMER
                ! Conformation tensor
                c11r(i,j,k) = real(c11(i,j,k))*fac
                c12r(i,j,k) = real(c12(i,j,k))*fac
                c13r(i,j,k) = real(c13(i,j,k))*fac
                c21r(i,j,k) = real(c21(i,j,k))*fac
                c22r(i,j,k) = real(c22(i,j,k))*fac
                c23r(i,j,k) = real(c23(i,j,k))*fac
                c31r(i,j,k) = real(c31(i,j,k))*fac
                c32r(i,j,k) = real(c32(i,j,k))*fac
                c33r(i,j,k) = real(c33(i,j,k))*fac
               
                ! Conformation tensor gradient
                dc111r(i,j,k) = real(dc111(i,j,k))*fac 
                dc112r(i,j,k) = real(dc112(i,j,k))*fac 
                dc113r(i,j,k) = real(dc113(i,j,k))*fac 
                dc121r(i,j,k) = real(dc121(i,j,k))*fac 
                dc122r(i,j,k) = real(dc122(i,j,k))*fac 
                dc123r(i,j,k) = real(dc123(i,j,k))*fac 
                dc131r(i,j,k) = real(dc131(i,j,k))*fac 
                dc132r(i,j,k) = real(dc132(i,j,k))*fac 
                dc133r(i,j,k) = real(dc133(i,j,k))*fac 

                dc211r(i,j,k) = real(dc211(i,j,k))*fac 
                dc212r(i,j,k) = real(dc212(i,j,k))*fac 
                dc213r(i,j,k) = real(dc213(i,j,k))*fac 
                dc221r(i,j,k) = real(dc221(i,j,k))*fac 
                dc222r(i,j,k) = real(dc222(i,j,k))*fac 
                dc223r(i,j,k) = real(dc223(i,j,k))*fac 
                dc231r(i,j,k) = real(dc231(i,j,k))*fac 
                dc232r(i,j,k) = real(dc232(i,j,k))*fac 
                dc233r(i,j,k) = real(dc233(i,j,k))*fac 

                dc311r(i,j,k) = real(dc311(i,j,k))*fac 
                dc312r(i,j,k) = real(dc312(i,j,k))*fac 
                dc313r(i,j,k) = real(dc313(i,j,k))*fac 
                dc321r(i,j,k) = real(dc321(i,j,k))*fac 
                dc322r(i,j,k) = real(dc322(i,j,k))*fac 
                dc323r(i,j,k) = real(dc323(i,j,k))*fac 
                dc331r(i,j,k) = real(dc331(i,j,k))*fac 
                dc332r(i,j,k) = real(dc332(i,j,k))*fac 
                dc333r(i,j,k) = real(dc333(i,j,k))*fac 

                ! Conformation tensor
                c11i(i,j,k) = aimag(c11(i,j,k))*fac
                c12i(i,j,k) = aimag(c12(i,j,k))*fac
                c13i(i,j,k) = aimag(c13(i,j,k))*fac
                c21i(i,j,k) = aimag(c21(i,j,k))*fac
                c22i(i,j,k) = aimag(c22(i,j,k))*fac
                c23i(i,j,k) = aimag(c23(i,j,k))*fac
                c31i(i,j,k) = aimag(c31(i,j,k))*fac
                c32i(i,j,k) = aimag(c32(i,j,k))*fac
                c33i(i,j,k) = aimag(c33(i,j,k))*fac

                ! Conformation tensor gradient
                dc111i(i,j,k) = aimag(dc111(i,j,k))*fac 
                dc112i(i,j,k) = aimag(dc112(i,j,k))*fac 
                dc113i(i,j,k) = aimag(dc113(i,j,k))*fac 
                dc121i(i,j,k) = aimag(dc121(i,j,k))*fac 
                dc122i(i,j,k) = aimag(dc122(i,j,k))*fac 
                dc123i(i,j,k) = aimag(dc123(i,j,k))*fac 
                dc131i(i,j,k) = aimag(dc131(i,j,k))*fac 
                dc132i(i,j,k) = aimag(dc132(i,j,k))*fac 
                dc133i(i,j,k) = aimag(dc133(i,j,k))*fac 

                dc211i(i,j,k) = aimag(dc211(i,j,k))*fac 
                dc212i(i,j,k) = aimag(dc212(i,j,k))*fac 
                dc213i(i,j,k) = aimag(dc213(i,j,k))*fac 
                dc221i(i,j,k) = aimag(dc221(i,j,k))*fac 
                dc222i(i,j,k) = aimag(dc222(i,j,k))*fac 
                dc223i(i,j,k) = aimag(dc223(i,j,k))*fac 
                dc231i(i,j,k) = aimag(dc231(i,j,k))*fac 
                dc232i(i,j,k) = aimag(dc232(i,j,k))*fac 
                dc233i(i,j,k) = aimag(dc233(i,j,k))*fac 

                dc311i(i,j,k) = aimag(dc311(i,j,k))*fac 
                dc312i(i,j,k) = aimag(dc312(i,j,k))*fac 
                dc313i(i,j,k) = aimag(dc313(i,j,k))*fac 
                dc321i(i,j,k) = aimag(dc321(i,j,k))*fac 
                dc322i(i,j,k) = aimag(dc322(i,j,k))*fac 
                dc323i(i,j,k) = aimag(dc323(i,j,k))*fac 
                dc331i(i,j,k) = aimag(dc331(i,j,k))*fac 
                dc332i(i,j,k) = aimag(dc332(i,j,k))*fac 
                dc333i(i,j,k) = aimag(dc333(i,j,k))*fac 
#ENDIF
            end do
        end do
    end do
    !$omp end parallel do

    ! Compute Real --> Real DCT-I on real and imaginary parts of spectral variables
    !$omp parallel do simd default(shared) private(j,k) collapse(2) schedule(simd: auto)
    do k = 1,nxh
        do j = 1,nz
            call fftw_execute_r2r(planY,ur(:,j,k),ur(:,j,k))
            call fftw_execute_r2r(planY,vr(:,j,k),vr(:,j,k))
            call fftw_execute_r2r(planY,wr(:,j,k),wr(:,j,k))
            
            call fftw_execute_r2r(planY,wxr(:,j,k),wxr(:,j,k))
            call fftw_execute_r2r(planY,wyr(:,j,k),wyr(:,j,k))
            call fftw_execute_r2r(planY,wzr(:,j,k),wzr(:,j,k))

            call fftw_execute_r2r(planY,u11r(:,j,k),u11r(:,j,k))
            call fftw_execute_r2r(planY,u12r(:,j,k),u12r(:,j,k))
            call fftw_execute_r2r(planY,u13r(:,j,k),u13r(:,j,k))
            call fftw_execute_r2r(planY,u21r(:,j,k),u21r(:,j,k))
            call fftw_execute_r2r(planY,u22r(:,j,k),u22r(:,j,k))
            call fftw_execute_r2r(planY,u23r(:,j,k),u23r(:,j,k))
            call fftw_execute_r2r(planY,u31r(:,j,k),u31r(:,j,k))
            call fftw_execute_r2r(planY,u32r(:,j,k),u32r(:,j,k))
            call fftw_execute_r2r(planY,u33r(:,j,k),u33r(:,j,k))

            call fftw_execute_r2r(planY,Lur(:,j,k),Lur(:,j,k))
            call fftw_execute_r2r(planY,Lvr(:,j,k),Lvr(:,j,k))
            call fftw_execute_r2r(planY,Lwr(:,j,k),Lwr(:,j,k))

            call fftw_execute_r2r(planY,  ui(:,j,k),  ui(:,j,k))
            call fftw_execute_r2r(planY,  vi(:,j,k),  vi(:,j,k))
            call fftw_execute_r2r(planY,  wi(:,j,k),  wi(:,j,k))
            
            call fftw_execute_r2r(planY, wxi(:,j,k), wxi(:,j,k))
            call fftw_execute_r2r(planY, wyi(:,j,k), wyi(:,j,k))
            call fftw_execute_r2r(planY, wzi(:,j,k), wzi(:,j,k))

            call fftw_execute_r2r(planY,u11i(:,j,k),u11i(:,j,k))
            call fftw_execute_r2r(planY,u12i(:,j,k),u12i(:,j,k))
            call fftw_execute_r2r(planY,u13i(:,j,k),u13i(:,j,k))
            call fftw_execute_r2r(planY,u21i(:,j,k),u21i(:,j,k))
            call fftw_execute_r2r(planY,u22i(:,j,k),u22i(:,j,k))
            call fftw_execute_r2r(planY,u23i(:,j,k),u23i(:,j,k))
            call fftw_execute_r2r(planY,u31i(:,j,k),u31i(:,j,k))
            call fftw_execute_r2r(planY,u32i(:,j,k),u32i(:,j,k))
            call fftw_execute_r2r(planY,u33i(:,j,k),u33i(:,j,k))

            call fftw_execute_r2r(planY, Lui(:,j,k), Lui(:,j,k))
            call fftw_execute_r2r(planY, Lvi(:,j,k), Lvi(:,j,k))
            call fftw_execute_r2r(planY, Lwi(:,j,k), Lwi(:,j,k))
#IFDEF SCALAR
            call fftw_execute_r2r(planY,scr(:,j,k),scr(:,j,k))
            call fftw_execute_r2r(planY,cxr(:,j,k),cxr(:,j,k))
            call fftw_execute_r2r(planY,cyr(:,j,k),cyr(:,j,k))
            call fftw_execute_r2r(planY,czr(:,j,k),czr(:,j,k))

            call fftw_execute_r2r(planY,sci(:,j,k),sci(:,j,k))
            call fftw_execute_r2r(planY,cxi(:,j,k),cxi(:,j,k))
            call fftw_execute_r2r(planY,cyi(:,j,k),cyi(:,j,k))
            call fftw_execute_r2r(planY,czi(:,j,k),czi(:,j,k))
#ENDIF
#IFDEF POLYMER
            if (it .ge. (src_start-1)) then
            ! Real part
            call fftw_execute_r2r(planY,c11r(:,j,k),c11r(:,j,k))
            call fftw_execute_r2r(planY,c12r(:,j,k),c12r(:,j,k))
            call fftw_execute_r2r(planY,c13r(:,j,k),c13r(:,j,k))
            call fftw_execute_r2r(planY,c21r(:,j,k),c21r(:,j,k))
            call fftw_execute_r2r(planY,c22r(:,j,k),c22r(:,j,k))
            call fftw_execute_r2r(planY,c23r(:,j,k),c23r(:,j,k))
            call fftw_execute_r2r(planY,c31r(:,j,k),c31r(:,j,k))
            call fftw_execute_r2r(planY,c32r(:,j,k),c32r(:,j,k))
            call fftw_execute_r2r(planY,c33r(:,j,k),c33r(:,j,k))

            call fftw_execute_r2r(planY,dc111r(:,j,k),dc111r(:,j,k))
            call fftw_execute_r2r(planY,dc112r(:,j,k),dc112r(:,j,k))
            call fftw_execute_r2r(planY,dc113r(:,j,k),dc113r(:,j,k))
            call fftw_execute_r2r(planY,dc121r(:,j,k),dc121r(:,j,k))
            call fftw_execute_r2r(planY,dc122r(:,j,k),dc122r(:,j,k))
            call fftw_execute_r2r(planY,dc123r(:,j,k),dc123r(:,j,k))
            call fftw_execute_r2r(planY,dc131r(:,j,k),dc131r(:,j,k))
            call fftw_execute_r2r(planY,dc132r(:,j,k),dc132r(:,j,k))
            call fftw_execute_r2r(planY,dc133r(:,j,k),dc133r(:,j,k))

            call fftw_execute_r2r(planY,dc211r(:,j,k),dc211r(:,j,k))
            call fftw_execute_r2r(planY,dc212r(:,j,k),dc212r(:,j,k))
            call fftw_execute_r2r(planY,dc213r(:,j,k),dc213r(:,j,k))
            call fftw_execute_r2r(planY,dc221r(:,j,k),dc221r(:,j,k))
            call fftw_execute_r2r(planY,dc222r(:,j,k),dc222r(:,j,k))
            call fftw_execute_r2r(planY,dc223r(:,j,k),dc223r(:,j,k))
            call fftw_execute_r2r(planY,dc231r(:,j,k),dc231r(:,j,k))
            call fftw_execute_r2r(planY,dc232r(:,j,k),dc232r(:,j,k))
            call fftw_execute_r2r(planY,dc233r(:,j,k),dc233r(:,j,k))

            call fftw_execute_r2r(planY,dc311r(:,j,k),dc311r(:,j,k))
            call fftw_execute_r2r(planY,dc312r(:,j,k),dc312r(:,j,k))
            call fftw_execute_r2r(planY,dc313r(:,j,k),dc313r(:,j,k))
            call fftw_execute_r2r(planY,dc321r(:,j,k),dc321r(:,j,k))
            call fftw_execute_r2r(planY,dc322r(:,j,k),dc322r(:,j,k))
            call fftw_execute_r2r(planY,dc323r(:,j,k),dc323r(:,j,k))
            call fftw_execute_r2r(planY,dc331r(:,j,k),dc331r(:,j,k))
            call fftw_execute_r2r(planY,dc332r(:,j,k),dc332r(:,j,k))
            call fftw_execute_r2r(planY,dc333r(:,j,k),dc333r(:,j,k))

            ! Imaginary part
            call fftw_execute_r2r(planY,c11i(:,j,k),c11i(:,j,k))
            call fftw_execute_r2r(planY,c12i(:,j,k),c12i(:,j,k))
            call fftw_execute_r2r(planY,c13i(:,j,k),c13i(:,j,k))
            call fftw_execute_r2r(planY,c21i(:,j,k),c21i(:,j,k))
            call fftw_execute_r2r(planY,c22i(:,j,k),c22i(:,j,k))
            call fftw_execute_r2r(planY,c23i(:,j,k),c23i(:,j,k))
            call fftw_execute_r2r(planY,c31i(:,j,k),c31i(:,j,k))
            call fftw_execute_r2r(planY,c32i(:,j,k),c32i(:,j,k))
            call fftw_execute_r2r(planY,c33i(:,j,k),c33i(:,j,k))

            call fftw_execute_r2r(planY,dc111i(:,j,k),dc111i(:,j,k))
            call fftw_execute_r2r(planY,dc112i(:,j,k),dc112i(:,j,k))
            call fftw_execute_r2r(planY,dc113i(:,j,k),dc113i(:,j,k))
            call fftw_execute_r2r(planY,dc121i(:,j,k),dc121i(:,j,k))
            call fftw_execute_r2r(planY,dc122i(:,j,k),dc122i(:,j,k))
            call fftw_execute_r2r(planY,dc123i(:,j,k),dc123i(:,j,k))
            call fftw_execute_r2r(planY,dc131i(:,j,k),dc131i(:,j,k))
            call fftw_execute_r2r(planY,dc132i(:,j,k),dc132i(:,j,k))
            call fftw_execute_r2r(planY,dc133i(:,j,k),dc133i(:,j,k))

            call fftw_execute_r2r(planY,dc211i(:,j,k),dc211i(:,j,k))
            call fftw_execute_r2r(planY,dc212i(:,j,k),dc212i(:,j,k))
            call fftw_execute_r2r(planY,dc213i(:,j,k),dc213i(:,j,k))
            call fftw_execute_r2r(planY,dc221i(:,j,k),dc221i(:,j,k))
            call fftw_execute_r2r(planY,dc222i(:,j,k),dc222i(:,j,k))
            call fftw_execute_r2r(planY,dc223i(:,j,k),dc223i(:,j,k))
            call fftw_execute_r2r(planY,dc231i(:,j,k),dc231i(:,j,k))
            call fftw_execute_r2r(planY,dc232i(:,j,k),dc232i(:,j,k))
            call fftw_execute_r2r(planY,dc233i(:,j,k),dc233i(:,j,k))

            call fftw_execute_r2r(planY,dc311i(:,j,k),dc311i(:,j,k))
            call fftw_execute_r2r(planY,dc312i(:,j,k),dc312i(:,j,k))
            call fftw_execute_r2r(planY,dc313i(:,j,k),dc313i(:,j,k))
            call fftw_execute_r2r(planY,dc321i(:,j,k),dc321i(:,j,k))
            call fftw_execute_r2r(planY,dc322i(:,j,k),dc322i(:,j,k))
            call fftw_execute_r2r(planY,dc323i(:,j,k),dc323i(:,j,k))
            call fftw_execute_r2r(planY,dc331i(:,j,k),dc331i(:,j,k))
            call fftw_execute_r2r(planY,dc332i(:,j,k),dc332i(:,j,k))
            call fftw_execute_r2r(planY,dc333i(:,j,k),dc333i(:,j,k))
            end if
#ENDIF
        end do
    end do
    !$omp end parallel do simd

    ! Loop over y-planes to do the Fourier transforms and compute nonlinear terms

    !$omp parallel do simd default(shared) private(us,vs,ws,wxs,wys,wzs,u11s,u12s,  &
    !$omp       u13s,u21s,u22s,u23s,u31s,u32s,u33s,up,vp,wp,wxp,wyp,wzp,u11p,  &
    !$omp       u12p,u13p,u21p,u22p,u23p,u31p,u32p,u33p,i,j,k,ii,jj,ipii,jpjj, &
    !$omp       idif,jdif,cflcheck,wdes,zj,rsq,argrad,fr,xi,argx,fx,segdrag,   &
    !$omp       Lus,Lvs,Lws,Lup,Lvp,Lwp, &
#IFDEF SCALAR
    !$omp       scs,cxs,cys,czs,scp,cxp,cyp,czp,vc,vcs,                        &
#ENDIF
#IFDEF POLYMER
    !$omp       c11s,c12s,c13s,c21s,c22s,c23s,c31s,c32s,c33s,                   &
    !$omp       c11p,c12p,c13p,c21p,c22p,c23p,c31p,c32p,c33p,                   &
    !$omp       dc111s,dc112s,dc113s,dc121s,dc122s,dc123s,dc131s,dc132s,dc133s, &
    !$omp       dc111p,dc112p,dc113p,dc121p,dc122p,dc123p,dc131p,dc132p,dc133p, &
    !$omp       dc211s,dc212s,dc213s,dc221s,dc222s,dc223s,dc231s,dc232s,dc233s, &
    !$omp       dc211p,dc212p,dc213p,dc221p,dc222p,dc223p,dc231p,dc232p,dc233p, &
    !$omp       dc311s,dc312s,dc313s,dc321s,dc322s,dc323s,dc331s,dc332s,dc333s, &
    !$omp       dc311p,dc312p,dc313p,dc321p,dc322p,dc323p,dc331p,dc332p,dc333p, &
    !$omp       c11np,c12np,c13np,c22np,c23np,c33np,                            &
    !$omp       c11ns,c12ns,c13ns,c22ns,c23ns,c33ns,                            &
    !$omp       str11np,str12np,str13np,str22np,str23np,str33np,                &
    !$omp       str11ns,str12ns,str13ns,str22ns,str23ns,str33ns,                &
    !$omp       qp11np,qp12np,qp13np,qp22np,qp23np,qp33np,                      &
    !$omp       qp11s,qp12s,qp13s,qp22s,qp23s,qp33s,                            &
    !$omp       trp,zbeta1,beta_poly,                                           &
#ENDIF
    !$omp       dragx,dragy,dragz,vwx,vwy,vwz,vwxs,vwys,vwzs) schedule(simd: auto)
    do i = 1,nyp

        ! Zero out 2D variables
        us = 0.0
        vs = 0.0
        ws = 0.0

        wxs = 0.0
        wys = 0.0
        wzs = 0.0

        u11s = 0.0
        u12s = 0.0
        u13s = 0.0
        u21s = 0.0
        u22s = 0.0
        u23s = 0.0
        u31s = 0.0
        u32s = 0.0
        u33s = 0.0

        Lus = 0.0
        Lvs = 0.0
        Lws = 0.0

#IFDEF SCALAR
        scs = 0.0
        cxs = 0.0
        cys = 0.0
        czs = 0.0
#ENDIF
#IFDEF POLYMER
        if (it .ge. (src_start-1)) then
        c11s = 0.0
        c12s = 0.0
        c13s = 0.0
        c21s = 0.0
        c22s = 0.0
        c23s = 0.0
        c31s = 0.0
        c32s = 0.0
        c33s = 0.0

        dc111s = 0.0
        dc112s = 0.0
        dc113s = 0.0
        dc121s = 0.0
        dc122s = 0.0
        dc123s = 0.0
        dc131s = 0.0
        dc132s = 0.0
        dc133s = 0.0

        dc211s = 0.0
        dc212s = 0.0
        dc213s = 0.0
        dc221s = 0.0
        dc222s = 0.0
        dc223s = 0.0
        dc231s = 0.0
        dc232s = 0.0
        dc233s = 0.0

        dc311s = 0.0
        dc312s = 0.0
        dc313s = 0.0
        dc321s = 0.0
        dc322s = 0.0
        dc323s = 0.0
        dc331s = 0.0
        dc332s = 0.0
        dc333s = 0.0
        end if
#ENDIF

        ! Copy data to 2D variables
        do k = 1,nxh
            do j = 1,nz
                if (j .le. nzh) jj = j
                if (j .gt. nzh) jj = (mz-nz) + j
                us(jj,k) = cmplx(ur(i,j,k),ui(i,j,k))
                vs(jj,k) = cmplx(vr(i,j,k),vi(i,j,k))
                ws(jj,k) = cmplx(wr(i,j,k),wi(i,j,k))

                wxs(jj,k) = cmplx(wxr(i,j,k),wxi(i,j,k))
                wys(jj,k) = cmplx(wyr(i,j,k),wyi(i,j,k))
                wzs(jj,k) = cmplx(wzr(i,j,k),wzi(i,j,k))

                u11s(jj,k) = cmplx(u11r(i,j,k),u11i(i,j,k))
                u12s(jj,k) = cmplx(u12r(i,j,k),u12i(i,j,k))
                u13s(jj,k) = cmplx(u13r(i,j,k),u13i(i,j,k))
                u21s(jj,k) = cmplx(u21r(i,j,k),u21i(i,j,k))
                u22s(jj,k) = cmplx(u22r(i,j,k),u22i(i,j,k))
                u23s(jj,k) = cmplx(u23r(i,j,k),u23i(i,j,k))
                u31s(jj,k) = cmplx(u31r(i,j,k),u31i(i,j,k))
                u32s(jj,k) = cmplx(u32r(i,j,k),u32i(i,j,k))
                u33s(jj,k) = cmplx(u33r(i,j,k),u33i(i,j,k))

                Lus(jj,k) = cmplx(Lur(i,j,k),Lui(i,j,k))
                Lvs(jj,k) = cmplx(Lvr(i,j,k),Lvi(i,j,k))
                Lws(jj,k) = cmplx(Lwr(i,j,k),Lwi(i,j,k))
#IFDEF SCALAR
                scs(jj,k) = cmplx(scr(i,j,k),sci(i,j,k))
                cxs(jj,k) = cmplx(cxr(i,j,k),cxi(i,j,k))
                cys(jj,k) = cmplx(cyr(i,j,k),cyi(i,j,k))
                czs(jj,k) = cmplx(czr(i,j,k),czi(i,j,k))
#ENDIF
#IFDEF POLYMER
                if (it .ge. (src_start-1)) then
                c11s(jj,k) = cmplx(c11r(i,j,k),c11i(i,j,k))
                c12s(jj,k) = cmplx(c12r(i,j,k),c12i(i,j,k))
                c13s(jj,k) = cmplx(c13r(i,j,k),c13i(i,j,k))
                c21s(jj,k) = cmplx(c21r(i,j,k),c21i(i,j,k))
                c22s(jj,k) = cmplx(c22r(i,j,k),c22i(i,j,k))
                c23s(jj,k) = cmplx(c23r(i,j,k),c23i(i,j,k))
                c31s(jj,k) = cmplx(c31r(i,j,k),c31i(i,j,k))
                c32s(jj,k) = cmplx(c32r(i,j,k),c32i(i,j,k))
                c33s(jj,k) = cmplx(c33r(i,j,k),c33i(i,j,k))

                dc111s(jj,k) = cmplx(dc111r(i,j,k),dc111i(i,j,k))
                dc112s(jj,k) = cmplx(dc112r(i,j,k),dc112i(i,j,k))
                dc113s(jj,k) = cmplx(dc113r(i,j,k),dc113i(i,j,k))
                dc121s(jj,k) = cmplx(dc121r(i,j,k),dc121i(i,j,k))
                dc122s(jj,k) = cmplx(dc122r(i,j,k),dc122i(i,j,k))
                dc123s(jj,k) = cmplx(dc123r(i,j,k),dc123i(i,j,k))
                dc131s(jj,k) = cmplx(dc131r(i,j,k),dc131i(i,j,k))
                dc132s(jj,k) = cmplx(dc132r(i,j,k),dc132i(i,j,k))
                dc133s(jj,k) = cmplx(dc133r(i,j,k),dc133i(i,j,k))

                dc211s(jj,k) = cmplx(dc211r(i,j,k),dc211i(i,j,k))
                dc212s(jj,k) = cmplx(dc212r(i,j,k),dc212i(i,j,k))
                dc213s(jj,k) = cmplx(dc213r(i,j,k),dc213i(i,j,k))
                dc221s(jj,k) = cmplx(dc221r(i,j,k),dc221i(i,j,k))
                dc222s(jj,k) = cmplx(dc222r(i,j,k),dc222i(i,j,k))
                dc223s(jj,k) = cmplx(dc223r(i,j,k),dc223i(i,j,k))
                dc231s(jj,k) = cmplx(dc231r(i,j,k),dc231i(i,j,k))
                dc232s(jj,k) = cmplx(dc232r(i,j,k),dc232i(i,j,k))
                dc233s(jj,k) = cmplx(dc233r(i,j,k),dc233i(i,j,k))

                dc311s(jj,k) = cmplx(dc311r(i,j,k),dc311i(i,j,k))
                dc312s(jj,k) = cmplx(dc312r(i,j,k),dc312i(i,j,k))
                dc313s(jj,k) = cmplx(dc313r(i,j,k),dc313i(i,j,k))
                dc321s(jj,k) = cmplx(dc321r(i,j,k),dc321i(i,j,k))
                dc322s(jj,k) = cmplx(dc322r(i,j,k),dc322i(i,j,k))
                dc323s(jj,k) = cmplx(dc323r(i,j,k),dc323i(i,j,k))
                dc331s(jj,k) = cmplx(dc331r(i,j,k),dc331i(i,j,k))
                dc332s(jj,k) = cmplx(dc332r(i,j,k),dc332i(i,j,k))
                dc333s(jj,k) = cmplx(dc333r(i,j,k),dc333i(i,j,k))
                end if
#ENDIF
            end do
        end do
        
        ! Complex --> Complex z-transform
        call fftw_execute_dft(planZb,  us,  us)
        call fftw_execute_dft(planZb,  vs,  vs)
        call fftw_execute_dft(planZb,  ws,  ws)

        call fftw_execute_dft(planZb, wxs, wxs)
        call fftw_execute_dft(planZb, wys, wys)
        call fftw_execute_dft(planZb, wzs, wzs)

        call fftw_execute_dft(planZb,u11s,u11s)
        call fftw_execute_dft(planZb,u12s,u12s)
        call fftw_execute_dft(planZb,u13s,u13s)
        call fftw_execute_dft(planZb,u21s,u21s)
        call fftw_execute_dft(planZb,u22s,u22s)
        call fftw_execute_dft(planZb,u23s,u23s)
        call fftw_execute_dft(planZb,u31s,u31s)
        call fftw_execute_dft(planZb,u32s,u32s)
        call fftw_execute_dft(planZb,u33s,u33s)

        call fftw_execute_dft(planZb, Lus, Lus)
        call fftw_execute_dft(planZb, Lvs, Lvs)
        call fftw_execute_dft(planZb, Lws, Lws)
#IFDEF SCALAR
        ! Scalar field and its gradient
        call fftw_execute_dft(planZb,scs,scs)
        call fftw_execute_dft(planZb,cxs,cxs)
        call fftw_execute_dft(planZb,cys,cys)
        call fftw_execute_dft(planZb,czs,czs)
            
#ENDIF
#IFDEF POLYMER
        if (it .ge. (src_start-1)) then
        ! Conformation tensor
        call fftw_execute_dft(planZb,c11s,c11s)
        call fftw_execute_dft(planZb,c12s,c12s)
        call fftw_execute_dft(planZb,c13s,c13s)
        call fftw_execute_dft(planZb,c21s,c21s)
        call fftw_execute_dft(planZb,c22s,c22s)
        call fftw_execute_dft(planZb,c23s,c23s)
        call fftw_execute_dft(planZb,c31s,c31s)
        call fftw_execute_dft(planZb,c32s,c32s)
        call fftw_execute_dft(planZb,c33s,c33s)

        ! Conformation tensor gradient
        call fftw_execute_dft(planZb,dc111s,dc111s)
        call fftw_execute_dft(planZb,dc112s,dc112s)
        call fftw_execute_dft(planZb,dc113s,dc113s)
        call fftw_execute_dft(planZb,dc121s,dc121s)
        call fftw_execute_dft(planZb,dc122s,dc122s)
        call fftw_execute_dft(planZb,dc123s,dc123s)
        call fftw_execute_dft(planZb,dc131s,dc131s)
        call fftw_execute_dft(planZb,dc132s,dc132s)
        call fftw_execute_dft(planZb,dc133s,dc133s)

        call fftw_execute_dft(planZb,dc211s,dc211s)
        call fftw_execute_dft(planZb,dc212s,dc212s)
        call fftw_execute_dft(planZb,dc213s,dc213s)
        call fftw_execute_dft(planZb,dc221s,dc221s)
        call fftw_execute_dft(planZb,dc222s,dc222s)
        call fftw_execute_dft(planZb,dc223s,dc223s)
        call fftw_execute_dft(planZb,dc231s,dc231s)
        call fftw_execute_dft(planZb,dc232s,dc232s)
        call fftw_execute_dft(planZb,dc233s,dc233s)

        call fftw_execute_dft(planZb,dc311s,dc311s)
        call fftw_execute_dft(planZb,dc312s,dc312s)
        call fftw_execute_dft(planZb,dc313s,dc313s)
        call fftw_execute_dft(planZb,dc321s,dc321s)
        call fftw_execute_dft(planZb,dc322s,dc322s)
        call fftw_execute_dft(planZb,dc323s,dc323s)
        call fftw_execute_dft(planZb,dc331s,dc331s)
        call fftw_execute_dft(planZb,dc332s,dc332s)
        call fftw_execute_dft(planZb,dc333s,dc333s)
        end if
#ENDIF

        ! Complex --> Real x-transform 
        call fftw_execute_dft_c2r(planXb,  us,  up)
        call fftw_execute_dft_c2r(planXb,  vs,  vp)
        call fftw_execute_dft_c2r(planXb,  ws,  wp)
                                                  
        call fftw_execute_dft_c2r(planXb, wxs, wxp)
        call fftw_execute_dft_c2r(planXb, wys, wyp)
        call fftw_execute_dft_c2r(planXb, wzs, wzp)
                                                  
        call fftw_execute_dft_c2r(planXb,u11s,u11p)
        call fftw_execute_dft_c2r(planXb,u12s,u12p)
        call fftw_execute_dft_c2r(planXb,u13s,u13p)
        call fftw_execute_dft_c2r(planXb,u21s,u21p)
        call fftw_execute_dft_c2r(planXb,u22s,u22p)
        call fftw_execute_dft_c2r(planXb,u23s,u23p)
        call fftw_execute_dft_c2r(planXb,u31s,u31p)
        call fftw_execute_dft_c2r(planXb,u32s,u32p)
        call fftw_execute_dft_c2r(planXb,u33s,u33p)
                                                  
        call fftw_execute_dft_c2r(planXb, Lus, Lup)
        call fftw_execute_dft_c2r(planXb, Lvs, Lvp)
        call fftw_execute_dft_c2r(planXb, Lws, Lwp)
#IFDEF SCALAR
        ! Scalar field and its gradient
        call fftw_execute_dft_c2r(planXb,scs,scp)
        call fftw_execute_dft_c2r(planXb,cxs,cxp)
        call fftw_execute_dft_c2r(planXb,cys,cyp)
        call fftw_execute_dft_c2r(planXb,czs,czp)
            
#ENDIF
#IFDEF POLYMER
        if (it .ge. (src_start-1)) then
        ! Conformation tensor
        call fftw_execute_dft_c2r(planXb,c11s,c11p)
        call fftw_execute_dft_c2r(planXb,c12s,c12p)
        call fftw_execute_dft_c2r(planXb,c13s,c13p)
        call fftw_execute_dft_c2r(planXb,c21s,c21p)
        call fftw_execute_dft_c2r(planXb,c22s,c22p)
        call fftw_execute_dft_c2r(planXb,c23s,c23p)
        call fftw_execute_dft_c2r(planXb,c31s,c31p)
        call fftw_execute_dft_c2r(planXb,c32s,c32p)
        call fftw_execute_dft_c2r(planXb,c33s,c33p)

        ! Conformation tensor gradient
        call fftw_execute_dft_c2r(planXb,dc111s,dc111p)
        call fftw_execute_dft_c2r(planXb,dc112s,dc112p)
        call fftw_execute_dft_c2r(planXb,dc113s,dc113p)
        call fftw_execute_dft_c2r(planXb,dc121s,dc121p)
        call fftw_execute_dft_c2r(planXb,dc122s,dc122p)
        call fftw_execute_dft_c2r(planXb,dc123s,dc123p)
        call fftw_execute_dft_c2r(planXb,dc131s,dc131p)
        call fftw_execute_dft_c2r(planXb,dc132s,dc132p)
        call fftw_execute_dft_c2r(planXb,dc133s,dc133p)

        call fftw_execute_dft_c2r(planXb,dc211s,dc211p)
        call fftw_execute_dft_c2r(planXb,dc212s,dc212p)
        call fftw_execute_dft_c2r(planXb,dc213s,dc213p)
        call fftw_execute_dft_c2r(planXb,dc221s,dc221p)
        call fftw_execute_dft_c2r(planXb,dc222s,dc222p)
        call fftw_execute_dft_c2r(planXb,dc223s,dc223p)
        call fftw_execute_dft_c2r(planXb,dc231s,dc231p)
        call fftw_execute_dft_c2r(planXb,dc232s,dc232p)
        call fftw_execute_dft_c2r(planXb,dc233s,dc233p)

        call fftw_execute_dft_c2r(planXb,dc311s,dc311p)
        call fftw_execute_dft_c2r(planXb,dc312s,dc312p)
        call fftw_execute_dft_c2r(planXb,dc313s,dc313p)
        call fftw_execute_dft_c2r(planXb,dc321s,dc321p)
        call fftw_execute_dft_c2r(planXb,dc322s,dc322p)
        call fftw_execute_dft_c2r(planXb,dc323s,dc323p)
        call fftw_execute_dft_c2r(planXb,dc331s,dc331p)
        call fftw_execute_dft_c2r(planXb,dc332s,dc332p)
        call fftw_execute_dft_c2r(planXb,dc333s,dc333p)
        end if
#ENDIF
        ! Compute nonlinear and IBF terms

        ! Initialize the force field array for IBF
        dragx = 0.0
        dragy = 0.0
        dragz = 0.0

            !---------------------------------------------------------------------!
            !        Compute the cross product of velocity and vorticity          !
            !---------------------------------------------------------------------!
   
        do k = 1,mx
            do j = 1,mz
                vwx(j,k) =  vp(j,k)*wzp(j,k) - wp(j,k)*wyp(j,k)
                vwy(j,k) = -up(j,k)*wzp(j,k) + wp(j,k)*wxp(j,k)
                vwz(j,k) =  up(j,k)*wyp(j,k) - vp(j,k)*wxp(j,k)
    
#IFDEF SCALAR
                 vc(j,k) = -(up(j,k)*cxp(j,k) + vp(j,k)*cyp(j,k) + wp(j,k)*czp(j,k))
#ENDIF
#IFDEF POLYMER    
                if (it .ge. src_start-1) then
                c11np(j,k)=c11p(j,k)*u11p(j,k)+c12p(j,k)*u12p(j,k)+c13p(j,k)*u13p(j,k)+  &
                             c11p(j,k)*u11p(j,k)+c21p(j,k)*u12p(j,k)+c31p(j,k)*u13p(j,k)-  &
                            (up(j,k)*dc111p(j,k)+vp(j,k)*dc112p(j,k)+wp(j,k)*dc113p(j,k))            
                  
                c12np(j,k)=c11p(j,k)*u21p(j,k)+c12p(j,k)*u22p(j,k)+c13p(j,k)*u23p(j,k)+   &
                             c12p(j,k)*u11p(j,k)+c22p(j,k)*u12p(j,k)+c32p(j,k)*u13p(j,k)-   &
                            (up(j,k)*dc121p(j,k)+vp(j,k)*dc122p(j,k)+wp(j,k)*dc123p(j,k))
                  
                c13np(j,k)=c11p(j,k)*u31p(j,k)+c12p(j,k)*u32p(j,k)+c13p(j,k)*u33p(j,k)+   &
                             c13p(j,k)*u11p(j,k)+c23p(j,k)*u12p(j,k)+c33p(j,k)*u13p(j,k)-   &
                             (up(j,k)*dc131p(j,k)+vp(j,k)*dc132p(j,k)+wp(j,k)*dc133p(j,k))
                  
                c22np(j,k)=c21p(j,k)*u21p(j,k)+c22p(j,k)*u22p(j,k)+c23p(j,k)*u23p(j,k)+   &
                             c12p(j,k)*u21p(j,k)+c22p(j,k)*u22p(j,k)+c32p(j,k)*u23p(j,k)-   &
                             (up(j,k)*dc221p(j,k)+vp(j,k)*dc222p(j,k)+wp(j,k)*dc223p(j,k))
                  
                c23np(j,k)=c21p(j,k)*u31p(j,k)+c22p(j,k)*u32p(j,k)+c23p(j,k)*u33p(j,k)+   &
                             c13p(j,k)*u21p(j,k)+c23p(j,k)*u22p(j,k)+c33p(j,k)*u23p(j,k)-   &
                             (up(j,k)*dc231p(j,k)+vp(j,k)*dc232p(j,k)+wp(j,k)*dc233p(j,k))

                c33np(j,k)=c31p(j,k)*u31p(j,k)+c32p(j,k)*u32p(j,k)+c33p(j,k)*u33p(j,k)+   &         
                             c13p(j,k)*u31p(j,k)+c23p(j,k)*u32p(j,k)+c33p(j,k)*u33p(j,k)-   &
                            (up(j,k)*dc331p(j,k)+vp(j,k)*dc332p(j,k)+wp(j,k)*dc333p(j,k))
     
                trp(j,k) = c11p(j,k) + c22p(j,k) + c33p(j,k)
     
                if (ipeter .eq. 0) then  ! peterlin function is 1.0
    
                    str11np(j,k)= c11p(j,k)
                    str12np(j,k)= c12p(j,k)
                    str13np(j,k)= c13p(j,k)
                    str22np(j,k)= c22p(j,k)
                    str23np(j,k)= c23p(j,k)
                    str33np(j,k)= c33p(j,k)
              
                else
               
                    str11np(j,k)=((zlmax**2 - 3.0)/(zlmax**2 - trp(j,k)))*c11p(j,k)
                    str12np(j,k)=((zlmax**2 - 3.0)/(zlmax**2 - trp(j,k)))*c12p(j,k)
                    str13np(j,k)=((zlmax**2 - 3.0)/(zlmax**2 - trp(j,k)))*c13p(j,k)
                    str22np(j,k)=((zlmax**2 - 3.0)/(zlmax**2 - trp(j,k)))*c22p(j,k)
                    str23np(j,k)=((zlmax**2 - 3.0)/(zlmax**2 - trp(j,k)))*c23p(j,k)
                    str33np(j,k)=((zlmax**2 - 3.0)/(zlmax**2 - trp(j,k)))*c33p(j,k)
               
                end if
               
                ! Add brownian motion terms
                str11np(j,k)=str11np(j,k)-1.0
                str22np(j,k)=str22np(j,k)-1.0
                str33np(j,k)=str33np(j,k)-1.0

                ! Polymer model
                if (itarget .eq. 0) then ! polymer ocean
                    beta_poly(j,k) = qbeta
                else if (itarget .eq. 1) then
                    ! Nonlinear model:
                    ! beta = exp(-alpha*gamma) --> beta_poly
                    ! alpha = 2.6e-03 PPM --> alpha_poly
                    ! gamma = scalar concentration (PPM) --> scp

                    if (alpha_poly*abs(scp(j,k)) .lt. 18.0) then
                        beta_poly(j,k) = exp(-alpha_poly*abs(scp(j,k)))
                    end if

                else if (itarget .eq. 2) then ! Linear polymer model
                    ! Linear model:
                    ! beta = (alpha*|gamma|)
                    beta_poly(j,k) = 1.0/(alpha_poly*abs(scp(j,k)) + 1.0)
                end if

                zbeta1 = (1.0 - beta_poly(j,k))/(re*beta_poly(j,k)*tpoly) ! = (nu_0 - nu_s)

                qp11np(j,k) = zbeta1*str11np(j,k)
                qp12np(j,k) = zbeta1*str12np(j,k)
                qp13np(j,k) = zbeta1*str13np(j,k)
                qp22np(j,k) = zbeta1*str22np(j,k)
                qp23np(j,k) = zbeta1*str23np(j,k)
                qp33np(j,k) = zbeta1*str33np(j,k)

                else
                    beta_poly(j,k) = 1.0 ! for printing
                end if ! it >= src_start-1
#ENDIF
                cflcheck = (abs(up(j,k))/delxm + abs(vp(j,k))/seght(i) + abs(wp(j,k))/delzm)*dt
                if(cflcheck .gt. cfl(i)) cfl(i) = cflcheck

                !---------------------------------------------------------------------!
                ! Compute immersed boundary force terms and add to nonlinear term     !
                ! The type of forcing is based on the value of imatrix at the grid    !
                ! point. Several of these cases probably require more coding and      !
                ! specification in the geometry file, but they have been unused since !
                ! I've used the code                                                  !
                !---------------------------------------------------------------------!
    
                if (imatrix(i,j,k) .eq. 1 .or. (imatrix(i,j,k) .eq. 4 .and. it .le. perturbtime)) then ! Generate solid surfaces
                    fxintg(i,j,k) = fxintg(i,j,k) + up(j,k)
                    fyintg(i,j,k) = fyintg(i,j,k) + vp(j,k)
                    fzintg(i,j,k) = fzintg(i,j,k) + wp(j,k)
                    dragx(j,k) = (-ugain*up(j,k)) - gain*fxintg(i,j,k)
                    dragy(j,k) = (-ugain*vp(j,k)) - gain*fyintg(i,j,k)
                    dragz(j,k) = (-ugain*wp(j,k)) - gain*fzintg(i,j,k)
        
                else if (imatrix(i,j,k) .eq. 3) then ! Spanwise-damping textures
                    fzintg(i,j,k) = fzintg(i,j,k) + wp(j,k)
                    dragz(j,k) = -ugain*wp(j,k) - gain*fzintg(i,j,k)
        
                else if (imatrix(i,j,k) .eq. 8) then ! Spanwise moving wall
                    if (k .le. (bfwidth*3 + 1)) then
                        wdes = (0.3*Uinf)/slopelength*(k - (bfwidth*3 + 1 - slopelength))
                    else
                        wdes = 0.3*Uinf
                    end if
    
                    fxintg(i,j,k) = fxintg(i,j,k) + up(j,k)
                    fyintg(i,j,k) = fyintg(i,j,k) + vp(j,k)
                    fzintg(i,j,k) = fzintg(i,j,k) + wp(j,k)-wdes
                    dragx(j,k)    = -ugain*up(j,k) - gain*fxintg(i,j,k)
                    dragy(j,k)    = -ugain*vp(j,k) - gain*fyintg(i,j,k)
                    dragz(j,k)    = -ugain*(wp(j,k)-wdes) - gain*fzintg(i,j,k)
                    
                else if (imatrix(i,j,k) .eq. 6) then ! Buffer zone
                    ! NOTE: initu and initw should technically be transformed first, but since
                    ! the buffer region has only been used for x- and z-constant flows, using
                    ! them as they are should be fine
                    fxintg(i,j,k) = fxintg(i,j,k) + up(j,k) - initu(i,j,k)
                    fyintg(i,j,k) = fyintg(i,j,k) + vp(j,k) - initv(i,j,k)
                    fzintg(i,j,k) = fzintg(i,j,k) + wp(j,k) - initw(i,j,k)
                    vwx(j,k) = vwx(j,k) + (-bfugain(k-bfhead+1)*(up(j,k) - initu(i,j,k)) - (bfgain(k-bfhead+1) * fxintg(i,j,k)))
                    vwy(j,k) = vwy(j,k) + (-bfugain(k-bfhead+1)*(vp(j,k) - initv(i,j,k)) - (bfgain(k-bfhead+1) * fyintg(i,j,k)))
                    vwz(j,k) = vwz(j,k) + (-bfugain(k-bfhead+1)*(wp(j,k) - initw(i,j,k)) - (bfgain(k-bfhead+1) * fzintg(i,j,k)))
    
                else if (imatrix(i,j,k) .eq. 7) then ! Suction region
                    fyintg(i,j,k) = fyintg(i,j,k) + (vp(j,k)-vdes(i))
                    vwy(j,k) = vwy(j,k) + (-gain*fyintg(i,j,k)) + (-ugain*(vp(j,k)-vdes(i)))
                    fxintg(i,j,k) = fxintg(i,j,k) + (up(j,k)-uinf)
                    vwx(j,k) = vwx(j,k) + (-gain*fxintg(i,j,k)) + (-ugain*(up(j,k)-uinf))
                    fzintg(i,j,k) = fzintg(i,j,k) + (wp(j,k)-initw(i,j,k))
                    vwz(j,k) = vwz(j,k) + (-gain*fzintg(i,j,k)) + (-ugain*(wp(j,k)-initw(i,j,k)))
     
                else if (imatrix(i,j,k) .eq. 11) then ! Constant x-force between wall and suction layer
                    dragx(j,k) = 1.0
    
                else if (imatrix(i,j,k) .eq. 0 .and. flow_select .eq. 1) then ! Constant pressure gradient
                    vwx(j,k) = vwx(j,k) - dPdx 
                
                else if (flow_select .eq. 5 .and. it .le. perturbtime) then ! Vortex ring
                    zj = float(j-1)*delzm
                    rsq = (ycoord(i) - ycenter)**2 + (zj - zcenter)**2
                    argrad = forbeta*(rad - sqrt(rsq))
                    fr = 0.5*(1.0 + tanh(argrad))
    
                    xi = float(k-1)*delxm
                    argx = forbeta*(L - abs(xi - xcenter))
                    fx = 0.5*(1.0 + tanh(argx))
  
                    vwx(j,k) = vwx(j,k) + bdyfx*fx*fr
    
                end if
#IFDEF SCALAR
                if (scl_flag .ge. 2) then
                    vc(j,k) = vc(j,k) + scsource(i,j,k)/dt
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
                            xsegdrag = segdrag*dragx(j,k)
                            ysegdrag = segdrag*dragy(j,k)
                            zsegdrag = segdrag*dragz(j,k)
                            vwx(jpjj,ipii) = vwx(jpjj,ipii) + xsegdrag
                            vwy(jpjj,ipii) = vwy(jpjj,ipii) + ysegdrag
                            vwz(jpjj,ipii) = vwz(jpjj,ipii) + zsegdrag
                        else if (j .eq. 1 .or. j .eq. 2) then 
                            jpjj = j + jj
                            if(jpjj .lt. 1 ) jpjj = jpjj + mz
                            xsegdrag = segdrag*dragx(j,k)
                            ysegdrag = segdrag*dragy(j,k)
                            zsegdrag = segdrag*dragz(j,k)
                            vwx(jpjj,ipii) = vwx(jpjj,ipii) + xsegdrag
                            vwy(jpjj,ipii) = vwy(jpjj,ipii) + ysegdrag
                            vwz(jpjj,ipii) = vwz(jpjj,ipii) + zsegdrag
                        else 
                            jpjj = j + jj
                            if(jpjj .gt. mz) jpjj = jpjj - mz
                            xsegdrag = segdrag*dragx(j,k)
                            ysegdrag = segdrag*dragy(j,k)
                            zsegdrag = segdrag*dragz(j,k)
                            vwx(jpjj,ipii) = vwx(jpjj,ipii) + xsegdrag
                            vwy(jpjj,ipii) = vwy(jpjj,ipii) + ysegdrag
                            vwz(jpjj,ipii) = vwz(jpjj,ipii) + zsegdrag
                        end if
                    end do ! jj
                end do ! ii
                end if
            end do ! k
        end do ! j

        ! Save certain variables for printing and/or particle tracking in 3D physical space
        up3d(i,:,:) = up
        vp3d(i,:,:) = vp
        wp3d(i,:,:) = wp
        wx3d(i,:,:) = wxp
        wy3d(i,:,:) = wyp
        wz3d(i,:,:) = wzp
        u11p3d(i,:,:) = u11p
        u12p3d(i,:,:) = u12p
        u13p3d(i,:,:) = u13p
        u21p3d(i,:,:) = u21p
        u22p3d(i,:,:) = u22p
        u23p3d(i,:,:) = u23p
        u31p3d(i,:,:) = u31p
        u32p3d(i,:,:) = u32p
        u33p3d(i,:,:) = u33p
    
        Lup3d(i,:,:) = Lup
        Lvp3d(i,:,:) = Lvp
        Lwp3d(i,:,:) = Lwp
#IFDEF SCALAR   
        scp3d(i,:,:) = scp
#ENDIF
#IFDEF POLYMER
        trC(i,:,:) = trp
        beta3d(i,:,:) = beta_poly
        p12(i,:,:) = qp12np
#ENDIF
        
        !---------------------------------------------------------------------!
        !              Now transform vxw back to spectral space               !
        !---------------------------------------------------------------------!
    

        ! Real --> Complex x-transform
        call fftw_execute_dft_r2c(planXf,vwx,vwxs)
        call fftw_execute_dft_r2c(planXf,vwy,vwys)
        call fftw_execute_dft_r2c(planXf,vwz,vwzs)

#IFDEF SCALAR
        call fftw_execute_dft_r2c(planXf,vc,vcs)
#ENDIF
#IFDEF POLYMER
        if (it .ge. (src_start-1)) then
        call fftw_execute_dft_r2c(planXf,c11np,c11ns)
        call fftw_execute_dft_r2c(planXf,c12np,c12ns)
        call fftw_execute_dft_r2c(planXf,c13np,c13ns)
        call fftw_execute_dft_r2c(planXf,c22np,c22ns)
        call fftw_execute_dft_r2c(planXf,c23np,c23ns)
        call fftw_execute_dft_r2c(planXf,c33np,c33ns)

        call fftw_execute_dft_r2c(planXf,str11np,str11ns)
        call fftw_execute_dft_r2c(planXf,str12np,str12ns)
        call fftw_execute_dft_r2c(planXf,str13np,str13ns)
        call fftw_execute_dft_r2c(planXf,str22np,str22ns)
        call fftw_execute_dft_r2c(planXf,str23np,str23ns)
        call fftw_execute_dft_r2c(planXf,str33np,str33ns)

        call fftw_execute_dft_r2c(planXf,qp11np,qp11s)
        call fftw_execute_dft_r2c(planXf,qp12np,qp12s)
        call fftw_execute_dft_r2c(planXf,qp13np,qp13s)
        call fftw_execute_dft_r2c(planXf,qp22np,qp22s)
        call fftw_execute_dft_r2c(planXf,qp23np,qp23s)
        call fftw_execute_dft_r2c(planXf,qp33np,qp33s)
        end if
#ENDIF

        ! Complex --> Complex z-transform
        call fftw_execute_dft(planZf,vwxs,vwxs)
        call fftw_execute_dft(planZf,vwys,vwys)
        call fftw_execute_dft(planZf,vwzs,vwzs)
        
#IFDEF SCALAR
        call fftw_execute_dft(planZf,vcs,vcs)
#ENDIF
#IFDEF POLYMER
        if (it .ge. (src_start-1)) then
        call fftw_execute_dft(planZf,c11ns,c11ns)
        call fftw_execute_dft(planZf,c12ns,c12ns)
        call fftw_execute_dft(planZf,c13ns,c13ns)
        call fftw_execute_dft(planZf,c22ns,c22ns)
        call fftw_execute_dft(planZf,c23ns,c23ns)
        call fftw_execute_dft(planZf,c33ns,c33ns)

        call fftw_execute_dft(planZf,str11ns,str11ns)
        call fftw_execute_dft(planZf,str12ns,str12ns)
        call fftw_execute_dft(planZf,str13ns,str13ns)
        call fftw_execute_dft(planZf,str22ns,str22ns)
        call fftw_execute_dft(planZf,str23ns,str23ns)
        call fftw_execute_dft(planZf,str33ns,str33ns)

        call fftw_execute_dft(planZf,qp11s,qp11s)
        call fftw_execute_dft(planZf,qp12s,qp12s)
        call fftw_execute_dft(planZf,qp13s,qp13s)
        call fftw_execute_dft(planZf,qp22s,qp22s)
        call fftw_execute_dft(planZf,qp23s,qp23s)
        call fftw_execute_dft(planZf,qp33s,qp33s)
        end if
#ENDIF

        ! Fill out complex spectral arrays
        do k = 1,nxh
            do j = 1,nz
                if (j .le. nzh) jj = j
                if (j .gt. nzh) jj = (mz-nz) + j
                gnr(i,j,k) = real(vwxs(jj,k))/float(mx*mz)
                gni(i,j,k) = aimag(vwxs(jj,k))/float(mx*mz)
                fnr(i,j,k) = real(vwys(jj,k))/float(mx*mz)
                fni(i,j,k) = aimag(vwys(jj,k))/float(mx*mz)
                omzr(i,j,k) = real(vwzs(jj,k))/float(mx*mz)
                omzi(i,j,k) = aimag(vwzs(jj,k))/float(mx*mz)
#IFDEF SCALAR
                scnr(i,j,k) = real(vcs(jj,k))/float(mx*mz)
                scni(i,j,k) = aimag(vcs(jj,k))/float(mx*mz)
#ENDIF
#IFDEF POLYMER
                if (it .ge. src_start-1) then
                c11nr(i,j,k)   = real(c11ns(jj,k))/float(mx*mz)
                c12nr(i,j,k)   = real(c12ns(jj,k))/float(mx*mz)
                c13nr(i,j,k)   = real(c13ns(jj,k))/float(mx*mz)
                c22nr(i,j,k)   = real(c22ns(jj,k))/float(mx*mz)
                c23nr(i,j,k)   = real(c23ns(jj,k))/float(mx*mz)
                c33nr(i,j,k)   = real(c33ns(jj,k))/float(mx*mz)

                str11nr(i,j,k) = real(str11ns(jj,k))/float(mx*mz)
                str12nr(i,j,k) = real(str12ns(jj,k))/float(mx*mz)
                str13nr(i,j,k) = real(str13ns(jj,k))/float(mx*mz)
                str22nr(i,j,k) = real(str22ns(jj,k))/float(mx*mz)
                str23nr(i,j,k) = real(str23ns(jj,k))/float(mx*mz)
                str33nr(i,j,k) = real(str33ns(jj,k))/float(mx*mz)

                qp11r(i,j,k)   = real(qp11s(jj,k))/float(mx*mz)
                qp12r(i,j,k)   = real(qp12s(jj,k))/float(mx*mz)
                qp13r(i,j,k)   = real(qp13s(jj,k))/float(mx*mz)
                qp22r(i,j,k)   = real(qp22s(jj,k))/float(mx*mz)
                qp23r(i,j,k)   = real(qp23s(jj,k))/float(mx*mz)
                qp33r(i,j,k)   = real(qp33s(jj,k))/float(mx*mz)

                c11ni(i,j,k)   = aimag(c11ns(jj,k))/float(mx*mz)
                c12ni(i,j,k)   = aimag(c12ns(jj,k))/float(mx*mz)
                c13ni(i,j,k)   = aimag(c13ns(jj,k))/float(mx*mz)
                c22ni(i,j,k)   = aimag(c22ns(jj,k))/float(mx*mz)
                c23ni(i,j,k)   = aimag(c23ns(jj,k))/float(mx*mz)
                c33ni(i,j,k)   = aimag(c33ns(jj,k))/float(mx*mz)

                str11ni(i,j,k) = aimag(str11ns(jj,k))/float(mx*mz)
                str12ni(i,j,k) = aimag(str12ns(jj,k))/float(mx*mz)
                str13ni(i,j,k) = aimag(str13ns(jj,k))/float(mx*mz)
                str22ni(i,j,k) = aimag(str22ns(jj,k))/float(mx*mz)
                str23ni(i,j,k) = aimag(str23ns(jj,k))/float(mx*mz)
                str33ni(i,j,k) = aimag(str33ns(jj,k))/float(mx*mz)

                qp11i(i,j,k)   = aimag(qp11s(jj,k))/float(mx*mz)
                qp12i(i,j,k)   = aimag(qp12s(jj,k))/float(mx*mz)
                qp13i(i,j,k)   = aimag(qp13s(jj,k))/float(mx*mz)
                qp22i(i,j,k)   = aimag(qp22s(jj,k))/float(mx*mz)
                qp23i(i,j,k)   = aimag(qp23s(jj,k))/float(mx*mz)
                qp33i(i,j,k)   = aimag(qp33s(jj,k))/float(mx*mz)
                end if
#ENDIF
            end do
        end do

    end do ! y-planes
    !$omp end parallel do simd

    ! Real --> Real DCT-I
    !$omp parallel do simd default(shared) private(j,k) collapse(2) schedule(simd: auto)
    do k = 1,nxh
        do j = 1,nz
            call fftw_execute_r2r(planY,fnr(:,j,k),fnr(:,j,k))
            call fftw_execute_r2r(planY,gnr(:,j,k),gnr(:,j,k))
            call fftw_execute_r2r(planY,omzr(:,j,k),omzr(:,j,k))
            call fftw_execute_r2r(planY,fni(:,j,k),fni(:,j,k))
            call fftw_execute_r2r(planY,gni(:,j,k),gni(:,j,k))
            call fftw_execute_r2r(planY,omzi(:,j,k),omzi(:,j,k))

            fn(:,j,k) = cmplx(fnr(:,j,k),fni(:,j,k))/float(ny)
            fn(1,j,k) = fn(1,j,k)/2.0
            fn(nyp,j,k) = fn(nyp,j,k)/2.0
            gn(:,j,k) = cmplx(gnr(:,j,k),gni(:,j,k))/float(ny)
            gn(1,j,k) = gn(1,j,k)/2.0
            gn(nyp,j,k) = gn(nyp,j,k)/2.0
            omz(:,j,k) = cmplx(omzr(:,j,k),omzi(:,j,k))/float(ny)
            omz(1,j,k) = omz(1,j,k)/2.0
            omz(nyp,j,k) = omz(nyp,j,k)/2.0

#IFDEF SCALAR
            call fftw_execute_r2r(planY,scnr(:,j,k),scnr(:,j,k))
            call fftw_execute_r2r(planY,scni(:,j,k),scni(:,j,k))

            scn(:,j,k) = cmplx(scnr(:,j,k),scni(:,j,k))/float(ny)
            scn(1,j,k) = scn(1,j,k)/2.0
            scn(nyp,j,k) = scn(nyp,j,k)/2.0
#ENDIF
#IFDEF POLYMER
            call fftw_execute_r2r(planY,c11nr(:,j,k),c11nr(:,j,k))
            call fftw_execute_r2r(planY,c12nr(:,j,k),c12nr(:,j,k))
            call fftw_execute_r2r(planY,c13nr(:,j,k),c13nr(:,j,k))
            call fftw_execute_r2r(planY,c22nr(:,j,k),c22nr(:,j,k))
            call fftw_execute_r2r(planY,c23nr(:,j,k),c23nr(:,j,k))
            call fftw_execute_r2r(planY,c33nr(:,j,k),c33nr(:,j,k))

            call fftw_execute_r2r(planY,str11nr(:,j,k),str11nr(:,j,k))
            call fftw_execute_r2r(planY,str12nr(:,j,k),str12nr(:,j,k))
            call fftw_execute_r2r(planY,str13nr(:,j,k),str13nr(:,j,k))
            call fftw_execute_r2r(planY,str22nr(:,j,k),str22nr(:,j,k))
            call fftw_execute_r2r(planY,str23nr(:,j,k),str23nr(:,j,k))
            call fftw_execute_r2r(planY,str33nr(:,j,k),str33nr(:,j,k))

            call fftw_execute_r2r(planY,qp11r(:,j,k),qp11r(:,j,k))
            call fftw_execute_r2r(planY,qp12r(:,j,k),qp12r(:,j,k))
            call fftw_execute_r2r(planY,qp13r(:,j,k),qp13r(:,j,k))
            call fftw_execute_r2r(planY,qp22r(:,j,k),qp22r(:,j,k))
            call fftw_execute_r2r(planY,qp23r(:,j,k),qp23r(:,j,k))
            call fftw_execute_r2r(planY,qp33r(:,j,k),qp33r(:,j,k))

            call fftw_execute_r2r(planY,c11ni(:,j,k),c11ni(:,j,k))
            call fftw_execute_r2r(planY,c12ni(:,j,k),c12ni(:,j,k))
            call fftw_execute_r2r(planY,c13ni(:,j,k),c13ni(:,j,k))
            call fftw_execute_r2r(planY,c22ni(:,j,k),c22ni(:,j,k))
            call fftw_execute_r2r(planY,c23ni(:,j,k),c23ni(:,j,k))
            call fftw_execute_r2r(planY,c33ni(:,j,k),c33ni(:,j,k))

            call fftw_execute_r2r(planY,str11ni(:,j,k),str11ni(:,j,k))
            call fftw_execute_r2r(planY,str12ni(:,j,k),str12ni(:,j,k))
            call fftw_execute_r2r(planY,str13ni(:,j,k),str13ni(:,j,k))
            call fftw_execute_r2r(planY,str22ni(:,j,k),str22ni(:,j,k))
            call fftw_execute_r2r(planY,str23ni(:,j,k),str23ni(:,j,k))
            call fftw_execute_r2r(planY,str33ni(:,j,k),str33ni(:,j,k))

            call fftw_execute_r2r(planY,qp11i(:,j,k),qp11i(:,j,k))
            call fftw_execute_r2r(planY,qp12i(:,j,k),qp12i(:,j,k))
            call fftw_execute_r2r(planY,qp13i(:,j,k),qp13i(:,j,k))
            call fftw_execute_r2r(planY,qp22i(:,j,k),qp22i(:,j,k))
            call fftw_execute_r2r(planY,qp23i(:,j,k),qp23i(:,j,k))
            call fftw_execute_r2r(planY,qp33i(:,j,k),qp33i(:,j,k))

            c11n(:,j,k) = cmplx(c11nr(:,j,k),c11ni(:,j,k))/float(ny)
            c12n(:,j,k) = cmplx(c12nr(:,j,k),c12ni(:,j,k))/float(ny)
            c13n(:,j,k) = cmplx(c13nr(:,j,k),c13ni(:,j,k))/float(ny)
            c22n(:,j,k) = cmplx(c22nr(:,j,k),c22ni(:,j,k))/float(ny)
            c23n(:,j,k) = cmplx(c23nr(:,j,k),c23ni(:,j,k))/float(ny)
            c33n(:,j,k) = cmplx(c33nr(:,j,k),c33ni(:,j,k))/float(ny)

            str11n(:,j,k) = cmplx(str11nr(:,j,k),str11ni(:,j,k))/float(ny)
            str12n(:,j,k) = cmplx(str12nr(:,j,k),str12ni(:,j,k))/float(ny)
            str13n(:,j,k) = cmplx(str13nr(:,j,k),str13ni(:,j,k))/float(ny)
            str22n(:,j,k) = cmplx(str22nr(:,j,k),str22ni(:,j,k))/float(ny)
            str23n(:,j,k) = cmplx(str23nr(:,j,k),str23ni(:,j,k))/float(ny)
            str33n(:,j,k) = cmplx(str33nr(:,j,k),str33ni(:,j,k))/float(ny)

            qp11(:,j,k) = cmplx(qp11r(:,j,k),qp11i(:,j,k))/float(ny)
            qp12(:,j,k) = cmplx(qp12r(:,j,k),qp12i(:,j,k))/float(ny)
            qp13(:,j,k) = cmplx(qp13r(:,j,k),qp13i(:,j,k))/float(ny)
            qp22(:,j,k) = cmplx(qp22r(:,j,k),qp22i(:,j,k))/float(ny)
            qp23(:,j,k) = cmplx(qp23r(:,j,k),qp23i(:,j,k))/float(ny)
            qp33(:,j,k) = cmplx(qp33r(:,j,k),qp33i(:,j,k))/float(ny)

            c11n(1,j,k) = c11n(1,j,k)/2.0
            c12n(1,j,k) = c12n(1,j,k)/2.0
            c13n(1,j,k) = c13n(1,j,k)/2.0
            c22n(1,j,k) = c22n(1,j,k)/2.0
            c23n(1,j,k) = c23n(1,j,k)/2.0
            c33n(1,j,k) = c33n(1,j,k)/2.0

            str11n(1,j,k) = str11n(1,j,k)/2.0
            str12n(1,j,k) = str12n(1,j,k)/2.0
            str13n(1,j,k) = str13n(1,j,k)/2.0
            str22n(1,j,k) = str22n(1,j,k)/2.0
            str23n(1,j,k) = str23n(1,j,k)/2.0
            str33n(1,j,k) = str33n(1,j,k)/2.0

            qp11(1,j,k) = qp11(1,j,k)/2.0
            qp12(1,j,k) = qp12(1,j,k)/2.0
            qp13(1,j,k) = qp13(1,j,k)/2.0
            qp22(1,j,k) = qp22(1,j,k)/2.0
            qp23(1,j,k) = qp23(1,j,k)/2.0
            qp33(1,j,k) = qp33(1,j,k)/2.0

            c11n(nyp,j,k) = c11n(nyp,j,k)/2.0
            c12n(nyp,j,k) = c12n(nyp,j,k)/2.0
            c13n(nyp,j,k) = c13n(nyp,j,k)/2.0
            c22n(nyp,j,k) = c22n(nyp,j,k)/2.0
            c23n(nyp,j,k) = c23n(nyp,j,k)/2.0
            c33n(nyp,j,k) = c33n(nyp,j,k)/2.0

            str11n(nyp,j,k) = str11n(nyp,j,k)/2.0
            str12n(nyp,j,k) = str12n(nyp,j,k)/2.0
            str13n(nyp,j,k) = str13n(nyp,j,k)/2.0
            str22n(nyp,j,k) = str22n(nyp,j,k)/2.0
            str23n(nyp,j,k) = str23n(nyp,j,k)/2.0
            str33n(nyp,j,k) = str33n(nyp,j,k)/2.0

            qp11(nyp,j,k) = qp11(nyp,j,k)/2.0
            qp12(nyp,j,k) = qp12(nyp,j,k)/2.0
            qp13(nyp,j,k) = qp13(nyp,j,k)/2.0
            qp22(nyp,j,k) = qp22(nyp,j,k)/2.0
            qp23(nyp,j,k) = qp23(nyp,j,k)/2.0
            qp33(nyp,j,k) = qp33(nyp,j,k)/2.0

#ENDIF            
        end do
    end do
    !$omp end parallel do simd

    !---------------------------------------------------------------------!
    !     Calculate swirl criterion and write data for visualization      !
    !---------------------------------------------------------------------!
    call calcswirl(u11p3d,u21p3d,u31p3d,u12p3d,u22p3d,u32p3d,u13p3d,u23p3d,u33p3d,swirl_3d)
!    call calcQ(u11p3d,u21p3d,u31p3d,u12p3d,u22p3d,u32p3d,u13p3d,u23p3d,u33p3d,swirl_3d)

    ! Process 3D variables (write outputs in physical space)
    if (print3d .ne. 0) then
    
        if ((mod(it,iprnfrq) .eq. 0 .and. it .ne. 0) .or. it .eq. 1) then

            print *,'Writing output data...'            
            ! Write output files
            if (print3d .eq. 1) then ! Write output in ASCII format
#IFDEF POLYMER
                call write_flowfield_ascii(up3d,vp3d,wp3d,wx3d,wy3d,wz3d,swirl_3d,beta3d)
#ELIF DEFINED SCALAR
                call write_flowfield_ascii(up3d,vp3d,wp3d,wx3d,wy3d,wz3d,swirl_3d,scp3d)
#ELSE
                call write_flowfield_ascii(up3d,vp3d,wp3d,wx3d,wy3d,wz3d,swirl_3d)
#ENDIF
#IFDEF OUTPUTFORM
            else if (print3d .eq. 3) then ! Write output in Tecplot binary (.szplt)
#IFDEF POLYMER
                call write_flowfield_plt(up3d,vp3d,wp3d,wx3d,wy3d,wz3d,swirl_3d,beta3d)
#ELIF DEFINED SCALAR                                         
                call write_flowfield_plt(up3d,vp3d,wp3d,wx3d,wy3d,wz3d,swirl_3d,scp3d)
#ELSE                                                        
                call write_flowfield_plt(up3d,vp3d,wp3d,wx3d,wy3d,wz3d,swirl_3d)
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
        call part_track(up3d,vp3d,wp3d,wx3d,wy3d,wz3d,u_old,v_old,w_old,  &
                        u11p3d,u12p3d,u13p3d,u21p3d,u22p3d,u23p3d,u31p3d, &
                        u32p3d,u33p3d,Lup3d,Lvp3d,Lwp3d,Lunm1,Lvnm1,Lwnm1)
    
        ! Save velocity and Laplacian for time derivatives inside particle integration
        u_old = up3d
        v_old = vp3d
        w_old = wp3d
    
        Lunm1 = Lup3d
        Lvnm1 = Lvp3d
        Lwnm1 = Lwp3d
    end if
   
    ! Calculate Q-criterion at particle locations
    do n = 1,npart
        call fluid_interp1(xpart(n),ypart(n),zpart(n),swirl_3d,swirl_part(n))
    end do
 
    !---------------------------------------------------------------------!
    ! Write Mean U Data, calculate output data and check polymer addition !
    !---------------------------------------------------------------------!
  
    if (flow_select .eq. 1 .or. flow_select .eq. 4) then ! Only relevant for wall-bounded turbulence
        ! Mean velocity
        do i = 1,nyp
            do k = 1,mx
                uzmean(k) = sum(up3d(i,:,k))/mz
            end do
            uxmean(i) = sum(uzmean)/mx
        end do
    end if
  

    !---------------------------------------------------------------------!
    !               Check for numerical instabilities                     !
    !---------------------------------------------------------------------!
    
    ! Check CFL condition 
    do i = 1,nyp
        if (cfl(i) > cflmax) cflmax = cfl(i)
    end do
   
    if (mod(it,cadence) .eq. 1) then 
    end if
    
    ! Check CFL condition
    if (cflmax .gt. 1.0) then
        write(*,*) 'CFL failure at time step ',it
        stop
    end if

    ! Calculate output data and write to outputs/<file>
    if (mod(it,cadence) .eq. 1) then ! only write every <cadence> time steps
        write(*,*) '    max CFL = ',cflmax
        write(*,*) '    Writing output data...'
        call writeoutputs(up3d,vp3d,wp3d,wx3d,wy3d,wz3d,u12p3d,uxmean,swirl_3d &
#IFDEF POLYMER
                          ,beta3d,p12,trC &
#ELIF DEFINED SCALAR
                          ,scp3d &
#ENDIF
                          )

        write(*,*) '    Done!'
    end if

#IFDEF POLYMER
    ! Check to see if we can stop adding polymer
    if (scl_flag .ge. 2 .and. it .le. src_stop .and. it .ge. src_start) then
        call calc_total_beta(it,delxm,delzm,scp3d,beta3d)
    end if
#ENDIF
    
    end subroutine vcw3dp
end module solvers


