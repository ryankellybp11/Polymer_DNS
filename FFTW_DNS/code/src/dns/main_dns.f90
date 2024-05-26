program dns

! ============================================================================ !
!                             Declare Modules                                  !
! ============================================================================ !
    use,intrinsic :: iso_c_binding
    use grid_size
    use omp_lib
    use solvers
    use derivs
    use helpers
! ---------------------------------------------------------------------------- !

! ============================================================================ !
!                      Declare all varibles explicitly                         !
! ============================================================================ !

    implicit none

    include 'fftw3.f03'

    ! FFTW Plan variables
    type(C_PTR) :: planZb,planXb,planY,planXf,planZf
    complex(C_DOUBLE_COMPLEX), dimension(nyp,mz,mx) :: aspec
    real(C_DOUBLE), dimension(nyp,mz,mx) :: aphys
    integer :: ierr

    ! Flow variables
    complex, dimension(nyp,nz,nxh) :: u,v,w,omx,omy,omz
    complex, dimension(nyp,nz,nxh) :: u11,u12,u13
    complex, dimension(nyp,nz,nxh) :: u21,u22,u23
    complex, dimension(nyp,nz,nxh) :: u31,u32,u33
    complex, dimension(nyp,nz,nxh) :: Lu,Lv,Lw
    real,    dimension(nyp,mz,mx)  :: Lu_old,Lv_old,Lw_old

#IFDEF SCALAR
    ! Scalar variables
    complex, dimension(nyp,nz,nxh) :: scalar,sclx,scly,sclz ! Concentration & derivatives
    complex, dimension(nyp,nz,nxh) :: scn,scnm1 ! Nonlinear terms
    complex, dimension(nz,nxh)     :: bcscltop = 0.0, bcsclbot = 0.0 ! Scalar BCs
    real :: atscl,btscl,abscl,bbscl

    real    :: sigmax,sigmay,sigmaz
    real    :: deltaT,diff
    integer :: scl_flag
    integer :: src_start,src_stop
#ENDIF
#IFDEF POLYMER
    ! Polymer variables
    complex, dimension(nyp,nz,nxh) :: c11,c12,c13,c21,c22,c23,c31,c32,c33
    complex, dimension(nyp,nz,nxh) :: c11n,c12n,c13n,c22n,c23n,c33n
    complex, dimension(nyp,nz,nxh) :: c11nm1,c12nm1,c13nm1,c22nm1,c23nm1,c33nm1
    complex, dimension(nyp,nz,nxh) :: c11NL,c12NL,c13NL,c22NL,c23NL,c33NL

    complex, dimension(nyp,nz,nxh) :: dc111,dc112,dc113,dc121,dc122,dc123,dc131,dc132,dc133
    complex, dimension(nyp,nz,nxh) :: dc211,dc212,dc213,dc221,dc222,dc223,dc231,dc232,dc233
    complex, dimension(nyp,nz,nxh) :: dc311,dc312,dc313,dc321,dc322,dc323,dc331,dc332,dc333

    complex, dimension(nyp,nz,nxh) :: str11n,str12n,str13n,str22n,str23n,str33n
    complex, dimension(nyp,nz,nxh) :: str11nm1,str12nm1,str13nm1,str22nm1,str23nm1,str33nm1
    complex, dimension(nyp,nz,nxh) :: qp11,qp12,qp13,qp22,qp23,qp33
    complex, dimension(nyp,nz,nxh) :: t1,t2,t3

    real :: alpha_poly,tpoly,zlmax,diffpoly,qbeta
    real :: atp,btp,abp,bbp
    real :: c11z,c22z,c33z
#ENDIF

    ! Flow setup variables
    real    :: re,Uinf,R_tau,dPdx

    ! Solver variables
    complex, dimension(nyp,nz,nxh) :: fn,fnm1,gn,gnm1,gf1,gf2
    complex, dimension(nz,nxh)     :: c1,c2,bctop,bcbot,dnv1b,dnv1t,dnv2b,dnv2t,denom
    real,    dimension(nyp)        :: u0,h1n,h1nm1
    real,    dimension(nyp)        :: w0,h3n,h3nm1
    real,    dimension(nxh)        :: wavx
    real,    dimension(nz)         :: wavz

    real :: gain,ugain,theta,alpha,beta,dyde

    ! BC variables
    integer :: nt,nb
    real    :: p1t,p1b,p2t,p2b
    real    :: atu,btu,gtu,abu,bbu,gbu

    ! Calculation variables
    complex, dimension(nyp,nz,nxh) :: wrkc, wrk1
#IFDEF POLYMER
    real,    dimension(nyp,nz,nx)  :: wrk11,wrk12,wrk13,wrk21,wrk22,wrk23,wrk31,wrk32,wrk33
#ENDIF
    real,    dimension(nyhp,nz,3)  :: a

    real    :: wn2,x,g,ysmth,zsmth
    integer :: i,j,k,ib,iyt,iyb,jj,uflag

    ! Simulation control variables
    integer :: irstrt,nsteps,iprnfrq,print3d,crstrt
    integer :: it
    real    :: dt

    ! Immersed boundary force variables
    real, dimension(nyp,mz,mx) :: fxintg,fyintg,fzintg
    real, dimension(mz2,mx2)   :: fspread

    ! temp
    real,dimension(nyp,nz,nx) :: initu,initv,initw
! ---------------------------------------------------------------------------- !


! ============================================================================ !
!                            Define common blocks                              !
! ============================================================================ !

    common/iocontrl/  irstrt,nsteps,iprnfrq,print3d,crstrt
    common/init/     initu,initv,initw
    common/waves/     wavx,wavz
    common/u0bcs/     atu,btu,gtu,abu,bbu,gbu
    common/grnfcn/    nt,nb,p1t,p1b,p2t,p2b
    common/solver/    gain,ugain,theta,alpha,beta,dyde
    common/itime/     it
    common/dtime/     dt
    common/ibforce/   fxintg,fyintg,fzintg,fspread
    common/flow/      re,Uinf,R_tau,dPdx
#IFDEF SCALAR
    common/scl_stuff/ sigmax,sigmay,sigmaz,deltaT,diff,scl_flag
    common/src_time/  src_start,src_stop
#ENDIF
#IFDEF POLYMER
    common/poly_var/  alpha_poly,tpoly,zlmax,diffpoly,qbeta
    common/c_init/    c11z,c22z,c33z
#ENDIF
! ---------------------------------------------------------------------------- !

! ============================================================================ !
!                          Initialization Routines                             !
! ============================================================================ !

    ! -------------------------------------------------------------------- !
    ! Initializes all preliminary variables before heading into the main   !
    ! time loop. Reads data from setup/dns.config, computes preliminary    !
    ! derivatives and Green's functions, sets up immersed boundary forces, !
    ! and uses a start-up routine to get a time step 0. The main loop uses !
    ! a second-order Adams-Bashforth integrator, so it needs the first     !
    ! time step in order to get start.                                     !
    ! -------------------------------------------------------------------- !
  
    call setstuff ! Reads setup file to initialize variables


    ! Plan FFT routines
    ierr = fftw_init_threads() ! Initialize threaded FFTW (OpenMP)
    if (ierr .eq. 0) then
        print *,'Error! FFTW could not initialize threads!'
        stop
    end if
    call fftw_plan_with_nthreads(OMP_GET_MAX_THREADS())

    ! NOTE:
    !   With the FFTW_PATIENT option, this process may take a few seconds
    !   but once the plans are created, they can be used indefinitely.
    !   The variables aspec and aphys are dummy variables that are
    !   overwritten during this process, so they can't be used after this.
    !   It also has the benefit of reducing the number of OpenMP threads
    !   used if it determines less is better.
    print *,'Planning FFT routines...'
    planZb = fftw_plan_dft_1d(mz,aspec,aspec,FFTW_BACKWARD,FFTW_PATIENT)
    planXb = fftw_plan_dft_c2r_1d(mx,aspec,aphys,FFTW_PATIENT)
    planY  = fftw_plan_r2r_1d(nyp,aphys,aphys,FFTW_REDFT00,FFTW_PATIENT)
    planXf = fftw_plan_dft_r2c_1d(mx,aphys,aspec,FFTW_PATIENT)
    planZf = fftw_plan_dft_1d(mz,aspec,aspec,FFTW_FORWARD,FFTW_PATIENT)

    print *,'Done!'
    
    ! Compute Green's functions for the two systems
    call gfcn(gf1,a,wrkc,wrk1,bctop,bcbot,dnv1b,dnv1t,nb,nt,p1b,p1t)
    call gfcn(gf2,a,wrkc,wrk1,bctop,bcbot,dnv2b,dnv2t,nb,nt,p2b,p2t)

    ! Calculate common denominator
    !$omp parallel do
    do k = 1,nxh
        do j = 1,nz
            denom(j,k) = dnv1t(j,k)*dnv2b(j,k) - dnv2t(j,k)*dnv1b(j,k)
        end do
    end do
    !$omp end parallel do

    call forcepre ! Initializes immersed boundary forces

    it = 0 ! Start time loop counter at 0

    ! Compute start-up condition to use Adams-Bashforth later
    call initial(u,u0,v,w,w0,omx,omy,omz,fn,fnm1,gn,gnm1,h1n,h1nm1,h3n,h3nm1, &
                 wrk1,wrkc,u11,u12,u13,u21,u22,u23,u31,u32,u33, &
#IFDEF SCALAR
                 scalar,sclx,scly,sclz,scn,scnm1,     &
#ENDIF
#IFDEF POLYMER
                 t1,t2,t3,c11,c12,c13,c21,c22,c23,c31,c32,c33,dc111,dc112,    &
                 dc113,dc121,dc122,dc123,dc131,dc132,dc133,dc211,dc212,dc213, &
                 dc221,dc222,dc223,dc231,dc232,dc233,dc311,dc312,dc313,dc321, &
                 dc322,dc323,dc331,dc332,dc333,c11n,c12n,c13n,c22n,c23n,c33n, &
                 c11nm1,c12nm1,c13nm1,c22nm1,c23nm1,c33nm1,str11n,str12n,     &
                 str13n,str22n,str23n,str33n,str11nm1,str12nm1,str13nm1,      &
                 str22nm1,str23nm1,str33nm1,qp11,qp12,qp13,qp22,qp23,qp33,    &
#ENDIF
                 Lu,Lv,Lw,Lu_old,Lv_old,Lw_old, &
                 planZb,planXb,planY,planXf,planZf) 
                 

! ---------------------------------------------------------------------------- !

! ============================================================================ !
!                                 Main Loop                                    !
! ============================================================================ !

    do it = 1,nsteps

        write(*,*) '   it = ',it

#IFDEF POLYMER
        ! -------------------------------------------------------------------- !
        ! Initialize conformation tensor on the time step before polymer is released
        if (it .eq. src_start - 1 .or. it .eq. irstrt .and. crstrt .eq. 0) then

            do k = 1,nx
                do j = 1,nz
                    do i = 1,nyp
                        wrk11(i,j,k) = c11z
                        wrk22(i,j,k) = c22z
                        wrk33(i,j,k) = c33z

                        wrk12(i,j,k) = 0.0
                        wrk13(i,j,k) = 0.0
                        wrk21(i,j,k) = 0.0
                        wrk23(i,j,k) = 0.0
                        wrk31(i,j,k) = 0.0
                        wrk32(i,j,k) = 0.0
                    end do
                end do
            end do

            call xyzfft(c11,wrk11,-1)
            call xyzfft(c12,wrk12,-1)
            call xyzfft(c13,wrk13,-1)
            call xyzfft(c21,wrk21,-1)
            call xyzfft(c22,wrk22,-1)
            call xyzfft(c23,wrk23,-1)
            call xyzfft(c31,wrk31,-1)
            call xyzfft(c32,wrk32,-1)
            call xyzfft(c33,wrk33,-1)
        end if

        ! -------------------------------------------------------------------- !
        ! Solve for conformation tensor with diffusion & Neumann BCs
        if (it .ge. src_start) then
            call polyNL(c11NL,c12NL,c13NL,c22NL,c23NL,c33NL,                   &
                        c11n,c12n,c13n,c22n,c23n,c33n,                         &
                        c11nm1,c12nm1,c13nm1,c22nm1,c23nm1,c33nm1,             &
                        str11n,str12n,str13n,str22n,str23n,str33n,             &
                        str11nm1,str12nm1,str13nm1,str22nm1,str23nm1,str33nm1)

            ! Note: The arguments are chosen such that I can still use phirhs for the C tensor solver
            !$omp parallel sections default(shared) private(wrkc,wrk1,bcbot,bctop)
            !$omp section
            call phirhs(c11,wrkc,wrk1,c11NL,2.0*c11NL,bcbot,bctop,diffpoly)
            !$omp section
            call phirhs(c12,wrkc,wrk1,c12NL,2.0*c12NL,bcbot,bctop,diffpoly)
            !$omp section
            call phirhs(c13,wrkc,wrk1,c13NL,2.0*c13NL,bcbot,bctop,diffpoly)
            !$omp section
            call phirhs(c22,wrkc,wrk1,c22NL,2.0*c22NL,bcbot,bctop,diffpoly)
            !$omp section
            call phirhs(c23,wrkc,wrk1,c23NL,2.0*c23NL,bcbot,bctop,diffpoly)
            !$omp section
            call phirhs(c33,wrkc,wrk1,c33NL,2.0*c33NL,bcbot,bctop,diffpoly)
            !$omp end parallel sections

            ! Copy n to nm1
            !$omp parallel do
            do k = 1,nxh
                do j = 1,nz
                    do i = 1,nyp
                        c11nm1(i,j,k) = c11n(i,j,k)
                        c12nm1(i,j,k) = c12n(i,j,k)
                        c13nm1(i,j,k) = c13n(i,j,k)
                        c22nm1(i,j,k) = c22n(i,j,k)
                        c23nm1(i,j,k) = c23n(i,j,k)
                        c33nm1(i,j,k) = c33n(i,j,k)

                        str11nm1(i,j,k) = str11n(i,j,k)
                        str12nm1(i,j,k) = str12n(i,j,k)
                        str13nm1(i,j,k) = str13n(i,j,k)
                        str22nm1(i,j,k) = str22n(i,j,k)
                        str23nm1(i,j,k) = str23n(i,j,k)
                        str33nm1(i,j,k) = str33n(i,j,k)
                    end do
                end do
            end do
            !$omp end parallel do

            g = (re*diffpoly)/(dt*theta)
            atp = 0.0
            btp = 1.0
            abp = 0.0
            bbp = 1.0

            do k = 1,nxh
                x = wavx(k)**2

                call penta(c11(1,1,k),x,g,dyde,ib,atp,btp,bctop(1,k),abp,bbp,bcbot(1,k),a,wrk1)
                call penta(c12(1,1,k),x,g,dyde,ib,atp,btp,bctop(1,k),abp,bbp,bcbot(1,k),a,wrk1)
                call penta(c13(1,1,k),x,g,dyde,ib,atp,btp,bctop(1,k),abp,bbp,bcbot(1,k),a,wrk1)
                call penta(c22(1,1,k),x,g,dyde,ib,atp,btp,bctop(1,k),abp,bbp,bcbot(1,k),a,wrk1)
                call penta(c23(1,1,k),x,g,dyde,ib,atp,btp,bctop(1,k),abp,bbp,bcbot(1,k),a,wrk1)
                call penta(c33(1,1,k),x,g,dyde,ib,atp,btp,bctop(1,k),abp,bbp,bcbot(1,k),a,wrk1)
            end do

            call norm(c11)
            call norm(c12)
            call norm(c13)
            call norm(c22)
            call norm(c23)
            call norm(c33)

            ! Use symmetry to get other components of C tensor
            !$omp parallel do
            do k = 1,nxh
                do j = 1,nz
                    do i = 1,nyp
                        c21(i,j,k) = c12(i,j,k)
                        c31(i,j,k) = c13(i,j,k)
                        c32(i,j,k) = c23(i,j,k)
                    end do
                end do
            end do
            !$omp end parallel do

        end if ! it .ge. src_start
#ENDIF
#IFDEF SCALAR
        ! -------------------------------------------------------------------- !
        ! Scalar calculated first
        call phirhs(scalar,wrkc,wrk1,scn,scnm1,bcsclbot,bcscltop,diff)

        ! Copy scn to scnm1
        !$omp parallel do
        do k = 1,nxh
            do j = 1,nz
                do i = 1,nyp
                    scnm1(i,j,k) = scn(i,j,k)
                end do
            end do
        end do
        !$omp end parallel do

        ! To use mixed BCs for scalar field, we need to use penta(...)
        ! Let's manually set zero flux conditions here for now. I don't know if/
        ! when we would ever want to change it, so I'll leave that for a later
        ! problem.
        x = 0.0
        g = (re*diff)/(dt*theta)
        ib = 0

        ! This should be atscl*scl + btscl*d(scl)/dy = bcscltop (for top, similar for bottom)
        atscl = 0.0
        btscl = 1.0
        abscl = 0.0
        bbscl = 1.0 

        do k = 1,nxh
            x = wavx(k)**2
            call penta(scalar(1,1,k),x,g,dyde,ib,atscl,btscl,bcscltop(1,k),abscl,bbscl,bcsclbot(1,k),a,wrk1)
        end do
#ENDIF

        ! -------------------------------------------------------------------- !
        ! Initial phi field is computed from the laplacian of the v field
        call cderiv(v,wrk1)
        call cderiv(wrk1,wrkc)

        !$omp parallel do
        do k = 1,nxh
            do j = 1,nz
                wn2 = wavx(k)**2 + wavz(j)**2
                do i = 1,nyp
                    v(i,j,k) = wrkc(i,j,k) - wn2*v(i,j,k)
                end do
            end do
        end do
        !$omp end parallel do


        ! Evaluate RHS of time discrete equation for phi
        call phirhs(v,wrkc,wrk1,fn,fnm1,bcbot,bctop,1.0)

        ! Copy fn to fnm1
        !$omp parallel do
        do k = 1,nxh
            do j = 1,nz
                do i = 1,nyp
                    fnm1(i,j,k) = fn(i,j,k)
                end do
            end do
        end do
        !$omp end parallel do

        ! -------------------------------------------------------------------- !
        g = re/(dt*theta) ! time step dependency parameter
        ib = 0

        do k = 1,nxh
            x = wavx(k)**2
            call solve(v(1,1,k),x,g,dyde,bctop(1,k),bcbot(1,k),ib,a,wrkc)
        end do

        ! Solve for vp (Laplace(vp) = phi)
        g = 0.0
        
        do k = 1,nxh
            x = wavx(k)**2
            call solve(v(1,1,k),x,g,dyde,bctop(1,k),bcbot(1,k),ib,a,wrkc)
        end do

        ! -------------------------------------------------------------------- !
        ! Evaluate RHS of time discrete equation for vorticity
        call phirhs(omy,wrkc,wrk1,gn,gnm1,bcbot,bctop,1.0)

        ! Copy gn to gnm1
        !$omp parallel do
        do k = 1,nxh
            do j = 1,nz
                do i = 1,nyp
                    gnm1(i,j,k) = gn(i,j,k)
                end do
            end do
        end do
        !$omp end parallel do

        ! ---------------------------------------------------------------------- !
        ! Evaluate derivatives necessary for constructing the composite solution !
        !                                                                        !
        ! v = vp + c1*gf1 + c2*gf2                                               !
        !                                                                        !
        ! evaluate d(nb)vp/dy(nb) at y = -1                                      !
        ! evaluate d(nt)vp/dy(nt) at y = +1                                      !
        ! ---------------------------------------------------------------------- !

        iyb = 2
        iyt = 4
        if (nb .eq. 1) then ! dv/dy
            call c1derbw(v,wrkc,iyb)
            if (nt .eq. 1) then
                call c1dertw(v,wrkc,iyt)
            else if (nt .eq. 2) then
                call cderiv(v,wrk1)
                call c1dertw(wrk1,wrkc,iyt)
            end if

        else if (nb .eq. 2) then ! d^2(v)/(dy)^2
            call cderiv(v,wrk1) 
            call c1derbw(wrk1,wrkc,iyb)
            if (nt .eq. 1) then
                call c1dertw(v,wrkc,iyt)
            else if (nt .eq. 2) then
                call cderiv(v,wrk1)
                call c1dertw(wrk1,wrkc,iyt)
            end if
        end if

        ! Compute c1 and c2
        do k = 1,nxh
            do j = 1,nz
                c1(j,k) = (wrkc(iyb,j,k)*dnv2t(j,k) - wrkc(iyt,j,k)*dnv2b(j,k))/denom(j,k)
                c2(j,k) = (wrkc(iyt,j,k)*dnv1b(j,k) - wrkc(iyb,j,k)*dnv1t(j,k))/denom(j,k)
            end do
        end do

        ! Form total v
        do k = 1,nxh
            do j = 1,nz
                do i = 1,nyp
                    v(i,j,k) = v(i,j,k) + c1(j,k)*gf1(i,j,k) + c2(j,k)*gf2(i,j,k)
                end do
            end do
        end do

        ! -------------------------------------------------------------------- !
        ! Now do vorticity field
        g = re/(dt*theta)

        do k = 1,nxh
            x = wavx(k)**2
            call penta(omy(1,1,k),x,g,dyde,ib,atu,btu,bctop(1,k),abu,bbu,bcbot(1,k),a,wrk1)
        end do

        ! Because of periodicity, the (0,0) mode for omy is always zero
        do i = 1,nyp
            omy(i,1,1) = (0.0,0.0)
        end do

        ! -------------------------------------------------------------------- !
        ! Solve for zeroth Fourier modes of streamwise and spanwise velocities
        uflag = 1
        call uzero(u0,h1n,h1nm1,uflag)
        uflag = 3
        call uzero(w0,h3n,h3nm1,uflag)

        do i = 1,nyp
            h1nm1(i) = h1n(i)
            h3nm1(i) = h3n(i)
        end do

        ! -------------------------------------------------------------------- !
        ! Calculate the remaining components of velocity from continuity and
        ! normal vorticity in spectral space:
        !
        ! U(n,kz,kx) = -((i*kx)(-dV/dy) + (i*kx)(omy))/(kx^2 + kz^2)
        ! if (kx^2 + kz^2) = 0, u0 above is the real part of u(i,kx=0)
        !
        ! W(n,kz,kx) = -((i*kz)(-dV/dy) - (i*kx)(omy))/(kx^2 + kz^2)
        ! If (kx^2 + kz^2) = 0, W(n,kz,kx) = 0.0
        !
        ! The imaginary part of the kx=kz=0 modes are used to store the real 
        ! part of U(i,kmax), W(i,kmax) which is set to 0 in the subroutine
        ! veloc(...) (? - original had "SUB. UVEL." which I'm guessing is veloc)
        !   For consistency with continuity, the subroutine veloc(...) also 
        !   resets the imaginary component of V(i,kx=0) (i.e., the real part 
        !   of the kmax mode) equal to 0
      
        call norm(v)
        call norm(omy)
        call veloc(u,v,w,u0,w0,omy,wrkc)
       
        ! Normalize to ensure conjugate symmetry
        call norm(u)
        call norm(w)

        ! Calculate the other vorticity components
        call vort(u,v,w,omx,omz)
        call norm(omx)
        call norm(omz)

#IFDEF SCALAR
        ! Calculate scalar gradient
        call norm(scalar)
        call gradscl(scalar,sclx,scly,sclz,wrkc)
        call norm(sclx)
        call norm(scly)
        call norm(sclz)
#ENDIF
 
        ! Calculate derivatives for velocity gradient tensor & Laplacian
        call gradv(u,v,w,u11,u12,u13,u21,u22,u23,u31,u32,u33,Lu,Lv,Lw)
        
        call norm(u11)
        call norm(u12)
        call norm(u13)
        call norm(u21)
        call norm(u22)
        call norm(u23)
        call norm(u31)
        call norm(u32)
        call norm(u33)
        
        call norm(Lu)
        call norm(Lv)
        call norm(Lw)

#IFDEF POLYMER
        ! Calculate derivatives of conformation tensor        
        call derivscji(c11,c12,c13,c21,c22,c23,c31,c32,c33,                   &
                       dc111,dc112,dc113,dc211,dc212,dc213,dc311,dc312,dc313, &
                       dc121,dc122,dc123,dc221,dc222,dc223,dc321,dc322,dc323, &
                       dc131,dc132,dc133,dc231,dc232,dc233,dc331,dc332,dc333)  
  
         call norm(c11)
         call norm(c12)
         call norm(c13)
         call norm(c21)
         call norm(c22)
         call norm(c23)
         call norm(c31)
         call norm(c32)
         call norm(c33)

         call norm(dc111)
         call norm(dc112)
         call norm(dc113)
         call norm(dc211)
         call norm(dc212)
         call norm(dc213)
         call norm(dc311)
         call norm(dc312)
         call norm(dc313)
         call norm(dc121)
         call norm(dc122)
         call norm(dc123)
         call norm(dc221)
         call norm(dc222)
         call norm(dc223)
         call norm(dc321)
         call norm(dc322)
         call norm(dc323)
         call norm(dc131)
         call norm(dc132)
         call norm(dc133)
         call norm(dc231)
         call norm(dc232)
         call norm(dc233)
         call norm(dc331)
         call norm(dc332) 
         call norm(dc333)                                                 
#ENDIF
       
        ! -------------------------------------------------------------------- !
        ! Compute v x omega in physical space - FFTs are performed inside subroutine
        ! Assumes all variables are in 3D spectral space (on input and output)
        call vcw3dp(u,v,w,omx,omy,omz,fn,gn,u11,u12,u13,u21,u22,u23,       &
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

        ! Normalizations
        call norm(v)
        call norm(omy)
        call norm(fn)
        call norm(gn)
        call norm(omz)
#IFDEF SCALAR
        call norm(scalar)
        call norm(scn)
#ENDIF
#IFDEF POLYMER
        if (it .ge. src_start-1) then
        call norm(c11) 
        call norm(c12)
        call norm(c13)
        call norm(c21)
        call norm(c22)
        call norm(c23)
        call norm(c31)
        call norm(c32)
        call norm(c33)

        call norm(c11n)
        call norm(c12n)
        call norm(c13n)
        call norm(c22n)
        call norm(c23n)
        call norm(c33n)

        call norm(str11n)
        call norm(str12n)
        call norm(str13n)
        call norm(str22n)
        call norm(str23n)
        call norm(str33n)

        call norm(qp11)
        call norm(qp12)
        call norm(qp13)
        call norm(qp22)
        call norm(qp23)
        call norm(qp33)
        end if

        ! -------------------------------------------------------------------- !
        ! Add polymer forces to the momentum equation
        ! qp is the polymer stress tensor, t = div(qp) which is the force on 
        ! the fluid due to polymer stress
        if (it .ge. src_start) then
            call polyforce(qp11,qp12,qp13,qp22,qp23,qp33,t1,t2,t3)
            call norm(t1)
            call norm(t2)
            call norm(t3)
            call subforce(gn,fn,omz,t1,t2,t3)


         if (it .eq. (src_start-1) .or. it .eq. irstrt .and. crstrt .eq. 0) then
         ! Set (n-1) terms equal to current terms for explicit time integration - Ryan 9/27/22
         !$omp parallel do
         do k=1,nxh
            do j=1,nz
                 do i=1,nyp
                    c11nm1(i,j,k)   = c11n(i,j,k)
                    c12nm1(i,j,k)   = c12n(i,j,k)
                    c13nm1(i,j,k)   = c13n(i,j,k)
                    c22nm1(i,j,k)   = c22n(i,j,k)
                    c23nm1(i,j,k)   = c23n(i,j,k)
                    c33nm1(i,j,k)   = c33n(i,j,k)
                    str11nm1(i,j,k) = str11n(i,j,k)
                    str12nm1(i,j,k) = str12n(i,j,k)
                    str13nm1(i,j,k) = str13n(i,j,k)
                    str22nm1(i,j,k) = str22n(i,j,k)
                    str23nm1(i,j,k) = str23n(i,j,k)
                    str33nm1(i,j,k) = str33n(i,j,k)
                 enddo
             enddo
         enddo
         !$omp end parallel do
         end if
        end if
#ENDIF

        ! -------------------------------------------------------------------- !
        ! Spectral smoothing, assuming force is in fn, gn, and omz
        
        ! y-direction smoothing
        !$omp parallel do
        do i = 1,nyp
            ysmth = exp(-3.0*((float(i-1)/float(ny))**10))
            do j = 1,nz
                do k = 1,nxh
                    omz(i,j,k) = omz(i,j,k)*ysmth
                    gn(i,j,k)  =  gn(i,j,k)*ysmth
                    fn(i,j,k)  =  fn(i,j,k)*ysmth
#IFDEF SCALAR
                    scn(i,j,k) = scn(i,j,k)*ysmth
#ENDIF
#IFDEF POLYMER
                    c11n(i,j,k) = c11n(i,j,k)*ysmth
                    c12n(i,j,k) = c12n(i,j,k)*ysmth
                    c13n(i,j,k) = c13n(i,j,k)*ysmth
                    c22n(i,j,k) = c22n(i,j,k)*ysmth
                    c23n(i,j,k) = c23n(i,j,k)*ysmth
                    c33n(i,j,k) = c33n(i,j,k)*ysmth
                    
                    str11n(i,j,k) = str11n(i,j,k)*ysmth
                    str12n(i,j,k) = str12n(i,j,k)*ysmth
                    str13n(i,j,k) = str13n(i,j,k)*ysmth
                    str22n(i,j,k) = str22n(i,j,k)*ysmth
                    str23n(i,j,k) = str23n(i,j,k)*ysmth
                    str33n(i,j,k) = str33n(i,j,k)*ysmth
#ENDIF
                end do
            end do
        end do
        !$omp end parallel do
        
        ! z-direction smoothing
        !$omp parallel do
        do j = 1,nz
            jj = j-1
            if (j .gt. nzh) jj = nz - j + 1
            zsmth = exp(-3.0*((float(jj)/float(nzh))**10))
            do i = 1,nyp
                do k = 1,nxh
                    omz(i,j,k) = omz(i,j,k)*zsmth
                    gn(i,j,k)  =  gn(i,j,k)*zsmth
                    fn(i,j,k)  =  fn(i,j,k)*zsmth
#IFDEF SCALAR
                    scn(i,j,k) = scn(i,j,k)*zsmth
#ENDIF
#IFDEF POLYMER
                    c11n(i,j,k) = c11n(i,j,k)*zsmth
                    c12n(i,j,k) = c12n(i,j,k)*zsmth
                    c13n(i,j,k) = c13n(i,j,k)*zsmth
                    c22n(i,j,k) = c22n(i,j,k)*zsmth
                    c23n(i,j,k) = c23n(i,j,k)*zsmth
                    c33n(i,j,k) = c33n(i,j,k)*zsmth
                    
                    str11n(i,j,k) = str11n(i,j,k)*zsmth
                    str12n(i,j,k) = str12n(i,j,k)*zsmth
                    str13n(i,j,k) = str13n(i,j,k)*zsmth
                    str22n(i,j,k) = str22n(i,j,k)*zsmth
                    str23n(i,j,k) = str23n(i,j,k)*zsmth
                    str33n(i,j,k) = str33n(i,j,k)*zsmth
#ENDIF
                end do
            end do
        end do
        !$omp end parallel do
        
        ! -------------------------------------------------------------------- !
        do i = 1,nyp
            h1n(i) = real(gn(i,1,1))
            h3n(i) = real(omz(i,1,1))
        end do
        
        ! On output: 
        !    fn = h2 = (u x om)_2
        !    gn = h1 = (u x om)_1
        !   omz = h3 = (u x om)_3
        !
        ! The subroutine nonlin(...) calculates the term fn in the fourth-
        ! order system for the vertical velocity, and the term gn for the 
        ! vertical velocity equation.
        
        call norm(fn)
        call norm(gn)
        call norm(omz)
        
        call nonlin(fn,omz,gn,wrkc,wrk1)
        
        call norm(fn)
        call norm(gn)
#IFDEF SCALAR
        call norm(scn)
#ENDIF
       
        ! -------------------------------------------------------------------- !
        ! Write restart file every iprnfrq time steps
        if (mod(it,iprnfrq) .eq. 0 .or. it .eq. nsteps) then
            open(122, file = 'outputs/last-restart',status='replace',form='unformatted')
            write(122) v,fn,fnm1,omy,gn,gnm1,u0,h1n,h1nm1,w0,h3n,h3nm1
            write(122) fxintg, fyintg, fzintg
            rewind(122)
            close(122)
#IF DEFINED SCALAR && DEFINED POLYMER
            if (it .ge. src_start) then
                open(123, file = 'outputs/c-last-restart',status='replace',form='unformatted')
                write(123) scalar,c11,c12,c13,c21,c22,c23,c31,c32,c33
                write(123) scn,c11n,c12n,c13n,c22n,c23n,c33n
                write(123) scnm1,c11nm1,c12nm1,c13nm1,c22nm1,c23nm1,c33nm1
                write(123) str11n,str12n,str13n,str22n,str23n,str33n
                write(123) str11nm1,str12nm1,str13nm1,str22nm1,str23nm1,str33nm1
                rewind(123)
                close(123)
            end if
#ENDIF
        end if

        ! -------------------------------------------------------------------- !
        ! -------------------------------------------------------------------- !
    end do ! it

! ---------------------------------------------------------------------------- !


! ============================================================================ !
!                               Post-Processing                                !
! ============================================================================ !

    call fftw_destroy_plan(planZb)
    call fftw_destroy_plan(planZf)
    call fftw_destroy_plan(planXb)
    call fftw_destroy_plan(planXf)
    call fftw_destroy_plan(planY)

! ---------------------------------------------------------------------------- !
end program dns

! ============================================================================ !
! ============================================================================ !
!                   Define setup and initial subroutines                       !
! ============================================================================ !
! ============================================================================ !

subroutine setstuff
!!----------------------------------------------------------------------------!! 
! This subroutine reads relevant input data from setup files and stores them   !
! in various common blocks for access elsewhere in the code                    !
!!----------------------------------------------------------------------------!! 

! ============================================================================ !
!                             Declare Modules                                  !
! ============================================================================ !
    use grid_size
! ---------------------------------------------------------------------------- !
 
! ============================================================================ !
!                      Declare all varibles explicitly                         !
! ============================================================================ !
    implicit none

    ! Temporary calculation variables
    real    :: pi
    integer :: i,j,k

    ! Solver variables
    real :: gain,ugain,theta,alpha,beta,dyde
    real :: wavx(nxh),wavz(nz) 
    real :: c(nyp)

    ! Geometry variables
    integer :: nyi,nzi,nxi,kwall,kmaxsurf
    real    :: xl,yl,zl
    integer, dimension(nyp,mz,mx) :: imatrix

    ! Buffer region
    real, dimension(1200) :: bfgain,bfugain
    integer :: bfhead,bftail,bfwidth
    real    :: vdes(mx)
    real    :: slopelength

    ! Flow setup variables
    real    :: re,Uinf,R_tau,dPdx
    integer :: geomtype,flow_select,perturbtime

    ! Simulation control variables
    integer :: irstrt,nsteps,iprnfrq,print3d,crstrt
    integer :: it
    real    :: dt

    ! Particle variables
    integer :: particle_flag,CD_switch
    real    :: ratio,ap,C_mu,gravity

    ! Vortex variables
    real    :: forbeta,xcenter,ycenter,zcenter,L,rad,bdyfx
    real    :: vortGamma,vortSigma,vortY,vortZ,vortSpace
    integer :: vNum
    
    ! Buffer region variables
    real    :: xstart

#IFDEF SCALAR
    ! Scalar variables
    real    :: sigmax,sigmay,sigmaz
    real    :: deltaT,diff
    integer :: scl_flag,scltarg
    integer :: src_start,src_stop
#ENDIF    
#IFDEF POLYMER
    ! Polymer variables
    integer :: ipolyflag,itarget,ipeter
    real    :: alpha_poly,tpoly,zlmax,diffpoly,qbeta,beta_min
    real    :: c11z,c22z,c33z
#ENDIF

! ---------------------------------------------------------------------------- !

! ============================================================================ !
!                            Define common blocks                              !
! ============================================================================ !

    common/iocontrl/   irstrt,nsteps,iprnfrq,print3d,crstrt
    common/solver/     gain,ugain,theta,alpha,beta,dyde
    common/waves/      wavx,wavz,c
    common/particles/  particle_flag,CD_switch,ratio,ap,C_mu,gravity
    common/setup/      geomtype,flow_select,perturbtime
    common/flow/       re,Uinf,R_tau,dPdx
    common/vortexring/ forbeta,xcenter,ycenter,zcenter,L,rad,bdyfx
    common/vort/       vortGamma,vortSigma,vortY,vortZ,vortSpace,vNum
    common/itime/      it
    common/dtime/      dt
    common/imat/       imatrix,kwall,kmaxsurf
    common/domain/     xl,yl,zl
    common/buffer/     bfgain,bfugain,vdes,bfhead,bftail,bfwidth
    common/slope/      slopelength
    common/buffx/      xstart
#IFDEF SCALAR
    common/scl_stuff/  sigmax,sigmay,sigmaz,deltaT,diff,scl_flag,scltarg
    common/src_time/   src_start,src_stop
#ENDIF
#IFDEF POLYMER
    common/poly_flgs/  ipolyflag,itarget,ipeter
    common/poly_var/   alpha_poly,tpoly,zlmax,diffpoly,qbeta
    common/c_init/     c11z,c22z,c33z
#ENDIF
! ---------------------------------------------------------------------------- !

! ============================================================================ !
!                        Read in data from dns.config                          !
! ============================================================================ !

    pi = 2.0*acos(0.0)

    open(110,file='setup/dns.config',status='old')

    ! Read in the data
    read(110,*) 
    read(110,*) 
    read(110,*) 
    read(110,*) nsteps
    read(110,*) iprnfrq
    read(110,*) dt    
    read(110,*) gain   
    read(110,*) ugain  
    read(110,*) theta  

    ! Read data characterizing geometry and fluid
    read(110,*) nyi   
    read(110,*) nzi    
    read(110,*) nxi     
    read(110,*) 
    read(110,*) yl       
    read(110,*) zl        
    read(110,*) xl         
    read(110,*) re          
    read(110,*) xstart
    read(110,*) uinf
    read(110,*) kwall
  
    ! Check consistency with parameter statements.
    if(nyi .ne. ny .or. nzi .ne. nz .or. nxi .ne. nx) then
         write(*,*) 'parameter statements inconsistent with ','prescribed grid distribution.'
         write(*,*) 'nyi = ',nyi
         write(*,*) 'ny = ',ny
         write(*,*) 'nzi = ',nzi
         write(*,*) 'nz = ',nz
         write(*,*) 'nxi = ',nxi
         write(*,*) 'nx = ',nx
         stop
    end if
  
    read(110,*) perturbtime
    read(110,*) R_tau
    read(110,*) gravity

    read(110,*) 
    read(110,*) 
    read(110,*) 
    read(110,*) 
    read(110,*) 

    ! Read in selection flags/switches 
    read(110,*) irstrt
    read(110,*) crstrt
    read(110,*) print3d
    read(110,*) geomtype
    read(110,*) particle_flag
    read(110,*) flow_select
#IFDEF POLYMER
    read(110,*) itarget
    read(110,*) ipeter
    read(110,*) ipolyflag
#ELSE
    read(110,*) 
    read(110,*) 
    read(110,*) 
#ENDIF
#IFDEF SCALAR
    read(110,*) scl_flag
#ELSE
    read(110,*)
#ENDIF

    read(110,*) 
    read(110,*) 
    read(110,*) 
    read(110,*) 
    read(110,*) 

    ! Read in data for particles
    read(110,*) 
    read(110,*) 
    read(110,*) ratio
    read(110,*) ap
    read(110,*) C_mu
    read(110,*) CD_switch

    read(110,*) 
    read(110,*) 
    read(110,*) 
    read(110,*) 
    read(110,*) 
#IFDEF POLYMER
    read(110,*) zlmax
    read(110,*) tpoly
    read(110,*) alpha_poly
    read(110,*) diffpoly
    read(110,*) c11z
    read(110,*) c22z
    read(110,*) c33z
    read(110,*) qbeta
    beta_min = qbeta
#ELSE
    read(110,*) 
    read(110,*) 
    read(110,*) 
    read(110,*) 
    read(110,*) 
    read(110,*) 
    read(110,*) 
    read(110,*) 
#ENDIF
#IFDEF SCALAR
    read(110,*) diff
    read(110,*) deltaT
    read(110,*) sigmax,sigmay,sigmaz
    read(110,*) src_start
    read(110,*) src_stop
    read(110,*) scltarg
#ENDIF
    close(110)

#IFDEF SCALAR
    if (src_stop .lt. src_start) then
        src_stop = nsteps
    end if
#ENDIF

! ============================================================================ !
!             Read in data for a vortex ring or streamwise vortex              !
! ============================================================================ !

    if (flow_select .eq. 5 .or. flow_select .ge. 10) then 
        open(111,file='setup/vort.config',status='old')
    
        read(111,*) 
        read(111,*) 
        read(111,*) 
        read(111,*) 
    
        ! Read in data for vortex ring
        read(111,*) bdyfx
        read(111,*) forbeta
        read(111,*) xcenter, ycenter, zcenter
        read(111,*) L, rad
    
        read(111,*) 
        read(111,*) 
        read(111,*) 
        read(111,*) 
        read(111,*) 
    
        read(111,*) vortGamma, vortSigma
        read(111,*) vortY, vortZ
        read(111,*) vNum 
        read(111,*) vortSpace
    
        close(111)
    end if

! ============================================================================ !
!                        Initialize Variables and FFTs                         !
! ============================================================================ !
    ! Set y coordinate vector for easy reference elsewhere
    ycoord(1) = yl/2.0
    do k = 2,nyp
        ycoord(k) = (cos(float(k-1)*pi/float(ny)))*(yl/2.0) ! [-YL/2,+YL/2]
        seght(k) = abs(ycoord(k) - ycoord(k-1))
    end do
    seght(1) = seght(nyp)
   
    ! Read in geometry matrix
    open(112,file='code/bin/geometry/geometry',form='unformatted')
    read(112) imatrix
    close(112)
    
    ! -------------------------------------------------------------------- !
    
    ! Calculate the resolvable wave numbers in x
    !   Assumes length XL has been normalized by length YL
    
    alpha = 2.0*pi/xl
    beta  = 2.0*pi/zl
    
    do k = 1,nxh
        wavx(k) = float(k-1)*alpha
    end do
    do j = 1,nzh
        wavz(j) = float(j-1)*beta
    end do
    do j = nzh+1,nz
        wavz(j) = float(j-nz-1)*beta
    end do
    
    ! Define Chebyshev coefficients
    do i = 1,nyp
        c(i) = 1.0
    end do
    c(1) = 2.0
    c(nyp) = 2.0

    ! -------------------------------------------------------------------- !
    ! Calculate dy/d(eta) scaling factor; see subroutine solve(...)
    !
    !   Assumes physical problem is solved on a space -yl < x2 < 0 and
    !   y = (2*x2 + yl)/yl. This scale factor is applied in the poisson
    !   solver routines and in the routines evaluating derivatives in the
    !   Chebyshev direction.
    
    dyde = 2.0/yl
    
    ! -------------------------------------------------------------------- !
    ! Set boundary conditions based on flow selection
    call setbcs

end subroutine setstuff

! ---------------------------------------------------------------------------- !

subroutine setbcs
!!----------------------------------------------------------------------------!! 
! This subroutine sets appropriate BC variables based on the flow selection    !
! switch instead of having to manually define it in the setup file each time.  !
! I added separate variables to handle u and w (I'm not sure why they were the !
! same by default in the first place).                                         !
!
! ---------------------------------------------------------------------------- !
! Set data for the calculation of the Green's functions phi1, phi2 in the      !
! subroutines gfcn1 & gfcn2. The composite solution ensures that v = 0 on both !
! boundaries normal to the Chebyshev direction by virtue of the BCs on vp, v1, !
! and v2. Specifying nt = 1 or 2 will ensure that either dv/dy(1) or           !
! d^2(v)/dy^2(1) vanishes at the top boundary, respectively. nb has the same   !
! effect for the bottom boundary, y = -1. The Dirichlet BCs for phi1 & phi2    !
! are arbitrary but must be independent.                                       !
!                                                                              !
! Set data for BCs on the zeroth Fourier mode of streamwise velocity. The      !
! general form of the BCs is:                                                  !
!                                                                              !
!     at*uo(1) + bt*d(uo)/dy(1) = gt(1)                                        !
!     ab*uo(-1) + bb*d(uo)/dy(-1) = gb(-1)                                     !
!                                                                              !
! The current formulation does not provide for ai,bi,gi to be time dependent.  !
! This capability can be included by providing recipes prior to the call of    !
! the subroutine uzero.                                                        !
!!----------------------------------------------------------------------------!! 

! ============================================================================ !
!                             Declare Modules                                  !
! ============================================================================ !
    use grid_size
! ---------------------------------------------------------------------------- !
 

! ============================================================================ !
!                      Declare all varibles explicitly                         !
! ============================================================================ !
    implicit none

    ! Boundary Condition Variables
    real    :: atu,btu,gtu,abu,bbu,gbu
    real    :: atw,btw,gtw,abw,bbw,gbw
    integer :: nt,nb
    real    :: p1t,p1b,p2t,p2b

    ! Flow setup variables
    real    :: re,Uinf,R_tau,dPdx

    ! Flow selection
    integer :: geomtype,flow_select

! ---------------------------------------------------------------------------- !


! ============================================================================ !
!                            Define common blocks                              !
! ============================================================================ !
    common/u0bcs/    atu,btu,gtu,abu,bbu,gbu
    common/w0bcs/    atw,btw,gtw,abw,bbw,gbw 
    common/grnfcn/   nt,nb,p1t,p1b,p2t,p2b  ! these won't change much
    common/setup/    geomtype,flow_select
    common/flow/     re,Uinf,R_tau,dPdx
! ---------------------------------------------------------------------------- !


! ============================================================================ !
!                            Set Boundary Conditions                           !
! ============================================================================ !

    print *,' Setting BCs...'
    print *,' flow_select = ',flow_select
    ! Setting BCs based on flow_select
    select case (flow_select)
        case (1,11) ! Channel flow
            print *,' BCs chosen for Channel flow:'
            write(*,*) ' u(0) = 0, u(yl) = 0'
            write(*,*) ' w(0) = 0, w(yl) = 0'
            print *,''
            ! These equate to no-slip BCs on top and bottom walls (y-direction)

            ! u=0  on y = 0 wall
            atu = 1.0
            btu = 0.0
            gtu = 0.0

            ! w=0 on y = 0 wall
            atw = 1.0
            btw = 0.0
            gtw = 0.0

            ! u=0 on y = yl wall
            abu = 1.0
            bbu = 0.0
            gbu = 0.0

            ! w=0 on y = yl wall
            abw = 1.0
            bbw = 0.0
            gbw = 0.0

            ! Green's Function BCs
            nt = 1.0
            nb = 1.0

        case (2,12) ! 2D Couette Flow
            print *,' BCs chosen for Couette flow:'
            write(*,*) ' u(0) = 0, u(yl) = ',Uinf
            write(*,*) ' w(0) = 0, w(yl) = 0'
            print *,''

            ! u = 0 on bottom wall, u = Uinf on top wall

            ! u=0 on y = 0 wall
            atu = 1.0
            btu = 0.0
            gtu = Uinf

            ! w=0 on y = 0 wall
            atw = 1.0
            btw = 0.0
            gtw = 0.0

            ! u=Uinf on y = yl wall
            abu = 1.0
            bbu = 0.0
            gbu = 0.0

            ! w=0 on y = yl wall
            abw = 1.0
            bbw = 0.0
            gbw = 0.0

            ! Green's Function BCs
            nt = 1.0
            nb = 1.0

        case (3) ! 3D Couette Flow
            print *,' BCs chosen for 3D Couette flow:'
            write(*,*) ' u(0) = 0, u(yl) = ',Uinf
            write(*,*) ' w(0) = 0, w(yl) = ',Uinf
            print *,''
            ! u = 0 on bottom wall, u = Uinf on top wall

            ! u=0 on y = 0 wall
            atu = 1.0
            btu = 0.0
            gtu = Uinf

            ! w=0 on y = 0 wall
            atw = 1.0
            btw = 0.0
            gtw = Uinf

            ! u=Uinf on y = yl wall
            abu = 1.0
            bbu = 0.0
            gbu = 0.0

            ! w=0 on y = yl wall
            abw = 1.0
            bbw = 0.0
            gbw = 0.0

            ! Green's Function BCs
            nt = 1.0
            nb = 1.0

        case (4) ! Blasius BL
            print *,' BCs chosen for Blasius BL:'
            write(*,*) ' u(0) = 0, du/dy(yl) = 0'
            write(*,*) ' w(0) = 0, w(yl) = 0'
            print *,''
            ! u = 0 on bottom wall, du/dy = 0 on top wall

            ! u=0 on y = 0 wall
            atu = 0.0
            btu = 1.0
            gtu = 0.0

            ! w=0 on y = 0 wall
            atw = 0.0
            btw = 1.0
            gtw = 0.0

            ! u=Uinf on y = yl wall
            abu = 1.0
            bbu = 0.0
            gbu = 0.0

            ! w=0 on y = yl wall
            abw = 1.0
            bbw = 0.0
            gbw = 0.0

            ! Green's Function BCs
            nt = 2.0
            nb = 1.0

        case (0,5,10) ! Vortex Ring or Still fluid
            print *,' BCs chosen for Still Fluid:'
            write(*,*) ' du/dy(0) = 0, du/dy(yl) = 0'
            write(*,*) ' dw/dy(0) = 0, dw/dy(yl) = 0'
            print *,''
            ! shear-free BCs on top and bottom

            ! du/dy=0 on y = 0 wall
            atu = 0.0
            btu = 1.0
            gtu = 0.0

            ! dw/dy=0 on y = 0 wall
            atw = 0.0
            btw = 1.0
            gtw = 0.0

            ! du/dy=0 on y = yl wall
            abu = 0.0
            bbu = 1.0
            gbu = 0.0

            ! dw/dy=0 on y = yl wall
            abw = 0.0
            bbw = 1.0
            gbw = 0.0

            ! Green's Function BCs
            nt = 2.0
            nb = 2.0

        case default
            print *,'Error: That is not a defined flow selection number. Please select a number from the following list:'
            print *,'       0,1,2,3,4,5,10,11,12'
            print *,''
    end select

    ! For all the cases so far, these numbers have not changed, so I'm defining
    ! them outside the selection loop
 
    p1b = 1.0
    p1t = -1.0
    p2b = -1.0
    p2t = 0.5 

end subroutine setbcs

! ---------------------------------------------------------------------------- !

subroutine forcepre

! ============================================================================ !
!                             Declare Modules                                  !
! ============================================================================ !
    use grid_size
! ---------------------------------------------------------------------------- !


! ============================================================================ !
!                      Declare all varibles explicitly                         !
! ============================================================================ !

    implicit none

    ! Flow variables
    real, dimension(nyp,nz,nx) :: initu,initv,initw
    real    :: re,Uinf,R_tau,dPdx

    ! Buffer region
    real, dimension(1200) :: bfgain,bfugain
    integer :: bfhead,bftail,bfwidth,bfwh,bfcntr
    real    :: vdes(mx), xsuction 
    real    :: xstart
    
    ! Indices
    integer :: i,j,k,ib
    
    ! Geometry
    integer, dimension(nyp,mz,mx) :: imatrix 
    integer :: kwall,kmaxsurf
    real    :: slopelength,ywall,ytop
    real    :: xl,yl,zl
    
    ! Calculation variables
    real :: delxm,x,vtemp,pi
    
    ! Immersed boundary force variables
    real, dimension(nyp,mz,mx) :: fxintg,fyintg,fzintg
    real, dimension(mz2,mx2)   :: fspread
    real :: gwidth,dsq
    
    ! Simulation control variables
    integer :: irstrt,nsteps,iprnfrq,print3d

! ---------------------------------------------------------------------------- !


! ============================================================================ !
!                            Define common blocks                              !
! ============================================================================ !

    common/iocontrl/ irstrt,nsteps,iprnfrq,print3d
    common/ibforce/  fxintg,fyintg,fzintg,fspread
    common/init/     initu,initv,initw
    common/imat/     imatrix,kwall,kmaxsurf
    common/buffer/   bfgain,bfugain,vdes,bfhead,bftail,bfwidth
    common/slope/    slopelength
    common/flow/     re,Uinf,R_tau,dPdx
    common/domain/   xl,yl,zl
    common/buffx/    xstart

! ---------------------------------------------------------------------------- !

! ============================================================================ !
!                   Calculate preliminary data for forcing terms               !
! ============================================================================ !

    delxm = xl/float(mx-1) ! x-spacing in physical space
    pi = 2.0*acos(0.0)

    ! -------------------------------------------------------------------- !
    ! Buffer region parameters
    bfhead  = 2
    bftail  = bftail_
    bfwidth = (bftail - bfhead)
    bfwh    = bfwidth/2
    bfcntr  = bfwh + 1
    bfgain  = 0.0
    bfugain = 0.0

    do i = 1,(bfwidth+1)
        ib = i - bfcntr
        if (iabs(ib) .le. bfwh/3) then
            bfgain(i) = 20.0*exp(-60.0*(float(ib)/float(bfwh))**2)
        end if
        bfugain(i) = 20.0*exp(-6.0*(float(ib)/float(bfwh))**2)
    end do

    slopelength = bfwidth/2.0
    if (kwall .ne. 0) then
        ywall = ycoord(kwall)
    else
        ywall = 0.0
    end if
    ytop = yl - ywall

    ! -------------------------------------------------------------------- !
    ! Initialize the flowfield

    call init_flow(initu,initv,initw) ! Sets initial flow profile for the domain

    ! -------------------------------------------------------------------- !
    ! Set geometry type for immersed boundary forces
    kmaxsurf = 0
    do k = 1,mx
        do j = 1,mz
            do i = 1,nyp
                if ((imatrix(i,j,k) .eq. 1) .or. (imatrix(i,j,k) .eq. 3) .or. (imatrix(i,j,k) .eq. 4)) then
                    if (k .ge. kmaxsurf) kmaxsurf = k
                end if
            end do
        end do
    end do

    xsuction = delxm*(bftail-1)

    do i = 1,mx
        x = delxm*(i-1)
        if (x .le. xsuction) then
            vdes(i) = 0.0
        else
            vdes(i) = (-0.8604*Uinf)/sqrt(re*Uinf*(x + xstart - xsuction))
        end if
    end do

    ! -------------------------------------------------------------------- !
    ! Calculate force terms
    ! Force spreading matrix
    gwidth = 2.0
    do i = 1,mx2
        do j = 1,mz2
            dsq = (float(i-1))**2 + (float(j-1))**2
            fspread(j,i) = 0.0
            if (dsq .lt. 13) fspread(j,i) = exp(-gwidth*dsq)
        end do
    end do

    ! Initialize integral forcing terms
    if (irstrt .eq. 0) then
        do k = 1,mx
            do j = 1,mz
                do i = 1,nyp
                    fxintg(i,j,k) = 0.0
                    fyintg(i,j,k) = 0.0
                    fzintg(i,j,k) = 0.0
                end do
            end do
        end do
    end if

    ! -------------------------------------------------------------------- !
    ! Low resolution restart case, commented because it's unused and untested
!   if (irstrt .eq. 2) then
!       numjstrips = mz/mzlowres
!       open(90,file='setup/low-resolution-forces',status='old',form='unformatted')
!       read(90) fxintLR
!       read(90) fyintLR
!       read(90) fzintLR
!       close(90)
!
!       do i = 1,nyp
!           do jstrip = 1,numjstrips
!               do jcount = 1,mzlowres
!                   j = ((jstrip-1)*mzlowres) + jcount
!                   do k = 1,mx
!                       fxintg(i,j,k) = fxintLR(i,jcount,k)
!                       fyintg(i,j,k) = fxintLR(i,jcount,k)
!                       fzintg(i,j,k) = fxintLR(i,jcount,k)
!                   end do
!               end do
!           end do
!       end do
!   end if


end subroutine forcepre

! ---------------------------------------------------------------------------- !

subroutine initial(u,u0,v,w,w0,omx,omy,omz,fn,fnm1,gn,gnm1,h1n,h1nm1,h3n,h3nm1, &
                   wrk1,wrkc,u11,u12,u13,u21,u22,u23,u31,u32,u33, &
#IFDEF SCALAR
                   scalar,sclx,scly,sclz,scn,scnm1,     &
#ENDIF
#IFDEF POLYMER
                   t1,t2,t3,c11,c12,c13,c21,c22,c23,c31,c32,c33,dc111,dc112,    &
                   dc113,dc121,dc122,dc123,dc131,dc132,dc133,dc211,dc212,dc213, &
                   dc221,dc222,dc223,dc231,dc232,dc233,dc311,dc312,dc313,dc321, &
                   dc322,dc323,dc331,dc332,dc333,c11n,c12n,c13n,c22n,c23n,c33n, &
                   c11nm1,c12nm1,c13nm1,c22nm1,c23nm1,c33nm1,str11n,str12n,     &
                   str13n,str22n,str23n,str33n,str11nm1,str12nm1,str13nm1,      &
                   str22nm1,str23nm1,str33nm1,qp11,qp12,qp13,qp22,qp23,qp33,    &
#ENDIF
                   Lu,Lv,Lw,Lu_old,Lv_old,Lw_old, &
                   planZb,planXb,planY,planXf,planZf) 

! ============================================================================ !
!                             Declare Modules                                  !
! ============================================================================ !
    use,intrinsic :: iso_c_binding
    use grid_size
    use solvers
    use derivs
    use helpers
! ---------------------------------------------------------------------------- !


! ============================================================================ !
!                      Declare all varibles explicitly                         !
! ============================================================================ !

    implicit none

    include 'fftw3.f03'

    ! FFTW variables
    type(C_PTR) :: planZb,planXb,planY,planXf,planZf

    ! Flow variables
    complex, dimension(nyp,nz,nxh) :: u,v,w,omx,omy,omz
    complex, dimension(nyp,nz,nxh) :: u11,u12,u13
    complex, dimension(nyp,nz,nxh) :: u21,u22,u23
    complex, dimension(nyp,nz,nxh) :: u31,u32,u33
    complex, dimension(nyp,nz,nxh) :: Lu,Lv,Lw
    real,    dimension(nyp,mz,mx)  :: Lu_old,Lv_old,Lw_old
    real,    dimension(nyp,nz,nx)  :: initu,initv,initw

#IFDEF SCALAR
    ! Scalar variables
    complex, dimension(nyp,nz,nxh) :: scalar,sclx,scly,sclz
    complex, dimension(nyp,nz,nxh) :: scn,scnm1,csource
    real,    dimension(nyp,nz,nx)  :: psource,scl
    integer :: src_start,src_stop
#ENDIF
#IFDEF POLYMER    
    ! Polymer variables
    complex, dimension(nyp,nz,nxh) :: c11,c12,c13,c21,c22,c23,c31,c32,c33
    complex, dimension(nyp,nz,nxh) :: c11n,c12n,c13n,c22n,c23n,c33n
    complex, dimension(nyp,nz,nxh) :: c11nm1,c12nm1,c13nm1,c22nm1,c23nm1,c33nm1

    complex, dimension(nyp,nz,nxh) :: dc111,dc112,dc113,dc121,dc122,dc123,dc131,dc132,dc133
    complex, dimension(nyp,nz,nxh) :: dc211,dc212,dc213,dc221,dc222,dc223,dc231,dc232,dc233
    complex, dimension(nyp,nz,nxh) :: dc311,dc312,dc313,dc321,dc322,dc323,dc331,dc332,dc333

    complex, dimension(nyp,nz,nxh) :: str11n,str12n,str13n,str22n,str23n,str33n
    complex, dimension(nyp,nz,nxh) :: str11nm1,str12nm1,str13nm1,str22nm1,str23nm1,str33nm1

    complex, dimension(nyp,nz,nxh) :: qp11,qp12,qp13,qp22,qp23,qp33

    complex, dimension(nyp,nz,nxh) :: t1,t2,t3

    integer :: ipolyflag,itarget,ipeter
    real    :: alpha_poly,tpoly,zlmax,diffpoly,qbeta
    real    :: c11z,c22z,c33z
#ENDIF

    ! Solver variables
    complex, dimension(nyp,nz,nxh) :: fn,fnm1,gn,gnm1
    real,    dimension(nyp)        :: u0,h1n,h1nm1
    real,    dimension(nyp)        :: w0,h3n,h3nm1
    real,    dimension(nxh)        :: wavx
    real,    dimension(nz)         :: wavz

    ! Temporary variables (multi-purpose calculation variables)
    complex, dimension(nyp,nz,nxh) :: wrkc, wrk1
#IFDEF POLYMER
    real,    dimension(nyp,nz,nx)  :: wrk11,wrk12,wrk13,wrk21,wrk22,wrk23,wrk31,wrk32,wrk33
#ENDIF
    integer :: i,j,k
    complex :: im
    
    
    ! Simulation control variables
    integer :: irstrt,nsteps,iprnfrq,print3d,crstrt
    real    :: dt
    integer :: it
    
    ! Immersed boundary force variables
    real, dimension(nyp,mz,mx) :: fxintg,fyintg,fzintg
    real, dimension(mz2,mx2)   :: fspread

! ---------------------------------------------------------------------------- !


! ============================================================================ !
!                            Define common blocks                              !
! ============================================================================ !

    common/iocontrl/  irstrt,nsteps,iprnfrq,print3d,crstrt
    common/waves/     wavx,wavz
    common/init/      initu,initv,initw
    common/ibforce/   fxintg,fyintg,fzintg,fspread
    common/itime/     it
    common/dtime/     dt
#IFDEF SCALAR
    common/src_time/  src_start,src_stop
#ENDIF
#IFDEF POLYMER
    common/poly_flgs/ ipolyflag,itarget,ipeter
    common/poly_var/  alpha_poly,tpoly,zlmax,diffpoly,qbeta
    common/c_init/    c11z,c22z,c33z
#ENDIF

! ---------------------------------------------------------------------------- !

! ============================================================================ !
!                          Set up initial conditions                           !
! ============================================================================ !

    ! -------------------------------------------------------------------- !
    ! Sets the initial condition or reads the restart file. The first time !
    ! step has to be handled differently than the rest of the code because !
    ! the explicit time integration requires an (n-1) value. This routine  !
    ! sets initial conditions in real space and then transforms to         !
    ! spectral space.                                                      !
    ! -------------------------------------------------------------------- !
   
    im = (0.0,1.0) ! imaginary number, i
    
    if (irstrt .eq. 1) then
        open(113,file='setup/restart',status='old',form='unformatted')
!        read(113)
        read(113) v,fn,fnm1,omy,gn,gnm1,u0,h1n,h1nm1,w0,h3n,h3nm1
        read(113) fxintg,fyintg,fzintg
        close(113)

#IF DEFINED SCALAR && !DEFINED POLYMER
    ! -------------------------------------------------------------------- !
    ! Initialize scalar field
        call setscl(scl,psource)
        call xyzfft(scalar,scl,-1)

        call xyzfft(csource,psource,-1)
#ELIF DEFINED SCALAR
        if (crstrt .eq. 0) then
            ! For now, just manually set Cij and scalar and leave nonlinear terms as 0
            call setpoly(scl,psource,wrk11,wrk12,wrk13,wrk21,wrk22,wrk23,wrk31,wrk32,wrk33)

            call xyzfft(scalar,scl,-1)
            call xyzfft(csource,psource,-1)
    
            call xyzfft(c11,wrk11,-1)
            call xyzfft(c12,wrk12,-1)
            call xyzfft(c13,wrk13,-1)
            call xyzfft(c21,wrk21,-1)
            call xyzfft(c22,wrk22,-1)
            call xyzfft(c23,wrk23,-1)
            call xyzfft(c31,wrk31,-1)
            call xyzfft(c32,wrk32,-1)
            call xyzfft(c33,wrk33,-1)
    
            ! Calculate scalar gradient
            call gradscl(scalar,sclx,scly,sclz,wrkc)
            call norm(sclx)
            call norm(scly)
            call norm(sclz)
    
            ! Calculate derivatives of conformation tensor        
            call derivscji(c11,c12,c13,c21,c22,c23,c31,c32,c33,                   &
                           dc111,dc112,dc113,dc211,dc212,dc213,dc311,dc312,dc313, &
                           dc121,dc122,dc123,dc221,dc222,dc223,dc321,dc322,dc323, &
                           dc131,dc132,dc133,dc231,dc232,dc233,dc331,dc332,dc333)  
      
            call norm(c11)
            call norm(c12)
            call norm(c13)
            call norm(c21)
            call norm(c22)
            call norm(c23)
            call norm(c31)
            call norm(c32)
            call norm(c33)
    
            call norm(dc111)
            call norm(dc112)
            call norm(dc113)
            call norm(dc211)
            call norm(dc212)
            call norm(dc213)
            call norm(dc311)
            call norm(dc312)
            call norm(dc313)
            call norm(dc121)
            call norm(dc122)
            call norm(dc123)
            call norm(dc221)
            call norm(dc222)
            call norm(dc223)
            call norm(dc321)
            call norm(dc322)
            call norm(dc323)
            call norm(dc131)
            call norm(dc132)
            call norm(dc133)
            call norm(dc231)
            call norm(dc232)
            call norm(dc233)
            call norm(dc331)
            call norm(dc332) 
            call norm(dc333)                                                 
    
            ! Set nonlinear terms to 0
            scn = 0.0
    
            c11n = 0.0
            c12n = 0.0
            c13n = 0.0
            c22n = 0.0
            c23n = 0.0
            c33n = 0.0
    
            str11n = 0.0
            str12n = 0.0
            str13n = 0.0
            str22n = 0.0
            str23n = 0.0
            str33n = 0.0
    
            qp11 = 0.0
            qp12 = 0.0
            qp13 = 0.0
            qp22 = 0.0
            qp23 = 0.0
            qp33 = 0.0
        else
            open(114,file='setup/c-restart',status='old',form='unformatted')
            read(114) scalar,c11,c12,c13,c21,c22,c23,c31,c32,c33
            read(114) scn,c11n,c12n,c13n,c22n,c23n,c33n
            read(114) scnm1,c11nm1,c12nm1,c13nm1,c22nm1,c23nm1,c33nm1
            read(114) str11n,str12n,str13n,str22n,str23n,str33n
            read(114) str11nm1,str12nm1,str13nm1,str22nm1,str23nm1,str33nm1
            close(114)
        end if ! Scalar/polymer restart
#ENDIF
    else ! assumes only options are 1 and 0
    ! -------------------------------------------------------------------- !
    ! Use initial flow values in real space to set u,v,w in spectral space

        call xyzfft(u,initu,-1)
        call xyzfft(v,initv,-1)
        call xyzfft(w,initw,-1)

#IF DEFINED SCALAR && !DEFINED POLYMER
    ! -------------------------------------------------------------------- !
    ! Initialize scalar field
        call setscl(scl,psource)

        call xyzfft(scalar,scl,-1)
        call xyzfft(csource,psource,-1)
#ELIF DEFINED SCALAR
        call setpoly(scl,psource,wrk11,wrk12,wrk13,wrk21,wrk22,wrk23,wrk31,wrk32,wrk33)

        call xyzfft(scalar,scl,-1)
        call xyzfft(csource,psource,-1)
    
        call xyzfft(c11,wrk11,-1)
        call xyzfft(c12,wrk12,-1)
        call xyzfft(c13,wrk13,-1)
        call xyzfft(c21,wrk21,-1)
        call xyzfft(c22,wrk22,-1)
        call xyzfft(c23,wrk23,-1)
        call xyzfft(c31,wrk31,-1)
        call xyzfft(c32,wrk32,-1)
        call xyzfft(c33,wrk33,-1)
#ENDIF

    ! -------------------------------------------------------------------- !
    ! Calculate other variables
        ! Calculate omy
        do k = 1,nxh
            do j = 1,nz
                do i = 1,nyp
                    omy(i,j,k) = im*wavz(j)*u(i,j,k) - im*wavx(k)*w(i,j,k)
                end do
            end do
        end do
    
        ! Because of periodicity, the (0,0) mode for omy is always 0
        do i = 1,nyp
            omy(i,1,1) = (0.0,0.0)
        end do
   
        ! Normalize variables
        call norm(u)
        call norm(v)
        call norm(w)
        call norm(omy)

#IFDEF SCALAR
        call norm(scalar)
        call norm(csource)
#ENDIF
#IFDEF POLYMER        
        call norm(c11)
        call norm(c12)
        call norm(c13)
        call norm(c21)
        call norm(c22)
        call norm(c23)
        call norm(c31)
        call norm(c32)
        call norm(c33)
#ENDIF

        ! Initialize u0/w0: 1st Fourier mode of streamwise/spanwise velocity
        do i = 1,nyp
            u0(i) = real(u(i,1,1))
            w0(i) = real(w(i,1,1))
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

        call norm(u)
        call norm(v)
        call norm(w)
        call norm(omy)
       
        ! Calculate remaining vorticity and nonlinear terms
        call vort(u,v,w,omx,omz)
        
        call norm(omx)
        call norm(omz)

#IFDEF SCALAR        
        ! Calculate scalar gradient
        call gradscl(scalar,sclx,scly,sclz,wrkc)
        call norm(sclx)
        call norm(scly)
        call norm(sclz)
#ENDIF

        ! Calculate derivatives for velocity gradient tensor & Laplacian
        call gradv(u,v,w,u11,u12,u13,u21,u22,u23,u31,u32,u33,Lu,Lv,Lw)
    
        call norm(u11)
        call norm(u12)
        call norm(u13)
        call norm(u21)
        call norm(u22)
        call norm(u23)
        call norm(u31)
        call norm(u32)
        call norm(u33)
    
        call norm(Lu)
        call norm(Lv)
        call norm(Lw)

#IFDEF POLYMER    
        ! Calculate derivatives of conformation tensor        
        call derivscji(c11,c12,c13,c21,c22,c23,c31,c32,c33,                   &
                       dc111,dc112,dc113,dc211,dc212,dc213,dc311,dc312,dc313, &
                       dc121,dc122,dc123,dc221,dc222,dc223,dc321,dc322,dc323, &
                       dc131,dc132,dc133,dc231,dc232,dc233,dc331,dc332,dc333)  
  
        call norm(c11)
        call norm(c12)
        call norm(c13)
        call norm(c21)
        call norm(c22)
        call norm(c23)
        call norm(c31)
        call norm(c32)
        call norm(c33)

        call norm(dc111)
        call norm(dc112)
        call norm(dc113)
        call norm(dc211)
        call norm(dc212)
        call norm(dc213)
        call norm(dc311)
        call norm(dc312)
        call norm(dc313)
        call norm(dc121)
        call norm(dc122)
        call norm(dc123)
        call norm(dc221)
        call norm(dc222)
        call norm(dc223)
        call norm(dc321)
        call norm(dc322)
        call norm(dc323)
        call norm(dc131)
        call norm(dc132)
        call norm(dc133)
        call norm(dc231)
        call norm(dc232)
        call norm(dc233)
        call norm(dc331)
        call norm(dc332) 
        call norm(dc333)                                                 
#ENDIF

!    write(1,*) u
    ! -------------------------------------------------------------------- !
    ! Compute v x omega in physical space - FFTs are performed inside subroutine
    ! Assumes all variables are in 3D spectral space (on input and output)
        call vcw3dp(u,v,w,omx,omy,omz,fn,gn,u11,u12,u13,u21,u22,u23,       &
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
   
        call norm(v)
        call norm(omy)
        call norm(fn)
        call norm(gn)
        call norm(omz)

#IFDEF POLYMER    
        ! Normalize variables
        call norm(c11)
        call norm(c12)
        call norm(c13)
        call norm(c21)
        call norm(c22)
        call norm(c23)
        call norm(c31)
        call norm(c32)
        call norm(c33)

        call norm(c11n)
        call norm(c12n)
        call norm(c13n)
        call norm(c22n)
        call norm(c23n)
        call norm(c33n)

        call norm(str11n)
        call norm(str12n)
        call norm(str13n)
        call norm(str22n)
        call norm(str23n)
        call norm(str33n)

        call norm(qp11) 
        call norm(qp12)
        call norm(qp13)
        call norm(qp22)
        call norm(qp23)
        call norm(qp33)

        ! -------------------------------------------------------------------- !
        ! Add polymer forces to the momentum equation
        ! qp is the polymer stress tensor, t = div(qp) which is the force on 
        ! the fluid due to polymer stress
        if (it .ge. src_start) then
            call polyforce(qp11,qp12,qp13,qp22,qp23,qp33,t1,t2,t3)
            call norm(t1)
            call norm(t2)
            call norm(t3)
            call subforce(gn,fn,omz,t1,t2,t3)

        end if
#ENDIF
#IFDEF SCALAR
        ! -------------------------------------------------------------------- !
        ! Add source to scalar term
        if (src_start .eq. 0) then
            do k = 1,nxh
                do j = 1,nz
                    do i = 1,nyp
                        scn(i,j,k) = scn(i,j,k) + csource(i,j,k)
                    end do
                end do
            end do
        end if
#ENDIF        
        ! -------------------------------------------------------------------- !
        do i=1,nyp
            h1n(i) = real(gn(i,1,1))
            h3n(i) = real(omz(i,1,1))
        end do
    
        call nonlin(fn,omz,gn,wrkc,wrk1)
    
        call norm(fn)
        call norm(gn)

#IFDEF SCALAR
        call norm(scn)
#ENDIF
    
        ! For the first time step, fnm1=fn & h1nm1=h1n
    
        do k=1,nxh
            do j=1,nz
                do i=1,nyp
                     fnm1(i,j,k) =  fn(i,j,k)
                     gnm1(i,j,k) =  gn(i,j,k)
#IFDEF SCALAR
                    scnm1(i,j,k) = scn(i,j,k)
#ENDIF
                end do
            end do
        end do
    
        do i=1,nyp
            h1nm1(i)=h1n(i)
            h3nm1(i)=h3n(i)
        end do
#IFDEF POLYMER   
        ! Add Cij terms 
        do k=1,nxh
            do j=1,nz
                do i=1,nyp
                    c11nm1(i,j,k) = c11n(i,j,k)
                    c12nm1(i,j,k) = c12n(i,j,k)
                    c13nm1(i,j,k) = c13n(i,j,k)
                    c22nm1(i,j,k) = c22n(i,j,k)
                    c23nm1(i,j,k) = c23n(i,j,k)
                    c33nm1(i,j,k) = c33n(i,j,k)
                    str11nm1(i,j,k) = str11n(i,j,k)
                    str12nm1(i,j,k) = str12n(i,j,k)
                    str13nm1(i,j,k) = str13n(i,j,k)
                    str22nm1(i,j,k) = str22n(i,j,k)
                    str23nm1(i,j,k) = str23n(i,j,k)
                    str33nm1(i,j,k) = str33n(i,j,k)
                end do
            end do
        end do
#ENDIF
    end if
end subroutine initial

! ---------------------------------------------------------------------------- !

#IFDEF POLYMER
subroutine setpoly(scl,psource,wrk11,wrk12,wrk13,wrk21,wrk22,wrk23,wrk31,wrk32,wrk33)

! ============================================================================ !
!                             Declare Modules                                  !
! ============================================================================ !
    use grid_size
! ---------------------------------------------------------------------------- !

! ============================================================================ !
!                      Declare all varibles explicitly                         !
! ============================================================================ !

    implicit none

    ! Passed variables
    real, dimension(nyp,nz,nx) :: scl,psource
    real, dimension(nyp,nz,nx) :: wrk11,wrk12,wrk13,wrk21,wrk22,wrk23,wrk31,wrk32,wrk33

    ! Scalar variables
    real    :: xshift,yshift,zshift,sigmax,sigmay,sigmaz
    real    :: deltaT,diff
    integer :: scl_flag

    ! Geometry variables
    real    :: xl,yl,zl

    ! Calculation variables
    integer :: i,j,k,n
    real    :: pi
    real    :: xc1,yc1,zc1,betax,betay,betaz
    real    :: xcor,zcor,xsq,ysq,zsq

    ! Polymer variables
    real    :: c11z,c22z,c33z
    integer :: src_start,src_stop

! ---------------------------------------------------------------------------- !

! ============================================================================ !
!                            Define common blocks                              !
! ============================================================================ !

    common/domain/    xl,yl,zl
    common/scl_stuff/ sigmax,sigmay,sigmaz,deltaT,diff,scl_flag
    common/c_init/    c11z,c22z,c33z
    common/src_time/  src_start,src_stop

! ---------------------------------------------------------------------------- !

! ============================================================================ !
!                              Begin Calculations                              !
! ============================================================================ !

    pi = 2.0*acos(0.0)

    ! Initialize scalar and source
    scl = 0.0
    psource = 0.0

    ! Use initial particle data for scalar/source location(s)
    do n = 1,npart
        open(95,file='setup/particles/particles.dat',status='old',action='read')
        do j = 1,n
            read(95,*)
        end do
        read(95,*) xc1,yc1,zc1
        close(95)

        do k = 1,nx
            xcor = xl*(float(k-1)/float(nx))
            xsq  = (xcor - xc1)**2
            betax = xsq/(2.0*sigmax**2)
            do j = 1,nz
                zcor = zl*(float(j-1)/float(nz))
                zsq  = (zcor - zc1)**2
                betaz = zsq/(2.0*sigmaz**2)
                do i = 1,nyp
                    ysq = (ycoord(i) - yc1)**2
                    betay = ysq/(2.0*sigmay**2)

                    ! Initialize Cij
                    wrk11(i,j,k) = c11z
                    wrk12(i,j,k) = 0.0
                    wrk13(i,j,k) = 0.0
                    wrk21(i,j,k) = 0.0
                    wrk22(i,j,k) = c22z
                    wrk23(i,j,k) = 0.0
                    wrk31(i,j,k) = 0.0
                    wrk32(i,j,k) = 0.0
                    wrk33(i,j,k) = c33z

                    ! Initialize scalar
                    if (scl_flag .eq. 1) then
                        scl(i,j,k) = deltaT*exp(-(betax + betay + betaz))
                    else if (scl_flag .eq. 2) then  
                        if (src_start .eq. 0) then
                            psource(i,j,k) = deltaT*exp(-(betax + betay + betaz))
                        end if
                    end if
                end do
            end do
        end do
    end do

    if (scl_flag .eq. 3) then
        scl(i,j,k) = 1.0
    end if

end subroutine setpoly
#ENDIF
! ---------------------------------------------------------------------------- !
#IF DEFINED SCALAR && !DEFINED POLYMER
subroutine setscl(scl,psource)

! ============================================================================ !
!                             Declare Modules                                  !
! ============================================================================ !
    use grid_size
! ---------------------------------------------------------------------------- !

! ============================================================================ !
!                      Declare all varibles explicitly                         !
! ============================================================================ !

    implicit none

    ! Passed variables
    real, dimension(nyp,nz,nx) :: scl,psource

    ! Scalar variables
    real    :: xshift,yshift,zshift,sigmax,sigmay,sigmaz
    real    :: deltaT,diff
    integer :: scl_flag

    ! Geometry variables
    real    :: xl,yl,zl

    ! Calculation variables
    integer :: i,j,k,n
    real    :: pi
    real    :: xc1,yc1,zc1,betax,betay,betaz
    real    :: xcor,zcor,xsq,ysq,zsq

    ! Source variables
    integer :: src_start,src_stop

! ---------------------------------------------------------------------------- !

! ============================================================================ !
!                            Define common blocks                              !
! ============================================================================ !

    common/domain/    xl,yl,zl
    common/scl_stuff/ sigmax,sigmay,sigmaz,deltaT,diff,scl_flag
    common/src_time/  src_start,src_stop

! ---------------------------------------------------------------------------- !

! ============================================================================ !
!                              Begin Calculations                              !
! ============================================================================ !

    pi = 2.0*acos(0.0)

    ! Initialize scalar and source
    scl = 0.0
    psource = 0.0

    ! Use initial particle data for scalar/source location(s)
    do n = 1,npart
        open(95,file='setup/particles/particles.dat',status='old',action='read')
        do j = 1,n
            read(95,*)
        end do
        read(95,*) xc1,yc1,zc1
        close(95)
        
        do k = 1,nx
            xcor = xl*(float(k-1)/float(nx))
            xsq  = (xcor - xc1)**2
            betax = xsq/(2.0*sigmax**2)
            do j = 1,nz
                zcor = zl*(float(j-1)/float(nz))
                zsq  = (zcor - zc1)**2
                betaz = zsq/(2.0*sigmaz**2)
                do i = 1,nyp
                    ysq = (ycoord(i) - yc1)**2
                    betay = ysq/(2.0*sigmay**2)
                    if (scl_flag .eq. 1) then
                        scl(i,j,k) = scl(i,j,k) + deltaT*exp(-(betax + betay + betaz))
                    else if (scl_flag .eq. 2) then  
                        if (src_start .eq. 0) then
                            psource(i,j,k) = psource(i,j,k) + deltaT*exp(-(betax + betay + betaz))
                        end if
                    end if
                end do
            end do
        end do
    end do

    if (scl_flag .eq. 3) then
       scl = 1.0
    end if 

end subroutine setscl
#ENDIF
