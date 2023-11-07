subroutine part_track(u,v,w,omx,omy,omz,u_old,v_old,w_old, & 
                      u11,u12,u13,u21,u22,u23,u31,u32,u33, &
                      Lu,Lv,Lw,Lu_old,Lv_old,Lw_old)

! -------------------------------------------------------- !
! This is an implementation of Lagrangian particle tracking 
! for passive inertial and tracer particles. This 
! subroutine uses flowfield data:
! -------------------------------------------------------- !
!                       INPUTS
!
!   u, v, w             : x-, y-, and z-components of fluid 
!                         velocity from DNS solver
!   omx, omy, omz       : x-, y-, and z-components of fluid
!                         vorticity from DNS solver
!   u_old, v_old, w_old : u,v,w from previous time step for
!                         computing Du/Dt
!   u11-u33             : Velocity gradient components from
!                         DNS solver (u11 = du/dx, 
!                         u12 = du/dy, u33 = dw/dz, etc.) 
!   Lu, Lv, Lw          : Laplacian of velocity field, used 
!                         to compute Faxen correction terms
! -------------------------------------------------------- !
!
! This subroutine does NOT change the flowfield data. The 
! only outputs are files written to the output directory.
! Each particle output file contains position, velocity,
! and swirl strength at the particle location for each
! particle initialized and read from a setup file called
! 'particles.dat'
!
! -------------------------------------------------------- !
!                       OUTPUTS                            
!
!   part-t-XXXXX.dat : Particle data output file, to be 
!                      read by Tecplot
!   swirlX           : Percentage of particles inside a
!                      vortex, identified by swirl > X
! -------------------------------------------------------- !
!
! -------------------------------------------------------- !
!                        METHOD
!
! The following is the general procedure followed by this 
! routine:
!
! 1. Initialize particles - Read in data and force initial
!    velocity to match that of fluid
!
! 2. Check if sub-step is needed (smaller time step for 
!    RK4 needed than fluid dt)
!
! 3. Integrate particle velocity to find new particle 
!    position. Uses RK4 integrator for inertial particles
!   
!    The following steps take place inside a loop, where 
!    each iteration is a unique particle trajectory
!
! 3a. Interpolate fluid velocity and vorticity at particle
!     location. 
! 3b. Compute Du/Dt at nearest neighbors and interpolate.
! 3c. Compute d/dt(laplacian(u)) at nearest neighbors and
!     interpolate.
! 3d. Call RKStage(...) to compute forces for each stage
!     based on MR equations
! 3e. Update particle position
!
! 4. Compute (interpolate) swirl strength at new position
!
! 5. Compute correlation percentage (swirlX file)
!
! 6. Write particle data file 
! -------------------------------------------------------- !
! -------------------------------------------------------- !

    ! Declare variables and modules

    use grid_size

    implicit none

    ! Common blocks
    real :: pi, dt
    common/data1/ pi,dt

    integer :: it
    common/pre5f/ it

    integer :: irstrt, nsteps, iprnfrq ! irstrt and nsteps unused
    common/iocntrl/ irstrt, nsteps, iprnfrq

    integer :: particle_flag
    common/particle/ particle_flag

    ! Flowfield variables
    real,dimension(nyp,mz,mx) :: u,v,w,u_old,v_old,w_old,omx,omy,omz
    real,dimension(nyp,mz,mx) :: u11,u12,u13,u21,u22,u23,u31,u32,u33
    real,dimension(nyp,mz,mx) :: Lu,Lv,Lw,Lu_old,Lv_old,Lw_old

    ! Particle variables
    real,save,dimension(npart) :: xpart,ypart,zpart,upart,vpart,wpart,swirl_part
    real :: xp,yp,zp,up,vp,wp
    common/part_traj/ xpart,ypart,zpart

    real :: p_cfl,umax1
    
    real :: g, ratio, a
    common/particle_params/ g,ratio,a ! g, ratio unused

    integer,save :: substep

    ! Loop integers
    integer :: j,jj

    ! File writing variables
    character*18  :: directory_part = 'outputs/particles/'
    character*18  :: filename_part
    character*20  :: filename_part_plt

    ! -------------------------------------------------------- !
    !                      Begin Calculations                  !
    ! -------------------------------------------------------- !


    ! At first time step, initialize particles
    if ((it .eq. 0 .and. irstrt .eq. 0) .or. (it .eq. 1 .and. irstrt .ne. 0)) then
        call init_part(xpart,ypart,zpart,upart,vpart,wpart,u,v,w)
        call substep_check(substep)
    else
      
        ! Change dt for substep
        dt = dt/float(substep)
 
        ! Compute new particle location
        do j = 1,npart
            if (particle_flag .eq. 0) then ! Tracer
                call tracer_update(xpart(j),ypart(j),zpart(j),upart(j),vpart(j),wpart(j),u,v,w)
            elseif (particle_flag .eq. 1) then ! Inertial passive particle
                do jj = 1,substep
                    ! Store in temporary variable to avoid memory conflict
                    xp = xpart(j)
                    yp = ypart(j)
                    zp = zpart(j)        
                    up = upart(j)
                    vp = vpart(j)
                    wp = wpart(j)        

                    call RK4(xp,yp,zp,up,vp,wp,u,v,w,u_old,v_old,w_old,omx, &
                             omy,omz,u11,u12,u13,u21,u22,u23,u31,u32,u33,   &
                             Lu,Lv,Lw,Lu_old,Lv_old,Lw_old,umax1)

                    ! Update particle data
                    xpart(j) = xp
                    ypart(j) = yp
                    zpart(j) = zp
                    upart(j) = up
                    vpart(j) = vp
                    wpart(j) = wp
                end do
            end if
        end do

        ! Change dt back to normal
        dt = dt*float(substep)

    end if

    ! Calculate percentage of particles in a vortex and write swirlX file
    call swirl_count(xpart,ypart,zpart,swirl_part,u11,u12,u13,u21,u22,u23,u31,u32,u33)

    ! Calculate concentration and particle entrapment as a function of y
!    call swirl_trap_y(ypart,swirl_part)

    ! Write particle data
    if ((mod(it,iprnfrq) .eq. 0 .or. it .eq. 1) .and. it .ne. 0) then
        if (npart .le. 10000) then ! write particle data to ASCII file
            write(filename_part,'("part-t-",i6.6,".dat")') it
            if (it .eq. 1) then
                open(31, file = directory_part//filename_part)
                write(31,*) 'variables = "x", "y", "z", "u", "v", "w", "swirl"' 
                write(31,*) 'zone f=point t="Particles", solutiontime=', it, new_line('a')
                do j = 1,npart
                    write(31, "(*(e14.6,1x))") xpart(j), ypart(j), zpart(j), upart(j), vpart(j), wpart(j), swirl_part(j)
                end do
                close(31)
            else if ((mod(it, iprnfrq) .eq. 0) .and. (it .gt. 1)) then
                open(31, file = directory_part//filename_part)
                write(31,*) 'variables = "x", "y", "z", "u", "v", "w", "swirl"' 
                write(31,*) 'zone f=point t="Particles", solutiontime=', it, new_line('a')
                do j = 1,npart
                    write(31, "(*(e14.6,1x))") xpart(j), ypart(j), zpart(j), upart(j), vpart(j), wpart(j), swirl_part(j)
                end do
                close(31)
            end if
        else ! Write particle data to Tecplot SZL file
            print *,'writing particle SZL file'
            call write_particles_szplt(xpart,ypart,zpart,upart,vpart,wpart,swirl_part)
        end if
    end if

    ! Calculate particle cfl number and check for explosions
    p_cfl = umax1*dt/a
    print *,'max relative velocity: ',umax1
    print *,'particle cfl: ',p_cfl
    if (p_cfl .ge. 10.0) then
        print *,'Particle cfl failure!'
        stop
    end if

end subroutine

! --------------------------------------------------------------------------------- !

subroutine init_part(xp,yp,zp,up,vp,wp,u,v,w)
! Initializes particle position and velocity. Reads in particle data from setup
! file generated in preconfigure command (before code compiles)

    use grid_size

    implicit none

    ! Input variables
    real,dimension(npart) :: xp,yp,zp,up,vp,wp,tp
    real,dimension(nyp,mz,mx) :: u,v,w

    ! Loop integer
    integer :: j

    ! Frozen turbulence flag
    logical :: frozen_flag
    common /frozen/ frozen_flag
 
    ! -------------------------------------------------------- !
    !                      Begin Calculations                  !
    ! -------------------------------------------------------- !

    ! Read in data from file for initial particle positions
    open(95, file = 'setup/particles/particles.dat', status = 'old', action = 'read')
    read(95,*)
    do j = 1,npart
        read(95,*) xp(j), yp(j), zp(j), up(j), vp(j), wp(j), tp(j)

        ! Force velocity to match flowfield
        call fluid_interp(xp(j),yp(j),zp(j),u,v,w,up(j),vp(j),wp(j))
        if (frozen_flag) then
            up(j) = 0.0
        end if
    end do
    close(95)

    ! At the moment, I have the particle release time (tp) in the setup file,
    ! but I'm not doing anything with it in the code - Ryan 5/26/23


end subroutine


! --------------------------------------------------------------------------------- !
subroutine substep_check(substep)
! Checks the time step and determines if we need a smaller time step just for
! integrating the particle. If not substep = 1 and it proceeds as normal.
! Otherwise, substep is the number of sub-iterations for each fluid time step to 
! sufficiently integrate the particle trajectory.

    implicit none

    ! Input variable
    integer :: substep

    ! Calculation variables
    real :: tau

    real :: g, ratio, a, C_mu, CD_switch, Pi_f, Pi_p
    common/particle_params/ g, ratio, a, C_mu, CD_switch, Pi_f, Pi_p

    real :: pi, dt
    common/data1/ pi,dt
 
    real :: re
    common/data2/ re

    ! -------------------------------------------------------- !
    !                      Begin Calculations                  !
    ! -------------------------------------------------------- !

    ! Calculate density ratios
    Pi_f = 1.0/(ratio + 0.5)
    Pi_p = ratio/(ratio + 0.5)

    ! Substep check - Ryan 4-16-22
    if (CD_switch .eq. 0) then
        tau = 1.0
    else if (CD_switch .eq. 1) then
        tau = re*(a**2)/(6.0*Pi_f)
    else 
        tau = 2.0*re*(a**2)/(9.0*Pi_f)
    end if

    substep = ceiling(max(1.0,2.0*dt/tau))
    if (substep .gt. 1) then
        print *,'NOTE: using sub-step because flow time step is too large for particle integration | substep = ',substep
    end if


end subroutine

! --------------------------------------------------------------------------------- !

subroutine tracer_update(xp,yp,zp,up,vp,wp,u,v,w)
! This subroutine updates the tracer particle position using the fluid velocity.
! Essentially, it uses an Explicit Euler integration, which is fine since fluid
! time step is already (necessarily) less than Kolmogorov time scale

    use grid_size
   
    implicit none

    ! Input variables
    real :: xp,yp,zp,up,vp,wp
    real,dimension(nyp,mz,mx) :: u,v,w

    ! Common blocks
    real :: pi, dt
    common/data1/ pi,dt

    real      :: g, ratio, a, C_mu, CD_switch, Pi_f, Pi_p
    common/particle_params/ g, ratio, a, C_mu, CD_switch, Pi_f, Pi_p

    integer :: kwall
    common /kwallpos/ kwall

    real :: re,alpha,beta,xl,zl,yl ! alpha and beta unused
    common/data2/ re,alpha,beta,xl,zl,yl

    ! Calculation variable
    real :: Height

    ! -------------------------------------------------------- !
    !                      Begin Calculations                  !
    ! -------------------------------------------------------- !

    ! Calculate height if there are walls
    Height = ycoord(nyp - kwall) - a

    ! Update particle position
    xp = mod(xl + xp + up*dt,xl)
    yp = max(min(yp + vp*dt,Height),a+ycoord(kwall))
    zp = mod(zl + zp + wp*dt,zl)
   
    ! Update particle velocity 
    call fluid_interp(xp,yp,zp,u,v,w,up,vp,wp)


end subroutine

! --------------------------------------------------------------------------------- !

subroutine fluid_interp(xp,yp,zp,u,v,w,uf,vf,wf)
! This subroutine interpolates the fluid velocity at the particle
! position and stores the values in uf, vf, and wf.

    use grid_size

    implicit none

    ! Common blocks
    real :: pi, dt
    common/data1/ pi,dt

    real :: re,alpha,beta,xl,zl,yl ! alpha and beta unused
    common/data2/ re,alpha,beta,xl,zl,yl

    ! Input variables
    real :: xp,yp,zp
    real,dimension(nyp,mz,mx) :: u,v,w

    ! Interpolation variables
    integer :: imin, imax, jmin, jmax, kmin, kmax
    real    :: xmin, xmax, ymin, ymax, zmin, zmax
    real    :: delxm, delzm, interpolate3

    ! Output variable
    real :: uf,vf,wf

    ! -------------------------------------------------------- !
    !                      Begin Calculations                  !
    ! -------------------------------------------------------- !


    ! Find nearest grid nodes with fluid data
    delxm = xl /float(mx-1)
    delzm = zl /float(mz-1)

    imin = mod(floor(xp/delxm) + mx, mx) + 1
    imax = mod(imin,mx) + 1
    jmin = mod(floor(zp/delzm) + mz, mz) + 1
    jmax = mod(jmin,mz) + 1
    kmin = ceiling(float(ny)/pi*dacos(1.0 - 2.0*yp/yl))
    if (kmin .eq. nyp) kmin = ny
    kmax = kmin + 1
    
    xmin = float(imin-1)*delxm
    xmax = float(imax-1)*delxm
    ymin = ycoord(kmin)
    ymax = ycoord(kmax)
    zmin = float(jmin-1)*delzm
    zmax = float(jmax-1)*delzm

if (imin .lt. 1 .or. jmin .lt. 1 .or. kmin .lt. 1) then
    print *,'imin/imax: ',imin,imax
    print *,'jmin/jmax: ',jmin,jmax
    print *,'kmin/kmax: ',kmin,kmax
    print *,'xp,yp,zp: ',xp,yp,zp
    stop
end if
    ! Interpolate fluid velocity at particle position
    uf = interpolate3(xp, yp, zp,       &
                  xmin, xmax, ymin, ymax, zmin, zmax,   &
                  u(kmin,jmin,imin), u(kmin,jmin,imax), &
                  u(kmin,jmax,imin), u(kmin,jmax,imax), &
                  u(kmax,jmin,imin), u(kmax,jmin,imax), &
                  u(kmax,jmax,imin), u(kmax,jmax,imax))
    vf = interpolate3(xp, yp, zp,           &
                  xmin, xmax, ymin, ymax, zmin, zmax,   &
                  v(kmin,jmin,imin), v(kmin,jmin,imax), &
                  v(kmin,jmax,imin), v(kmin,jmax,imax), &
                  v(kmax,jmin,imin), v(kmax,jmin,imax), &
                  v(kmax,jmax,imin), v(kmax,jmax,imax))
    wf = interpolate3(xp, yp, zp,           &
                  xmin, xmax, ymin, ymax, zmin, zmax,   &
                  w(kmin,jmin,imin), w(kmin,jmin,imax), &
                  w(kmin,jmax,imin), w(kmin,jmax,imax), &
                  w(kmax,jmin,imin), w(kmax,jmin,imax), &
                  w(kmax,jmax,imin), w(kmax,jmax,imax))

end subroutine

! --------------------------------------------------------------------------------- !

subroutine RK4(xpart,ypart,zpart,up,vp,wp,u,v,w,u_old,v_old,w_old,omx,omy,omz, &
               u11,u12,u13,u21,u22,u23,u31,u32,u33,Lu,Lv,Lw,Lu_old,Lv_old,Lw_old,umax1)
! Integrates a single particle at position (xp,yp,zp) with velocity (up,vp,wp)
! and the given fluid variables using RK4 to integrate MR equations.

    use grid_size

    implicit none

    ! Input variables
    real :: xpart,ypart,zpart,up,vp,wp
    real,dimension(nyp,mz,mx) :: u,v,w,u_old,v_old,w_old,omx,omy,omz
    real,dimension(nyp,mz,mx) :: u11,u12,u13,u21,u22,u23,u31,u32,u33
    real,dimension(nyp,mz,mx) :: Lu,Lv,Lw,Lu_old,Lv_old,Lw_old
    real :: umax1

    ! Common Blocks
    integer :: kwall
    common /kwallpos/ kwall

    real :: pi, dt
    common/data1/ pi,dt

    real :: re,alpha,beta,xl,zl,yl ! alpha and beta unused
    common/data2/ re,alpha,beta,xl,zl,yl

    real :: g, ratio, a
    common/particle_params/ g, ratio, a

    ! Calculation variables
    real :: Height,umax

    real :: xp,yp,zp,vp_new

    real :: Ku1,Kv1,Kw1
    real :: Ku2,Kv2,Kw2
    real :: Ku3,Kv3,Kw3
    real :: Ku4,Kv4,Kw4

    real :: u1,u2,u3,u4
    real :: v1,v2,v3,v4
    real :: w1,w2,w3,w4
    
    real :: un,vn,wn,Ku,Kv,Kw
    real :: ufluid,vfluid,wfluid
    
    ! -------------------------------------------------------- !
    !                      Begin Calculations                  !
    ! -------------------------------------------------------- !

    ! Set height of top wall
    Height = ycoord(nyp - kwall) - a


    xp = xpart
    yp = ypart
    zp = zpart
    un = up
    vn = vp
    wn = wp 
    call  RKstage(1,un,vn,wn,xp,yp,zp,Ku,Kv,Kw,u,v,w,  &
                  u_old,v_old,w_old,omx,omy,omz,       &
                  u11,u12,u13,u21,u22,u23,u31,u32,u33, &
                  Lu,Lv,Lw,Lu_old,Lv_old,Lw_old,       &
                  xpart,ypart,zpart,up,vp,wp)
     u2 = un
     v2 = vn
     w2 = wn
    Ku1 = Ku
    Kv1 = Kv
    Kw1 = Kw
    call  RKstage(2,un,vn,wn,xp,yp,zp,Ku,Kv,Kw,u,v,w,  &
                  u_old,v_old,w_old,omx,omy,omz,       &
                  u11,u12,u13,u21,u22,u23,u31,u32,u33, &
                  Lu,Lv,Lw,Lu_old,Lv_old,Lw_old,       &
                  xpart,ypart,zpart,up,vp,wp)
     u3 = un
     v3 = vn
     w3 = wn
    Ku2 = Ku
    Kv2 = Kv
    Kw2 = Kw
    call  RKstage(3,un,vn,wn,xp,yp,zp,Ku,Kv,Kw,u,v,w,  &
                  u_old,v_old,w_old,omx,omy,omz,       &
                  u11,u12,u13,u21,u22,u23,u31,u32,u33, &
                  Lu,Lv,Lw,Lu_old,Lv_old,Lw_old,       &
                  xpart,ypart,zpart,up,vp,wp)
     u4 = un
     v4 = vn
     w4 = wn
    Ku3 = Ku
    Kv3 = Kv
    Kw3 = Kw
    call  RKstage(4,un,vn,wn,xp,yp,zp,Ku,Kv,Kw,u,v,w,  &
                  u_old,v_old,w_old,omx,omy,omz,       &
                  u11,u12,u13,u21,u22,u23,u31,u32,u33, &
                  Lu,Lv,Lw,Lu_old,Lv_old,Lw_old,       &
                  xpart,ypart,zpart,up,vp,wp)
    Ku4 = Ku
    Kv4 = Kv
    Kw4 = Kw
 

    ! Update position
    xpart = mod(xl + xpart + dt*(up/6 + u2/3 + u3/3 + u4/6),xl)
    ypart = max(min(ypart + dt*(vp/6 + v2/3 + v3/3 + v4/6),Height),a+ycoord(kwall))
    zpart = mod(zl + zpart + dt*(wp/6 + w2/3 + w3/3 + w4/6),zl)

    ! Update velocity 
    up = up + dt*(Ku1/6 + Ku2/3 + Ku3/3 + Ku4/6)
    vp_new = vp + dt*(Kv1/6 + Kv2/3 + Kv3/3 + Kv4/6)
    wp = wp + dt*(Kw1/6 + Kw2/3 + Kw3/3 + Kw4/6)

    ! Handle particles near walls
    if ((ypart .eq. Height .and. vp_new .gt. 0.0) .or. (ypart .eq. a+ycoord(kwall) .and. vp_new .lt. 0.0)) then
        vp = -vp_new ! Bounce elastically off the walls  
    else
        vp = vp_new
    end if


    ! Find max particle relative velocity
    call fluid_interp(xpart,ypart,zpart,u,v,w,ufluid,vfluid,wfluid)

    umax = sqrt((ufluid-up)**2 + (vfluid-vp)**2 + (wfluid-wp)**2)
    if (umax .gt. umax1) then
        umax1 = umax
    end if

end subroutine

! --------------------------------------------------------------------------------- !

subroutine RKstage(sn,un,vn,wn,xp,yp,zp,Ku,Kv,Kw,u,v,w, &
                   u_old,v_old,w_old,omx,omy,omz,       &
                   u11,u12,u13,u21,u22,u23,u31,u32,u33, &
                   Lu,Lv,Lw,Lu_old,Lv_old,Lw_old,       &
                   xp_old,yp_old,zp_old,up_old,vp_old,wp_old)

! Used for each stage (determined by sn) to compute updated velocity and
! acceleration based on MR equations

    use grid_size
    implicit none
 
    ! Common blocks
    real :: pi, dt
    common/data1/ pi, dt

    real :: re, alpha, beta, xl, yl, zl
    common/data2/ re, alpha, beta, xl, zl, yl

    integer :: kwall
    common /kwallpos/ kwall

    real :: g, ratio, a, C_mu, CD_switch, Pi_f, Pi_p
    common/particle_params/ g, ratio, a, C_mu, CD_switch, Pi_f, Pi_p

    logical :: frozen_flag
    common/frozen/ frozen_flag

    ! Input variables 
    real, dimension(nyp,mz,mx) :: u,v,w,omx,omy,omz,u_old,v_old,w_old
    real, dimension(nyp,mz,mx) :: u11,u12,u13,u21,u22,u23,u31,u32,u33
    real, dimension(nyp,mz,mx) :: Lu,Lv,Lw,Lu_old,Lv_old,Lw_old
    real :: xp,yp,zp,un,vn,wn,Ku,Kv,Kw
    integer :: sn

    ! Calculation variables
    real :: C_D, Re_p, C_L, Height
    real :: h, umax
    real :: subDufDt, subDvfDt, subDwfDt
    real :: ufluid, vfluid, wfluid, wxf, wyf, wzf, Luf, Lvf, Lwf
    real :: PGX, DX, LX, PGY, DY, LY, PGZ, DZ, LZ, FX, FY, FZ
    real :: up_old, vp_old, wp_old, xp_old, yp_old, zp_old
    real :: dLufdt, dLvfdt, dLwfdt

    ! -------------------------------------------------------- !
    !                      Begin Calculations                  !
    ! -------------------------------------------------------- !

    ! For now, set C_L to 1/2
    C_L = 0.5
 
    ! Set height of top wall
    Height = ycoord(nyp - kwall) - a

    ! Calculate fluid velocity and vorticity at particle location
    call fluid_interp(xp,yp,zp,u,v,w,ufluid,vfluid,wfluid)
    call fluid_interp(xp,yp,zp,omx,omy,omz,wxf,wyf,wzf)
    call fluid_interp(xp,yp,zp,Lu,Lv,Lw,Luf,Lvf,Lwf)

    ! Calculate substantial derivative at particle location
    call SubDerivative(u,v,w,u_old,v_old,w_old,ufluid,vfluid,wfluid,subDufDt,subDvfDt,subDwfDt, &
                       xp,yp,zp,u11,u12,u13,u21,u22,u23,u31,u32,u33,Lu,Lv,Lw,Lu_old,Lv_old,     &
                       Lw_old,dLufdt,dLvfdt,dLwfdt)
    
    Re_p = 2*a*re*sqrt((ufluid - un)**2 + (vfluid - vn)**2 + (wfluid - wn)**2)

    ! Drag model
    if (CD_switch .eq. 0) then
        C_D = 0.0
    elseif (CD_switch .eq. 1) then
        C_D = 16/(2*a*re)*C_mu
    elseif (CD_switch .eq. 2) then
        C_D = 24/(2*a*re)*(1 + 0.15*Re_p**(0.687))
    elseif (CD_switch .eq. 3) then
       if (Re_p .gt. 1.0) then
           C_D = 16.31*Re_p**(-0.8331)
       else
           C_D = 17.28*Re_p**(-0.9711)
       end if
    end if

    ! Calculate each force component
    if (frozen_flag .eq. .true.) then
        un = ufluid
    end if

    ! x forces
    PGX = 1.5*subDufDt
    DX = C_D*3/(8*a)*(ufluid - un + a**2/6.0*Luf)
    LX = C_L*((vfluid - vn)*wzf - (wfluid - wn)*wyf)
    FX = -a**2/20.0*dLufdt ! Faxen correction term

    ! y forces    
    PGY = 1.5*subDvfDt
    DY = C_D*3/(8*a)*(vfluid - vn + a**2/6.0*Lvf)
    LY = C_L*((wfluid - wn)*wxf - (ufluid - un)*wzf)
    FY = -a**2/20.0*dLvfdt ! Faxen correction term

    ! z forces
    PGZ = 1.5*subDwfDt
    DZ = C_D*3/(8*a)*(wfluid - wn + a**2/6.0*Lwf)
    LZ = C_L*((ufluid - un)*wyf - (vfluid - vn)*wxf)
    FZ = -a**2/20.0*dLwfdt ! Faxen correction term

    ! Combine to find acceleration component
    Ku = (Pi_f*(PGX + DX + LX + FX))
    Kv = (Pi_f*(PGY + DY + LY + FY) - g*(Pi_f - Pi_p))
    Kw = (Pi_f*(PGZ + DZ + LZ + FZ))

    if (sn .eq. 3) then
        h = dt
    else 
        h = dt/2.0
    end if

  
    un = up_old + h*Ku
    vn = vp_old + h*Kv
    wn = wp_old + h*Kw
    
    xp = mod(xl + xp_old + h*un,xl)
    yp = max(min(yp_old + h*vn, Height),a+ycoord(kwall))
    zp = mod(zl + zp_old + h*wn,zl)

!    print *,'yp_old = ',yp_old
!    print *,'vn = ',vn
!    print *,'h = ',h
!    print *,'Height = ',Height
!    print *,'yp = ',yp
!    print *,''
!    print *,'yp = max(min(yp_old + h*vn, Height),a+ycoord(kwall))'

end subroutine 

! --------------------------------------------------------------------------------- !

subroutine SubDerivative(u,v,w,u_old,v_old,w_old,ufluid,vfluid,wfluid,subDufDt,subDvfDt,subDwfDt, &
                         xpart,ypart,zpart,u11,u12,u13,u21,u22,u23,u31,u32,u33,Lu,Lv,Lw,Lu_old,   &
                         Lv_old,Lw_old,dLufdt,dLvfdt,dLwfdt)
    use grid_size

    implicit none

    ! Input variables
    real,dimension(nyp,mz,mx) :: u,v,w,u_old,v_old,w_old
    real,dimension(nyp,mz,mx) :: u11,u12,u13,u21,u22,u23,u31,u32,u33
    real,dimension(nyp,mz,mx) :: Lu,Lv,Lw,Lu_old,Lv_old,Lw_old

    real :: ufluid,vfluid,wfluid
    real :: subDufDt,subDvfDt,subDwfDt
    real :: dLufdt,dLvfdt,dLwfdt

    real :: xpart,ypart,zpart

    ! Common blocks
    real :: pi, dt
    common/data1/ pi, dt

    real :: re, alpha, beta, xl, yl, zl
    common/data2/ re, alpha, beta, xl, zl, yl

    ! Interpolation variables
    integer :: imin, imax, jmin, jmax, kmin, kmax
    real    :: xmin, xmax, ymin, ymax, zmin, zmax
    real    :: delxm, delzm, interpolate3

    ! Loop integers
    integer :: n_grid

    ! Calculation variables
    integer,dimension(8) :: gridi, gridj, gridk
    real,dimension(3,8)  :: dudt,dLudt
    real,dimension(3)    :: dufdt

    real :: dufdx, dvfdx, dwfdx
    real :: dufdy, dvfdy, dwfdy
    real :: dufdz, dvfdz, dwfdz


    ! -------------------------------------------------------- !
    !                      Begin Calculations                  !
    ! -------------------------------------------------------- !

    delxm = xl /float(mx-1)
    delzm = zl /float(mz-1)
 
    imin = mod(floor(xpart/delxm) + mx, mx) + 1
    imax = mod(imin,mx) + 1
    jmin = mod(floor(zpart/delzm) + mz, mz) + 1
    jmax = mod(jmin,mz) + 1
    kmin = ceiling(float(ny)/pi*dacos(1.0 - 2.0*ypart/yl))
    if (kmin .eq. nyp) kmin = ny
    kmax = kmin + 1
    
    xmin = float(imin-1)*delxm
    xmax = float(imax-1)*delxm
    ymin = ycoord(kmin)
    ymax = ycoord(kmax)
    zmin = float(jmin-1)*delzm
    zmax = float(jmax-1)*delzm


    gridk = [kmin, kmin, kmin, kmin, kmax, kmax, kmax, kmax] ! y index
    gridj = [jmin, jmin, jmax, jmax, jmin, jmin, jmax, jmax] ! z index
    gridi = [imin, imax, imin, imax, imin, imax, imin, imax] ! x index

    ! Time derivative
    do n_grid = 1,8
        dudt(1,n_grid) = (u(gridk(n_grid),gridj(n_grid),gridi(n_grid)) - u_old(gridk(n_grid),gridj(n_grid),gridi(n_grid)))/dt 
        dudt(2,n_grid) = (v(gridk(n_grid),gridj(n_grid),gridi(n_grid)) - v_old(gridk(n_grid),gridj(n_grid),gridi(n_grid)))/dt 
        dudt(3,n_grid) = (w(gridk(n_grid),gridj(n_grid),gridi(n_grid)) - w_old(gridk(n_grid),gridj(n_grid),gridi(n_grid)))/dt 

        dLudt(1,n_grid) = (Lu(gridk(n_grid),gridj(n_grid),gridi(n_grid)) - Lu_old(gridk(n_grid),gridj(n_grid),gridi(n_grid)))/dt 
        dLudt(2,n_grid) = (Lv(gridk(n_grid),gridj(n_grid),gridi(n_grid)) - Lv_old(gridk(n_grid),gridj(n_grid),gridi(n_grid)))/dt 
        dLudt(3,n_grid) = (Lw(gridk(n_grid),gridj(n_grid),gridi(n_grid)) - Lw_old(gridk(n_grid),gridj(n_grid),gridi(n_grid)))/dt 
    end do

    ! Interpolating substantial derivative at particle position
    dufdx = interpolate3(xpart,ypart,zpart,xmin,xmax,ymin,ymax,zmin,zmax, &
                         u11(kmin,jmin,imin), u11(kmin,jmin,imax), &
                         u11(kmin,jmax,imin), u11(kmin,jmax,imax), &
                         u11(kmax,jmin,imin), u11(kmax,jmin,imax), &
                         u11(kmax,jmax,imin), u11(kmax,jmax,imax))
    dufdy = interpolate3(xpart,ypart,zpart,xmin,xmax,ymin,ymax,zmin,zmax, &
                         u12(kmin,jmin,imin), u12(kmin,jmin,imax), &
                         u12(kmin,jmax,imin), u12(kmin,jmax,imax), &
                         u12(kmax,jmin,imin), u12(kmax,jmin,imax), &
                         u12(kmax,jmax,imin), u12(kmax,jmax,imax))
    dufdz = interpolate3(xpart,ypart,zpart,xmin,xmax,ymin,ymax,zmin,zmax, &
                         u13(kmin,jmin,imin), u13(kmin,jmin,imax), &
                         u13(kmin,jmax,imin), u13(kmin,jmax,imax), &
                         u13(kmax,jmin,imin), u13(kmax,jmin,imax), &
                         u13(kmax,jmax,imin), u13(kmax,jmax,imax))
    dvfdx = interpolate3(xpart,ypart,zpart,xmin,xmax,ymin,ymax,zmin,zmax, &
                         u21(kmin,jmin,imin), u21(kmin,jmin,imax), &
                         u21(kmin,jmax,imin), u21(kmin,jmax,imax), &
                         u21(kmax,jmin,imin), u21(kmax,jmin,imax), &
                         u21(kmax,jmax,imin), u21(kmax,jmax,imax))
    dvfdy = interpolate3(xpart,ypart,zpart,xmin,xmax,ymin,ymax,zmin,zmax, &
                         u22(kmin,jmin,imin), u22(kmin,jmin,imax), &
                         u22(kmin,jmax,imin), u22(kmin,jmax,imax), &
                         u22(kmax,jmin,imin), u22(kmax,jmin,imax), &
                         u22(kmax,jmax,imin), u22(kmax,jmax,imax))
    dvfdz = interpolate3(xpart,ypart,zpart,xmin,xmax,ymin,ymax,zmin,zmax, &
                         u23(kmin,jmin,imin), u23(kmin,jmin,imax), &
                         u23(kmin,jmax,imin), u23(kmin,jmax,imax), &
                         u23(kmax,jmin,imin), u23(kmax,jmin,imax), &
                         u23(kmax,jmax,imin), u23(kmax,jmax,imax))
    dwfdx = interpolate3(xpart,ypart,zpart,xmin,xmax,ymin,ymax,zmin,zmax, &
                         u31(kmin,jmin,imin), u31(kmin,jmin,imax), &
                         u31(kmin,jmax,imin), u31(kmin,jmax,imax), &
                         u31(kmax,jmin,imin), u31(kmax,jmin,imax), &
                         u31(kmax,jmax,imin), u31(kmax,jmax,imax))
    dwfdy = interpolate3(xpart,ypart,zpart,xmin,xmax,ymin,ymax,zmin,zmax, &
                         u32(kmin,jmin,imin), u32(kmin,jmin,imax), &
                         u32(kmin,jmax,imin), u32(kmin,jmax,imax), &
                         u32(kmax,jmin,imin), u32(kmax,jmin,imax), &
                         u32(kmax,jmax,imin), u32(kmax,jmax,imax))
    dwfdz = interpolate3(xpart,ypart,zpart,xmin,xmax,ymin,ymax,zmin,zmax, &
                         u33(kmin,jmin,imin), u33(kmin,jmin,imax), &
                         u33(kmin,jmax,imin), u33(kmin,jmax,imax), &
                         u33(kmax,jmin,imin), u33(kmax,jmin,imax), &
                         u33(kmax,jmax,imin), u33(kmax,jmax,imax))

    ! Interpolating time derivatives at particle position
    dufdt(1) = interpolate3(xpart, ypart, zpart,               &
                            xmin, xmax, ymin, ymax, zmin, zmax,         &
                            dudt(1,1), dudt(1,2), dudt(1,3), dudt(1,4), &
                            dudt(1,5), dudt(1,6), dudt(1,7), dudt(1,8))
 
    dufdt(2) = interpolate3(xpart, ypart, zpart,               &
                            xmin, xmax, ymin, ymax, zmin, zmax,         &
                            dudt(2,1), dudt(2,2), dudt(2,3), dudt(2,4), &
                            dudt(2,5), dudt(2,6), dudt(2,7), dudt(2,8))
    
    dufdt(3) = interpolate3(xpart, ypart, zpart,               &
                            xmin, xmax, ymin, ymax, zmin, zmax,         &
                            dudt(3,1), dudt(3,2), dudt(3,3), dudt(3,4), &
                            dudt(3,5), dudt(3,6), dudt(3,7), dudt(3,8))
             
    dLufdt   = interpolate3(xpart, ypart, zpart,               &
                            xmin, xmax, ymin, ymax, zmin, zmax,         &
                            dLudt(1,1), dLudt(1,2), dLudt(1,3), dLudt(1,4), &
                            dLudt(1,5), dLudt(1,6), dLudt(1,7), dLudt(1,8))
 
    dLvfdt   = interpolate3(xpart, ypart, zpart,               &
                            xmin, xmax, ymin, ymax, zmin, zmax,         &
                            dLudt(2,1), dLudt(2,2), dLudt(2,3), dLudt(2,4), &
                            dLudt(2,5), dLudt(2,6), dLudt(2,7), dLudt(2,8))
    
    dLwfdt   = interpolate3(xpart, ypart, zpart,               &
                            xmin, xmax, ymin, ymax, zmin, zmax,         &
                            dLudt(3,1), dLudt(3,2), dLudt(3,3), dLudt(3,4), &
                            dLudt(3,5), dLudt(3,6), dLudt(3,7), dLudt(3,8))

    ! Calculate substantial derivative of fluid
    subDufDt = dufdt(1) + ufluid*dufdx + vfluid*dufdy + wfluid*dufdz
    subDvfDt = dufdt(2) + ufluid*dvfdx + vfluid*dvfdy + wfluid*dvfdz
    subDwfDt = dufdt(3) + ufluid*dwfdx + vfluid*dwfdy + wfluid*dwfdz

end subroutine

! --------------------------------------------------------------------------------- !

subroutine swirl_count(xpart,ypart,zpart,swirl_part,u11,u12,u13,u21,u22,u23,u31,u32,u33)
! Computes swirl at each particle position and then computes percentage of
! particles inside a vortex using a given threshold value

    use grid_size

    implicit none

    ! Input variables
    real,dimension(npart) :: xpart,ypart,zpart,swirl_part
    real,dimension(nyp,mz,mx) :: u11,u12,u13,u21,u22,u23,u31,u32,u33

    ! Common blocks
    integer :: it
    common/pre5f/ it

    ! Calculation variables
    integer :: swirl_cnt, j

    ! -------------------------------------------------------- !
    !                      Begin Calculations                  !
    ! -------------------------------------------------------- !

    ! Compute percentage of particles in vortex
    swirl_cnt = 0;
    do j = 1,npart
        ! Compute swirl at particle locations
        call interp_swirl(xpart(j),ypart(j),zpart(j),swirl_part(j),u11,u12,u13,u21,u22,u23,u31,u32,u33)
        if (swirl_part(j) .gt. 10) then
            swirl_cnt = swirl_cnt + 1
        end if
    end do

    ! Write output to data file
    if (it .eq. 1) then
        open(30, file = "outputs/swirl10")
        write(30,*) float(swirl_cnt)/float(npart)*100
    else
        open(30, file = "outputs/swirl10", position = "append")
        write(30,*) float(swirl_cnt)/float(npart)*100
    end if
    close(30)

end subroutine

! --------------------------------------------------------------------------------- !

! Trilinear interpolation function
real function interpolate3(x, y, z, x0, x1, y0, y1, z0, z1, c000, c100, c001, c101, c010, c110, c011, c111)
    ! This function uses trilinear interpolation to estimate the value
    ! of a function f at a given point (x,y,z)
    ! Reference: http://en.wikipedia.org/wiki/Trilinear_interpolation
    implicit none
    real :: x, y, z, x0, x1, y0, y1, z0, z1
    real :: xd, yd, zd
    real :: c000, c100, c001, c101, c010, c110, c011, c111
    real :: c00, c01, c10, c11
    real :: c0, c1
 
    if (x1 .gt. x0) then
        xd = (x - x0)/(x1 - x0)
    else
        xd = 0.0
    end if

    if (y1 .gt. y0) then
        yd = (y - y0)/(y1 - y0)
    else
        yd = 0.0
    end if

    if (z1 .gt. z0) then
        zd = (z - z0)/(z1 - z0)
    else
        zd = 0.0
    end if

    c00 = c000*(1.0 - xd) + c100*xd
    c01 = c001*(1.0 - xd) + c101*xd
    c10 = c010*(1.0 - xd) + c110*xd
    c11 = c011*(1.0 - xd) + c111*xd

    c0 = c00*(1.0 - yd) + c10*yd
    c1 = c01*(1.0 - yd) + c11*yd
    
    interpolate3 = c0*(1.0 - zd) + c1*zd
end function interpolate3


! --------------------------------------------------------------------------------- !


subroutine interp_swirl(xpart,ypart,zpart,swirl,u11,u12,u13,u21,u22,u23,u31,u32,u33)

    use grid_size

    implicit none


    ! Input variables
    real :: xpart, ypart, zpart
    real, dimension(nyp,mz,mx) :: u,v,w,u11,u12,u13,u21,u22,u23,u31,u32,u33

    ! Common blocks:
    real :: pi, dt
    common/data1/ pi, dt

    real :: re, alpha, beta, xl, yl, zl
    common/data2/ re, alpha, beta, xl, zl, yl

    integer :: it
    common/pre5f/ it

    ! Interoplation variables
    integer :: imin, imax, jmin, jmax, kmin, kmax
    real    :: xmin, xmax, ymin, ymax, zmin, zmax
    real    :: delxm, delzm, interpolate3

    real    :: swirl
    real    :: s1,s2,s3,s4,s5,s6,s7,s8


    ! -------------------------------------------------------- !
    !                      Begin Calculations                  !
    ! -------------------------------------------------------- !


    delxm = xl /float(mx-1)
    delzm = zl /float(mz-1)
 

    ! Normal interpolation
    imin = mod(floor(xpart/delxm) + mx, mx) + 1
    imax = mod(imin,mx) + 1
    jmin = mod(floor(zpart/delzm) + mz, mz) + 1
    jmax = mod(jmin,mz) + 1
    kmin = ceiling(float(ny)/pi*dacos(1.0 - 2.0*ypart/yl))
    if (kmin .eq. nyp) kmin = ny
    kmax = kmin + 1
    
    xmin = float(imin-1)*delxm
    xmax = float(imax-1)*delxm
    ymin = ycoord(kmin)
    ymax = ycoord(kmax)
    zmin = float(jmin-1)*delzm
    zmax = float(jmax-1)*delzm

    ! Calculate swirl at nearest points
    call calcswirl(u11(kmin,jmin,imin),u21(kmin,jmin,imin),u31(kmin,jmin,imin), &
                   u12(kmin,jmin,imin),u22(kmin,jmin,imin),u32(kmin,jmin,imin), &
                   u13(kmin,jmin,imin),u23(kmin,jmin,imin),u33(kmin,jmin,imin),s1)
    call calcswirl(u11(kmin,jmin,imax),u21(kmin,jmin,imax),u31(kmin,jmin,imax), &
                   u12(kmin,jmin,imax),u22(kmin,jmin,imax),u32(kmin,jmin,imax), &
                   u13(kmin,jmin,imax),u23(kmin,jmin,imax),u33(kmin,jmin,imax),s2)
    call calcswirl(u11(kmin,jmax,imin),u21(kmin,jmax,imin),u31(kmin,jmax,imin), &
                   u12(kmin,jmax,imin),u22(kmin,jmax,imin),u32(kmin,jmax,imin), &
                   u13(kmin,jmax,imin),u23(kmin,jmax,imin),u33(kmin,jmax,imin),s3)
    call calcswirl(u11(kmin,jmax,imax),u21(kmin,jmax,imax),u31(kmin,jmax,imax), &
                   u12(kmin,jmax,imax),u22(kmin,jmax,imax),u32(kmin,jmax,imax), &
                   u13(kmin,jmax,imax),u23(kmin,jmax,imax),u33(kmin,jmax,imax),s4)
    call calcswirl(u11(kmax,jmin,imin),u21(kmax,jmin,imin),u31(kmax,jmin,imin), &
                   u12(kmax,jmin,imin),u22(kmax,jmin,imin),u32(kmax,jmin,imin), &
                   u13(kmax,jmin,imin),u23(kmax,jmin,imin),u33(kmax,jmin,imin),s5)
    call calcswirl(u11(kmax,jmin,imax),u21(kmax,jmin,imax),u31(kmax,jmin,imax), &
                   u12(kmax,jmin,imax),u22(kmax,jmin,imax),u32(kmax,jmin,imax), &
                   u13(kmax,jmin,imax),u23(kmax,jmin,imax),u33(kmax,jmin,imax),s6)
    call calcswirl(u11(kmax,jmax,imin),u21(kmax,jmax,imin),u31(kmax,jmax,imin), &
                   u12(kmax,jmax,imin),u22(kmax,jmax,imin),u32(kmax,jmax,imin), &
                   u13(kmax,jmax,imin),u23(kmax,jmax,imin),u33(kmax,jmax,imin),s7)
    call calcswirl(u11(kmax,jmax,imax),u21(kmax,jmax,imax),u31(kmax,jmax,imax), &
                   u12(kmax,jmax,imax),u22(kmax,jmax,imax),u32(kmax,jmax,imax), &
                   u13(kmax,jmax,imax),u23(kmax,jmax,imax),u33(kmax,jmax,imax),s8)

    ! Interpolating swirl strength at particle position
    swirl = interpolate3(xpart,ypart,zpart,xmin,xmax,ymin,ymax,zmin,zmax,s1,s2,s3,s4,s5,s6,s7,s8)


end subroutine 


! --------------------------------------------------------------------------------- !


subroutine swirl_trap_y(ypart,swirl_part)

    use grid_size

    implicit none

    ! Input variable
    real,dimension(npart) :: ypart,swirl_part

    ! Common blocks
    integer :: it
    common/pre5f/ it

    integer :: irstrt, nsteps, iprnfrq ! irstrt and nsteps unused
    common/iocntrl/ irstrt, nsteps, iprnfrq

    real :: re, alpha, beta, xl, yl, zl
    common/data2/ re, alpha, beta, xl, zl, yl

    real :: dPdx
    common/pressure/ dPdx

    ! Calculation variables
    integer :: j,k
    integer,dimension(18) :: C,D
    real :: u_tau,yp_plus


    ! -------------------------------------------------------- !
    !                      Begin Calculations                  !
    ! -------------------------------------------------------- !


    ! Define friction velocity
    u_tau = sqrt(-dPdx) ! Technically this should be sqrt(-dPdx/rho), but we don't have rho in the code, and it's just 1 for water

    ! Initialize concentration vector
    C(:) = 0 ! # of particles total
    D(:) = 0 ! # of particles inside vortex

    ! Count number of particles in each y-plane 
    do j = 1,npart
        if (ypart(j) .le. yl/2.0) then
            yp_plus = ypart(j)*u_tau*re
        else
            yp_plus = (yl - ypart(j))*u_tau*re
        end if

        do k = 1,18
            if (yp_plus .gt. float(k-1)*10.0 .and. yp_plus .le. float(k)*10.0) then
                C(k) = C(k) + 1
                if (swirl_part(j) .gt. 10.0) then
                    D(k) = D(k) + 1
                end if
            end if
        end do

!        if (yp_plus .le. 5) then
!            C(1) = C(1) + 1 ! viscous sublayer
!            if (swirl_part(j) .gt. 10.0) then
!                D(1) = D(1) + 1
!            end if
!        else if (yp_plus .gt. 5 .and. yp_plus .le. 30) then
!            C(2) = C(2) + 1 ! buffer region
!            if (swirl_part(j) .gt. 10.0) then
!                D(2) = D(2) + 1
!            end if
!        else if (yp_plus .gt. 30 .and. yp_plus .le. 100) then
!            C(3) = C(3) + 1 ! inner layer
!            if (swirl_part(j) .gt. 10.0) then
!                D(3) = D(3) + 1
!            end if
!        else if (yp_plus .gt. 100) then
!            C(4) = C(4) + 1 ! log-law
!            if (swirl_part(j) .gt. 10.0) then
!                D(4) = D(4) + 1
!            end if
!        end if
    end do


    ! Write concentration to a file each time step
    if (it .eq. 1) then
        open(300, file = "outputs/part_concentration")
        write(300,"(*(f8.4,1x))") float(C)/float(npart)*100
    else if (mod(it,iprnfrq) .eq. 0 .and. it .ne. 0) then
        open(300, file = "outputs/part_concentration", position = "append")
        write(300,"(*(f8.4,1x))") float(C)/float(npart)*100
    end if
    close(300)

    ! Write trap percentage to a file each time step
    if (it .eq. 1) then
        open(400, file = "outputs/trap_y")
        write(400,"(*(f8.4,1x))") float(D)/float(C)*100
    else if (mod(it,iprnfrq) .eq. 0 .and. it .ne. 0) then
        open(400, file = "outputs/trap_y", position = "append")
        write(400,"(*(f8.4,1x))") float(D)/float(C)*100
    end if
    close(400)

end subroutine


! --------------------------------------------------------------------------------- !


#ifdef OUTPUTFORM
subroutine write_particles_szplt(xp,yp,zp,up,vp,wp,swirl)

    use iso_c_binding 
    use grid_size
    implicit none

    include 'tecio.f90'

    integer :: it, i
    real, dimension(npart) :: xp,yp,zp,up,vp,wp,swirl
    common/pre5f/ it

    character*1 :: NULLCHR
    character (37) :: filename

    integer(c_int64_t) :: numValues
    integer(c_int32_t) :: fileFormat = 1 ! 0 = .plt | 1 = .szplt
    integer(c_int32_t) :: gridFileType = 1 ! GRID
    integer(c_int32_t) :: outputFileType = 0 ! Grid + Solution
    integer(c_int32_t) :: gridZone = 1
    integer(c_int32_t) :: outputZone = 1
    integer(c_int32_t) :: zoneType = 0 ! ZoneType: 0-ordered
    integer(c_int32_t), allocatable :: varTypes(:)
    integer(c_int32_t), allocatable :: shareVarFromZone(:)
    integer(c_int32_t), allocatable :: valueLocation(:)
    integer(c_int32_t), allocatable :: passiveVarList(:)
    integer(c_int32_t) :: shareFaceNeighborsFromZone = 0
    integer(c_int32_t) :: numFaceConnections = 0
    integer(c_int32_t) :: faceNeighborMode = 0
    integer(c_int32_t) :: totalNumFaceNodes = 0
    integer(c_int32_t) :: numConnectedBoundaryFaces = 0
    integer(c_int32_t) :: totalNumBoundaryConnections = 0
    integer(c_int32_t) :: shareConnectivityFromZone = 0
    integer(c_int32_t) :: isDouble = 0
    real(c_float), allocatable :: floatValues(:)
    real(c_float), allocatable :: coordinates(:,:,:,:)
    type(c_ptr) :: gridFileHandle = C_NULL_PTR
    type(c_ptr) :: outputFileHandle = C_NULL_PTR

    NULLCHR = CHAR(0)

    ! Initialize helper arrays
    numValues = npart
    allocate(floatValues(numValues))
    allocate(varTypes(numValues))
    allocate(shareVarFromZone(numValues))
    allocate(valueLocation(numValues))
    allocate(passiveVarList(numValues))

    varTypes = 1
    shareVarFromZone = 0
    valueLocation = 1
    passiveVarList = 0

    ! If this is the first time step, write the grid solution
    ! Open files
    write(filename,'("outputs/particles/part-t-",i6.6,".szplt")')it
    i = tecini142("Solution", "x y z u v w swirl", filename // NULLCHR, ".", &
        fileFormat, outputFileType, 1, 0)

    ! Create zone
    i = teczne142( &
        "Particles", &
        zoneType, &
        npart, &
        1, &
        1, &
        0, 0, 0, &
        1.0 * it, &
        1, & ! StrandID
        0, &
        1, &
        numFaceConnections, &
        faceNeighborMode, &
        totalNumFaceNodes, &
        numConnectedBoundaryFaces, &
        totalNumBoundaryConnections, &
        passiveVarList, &
        valueLocation, &
        shareVarFromZone, &
        shareConnectivityFromZone &
    )


    ! Write x
    floatValues = pack(xp, .true.)
    i = tecdat142(numValues, floatValues, isDouble)

    ! Write y
    floatValues = pack(yp, .true.)
    i = tecdat142(numValues, floatValues, isDouble)

    ! Write z
    floatValues = pack(zp, .true.)
    i = tecdat142(numValues, floatValues, isDouble)

    ! Write u
    floatValues = pack(up, .true.)
    i = tecdat142(numValues, floatValues, isDouble)

    ! Write v
    floatValues = pack(vp, .true.)
    i = tecdat142(numValues, floatValues, isDouble)

    ! Write w
    floatValues = pack(wp, .true.)
    i = tecdat142(numValues, floatValues, isDouble)

    ! Write swirl
    floatValues = pack(swirl, .true.)
    i = tecdat142(numValues, floatValues, isDouble)

    ! Close file
    i = tecend142()


end subroutine

#else

subroutine write_particles_szplt(xp,yp,zp,up,vp,wp,swirl)

    use grid_size
    implicit none

    real,dimension(npart) :: xp,yp,zp,up,vp,wp,swirl

end subroutine
#endif
