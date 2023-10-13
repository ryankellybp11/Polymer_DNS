!      *** PURPOSE ***
! 
!      THE PURPOSE OF THIS PROGRAM IS TO SOLVE A 3D VERSION OF THE
!      NAVIER-STOKES EQNS. AS RECAST IN A 4TH ORDER EQN. SYSTEM FOR
!      THE VERTICAL VELOCITY. THIS SYSTEM IS DESCRIBED IN JFM V177
!      PP 138-140.  THE CURRENT APPLICATIONS ARE TO FLOW GEOMETRIES
!      WHEREIN THE VERTICAL VELOCITY IS ZERO AT BOTH BOUNDARIES NORMAL
!      TO THE CHEB. DIRECTION.  THE STREAMWISE VEL. MUST EITHER BE
!      CONSTANT AT THESE BOUNDARIES OR ITS NORMAL GRADIENT MUST BE
!      CONSTANT (CONSTANT SHEAR).  BY CONTINUITY THESE CONDITIONS
!      IMPLY THE FOLLOWING BEHAVIOR FOR NORMAL DERIVATIVES OF THE
!      VERTICAL VELOCITY AT THESE BOUNDARIES.
! 
!         dV/dy = 0 : U = constant on boundary all X
! 
!         d(dV/dy)/dy = 0 : dU/dy = constant in X at boundary
! 
!      WITH MINOR MODIFICATIONS THE CURRENT FORMULATION CAN BE MADE
!      TO SATISFY CONDITIONS ON THE STREAMWISE VELOCITY SUCH THAT
!      THESE nth ORDER V DERIVATIVES ARE NON-ZERO AND FUNCTIONS OF
!      X.  IF THE FLOW GEOMETRY REQUIRES NON-ZERO VERTICAL VELOCITY
!      AT EITHER OF THE CHEB. BNDRYS THE GREEN'S FCT. FORMULATION
!      MUST BE RECODED.
!      
!      THE VERTICAL VELOCITY IS CONSTRUCTED FROM A COMPOSITE SOLUTION
!     
!         V = Vp + C1*V1 + C2*V2
! 
!      WHERE Vp IS THE SOLUTION TO THE NON-HOMOGENEOUS AND NON-LINEAR
!      PART OF THE DISCRETE MODEL (TIME-STEPPED) AND V1,V2 ARE DETERMINED
!      FROM THE LINEAR HOMOGENEOUS PART OF THE MODEL DURING THE INITIALIZATION.
!      VANISHING OF THE VERTICAL VEL. AT BNDRYS IS REQUIRED BY THE IMPOSITION
!      OF HOMOGENEOUS DIRICLET BCS ON THESE VARIABLES IN THE MAIN DRIVE LOOP
!      AND IN SBR'S GFNC1 & GFNC2. 
! 
!      THE DERIVATIVE BCS ON V ARE SATISFIED BY EVALUATING THE (nt)th AND (nb)th
!      DERIVATIVES OF THE ABOVE EQUATION AT THE TOP AND BOTTOM BNDRYS
!      RESPECTIVELY AND INVERTING FOR C1,C2.
! 
!         C1*D= d(nb)V2/dy(nb)*[+d(nt)V/dy(nt)-d(nt)Vp/dy(nt)]
!              +d(nt)V2/dy(nt)*[-d(nb)V/dy(nb)+d(nb)Vp/dy(nb)]
! 
!         C2*D=-d(nb)V1/dy(nb)*[+d(nt)V/dy(nt)-d(nt)Vp/dy(nt)]
!              +d(nt)V1/dy(nt)*[-d(nb)V/dy(nb)+d(nb)Vp/dy(nb)]
! 
!            D= d(nt)V1/dy(nt)*d(nb)V2/dy(nb)-d(nt)V2/dy(nt)*d(nb)V1/dy(nb)
! 
!      THE ABOVE DETERMINATION IS REQUIRED AT EACH WAVE NO. IN X.
!      THE MAIN DRIVE LOOP BELOW ASSUMES THE DERIVATIVE BCS ON
!      V AT BOUNDARIES ARE HOMOGENEOUS NEUMANN.  IF THE FLOW
!      PHYSICS REQUIRE NON-VANISHING DERIVATIVES THE DO 20 LOOP
!      BELOW SHOULD BE ALTERED TO INCLUDE THE TERMS OMMITED FROM
!      THE ABOVE EQNS.  FURTHER MODS. TO INSURE THE CONSISTENT
!      Uo AND P CALCULATIONS WILL BE NECESSARY IN SBR'S UZERO AND
!      PRESS.
!      Domain has x as the streamwise direction, z as the spanwise
!      direction and y as the height.
! -----------------------------------------------------------------------

program threed
    use omp_lib
    use grid_size
  
    implicit none
  
    complex u(nyp,nz,nxh), v(nyp,nz,nxh), w(nyp,nz,nxh)
    complex omx(nyp,nz,nxh), omy(nyp,nz,nxh), omz(nyp,nz,nxh)
    complex fn(nyp,nz,nxh), fnm1(nyp,nz,nxh)
    complex gn(nyp,nz,nxh), gnm1(nyp,nz,nxh)
    complex wrkc(nyp,nz,nxh)
    complex im
    complex gf1(nyp,nz,nxh),gf2(nyp,nz,nxh)
    complex c1(nz,nxh),c2(nz,nxh)
    complex bctop(nz,nxh),bcbot(nz,nxh)
    complex dnv1b(nz,nxh),dnv1t(nz,nxh)
    complex dnv2b(nz,nxh),dnv2t(nz,nxh),denom(nz,nxh)
    
    ! Velocity gradient tensor
    complex u11(nyp,nz,nxh), u12(nyp,nz,nxh), u13(nyp,nz,nxh)
    complex u21(nyp,nz,nxh), u22(nyp,nz,nxh), u23(nyp,nz,nxh)
    complex u31(nyp,nz,nxh), u32(nyp,nz,nxh), u33(nyp,nz,nxh)
    complex,dimension(nyp,nz,nxh) :: Lu,Lv,Lw ! Laplacian terms - Ryan 7/24/23

    real ampbdy1,ampbdy2,c11amp,tfrac
    real u0(nyp),h1n(nyp),h1nm1(nyp),f
    real w0(nyp),h3n(nyp),h3nm1(nyp)
    real umn(nyp),urms(nyp),vrms(nyp),wrms(nyp),prms(nyp)
    real a(nyhp,nz,3),t(nyp)
    real wrk1(nyp,nz,nx)
    real wfft1(nmax), wfft2(nmax)
    real wfft3(nmax), wfft4(nmax)
    real w1(mx,mz), w2(mx,mz)
    real uchbeng(nyp), vchbeng(nyp), wchbeng(nyp)
    real trigx(2*nx), trigz(2*nz), trigy(2*ny), sine(ny), cosine(ny)
    real pi,dt,theta,wavz(nz),wavx(nxh),c(nyp),yc(nyp)
    real re,alpha,beta,xl,zl,yl,time,dyde,xstart,uinf
    real e1(nz,nxh),ckx(nz,nxh),tm1
    real seght(nyp),spread(mz2,mx2),fxintg(nyp,mz,mx),fzintg(nyp,mz,mx), fyintg(nyp,mz,mx)
    real gain,ugain
    real wn2,x,g,ysmth,zsmth,sums,p1t,p1b,p2t,p2b
    ! new BC variables - Ryan 6/8/23
    real :: atu,btu,gtu,abu,bbu,gbu
    real :: atw,btw,gtw,abw,bbw,gbw
  
    integer irstrt,nsteps,iprnfrq,nt,nb,ixfax(19),izfax(19),iyfax(19)
    integer ib,it,i,j,k,iyt,iyb,is,jj
    integer print3d

    common/preforce/ seght
    common/pre2forc/ spread,fxintg
    common/pre4f/ fzintg, fyintg
    common/pre5f/ it
    common/params/ gain, ugain
    common/data1/ pi,dt,theta,wavz,wavx,c,yc
    common/data2/ re,alpha,beta,xl,zl,yl,time,dyde,xstart,uinf
    common/energ/ e1,ckx,tm1
    common/iocntrl/ irstrt,nsteps,iprnfrq,print3d
    common/u0bcs/ atu,btu,gtu,abu,bbu,gbu
    common/w0bcs/ atw,btw,gtw,abw,bbw,gbw ! adding separate BCs for w
    common/grnfct/ nt,nb,p1t,p1b,p2t,p2b
    common/trigx/ trigx, ixfax
    common/trigz/ trigz, izfax
    common/trigy/ trigy,sine,cosine,iyfax
  
    integer bfhead            ! index # for the buffer zone head face.
    integer bftail            ! index # for the buffer zone tail face.
    real    bfgain(1200)          ! buffer zone stiffness distribution array.
    real    bfugain(1200)         ! buffer zone damping coef. distribution array.
    integer bfwidth
    common /buffer/ bfgain, bfugain, bfhead, bftail, bfwidth
  
    real,dimension(nyp,nz,nx) :: initu, initv, initw
    common /init/ initu, initv, initw
  
    integer kwall, kmaxsurf
    common /kwallpos/ kwall, kmaxsurf
  
    integer geomtype, flow_select
    common /geometry/ geomtype, flow_select
  !............................................................................
  
    ib = 0
   
  !      ib    = switch to implement various type of bc's; 
  !              ib=0 then dirichlet and t(i)=1
  !              ib=1 then neumann and t(i)=k*k
  !               ip = 0, even;ip = 1, odd.
  !      see subroutine  solve()
  
  
  !      call sub. setstuf to fill common blocks data1,data2,iocntrl
  !      and idiag data described in setstuf 
    call setstuf

  !  if(nmax .lt. (4*nx*nyp) .or. nmax .lt. (4*nx*nz) .or. nmax .lt. (4*nyp*nz)) then
  !     write(*,*) 'error error error, nmax too small'
  !     stop
  !  end if
  
  !   zero out some stuff
   
    do i=1,nyp
       umn(i)  = 0.0
       urms(i) = 0.0
       vrms(i) = 0.0
       wrms(i) = 0.0
       prms(i) = 0.0
    end do
   
   
   
  !     determine the sol'n to the first homogeneous system. sub. gfcn1 computes
  !     v1 and returns the soln in gf1.  the derivatives necessary for computation
  !     of the constant c1 are computed.
  !     subroutine gfcn2 performs the same for the second homogeneous system.

    call gfcn (gf1,a,t,wrkc,wrk1,bctop,bcbot,dnv1b,dnv1t,nb,nt,p1b,p1t)   
    call gfcn (gf2,a,t,wrkc,wrk1,bctop,bcbot,dnv2b,dnv2t,nb,nt,p2b,p2t)
  
  
  !     calculate the common denominator for the constants c1(kz,kx),c2(kz,kx)
  !$omp parallel do 
   do k=1,nxh
       do j=1,nz
          denom(j,k)=dnv1t(j,k)*dnv2b(j,k)-dnv2t(j,k)*dnv1b(j,k)
       end do
    end do
  !$omp end parallel do
   
   
  ! --  initialize the fluid. sub. initial supplies initial values
  !     of phi,fn,fnm1,u0,h1n,h1nm1.  it also supplies the forcing
  !     term f.
    call forcepre
   
    it = 0

    call initial(u,u0,v,w,w0,omx,omy,omz,fn,fnm1,gn,gnm1,h1n,h1nm1,h3n,h3nm1,  &
                 wrk1,wrkc,wfft1,wfft2,wfft3,wfft4,w1,w2,bctop,bcbot,a,t,      &
                 u11,u12,u13,u21,u22,u23,u31,u32,u33,Lu,Lv,Lw)

  ! -----------------------------------------------------------------
  ! ------------------------  main loop  ----------------------------
  ! -----------------------------------------------------------------
   
    do it = 1, nsteps
  
       
       write(*,22) it
  22   format('  it = ',i7)


!     initial phi field is computed from the laplacian of the v field.
       
       call cderiv(v,wrk1)
       call cderiv(wrk1,wrkc)
  !$omp parallel do     
       do k=1,nxh
          do j=1,nz
             wn2 = wavx(k)**2 +wavz(j)**2
             do i=1,nyp
                v(i,j,k)=wrkc(i,j,k) -wn2*v(i,j,k)
             end do
          end do
       end do
  !$omp end parallel do
    
  ! --  evaluate rhs of time discrete equation for particular phi
  ! --  phirhs also sets bcbot=bctop=(0,0).
       
       call phirhs(v,wrkc,wrk1,fn,fnm1,bcbot,bctop)
   
  !   copy fn to fnm1
   
       do k=1,nxh
          do j=1,nz
             do i=1,nyp
                fnm1(i,j,k) = fn(i,j,k)
             end do
          end do
       end do
        
       ib = 0
       x = 0.0              ! parameter for x (or kx) dependency
       g = re/(dt*theta)    ! parameter for timestep dep., theta=.5 cn method
  
  !this loop cannot be parallelized
       do k=1,nxh
          x = wavx(k)**2
          call solve(  v(1,1,k), x,g,dyde,bctop(1,k),bcbot(1,k),ib,wavz,c,a,t,wrkc)
       end do
  
  !     the rhs of the poisson eqn for vp is phi
   
  !      solve for vp
  !      the time stepping term is g=0
         
       g = 0.0
  
       ib = 0
  
  !this loop cannot be parallelized
       do k=1,nxh
          x = wavx(k)**2
          call solve(  v(1,1,k),x,g,dyde,bctop(1,k),bcbot(1,k),ib,wavz,c,a,t,wrkc)
       end do
  
  ! --  evaluate rhs of time discrete equation for vorticity
  ! --  omyrhs also sets the boundary conditions which may have 
  ! --  been altered during pressure sol'n.
       
       call omyrhs(omy,wrkc,wrk1,gn,gnm1,bcbot,bctop)
   
  !   copy gn to gnm1
       do k=1,nxh
          do j=1,nz
             do i=1,nyp
                gnm1(i,j,k) = gn(i,j,k)
             end do
          end do
       end do
   
  !      evaluate the derivs. necessary for construction of the composite
  !      solution. all derivatives of the green's fcts have already been
  !      calculated in sbr's gfcn1 & gfcn2 & the commom denominator has
  !      been formed above in the "do 5" loop.
  ! 
  !          v=vp+c1*gf1+c2*gf2
  ! 
  !          
  !      evaluate d(nb)vp/dy(nb) at y=-1
         
       iyb=2
       iyt=4
       if (nb .eq. 1) then
          call c1derbw(v,wrkc,iyb)
          if (nt .eq. 1) then
             call c1dertw(v,wrkc,iyt)
          else if (nt .eq. 2) then
             call cderiv(v,wrk1)
             call c1dertw(wrk1,wrkc,iyt)
          end if
  
       else if (nb .eq. 2) then
          call cderiv(v,wrk1)
          call c1derbw(wrk1,wrkc,iyb)
          if (nt .eq. 1) then
             call c1dertw(v,wrkc,iyt)
          else if (nt .eq. 2) then
             call cderiv(v,wrk1)
             call c1dertw(wrk1,wrkc,iyt)
          end if
       end if
   
  !      form total v
  !$omp parallel do 
       do k=1,nxh
          do j=1,nz
             c1(j,k)=(wrkc(iyb,j,k)*dnv2t(j,k)-wrkc(iyt,j,k)*dnv2b(j,k))/denom(j,k)
             c2(j,k)=(wrkc(iyt,j,k)*dnv1b(j,k)-wrkc(iyb,j,k)*dnv1t(j,k))/denom(j,k)
          end do
       end do
  !$omp end parallel do
  
  !$omp parallel do 
       do k=1,nxh
          do j=1,nz
             do i=1,nyp
                v(i,j,k)=v(i,j,k)+c1(j,k)*gf1(i,j,k)+c2(j,k)*gf2(i,j,k)
             end do
          end do
       end do
  !$omp end parallel do
  
  !    now do vorticity field
  
       x = 0.0              ! parameter for x (or kx) dependency
       g = re/(dt*theta)    ! parameter for timestep dep., theta=.5 cn method
       ib = 0
  
  !not parallizable
       do k=1,nxh
          x = wavx(k)**2
          call penta(omy(1,1,k),x,g,dyde,ib,wavz,c,atu,btu,bctop(1,k),abu,bbu,bcbot(1,k),a,t,wrk1)
       end do
  
  !  because of periodicity, the 0,0 mode for omy is always 
  !  zero.
  
       do i=1,nyp
          omy(i,1,1) = (0.0,0.0)
       end do
  
  !   solve for the zeroth fourier mode of streamwise velocity
  
       call uzero(u0,h1n,h1nm1,a,t,wrk1)
       call wzero(w0,h3n,h3nm1,a,t,wrk1)
  
       do i=1,nyp
          h1nm1(i) = h1n(i)
          h3nm1(i) = h3n(i)
       end do
  
  !     CALCULATE THE REMAINING COMPONENTS OF VELOCITY FROM CONTINUITY 
  !     AND NORMAL VORTICITY IN SPECTRAL SPACE
  !    
  !     U(N,Kz,Kx) = 
  !          -((i*Kx)(-dV/dy)+(i*Kz)(Omy))/(Kx**2+Kz**2) : Kx**2+Kz**2 > 0 
  !     for Kx**2 + Kz**2 =0 HAVE U0 ABOVE AS REAL PART OF U(I,KX=0).
  !
  !     W(N<Kz,Kx) 
  !          -((i*Kz)(-dV/dy)-(i*Kx)(Omy))/(Kx**2+Kz**2) : Kx**2+Kz**2 > 0 
  !     For Kx**2 + Kz**2 = 0, W(n,Kz,Kx) = 0.0
  !
  !     THE  IMAGINARY PART OF THE Kx=Kz=0  MODES
  !     ARE USED TO STORE THE REAL PART OF U(I,KMAX),W(I,KMAX) WHICH IS SET
  !     TO ZERO IN SUB. UVEL.  
  !     FOR CONSISTENCY WITH CONTINUITY SUBROUTINE  VELOC ALSO RESETS
  !     THE IMAG. COMPONENT OF V(I,KX=0) (I.E. THE REAL PART OF KMAX MODE) 
  !     EQUAL ZERO    
      
       call norm(v) 
       call norm(omy)
       call veloc(u,v,w,u0,w0,omy,wrkc)
 

  !    normalize to ensure conjg sym.
      
       call norm(u)      
       call norm(w)
 
  !    calculate the other vorticity components.
      
       call vort(u,v,w,omx,omz,wrkc)
       call norm(omx)
       call norm(omz)

! Calculate derivates needed for conformation tensor sources
! Keep this for velocity gradient tensor
        call derivscji(u,v,w,u11,u12,u13,u21,u22,u23,u31,u32,u33,Lu,Lv,Lw)

        call norm(u11)
        call norm(u12)
        call norm(u13)
        call norm(u21)
        call norm(u22)
        call norm(u23)
        call norm(u31)
        call norm(u32)
        call norm(u33)
    
        ! Laplacian terms - Ryan 7/24/23
        call norm(Lu)
        call norm(Lv)
        call norm(Lw)

  !    at frequencies specified by iprnfrq calculate the pressure
  !    and output u,v,w, and p and determine statistics.    
  
   
  !  transform all data into y-physical
       
       is = 1


  !$omp parallel sections default(shared) private(wfft1,wfft2,wfft3,wfft4)      
  !$omp section      
       call yfft(u,wfft1,wfft2,wfft3,wfft4,is)
  !$omp section
       call yfft(v,wfft1,wfft2,wfft3,wfft4,is)
  !$omp section   
       call yfft(w,wfft1,wfft2,wfft3,wfft4,is)
  !$omp section   
       call yfft(omx,wfft1,wfft2,wfft3,wfft4,is)
  !$omp section   
       call yfft(omy,wfft1,wfft2,wfft3,wfft4,is)
  !$omp section   
       call yfft(omz,wfft1,wfft2,wfft3,wfft4,is)
 
! -------------------------------------------------------------------------------- !
!                   Adding the other ffts - Ryan 2-18-22                           !
! -------------------------------------------------------------------------------- !

  !$omp section   
      call yfft(u11,wfft1,wfft2,wfft3,wfft4,is)    
  !$omp section   
      call yfft(u12,wfft1,wfft2,wfft3,wfft4,is)
  !$omp section   
      call yfft(u13,wfft1,wfft2,wfft3,wfft4,is)
  !$omp section   
      call yfft(u21,wfft1,wfft2,wfft3,wfft4,is)
  !$omp section   
      call yfft(u22,wfft1,wfft2,wfft3,wfft4,is)
  !$omp section   
      call yfft(u23,wfft1,wfft2,wfft3,wfft4,is)
  !$omp section   
      call yfft(u31,wfft1,wfft2,wfft3,wfft4,is)
  !$omp section   
      call yfft(u32,wfft1,wfft2,wfft3,wfft4,is)
  !$omp section   
      call yfft(u33,wfft1,wfft2,wfft3,wfft4,is)
    
  ! Laplacian terms - Ryan 7/24/23
  !$omp section   
      call yfft(Lu,wfft1,wfft2,wfft3,wfft4,is)
  !$omp section   
      call yfft(Lv,wfft1,wfft2,wfft3,wfft4,is)
  !$omp section   
      call yfft(Lw,wfft1,wfft2,wfft3,wfft4,is)
  !$omp end parallel sections


      call vcw3d(u,v,w,omx,omy,omz,fn,gn,                 &
                 u11,u12,u13,u21,u22,u23,u31,u32,u33,Lu,Lv,Lw)

  
  write(*,*) '1'
  !   transform all data into y-spectral
       
       is = -1
 

  !$omp parallel sections default(shared) private(wfft1,wfft2,wfft3,wfft4)   
  !$omp section     
       call yfft(v,wfft1,wfft2,wfft3,wfft4,is)
  !$omp section     
       call yfft(omy,wfft1,wfft2,wfft3,wfft4,is)
  !$omp section     
       call yfft(fn,wfft1,wfft2,wfft3,wfft4,is)
  !$omp section     
       call yfft(gn,wfft1,wfft2,wfft3,wfft4,is)
  !$omp section     
       call yfft(omz,wfft1,wfft2,wfft3,wfft4,is)
  !$omp end parallel sections

       call norm(v)
       call norm(omy)
       call norm(fn)
       call norm(gn)
       call norm(omz)


  write(*,*) '2'
  !......................................................................
  !
  ! use the smoothing prof gottlieb suggested, assuming the force is in
  ! fn, gn and omz.  * Should I add something for the scalar here? Ryan 3-18-22

  if (flow_select .ge. 0) then  
  !.... spectral smoothing in y-direction
  !$omp parallel do         
       do k = 1, nyp
          ysmth = exp(-3.0 * ((float(k-1)/float(ny))**10))
          do  j = 1, nz
             do i = 1, nxh
                omz(k,j,i) = omz(k,j,i) * ysmth
                gn(k,j,i)  = gn(k,j,i) * ysmth
                fn(k,j,i)  = fn(k,j,i) * ysmth
             end do
          end do
       end do
  !$omp end parallel do
  !.... spectral smoothing in z-direction
  !$omp parallel do         
       do j = 1, nz
          jj = j - 1
          if(j .gt. nzh) jj = nz - j + 1
          zsmth = exp( -3.0 * ((float(jj)/float(nzh))**10))
          do k = 1, nyp
             do i = 1, nxh
                omz(k,j,i) = omz(k,j,i) * zsmth
                gn(k,j,i)  = gn(k,j,i) * zsmth
                fn(k,j,i)  = fn(k,j,i) * zsmth
             end do
          end do
       end do
  !$omp end parallel do
  end if

  write(*,*) '3'
  !......................................................................
  !
  !        
       do i=1,nyp
          h1n(i)=real(gn(i,1,1))
          h3n(i)=real(omz(i,1,1))
       end do
  
  ! on output fn = h2 = (u x om)2
  !           gn = h1 = (u x om)1
  !           omz= h3 = (u x om)3
  
  !     subroutine  nonlin calculates the term fn in the fourth order system 
  !     for the vertical velocity, and the term gn for the vertical vorticity 
  !     equation.
  
       call norm(fn)
       call norm(gn)
       call norm(omz)
  
       call nonlin(fn,omz,gn,wrkc,wrk1)
  
       call norm(fn)
       call norm(gn)
  

! ------------------------------------------------------------------------ !

  !     update the running time
  
       time=time+dt
  write(*,*) '4',it,nsteps
       if(mod(it,iprnfrq) .eq. 0) then
  
          open(122, file = 'outputs/last-restart', status = 'replace', form = 'unformatted')   ! alex 09/14/2018
          write(122) time
          write(122) v,fn,fnm1,omy,gn,gnm1,u0,h1n,h1nm1,w0,h3n,h3nm1
          write(122) fxintg, fyintg, fzintg
          rewind(122)
          close(122)

       end if
      write(*,*) 'Time step: ', it, ' of ', nsteps
    end do
  
  !-----------------   end main loop -------------------------------
       
  !    calculate averages
  
    sums=float(nsteps)/(float(iprnfrq))
    
!   open(112, file = 'restart', form = 'unformatted')   ! alex 09/14/2018
    write(*,112) nsteps,iprnfrq, sums
  112 format(2x,2i5,e16.9)
!   close(112)                           ! alex 09/14/2018
  
  
  !$omp parallel do
    do i=1,nyp
       umn(i)=umn(i)/sums
       urms(i)=sqrt(urms(i))/sums
       prms(i)=sqrt(prms(i))/sums
       vrms(i)=sqrt(vrms(i))/sums
       wrms(i)=sqrt(wrms(i))/sums
    end do
  !$omp end parallel do

end program threed

!---------------------------------------------------------------
!---------------------------------------------------------------

subroutine setstuf 
    use grid_size

    implicit none
  
    real trigx(2*nx), trigz(2*nz), trigy(2*ny), sine(ny), cosine(ny)
    real pi,dt,theta,wavz(nz),wavx(nxh),c(nyp),yc(nyp)
    real re,alpha,beta,xl,zl,yl,time,dyde,xstart,uinf
    real e1(nz,nxh),ckx(nz,nxh),tm1
    real seght(nyp),spread(mz2,mx2),fxintg(nyp,mz,mx),fzintg(nyp,mz,mx), fyintg(nyp,mz,mx)
    real gain,ugain
    real wn2,x,g,ysmth,zsmth,sums,p1t,p1b,p2t,p2b
    ! New BC variables - Ryan 6/8/23
    real :: atu,btu,gtu,abu,bbu,gbu
    real :: atw,btw,gtw,abw,bbw,gbw
    real trigx32(16000),trigz32(4000)
    real bfgain, bfugain
    real rn,rj,yj

    integer irstrt,nsteps,iprnfrq,nt,nb,ixfax(19),izfax(19),iyfax(19),ixfax32(19),izfax32(19)
    integer ib,it,i,j,k,iyt,iyb,is,jj, imatrix(nyp,mz,mx), bfhead, bftail, bfwidth, nyi, nzi, nxi
    integer print3d
    common/preforce/ seght
    common/pre2forc/ spread,fxintg
    common/pre4f/ fzintg, fyintg
    common/pre5f/ it
    common/params/ gain, ugain
    common/data1/ pi,dt,theta,wavz,wavx,c,yc
    common/data2/ re,alpha,beta,xl,zl,yl,time,dyde,xstart,uinf
    common/energ/ e1,ckx,tm1
    common/iocntrl/ irstrt,nsteps,iprnfrq,print3d
    common/u0bcs/ atu,btu,gtu,abu,bbu,gbu
    common/w0bcs/ atw,btw,gtw,abw,bbw,gbw ! adding separate BCs for w
    common/grnfct/ nt,nb,p1t,p1b,p2t,p2b
    common/trigx/ trigx, ixfax
    common/trigz/ trigz, izfax
    common/trigy/ trigy,sine,cosine,iyfax
    common/trigx32/trigx32,ixfax32  !use in 3/2 rule calcs
    common/trigz32/trigz32,izfax32  !use in 3/2 rule calcs
    common/geommat/ imatrix
  
    common /buffer/ bfgain, bfugain, bfhead, bftail, bfwidth
  
    integer perturbtime
    common /perturbation/ perturbtime
  
    integer kwall, kmaxsurf
    common /kwallpos/ kwall, kmaxsurf
  
    integer geomtype, flow_select
    common /geometry/ geomtype, flow_select
  
    real geomwidth, geomheight
    common /geomprops/ geomwidth, geomheight
  
    integer readvdes
    real vdes(mx), xsuction
    common /suction/ vdes, xsuction, readvdes
  
    integer particle_flag
    common /particle/ particle_flag
    real gravity, ratio, a, C_mu, CD_switch
    common/particle_params/ gravity, ratio, a, C_mu, CD_switch

    real bdyfx,L,rad,xcenter,ycenter,zcenter,ampbdy1,ampbdy2,c11amp,tfrac 
    common/vortRing/ bdyfx,L,rad,xcenter,ycenter,zcenter,ampbdy1,ampbdy2,c11amp,tfrac 

    real forbeta
    common/data22/ forbeta
  
    real dPdx, R_tau
    common/pressure/ dPdx, R_tau

    integer nlines, io

    real yfrac, a_spacing, r_spacing
  
  !    first set some stuff - note that grid densities must be coded into 
  !    the parameter statements.
  !
  !    set pi to something - suggest approximately---
  !    pi=3.14159 26535 89793 23846 26433 832795 028841 971693 9937511
  !    unless it proves difficult to achieve good results with this value
  
    pi= 2.0*acos(0.0)
  
  !    set the running time to 0. this will be changed as the fluid
  !    is initialized if a resstart package is used.
    
    time = 0.0
  
  
  
    open(110, file='setup/dns.config', status='old')     
  
    !    read the dns setup data
    read(110,*) 
    read(110,*) 
    read(110,*) 
    read(110,*) nsteps    !--- nsteps -- the total no. of time steps to integrate eqns. during run
    read(110,*) iprnfrq   !--- iprnfrq -- freq. to dump solns for analysis. freq. to calc. pressure
    read(110,*) dt        !--- dt -- the time step. time = time + dt at each step.
    read(110,*) gain   
    read(110,*) ugain  
    read(110,*) theta     !--- theta -- implicitness factor for 2nd ord. time accuracy. (.5 for cn)
  
    !   read data characterizing geometry and fluid
    read(110,*) nyi         !--- nyi -- order of the cheb poly in the nonhom. direction: npts=nyi+1
    read(110,*) nzi         !--- nzi -- no. of grid points in periodic direction z
    read(110,*) nxi         !--- nxi -- no. of grid points in periodic direction x
    read(110,*)
    read(110,*) yl          !--- yl -- half-length in nonhomogeneous direction
    read(110,*) zl          !--- zl -- length in periodic direction
    read(110,*) xl          !--- xl -- length in periodic direction
    read(110,*) re          !--- re -- reynolds number
    read(110,*) xstart
    read(110,*) uinf
    read(110,*) kwall
  
    write(*,*) "kwall is equal to ", kwall
    write(*,*) "xstart is equal to ", xstart
  
    !    check consistency with parameter statements.
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
    read(110,*) 
    read(110,*) 
    read(110,*) gravity
    read(110,*) 

    read(110,*) 
    read(110,*) 
    read(110,*) 
    read(110,*) 
    read(110,*) 

    ! Read in selection flags/switches 
    read(110,*) irstrt    !--- irstrt -- ne 0 means a restart file is used to initialize this run
    read(110,*) 
    read(110,*) print3d   !--- print3d -- switch how we print the flowfield
    read(110,*) geomtype
    read(110,*) readvdes
    read(110,*) particle_flag
    read(110,*) flow_select
    read(110,*) 
    read(110,*) 
    read(110,*) 
    read(110,*) 

    read(110,*) 
    read(110,*) 
    read(110,*) 
    read(110,*) 
    read(110,*) 

    ! Read in data for particles
    read(110,*) 
    read(110,*) 
    read(110,*) ratio
    read(110,*) a
    read(110,*) C_mu
    read(110,*) CD_switch

    !gravity = grav ! Ryan 7/14/22 - There's probably a better way to do this, but this was easiest to have the same gravity in two pre-defined common blocks

    read(110,*) 
    read(110,*) 
    read(110,*) 
    read(110,*) 
    read(110,*) 

    ! Read in data for scalar and polymer
    read(110,*) 
    read(110,*) 
    read(110,*) 
    read(110,*) 
    read(110,*) 
    read(110,*) 
    read(110,*) 
    read(110,*) 
    read(110,*) 
    read(110,*) 
    read(110,*) 
    read(110,*) 
    read(110,*) 
    read(110,*) 
    read(110,*) 

!    if (src_stop .lt. 1) then
!        src_stop = nsteps
!    end if

    read(110,*) 
    read(110,*) 
    read(110,*) 
    read(110,*) 
    read(110,*) 

    ! Read in data for body force terms (used for vortex ring generation)
    read(110,*) ampbdy1, ampbdy2, c11amp
    read(110,*) tfrac
    read(110,*) forbeta
    read(110,*) xcenter, ycenter, zcenter
    read(110,*) L, rad
    close(110)

    bdyfx = ampbdy1 ! Ryan 7/14/22

    ! Set y coordinate
    do k = 1, nyp
        ycoord(k) = (1. - cos(float(k - 1)*pi/float(ny)))*(yl / 2.0)  ! [0,YL]
    end do
  
  !    initialize the fft package
  
    call rcsexp(nx,ixfax,trigx)
    call ccosexp(ny,sine,cosine,iyfax,trigy)
    call ccexp(nz,izfax,trigz)
    call rcsexp(nx32,ixfax32,trigx32)
    call ccexp(nz32,izfax32,trigz32)
  
    ! Read in the geometry matrix.
    open(111, file = 'code/bin/geometry/geometry', form = 'unformatted')  ! alex 09/14/2018
    read(111) imatrix
    close(111)
  
  !---------------------------------------------------------------
  
  !    calculate the resolvable wave nos. in x: assumes length xl has
  !    been non-dimensionalized with the length yl.
  
    alpha=2.*pi/xl
    beta=2.*pi/zl
  
    do k=1,nxh
       wavx(k) = float(k-1)*alpha   
    end do
    do j=1,nz/2
       wavz(j) = float(j-1)*beta 
    end do
    do j=nz/2+1,nz
       wavz(j) = float(j-nz-1)*beta
    end do
  
  !    define coefficients ck for cheb. expansions
  
    do i = 1,nyp
       c(i) = 1.
    end do
    c(1) = 2.0
    c(nyp) = 2.0
  
  
  !    calculate dy/d(eta) scaling factor: see sub. solve: assumes
  !    physical problem is solved on a space 0>x2>-yl and 
  !    y=(2*x2+yl)/yl.  this scale factor is applied in the
  !    poisson solver routines and in the routines evaluating derivatives
  !    in the chebyshev direction.
  
    dyde = 2.0/yl
  
  !    calculate phys. coords. in cheb. dir. for later use.
  
    rn = float(ny)
    do j=1,nyp
       rj=float(j-1)
       yj=cos(pi*rj/rn)
       yc(j)=(yj-1.)*yl/2.
    end do
 

    ! Set boundary conditions based on flow - Ryan 6/8/23
    call setbcs ! Only manipulates common block variables so no need for args 

end subroutine setstuf

! Ryan 6/8/23
subroutine setbcs
! This subroutine sets appropriate BC variables based on the flow selection
! switch instead of having to manually define it in the setup file each time. I
! added separate variables to handle u and v (I'm not sure why they were the
! same by default in the first place).
!
! One odd note is that all the "top" BCs correspond to our output of y = 0 and
! all the "bottom" BCs correspond to y = yl. This is a weird thing with how we
! apparently write our 3D data in real space after the transform, but I'm just
! rolling with it

! From the original place this was called (for context):
    !     read data for the calculation of the greens fcts phi1,phi2 in sbr's
    !     gfcn1,gfcn2.  the composite sol'n insures that v=0 on both boundaries
    !     normal to the cheb dir. by virtue of the bc's on vp,v1,v2.  specifying
    !     nt=1 or 2 will insure that either dv/dy(1) or d(dv/dy)/dy(1) vanishes
    !     at the top boundary respectively.  nb has the same effect for the bottom 
    !     boundary y=-1.  the dirichlet bc's for phi1 & phi2 are arbitrary but
    !     must be independent. guess what p1b & p1t etc. denote.
  
    !     read data for bc's on the zeroth fourier mode of streamwise velocity.
    !     general form of bc's are
  
    !     at*uo(1)+bt*duo/dy(1)=gt(1)
    !     ab*uo(-1)+bb*duo/dy(-1)=gb(-1)
  
    !     current formulation does not provide for ai,bi,gi to be time dependent
    !     this capability easily included by providing recipes prior to call to
    !     sbr uzero
  
    use grid_size

    implicit none


    ! Boundary Condition Variables
    real :: atu,btu,gtu,abu,bbu,gbu
    real :: atw,btw,gtw,abw,bbw,gbw
    integer :: nt,nb
    real :: p1t,p1b,p2t,p2b
    real :: re,alpha,beta,xl,zl,yl,time,dyde,xstart,Uinf

    common/u0bcs/ atu,btu,gtu,abu,bbu,gbu
    common/w0bcs/ atw,btw,gtw,abw,bbw,gbw ! adding separate BCs for w
    common/grnfct/ nt,nb,p1t,p1b,p2t,p2b  ! these won't change much
    common/data2/ re,alpha,beta,xl,zl,yl,time,dyde,xstart,Uinf

    ! Flow selector
    integer :: geomtype,flow_select
    common/geometry/ geomtype,flow_select


    print *,' Setting BCs...'
    print *,' flow_select = ',flow_select
    ! Setting BCs based on flow_select
    select case (flow_select)
        case (1,11) ! Channel flow
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
            gtu = 0.0

            ! w=0 on y = 0 wall
            atw = 1.0
            btw = 0.0
            gtw = 0.0

            ! u=Uinf on y = yl wall
            abu = 1.0
            bbu = 0.0
            gbu = Uinf

            ! w=0 on y = yl wall
            abw = 1.0
            bbw = 0.0
            gbw = 0.0

            ! Green's Function BCs
            nt = 1.0
            nb = 1.0

        case (3) ! 3D Couette Flow
            ! u = 0 on bottom wall, u = Uinf on top wall

            ! u=0 on y = 0 wall
            atu = 1.0
            btu = 0.0
            gtu = 0.0

            ! w=0 on y = 0 wall
            atw = 1.0
            btw = 0.0
            gtw = 0.0

            ! u=Uinf on y = yl wall
            abu = 1.0
            bbu = 0.0
            gbu = Uinf

            ! w=0 on y = yl wall
            abw = 1.0
            bbw = 0.0
            gbw = Uinf

            ! Green's Function BCs
            nt = 1.0
            nb = 1.0

        case (4) ! Blasius BL
            ! u = 0 on bottom wall, du/dy = 0 on top wall

            ! u=0 on y = 0 wall
            atu = 1.0
            btu = 0.0
            gtu = 0.0

            ! w=0 on y = 0 wall
            atw = 1.0
            btw = 0.0
            gtw = 0.0

            ! u=Uinf on y = yl wall
            abu = 0.0
            bbu = 1.0
            gbu = 0.0

            ! w=0 on y = yl wall
            abw = 0.0
            bbw = 1.0
            gbw = 0.0

            ! Green's Function BCs
            nt = 1.0
            nb = 2.0

        case (0,10,5) ! Vortex Ring or Still fluid
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
            stop
    end select

    ! For all the cases so far, these numbers have not changed, so I'm defining
    ! them outside the selection loop
 
    p1b = 1.0
    p1t = -1.0
    p2b = -1.0
    p2t = 0.5 

end subroutine setbcs

subroutine initial(u,u0,v,w,w0,omx,omy,omz,fn,fnm1,gn,gnm1,h1n,h1nm1,h3n,h3nm1,  &
                   wrk1,wrkc,wfft1,wfft2,wfft3,wfft4,w1,w2,bctop,bcbot,a,t,      &
                   u11,u12,u13,u21,u22,u23,u31,u32,u33,Lu,Lv,Lw)

  use grid_size

  complex u(nyp,nz,nxh), v(nyp,nz,nxh), w(nyp,nz,nxh)
  complex omx(nyp,nz,nxh), omy(nyp,nz,nxh), omz(nyp,nz,nxh)
  complex fn(nyp,nz,nxh), fnm1(nyp,nz,nxh)
  complex gn(nyp,nz,nxh), gnm1(nyp,nz,nxh)
  complex wrkc(nyp,nz,nxh)
  complex bctop(nz,nxh),bcbot(nz,nxh)
  complex im
  real u0(nyp),h1n(nyp),h1nm1(nyp)
  real w0(nyp),h3n(nyp),h3nm1(nyp)
  real wrk1(nyp,nz,nx)
  real  w1(1), w2(1), wfft1(1), wfft2(1)
  real wfft3(1), wfft4(1)
  real  w5(1), w6(1), w7(1), w8(1), w9(1)
  real a(nyhp,nz,3),t(nyp)

!     velocity gradient tensor 
  complex u11(nyp,nz,nxh),u12(nyp,nz,nxh),u13(nyp,nz,nxh)
  complex u21(nyp,nz,nxh),u22(nyp,nz,nxh),u23(nyp,nz,nxh) 
  complex u31(nyp,nz,nxh),u32(nyp,nz,nxh),u33(nyp,nz,nxh)
  complex,dimension(nyp,nz,nxh) :: Lu,Lv,Lw ! Laplacian terms - Ryan 7/24/23

  common/data1/ pi,dt,theta,wavz(nz),wavx(nxh),c(nyp),yc(nyp)
  common/data2/ re,alpha,beta,xl,zl,yl,time,dyde,xstart,uinf
  common/energ/ e1(nz,nxh),ckx(nz,nxh),tm1
  common/grnfct/nt,nb,p1t,p1b,p2t,p2b
  common/iocntrl/irstrt,nsteps,iprnfrq,print3d
  common/pre2forc/ spread(mz2,mx2),fxintg(nyp,mz,mx)
  common/pre4f/ fzintg(nyp,mz,mx),fyintg(nyp,mz,mx)
  common/pre5f/ it
  common/params/ gain, ugain

  integer bfhead   ! index # for the buffer zone head face.
  integer bftail   ! index # for the buffer zone tail face.
  real bfgain(1200)  ! buffer zone stiffness distribution array.
  real bfugain(1200) ! buffer zone damping coef. distribution array.
  integer bfwidth
  common /buffer/ bfgain, bfugain, bfhead, bftail, bfwidth

  real,dimension(nyp,nz,nx) :: initu, initv, initw
  common /init/ initu, initv, initw

  integer geomtype, flow_select
  common /geometry/ geomtype, flow_select

  im = (0.0,1.0)


!   upon exit from this routine the following initial conditions
!   must have been established.

!   phi,fn,fnm1: fourier(x) and chebyshev(y) transformed (complex)
!    see sub. phirhs

!   u0,h1,h1nm1: fourier and chebyshev tranformed (real)
!    see sub. uzero 

!   if this run is initialized from a restart package, read in the
!   old data & exit

  
  if(irstrt.ne.1) then


!    set up initial cond in real space
!    and then transform

!    subroutine readuv is supplied by the user. please give me back
!    u, v, and w stored in omx, omy, omz
!    respectively. 
     
     call readuv(omx,omy,omz)

     call scram(omy,v)
     call xyzfft(v,w1,w2,wfft1,wfft2,wfft3,wfft4,-1)

     call scram(omx,u)
     call xyzfft(u,w1,w2,wfft1,wfft2,wfft3,wfft4,-1)

     call scram(omz,w)
     call xyzfft(w,w1,w2,wfft1,wfft2,wfft3,wfft4,-1)

!   calc omy

     do k=1,nxh
        do j=1,nz
           do i=1,nyp
              omy(i,j,k) = im*wavz(j)*u(i,j,k) - im*wavx(k)*w(i,j,k)
           end do
        end do
     end do

!    because of periodicity, the 0,0 mode for omy is always 
!    zero.
     
     do i=1,nyp
        omy(i,1,1) = (0.0,0.0)
     end do

     call norm(u)
     call norm(v)
     call norm(w)
     call norm(omy)

!    initialize uo: 1st fou. mode of streamwise vel.
     
     do i=1,nyp
        u0(i) = real(u(i,1,1))
        w0(i) = real(w(i,1,1))
     end do
    
!    calc remaining components of velocity
     
    call veloc(u,v,w,u0,w0,omy,wrkc)

    call norm(u)
    call norm(v)
    call norm(w)
    call norm(omy)

!    calculate the vorticity and the nonlinear terms.
  
    call vort(u,v,w,omx,omz,wrkc)

    call norm(omx)
    call norm(omz)

! Calculate derivates needed for conformation tensor sources

    call derivscji(u,v,w,u11,u12,u13,u21,u22,u23,u31,u32,u33,Lu,Lv,Lw)

    call norm(u11) ! asdf
    call norm(u12)
    call norm(u13)
    call norm(u21)
    call norm(u22)
    call norm(u23)
    call norm(u31)
    call norm(u32)
    call norm(u33)

    ! Laplacian terms - Ryan 7/24/23
    call norm(Lu)
    call norm(Lv)
    call norm(Lw)

!    transform all data into y-physical
  
    is = 1

    call yfft(u,wfft1,wfft2,wfft3,wfft4,is)
    call yfft(v,wfft1,wfft2,wfft3,wfft4,is)
    call yfft(w,wfft1,wfft2,wfft3,wfft4,is)
    call yfft(omx,wfft1,wfft2,wfft3,wfft4,is)
    call yfft(omy,wfft1,wfft2,wfft3,wfft4,is)
    call yfft(omz,wfft1,wfft2,wfft3,wfft4,is)

    call yfft(u11,wfft1,wfft2,wfft3,wfft4,is)    
    call yfft(u12,wfft1,wfft2,wfft3,wfft4,is)
    call yfft(u13,wfft1,wfft2,wfft3,wfft4,is)
    call yfft(u21,wfft1,wfft2,wfft3,wfft4,is)
    call yfft(u22,wfft1,wfft2,wfft3,wfft4,is)
    call yfft(u23,wfft1,wfft2,wfft3,wfft4,is)
    call yfft(u31,wfft1,wfft2,wfft3,wfft4,is)
    call yfft(u32,wfft1,wfft2,wfft3,wfft4,is)
    call yfft(u33,wfft1,wfft2,wfft3,wfft4,is)

    ! Laplacian terms - Ryan 7/24/23
    call yfft(Lu,wfft1,wfft2,wfft3,wfft4,is)
    call yfft(Lv,wfft1,wfft2,wfft3,wfft4,is)
    call yfft(Lw,wfft1,wfft2,wfft3,wfft4,is)

! New vcw3d (includes vcw3dp stuff) - Ryan 6/28/22
      call vcw3d(u,v,w,omx,omy,omz,fn,gn,               &
                 u11,u12,u13,u21,u22,u23,u31,u32,u33,Lu,Lv,Lw)

!    transform all data into y-spectral
  
     is = -1

     call yfft(v,wfft1,wfft2,wfft3,wfft4,is)
     call yfft(omy,wfft1,wfft2,wfft3,wfft4,is)
     call yfft(fn,wfft1,wfft2,wfft3,wfft4,is)
     call yfft(gn,wfft1,wfft2,wfft3,wfft4,is)
     call yfft(omz,wfft1,wfft2,wfft3,wfft4,is)

     call norm(v)
     call norm(omy)
     call norm(fn)
     call norm(gn)
     call norm(omz)

     do i=1,nyp
        h1n(i)=real(gn(i,1,1))
        h3n(i)=real(omz(i,1,1))
     end do

     call nonlin(fn,omz,gn,wrkc,wrk1)

     call norm(fn)
     call norm(gn)

!    for the first time step fnm1=fn & h1nm1=h1n

     do k=1,nxh
        do j=1,nz
           do i=1,nyp
              fnm1(i,j,k)=fn(i,j,k)
              gnm1(i,j,k)=gn(i,j,k)
           end do
        end do
     end do
     do i=1,nyp
        h1nm1(i)=h1n(i)
        h3nm1(i)=h3n(i)
     end do

  else

!    read restart file
     open(112, file = 'setup/restart', status = 'old', form = 'unformatted')   ! alex 09/14/2018
     read(112) time
     read(112) v,fn,fnm1,omy,gn,gnm1,u0,h1n,h1nm1,w0,h3n,h3nm1
     read(112) fxintg, fyintg, fzintg
     close(112)                                           ! alex 09/14/2018

  end if

end subroutine initial

!---------------------------------------------------------------
!---------------------------------------------------------------

subroutine phirhs(v,wrkc,wrk1,fn,fnm1,bcbot,bctop)
  use omp_lib
  use grid_size
!    if theta .le. 0.99 crank-nicholson,
!    works on both odd and even parts simultaneously
  
  complex v(nyp,nz,nxh),wrkc(nyp,nz,nxh),wrk1(nyp,nz,nxh)
  complex fn(nyp,nz,nxh),fnm1(nyp,nz,nxh)
  complex bcbot(nz,nxh),bctop(nz,nxh)
  common/data1/ pi,dt,theta,wavz(nz),wavx(nxh),c(nyp),yc(nyp)
  common/data2/ re,alpha,beta,xl,zl,yl,time,dyde,xstart,uinf

!**********************************************************************
!                                                                     *
!     THIS SUBROUTINE EVALUATES THE RHS OF THE POISSON EQN. FOR THE   *
!     PARTICULAR PHI.                                                 *
!                                                                     *
!     SEE P. 139 JFM V.177 (1987) KMM                                 *
!                                                                     *
!     RHS=-(1.-THTA)(PHI"(N)-(KX*KX+KZ*KZ)*PHI(N))/THTA               *
!         -RE*PHI(N)/(DT*THTA)                                        *
!         -(3FN-FNM1)*RE/(2*THTA)                                     *
!                                                                     *
!**********************************************************************

!    first reset bc's to 0 since they may have been altered in 
!    pressure sol'n.
 
  do k=1,nxh
     do j=1,nz
        bctop(j,k)=(0.0,0.0)
        bcbot(j,k)=(0.0,0.0)  
     end do
  end do
!$omp parallel do
  do k=1,nxh
     do j=1,nz
        do i=1,nyp
           wrkc(i,j,k)= -re/(dt*theta)*v(i,j,k)-re/(2.*theta)*(3.*fn(i,j,k) -fnm1(i,j,k))
        end do
     end do
  end do

  if(theta.gt..99) then
!$omp parallel do
     do k=1,nxh
        do j=1,nz
           do i=1,nyp
              v(i,j,k) = wrkc(i,j,k)
           end do
        end do
     end do
  end if
!$omp parallel do
  do k=1,nxh
     do j=1,nz
        w2 = wavz(j)**2 +wavx(k)**2
        do i=1,nyp
           wrkc(i,j,k) = wrkc(i,j,k) +(1.-theta)*w2*v(i,j,k)/theta
        end do
     end do
  end do

!    first deriv.
  call cderiv(v,wrk1)

!    second deriv.
  call cderiv(wrk1,v)
!$omp parallel do
  do k=1,nxh
     do j=1,nz
        do i=1,nyp
           v(i,j,k)=wrkc(i,j,k) -(1.-theta)*v(i,j,k)/theta
        end do
     end do
  end do

end subroutine phirhs

!---------------------------------------------------------------

subroutine omyrhs(omy,wrkc,wrk1,gn,gnm1,bcbot,bctop)
use omp_lib
use grid_size
!    if theta .le. 0.99 crank-nickolson,
!    works on both odd and even parts simultaneously
  
  complex omy(nyp,nz,nxh),wrkc(nyp,nz,nxh),wrk1(nyp,nz,nxh)
  complex gn(nyp,nz,nxh),gnm1(nyp,nz,nxh)
  complex bcbot(nz,nxh),bctop(nz,nxh)
  common/data1/ pi,dt,theta,wavz(nz),wavx(nxh),c(nyp),yc(nyp)
  common/data2/ re,alpha,beta,xl,zl,yl,time,dyde,xstart,uinf

!**********************************************************************
!                                                                     *
!     THIS SUBROUTINE EVALUATES THE RHS OF THE POISSON EQN. FOR THE   *
!     PARTICULAR OMY.                                                 *
!                                                                     *
!     SEE P. 139 JFM V.177 (1987) KMM                                 *
!                                                                     *
!     RHS=-(1.-THTA)(OMY"(N)-(KX*KX+KZ*KZ)*OMY(N))/THTA               *
!         -RE*OMY(N)/(DT*THTA)                                        *
!         -(3FN-FNM1)*RE/(2*THTA)                                     *
!                                                                     *
!**********************************************************************
!
!    FIRST RESET BC'S TO 0 SINCE THEY MAY HAVE BEEN ALTERED IN 
!    PRESSURE SOL'N.

  do k=1,nxh
     do j=1,nz
        bctop(j,k)=(0.0,0.0)
        bcbot(j,k)=(0.0,0.0)  
     end do
  end do

!$omp parallel do
  do k=1,nxh
     do j=1,nz
        do i=1,nyp
           wrkc(i,j,k)= -re/(dt*theta)*omy(i,j,k)-re/(2.*theta)*(3.*gn(i,j,k) -gnm1(i,j,k))
        end do
     end do
  end do
  
  if(theta.gt..99) then
!$omp parallel do
     do k=1,nxh
        do j=1,nz
           do i=1,nyp
              omy(i,j,k) = wrkc(i,j,k)
           end do
        end do
     end do
  end if

!$omp parallel do
  do k=1,nxh
     do j=1,nz
        w2 = wavz(j)**2 +wavx(k)**2
        do i=1,nyp
           wrkc(i,j,k) = wrkc(i,j,k) +(1.-theta)*w2*omy(i,j,k)/theta
        end do
     end do
  end do

!    first deriv.
  call cderiv(omy,wrk1)

!    second deriv.
  call cderiv(wrk1,omy)
!$omp parallel do
  do k=1,nxh
     do j=1,nz
        do i=1,nyp
           omy(i,j,k)=wrkc(i,j,k) -(1.-theta)*omy(i,j,k)/theta
        end do
     end do
  end do
  
end subroutine omyrhs

!----------------------------------------------------------

subroutine uzero(u0,h1,h1l,a,t,wrk1)
  use omp_lib
  use grid_size
  
  real u0(nyp), f(nyp), h1(nyp),h1l(nyp)
  real a(nyhp,nz,3),t(nyp),wrk1(nyp,nz)
    ! New BC variables - Ryan 6/8/23
    real :: atu,btu,gtu,abu,bbu,gbu
    real :: atw,btw,gtw,abw,bbw,gbw
  common/data1/ pi,dt,theta,wavz(nz),wavx(nxh),c(nyp),yc(nyp)
    common/u0bcs/ atu,btu,gtu,abu,bbu,gbu
    common/w0bcs/ atw,btw,gtw,abw,bbw,gbw ! adding separate BCs for w
  common/data2/ re,alpha,beta,xl,zl,yl,time,dyde,xstart,uinf

!**********************************************************************
!                                                                     *
!    THIS SUBROUTINE SETS UP TO SOLVE THE X-MOMENTUM EQN. FOR KX=KZ=0,*
!    THE ZEROTH FOURIER MODE OF STREAMWISE VELOCITY.                  *
!                                                                     *
!      dUo/dt = H1 + F + (1/Re)*d(dUo/dy)/dy                          *
!                                                                     *
!    DISCRETE FORM:                                                   *
!                                                                     *
!    Uo"(n+1)-(Re/thta*Dt)Uo(n+1)=-(Re/thta*Dt)Uo(n)-Re*F/thta        *
!             -(Re/2*thta)(3H1(n)-H1(n-1))-(1-thta)Uo"(n)/thta        *
!                                                                     *
!                                                                     *
!                                                                     *
!**********************************************************************

  do i = 1,nyp
     f(i) = 0.0
  end do
!$omp parallel do    
  do i=1,nyp
     wrk1(i,1)=-re*u0(i)/(dt*theta)-re*(3.*h1(i)-h1l(i))/(2.*theta)-re*f(i)/theta
  end do

  if(theta .gt. 0.99) then
     do i=1,nyp
        u0(i) = wrk1(i,1)
     end do
  else

!    calc. 2nd deriv of u0
     call cdrv1d(u0,t)  !first deriv. 
     call cdrv1d(t,u0)  !second deriv.

     do i=1,nyp
        u0(i)=wrk1(i,1) -(1.-theta)*u0(i)/theta
     end do
  end if

!      call 1d penta-diag. solver
!      rhs on input is u0, sol'n is returned in u0.
!      solver assumes mixed bc's are of form:

!      at*u0 + bt*du0/dy =gt  at y=1
!      ab*u0 + bb*du0/dy =gb  at y=-1

!      bc's have been set in subroutine setstuf or by other means
!      set time-step term (g)
  
  g=re/(theta*dt)

  call psolv1d(u0,g,dyde,atu,btu,gtu,abu,bbu,gbu,c,a,t,wrk1)

end subroutine uzero

! *************************************************

subroutine wzero(u0,h1,h1l,a,t,wrk1)
  use omp_lib
  use grid_size

  real u0(nyp), f(nyp), h1(nyp),h1l(nyp)
  real a(nyhp,nz,3),t(nyp),wrk1(nyp,nz)
    ! New BC variables - Ryan 6/8/23
    real :: atu,btu,gtu,abu,bbu,gbu
    real :: atw,btw,gtw,abw,bbw,gbw
  common/data1/ pi,dt,theta,wavz(nz),wavx(nxh),c(nyp),yc(nyp)
    common/u0bcs/ atu,btu,gtu,abu,bbu,gbu
    common/w0bcs/ atw,btw,gtw,abw,bbw,gbw ! adding separate BCs for w
  common/data2/ re,alpha,beta,xl,zl,yl,time,dyde,xstart,uinf

!**********************************************************************
!                                                                     *
!    THIS SUBROUTINE SETS UP TO SOLVE THE X-MOMENTUM EQN. FOR KX=KZ=0,*
!    THE ZEROTH FOURIER MODE OF STREAMWISE VELOCITY.                  *
!                                                                     *
!      dUo/dt = H1 + F + (1/Re)*d(dUo/dy)/dy                          *
!                                                                     *
!    DISCRETE FORM:                                                   *
!                                                                     *
!    Uo"(n+1)-(Re/thta*Dt)Uo(n+1)=-(Re/thta*Dt)Uo(n)-Re*F/thta        *
!             -(Re/2*thta)(3H1(n)-H1(n-1))-(1-thta)Uo"(n)/thta        *
!                                                                     *
!                                                                     *
!                                                                     *
!**********************************************************************
!    SET UP RHS: NOTE THE FORCING TERM F WILL ONLY EFFECT 
!       THE ZEROTH CHEB MODE


  do i = 1,nyp
     f(i) = 0.0
  end do

!$omp parallel do    
  do i=1,nyp
     wrk1(i,1)=-re*u0(i)/(dt*theta)-re*(3.*h1(i)-h1l(i))/(2.*theta)-re*f(i)/theta
  end do

  if(theta.gt..99) then
     do i=1,nyp
        u0(i) = wrk1(i,1)
     end do
  else
!  calc. 2nd deriv of u0
!  first deriv.
     call cdrv1d(u0,t)
!  second deriv.
     call cdrv1d(t,u0)
     
     do i=1,nyp
        u0(i)=wrk1(i,1) -(1.-theta)*u0(i)/theta
     end do
  end if


!    CALL 1D PENTA-DIAG. SOLVER
!    RHS ON INPUT IS U0, SOL'N IS RETURNED IN U0.
!    SOLVER ASSUMES MIXED BC'S ARE OF FORM:
!    AT*U0 + BT*dU0/dy =GT  at y=1
!    AB*U0 + BB*dU0/dy =GB  at y=-1
!    BC'S HAVE BEE SET IN SBR SETSTUF OR BY OTHER MEANS
!    SET TIME-STEP TERM (G)

  g=re/(theta*dt)

  !call psolv1d(u0,g,dyde,at,bt,gt,ab,bb,gb,c,a,t,wrk1)
  call psolv1d(u0,g,dyde,atw,btw,gtw,abw,bbw,gbw,c,a,t,wrk1) ! Changing for Couette flow


end subroutine wzero

!----------------------------------------------------------------------

subroutine veloc(u,v,w,u0,w0,omy,wrkc)
  use omp_lib
  use grid_size

  complex u(nyp,nz,nxh),v(nyp,nz,nxh),w(nyp,nz,nxh)
  complex omy(nyp,nz,nxh)
  complex wrkc(nyp,nz,nxh)
  real u0(nyp), w0(nyp)
  complex im
  common/data1/ pi,dt,theta,wavz(nz),wavx(nxh),c(nyp),yc(nyp)
  common/data2/ re,alpha,beta,xl,zl,yl,time,dyde,xstart,uinf

!**********************************************************************
!
!     CALCULATE THE REMAINING COMPONENTS OF VELOCITY FROM CONTINUITY 
!     AND NORMAL VORTICITY IN SPECTRAL SPACE
!    
!     U(N,Kz,Kx) = 
!          -((i*Kx)(-dV/dy)+(i*Kz)(Omy))/(Kx**2+Kz**2) : Kx**2+Kz**2 > 0 
!     for Kx**2 + Kz**2 =0 HAVE U0 ABOVE AS REAL PART OF U(I,KX=0).
!
!     W(N<Kz,Kx) =
!          -((i*Kz)(-dV/dy)-(i*Kx)(Omy))/(Kx**2+Kz**2) : Kx**2+Kz**2 > 0 
!     For Kx**2 + Kz**2 = 0, W(n,Kz,Kx) = 0.0
!
!
!**********************************************************************
  
  im = (0.0,1.0)
  call cderiv(v,wrkc)

!********************* calc. for all but first mode (k,j=1)
!$omp parallel do    
  do k=2,nxh
     do j=1,nz
        w2 = wavx(k)**2 + wavz(j)**2
        do i=1,nyp
           u(i,j,k)=-(-im*wavx(k)*wrkc(i,j,k)+im*wavz(j)*omy(i,j,k))/w2
           w(i,j,k)=-(-im*wavz(j)*wrkc(i,j,k)-im*wavx(k)*omy(i,j,k))/w2
        end do
     end do
  end do

  k = 1          ! skip only j=k=1 mode
!$omp parallel do
  do j=2,nz
     w2 =  wavz(j)**2
     do i=1,nyp
        u(i,j,k)=-(im*wavz(j)*omy(i,j,k))/w2
        w(i,j,k)=-(-im*wavz(j)*wrkc(i,j,k))/w2
     end do
  end do

!********************* calc. for first mode

!  update the first (kz=0) mode of streamwise velocity.  the imag. part
!  of the first mode has been used to store the real part of the last
!  fourier mode (nz/2 + 1). 
  
  do i=1,nyp
     u(i,1,1) = cmplx(u0(i),0.0)
     v(i,1,1) = 0.0
     w(i,1,1) = cmplx(w0(i),0.0)
  end do

end subroutine veloc

!----------------------------------------------------------

subroutine norm(a)
use omp_lib
use grid_size
  complex a(nyp,nz,nxh)

!**********************************************************************

!  'normalize' fields, ensure symmetries and reality

!!  k=j=1 mode
  do i = 1,2
  end do

!$omp parallel do
  do i=1,nyp
     a(i,1,1)= real(a(i,1,1))
  end do
!$omp end parallel do
!  k=1, j>1 modes ensure conjg sym.

!$omp parallel do
  do j=2,nz/2
     jp = nz +2 -j
     do i=1,nyp
        a(i,j,1) = 0.5*( a(i,j,1) +conjg(a(i,jp,1)) )
     end do
  end do
!$omp end parallel do

!$omp parallel do
  do j=2,nz/2
     jp = nz +2 -j
     do i=1,nyp
        a(i,jp,1) = conjg(a(i,j,1))
     end do
  end do
!$omp end parallel do

!  zero out highest mode in y dir, for reality
  !nzhp = nz/2+1

!$omp parallel do
  do k = 1,nxh
     do i = 1,nyp
        a(i,nzhp,k) = 0.0
     end do
  end do
!$omp end parallel do

end subroutine norm

!----------------------------------------------------------------------

subroutine vort(u,v,w,omx,omz,wrkc)
  use grid_size
  complex u(nyp,nz,nxh), v(nyp,nz,nxh), w(nyp,nz,nxh)
  complex omx(nyp,nz,nxh), omz(nyp,nz,nxh), wrkc(nyp,nz,nxh)
  complex im
  common/data1/ pi,dt,theta,wavz(nz),wavx(nxh),c(nyp),yc(nyp)

!**********************************************************************
!     omx =  dw/dy-dv/dz                                              *
!     omy =  du/dz-dw/dx   calculated above                           *
!     omz =  dv/dx-du/dy                                              *
!     all quantities are spectral                                     *
!**********************************************************************

  im = (0.0,1.0)

  call cderiv (u,omz)

!$omp parallel do
  do k=1,nxh
     do j=1,nz
        do i=1,nyp
           omz(i,j,k) = im*wavx(k)*v(i,j,k)  - omz(i,j,k)
        end do
     end do
  end do


  call cderiv (w,omx)

!$omp parallel do 
  do k=1,nxh
     do j=1,nz
        do i=1,nyp
           omx(i,j,k) = omx(i,j,k) - im*wavz(j)*v(i,j,k)
        end do
     end do
  end do
 

end subroutine vort

!----------------------------------------------------------------------
!----------------------------------------------------------------------

subroutine nonlin(fn,omz,gn,wrkc,wrk1)
  use grid_size
  complex fn(nyp,nz,nxh), gn(nyp,nz,nxh), omz(nyp,nz,nxh)
  complex wrkc(nyp,nz,nxh), wrk1(nyp,nz,nxh)
  complex im
  common/data1/ pi,dt,theta,wavz(nz),wavx(nxh),c(nyp),yc(nyp)
  common/data2/ re,alpha,beta,xl,zl,yl,time,dyde,xstart,uinf

!**********************************************************************
!     THIS SUBROUTINE CALCULATES THE NONLINEAR TERM IN THE 4th ORDER  *
!     SYSTEM.                                                         *
!                                                                     *
!                                                                     *
!      FN=-d(dH1/dx + dH3/dz)/dy+d(dH2/dx)/dx +d(dH2/dz)/dz           *
!                                                                     *
!      GN=  dH1/dz - dH3/dx                                           *
!**********************************************************************     

  im = (0.0,1.0)

  do k=1,nxh
     do j=1,nz
        do i=1,nyp
           wrkc(i,j,k) = im*( wavz(j)*(omz(i,j,k)) + wavx(k)*(gn(i,j,k)) )
        end do
     end do
  end do

  call cderiv(wrkc,wrk1)

  do k=1,nxh
     do j=1,nz
        w2 = -wavx(k)**2 -wavz(j)**2
        do i=1,nyp
           fn(i,j,k) = w2*(fn(i,j,k) ) - wrk1(i,j,k)           ! note fnorm addition 
        end do
     end do
  end do

!   calculate gn
  do k=1,nxh
     do j=1,nz
        do i=1,nyp
           gn(i,j,k) = im*( wavz(j)*( gn(i,j,k) ) - wavx(k)*(omz(i,j,k)) )
        end do
     end do
  end do

end subroutine nonlin


subroutine gfcn(gfns,a,t,wrkc,wrk1,bct,bcb,dnvb,dnvt,nb,nt,pb,pt)

!**********************************************************************
!*                                                                    *
!*    This subroutine performs the homogeneous solution of the        *
!*                                                                    *
!*    system:                                                         *
!*                                                                    *
!*     d  d(phi)   d  d(phi)   d  d(phi)      Re                      *
!*     -- ------ + -- ------ + -- ------ - -------- phi = 0           *
!*     dx   dx     dy   dy     dz   dz     theta*dt                   *
!*                                                                    *
!*     d  d(Vm)   d  d(Vm)   d  d(Vm)                                 *
!*     -- ----- + -- ----- + -- ----- =  phi                          *
!*     dx  dx     dy  dy     dz  dz                                   * 
!*                                                                    *
!*    with bc's                                                       *
!*                                                                    *
!*        phi(1) = pt, phi(-1) = pb, vm(1) = vm(-1) = 0               *
!*                                                                    *
!**********************************************************************
 
  use grid_size

  complex wrkc(nyp,nz,NXH), wrk1(nyp,nz,NXH)
  complex gfns(nyp,nz,NXH)
  complex bct(nz,NXH), bcb(nz,NXH)
  complex dnvb(nz,NXH), dnvt(nz,NXH)

  real a(nyhp,nz,3), t(nyp)

  common/data1/ pi, dt, theta, wavz(nz), wavx(NXH), c(nyp), yc(nyp)
  common/data2/ Re, alpha, beta, xl, zl, yl, time, dyde, xstart, uinf

!... Set RHS & bc's for phi eqn.

  g = Re / (theta * dt)

  do k = 1, NXH
     do j = 1, nz
        bct(j,k) = cmplx(pt, pt)
        bcb(j,k) = cmplx(pb, pb)
     end do
  end do

  do j = nz/2+1, nz
     bct(j,1) = conjg(bct(j,1))
     bcb(j,1) = conjg(bcb(j,1))
  end do

  do k=1,NXH
     do j=1,nz
        do i=1,nyp
           gfns(i,j,k) = (0.0,0.0)
        end do
     end do
  end do

  bct(1,1) = cmplx(pt, 0.0)
  bcb(1,1) = cmplx(pb, 0.0)

!    Dirichlet bc's -- ib=0

  ib = 0

!... Solver for phi by the  Poisson solver, "solve,"
!... sol'n returned in gfns(i,j).

  do k = 1, NXH
     x = wavx(k)**2
     call solve (gfns(1,1,k), x, g, dyde, bct(1,k), bcb(1,k), ib, wavz, c, a, t, wrkc)
  end do
      
!----------------------------------------------------------------------
!... Set RHS & bc's for Vm eqn,
!...     RHS is gfns returned from above call

!... NOTE:  no dt term in Vm eqn.
     
  g = 0.0

  do k = 1,NXH
     do j = 1,nz
        bct(j,k) = (0.0, 0.0)
        bcb(j,k) = (0.0, 0.0)
     end do
  end do

  do k = 1, NXH
     x = wavx(k)**2
     call solve (gfns(1,1,k), x, g, dyde, bct(1,k), bcb(1,k), ib, wavz, c, a, t, wrkc)
  end do

  do i = 1, nyp
     gfns(i, 1, 1) = real(gfns(i, 1, 1))
  end do

!  k=1, j>1 modes ensure conjg sym.
  
  do j = 2, nz/2
     jp = nz + 2 - j
     do i = 1, nyp
        gfns(i,j,1) = 0.5 * (gfns(i,j,1) + conjg(gfns(i,jp,1)))
     end do
  end do

  do j = 2, nz/2
     jp = nz + 2 - j
     do i = 1, nyp
        gfns(i,jp,1) = conjg(gfns(i,j,1))
     end do
  end do

!... Calculate the required wall derivs. of gfns & store in dnvt,dnvb
  
  iy1 = 1
  iy2 = 2

  if (nb .eq. 1) then
     call c1derbw(gfns, wrkc, iy1)
  else if (nb .eq. 2) then
     call cderiv(gfns, wrk1)
     call c1derbw(wrk1, wrkc, iy1) 
  end if

  if (nt .eq. 1) then
     call c1dertw(gfns, wrkc, iy2)        
  else if (nt .eq. 2) then
     call cderiv(gfns, wrk1)
     call c1dertw(wrk1, wrkc, iy2)        
  end if

!.... Store the wall derivatives in "dnvt(j,k)" and "dnvb(j,k)".
     
  do k = 1, NXH
     do j = 1, nz
        dnvt(j, k) = wrkc(iy2, j, k)
        dnvb(j, k) = wrkc(iy1, j, k)
     end do
  end do

end subroutine gfcn

!---------------------------------------------------------------

subroutine readuv(u1,u2,u3)
  use grid_size
  implicit none

  integer, parameter :: nzlowres = 4, mzlowres = 3*nzlowres/2

  real,dimension(nyp,nz,nx) :: u1,u2,u3,initu,initv,initw,sc! added sc - Ryan 3-2-22

  real re,alpha,beta,xl,zl,yl,time,dyde,xstart,uinf
  integer i, j, k, irstrt, nsteps, iprnfrq, numjstrips, jstrip, jcount

  common /init/ initu, initv, initw
  common/iocntrl/ irstrt, nsteps, iprnfrq
  common/data2/ re,alpha,beta,xl,zl,yl,time,dyde,xstart,uinf

  real u_lowres(nyp,nzlowres,nx), v_lowres(nyp,nzlowres,nx), w_lowres(nyp,nzlowres,nx)
  real pi


  pi = 2.0 * acos(0.0)

  if (irstrt .eq. 0) then
    do  k = 1,nyp
       do  j = 1,nz       
          do  i = 1,nx

            ! Note (k,j,i) order!
             u1(k,j,i) = initu(k,j,i)
             u2(k,j,i) = initv(k,j,i)
             u3(k,j,i) = initw(k,j,i)

          end do
       end do
    end do
  end if

  numjstrips = nz/nzlowres

  if (irstrt .eq. 2) then
     write(*,*) "In readuv nzlowres = ", nzlowres, "and numjstrips = ", numjstrips
     open(80, file = 'setup/low-resolution-velocities', status = 'old', form = 'unformatted')   ! alex 09/14/2018
     read(80) u_lowres
     read(80) v_lowres
     read(80) w_lowres
     close(80)
     do  k = 1,nyp
        do jstrip = 1,numjstrips
           do jcount = 1,nzlowres
              j = ((jstrip-1)*nzlowres) + jcount    
              do  i = 1,NX
                 u1(k,j,i) = u_lowres(k,jcount,i)
                 u2(k,j,i) = v_lowres(k,jcount,i)
                 u3(k,j,i) = w_lowres(k,jcount,i)
              end do
           end do
        end do
     end do
  end if

end subroutine readuv

subroutine forcepre
!***********************************************************************
!*            Calculate preliminary data for forcing  terms            *
!***********************************************************************
  use grid_size
  implicit None

  integer, parameter :: nzlowres = 4, mzlowres = 3*nzlowres/2

  real spread(mz2,mx2), fxintg(nyp,mz,mx), fzintg(nyp,mz,mx), fyintg(nyp,mz,mx), seght(nyp)
  real re, alpha, beta, xl, yl, zl, time, dyde, xstart, uinf, gain, ugain
  real ywall, ytop, y, ytemp, u_temp, ym1, yp1, vtemp, gwidth, dsq

  integer imatrix(nyp,mz,mx)
  integer it, irstrt, iprnfrq, nsteps, i, j, k, numjstrips, jstrip, jcount


  common/pre2forc/ spread,fxintg
  common/pre4f/ fzintg,fyintg
  common/pre5f/ it
  common/geommat/ imatrix

  common/data2/ re,alpha,beta,xl,zl,yl,time,dyde,xstart,uinf
  common/params/ gain, ugain
  common/preforce/ seght
  common/iocntrl/irstrt,nsteps,iprnfrq

  integer bfhead            ! Index # for the Buffer zone head face.
  integer bftail            ! Index # for the Buffer zone tail face.
  real bfgain(1200)          ! Buffer zone stiffness distribution array.
  real bfugain(1200)         ! Buffer zone damping coef. distribution array.
  integer bfwidth, bfhw, bfcntr
  common /buffer/ bfgain, bfugain, bfhead, bftail, bfwidth

  real,dimension(nyp,nz,nx) :: initu, initv, initw
  real ycrd(nyp)
  common/init/initu, initv, initw, ycrd

  integer slopelength, slopepos, ib
  integer ribvpos(mzp), ribtempvpos(mzp)
  real c1, c2, c3, c4, c5, c6, ydelta, bconst, eta, delxm, x, pi

  integer readvdes
  real vdes(mx), xsuction
  common /suction/ vdes, xsuction, readvdes

  integer kwall, kmaxsurf
  common /kwallpos/ kwall, kmaxsurf

  integer geomtype, flow_select
  common /geometry/ geomtype, flow_select

  real geomwidth, geomheight
  integer numgeoms, ribhpos, geomspacing
  common /geomprops/ geomwidth, geomheight
  
  real fxintLR(nyp,mzlowres,mx), fyintLR(nyp,mzlowres,mx), fzintLR(nyp,mzlowres,mx)
  real dPdx, R_tau
  common/pressure/ dPdx, R_tau


  
  delxm = xl /float(mx-1)
  
  pi = 2.0 * acos(0.0)

  bfhead  = 2
  bftail  = bftail_              ! alex 09/20/2018 - Check if correct
! bftail  = (30*nx/128)+1        ! alex 09/20/2018 - Check if correct
  bfwidth = (bftail - bfhead)
  bfhw    = bfwidth / 2
  bfcntr  = bfhw + 1
  bfgain  = 0.0
  bfugain = 0.0
  do i = 1, (bfwidth + 1)
     ib = i - bfcntr
     if (iabs(ib) .le. bfhw/3) then
        bfgain(i) = 20.0 * exp(-60.*(float(ib)/float(bfhw))**2)
     end if
     bfugain(i)  =  20.0 * exp(-6.0*(float(ib)/float(bfhw))**2)
  end do

  slopelength = bfwidth/2.0
  ywall = ycoord(kwall)
  ytop = yl - ywall

  call init_flow(initu,initv,initw)

  write(*,*) "Geomtype equals ", geomtype
 
  kmaxsurf = 0
  do k = 1,nyp
     do j = 1,mz
        do i = 1,mx
           if ((imatrix(k,j,i) .eq. 1) .or. (imatrix(k,j,i) .eq. 3) .or. (imatrix(k,j,i) .eq. 4)) then
              if (k .ge. kmaxsurf) kmaxsurf = k
           end if
        end do
     end do
  end do
           
  write(*,*) "Kmaxsurf equals ", kmaxsurf

  xsuction = delxm * (bftail - 1)

  do i = 1,mx
     x = delxm * (i-1)
     if (x .le. xsuction) then
        vdes(i) = 0.0
     else
        if (readvdes .eq. 0) then
           vdes(i) = (-0.8604*Uinf) / sqrt(re*Uinf*(x+xstart-xsuction))
        else if (readvdes .eq. 1) then
           read(699,*) vtemp
           vdes(i) = -vtemp
        end if
     end if
  end do
   
!!!.... Setup y-coordiante and free stream velocity.
!!  do i = 1, nyp
!!     ycrd(i) = (cos(float(i-1) * pi / float(ny)) + 1.) / 2.
!!     ycrd(i) = ycrd(i) * yl
!!  end do

!.... Compute the segment height between neighboring y-mesh points.
  do i = 1, ny
     seght(i) = abs(ycoord(i+1) - ycoord(i))
  end do
  seght(nyp) = seght(ny)


!.... Calculate the force spreading matrix.
  gwidth = 2.0
  do i = 1, mx2
     do j = 1, mz2
        dsq = (float(i-1))**2 + (float(j-1))**2
        spread(j,i) = 0.
        if(dsq .lt. 13) spread(j,i) = exp(-gwidth * dsq)
     end do
  end do
     
!.... Initialize the integral forcing terms.
  if (irstrt .eq. 0) then
     do i = 1, mx
        do j = 1, mz
           do k = 1, nyp
              fxintg(k,j,i) = 0.0
              fyintg(k,j,i) = 0.0
              fzintg(k,j,i) = 0.0
           end do
        end do
     end do
  end if

  numjstrips = mz/mzlowres

 ! write(*,*) "In forcpre nzlowres = ", nzlowres, "and numjstrips = ", numjstrips

  if (irstrt .eq. 2) then
     open(90, file = 'setup/low-resolution-forces', status = 'old', form = 'unformatted')   ! alex 09/14/2018
     read(90) fxintLR
     read(90) fyintLR
     read(90) fzintLR
     close(90)
     do  k = 1,nyp
        do jstrip = 1,numjstrips
           do jcount = 1,mzlowres
              j = ((jstrip-1)*mzlowres) + jcount    
              do  i = 1,mx
                 fxintg(k,j,i) = fxintLR(k,jcount,i)
                 fyintg(k,j,i) = fyintLR(k,jcount,i)
                 fzintg(k,j,i) = fzintLR(k,jcount,i)
              end do
           end do
        end do
     end do
  end if
  

end subroutine forcepre


!---------------------------------------------------------------


subroutine output(u,v,w,urms,vrms,wrms,umn,wrk1,wrkc,wfft1,wfft2,wfft3,wfft4,w1,w2)

  use grid_size
  implicit none 

  complex u(nyp,nz,nxh), v(nyp,nz,nxh), w(nyp,nz,nxh)
  complex wrkc(nyp,nz,nxh)
  complex im
  real wrk1(nyp,nz,nx), arms(nyp), amn(nyp)
  real wfft1(1),wfft2(1),w1(1), w2(1)
  real wfft3(1),wfft4(1)
  real urms(nyp), vrms(nyp), wrms(nyp), umn(nyp)
  real up(nyp,nz,nx), vp(nyp,nz,nx), wp(nyp,nz,nx)
  real wavz(nz), wavx(nxh), c(nyp), yc(nyp), spread(mz2,mx2)
  real fxintg(nyp,mz,mx), fzintg(nyp,mz,mx), fyintg(nyp,mz,mx)
  real pi, dt, theta, re, alpha, beta, xl, zl, yl, time, dyde, xstart, uinf
  real x(nx), y(nyp), z(nz)

  integer irstrt, nsteps, iprnfrq, i, j, k, it, filenumber

  character*16 :: filename

  common/data1/ pi,dt,theta,wavz,wavx,c,yc
  common/data2/ re,alpha,beta,xl,zl,yl,time,dyde,xstart,uinf
  common/iocntrl/irstrt,nsteps,iprnfrq
  common/pre2forc/ spread,fxintg
  common/pre4f/ fzintg,fyintg
  common/pre5f/ it

  open(60, file="GMU_out/fort.60")
  call ccopy (nz*nxh*nyp,u,1,wrkc,1)
  call xyzfft(wrkc,w1,w2,wfft1,wfft2,wfft3,wfft4,1)
  call unscram(wrkc,wrk1)
  write(60,*) (((wrk1(i,j,k),i=1,nyp),j=1,nz),k=1,nx)

  call ccopy (nz*nxh*nyp,v,1,wrkc,1)
  call xyzfft(wrkc,w1,w2,wfft1,wfft2,wfft3,wfft4,1)
  call unscram(wrkc,wrk1)
  write(60,*) (((wrk1(i,j,k),i=1,nyp),j=1,nz),k=1,nx)

  call ccopy (nz*nxh*nyp,w,1,wrkc,1)
  call xyzfft(wrkc,w1,w2,wfft1,wfft2,wfft3,wfft4,1)
  call unscram(wrkc,wrk1)
  write(60,*) (((wrk1(i,j,k),i=1,nyp),j=1,nz),k=1,nx)
 
end subroutine output

!---------------------------------------------------------------

subroutine cderiv(f,df)
  use omp_lib
  use grid_size
  complex f(nyp,nz,nxh),df(nyp,nz,nxh)
  common/data2/ re,alpha,beta,xl,zl,yl,time,dyde,xstart,uinf

  do i=1,2 !not useful loop to fool tau. without this, tau instrumented code doesn't work
  end do

!$omp parallel do
  do k=1,nxh
     do j=1,nz
        df(nyp,j,k)=0.
        df(ny,j,k)=2.*float(ny)*f(nyp,j,k)*dyde
     end do
  end do
!$omp end parallel do

!$omp parallel do
  do j = 1, nz
     do i=nym,1,-1
        do k=1,nxh
              df(i,j,k)=df(i+2,j,k)+2.*float(i)*f(i+1,j,k)*dyde
        end do
     end do
  end do
!$omp end parallel do


!$omp parallel do
  do k=1,nxh
     do j=1,nz
        df(1,j,k)=.5*df(1,j,k)
     end do
  end do
!$omp end parallel do




end subroutine cderiv

!-------------------------------------------------------------------

subroutine c0derbw(wrkc,wrk1,iy)
  use grid_size
  complex wrkc(nyp,nz,nxh),wrk1(nyp,nz,nxh)
  common/data2/ re,alpha,beta,xl,zl,yl,time,dyde,xstart,uinf

!***********************************************************************
!     this subroutine evaluates the function at the y=-1 bndry*
!     it assumes wrkc is spectral in both z and y (the cheb. direction)*
!     wrk1 is the normal deriv. at y=-1 & is physical in y.             *
!***********************************************************************

  do k=1,nxh
     do j=1,nz
        wrk1(iy,j,k)=0.
     end do
  end do
  sgn=-1.
  do i=1,nyp
     sgn=-sgn
     do k=1,nxh
        do j=1,nz
           wrk1(iy,j,k)=wrk1(iy,j,k)+sgn*wrkc(i,j,k)
        end do
     end do
  end do
  
end subroutine c0derbw

!-------------------------------------------------------------------

subroutine c1derbw(wrkc,wrk1,iy)
  use grid_size
  complex wrkc(nyp,nz,nxh),wrk1(nyp,nz,nxh)
  common/data2/ re,alpha,beta,xl,zl,yl,time,dyde,xstart,uinf

!***********************************************************************
!     this subroutine evaluates the normal derivative at the y=-1 bndry*
!     it assumes wrkc is spectral in both z and y (the cheb. direction)*
!     wrk1 is the normal deriv. at y=-1 & is physical in y.             *
!***********************************************************************

  do k=1,nxh
     do j=1,nz
        wrk1(iy,j,k)=0.
     end do
  end do
  sgn=1.
  do i=1,nyp
     sgn=-sgn
     rp=float(i-1)
     do k=1,nxh
        do j=1,nz
           wrk1(iy,j,k)=wrk1(iy,j,k)+rp*rp*sgn*wrkc(i,j,k)*dyde
        end do
     end do
  end do

end subroutine c1derbw

!----------------------------------------------------------

subroutine c2derbw(wrkc,wrk1,iy)
  use grid_size
  complex wrkc(nyp,nz,nxh),wrk1(nyp,nz,nxh)
  common/data2/ re,alpha,beta,xl,zl,yl,time,dyde,xstart,uinf

!***********************************************************************
!     this subroutine evaluates the second derivative at the y=-1 bndry*
!     it assumes wrkc is spectral in both z and y (the cheb. direction)*
!     wrk1 is the normal deriv. at y=-1 & is physical in y.            *
!***********************************************************************

  do k=1,nxh
     do j=1,nz
        wrk1(iy,j,k)=0.
     end do
  end do
  sgn=-1.
  do i=1,nyp
     sgn=-sgn
     rp=float(i-1)
     do k=1,nxh
        do j=1,nz
           wrk1(iy,j,k)=wrk1(iy,j,k)+sgn*rp*rp*(rp*rp-1.)*wrkc(i,j,k)*dyde*dyde/3.
        end do
     end do
  end do

end subroutine c2derbw

!-------------------------------------------------------------------

subroutine c0dertw(wrkc,wrk1,iy)
  use grid_size
  complex wrkc(nyp,nz,nxh),wrk1(nyp,nz,nxh)
  common/data2/ re,alpha,beta,xl,zl,yl,time,dyde,xstart,uinf

!***********************************************************************
!     this subroutine evaluates the function at the y=-1 bndry*
!     it assumes wrkc is spectral in both z and y (the cheb. direction)*
!     wrk1 is the normal deriv. at y=-1 & is physical in y.             *
!***********************************************************************

  do k=1,nxh
     do j=1,nz
        wrk1(iy,j,k)=0.0
     end do
  end do
  do i=1,nyp
     do k=1,nxh
        do j=1,nz
           wrk1(iy,j,k)=wrk1(iy,j,k)+wrkc(i,j,k)
        end do
     end do
  end do

end subroutine c0dertw

!----------------------------------------------------------------------

subroutine c1dertw(wrkc,wrk1,iy)
  use grid_size
  complex wrkc(nyp,nz,nxh),wrk1(nyp,nz,nxh)
  common/data2/ re,alpha,beta,xl,zl,yl,time,dyde,xstart,uinf

!**********************************************************************!
!     this subroutine evaluates the normal derivative at the y=1 bndry !
!     it assumes wrkc is spectral in both z and y (the cheb. direction)!
!     wrk is the normal deriv. at y=1 & is physical in y.              !
!**********************************************************************!

  do k=1,nxh
     do j=1,nz
        wrk1(iy,j,k)=0.
     end do
  end do
  do i=1,nyp
     rp=float(i-1)
     do k=1,nxh
        do j=1,nz
           wrk1(iy,j,k)=wrk1(iy,j,k)+rp*rp*wrkc(i,j,k)*dyde
        end do
     end do
  end do

end subroutine c1dertw

!----------------------------------------------------------------------

subroutine c2dertw(wrkc,wrk1,iy)
  use grid_size
  complex wrkc(nyp,nz,nxh),wrk1(nyp,nz,nxh)
  common/data2/ re,alpha,beta,xl,zl,yl,time,dyde,xstart,uinf

!***********************************************************************
!     this subroutine evaluates the normal derivative at the y=1 bndry *
!     it assumes wrkc is spectral in both z and y (the cheb. direction)*
!     wrk is the normal deriv. at y=1 & is physical in y.              *
!***********************************************************************

  do k=1,nxh
     do j=1,nz
        wrk1(iy,j,k)=0.
     end do
  end do
  do i=1,nyp
     rp=float(i-1)
     do k=1,nxh
        do j=1,nz
           wrk1(iy,j,k)=wrk1(iy,j,k)+rp*rp*(rp*rp-1.)*wrkc(i,j,k)*dyde*dyde/3.0
        end do
     end do
  end do

end subroutine c2dertw

!----------------------------------------------------------------------

subroutine cdrv1d(f,df)
  use grid_size

  real f(nyp),df(nyp)
  common/data2/ re,alpha,beta,xl,zl,yl,time,dyde,xstart,uinf

  df(nyp)=0.
  df(ny)=2.*float(ny)*f(nyp)*dyde
  do i=nym,1,-1
     df(i)=df(i+2)+2.*float(i)*f(i+1)*dyde
  end do
  df(1)=.5*df(1)

end subroutine cdrv1d


!----------------------------------------------------------------------


subroutine xyzfft(a,c,d,wfft1,wfft2,wfft3,wfft4,is)

!-- fast fourier  trans into spectral space
!-    if is = -1 physical -> spectral
!-    if is = 1  spectral -> physical
!-    else stop
  
  use grid_size

  real wfft3(1), wfft4(1)
  complex a(nyp,nz,nxh)
  real wfft1(1), wfft2(1)
  real sumre(nz),sumim(nz)
  complex c(nxh,nz), d(nz,nxh) 
  common/trigx/ trigx(2*nx), ixfax(19)
  common/trigz/ trigz(2*nz), izfax(19)
  common/trigy/ trigy(2*ny), sine(ny), cosine(ny), iyfax(19)

!--  physical to spectral

  if(is.eq.-1) then

!--  do x-z  planes, looping on y
     do j=1,nyp

        do k=1,nz
           do i=1,nxh
              c(i,k) = a(j,k,i) 
           end do
        end do

!--   x transform

        call rcs(c,wfft1,wfft2,nx,nz,ixfax,trigx)
        do i = 1,nxh
           do k = 1,nz
              d(k,i) = c(i,k)
           end do
        end do

!--   z transform

        call ccfft(d,wfft3,wfft4,wfft2,nz,nxh,-1,izfax,trigz)

        fac=1./(2.*float(ny)*float(nz))     ! needed for ccheb
        do i=1,nxh
           do k=1,nz
              a(j,k,i) = d(k,i)*fac
           end do
        end do

     end do

!--  finally y transform

     do i=1,nxh
        call ccheb(a(1,1,i),wfft1,wfft3,wfft4,wfft2,sumre,sumim,ny,nz,-1,iyfax,trigy,sine,cosine)
     end do
  else if(is.eq.1) then   !spectral to real space

     do i=1,nxh
        call ccheb(a(1,1,i),wfft1,wfft3,wfft4,wfft2,sumre,sumim,ny,nz,1,iyfax,trigy,sine,cosine)
     end do

     do j=1,nyp

        fac = 2.*float(ny)*float(nz)

        do i=1,nxh
           do k=1,nz
              d(k,i) = a(j,k,i)*fac
           end do
        end do

        call ccfft(d,wfft3,wfft4,wfft2,nz,nxh,1,izfax,trigz)

        do i = 1,nxh
           do k = 1,nz
              c(i,k) = d(k,i)
           end do
        end do

        call csr(c,wfft1,wfft2,nx,nz,ixfax,trigx)

        do i = 1,nxh
           do k = 1,nz
              a(j,k,i) = c(i,k)
           end do
        end do

     end do

  else
     write(*,110) is
110  format(' error: is = ',i5)
     stop
  end if

!--  results are in complex array a 
!-       & are properly ordered and scaled ( i think).

end subroutine xyzfft

!----------------------------------------------------------

subroutine yfft(a,wfft1,wfft2,wfft3,wfft4,is)

!-- fast fourier  trans into spectral space
!-    if is = -1 physical -> spectral
!-    if is = 1  spectral -> physical
!-    else stop

  use grid_size
  
  real wfft3(1), wfft4(1)

  complex a(nyp,nz,nxh)
  real wfft1(1), wfft2(1)
  real sumre(nz),sumim(nz)
  common/trigy/ trigy(2*ny), sine(ny), cosine(ny), iyfax(19)

!--  physical to spectral

  if(is.eq.-1) then

     fac=1./(2.*float(ny))     ! needed for ccheb

     do i=1,nxh
        do k=1,nz
           do j=1,nyp
              a(j,k,i) = a(j,k,i)*fac
           end do
        end do
     end do

!--  y transform

     do i=1,nxh
        call ccheb(a(1,1,i),wfft1,wfft3,wfft4,wfft2,sumre,sumim,ny,nz,-1,iyfax,trigy,sine,cosine)
     end do

  else if(is.eq.1) then   !spectral to real space

     do i=1,nxh
        call ccheb(a(1,1,i),wfft1,wfft3,wfft4,wfft2,sumre,sumim,ny,nz,1,iyfax,trigy,sine,cosine)
     end do

     fac = 2.*float(ny)

     do i=1,nxh
        do k=1,nz
           do j=1,nyp
              a(j,k,i) = a(j,k,i)*fac
           end do
        end do
     end do

  else
     write(*,110) is 
110  format(' error: is = ',i5)
     stop
  end if

!--  results are in complex array a 
!-       & are properly ordered and scaled ( i think).

end subroutine yfft

!----------------------------------------------------------

subroutine scram(a,b)
  use grid_size
  real a(nyp,nz,nx)
  real b(nyp2,nz,nxh)

  do j=1,nxh
     j1 = 2*j-1
     j2 = 2*j
     do i=1,nyp
        i1 = 2*i-1
        i2 = 2*i
        do k=1,nz
           b(i1,k,j) = a(i,k,j1)
           b(i2,k,j) = a(i,k,j2)
        end do
     end do
  end do

end subroutine scram

!----------------------------------------------------------

subroutine unscram(b,a)
  use grid_size
  real a(nyp,nz,nx)
  real b(nyp2,nz,nxh)

  do j=1,nxh
     j1 = 2*j-1
     j2 = 2*j
     do i=1,nyp
        i1 = 2*i-1
        i2 = 2*i
        do k=1,nz
           a(i,k,j1) = b(i1,k,j) 
           a(i,k,j2) = b(i2,k,j) 
        end do
     end do
  end do

end subroutine unscram


!----------------------------------------------------------
!--  subroutine solve()  version 1.2
!----------------------------------------------------------
!
!--  history: this code started as a testbed for some spectral-tau
!      algorithms for the solution of second order problems using
!      chebeshev polynomials. originally written by tom swean and
!      extensively modified by richard leighton. the original code 
!      included a pentadiagonal solver for robin boundary conditions
!      as well as a tridiagonal solver for dirichlet and neumann bc's.
!      the modifications were mainly to limit the scope to only the
!      tridiagonal systems, and to solve two dimensional problems
!
!----------------------------------------------------------
!
!--   the problem solved: example two dimensional (x parameterically)
!     time dependent heat equation, fourier transformed in the z-direction
!
!       du/dt-1/r*(d2 /dy2 -kz**2 -x) u = f(y,kz,x).
!   
!     with crank-nickolson timestepping yields
!
!       (d2 /dy2 -kz**2 -x -2*r/dt) u(n+1) = f(y,kz,x,n,n-1)
!
!     note: the coefficient to the second derivative term is one.
!
!----------------------------------------------------------
!
!--  the following are modifications to allow the rescaling of
!    the y-axis by a linear function. that is:
!
!              y = s(eta) = (dy/d(eta))*eta + const.
!
!       or     eta = t(y) = (d(eta)/dy)*y + const.
!
!      with dy/d(eta) and d(eta)/dy constant.
!
!    let v(eta) = v(t(y)) = u(y),
!    
!    then  dv/d(eta) = du/dy*dy/d(eta)
!    and   d2v/d(eta)2 = d2u/dy2*(dy/d(eta))**2.
!
!
!    therefore if we are interested in solving
!
!       dv/dt-1/r*(d2 /d(eta)2 -kz**2 -x) v(eta,...) = g(eta,kz,x),
!
!       where   eta0 <= eta <= eta1,
!
!    then a change of variables yields,
!
!      du/dt-1/r*((dy/d(eta)**2)*d2/dy2 -kz**2 -x) u(y,..) = f(y,kz,x),
!
!      where f(y,..) = g(t(y),..) = g(eta,..).
!   
!
!    therefore first derivatives must be multplied by the variable
!    dyde (= dy/d(eta) ), including those in the neumann boundary 
!    conditions.
!
!----------------------------------------------------------
!
!     --  subr. belonging to second order solver only --
!
!     solve( u,                            !input/output
!            x,g,dyde,                     !input
!            bctop,bcbot,ib,wavz,c,        !input
!            a,t,wrk)                      !input
!
!
!  general notes:
!
!     variables:
!        u:    martix containing rhs of discrete eqns before tri-diag on input
!              resulting solution in spectral space on output
!        x:    represtnts parameteric dependence on x, for 2-d x = 0,
!        g:    represents time dependency form time stepping routine
!        dyde: a scaling variable (y = dyde*eta + const ).
!               dyde = 1 for std domain ( -1<+y<=+1 ).
!        bctop:matrix containing bcs in fourier space at top boundary
!        bcbot:   "       "       "   "     "     "    " bottom boundary
!        ib:   integer flag for boundary conditions:
!               ib.eq.0:  dirichlet
!                   u(y=+1,kz) = bctop(kz), u(-1,kz)=bcbot(kz)
!               ib.ne.0:  neumann
!                   du/dy(y=1,kz) = bctop(kz), du/dy(y=-1,kz) = bcbot(kz)
!        wavz: vector containing wave numbers in z.
!        c:    chebychev "c" vector. c(1) = c(nyp) = 2, otherwise c() = 1
!        a:    matrix containing tri-diagonal part of tau d2/dy2 operator
!        t:    vector containing  elements of tri-diag repr boundary
!        wrkc: complex work matrix
!
!
!     u is both input, containing the source terms from sbr.rhs, and the 
!     output, containing the updated velocity field in cheb. space.
!
!
!     proceedure:
!
!         solves even and odd cheb coeffs. separately. 
!
!            evnc(a,x) and oddc(a,x) setup tridiagonal matrix which is indep. 
!              of time step, but are wavenumber dependent.
!            evns(u,wrkc) and odds(u,wrkc) use results from timestepping 
!              (u = source) and return even and odd portions of rhs 
!              of matrix problem in matrix wrkc. 
!            dircbc(t,nrank) sets up matrix elements of dirichlet b.c. equation.
!              dirichlet conditions are used and independent of wavenumber. 
!            neumbc(t,nrank,ip) set up matrix elements of neumann b.c. equation.
!              if ip = 0, for even part
!              if ip = 0, for odd part
!            eqsolv uses modified gaussian elimination to solve matrix problem
!
!
!     subroutines:
!
!       evnc/oddc- fills tri-diagonal part of the matrix 'a' in the eqn au=f. 
!       note that 'a' also contains the bc vector t, calculated in sbr bndry.
!
!       odds/evns- the rhs vector u has been passed to te solver, and these
!       combine functions of u(k+2),u(k), & u(k-2) in proper manner & puts 
!       result in wrkc.
!
!       eqsolv- solves system by gauss elimination. this solves nz equations
!       in parallel
!
!       dircbc/neumbc- fills the bc vector depending on the switch "ib". 
!       the boundary.
!       terms occuring on rhs of matrix problem are set in sbr setstuf
!
!          ---   note: this version solves two dim problem.
!                therefore nz physical locations become nz/2
!                complex fourier coefficients. in the 3-d version, 
!                real to conjg symm. transforms are used in x direction
!                so that nz spanwise locations become nz complex
!                fourier coefficients. ---
!
!----------------------------------------------------------

subroutine solve(u,x,g,dyde,bctop,bcbot,ib,wavz,c,a,t,wrkc)
  use grid_size
  complex u(nyp,nz),wrkc(nyp,nz)
  complex bctop(nz), bcbot(nz)
  real a(nyhp,nz,3),t(nyp)
  real wavz(nz),c(nyp),wn(nz)

!-----------------------------------------------------------------

  if(ib.eq.0) call dircbc(t,nyp)

  do k=1,nz
     wn(k) = wavz(k)**2 +x +g
  end do

  kzero = 0
  if(ib.ne.0) then         ! neumann
     do k=1,nz
        if(wn(k).le.1.0e-20) then
           kzero = k
           wn(k) = 1.0
        end if
     end do
  end if

!     even problem

  ip = 0
  nrank = nyhp
  call evnc(a,wn,dyde)
  call evns(u,bctop,bcbot,ib,c,wrkc)
  if(ib.ne.0) call neumbc(t,dyde,nrank,ip)
  call eqsolv(a,t,wrkc,nrank)

!     update even values of u

  do j=1,nyp,2
     i=(j+1)/2
     do k=1,nz
        u(j,k)=wrkc(i,k)
     end do
  end do

!      odd problem

  ip = 1
  nrank = nyh
  call oddc(a,wn,dyde)
  call odds(u,bctop,bcbot,ib,c,wrkc)
  if(ib.ne.0) call neumbc(t,dyde,nrank,ip)
  call eqsolv(a,t,wrkc,nrank)

!--  update odd values of u

  do j=2,ny,2
     i=j/2
     do k=1,nz
        u(j,k)=wrkc(i,k)
     end do
  end do

  if((ib.ne.0).and.(kzero.ne.0)) then
     do j=1,nyp
        u(j,kzero) = 0.0
     end do
  end if

end subroutine solve

!----------------------------------------------------------

subroutine dircbc(t,n)

  real t(n)

  do j=1,n
     t(j)=1.
  end do

end subroutine dircbc

!---------------------------------------------------------

subroutine neumbc(t,dyde,n,ip)

  real t(n)

  if(ip.eq.0) then
     do j=1,n                      !k=0
        t(j)=dyde*(2*(j-1))**2           !k*k
     end do                        !k=k+2
  else
     do j=1,n                      !k=1
        t(j)= dyde*(2*(j-1) +1)**2       !k*k
     end do                        !k=k+2
  end if

end subroutine neumbc

!----------------------------------------------------------

subroutine evnc (a,wn,dyde)

!-- sets up 'a' matrix  even (ns=2) problems
!   nabs is the c(i) factor, c(n) = 2,c(i) =1 otherwise
!  'a' does not include the boundary condition
  
  use grid_size
  real a(nyhp,nz,3)
  real wn(nz)

  do i=1,nyh
     j = 2*i
     rj1=1./float(4*j*(j-1))
     rj2=1./float(2*(j*j-1))
     rj3=1./float(4*j*(j+1))
     do k = 1,nz
        a(i,k,1)= wn(k)*rj1
        a(i,k,2)=-1.*dyde*dyde-wn(k)*rj2
        a(i,k,3)= wn(k)*rj3
     end do
  end do

  do k = 1,nz
     a(1,k,1) = 2.*a(1,k,1)
     a(nyhm,k,3)=0.0
     a(nyh,k,2)=-1.0*dyde*dyde
     a(nyh,k,3)=0.0
  end do

end subroutine evnc

!----------------------------------------------------------

subroutine oddc (a,wn,dyde)

!-- sets up 'a' matrix for odd (ns=3) problems
!  'a' does not include the boundary condition

  use grid_size
  real a(nyhp,nz,3)
  real wn(nz)

  do i = 1,nyhm
     j=2*i+1
     rj1=1./float(4*j*(j-1))
     rj2=1./float(2*(j*j-1))
     rj3=1./float(4*j*(j+1))
     do k=1,nz
        a(i,k,1)=wn(k)*rj1
        a(i,k,2)=-1.*dyde*dyde-wn(k)*rj2
        a(i,k,3)=wn(k)*rj3
     end do
  end do

  do k=1,nz
     a(nyhm-1,k,3)=0.0
     a(nyhm,k,2)=  -1.0*dyde*dyde
     a(nyhm,k,3)=  0.0
  end do

end subroutine oddc

!----------------------------------------------------------

subroutine evns (s,bctop,bcbot,ib,c,wrk)

!--  sets up rhs for even tri-diag problem 
!-    given results of sbr rhs
!-   *** note use of zero based arrays ***

  use grid_size
  complex s(0:ny,nz),wrk(nyp,nz)
  complex bctop(nz), bcbot(nz)
  real c(0:ny)

  do i = 2,nyhm
     j = 2*i -2
     jp2 = j+2
     jm2 = j-2
     rj1=1./float(4*j*(j-1))
     rj2=1./float(2*(j*j-1))
     rj3=1./float(4*j*(j+1))
     do k = 1,nz
        wrk(i,k) = - c(jm2)*s(jm2,k)*rj1 + s(j,k)*rj2 - s(jp2,k)*rj3
     end do
  end do

  j = ny -2
  jm2 = ny -4
  rj1=1./float(4*j*(j-1))
  rj2=1./float(2*(j*j-1))
  do k = 1,nz
     wrk(nyh,k) = -s(jm2,k)*rj1+s(j,k)*rj2
  end do

  j = ny
  jm2 = ny -2
  rj1=1./float(4*j*(j-1))
  do k = 1,nz
     wrk(nyhp,k)= -s(jm2,k)*rj1
  end do

! bc's
  
  if(ib.eq.0) then
     do k=1,nz
        wrk(1,k) = 0.5*(bctop(k)+bcbot(k))
     end do
  else
     do k=1,nz
        wrk(1,k) = 0.5*(bctop(k)-bcbot(k))
     end do
  end if
  
end subroutine evns

!----------------------------------------------------------

subroutine odds (s,bctop,bcbot,ib,c,wrk)

!--  sets up rhs for odd tri-diag problem 
!-    given results of sbr rhs

  use grid_size
  complex s(0:ny,nz), wrk(nyp,nz)
  complex bctop(nz), bcbot(nz)
  real c(nyp)

  do i = 2,nyh-2
     j = 2*i -1
     jp2 = j+2
     jm2 = j-2
     rj1=1./float(4*j*(j-1))
     rj2=1./float(2*(j*j-1))
     rj3=1./float(4*j*(j+1))
     do k = 1,nz
        wrk(i,k) =  - s(jm2,k)*rj1 + s(j,k)*rj2 - s(jp2,k)*rj3
     end do
  end do

  i = nyhm
  j = ny -3
  jm2 = ny -5
  rj1=1./float(4*j*(j-1))
  rj2=1./float(2*(j*j-1))
  do k = 1,nz
     wrk(i,k) = - s(jm2,k)*rj1 + s(j,k)*rj2
  end do
  i = nyh
  j = ny-1
  jm2 = ny -3
  rj1=1./float(4*j*(j-1))
  do k = 1,nz
     wrk(i,k)= - s(jm2,k)*rj1
  end do

! bc's

  if(ib.eq.0) then
     do k=1,nz
        wrk(1,k) = 0.5*(bctop(k)-bcbot(k))
     end do
  else
     do k=1,nz
        wrk(1,k) = 0.5*(bctop(k)+bcbot(k))
     end do
  end if

end subroutine odds

!----------------------------------------------------------

subroutine eqsolv (a,t,wrk,nrk)

!--  solves odd or even problem by gaussian elimination
  
  use grid_size
  complex wrk(nyp,nz)
  real d(nyp,nz),f(nz)
  real a(nyhp,nz,3),t(nyp)

  m=nrk-1
  mm1=m-1
  mm2=m-2

  do k = 1,nz
     d(m,k)=a(m,k,2)
  end do
  do i=mm1,1,-1
     do k = 1,nz
        d(i,k)=a(i,k,2)-a(i,k,3)*a(i+1,k,1)/d(i+1,k)
        wrk(i+1,k)=wrk(i+1,k)-a(i,k,3)*wrk(i+2,k)/d(i+1,k)
     end do
  end do

!     eliminate top row

  do k = 1,nz
     f(k)=t(m)-t(nrk)*a(m,k,1)/d(m,k)
     wrk(1,k)=wrk(1,k)-t(nrk)*wrk(nrk,k)/d(m,k)
  end do
  do i=mm1,1,-1
     do k = 1,nz
        wrk(1,k)=wrk(1,k)-f(k)*wrk(i+1,k)/d(i,k)
        f(k)=t(i)-f(k)*a(i,k,1)/d(i,k)
     end do
  end do

!     back substitute: put sol'n in wrk

  do k = 1,nz
     wrk(1,k)=wrk(1,k)/f(k)
  end do
  do i=2,nrk
     do k = 1,nz
        wrk(i,k)=(wrk(i,k)-a(i-1,k,1)*wrk(i-1,k))/d(i-1,k)
     end do
  end do
  
end subroutine eqsolv

!----------------------------------------------------------
!
!----------------------------------------------------------
!-
!--  subroutine penta()  version 1.2
!-
!----------------------------------------------------------
!-
!--  history: this code started as a testbed for some spectral-tau
!      algorithms for the solution of second order problems using
!      chebeshev polynomials. originally written by tom swean and
!      extensively modified by richard leighton. the original code 
!      included a pentadiagonal solver for robin boundary conditions
!      as well as a tridiagonal solver for dirichlet and neumann bc's.
!      the modifications were mainly to limit the scope to only the
!      pentadiagonal systems, and to solve two dimensional problems
!
!----------------------------------------------------------
!
!--   the problem solved: example two diminsional (x parameterically)
!     time dependent heat equation, fourier transformed in the z-direction
!
!       du/dt-1/r*(d2 /dy2 -kz**2 -x) u = f(y,kz,x).
!   
!     with crank-nickolson timestepping yields
!
!       (d2 /dy2 -kz**2 -x -2*r/dt) u(n+1) = f(y,kz,x,n,n-1)
!
!            g = 2*r/dt
!
!     note: the coefficient to the second derivative term is one.
!
!----------------------------------------------------------
!
!  the boundary condition being used is :
!
!         at*u(y=1,kz)  + bt*du/dy(y=1,kz)  = gtop(kz)
!
!         ab*u(y=-1,kz) + bb*du/dy(y=-1,kz) = gbot(kz)
!
!----------------------------------------------------------
!
!--  the following are modifications to allow the rescaling of
!    the y-axis by a linear function. that is:
!
!              y = s(eta) = (dy/d(eta))*eta + const.
!
!       or     eta = t(y) = (d(eta)/dy)*y + const.
!
!      with dy/d(eta) and d(eta)/dy constant.
!
!    let v(eta) = v(t(y)) = u(y),
!    
!    then  dv/d(eta) = du/dy*dy/d(eta)
!    and   d2v/d(eta)2 = d2u/dy2*(dy/d(eta))**2.
!
!
!    therefore if we are interested in solving
!
!       dv/dt-1/r*(d2 /d(eta)2 -kz**2 -x) v(eta,...) = g(eta,kz,x),
!
!       where   eta0 <= eta <= eta1,
!
!    then a change of variables yields,
!
!   du/dt-1/r*(( (dy/d(eta))**2 )*d2/dy2 -kz**2 -x) u(y,..) = f(y,kz,x),
!
!      where f(y,..) = g(t(y),..) = g(eta,..).
!   
!
!    therefore first derivatives must be multplied by the variable
!    dyde (= dy/d(eta) ), including those in the robin boundary 
!    conditions.
!
!----------------------------------------------------------
!
!      subroutine pent(   u,                       ! output/input
!     1                   x,g,dyde,                ! input
!     2                   ib,wavz,c,               ! input
!     3                   at,bt,gtop,              ! input
!     4                   ab,bb,gbot,              ! input
!     5                   a,t,wrkc)                ! work space
!
!----------------------------------------------------------
!
!  general notes:
!
!     variables:
!        u:    martix containing rhs of discrete eqns before penta-diag 
!              on input resulting solution in spectral space on output
!        x:    represents parameteric dependence on x, for 2-d x = 0,
!        g:    represents time dependency form time stepping routine
!        dyde: a scaling factor in the y-direction. is the usual scale is
!              used (-1<=y<=1 ), dyde = 1.
!        ib:   ib.eq.0 for for general robin conditions, excluding neumann
!              ib.en.0 for neumann condtions
!        wavz: vector of wavenumbers in z direction
!        c:    chebychev "c" vector
!        at,bt,gtop,ab,bb,gbot: are boundary condition terms in fourier 
!              space. see above.
!        a:    matrix containing tri-diagonal part of tau d2/dy2 operator
!        t:    vector containing  elements of tri-diag repr boundary
!        wrkc: complex work martix
!
!
!     u is both input, containing the source terms from sbr.rhs, and the 
!     output, containing the updated velocity field in cheb. space.
!
!     proceedure:
!
!         solves even and odd cheb coeffs. separately, combining them
!          in the second call to pntsolv
!
!            pntevnc(a,x) and pntoddc(a,x) setup tridiagonal matrix which 
!              is indep. of time step, but are wavenumber dependent.
!            pntevns(u,wrkc) and pntodds(u,wrkc) use results from timestepping 
!              (u = source) and return even and odd portions of rhs 
!              of matrix problem in matrix wrkc. 
!            mixedbc(  ) sets up all elements of the b.c.
!            fillt ( )  sets up the t-matrix depending on parameters ip and
!              ipass
!            pntsolv uses modified gaussian elimination to solve matrix problem
!
!
!       subroutines:
!
!       pntevnc/pntoddc- fills tri-diagonal part of the matrix 'a' in the eqn 
!       au=f. note that 'a' also contains the bc vector t, calculated in 
!       sbrs fillt and mixedbc.
!
!       pntodds/pntevns- the rhs vector u has been passed to the solver,
!       and these combine functions of u(k+2),u(k), & u(k-2) in proper 
!       manner & puts result in wrkc.
!
!       pntsolv- solves system by gauss elimination. this solves nz equations
!       in parallel
!
!       terms occuring on rhs of matrix problem are set in sbr setstuf
!
!          ---   note: this version solves two dim problem.
!                therefore nz physical locations become nz/2
!                complex fourier coefficients. in the 3-d version, 
!                real to conjg symm. transforms are used in x direction
!                so that nz spanwise locations become nz complex
!                fourier coefficients. ---
!
!
!
!-
!----------------------------------------------------------

subroutine penta(u,x,g,dyde,ib,wavz,c,at,bt,gtop,ab,bb,gbot,a,t,wrkc)

  use grid_size
  complex u(nyp,nz),wrkc(nyp,nz)
  complex gtop(nz), gbot(nz)
  real a(nyhp,nz,3),t(nyp)
  real wavz(nz),c(nyp),wn(nz)

!-----------------------------------------------------------------

  call mixedbc(at,ab,bt,bb,dyde)

  do k=1,nz
     wn(k) = wavz(k)**2 +x +g
  end do

  kzero = 0
  if(ib.ne.0) then       ! strictly neumann
     do k=1,nz
        if(wn(k).le.1.0e-20) then
           kzero = k
           wn(k) = 1.0
        end if
     end do
  end if

!     even problem
  
  ip = 1
  nrank = nyhp
  call pntevnc(a,wn,dyde)
  call pntevns(u,ib,c,wrkc)
  call pntsolv(a,t,wrkc,gtop,gbot,ip,nrank)

!      odd problem
  
  ip = 2
  nrank = nyh
  call pntoddc(a,wn,dyde)
  call pntodds(u,ib,c,wrkc)
  call pntsolv(a,t,wrkc,gtop,gbot,ip,nrank)

!--  update u

  do k=1,nz
     do j=1,nyp
        u(j,k)=wrkc(j,k)
     end do
  end do
  if((ib.ne.0).and.(kzero.ne.0)) then
     do j=1,nyp
        u(j,kzero) = 0.0
     end do
  end if

end subroutine penta

!----------------------------------------------------------

subroutine pntevnc (a,wn,dyde)

!-- sets up 'a' matrix  even (ns=2) problems
!   nabs is the c(i) factor, c(n) = 2,c(i) =1 otherwise
!  'a' does not include the boundary condition
  
  use grid_size
  real a(nyhp,nz,3)
  real wn(nz)

  do i=1,nyh
     j = 2*i
     rj1=1./float(4*j*(j-1))
     rj2=1./float(2*(j*j-1))
     rj3=1./float(4*j*(j+1))
     do k = 1,nz
        a(i,k,1)= wn(k)*rj1
        a(i,k,2)=-1.*dyde*dyde-wn(k)*rj2
        a(i,k,3)= wn(k)*rj3
     end do
  end do

  do k = 1,nz
     a(1,k,1) = 2.*a(1,k,1)
     a(nyhm,k,3)=0.0
     a(nyh,k,2)=-1.0*dyde*dyde
     a(nyh,k,3)=0.0
  end do

end subroutine pntevnc

!----------------------------------------------------------

subroutine pntoddc (a,wn,dyde)

!-- sets up 'a' matrix for odd (ns=3) problems
!  'a' does not include the boundary condition

  use grid_size
  real a(nyhp,nz,3)
  real wn(nz)

  do i = 1,nyhm
     j=2*i+1
     rj1=1./float(4*j*(j-1))
     rj2=1./float(2*(j*j-1))
     rj3=1./float(4*j*(j+1))
     do k=1,nz
        a(i,k,1)=wn(k)*rj1
        a(i,k,2)=-1.*dyde*dyde-wn(k)*rj2
        a(i,k,3)=wn(k)*rj3
     end do
  end do

  do k=1,nz
     a(nyhm-1,k,3)=0.0
     a(nyhm,k,2)=  -1.0*dyde*dyde
     a(nyhm,k,3)=  0.0
  end do

end subroutine pntoddc

!----------------------------------------------------------

subroutine pntevns (s,ib,c,wrk)

!--  sets up rhs for even tri-diag problem 
!-    given results of sbr rhs
!-   *** note use of zero based arrays ***

  use grid_size
  complex s(0:ny,nz),wrk(nyp,nz)
  real c(0:ny)

  do i = 2,nyhm
     j = 2*i -2
     jp2 = j+2
     jm2 = j-2
     rj1=1./float(4*j*(j-1))
     rj2=1./float(2*(j*j-1))
     rj3=1./float(4*j*(j+1))
     do k = 1,nz
        wrk(i,k) = - c(jm2)*s(jm2,k)*rj1 + s(j,k)*rj2 - s(jp2,k)*rj3
     end do
  end do

  j = ny -2
  jm2 = ny -4
  rj1=1./float(4*j*(j-1))
  rj2=1./float(2*(j*j-1))
  do k = 1,nz
     wrk(nyh,k) = -s(jm2,k)*rj1+s(j,k)*rj2
  end do

  j = ny
  jm2 = ny -2
  rj1=1./float(4*j*(j-1))
  do k = 1,nz
     wrk(nyhp,k)= -s(jm2,k)*rj1
  end do

end subroutine pntevns

!----------------------------------------------------------

subroutine pntodds (s,ib,c,wrk)

!-    sets up rhs for odd tri-diag problem 
!-    given results of sbr rhs

  use grid_size
  complex s(0:ny,nz), wrk(nyp,nz)
  real c(nyp)

  do i = 2,nyh-2
     j = 2*i -1
     jp2 = j+2
     jm2 = j-2
     rj1=1./float(4*j*(j-1))
     rj2=1./float(2*(j*j-1))
     rj3=1./float(4*j*(j+1))
     do k = 1,nz
        wrk(i,k) =  - s(jm2,k)*rj1 + s(j,k)*rj2 - s(jp2,k)*rj3
     end do
  end do

  i = nyhm
  j = ny -3
  jm2 = ny -5
  rj1=1./float(4*j*(j-1))
  rj2=1./float(2*(j*j-1))
  do k = 1,nz
     wrk(i,k) = - s(jm2,k)*rj1 + s(j,k)*rj2
  end do

  i = nyh
  j = ny-1
  jm2 = ny -3
  rj1=1./float(4*j*(j-1))
  do k = 1,nz
     wrk(i,k)= - s(jm2,k)*rj1
  end do

end subroutine pntodds

!----------------------------------------------------------

subroutine mixedbc(at,ab,bt,bb,dyde)
  use grid_size
  real at,ab,bt,bb
  common /robbin/ tt(nyp),tb(nyp)

  do j=1,nyp
     j1=j-1
     tt(j)=at+dyde*bt*float(j1*j1)
     tb(j)=ab-dyde*bb*float(j1*j1)
  end do

  do j=2,nyp,2
     tb(j)=-tb(j)
  end do

end subroutine mixedbc

!----------------------------------------------------------

subroutine fillt(nrk,ip,ipass,t)
  use grid_size
  common /robbin/ tt(nyp),tb(nyp)
  dimension t(nyp)

  if(ip.eq.1.and.ipass.eq.1) then        ! even, first pass
     do j=1,nrk
        k=2*j-1
        t(j)=tb(k)
     end do
  else if(ip.eq.1.and.ipass.eq.2) then    ! even, second pass
     do j=1,nrk
        k=2*j-1
        t(j)=tt(k)
     end do
    else if(ip.eq.2.and.ipass.eq.1) then    ! odd, first pass
     do j=1,nrk
        k=2*j
        t(j)=tb(k)
     end do
  else                                   ! odd, second pass
     do j=1,nrk
        k=2*j
        t(j)=tt(k)
     end do
  end if

end subroutine fillt

!----------------------------------------------------------

subroutine pntsolv (a,t,wrk,gtop,gbot,ip,nrk)
  use grid_size
  complex wrk(nyp,nz), gtop(nz), gbot(nz)
  real a(nyhp,nz,3), t(nyp)
  complex s(2,nz), gs(nyhp,nz,2)
  real cs(nyhp,nz,2),f(2,nz,2),ds(nyhp,nz,2)
  save s, gs, cs, f, ds

  ipass=1
  if(ip.eq.1) then
     do k=1,nz
        s(1,k)=gtop(k)
        s(2,k)=gbot(k)
     end do
  end if

  js=2
  n = nrk
  m=n-1
  mm1=m-1
  mm2=m-2

  do i=1,m
     do k=1,nz
        gs(i,k,ip)=wrk(i,k)
        cs(i,k,ip)=a(i,k,1)
     end do
  end do

  do k=1,nz
     ds(m,k,ip)=a(m,k,2)
     gs(n,k,ip)=wrk(n,k)
  end do

  do i=mm1,1,-1
     do k=1,nz
        ds(i,k,ip)=a(i,k,2)-a(i,k,3)*cs(i+1,k,ip)/ds(i+1,k,ip)
        gs(i+1,k,ip)=gs(i+1,k,ip)-a(i,k,3)*gs(i+2,k,ip)/ds(i+1,k,ip)
     end do
  end do

!     eliminate 2nd row on ipass=1 & 1st row on ipass=2

2 continue

  call fillt(n,ip,ipass,t)

  do k=1,nz
     f(js,k,ip)=t(m)-t(n)*cs(m,k,ip)/ds(m,k,ip)
     s(js,k)=s(js,k)-t(n)*gs(n,k,ip)/ds(m,k,ip)
  end do

  do i=mm1,1,-1
     do k=1,nz
        s(js,k)=s(js,k)-f(js,k,ip)*gs(i+1,k,ip)/ds(i,k,ip)
        f(js,k,ip)=t(i)-f(js,k,ip)*cs(i,k,ip)/ds(i,k,ip)
     end do
  end do

  if(ipass.eq.1) then
     ipass=ipass+1
     js=js-1
     goto 2
  end if

!     solve the 2 eqns. f(1,k,1)*u0+f(1,k,2)*u1=s1; f(2,k,1)*u0+f(2,k,2)*u1=s2
!     put the sol'n in wrk
      
  if(ip.eq.2.and.ipass.eq.2) then

     do k=1,nz
        wrk(1,k)=(s(1,k)*f(2,k,2)-s(2,k)*f(1,k,2))/(f(1,k,1)*f(2,k,2)-f(2,k,1)*f(1,k,2))
        wrk(2,k)=(s(2,k)*f(1,k,1)-s(1,k)*f(2,k,1))/(f(1,k,1)*f(2,k,2)-f(2,k,1)*f(1,k,2))
     end do

!     even u
         
     ne=n+1
     do i=2,ne
        ii=2*i -1
        do k=1,nz
           wrk(ii,k)=(gs(i,k,1)-cs(i-1,k,1)*wrk(ii-2,k))/ds(i-1,k,1)
        end do
     end do

!     odd u

     do i=2,n
        ii=2*i
        do k=1,nz
           wrk(ii,k)=(gs(i,k,2)-cs(i-1,k,2)*wrk(ii-2,k))/ds(i-1,k,2)
        end do
     end do

  end if
  
end subroutine pntsolv

! ----------------------------------------------------------
! -
! -- subroutine penta1d()  version 1.1
! -
! ----------------------------------------------------------
! -
! --   this is a one dimensional version of penta.sbr, it is also 
!      entirely real
! ------------------------------------------------------------
! 
! --  history: this code started as a testbed for some spectral-tau
!       algorithms for the solution of second order problems using
!       chebeshev polynomials. originally written by tom swean and
!       extensively modified by richard leighton. the original code 
!       included a pentadiagonal solver for robin boundary conditions
!       as well as a tridiagonal solver for dirichlet and neumann bc's.
! 
! --  note: the algorithm and subroutine names are simlar to those in 
!       penta.sbr and solve.sbr, but they are not the same. 
!       they are not interchangabler.
! ----------------------------------------------------------
! 
! --   the problem solved: example one dimensional (x parameterically)
!      time dependent heat equation, fourier transformed in the z-direction
! 
!        du/dt-1/r*(d2 /dy2 -x) u = f(y,kz=0,x).
!    
!      with crank-nickolson timestepping yields
! 
!        (d2 /dy2 -2*r/dt -x) u0(n+1) = f(y,kz=0,x,n,n-1)
! 
!          g = 2*r/dt
! 
!      note: the coefficient to the second derivative term is one.
! 
! ----------------------------------------------------------
! 
!   the boundary condition being used is :c
!          at*u(y=1)  + bt*du/dy(y=1)  = gtop(kz=0)
! 
!          ab*u(y=-1) + bb*du/dy(y=-1) = gbot(kz=0)
! 
! ----------------------------------------------------------
! 
! --  the following are modifications to allow the rescaling of
!     the y-axis by a linear function. that is:
! 
!               y = s(eta) = (dy/d(eta))*eta + const.
! 
!        or     eta = t(y) = (d(eta)/dy)*y + const.
! 
!       with dy/d(eta) and d(eta)/dy constant.
! 
!     let v(eta) = v(t(y)) = u(y),
!     
!     then  dv/d(eta) = du/dy*dy/d(eta)
!     and   d2v/d(eta)2 = d2u/dy2*(dy/d(eta))**2.
! 
! 
!     therefore if we are interested in solving
! 
!        dv/dt-1/r*(d2 /d(eta)2 -kz**2 -x) v(eta,...) = g(eta,kz,x),
! 
!        where   eta0 <= eta <= eta1,
! 
!     then a change of variables yields,
! 
!       du/dt-1/r*((dy/d(eta)**2)*d2/dy2 -kz**2 -x) u(y,..) = f(y,kz,x),
! 
!       where f(y,..) = g(t(y),..) = g(eta,..).
!    
! 
!     therefore first derivatives must be multplied by the variable
!     dyde (= dy/d(eta) ), including those in the robin boundary 
!     conditions.
! 
! ----------------------------------------------------------
! 
! 
!      psolv1d( u,                          !input/output
!             g,dyde,                       !input
!             at, bt, gtop,                 !input
!             ab, bb, gbot,                 !input
!             c,a,t,wrk)                    !work space
! 
! 
!   general notes:
! 
!      variables:
!         u:    martix containing rhs of discrete eqns before penta-diag 
!               on input resulting solution in spectral space on output
!         x:    represtnts parameteric dependence on x, for 2-d x = 0,
!         g:    represents time dependency form time stepping routine
!         dyde: is a scalefactor in the y-direction.
!         at,bt,gtop,ab,bb,gbot: bc terms, see above.
!         wavz: spanwise wavenumber, left internal to code, but not used.
!         a:    matrix containing tri-diagonal part of tau d2/dy2 operator
!         t:    vector containing  elements of tri-diag repr boundary
!         wrkc: real work martix
! 
!      u is both input, containing the source terms from sbr.rhs, and the 
!      output, containing the updated velocity field in cheb. space.
! 
!      proceedure:
! 
!          solves even and odd cheb coeffs. separately. 
! 
!             p1devnc(a,x) and p1doddc(a,x) setup tridiagonal matrix which is indep. 
!               of time step, but are wavenumber dependent.
!             p1devns(u,wrkc) and p1dodds(u,wrkc) use results from timestepping 
!               (u = source) and return even and odd portions of rhs 
!               of matrix problem in matrix wrkc. 
!             mxdbc1d(  ) sets up all elements of the b.c.
!             fillt1d( )  sets up the t-matrix depending on parameters ip and
!               ipass
!             pntsolv uses modified gaussian elimination to solve matrix problem
! 
! 
!        subroutines
! 
!        evnc/oddc- fills tri-diagonal part of the matrix 'a' in the eqn au=f. 
!        note that 'a' also contains the bc vector t, calculated in sbr bndry.
! 
!        odds/evns- the rhs vector u has been passed to te solver, and these
!        combine functions of u(k+2),u(k), & u(k-2) in proper manner & puts 
!        result in wrkc.
! 
!        pntsolv- solves system by gauss elimination.
! 
!        terms occuring on rhs of matrix problem are set in sbr setstu
! 
!           ---   note: this version solves one dim problem only
!                 and is entirely real. ---
! 
! -
! ----------------------------------------------------------
! -

subroutine psolv1d ( u,g,dyde,at, bt, gtop,ab, bb, gbot,c,a,t,wrk)
  use grid_size
  real u(nyp),wrkc(nyp),wrk(nyp)
  real gtop, gbot
  real a(nyhp,3),t(nyp)
  real wavz,c(nyp)
  data wavz / 0.0 /

!-----------------------------------------------------------------
  
  x = 0.0
  call mxdbc1d(at,ab,bt,bb,dyde)

!     even problem
  
  ip = 1
  nrank = nyhp
  call p1devnc(a,x,g,wavz,dyde)
  call p1devns(u,c,wrkc)
  call p1dsolv(a,t,wrkc,gtop,gbot,ip,nrank)

!     odd problem
  
  ip = 2
  nrank = nyh
  call p1doddc(a,x,g,wavz,dyde)
  call p1dodds(u,c,wrkc)
  call p1dsolv(a,t,wrkc,gtop,gbot,ip,nrank)

!     update u
  
  do j=1,nyp
     u(j)=wrkc(j)
  end do
  
end subroutine psolv1d

!----------------------------------------------------------

subroutine p1devnc (a,x,g,wavz,dyde)

!-- sets up 'a' matrix  even (ns=2) problems
!   nabs is the c(i) factor, c(n) = 2,c(i) =1 otherwise
!  'a' does not include the boundary condition

  use grid_size
  real a(nyhp,3)
  real wavz

  do i=1,nyh
     j = 2*i
     rj1=1./float(4*j*(j-1))
     rj2=1./float(2*(j*j-1))
     rj3=1./float(4*j*(j+1))
     wn = wavz**2 +x +g
     a(i,1)= wn*rj1
     a(i,2)=-1.*dyde*dyde-wn*rj2
     a(i,3)= wn*rj3
  end do

  a(1,1) = 2.*a(1,1)
  a(nyhm,3)=0.0
  a(nyh,2)=-1.0*dyde*dyde
  a(nyh,3)=0.0
  
end subroutine p1devnc

!----------------------------------------------------------

subroutine p1doddc (a,x,g,wavz,dyde)

!-- sets up 'a' matrix for odd (ns=3) problems
!  'a' does not include the boundary condition
  
  use grid_size
  real a(nyhp,3)
  real wavz

  do i = 1,nyhm
     j=2*i+1
     rj1=1./float(4*j*(j-1))
     rj2=1./float(2*(j*j-1))
     rj3=1./float(4*j*(j+1))
     wn = wavz**2 +x +g
     a(i,1)=wn*rj1
     a(i,2)=-1.*dyde*dyde-wn*rj2
     a(i,3)=wn*rj3
  end do

  a(nyhm-1,3)=0.0
  a(nyhm,2)=  -1.0*dyde*dyde
  a(nyhm,3)=  0.0
  
end subroutine p1doddc

!----------------------------------------------------------

subroutine p1devns (s,c,wrk)

!--  sets up rhs for even tri-diag problem 
!-    given results of sbr rhs
!-   *** note use of zero based arrays ***
  
  use grid_size
  real s(0:ny),wrk(nyp)
  real c(0:ny)

  do i = 2,nyhm
     j = 2*i -2
     jp2 = j+2
     jm2 = j-2
     rj1=1./float(4*j*(j-1))
     rj2=1./float(2*(j*j-1))
     rj3=1./float(4*j*(j+1))
     wrk(i) = - c(jm2)*s(jm2)*rj1 + s(j)*rj2  - s(jp2)*rj3
  end do

  j = ny-2
  jm2 = ny-4
  rj1=1./float(4*j*(j-1))
  rj2=1./float(2*(j*j-1))
  wrk(nyh) = -s(jm2)*rj1 + s(j)*rj2

  j = ny
  jm2 = ny-2
  rj1=1./float(4*j*(j-1))
  wrk(nyhp)= -s(jm2)*rj1

end subroutine p1devns

!----------------------------------------------------------

subroutine p1dodds (s,c,wrk)

!-    sets up rhs for odd tri-diag problem 
!-    given results of sbr rhs

  use grid_size
  real s(0:ny), wrk(nyp)
  real c(nyp)

  do i = 2,nyh-2
     j = 2*i -1
     jp2 = j+2
     jm2 = j-2
     rj1=1./float(4*j*(j-1))
     rj2=1./float(2*(j*j-1))
     rj3=1./float(4*j*(j+1))
     wrk(i) =  - s(jm2)*rj1 + s(j)*rj2 - s(jp2)*rj3
  end do

  i = nyhm
  j = ny -3
  jm2 = ny -5
  rj1=1./float(4*j*(j-1))
  rj2=1./float(2*(j*j-1))
  wrk(i) = - s(jm2)*rj1 + s(j)*rj2

  i = nyh
  j = ny-1
  jm2 = ny -3
  rj1=1./float(4*j*(j-1))
  wrk(i)= - s(jm2)*rj1
  
end subroutine p1dodds

!----------------------------------------------------------

subroutine mxdbc1d (at,ab,bt,bb,dyde)

  use grid_size
  real at,ab,bt,bb
  common /robbin/ tt(nyp),tb(nyp)

  do j=1,nyp
     j1=j-1
     tt(j)=at+dyde*bt*float(j1*j1)
     tb(j)=ab-dyde*bb*float(j1*j1)
  end do

  do j=2,nyp,2
     tb(j)=-tb(j)
  end do

end subroutine mxdbc1d

!----------------------------------------------------------

subroutine fillt1d(nrk,ip,ipass,t)
  use grid_size
  common /robbin/ tt(nyp),tb(nyp)
  dimension t(nyp)

  if(ip.eq.1.and.ipass.eq.1) then        ! even, first pass
     do j=1,nrk
        k=2*j-1
        t(j)=tb(k)
     end do
  else if(ip.eq.1.and.ipass.eq.2) then    ! even, second pass
     do j=1,nrk
        k=2*j-1
        t(j)=tt(k)
     end do
  else if(ip.eq.2.and.ipass.eq.1) then    ! odd, first pass
     do j=1,nrk
        k=2*j
        t(j)=tb(k)
     end do
  else                                   ! odd, second pass
     do j=1,nrk
        k=2*j
        t(j)=tt(k)
     end do
  end if

end subroutine fillt1d

!----------------------------------------------------------

subroutine p1dsolv (a,t,wrk,gamt,gamb,ip,nrk)
  use grid_size
  real wrk(nyp), gamb, gamt
  real a(nyhp,3), t(nyp)
  real s(2), gs(nyhp,2)
  real cs(nyhp,2),f(2,2),ds(nyhp,2)
  save s, gs, cs, f, ds

  ipass=1
  if(ip.eq.1) then
     s(1)=gamt
     s(2)=gamb
  end if

  js=2
  n = nrk
  m=n-1
  mm1=m-1
  mm2=m-2

  do i=1,m
     gs(i,ip)=wrk(i)
     cs(i,ip)=a(i,1)
  end do

  ds(m,ip)=a(m,2)
  gs(n,ip)=wrk(n)

  do i=mm1,1,-1
     ds(i,ip)=a(i,2)-a(i,3)*cs(i+1,ip)/ds(i+1,ip)
     gs(i+1,ip)=gs(i+1,ip)-a(i,3)*gs(i+2,ip)/ds(i+1,ip)
  end do

!     eliminate 2nd row on ipass=1 & 1st row on ipass=2

2 continue

  call fillt1d(n,ip,ipass,t)

  f(js,ip) = t(m)-t(n)*cs(m,ip)/ds(m,ip)
  s(js) = s(js)-t(n)*gs(n,ip)/ds(m,ip)

  do i=mm1,1,-1
     s(js)=s(js)-f(js,ip)*gs(i+1,ip)/ds(i,ip)
     f(js,ip)=t(i)-f(js,ip)*cs(i,ip)/ds(i,ip)
  end do

  if(ipass.eq.1) then
     ipass=ipass+1
     js=js-1
     goto 2
  end if

!     solve the 2 eqns. f(1,1)*u0+f(1,2)*u1=s1; f(2,1)*u0+f(2,2)*u1=s2
!     put the sol'n in wrk
  
  if(ip.eq.2.and.ipass.eq.2) then

     wrk(1)=(s(1)*f(2,2)-s(2)*f(1,2))/(f(1,1)*f(2,2)-f(2,1)*f(1,2))
     wrk(2)=(s(2)*f(1,1)-s(1)*f(2,1))/(f(1,1)*f(2,2)-f(2,1)*f(1,2))

!     even u

     ne=n+1
     do i=2,ne
        ii=2*i -1
        wrk(ii)=(gs(i,1)-cs(i-1,1)*wrk(ii-2))/ds(i-1,1)
     end do

!     odd u

     do i=2,n
        ii=2*i
        wrk(ii)=(gs(i,2)-cs(i-1,2)*wrk(ii-2))/ds(i-1,2)
     end do

  end if
  
end subroutine p1dsolv
     

subroutine vcw3d(u,v,w,omx,omy,omz,fn,gn,           &
           u11,u12,u13,u21,u22,u23,u31,u32,u33,Lu,Lv,Lw)

    use omp_lib
    use grid_size
!**********************************************************************
!     this subroutine calculates the nonlinear term in the rotational *
!     form of the nse's. (3d-calculation)                             *
!                                                                     *
!      h1=(u x om)x
!                                                                     *
!      h2=(u x om)y
!
!      h3=(u x om)z 
!                                                                     *
!**********************************************************************     
!     the velocity and vorticity fields coming in are assumed to be 
!     in real space in y and fourier in x and z. 
!     before this code is called we must do a ytrans on six 
!     arrays containing the velocity and vorticity fields. 
!     also , we assume that the appropriate modes in all fields 
!     (i.e. first and last modes) are set to zero before we 
!     get to this point.
!      
!......................................................................

    !    parameters:      
    implicit none
    
    real xpg,mpg,xcpg,sigmapg ! for adverse pressure gradient
    real wavz(nz), wavx(nxh), c(nyp), yc(nyp), seght(nyp), trigx32(16000), trigz32(4000), spread(mz2,mx2)
    real fxintg(nyp,mz,mx), fyintg(nyp,mz,mx), fzintg(nyp,mz,mx)
    real pi, dt, theta, re, alpha, beta, xl, yl, zl, time, dyde, xstart, uinf
    real gain, ugain, delx, delz, delxm, delzm, segdrag, xsegdrag, ysegdrag, zsegdrag
  
    integer imatrix(nyp,mz,mx), ixfax32(19), izfax32(19)
    integer irstrt, nsteps, iprnfrq, i, j, k, idif, jdif, ipii, jpjj, inc, isign, lot, jump, ii, i1, i2, jj, it, kk
    integer print3d
    !    common blocks:
    common/data1/ pi, dt, theta, wavz, wavx, c, yc
    common/data2/ re, alpha, beta, xl, zl, yl, time, dyde, xstart, uinf
    common/iocntrl/ irstrt, nsteps, iprnfrq, print3d
        
    !    fft utilities arrays -- used in 3/2 rule calcs
    common/trigx32/ trigx32, ixfax32
    common/trigz32/ trigz32, izfax32
  
    common/preforce/ seght
    common/pre2forc/ spread,fxintg
    common/pre4f/ fzintg, fyintg
    common/pre5f/ it
    common/params/ gain, ugain
    common/geommat/ imatrix
        
    integer bfhead            ! index # for the buffer zone head face.
    integer bftail            ! index # for the buffer zone tail face.
    real bfgain(1200)          ! buffer zone stiffness distribution array.
    real bfugain(1200)         ! buffer zone damping coef. distribution array.
    integer bfwidth
    common /buffer/ bfgain, bfugain, bfhead, bftail, bfwidth
  
    real,dimension(nyp,nz,nx) :: initu, initv, initw
    real,dimension(nyp) :: ycrd
    common /init/ initu, initv, initw, ycrd
  
    integer perturbtime
    common /perturbation/ perturbtime
  
    integer kwall, kmaxsurf
    common /kwallpos/ kwall, kmaxsurf
  
    integer readvdes
    real vdes(mx), xsuction
    common /suction/ vdes, xsuction, readvdes

    integer particle_flag
    common /particle/ particle_flag


      
  !    passed variables:
    complex u(nyp,nz,nxh), v(nyp,nz,nxh), w(nyp,nz,nxh)
    complex omx(nyp,nz,nxh), omy(nyp,nz,nxh), omz(nyp,nz,nxh)
    complex fn(nyp,nz,nxh), gn(nyp,nz,nxh)
  
  !    local variables:
    real up(mzp,mxp2), vp(mzp,mxp2), wp(mzp,mxp2)
    real wx(mzp,mxp2), wy(mzp,mxp2), wz(mzp,mxp2)
    real vwx(mzp,mxp2), vwy(mzp,mxp2), vwz(mzp,mxp2)
    real wrk(nwrk)
  
    real u11p(mzp,mxp2),u12p(mzp,mxp2),u13p(mzp,mxp2)
    real u21p(mzp,mxp2),u22p(mzp,mxp2),u23p(mzp,mxp2) 
    real u31p(mzp,mxp2),u32p(mzp,mxp2),u33p(mzp,mxp2)
    real,dimension(mzp,mxp2) :: Lup,Lvp,Lwp ! Laplacian terms - Ryan 7/24/23
    
!      velocity gradient tensor 
    complex u11(nyp,nz,nxh),u12(nyp,nz,nxh),u13(nyp,nz,nxh)
    complex u21(nyp,nz,nxh),u22(nyp,nz,nxh),u23(nyp,nz,nxh) 
    complex u31(nyp,nz,nxh),u32(nyp,nz,nxh),u33(nyp,nz,nxh)
    complex,dimension(nyp,nz,nxh) :: Lu,Lv,Lw  ! Laplacian terms - Ryan 7/24/23

    real bdyfx, L, rad, xcenter, ycenter, zcenter,ampbdy1,ampbdy2,c11amp,tfrac,polyrate 
    common/vortRing/ bdyfx, L, rad, xcenter, ycenter, zcenter,ampbdy1,ampbdy2,c11amp,tfrac,polyrate  

    real dragx(mz,mx), dragy(mz,mx), dragz(mz,mx)
    real cfl(nyp)
  
    real cflcheck, cflmax, x, y, z

    ! Extra variables for vortex ring
    real xi, zj, argx, argrad, fx, fr, rsq
    real sumens, KE, argy
    real forbeta
    common/data22/ forbeta
  
    real vvv
  
    integer filenumber, filenumber3,filenumber2
      
    real wdes
    real slopelength  
      
    integer iiter, jiter, kiter   
    integer numxout, numyout, numzout
    
    ! Frozen turbulence variables
    real,save,dimension(nyp,mz,mx) :: u_old,v_old,w_old,up1,vp1,wp1,wx1,wy1,wz1
    real,save,dimension(nyp,mz,mx) :: u11p1,u12p1,u13p1,u21p1,u22p1,u23p1,u31p1,u32p1,u33p1,Lup1,Lvp1,Lwp1 ! Laplacian terms - Ryan 7/24/23
    logical :: frozen_flag
    common /frozen/ frozen_flag

    real,dimension(nyp,mz,mx) :: up3d,vp3d,wp3d,wx3d,wy3d,wz3d,scp3d,beta3d
    real,dimension(nyp,mz,mx) :: u11p3d,u12p3d,u13p3d,u21p3d,u22p3d,u23p3d,u31p3d,u32p3d,u33p3d,Lup3d,Lvp3d,Lwp3d ! Laplacian terms - Ryan 7/24/23
    real, dimension(mz) :: uxmean
    real :: dudx, dvdx, dwdx, dudy, dvdy, dwdy, dudz, dvdz, dwdz
    real,dimension(nyp) :: uzmean
    real dely,massFlux
    real swirl
    real, dimension(nyp,mz,mx) :: lamb2_3d, swirl_3d
    integer num_threads 
    integer maxrec

    real :: dPdx, R_tau
    common /pressure/ dPdx, R_tau
   
    real :: xp, yp, zp, u_interp, ureal, perr
  
    integer geomtype, flow_select
    common /geometry/ geomtype, flow_select


    ! Gain smoothing
    real dyc, dyhalf, dxc, dxhalf, dzc, dzhalf, dugain, digain
  
  
    ! Output grid data
    integer imin, imax, jmin, jmax, kmin, kmax, istep, kstep
  
  
  !++++++++++++++++++  
  
    bfhead  = 2
    bftail  = bftail_             
    bfwidth = (bftail - bfhead)
    slopelength =bfwidth/2
  !++++++++++++++++++  
  
    filenumber = it/iprnfrq + 200
    filenumber2 = it/iprnfrq + 1200
    
    pi = 2.0 * acos(0.0)
    delx = xl / float(nx-1)
    delz = zl / float(nz-1)
    delxm = xl /float(mx-1)
    delzm = zl /float(mz-1)
    cfl = 0.0
    cflmax = 0.0
    
  !++++++++++++++++++  
    kiter=0
  
    
   
    
   !    calculate the v cross omega term plane by plane
  
  numxout=mxp
  numzout=mzp
  numyout=ny/4   

!-------------------------------------------------------------------------------!
  !$omp parallel do default(shared) private(k,j,i,jj,i1,i2,up,vp,wp,     &
  !$omp   wx,wy,wz,vwx,vwy,vwz,dragx,dragy,dragz,wrk,inc,isign,          &
  !$omp   u11p,u12p,u13p,u21p,u22p,u23p,u31p,u32p,u33p,Lup,Lvp,Lwp,      &
  !$omp   jump,lot,ii,ipii,idif,jdif,segdrag,jpjj,                       &
  !$omp   xsegdrag,ysegdrag,zsegdrag,cflcheck,y,x,z,iiter,jiter) schedule(dynamic)

!-------------------------------------------------------------------------------!
  
  do k = 1, nyp
      iiter=0
      jiter=0
  !    initialize the velocities and vorticities on the 
  !    k-th x-z plane with zeros, we need to pad the
  !    extra 1/2 modes with zeros by 3/2 rule.   
           
       do j = 1, mzp
          do i = 1, mxp2
  
             up(j,i) = 0.0
             vp(j,i) = 0.0
             wp(j,i) = 0.0
  
             wx(j,i) = 0.0
             wy(j,i) = 0.0
             wz(j,i) = 0.0
  
             vwx(j,i) = 0.0
             vwy(j,i) = 0.0
             vwz(j,i) = 0.0

             ! use in particle tracking
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
          end do
       end do
       
       do j = 1,nz
          if(j.le.nzh) jj = j
          if(j.gt.nzh) jj=(mz-nz) +j
  
          do i = 1,nxh
             i1 = 2*(i-1) +1
  
             up(jj,i1) = real(u(k,j,i))
             vp(jj,i1) = real(v(k,j,i))
             wp(jj,i1) = real(w(k,j,i))
  
             wx(jj,i1) = real(omx(k,j,i))
             wy(jj,i1) = real(omy(k,j,i))
             wz(jj,i1) = real(omz(k,j,i))
 
             u11p(jj,i1) = real(u11(k,j,i))
             u12p(jj,i1) = real(u12(k,j,i))
             u13p(jj,i1) = real(u13(k,j,i))
             u21p(jj,i1) = real(u21(k,j,i))
             u22p(jj,i1) = real(u22(k,j,i))
             u23p(jj,i1) = real(u23(k,j,i))
             u31p(jj,i1) = real(u31(k,j,i))
             u32p(jj,i1) = real(u32(k,j,i))
             u33p(jj,i1) = real(u33(k,j,i))
            
             Lup(jj,i1) = real(Lu(k,j,i))
             Lvp(jj,i1) = real(Lv(k,j,i))
             Lwp(jj,i1) = real(Lw(k,j,i))
          end do
          do i = 1,nxh
             i2 = 2*i
  
             up(jj,i2) = aimag(u(k,j,i))
             vp(jj,i2) = aimag(v(k,j,i))
             wp(jj,i2) = aimag(w(k,j,i))
  
             wx(jj,i2) = aimag(omx(k,j,i))
             wy(jj,i2) = aimag(omy(k,j,i))
             wz(jj,i2) = aimag(omz(k,j,i))
  
             u11p(jj,i2) = aimag(u11(k,j,i))
             u12p(jj,i2) = aimag(u12(k,j,i))
             u13p(jj,i2) = aimag(u13(k,j,i))
             u21p(jj,i2) = aimag(u21(k,j,i))
             u22p(jj,i2) = aimag(u22(k,j,i))
             u23p(jj,i2) = aimag(u23(k,j,i))
             u31p(jj,i2) = aimag(u31(k,j,i))
             u32p(jj,i2) = aimag(u32(k,j,i))
             u33p(jj,i2) = aimag(u33(k,j,i))
           
             Lup(jj,i2) = aimag(Lu(k,j,i))
             Lvp(jj,i2) = aimag(Lv(k,j,i))
             Lwp(jj,i2) = aimag(Lw(k,j,i))
          end do
       end do
  
  
  !......................................................................
  !
  !    transform (interpolate) into 3/2 grid physical space.
           
       inc   = 1
       isign = 1
       jump  = 2*mzp
       lot   = nx/2
  
       call cfftmlt(up(1,1),up(1,2),wrk,trigz32,izfax32,inc,jump,mz,lot,isign)
       call cfftmlt(vp(1,1),vp(1,2),wrk,trigz32,izfax32,inc,jump,mz,lot,isign)
       call cfftmlt(wp(1,1),wp(1,2),wrk,trigz32,izfax32,inc,jump,mz,lot,isign)
       call cfftmlt(wx(1,1),wx(1,2),wrk,trigz32,izfax32,inc,jump,mz,lot,isign)
       call cfftmlt(wy(1,1),wy(1,2),wrk,trigz32,izfax32,inc,jump,mz,lot,isign)
       call cfftmlt(wz(1,1),wz(1,2),wrk,trigz32,izfax32,inc,jump,mz,lot,isign)
! New variables - Ryan 2-24-22
       call cfftmlt(u11p(1,1),u11p(1,2),wrk,trigz32,izfax32,inc,jump,mz,lot,isign)
       call cfftmlt(u12p(1,1),u12p(1,2),wrk,trigz32,izfax32,inc,jump,mz,lot,isign)
       call cfftmlt(u13p(1,1),u13p(1,2),wrk,trigz32,izfax32,inc,jump,mz,lot,isign)
       call cfftmlt(u21p(1,1),u21p(1,2),wrk,trigz32,izfax32,inc,jump,mz,lot,isign)
       call cfftmlt(u22p(1,1),u22p(1,2),wrk,trigz32,izfax32,inc,jump,mz,lot,isign)
       call cfftmlt(u23p(1,1),u23p(1,2),wrk,trigz32,izfax32,inc,jump,mz,lot,isign)
       call cfftmlt(u31p(1,1),u31p(1,2),wrk,trigz32,izfax32,inc,jump,mz,lot,isign)
       call cfftmlt(u32p(1,1),u32p(1,2),wrk,trigz32,izfax32,inc,jump,mz,lot,isign)
       call cfftmlt(u33p(1,1),u33p(1,2),wrk,trigz32,izfax32,inc,jump,mz,lot,isign)
     
       call cfftmlt(Lup(1,1),Lup(1,2),wrk,trigz32,izfax32,inc,jump,mz,lot,isign)
       call cfftmlt(Lvp(1,1),Lvp(1,2),wrk,trigz32,izfax32,inc,jump,mz,lot,isign)
       call cfftmlt(Lwp(1,1),Lwp(1,2),wrk,trigz32,izfax32,inc,jump,mz,lot,isign)

       do j = 1, mz
          up(j,nxp) = up(j,2)
          vp(j,nxp) = vp(j,2)
          wp(j,nxp) = wp(j,2)
          wx(j,nxp) = wx(j,2)
          wy(j,nxp) = wy(j,2)
          wz(j,nxp) = wz(j,2)
       end do
  
       do j = 1, mz
          up(j,nxp2) = 0.0 
          vp(j,nxp2) = 0.0
          wp(j,nxp2) = 0.0
          wx(j,nxp2) = 0.0
          wy(j,nxp2) = 0.0
          wz(j,nxp2) = 0.0
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
       end do
  
       isign = 1
       inc   = mzp
       jump  = 1
       lot   = mz
  
       call rfftmlt(up,wrk,trigx32,ixfax32,inc,jump,mx,lot,isign)
       call rfftmlt(vp,wrk,trigx32,ixfax32,inc,jump,mx,lot,isign)
       call rfftmlt(wp,wrk,trigx32,ixfax32,inc,jump,mx,lot,isign)
       call rfftmlt(wx,wrk,trigx32,ixfax32,inc,jump,mx,lot,isign)
       call rfftmlt(wy,wrk,trigx32,ixfax32,inc,jump,mx,lot,isign)
       call rfftmlt(wz,wrk,trigx32,ixfax32,inc,jump,mx,lot,isign)

       call rfftmlt(u11p,wrk,trigx32,ixfax32,inc,jump,mx,lot,isign)
       call rfftmlt(u12p,wrk,trigx32,ixfax32,inc,jump,mx,lot,isign)
       call rfftmlt(u13p,wrk,trigx32,ixfax32,inc,jump,mx,lot,isign)
       call rfftmlt(u21p,wrk,trigx32,ixfax32,inc,jump,mx,lot,isign)
       call rfftmlt(u22p,wrk,trigx32,ixfax32,inc,jump,mx,lot,isign)
       call rfftmlt(u23p,wrk,trigx32,ixfax32,inc,jump,mx,lot,isign)
       call rfftmlt(u31p,wrk,trigx32,ixfax32,inc,jump,mx,lot,isign)
       call rfftmlt(u32p,wrk,trigx32,ixfax32,inc,jump,mx,lot,isign)
       call rfftmlt(u33p,wrk,trigx32,ixfax32,inc,jump,mx,lot,isign)
     
       call rfftmlt(Lup,wrk,trigx32,ixfax32,inc,jump,mx,lot,isign)
       call rfftmlt(Lvp,wrk,trigx32,ixfax32,inc,jump,mx,lot,isign)
       call rfftmlt(Lwp,wrk,trigx32,ixfax32,inc,jump,mx,lot,isign)
  !......................................................................
  
  !    compute the cross product of velocity and vorticity.
  
       do j = 1, mz
          do i = 1, mx
             vwx(j,i) =  vp(j,i) * wz(j,i) - wp(j,i) * wy(j,i)
             vwy(j,i) = -up(j,i) * wz(j,i) + wp(j,i) * wx(j,i)
             vwz(j,i) =  up(j,i) * wy(j,i) - vp(j,i) * wx(j,i)

             cflcheck = (abs(up(j,i))/delx + abs(vp(j,i))/seght(k) + abs(wp(j,i))/delz) * dt
             
             if(cflcheck .gt. cfl(k)) cfl(k) = cflcheck
          end do
       end do
  
  !......................................................................
  
  if (flow_select .ge. 0) then    
  !    initialize the force field array.
  
       dragx = 0.0
       dragy = 0.0
       dragz = 0.0      
  
       do j = 1,mz
          do i = 1,mx 
  
  !    apply the forces to generate solid surfaces.
             if (imatrix(k,j,i) .eq. 1) then 
                 fxintg(k,j,i) = fxintg(k,j,i) + up(j,i)
                 fyintg(k,j,i) = fyintg(k,j,i) + vp(j,i)
                 fzintg(k,j,i) = fzintg(k,j,i) + wp(j,i)
                 dragx(j,i) = (-ugain*up(j,i) - gain*fxintg(k,j,i))
                 dragy(j,i) = (-ugain*vp(j,i) - gain*fyintg(k,j,i))
                 dragz(j,i) = (-ugain*wp(j,i) - gain*fzintg(k,j,i))
  
  !    apply the forces to generate spanwise-damping textures.
            else if (imatrix(k,j,i) .eq. 3) then 
                fzintg(k,j,i) = fzintg(k,j,i) + wp(j,i)
                dragz(j,i) = (-ugain*wp(j,i) - gain*fzintg(k,j,i))
                
  !    apply the forces to generate spanwise moving wall
            else if (imatrix(k,j,i) .eq. 8) then            
                if (i .le. bfwidth*3+1) then
                   wdes = (0.3*uinf)/slopelength*(i-(bfwidth*3+1-slopelength))
                else
                   wdes = 0.3*uinf
                end if
  
                fxintg(k,j,i) = fxintg(k,j,i) + up(j,i)
                fyintg(k,j,i) = fyintg(k,j,i) + vp(j,i)
                fzintg(k,j,i) = fzintg(k,j,i) + wp(j,i)-wdes
                dragx(j,i)    = (-ugain*up(j,i) - gain*fxintg(k,j,i))
                dragy(j,i)    = (-ugain*vp(j,i) - gain*fyintg(k,j,i))
                dragz(j,i)    = (-ugain*(wp(j,i)-wdes) - gain*fzintg(k,j,i))
                
  !    apply the buffer zone forces.
            else if (imatrix(k,j,i) .eq. 6) then
                fxintg(k,j,i) = fxintg(k,j,i) + (up(j,i)-initu(k,j,i)) ! update buffer zone integrals. 
                fyintg(k,j,i) = fyintg(k,j,i) + (vp(j,i)-initv(k,j,i))
                fzintg(k,j,i) = fzintg(k,j,i) + (wp(j,i)-initw(k,j,i))
                vwx(j,i) = vwx(j,i) + (-bfugain(i-bfhead+1) * (up(j,i) - initu(k,j,i)) - (bfgain(i-bfhead+1) * fxintg(k,j,i)))
                vwy(j,i) = vwy(j,i) + (-bfugain(i-bfhead+1) * (vp(j,i) - initv(k,j,i)) - (bfgain(i-bfhead+1) * fyintg(k,j,i)))
                vwz(j,i) = vwz(j,i) + (-bfugain(i-bfhead+1) * (wp(j,i) - initw(k,j,i)) - (bfgain(i-bfhead+1) * fzintg(k,j,i)))

  !    apply the alternative buffer zone forces (for flow at an angle).
            else if (imatrix(k,j,i) .eq. 26) then
                vwx(j,i) = vwx(j,i) - ugain * (up(j,i) - uinf * cos(10. * pi / 180.))
                vwy(j,i) = vwy(j,i) - ugain *  vp(j,i)
                vwz(j,i) = vwz(j,i) - ugain * (wp(j,i) - uinf * sin(10. * pi / 180.))

            else if (imatrix(k,j,i) .eq. 120) then
                vwx(j,i) = vwx(j,i) - ugain * (up(j,i) - uinf * cos(20. * pi / 180.))
                vwy(j,i) = vwy(j,i) - ugain *  vp(j,i)
                vwz(j,i) = vwz(j,i) - ugain * (wp(j,i) - uinf * sin(20. * pi / 180.))

  !    apply the suction force for the vertical velocity over the boundary layer.
            else if (imatrix(k,j,i) .eq. 7) then
                fyintg(k,j,i) = fyintg(k,j,i) + (vp(j,i)-vdes(i))
                vwy(j,i) = vwy(j,i) + (-gain*fyintg(k,j,i)) + (-ugain*(vp(j,i)-vdes(i)))
                fxintg(k,j,i) = fxintg(k,j,i) + (up(j,i)-uinf)
                vwx(j,i) = vwx(j,i) + (-gain*fxintg(k,j,i)) + (-ugain*(up(j,i)-uinf))
                fzintg(k,j,i) = fzintg(k,j,i) + (wp(j,i)-initw(k,j,i))
                vwz(j,i) = vwz(j,i) + (-gain*fzintg(k,j,i)) + (-ugain*(wp(j,i)-initw(k,j,i)))
  
  !    apply a constant x-force in a region between the wall and the suction layer.
            else if (imatrix(k,j,i) .eq. 11) then
                dragx(j,i) = 1.0
  
  !    apply the forces to generate the perturbation which triggers the turbulent spot.
            else if (imatrix(k,j,i) .eq. 4) then
                if (it .le. perturbtime) then
                   fxintg(k,j,i) = fxintg(k,j,i) + up(j,i)
                   fyintg(k,j,i) = fyintg(k,j,i) + vp(j,i)
                   fzintg(k,j,i) = fzintg(k,j,i) + wp(j,i)
                   dragx(j,i) = (-ugain*up(j,i) - gain*fxintg(k,j,i))
                   dragy(j,i) = (-ugain*vp(j,i) - gain*fyintg(k,j,i))
                   dragz(j,i) = (-ugain*wp(j,i) - gain*fzintg(k,j,i))
                end if

  !    apply a force to simulate a pressure gradient - Ryan 7/8/21
            else if (imatrix(k,j,i) .eq. 0 .and. flow_select .eq. 1) then
                   vwx(j,i) = vwx(j,i) - dPdx

  !    apply a force to generate a vortex ring - Ryan 7/11/22
            else if (flow_select .eq. 5 .and. time .lt. tfrac*nsteps*dt) then
                zj = float(j-1)*delzm
                rsq = (ycoord(k) - ycenter)**2 + (zj - zcenter)**2
                argrad = forbeta*(rad - sqrt(rsq))
                fr = 0.5*(1.0 + tanh(argrad))

                xi = float(i-1)*delxm
                argx = forbeta*(L - abs(xi - xcenter))
                fx = 0.5*(1.0 + tanh(argx))

                vwx(j,i) = vwx(j,i) + bdyfx*fx*fr

            end if
           
          end do
      end do 
    end if

  !......................................................................

         
  !    apply the force field to the solid surface
       if(k .le. kmaxsurf) then
          do i = 1,mx
             do ii = -2, 2
                ipii = i + ii
                idif = 1 + iabs(ii)        
                if(ipii .lt. 1 ) ipii = ipii + mx         
                if(ipii .gt. mx) ipii = ipii - mx
                do jj = -2, 2
                   jdif = 1 + iabs(jj)
                   segdrag = spread(jdif,idif)
                   do j = 3, mzm2
                      jpjj = j + jj
                      xsegdrag = segdrag * dragx(j,i)
                      ysegdrag = segdrag * dragy(j,i)
                      zsegdrag = segdrag * dragz(j,i)
                      vwx(jpjj,ipii) = vwx(jpjj,ipii) + xsegdrag
                      vwy(jpjj,ipii) = vwy(jpjj,ipii) + ysegdrag
                      vwz(jpjj,ipii) = vwz(jpjj,ipii) + zsegdrag             
                   end do
                   do j = 1, 2
                      jpjj = j + jj
                      if(jpjj .lt. 1 ) jpjj = jpjj + mz
                      xsegdrag = segdrag * dragx(j,i)
                      ysegdrag = segdrag * dragy(j,i)
                      zsegdrag = segdrag * dragz(j,i)
                      vwx(jpjj,ipii) = vwx(jpjj,ipii) + xsegdrag
                      vwy(jpjj,ipii) = vwy(jpjj,ipii) + ysegdrag
                      vwz(jpjj,ipii) = vwz(jpjj,ipii) + zsegdrag
                   end do
                   do j = mzm1, mz
                      jpjj = j + jj
                      if(jpjj .gt. mz) jpjj = jpjj - mz
                      xsegdrag = segdrag * dragx(j,i)
                      ysegdrag = segdrag * dragy(j,i)
                      zsegdrag = segdrag * dragz(j,i)
                      vwx(jpjj,ipii) = vwx(jpjj,ipii) + xsegdrag
                      vwy(jpjj,ipii) = vwy(jpjj,ipii) + ysegdrag
                      vwz(jpjj,ipii) = vwz(jpjj,ipii) + zsegdrag
                   end do
                end do
             end do
          end do
       end if
       
  
       ! Saving variables to print  
       ! Flip the sign of wx and wz due to the upside down nature wall normal direction
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
       Lvp3d(k,1:mz,1:mx) = Lvp(1:mz,1:mx)  ! I'm not sure if I need to flip signs of any of these...
       Lwp3d(k,1:mz,1:mx) = Lwp(1:mz,1:mx) 
  !......................................................................
  !.... now transform vxw back to spectral space.
       
       isign = -1
       inc=mzp
       jump = 1
       lot = mz
  
       call rfftmlt(vwx,wrk,trigx32,ixfax32,inc,jump,mx,lot,isign)
       call rfftmlt(vwy,wrk,trigx32,ixfax32,inc,jump,mx,lot,isign)
       call rfftmlt(vwz,wrk,trigx32,ixfax32,inc,jump,mx,lot,isign)

  !now wiping out unwanted last mode.
       do j=1,mz
          vwx(j,2) = 0.0          
          vwy(j,2) = 0.0          
          vwz(j,2) = 0.0
       end do
  
       inc = 1
       isign = -1
       jump = 2*mzp
       lot = nx/2
  
       call cfftmlt(vwx(1,1),vwx(1,2),wrk,trigz32,izfax32,inc,jump,mz,lot,isign)
       call cfftmlt(vwy(1,1),vwy(1,2),wrk,trigz32,izfax32,inc,jump,mz,lot,isign)
       call cfftmlt(vwz(1,1),vwz(1,2),wrk,trigz32,izfax32,inc,jump,mz,lot,isign)

       do j = 1, nz
          if(j.le.nzh) then
             jj = j
          else if(j.gt.nzh) then
             jj=(mz-nz) + j
          end if
          do i = 1, nxh
             i1 = 2*(i-1) +1
             i2 = 2*i
             gn(k,j,i)  = cmplx(vwx(jj,i1)*rmz, vwx(jj,i2)*rmz)
             fn(k,j,i)  = cmplx(vwy(jj,i1)*rmz, vwy(jj,i2)*rmz)
             omz(k,j,i) = cmplx(vwz(jj,i1)*rmz, vwz(jj,i2)*rmz)
          end do
       end do
         
  end do
  !+++++++++++++++++++++++++++
  !$omp end parallel do
  
  !                            -----  end of loop over normal(y) planes  -----
  
  
!*** Removed lambda2 calculations in favor of swirl calculations (I think it looks nicer) - Ryan 3/10/23 ***!

      ! Calculate swirl and print flowfield 
      if (print3d .ne. 0) then 


        if((mod(it,iprnfrq) .eq. 0 .and. it .ne. 0) .or. it .eq. 1) then
            ! Calculate swirl
            write(*,*) 'Start calculating swirl'
            !$omp parallel do private(i,j,k,dudx,dvdx,dwdx,dudy,dvdy,dwdy,dudz,dvdz,dwdz,swirl)     
                do k = 1,nyp
                    do j = 1,mz
                        do i = 1,mx
                            call calcswirl(u11p3d(k,j,i),u21p3d(k,j,i),u31p3d(k,j,i),u12p3d(k,j,i),u22p3d(k,j,i), &
                                           u32p3d(k,j,i),u13p3d(k,j,i),u23p3d(k,j,i),u33p3d(k,j,i),swirl)

                            swirl_3d(k,j,i) = swirl
                        end do
                    end do
                end do
            !$omp end parallel do
            write(*,*) 'Finished calculating swirl'
            ! End swirl calculation
    
            ! Write output files
            if (print3d .eq. 1 .or. print3d .eq. 2) then
                call write_flowfield_ascii(it, up3d, vp3d, wp3d, wx3d, wy3d, wz3d, swirl_3d, real(imatrix), xl, zl, print3d) ! Changed scp3d to beta3d - Ryan 9/23/22
            else 
                call write_flowfield_plt(it, up3d, vp3d, wp3d, wx3d, wy3d, wz3d, swirl_3d, real(imatrix), xl, zl)
            end if
        end if
    end if
    

    !===================================================================================!
    !                              Interpolation Check                                  !    
    !===================================================================================!
!    if (it .eq. 1) then
!        open(70, file = "outputs/u_interp.dat")
! 
!        do i = 1,ny
!            xp = float(mx/2)*delxm  - delxm/2 ! Arbitrary point between x grid nodes
!            yp = (ycoord(i) + ycoord(i+1))/2.0
!            zp = float(mz/2)*delzm - delzm/2 ! Arbitrary point betwen z grid nodes
!            ureal = yp
!            !u_interp = (up3d(i,mx/2,mz/2) + up3d(i+1,mx/2,mz/2))/2.0
!            call spec(u_interp,up3d,xp,yp,zp)
!            
!            perr = abs((u_interp - ureal)/ureal)*100.0
!            write(70,"(*(e14.6,1x))") yp, ureal, u_interp, perr
!        end do
!        close(70)
!    end if    

!    if (it .eq. 2) then
!        write(797,*) u_old
!        write(798,*) v_old
!        write(799,*) w_old
!        write(897,*) wx3d
!        write(898,*) wy3d
!        write(899,*) wz3d
!        write(997,*) up3d
!        write(998,*) vp3d
!        write(999,*) wp3d
!    end if

    !-----------------------------------------------------------------------------------!


    !===================================================================================!
    !                           Particle Tracking Implentation                          !    
    !===================================================================================!

    ! Just change this flag to turn frozen flow for the particles on or off
    frozen_flag = .false.

    if (npart .ne. 0) then
       if (it .le. 2 .or. frozen_flag .eq. .false.) then  ! Frozen turbulence - Ryan 3-10-22
          call part_track(up3d,vp3d,wp3d,wx3d,wy3d,wz3d,u_old,v_old,w_old,  &
                          u11p3d,u12p3d,u13p3d,u21p3d,u22p3d,u23p3d,u31p3d, &
                          u32p3d,u33p3d,Lup3d,Lvp3d,Lwp3d,Lup1,Lvp1,Lwp1)
       else
           call part_track(up1,vp1,wp1,wx1,wy1,wz1,u_old,v_old,w_old, &
                          u11p1,u12p1,u13p1,u21p1,u22p1,u23p1,u31p1,u32p1,u33p1)
       end if
    end if

    if (it .eq. 1 .or. frozen_flag .eq. .false.) then
       u_old = up3d
       v_old = vp3d
       w_old = wp3d

       Lup1 = Lup3d
       Lvp1 = Lvp3d
       Lwp1 = Lwp3d
    elseif (it .eq. 2) then
       up1 = up3d
       vp1 = vp3d
       wp1 = wp3d
       wx1 = wx3d
       wy1 = wy3d
       wz1 = wz3d

       u11p1 = u11p3d
       u12p1 = u12p3d
       u13p1 = u13p3d
       u21p1 = u21p3d
       u22p1 = u22p3d
       u23p1 = u23p3d
       u31p1 = u31p3d
       u32p1 = u32p3d
       u33p1 = u33p3d
    end if 
    !-----------------------------------------------------------------------------------!


    !===================================================================================!
    !                               Writing Mean U Data                                 !    
    !===================================================================================!

    if (flow_select .eq. 1) then ! This is only relevant for channel flow
    print *,'Writing mean U data...'
        ! Mean velocity
        if (it .eq. 1) then
            open(71, file = "outputs/mean_u_data.dat")
        else
            open(71, file = "outputs/mean_u_data.dat", position = "append")
        end if
    
        do k = 1,nyp
            do j = 1,mz
                uxmean(j) = sum(up3d(k,j,:))/mx
            end do
            uzmean(k) = sum(uxmean)/mz
            write(71,"(*(e14.6,1x))") uzmean ! up3d(k,floor(mz/2.0),floor(mx/2.0))
        end do
        close(71)
    print *,'Finished writing mean U data'
    end if

    ! Computing Mass Flux: integrate <U> over y - Ryan 8/2/23
    if (irstrt .eq. it) then
        open(73,file="outputs/mass_flux")
        open(74,file="outputs/C_f")
    else
        open(73,file="outputs/mass_flux",position="append")
        open(74,file="outputs/C_f",position="append")
    end if

    massFlux = 0.0
    do i = 1,ny 
        massFlux = massFlux + 1.0/yl*0.5*(uzmean(i+1) + uzmean(i))*(ycoord(i+1) - ycoord(i)) ! bulk velocity
    end do

    write(73,*) massFlux
    close(73)

    !-----------------------------------------------------------------------------------!
   

    ! Computing Average KE and Enstrophy for VR
    if (flow_select .eq. -1 .or. flow_select .eq. 5 .and. print3d .ge. 0) then
        sumens = 0.0
        KE = 0.0

        do i = 1,nyp
          argy = pi*float(i-1)/float(ny)
          do j = 1,mz
            do k = 1,mx
              sumens = sumens + (wx3d(i,j,k)**2 + wy3d(i,j,k)**2 + wz3d(i,j,k)**2)*sin(argy)
              KE = KE + 0.5*(up3d(i,j,k)**2 + vp3d(i,j,k)**2 + wp3d(i,j,i)**2)
            end do
          end do
        end do

        if (it .eq. 0) then
            call system('rm outputs/enstrophy')
            call system('rm outputs/KE')
        end if
        open(45,file='outputs/enstrophy',position="append")
        write(45,*) sumens
        close(45)

        open(46,file='outputs/KE',position="append")
        write(46,*) KE
        close(46)
    end if
                     
 
      do k = 1,nyp
         if (cfl(k) > cflmax) cflmax = cfl(k)
      end do
    
      write(*,*) 'cflmax = ', cflmax
          
    !    check cfl number      
      if(cflmax .gt. 1) then
         write(*,*) 'cfl failure, time step', it
         stop
      end if
   
end subroutine vcw3d

!----------------------------------------------------------

!subroutine calclambda2(dudx, dvdx, dwdx, dudy, dvdy, dwdy, dudz, dvdz, dwdz, lambda2)
!
!    real :: dudx, dvdx, dwdx, dudy, dvdy, dwdy, dudz, dvdz, dwdz, lambda2
!    real :: s11, s12, s13, s22, s23, s33
!    real :: omga12, omga13, omga23
!    character (len=1) :: uplo,jobz
!    integer :: ia, info, iz, ja, jz, lwork, n
!    real :: s2o2(3,3)
!    real :: w(3),work(9)
!    
!    
!    uplo = 'u' !upper or lower triangular
!    n = 3 !order of matrix
!    lda = 3
!    ia = 1
!    lwork = 9
!    ja = 1
!    jobz = 'n' !calculate eigenvalue only
!    
!    s11 = dudx
!    s12 = 0.5*(dudy+dvdx)
!    s13 = 0.5*(dudz+dwdx)
!    s22 = dvdy
!    s23 = 0.5*(dvdz+dwdy)
!    s33 = dwdz
!        
!    omga12 = 0.5*(dudy-dvdx)
!    omga13 = 0.5*(dudz-dwdx)
!    omga23 = 0.5*(dvdz-dwdy)
!    
!    s2o2(1,1) = s11**2 + s12**2 + s13**2 - omga12**2 - omga13**2
!    s2o2(1,2) = s11*s12 + s12*s22 + s13*s23 - omga13*omga23
!    s2o2(1,3) = s11*s13 + s12*s23 + s13*s33 -omga12*omga23
!    s2o2(2,2) = s12**2 + s22**2 + s23**2 - omga12**2 - omga23**2
!    s2o2(2,3) = s12*s13 + s22*s23 + s23*s33 - omga12*omga13
!    s2o2(3,3) = s13**2 + s23**2 + s33**2 - omga13**2 - omga**2
!    s2o2(2,1) = s2o2(1,2)
!    s2o2(3,1) = s2o2(1,3)
!    s2o2(3,2) = s2o2(2,3)
!
!    call dsyev(jobz,uplo,n,s2o2,lda,w,work,lwork,info)
!
!    lambda2=w(2)
!
!end subroutine

! Updating flowfield writing 7/8/22
subroutine write_flowfield_ascii(it, u, v, w, wx, wy, wz, swirl, ctp, xl, zl, print3d)

    use grid_size
    integer :: it, print3d, mz_copy,mx_copy
    integer :: irstrt, nsteps,iprnfrq
    common/iocntrl/ irstrt,nsteps,iprnfrq
    real, dimension(nyp,mz,mx) :: u, v, w, wx, wy, wz, ctp, swirl
    character (33) :: filename

    real :: x, y, z, delxm, delzm
    integer :: i, j, k

    real xl, yl, zl

    ! Printing in 2D workaround
    if (nz .gt. 10) then
       mz_copy = mz
    else
        mz_copy = 1
    end if
    if (nx .gt. 10) then
       mx_copy = mx
    else
        mx_copy = 1
    end if

    delxm = xl /float(mx-1)
    delzm = zl /float(mz-1)

    print *, 'Writing data at time step ', it
    
    write(filename,'("outputs/flowfield/time-",i6.6,".dat")')it
    
    ! Write 3d-data 
    open(6, file = filename)
    write(6,9711) 'filetype = solution, variables = "swirl", "u", "v", "w", "wx", "wy", "wz", "ctp"'
    write(6,9712) 'zone f=point t="Field", solutiontime=', it/iprnfrq,',i=',mx, 'j=',1, 'k=', nyp, new_line('a')
    
    do k = 1,nyp
        y = ycoord(k)
        do j = 1,mz_copy
            z = float(j-1)*delzm 
            do i = 1,mx
                x = float(i-1)*delxm
                write(6,8611) swirl(k,j,i),u(k,j,i),v(k,j,i),w(k,j,i),wx(k,j,i),wy(k,j,i),wz(k,j,i),ctp(k,j,i)
             end do
        end do
    end do 
    
    write(*,*) 'Finished writing flowfield'
    close(6)

    8611 format(8(e14.6,1x))
    9711 format(a93)
    9712 format(a40,i3,a5,i4,a5,i4,a5,i4,a2)
   
    if (print3d .eq. 2) then 
        ! Convert ascii to tecplot binary
        call system('preplot-bin ' // filename // ' && rm ' // filename)
    end if
    

    ! Write grid file if it == 1
    if (it .eq. 1) then
        open(5, file = "outputs/flowfield/grid.dat", status = "replace")
    
        write(5,*) 'filetype = grid, variables = "x", "y", "z"'
        write(5,*) 'zone f=point t="Grid", i=', mx, ' j=', 1, 'k=', nyp
        do k = 1,nyp
            y = ycoord(k)
            do j = 1,mz_copy
                z = float(j-1)*delzm 
                do i = 1,mx
                    x = float(i-1)*delxm
                    write(5,9400) x,y,z
                 end do
            end do
        end do 
    
        close(5)
    
        9400 format(12(e14.6,1x))  

        if (print3d .eq. 2) then
            ! Convert ascii to tecplot binary
            call system('preplot-bin outputs/flowfield/grid.dat && rm outputs/flowfield/grid.dat')
        end if

    end if
end subroutine


#ifdef OUTPUTFORM

subroutine write_flowfield_plt(it, u, v, w, wx, wy, wz, swirl, ctp, xl, zl)
    use iso_c_binding 
    use grid_size
    implicit none

    include 'tecio.f90'

    integer :: it, mz_copy,mx_copy
    integer :: irstrt, nsteps,iprnfrq
    common/iocntrl/ irstrt,nsteps,iprnfrq
    real, dimension(nyp,mz,mx) :: u, v, w, wx, wy, wz, ctp
    real, dimension(nyp,mz,mx) :: swirl
    character*1 :: NULLCHR
    character (35) :: filename
    

    real :: x, y, z, delxm, delzm
    integer :: i, j, k, ii, jj, kk

    real xl, yl, zl

    integer(c_int64_t) :: numValues
    integer(c_int32_t) :: fileFormat = 1 ! 0 = .plt | 1 = .szplt
    integer(c_int32_t) :: gridFileType = 1 ! GRID
    integer(c_int32_t) :: outputFileType = 2 ! SOLUTION
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

    ! Printing in 2D workaround
    if (nz .gt. 10) then
       mz_copy = mz
    else
        mz_copy = 1
    end if

    if (nx .gt. 10) then
       mx_copy = mx
    else
        mx_copy = 1
    end if
    ! Initialize helper arrays
    numValues = (nyp) * (mz_copy) * (mx_copy)
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
    if (it .eq. 1) then
        ! Open files
        i = tecini142("Grid", "x y z", "outputs/flowfield/grid.szplt"//NULLCHR, ".", &
            fileFormat, gridFileType, 1, 0)

        ! Create zone
        i = teczne142( &
            "Flowfield", &
            zoneType, &
            (nyp), &
            (mz_copy), &
            (mx_copy), &
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

        ! Compute Coordinates
        allocate(coordinates(nyp,mz_copy,mx_copy, 3))
        delxm = xl /float(mx-1)
        delzm = zl /float(mz-1)
        do k = 1,nyp
            y = ycoord(k)
            do j = 1,mz_copy
                z = float(j-1) * delzm 
                do i = 1,mx_copy
                    x = float(i-1) * delxm
                    coordinates(k,j,i,1) = x
                    coordinates(k,j,i,2) = y
                    coordinates(k,j,i,3) = z
                end do
            end do
        end do 

        ! Write x
        floatValues = pack(coordinates(:,:,:,1), .true.)
        i = tecdat142(numValues, floatValues, isDouble)

        ! Write y
        floatValues = pack(coordinates(:,:,:,2), .true.)
        i = tecdat142(numValues, floatValues, isDouble)

        ! Write z
        floatValues = pack(coordinates(:,:,:,3), .true.)
        i = tecdat142(numValues, floatValues, isDouble)

        ! Close file
        i = tecend142()

    end if

    ! Now write the flowfield
        
    write(filename,'("outputs/flowfield/time-",i6.6,".szplt")')it
    
    i = tecini142("Solution", "swirl u v w wx wy wz ctp", filename//NULLCHR, ".", &
        fileFormat, outputFileType, 1, 0)
    
    varTypes = 1
    shareVarFromZone = 0
    valueLocation = 1
    passiveVarList = 0
    
    ! Create zone
    i = teczne142( &
        "Flowfield", &
        zoneType, & 
        nyp, &
        mz_copy, &
        mx_copy, &
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
    
    ! Write swirl
    floatValues = pack(swirl(:,1:mz_copy,:), .true.)
    i = tecdat142(numValues, floatValues, isDouble)

    ! Write u
    floatValues = pack(u(:,1:mz_copy,:), .true.)
    i = tecdat142(numValues, floatValues, isDouble)
    
    ! Write v
    floatValues = pack(v(:,1:mz_copy,:), .true.)
    i = tecdat142(numValues, floatValues, isDouble)
    
    ! Write w
    floatValues = pack(w(:,1:mz_copy,:), .true.)
    i = tecdat142(numValues, floatValues, isDouble)
    
    ! Write wx
    floatValues = pack(wx(:,1:mz_copy,:), .true.)
    i = tecdat142(numValues, floatValues, isDouble)
    
    ! Write wy
    floatValues = pack(wy(:,1:mz_copy,:), .true.)
    i = tecdat142(numValues, floatValues, isDouble)
    
    ! Write wz
    floatValues = pack(wz(:,1:mz_copy,:), .true.)
    i = tecdat142(numValues, floatValues, isDouble)
    
    ! Write ctp
    floatValues = pack(ctp(:,1:mz_copy,:), .true.)
    i = tecdat142(numValues, floatValues, isDouble)
    
    ! Close file
    i = tecend142()

end subroutine

#else

! The following function is empty - compiled when tecio is not available
subroutine write_flowfield_plt(it, u, v, w, wx, wy, wz, swirl, ctp, &
        istart, iend, jstart, jend, kstart, kend, xl, zl)

    use grid_size
    implicit none

    integer :: it, istart, iend, jstart, jend, kstart, kend
    real, dimension(nyp,mz,mx) :: u, v, w, wx, wy, wz, ctp
    real, dimension(kend-kstart+1,jend-jstart+1,iend-istart+1) :: swirl 
    character (32) :: filename

    real :: x, y, z, delxm, delzm
    integer :: i, j, k, ii, jj, kk

    real xl, yl, zl

    ! THIS DOES NOTHING - IT IS ONLY A PLACEHOLDER

end subroutine

#endif


subroutine ubary(u,uavy,wfft1,wfft2,wfft3,wfft4) 
      use grid_size
      complex u(nyp,nz,nxh), ay(nyp)
      real uavy(nyp)
      real wfft1(1), wfft2(1), wfft3(1), wfft4(1)
      real sumre(512), sumim(512)
      common/data1/ pi,dt,theta,wavz,wavx,c,yc
      common/data2/ re,alpha,beta,xl,zl,yl,time,dyde,xstart,uinf
      common/trigx/ trigx(2*nx), ixfax(19)
      common/trigz/ trigz(2*nz), izfax(19)
      common/trigy/ sine(ny), cosine(ny), trigy(2*ny), iyfax(19)

      faccheb = 2.*float(ny)

!cccccc   compute average velocity profile

         do i =1,nyp
             ay(i)  =   u(i,1,1)
         end do

!cccccccccccc  take y-transform ccccccccccccccccc

      call ccheb(ay,wfft1,wfft2,wfft3,wfft4,sumre,sumim,ny,1,1,iyfax,trigy,sine,cosine)


        do i =1,nyp
            uavy(i) =  faccheb*real(ay(i))
        end do

end


subroutine ubulk(u,ub)
    use grid_size
      complex u(nyp,nz,nxh)
      real ub  

         ub =0.0

 
         do i =1,nyp,2
               im1 = i-1
               fac = float(1-im1**2)
               fac1 = 1./fac
             ub   =  ub  + real(u(i,1,1)*fac1)
         end do
end

! Below I'm implementing a solver to determine eigenvalues to calculate the
! swirl strength - Ryan 10/26/22
subroutine eig(u11,u12,u13,u21,u22,u23,u31,u32,u33, &
               lambda1,lambda2,lambda3)

! ==================================================== !
! This subroutine is used to get the eigenvalues of a  !
! 3x3 matrix, using a pre-determined expression for    !
! the values. Because of how it will be used in the    !
! DNS code, I have the matrix elements as separate     !
! inputs (u11 - u33).                                  !
!                                                      !
! The outputs are lambda1, lambda2, and lambda3, which !
! are the eigenvalues of the matrix U.                 !
!                                                      !
! Created 10/26/2022 by Ryan Kelly                     !
! ==================================================== !

    implicit none

    real u11,u12,u13,u21,u22,u23,u31,u32,u33 ! velocity gradients
    complex lambda1, lambda2, lambda3 ! Eigenvalues
    complex im

    im = (0.0,1.0) ! imaginary number i

    lambda1 = u11/3. + u22/3. + u33/3. + ((u11 + u22 + u33)**2./9. - (u11*u22)/3. +                   &
              (u12*u21)/3 - (u11*u33)/3 + (u13*u31)/3 - (u22*u33)/3 + (u23*u32)/3)/((u11 + u22 &
              + u33)**3./27 - ((u11 + u22 + u33)*(u11*u22 - u12*u21 + u11*u33 - u13*u31 +        &
              u22*u33 - u23*u32))/6 + (((u11 + u22 + u33)**3/27 - ((u11 + u22 + u33)*(u11*u22 - &
              u12*u21 + u11*u33 - u13*u31 + u22*u33 - u23*u32))/6 + (u11*u22*u33)/2 -          &
              (u11*u23*u32)/2 - (u12*u21*u33)/2 + (u12*u23*u31)/2 + (u13*u21*u32)/2 -          &
              (u13*u22*u31)/2)**2 - ((u11 + u22 + u33)**2/9 - (u11*u22)/3 + (u12*u21)/3 -        &
              (u11*u33)/3 + (u13*u31)/3 - (u22*u33)/3 + (u23*u32)/3)**3)**(1./2.) +                &
              (u11*u22*u33)/2 - (u11*u23*u32)/2 - (u12*u21*u33)/2 + (u12*u23*u31)/2 +          &
              (u13*u21*u32)/2 - (u13*u22*u31)/2)**(1./3.) + ((u11 + u22 + u33)**3/27 - ((u11 + u22 &
              + u33)*(u11*u22 - u12*u21 + u11*u33 - u13*u31 + u22*u33 - u23*u32))/6 + (((u11 + &
              u22 + u33)**3/27 - ((u11 + u22 + u33)*(u11*u22 - u12*u21 + u11*u33 - u13*u31 +    &
              u22*u33 - u23*u32))/6 + (u11*u22*u33)/2 - (u11*u23*u32)/2 - (u12*u21*u33)/2 +    &
              (u12*u23*u31)/2 + (u13*u21*u32)/2 - (u13*u22*u31)/2)**2 - ((u11 + u22 + u33)**2/9  &
              - (u11*u22)/3 + (u12*u21)/3 - (u11*u33)/3 + (u13*u31)/3 - (u22*u33)/3 +          &
              (u23*u32)/3)**3)**(1./2.) + (u11*u22*u33)/2 - (u11*u23*u32)/2 - (u12*u21*u33)/2 +    &
              (u12*u23*u31)/2 + (u13*u21*u32)/2 - (u13*u22*u31)/2)**(1./3.) 

    lambda2 = u11/3 + u22/3 + u33/3 - (3**(1./2.)*(((u11 + u22 + u33)**2/9 -                       &
              (u11*u22)/3 + (u12*u21)/3 - (u11*u33)/3 + (u13*u31)/3 - (u22*u33)/3 +            &
              (u23*u32)/3)/((u11 + u22 + u33)**3/27 - ((u11 + u22 + u33)*(u11*u22 - u12*u21 +   &
              u11*u33 - u13*u31 + u22*u33 - u23*u32))/6 + (((u11 + u22 + u33)**3/27 - ((u11 +   &
              u22 + u33)*(u11*u22 - u12*u21 + u11*u33 - u13*u31 + u22*u33 - u23*u32))/6 +      &
              (u11*u22*u33)/2 - (u11*u23*u32)/2 - (u12*u21*u33)/2 + (u12*u23*u31)/2 +          &
              (u13*u21*u32)/2 - (u13*u22*u31)/2)**2 - ((u11 + u22 + u33)**2/9 - (u11*u22)/3 +    &
              (u12*u21)/3 - (u11*u33)/3 + (u13*u31)/3 - (u22*u33)/3 + (u23*u32)/3)**3)**(1./2.) +  &
              (u11*u22*u33)/2 - (u11*u23*u32)/2 - (u12*u21*u33)/2 + (u12*u23*u31)/2 +          &
              (u13*u21*u32)/2 - (u13*u22*u31)/2)**(1./3.) - ((u11 + u22 + u33)**3/27 - ((u11 + u22 &
              + u33)*(u11*u22 - u12*u21 + u11*u33 - u13*u31 + u22*u33 - u23*u32))/6 + (((u11 + &
              u22 + u33)**3/27 - ((u11 + u22 + u33)*(u11*u22 - u12*u21 + u11*u33 - u13*u31 +    &
              u22*u33 - u23*u32))/6 + (u11*u22*u33)/2 - (u11*u23*u32)/2 - (u12*u21*u33)/2 +    &
              (u12*u23*u31)/2 + (u13*u21*u32)/2 - (u13*u22*u31)/2)**2 - ((u11 + u22 + u33)**2/9  &
              - (u11*u22)/3 + (u12*u21)/3 - (u11*u33)/3 + (u13*u31)/3 - (u22*u33)/3 +          &
              (u23*u32)/3)**3)**(1./2.) + (u11*u22*u33)/2 - (u11*u23*u32)/2 - (u12*u21*u33)/2 +    &
              (u12*u23*u31)/2 + (u13*u21*u32)/2 - (u13*u22*u31)/2)**(1./3.))*im)/2 - ((u11 + u22  &
              + u33)**2/9 - (u11*u22)/3 + (u12*u21)/3 - (u11*u33)/3 + (u13*u31)/3 - (u22*u33)/3 &
              + (u23*u32)/3)/(2*((u11 + u22 + u33)**3/27 - ((u11 + u22 + u33)*(u11*u22 -        &
              u12*u21 + u11*u33 - u13*u31 + u22*u33 - u23*u32))/6 + (((u11 + u22 + u33)**3/27 - &
              ((u11 + u22 + u33)*(u11*u22 - u12*u21 + u11*u33 - u13*u31 + u22*u33 -            &
              u23*u32))/6 + (u11*u22*u33)/2 - (u11*u23*u32)/2 - (u12*u21*u33)/2 +              &
              (u12*u23*u31)/2 + (u13*u21*u32)/2 - (u13*u22*u31)/2)**2 - ((u11 + u22 + u33)**2/9  &
              - (u11*u22)/3 + (u12*u21)/3 - (u11*u33)/3 + (u13*u31)/3 - (u22*u33)/3 +          &
              (u23*u32)/3)**3)**(1./2.) + (u11*u22*u33)/2 - (u11*u23*u32)/2 - (u12*u21*u33)/2 +    & 
              (u12*u23*u31)/2 + (u13*u21*u32)/2 - (u13*u22*u31)/2)**(1./3.)) - ((u11 + u22 +      &
              u33)**3/27 - ((u11 + u22 + u33)*(u11*u22 - u12*u21 + u11*u33 - u13*u31 + u22*u33  &
              - u23*u32))/6 + (((u11 + u22 + u33)**3/27 - ((u11 + u22 + u33)*(u11*u22 - u12*u21 &
              + u11*u33 - u13*u31 + u22*u33 - u23*u32))/6 + (u11*u22*u33)/2 - (u11*u23*u32)/2  &
              - (u12*u21*u33)/2 + (u12*u23*u31)/2 + (u13*u21*u32)/2 - (u13*u22*u31)/2)**2 -     &
              ((u11 + u22 + u33)**2/9 - (u11*u22)/3 + (u12*u21)/3 - (u11*u33)/3 + (u13*u31)/3 - &
              (u22*u33)/3 + (u23*u32)/3)**3)**(1./2.) + (u11*u22*u33)/2 - (u11*u23*u32)/2 -        &
              (u12*u21*u33)/2 + (u12*u23*u31)/2 + (u13*u21*u32)/2 - (u13*u22*u31)/2)**(1./3.)/2

    lambda3 = u11/3 + u22/3 + u33/3 + (3**(1./2.)*(((u11 + u22 + u33)**2/9 - (u11*u22)/3           &
              + (u12*u21)/3 - (u11*u33)/3 + (u13*u31)/3 - (u22*u33)/3 + (u23*u32)/3)/((u11 +   &
              u22 + u33)**3/27 - ((u11 + u22 + u33)*(u11*u22 - u12*u21 + u11*u33 - u13*u31 +    &
              u22*u33 - u23*u32))/6 + (((u11 + u22 + u33)**3/27 - ((u11 + u22 + u33)*(u11*u22 - &
              u12*u21 + u11*u33 - u13*u31 + u22*u33 - u23*u32))/6 + (u11*u22*u33)/2 -          &
              (u11*u23*u32)/2 - (u12*u21*u33)/2 + (u12*u23*u31)/2 + (u13*u21*u32)/2 -          &
              (u13*u22*u31)/2)**2 - ((u11 + u22 + u33)**2/9 - (u11*u22)/3 + (u12*u21)/3 -        &
              (u11*u33)/3 + (u13*u31)/3 - (u22*u33)/3 + (u23*u32)/3)**3)**(1./2.) +                &
              (u11*u22*u33)/2 - (u11*u23*u32)/2 - (u12*u21*u33)/2 + (u12*u23*u31)/2 +          &
              (u13*u21*u32)/2 - (u13*u22*u31)/2)**(1./3.) - ((u11 + u22 + u33)**3/27 - ((u11 + u22 &
              + u33)*(u11*u22 - u12*u21 + u11*u33 - u13*u31 + u22*u33 - u23*u32))/6 + (((u11 + &
              u22 + u33)**3/27 - ((u11 + u22 + u33)*(u11*u22 - u12*u21 + u11*u33 - u13*u31 +    &
              u22*u33 - u23*u32))/6 + (u11*u22*u33)/2 - (u11*u23*u32)/2 - (u12*u21*u33)/2 +    &
              (u12*u23*u31)/2 + (u13*u21*u32)/2 - (u13*u22*u31)/2)**2 - ((u11 + u22 + u33)**2/9  &
              - (u11*u22)/3 + (u12*u21)/3 - (u11*u33)/3 + (u13*u31)/3 - (u22*u33)/3 +          &
              (u23*u32)/3)**3)**(1./2.) + (u11*u22*u33)/2 - (u11*u23*u32)/2 - (u12*u21*u33)/2 +    &
              (u12*u23*u31)/2 + (u13*u21*u32)/2 - (u13*u22*u31)/2)**(1./3.))*im)/2 - ((u11 + u22  &
              + u33)**2/9 - (u11*u22)/3 + (u12*u21)/3 - (u11*u33)/3 + (u13*u31)/3 - (u22*u33)/3 &
              + (u23*u32)/3)/(2*((u11 + u22 + u33)**3/27 - ((u11 + u22 + u33)*(u11*u22 -        &
              u12*u21 + u11*u33 - u13*u31 + u22*u33 - u23*u32))/6 + (((u11 + u22 + u33)**3/27 - &
              ((u11 + u22 + u33)*(u11*u22 - u12*u21 + u11*u33 - u13*u31 + u22*u33 -            &
              u23*u32))/6 + (u11*u22*u33)/2 - (u11*u23*u32)/2 - (u12*u21*u33)/2 +              &
              (u12*u23*u31)/2 + (u13*u21*u32)/2 - (u13*u22*u31)/2)**2 - ((u11 + u22 + u33)**2/9  &
              - (u11*u22)/3 + (u12*u21)/3 - (u11*u33)/3 + (u13*u31)/3 - (u22*u33)/3 +          & 
              (u23*u32)/3)**3)**(1./2.) + (u11*u22*u33)/2 - (u11*u23*u32)/2 - (u12*u21*u33)/2 +    &
              (u12*u23*u31)/2 + (u13*u21*u32)/2 - (u13*u22*u31)/2)**(1./3.)) - ((u11 + u22 +      &
              u33)**3/27 - ((u11 + u22 + u33)*(u11*u22 - u12*u21 + u11*u33 - u13*u31 + u22*u33  &
              - u23*u32))/6 + (((u11 + u22 + u33)**3/27 - ((u11 + u22 + u33)*(u11*u22 - u12*u21 & 
              + u11*u33 - u13*u31 + u22*u33 - u23*u32))/6 + (u11*u22*u33)/2 - (u11*u23*u32)/2  &
              - (u12*u21*u33)/2 + (u12*u23*u31)/2 + (u13*u21*u32)/2 - (u13*u22*u31)/2)**2 -     &
              ((u11 + u22 + u33)**2/9 - (u11*u22)/3 + (u12*u21)/3 - (u11*u33)/3 + (u13*u31)/3 - &
              (u22*u33)/3 + (u23*u32)/3)**3)**(1./2.) + (u11*u22*u33)/2 - (u11*u23*u32)/2 -        &
              (u12*u21*u33)/2 + (u12*u23*u31)/2 + (u13*u21*u32)/2 - (u13*u22*u31)/2)**(1./3.)/2 
end subroutine

subroutine calcswirl(u11,u21,u31,u12,u22,u32,u13,u23,u33,swirl)
! This subroutine calculates the swirl strength, defined as the magnitude of the
! complex conjugate of the eigenvalues of the velocity gradient tensor

! Calls subroutine eig(...) to calculate eigenvalues of velocity gradient tensor

    use grid_size

    implicit none

    real :: u11,u12,u13,u21,u22,u23,u31,u32,u33
    complex :: l1, l2, l3
    real :: swirl, l1_im, l2_im, l3_im

    call eig(u11,u12,u13,u21,u22,u23,u31,u32,u33,l1,l2,l3)

    l1_im = aimag(l1)
    l2_im = aimag(l2)
    l3_im = aimag(l3)

    swirl = sqrt(l1_im**2 + l2_im**2 + l3_im**2) 

    ! Certain velocity gradient tensors produce NaNs using the equations in
    ! eig(...). When this happens, the actual eigenvalues are all 0. This is a
    ! quick check to see if the swirl is a NaN. If it's any real value,
    ! floor(swirl/swirl) will always be 1; otherwise it's a NaN (or Inf, but
    ! that doesn't happen in this case). - Ryan 10/28/22
    if (floor(swirl/swirl) .ne. 1) then
        swirl = 0
    end if

end subroutine


subroutine derivscji(u,v,w,u11,u12,u13,u21,u22,u23,u31,u32,u33,Lu,Lv,Lw)

!***************************************************************
!  calculate velocity gradient tensor and derivs of conformation tensor  
!  tj and rah 6/16/2014
!
!  edited for only velocity gradient - Ryan 6/1/23
!  edited for laplacian - Ryan 7/24/23
!
!************************************************************************
    use grid_size

    complex u(nyp,nz,nxh),v(nyp,nz,nxh),w(nyp,nz,nxh)

    ! velocity gradient tensor 
    complex u11(nyp,nz,nxh),u12(nyp,nz,nxh),u13(nyp,nz,nxh)
    complex u21(nyp,nz,nxh),u22(nyp,nz,nxh),u23(nyp,nz,nxh) 
    complex u31(nyp,nz,nxh),u32(nyp,nz,nxh),u33(nyp,nz,nxh)

    ! Laplacian
    complex,dimension(nyp,nz,nxh) :: Lu,Lv,Lw


    complex wrk1(nyp,nz,nxh),wrk2(nyp,nz,nxh),wrk3(nyp,nz,nxh)   

    complex im
    common/data1/ pi,dt,theta,wavz(nz),wavx(nxh),c(nyp),yc(nyp)
   
         im = (0.0,1.0)

    ! Calculate the velocity gradient tensor

    call cderiv(u,u12) !du/dy
    call cderiv(v,u22) !dv/dy
    call cderiv(w,u32) !dw/dy

    ! Laplacian terms
    call cderiv(u12,wrk1) ! d^2 u / dy^2
    call cderiv(u22,wrk2) ! d^2 v / dy^2
    call cderiv(u32,wrk3) ! d^2 w / dy^2

    do k=1,nxh
        do j=1,nz
            do i=1,nyp
                u11(i,j,k) = im*wavx(k)*u(i,j,k)
                u13(i,j,k) = im*wavz(j)*u(i,j,k)
                u21(i,j,k) = im*wavx(k)*v(i,j,k)
                u23(i,j,k) = im*wavz(j)*v(i,j,k)
                u31(i,j,k) = im*wavx(k)*w(i,j,k)
                u33(i,j,k) = im*wavz(j)*w(i,j,k)
    
                ! Laplacian terms
                Lu(i,j,k) = im*wavx(k)*u11(i,j,k) + wrk1(i,j,k) + im*wavz(j)*u13(i,j,k)
                Lv(i,j,k) = im*wavx(k)*u21(i,j,k) + wrk2(i,j,k) + im*wavz(j)*u23(i,j,k)
                Lw(i,j,k) = im*wavx(k)*u31(i,j,k) + wrk3(i,j,k) + im*wavz(j)*u33(i,j,k)

            end do
        end do
    end do

end

