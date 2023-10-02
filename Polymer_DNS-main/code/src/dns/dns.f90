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
    
! ------------------------------------------------------------------------- !
!                Adding in polymer stuff - Ryan 2-18-22                     ! 
! ------------------------------------------------------------------------- !

    ! General (non-specified) variables
    complex scalar(nyp,nz,nxh) ! Passive scalar
    complex sclx(nyp,nz,nxh), scly(nyp,nz,nxh), sclz(nyp,nz,nxh)
    complex csource(nyp,nz,nxh), psource(nyp,nz,nx) ! Scalar source added 8/3/22
    complex scn(nyp,nz,nxh), scnm1(nyp,nz,nxh) ! Scalar non-linear sources
    complex c11n(nyp,nz,nxh), c12n(nyp,nz,nxh), c13n(nyp,nz,nxh)
    complex c23n(nyp,nz,nxh), c22n(nyp,nz,nxh), c33n(nyp,nz,nxh) ! nonlinear terms for conformation
    complex c11nm1(nyp,nz,nxh), c12nm1(nyp,nz,nxh), c13nm1(nyp,nz,nxh)
    complex c23nm1(nyp,nz,nxh), c22nm1(nyp,nz,nxh), c33nm1(nyp,nz,nxh) ! nonlinear terms for conformation
    complex str11n(nyp,nz,nxh), str12n(nyp,nz,nxh), str13n(nyp,nz,nxh)
    complex str23n(nyp,nz,nxh), str22n(nyp,nz,nxh), str33n(nyp,nz,nxh) ! nonlinear terms for conformation
    complex str11nm1(nyp,nz,nxh), str12nm1(nyp,nz,nxh), str13nm1(nyp,nz,nxh)
    complex str23nm1(nyp,nz,nxh), str22nm1(nyp,nz,nxh), str33nm1(nyp,nz,nxh) ! nonlinear terms for conformation
    complex bcscltop(nz,nxh),bcsclbot(nz,nxh) ! BCs for scalar

    real ur(nyp,nz,nx), vr(nyp,nz,nx)

    ! Body forces - check to see if these need to be changed for IB
    complex forcebdy1(nyp,nz,nxh),forcebdy2(nyp,nz,nxh)

    ! GMU output variables
    real uavy(nyp),vavy(nyp),wavy(nyp)
    real u1b(nyp),u2b(nyp),u3b(nyp)
    real uvar(nyp),vvar(nyp),wvar(nyp)
    real scalaravy(nyp),sclvar(nyp)
    real c11avy(nyp),c22avy(nyp),c33avy(nyp),c12avy(nyp),c13avy(nyp),c23avy(nyp)
    real con, trc, dtdybw, dudybw, c11b, c22b, c33b

    ! LES arrays
    complex T1(nyp,nz,nxh), T2(nyp,nz,nxh), T3(nyp,nz,nxh) 
    real TKE(nyp,nz,nx), TKE1(nyp,nz,nx)
    complex znuc(nyp,nz,nxh)
    real znub, znub1, ub

    ! Conformation tensor
    complex c11(nyp,nz,nxh), c12(nyp,nz,nxh), c13(nyp,nz,nxh)
    complex c21(nyp,nz,nxh), c22(nyp,nz,nxh), c23(nyp,nz,nxh)
    complex c31(nyp,nz,nxh), c32(nyp,nz,nxh), c33(nyp,nz,nxh)

    ! Velocity gradient tensor
    complex u11(nyp,nz,nxh), u12(nyp,nz,nxh), u13(nyp,nz,nxh)
    complex u21(nyp,nz,nxh), u22(nyp,nz,nxh), u23(nyp,nz,nxh)
    complex u31(nyp,nz,nxh), u32(nyp,nz,nxh), u33(nyp,nz,nxh)
    complex,dimension(nyp,nz,nxh) :: Lu,Lv,Lw ! Laplacian terms - Ryan 7/24/23

    ! Derivatives of conformation tensor
    complex dc111(nyp,nz,nxh), dc112(nyp,nz,nxh), dc113(nyp,nz,nxh)
    complex dc211(nyp,nz,nxh), dc212(nyp,nz,nxh), dc213(nyp,nz,nxh)
    complex dc311(nyp,nz,nxh), dc312(nyp,nz,nxh), dc313(nyp,nz,nxh)
    
    complex dc121(nyp,nz,nxh), dc122(nyp,nz,nxh), dc123(nyp,nz,nxh)
    complex dc221(nyp,nz,nxh), dc222(nyp,nz,nxh), dc223(nyp,nz,nxh)
    complex dc321(nyp,nz,nxh), dc322(nyp,nz,nxh), dc323(nyp,nz,nxh)

    complex dc131(nyp,nz,nxh), dc132(nyp,nz,nxh), dc133(nyp,nz,nxh)
    complex dc231(nyp,nz,nxh), dc232(nyp,nz,nxh), dc233(nyp,nz,nxh)
    complex dc331(nyp,nz,nxh), dc332(nyp,nz,nxh), dc333(nyp,nz,nxh)

    ! Nonlinear terms for conformation tensor
    complex c11NL(nyp,nz,nxh), c12NL(nyp,nz,nxh), c13NL(nyp,nz,nxh)
    complex c22NL(nyp,nz,nxh), c23NL(nyp,nz,nxh)
    complex c33NL(nyp,nz,nxh)

    ! Polymer stresses used for targeting
    complex qp11(nyp,nz,nxh), qp12(nyp,nz,nxh), qp13(nyp,nz,nxh)
    complex qp22(nyp,nz,nxh), qp23(nyp,nz,nxh)
    complex qp33(nyp,nz,nxh)

    ! Adding input variables that are not declared in the original code
    real omearth,rotang,grav,virtual
    real zlmax,tpoly,qbeta,diff,diffpoly,deltat,c11z,c22z,c33z,alpha_poly,mpoly
    real ampbdy1,ampbdy2,c11amp,tfrac
    real att,btt,abb,bbb ! BCs for scalar ?
    real polyrate
    
    integer ipolyflag,itarget,ipeter,scl_flag,src_start,src_stop

    ! These are all the variables from the polymer code
! ====================================================================== !

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

    real wrk11(nyp,nz,nx),wrk12(nyp,nz,nx),wrk13(nyp,nz,nx)
    real wrk21(nyp,nz,nx),wrk22(nyp,nz,nx),wrk23(nyp,nz,nx)
    real wrk31(nyp,nz,nx),wrk32(nyp,nz,nx),wrk33(nyp,nz,nx)
  
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
! Adding in common blocks from polymer code - Ryan 2-18-22
    common/thermcor/ omearth,rotang,grav,virtual
    common/polymer1/ zlmax,tpoly,qbeta,diff,diffpoly,deltat,c11z,c22z,c33z,alpha_poly,mpoly
    common/polymer2/ ampbdy1,ampbdy2,c11amp,tfrac,polyrate,src_start,src_stop
    common/polymer3/ ipolyflag,itarget,ipeter,scl_flag
    common/fbdy/ forcebdy1, forcebdy2
! Done adding common blocks (I regrouped some, but that may change as I see how they're used)
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
   
    ! Adding the new variables to this subroutine - Ryan 2-18-22     
    it = 0

    call initial(u,u0,v,w,w0,omx,omy,omz,fn,fnm1,gn,gnm1,h1n,h1nm1,h3n,h3nm1,      &
     wrk1,wrkc,wfft1,wfft2,wfft3,wfft4,w1,w2,bctop,bcbot,a,t,                      &
     scalar,sclx,scly,sclz,scn,scnm1,T1,T2,T3,c11,c12,c13,c21,c22,c23,c31,c32,c33, &
     u11,u12,u13,u21,u22,u23,u31,u32,u33,dc111,dc112,dc113,dc211,dc212,dc213,      &
     dc311,dc312,dc313,dc121,dc122,dc123,dc221,dc222,dc223,dc321,dc322,dc323,      &
     dc131,dc132,dc133,dc231,dc232,dc233,dc331,dc332,dc333,c11n,c12n,c13n,         &
     c22n,c23n,c33n,c11nm1,c12nm1,c13nm1,c22nm1,c23nm1,c33nm1,str11n,str12n,       &
     str13n,str22n,str23n,str33n,str11nm1,str12nm1,str13nm1,str22nm1,str23nm1,     &
     str33nm1,qp11,qp12,qp13,qp22,qp23,qp33,csource,Lu,Lv,Lw)


  ! -----------------------------------------------------------------
  ! ------------------------  main loop  ----------------------------
  ! -----------------------------------------------------------------
   
    do it = 1, nsteps
  
       
       write(*,22) it
  22   format('  it = ',i7)

! --------------------------------------------------------------------------- !
!                   Solving the conformation tensor - Ryan 2-18-22            !
! --------------------------------------------------------------------------- !

    ! Reset conformation tensor to unstretched state when source starts - Ryan 9/27/22
    if (it .eq. (src_start-1)) then 

    ! reset just before src_start so it evolves one time step. This allows us to subsequently integrate in time

      do i = 1,nyp
        do j = 1,nz
          do k = 1,nx
            wrk11(i,j,k) = 1.0
            wrk22(i,j,k) = 1.0
            wrk33(i,j,k) = 1.0
    
            wrk12(i,j,k) = 0.0
            wrk13(i,j,k) = 0.0
            wrk21(i,j,k) = 0.0
            wrk23(i,j,k) = 0.0
            wrk31(i,j,k) = 0.0
            wrk32(i,j,k) = 0.0
          end do
        end do
      end do


      call scram(wrk11,c11)
      call xyzfft(c11,w1,w2,wfft1,wfft2,wfft3,wfft4,-1)
      call scram(wrk12,c12)
      call xyzfft(c12,w1,w2,wfft1,wfft2,wfft3,wfft4,-1)
      call scram(wrk13,c13)
      call xyzfft(c13,w1,w2,wfft1,wfft2,wfft3,wfft4,-1)
      call scram(wrk21,c21)
      call xyzfft(c21,w1,w2,wfft1,wfft2,wfft3,wfft4,-1)
      call scram(wrk22,c22)
      call xyzfft(c22,w1,w2,wfft1,wfft2,wfft3,wfft4,-1)
      call scram(wrk23,c23)
      call xyzfft(c23,w1,w2,wfft1,wfft2,wfft3,wfft4,-1)
      call scram(wrk31,c31)
      call xyzfft(c31,w1,w2,wfft1,wfft2,wfft3,wfft4,-1)
      call scram(wrk32,c32)
      call xyzfft(c32,w1,w2,wfft1,wfft2,wfft3,wfft4,-1)
      call scram(wrk33,c33)
      call xyzfft(c33,w1,w2,wfft1,wfft2,wfft3,wfft4,-1)

    end if

    ! Throughout the code, I have this logic gate, but for now it's only
    ! relevant for scl_flag = 2 (polymer emitted from a particle)  
    if (ipolyflag .eq. 1 .and. it .ge. src_start) then ! Note this is only called after C_ij has been initialized


!    SOLVE FOR CONFORMATION TENSOR WITH DIFFUSION AND NEUMANN BCS

! parallel sections?
      call polynl(c11nl,c12nl,c13nl,c22nl,c23nl,c33nl,                    &
                  c11n,c12n,c13n,c22n,c23n,c33n,                          &
                  c11nm1,c12nm1,c13nm1,c22nm1,c23nm1,c33nm1,              &
                  str11n,str12n,str13n,str22n,str23n,str33n,              &
                  str11nm1,str12nm1,str13nm1,str22nm1,str23nm1,str33nm1)

    !$omp parallel sections default(shared) private(wrkc,wrk1,bcbot,bctop)

    !$omp section
      call polyrhs(c11,wrkc,wrk1,c11nl,bcbot,bctop)
    !$omp section
      call polyrhs(c22,wrkc,wrk1,c22nl,bcbot,bctop)
    !$omp section
      call polyrhs(c33,wrkc,wrk1,c33nl,bcbot,bctop)
    !$omp section
      call polyrhs(c12,wrkc,wrk1,c12nl,bcbot,bctop)
    !$omp section
      call polyrhs(c13,wrkc,wrk1,c13nl,bcbot,bctop)
    !$omp section
      call polyrhs(c23,wrkc,wrk1,c23nl,bcbot,bctop)

    !$omp end parallel sections

!      copy n  to nm1

      !$omp parallel do
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
      !$omp end parallel do

      g = (re*diffpoly)/(dt*theta)    ! parameter for timestep dep., theta=.5 cn method
      att = 0.0  
      btt = 1.0  ! no flux 
      abb = 0.0
      bbb = 1.0  !no flux rah 4/1/2016

! Parallelize?
      do k=1,nxh
          x = wavx(k)**2
    
          call penta(c11(1,1,k),x,g,dyde,ib,wavz,c, &
                     att,btt,bctop(1,k),abb,bbb,    &
                     bcbot(1,k),a,t,wrk1)
          call penta(c22(1,1,k),x,g,dyde,ib,wavz,c, &
                     att,btt,bctop(1,k),abb,bbb,    &
                     bcbot(1,k),a,t,wrk1)
          call penta(c33(1,1,k),x,g,dyde,ib,wavz,c, &
                     att,btt,bctop(1,k),abb,bbb,    &
                     bcbot(1,k),a,t,wrk1)
          call penta(c12(1,1,k),x,g,dyde,ib,wavz,c, &
                     att,btt,bctop(1,k),abb,bbb,    &
                     bcbot(1,k),a,t,wrk1)
          call penta(c13(1,1,k),x,g,dyde,ib,wavz,c, &
                     att,btt,bctop(1,k),abb,bbb,    &
                     bcbot(1,k),a,t,wrk1)
          call penta(c23(1,1,k),x,g,dyde,ib,wavz,c, &
                     att,btt,bctop(1,k),abb,bbb,    &
                     bcbot(1,k),a,t,wrk1)
    
      end do

!   Now use symmetry to get other components of C tensor

! Parallelize?
      !$omp parallel do
      do k=1,nxh
          do j=1,nz
              do i=1,nyp
                  c21(i,j,k)=c12(i,j,k)
                  c31(i,j,k)=c13(i,j,k)
                  c32(i,j,k)=c23(i,j,k)
              end do
          end do
      end do
      !$omp end parallel do

!      Norms

       call norm(c11)
       call norm(c12)
       call norm(c13)
       call norm(c21)
       call norm(c22)
       call norm(c23)
       call norm(c31)
       call norm(c32)
       call norm(c33)
     end if ! ipolyflag = 1
!  End of conformation tensor solutions

      if (scl_flag .ge. 1) then
!--  Evaluate rhs of time discrete equation for the passive scalar
!--  Also sets bcsclbot=bcscltop=(k,0) or whatever is desired.
      call sclrhs(scalar,wrkc,wrk1,scn,scnm1,bcsclbot,bcscltop)
      call temperbc(bcscltop,bcsclbot,w1,w2,wfft1,wfft2,wfft3,wfft4)

! Copy scn to scnm1
! Parallelize?
      !$omp parallel do
      do k=1,nxh
          do j=1,nz
             do i=1,nyp
                scnm1(i,j,k) = scn(i,j,k)
             end do
          end do
      end do
      !$omp end parallel do

!       To now used mixed bc's for the temperature field we need to use PENTA 
!       Here we will use constant flux conditions at the free surface and 
!       temperature conditions on the bottom (RAH  7/11/96)


      x = 0.0              ! parameter for x (or kx) dependency
      g = (re*diff)/(dt*theta)    ! parameter for timestep dep., theta=.5 cn method
      ib = 0
      att = 1.0 ! zero temp on top 1/11/10
      btt = 0.0
      abb = 1.0
      bbb = 0.0 ! zero temp on bot 7/12/07 rah 

! I don't think this is parallelizable
      do k=1,nxh
          x = wavx(k)**2

          call penta(scalar(1,1,k),x,g,dyde,ib,wavz,c,att,btt,bcscltop(1,k),abb,bbb,bcsclbot(1,k),a,t,wrk1)
      end do
! End of C tensor calculations
! ============================================================================= !
      end if ! scl_flag = 1

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
 
! ------------------------------------------------------------------------------------- !
!                           Adding GMU Printing Variables (Ryan 7/15/22)                !
! ------------------------------------------------------------------------------------- !
    

    if (print3d .eq. -1) then
      if (it .eq. 1) then
          call system('rm -rf GMU_out') ! Delete old output directory
          call system('mkdir GMU_out')  ! Re-create a clean directory
      end if

      if(mod((it),iprnfrq).eq.0) then
         call output(u,v,w,scalar,c11,c12,c22,wrk1,wrkc,wfft1,wfft2,wfft3,wfft4,w1,w2)
      end if


!    at frequencies specified by idmpfrq we compute some 
!    global quantities on the fly and dump them out

          open(unit=40,file='GMU_out/heat_flux',form='formatted',status='unknown')
          open(unit=41,file='GMU_out/shear_stress',form='formatted',status='unknown')
          open(unit=42,file='GMU_out/mass_flux',form='formatted',status='unknown')
          open(unit=43,file='GMU_out/kinetic_energy',form='formatted',status='unknown')
          open(unit=70,file='GMU_out/kinetic_energy_total',form='formatted',status='unknown')
          open(unit=44,file='GMU_out/trace_conformation',form='formatted',status='unknown')

!----------------compute global parameters---------------------------
 
      if(mod((it-1),iprnfrq).eq.0) then 

         call c1derbw(scalar,wrkc,1)
         call c1derbw(u,wrkc,2)
         call ubulk(u,ub)
         call ubulk(c11,c11b)
         call ubulk(c22,c22b)
         call ubulk(c33,c33b)
         trc = c11b+c22b+c33b
         dtdybw = real(wrkc(1,1,1))     !heat flux
         dudybw = real(wrkc(2,1,1))     !shear stress
         write(40,1234) time,dtdybw
         write(41,1234) time,dudybw
         write(42,1234) time,ub         !mass flux
         write(44,1234) time,trc        !trace
           endif

 1234        format(e20.13,1x,e20.13)

!-------------------compute mean and rms profies-----------------------

      if(mod((it-1),iprnfrq).eq.0) then
           
          con = float(nx*nz)

          call ubary(u,uavy,wfft1,wfft2,wfft3,wfft4) 
          call ubary(v,vavy,wfft1,wfft2,wfft3,wfft4)
          call ubary(w,wavy,wfft1,wfft2,wfft3,wfft4)
          call ubary(scalar,scalaravy,wfft1,wfft2,wfft3,wfft4)

!         conformation tensor results 
          call ubary(c11,c11avy,wfft1,wfft2,wfft3,wfft4)
          call ubary(c22,c22avy,wfft1,wfft2,wfft3,wfft4)
          call ubary(c33,c33avy,wfft1,wfft2,wfft3,wfft4)
          call ubary(c12,c12avy,wfft1,wfft2,wfft3,wfft4)
          call ubary(c13,c13avy,wfft1,wfft2,wfft3,wfft4)
          call ubary(c23,c23avy,wfft1,wfft2,wfft3,wfft4)

  
          call ccopy (nz*nxh*nyp,u,1,wrkc,1)
          call xyzfft(wrkc,w1,w2,wfft1,wfft2,wfft3,wfft4,1)
          call unscram(wrkc,ur)

          call ccopy (nz*nxh*nyp,v,1,wrkc,1)
          call xyzfft(wrkc,w1,w2,wfft1,wfft2,wfft3,wfft4,1)
          call unscram(wrkc,vr)

          call ccopy (nz*nxh*nyp,w,1,wrkc,1)
          call xyzfft(wrkc,w1,w2,wfft1,wfft2,wfft3,wfft4,1)
          call unscram(wrkc,wrk1)

!         calculate tke

          do i=1,nyp
          do j=1,nz
          do k=1,nx
               tke(i,j,k)=(ur(i,j,k)-uavy(i))**2+(vr(i,j,k)-vavy(i))**2+(wrk1(i,j,k)-wavy(i))**2
               tke1(i,j,k)=(ur(i,j,k))**2+(vr(i,j,k))**2+(wrk1(i,j,k))**2
               tke(i,j,k)=tke(i,j,k)/2.0
               tke1(i,j,k)=tke1(i,j,k)/2.0
          enddo
          enddo
          enddo

          call scram(tke,znuc)
          call xyzfft(znuc,w1,w2,wfft1,wfft2,wfft3,wfft4,-1)
          call ubulk(znuc,znub)
          call scram(tke1,znuc)
          call xyzfft(znuc,w1,w2,wfft1,wfft2,wfft3,wfft4,-1)
          call ubulk(znuc,znub1)


          call ccopy (nz*nxh*nyp,scalar,1,wrkc,1)
          call xyzfft(wrkc,w1,w2,wfft1,wfft2,wfft3,wfft4,1)
          call unscram(wrkc,wrk1)

          do i=1,nyp
             sclvar(i) = 0.0 
          do j=1,nz
          do k=1,nx
            sclvar(i) = sclvar(i) + (wrk1(i,j,k)-scalaravy(i))**2
          enddo
          enddo
          enddo


          do i=1,nyp
             uavy(i) = uavy(i)
             vavy(i) = vavy(i)
             wavy(i) = wavy(i)
             scalaravy(i) = scalaravy(i)
             uvar(i) = sqrt(uvar(i)/con)
             vvar(i) = sqrt(vvar(i)/con)
             wvar(i) = sqrt(wvar(i)/con)
             sclvar(i) = sqrt(sclvar(i)/con)
          end do


          write(43,1235) time,znub ! tke dumped in fort.35 !
          write(70,1235) time,znub1
          open(36,file="GMU_out/fort.36")
            write(36,*) 'time ='
          write(36,1235) time
             write(36,*) 'average u,v,w,temp'
          write(36,1236)(uavy(i),i=1,nyp)
          write(36,1236)(vavy(i),i=1,nyp)
          write(36,1236)(wavy(i),i=1,nyp)
          write(36,1236)(scalaravy(i),i=1,nyp)
             write(36,*) 'conformation' 
          write(36,1236)(c11avy(i),i=1,nyp)
          write(36,1236)(c22avy(i),i=1,nyp)
          write(36,1236)(c33avy(i),i=1,nyp)
             write(36,*) 'rms u,v,w,temp'
          write(36,1236)(uvar(i),i=1,nyp)
          write(36,1236)(vvar(i),i=1,nyp)
          write(36,1236)(wvar(i),i=1,nyp)
          write(36,1236)(sclvar(i),i=1,nyp)
          close(36)
      endif 

   
 1235       format(e20.13,1x,e20.13)
 1236       format(4e20.13)         
 1237       format(4e12.5) 
 1238       format(4e20.13)
    end if
! ------------------------------------------------------------------------------------- ! 
! ------------------------------------------------------------------------------------- !

  !    normalize to ensure conjg sym.
      
       call norm(u)      
       call norm(w)
       call norm(scalar)
 
  !    calculate the other vorticity components.
      
       call vort(u,v,w,scalar,omx,omz,sclx,scly,sclz,wrkc) ! more arguments - Ryan 6/24/22
       call norm(omx)
       call norm(omz)

! ------------------------------------------------------------------------------------- !
!                   Adding main function calls for C tensor - Ryan                      !
! ------------------------------------------------------------------------------------- !
        call norm(omy)
        call norm(sclx)
        call norm(scly)
        call norm(sclz)

! Calculate derivates needed for conformation tensor sources
        call derivscji(u,v,w,c11,c12,c13,c21,c22,c23,c31,c32,c33,                &
                       u11,u12,u13,u21,u22,u23,u31,u32,u33,Lu,Lv,Lw,             &
                       dc111,dc112,dc113,dc211,dc212,dc213,dc311,dc312,dc313,    &
                       dc121,dc122,dc123,dc221,dc222,dc223,dc321,dc322,dc323,    &
                       dc131,dc132,dc133,dc231,dc232,dc233,dc331,dc332,dc333)    

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

    if (ipolyflag .eq. 1 .and. it .ge. src_start-1) then

! Need to normalize all arrays 

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

! ===================================================================================== !
    end if
  !    at frequencies specified by iprnfrq calculate the pressure
  !    and output u,v,w, and p and determine statistics.    

! Adding additional arguments - Ryan 2-24-22
!       if(mod(it,nsteps).eq.0) then
!          call output(u,v,w,urms,vrms,wrms,umn,wrk1,wrkc,wfft1,wfft2,wfft3,wfft4,w1,w2,scalar)
!       end if
  
!  !   print out the chebyshev spectrum of the velocity.
!        
!       if(mod((it-1),iprnfrq).eq.0) then
!          do i = 1, nyp
!             uchbeng(i) = 0.
!             vchbeng(i) = 0.
!             wchbeng(i) = 0.
! !            write(*,771) float(i), uchbeng(i),vchbeng(i), wchbeng(i) ! DEBUG
!          end do
!  
!  !$omp parallel do 
!          do i = 1,nyp
!             do j = 1, nz
!                do  k = 1, nxh
!                   uchbeng(i) = uchbeng(i) + u(i,j,k)*conjg(u(i,j,k))
!                   vchbeng(i) = vchbeng(i) + v(i,j,k)*conjg(v(i,j,k))
!                   wchbeng(i) = wchbeng(i) + w(i,j,k)*conjg(w(i,j,k))
!  !               write(*,*) 'one ', float(i), u(i,j,k), v(i,j,k), w(i,j,k) ! DEBUG
!  !               write(*,*) float(i), conjg(u(i,j,k)), conjg(v(i,j,k)), conjg(w(i,j,k)) ! DEBUG
!  !               write(*,771) float(i), uchbeng(i),vchbeng(i), wchbeng(i) ! DEBUG
!                end do
!             end do
!          end do
!  !$omp end parallel do
!  
! !     do i = 1, nyp
! !        write(*,771) float(i), uchbeng(i)/float(nz*nxh),vchbeng(i)/float(nz*nxh), wchbeng(i)/float(nz*nxh)
! !     end do
!      end if
 ! 771  format('i,u,v,w chbengy(i) = ',e10.3,e11.3,e11.3,e11.3)
  
   
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
! Perhaps I should add an 'if' statement here, but I may need to be careful if
! I'm using it in the middle of paralle sections - Ryan 6/24/22 

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

  if (scl_flag .ge. 1) then 
  !$omp parallel sections default(shared) private(wfft1,wfft2,wfft3,wfft4)      
  !$omp section   
      call yfft(scalar,wfft1,wfft2,wfft3,wfft4,is) 
  !$omp section   
      call yfft(sclx,wfft1,wfft2,wfft3,wfft4,is)
  !$omp section   
      call yfft(scly,wfft1,wfft2,wfft3,wfft4,is)
  !$omp section   
      call yfft(sclz,wfft1,wfft2,wfft3,wfft4,is)
  !$omp end parallel sections
  end if

  if (ipolyflag .eq. 1 .and. it .ge. (src_start-1)) then 
  !$omp parallel sections default(shared) private(wfft1,wfft2,wfft3,wfft4)      
  !$omp section   
      call yfft(c11,wfft1,wfft2,wfft3,wfft4,is)    
  !$omp section   
      call yfft(c12,wfft1,wfft2,wfft3,wfft4,is)
  !$omp section   
      call yfft(c13,wfft1,wfft2,wfft3,wfft4,is)
  !$omp section   
      call yfft(c21,wfft1,wfft2,wfft3,wfft4,is)
  !$omp section   
      call yfft(c22,wfft1,wfft2,wfft3,wfft4,is)
  !$omp section   
      call yfft(c23,wfft1,wfft2,wfft3,wfft4,is)
  !$omp section   
      call yfft(c31,wfft1,wfft2,wfft3,wfft4,is)
  !$omp section   
      call yfft(c32,wfft1,wfft2,wfft3,wfft4,is)
  !$omp section   
      call yfft(c33,wfft1,wfft2,wfft3,wfft4,is)

  !$omp section   
      call yfft(dc111,wfft1,wfft2,wfft3,wfft4,is)    
  !$omp section   
      call yfft(dc112,wfft1,wfft2,wfft3,wfft4,is)
  !$omp section   
      call yfft(dc113,wfft1,wfft2,wfft3,wfft4,is)
  !$omp section   
      call yfft(dc211,wfft1,wfft2,wfft3,wfft4,is)
  !$omp section   
      call yfft(dc212,wfft1,wfft2,wfft3,wfft4,is)
  !$omp section   
      call yfft(dc213,wfft1,wfft2,wfft3,wfft4,is)
  !$omp section   
      call yfft(dc311,wfft1,wfft2,wfft3,wfft4,is)
  !$omp section   
      call yfft(dc312,wfft1,wfft2,wfft3,wfft4,is)
  !$omp section   
      call yfft(dc313,wfft1,wfft2,wfft3,wfft4,is)

  !$omp section   
      call yfft(dc121,wfft1,wfft2,wfft3,wfft4,is)    
  !$omp section   
      call yfft(dc122,wfft1,wfft2,wfft3,wfft4,is)
  !$omp section   
      call yfft(dc123,wfft1,wfft2,wfft3,wfft4,is)
  !$omp section   
      call yfft(dc221,wfft1,wfft2,wfft3,wfft4,is)
  !$omp section   
      call yfft(dc222,wfft1,wfft2,wfft3,wfft4,is)
  !$omp section   
      call yfft(dc223,wfft1,wfft2,wfft3,wfft4,is)
  !$omp section   
      call yfft(dc321,wfft1,wfft2,wfft3,wfft4,is)
  !$omp section   
      call yfft(dc322,wfft1,wfft2,wfft3,wfft4,is)
  !$omp section   
      call yfft(dc323,wfft1,wfft2,wfft3,wfft4,is)

  !$omp section   
      call yfft(dc131,wfft1,wfft2,wfft3,wfft4,is)    
  !$omp section   
      call yfft(dc132,wfft1,wfft2,wfft3,wfft4,is)
  !$omp section   
      call yfft(dc133,wfft1,wfft2,wfft3,wfft4,is)
  !$omp section   
      call yfft(dc231,wfft1,wfft2,wfft3,wfft4,is)
  !$omp section   
      call yfft(dc232,wfft1,wfft2,wfft3,wfft4,is)
  !$omp section   
      call yfft(dc233,wfft1,wfft2,wfft3,wfft4,is)
  !$omp section   
      call yfft(dc331,wfft1,wfft2,wfft3,wfft4,is)
  !$omp section   
      call yfft(dc332,wfft1,wfft2,wfft3,wfft4,is)
  !$omp section   
      call yfft(dc333,wfft1,wfft2,wfft3,wfft4,is)
  !$omp end parallel sections   
  end if

! End of extra ffts 
! ================================================================================ !     

! New vcw3d (includes vcw3dp stuff) - Ryan 6/28/22
      call vcw3d(u,v,w,omx,omy,omz,fn,gn,scalar,sclx,scly,sclz,scn,     &
                 c11,c12,c13,c21,c22,c23,c31,c32,c33,                   &
                 u11,u12,u13,u21,u22,u23,u31,u32,u33,Lu,Lv,Lw,          &
                 dc111,dc112,dc113,dc211,dc212,dc213,dc311,dc312,dc313, &
                 dc121,dc122,dc123,dc221,dc222,dc223,dc321,dc322,dc323, &
                 dc131,dc132,dc133,dc231,dc232,dc233,dc331,dc332,dc333, &
                 c11n,c12n,c13n,c22n,c23n,c33n,                         &
                 str11n,str12n,str13n,str22n,str23n,str33n,             &
                 qp11,qp12,qp13,qp22,qp23,qp33) 

  
  write(*,*) '1'
  !   transform all data into y-spectral
       
       is = -1
 
! Adding u and w yffts for GMU test because I think whatever is happening here
! is covered by the immersed boundary method - Ryan 7/7/22
       if (flow_select .lt. 0) then
          call yfft(u,wfft1,wfft2,wfft3,wfft4,is) 
          call yfft(w,wfft1,wfft2,wfft3,wfft4,is) 
       end if

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

! -------------------------------------------------------------------------------- !
!                   Adding the other ffts - Ryan 2-18-22                           !
! -------------------------------------------------------------------------------- !

  ! Leaving this since it's less than 8 ffts and we use at least 8 cores in
  ! parallel (i.e., this calculation does not slow down the code) - Ryan 7/13/22
  !$omp section     
      call yfft(scalar,wfft1,wfft2,wfft3,wfft4,is) !added scalar 1/15/2020
  !$omp section     
      call yfft(scn,wfft1,wfft2,wfft3,wfft4,is)
  !$omp end parallel sections

  if (ipolyflag .eq. 1 .and. it .ge. (src_start-1)) then
  !$omp parallel sections default(shared) private(wfft1,wfft2,wfft3,wfft4)   
  !$omp section     
      call yfft(c11,wfft1,wfft2,wfft3,wfft4,is)    
  !$omp section     
      call yfft(c12,wfft1,wfft2,wfft3,wfft4,is)
  !$omp section     
      call yfft(c13,wfft1,wfft2,wfft3,wfft4,is)
  !$omp section     
      call yfft(c21,wfft1,wfft2,wfft3,wfft4,is)
  !$omp section     
      call yfft(c22,wfft1,wfft2,wfft3,wfft4,is)
  !$omp section     
      call yfft(c23,wfft1,wfft2,wfft3,wfft4,is)
  !$omp section     
      call yfft(c31,wfft1,wfft2,wfft3,wfft4,is)
  !$omp section     
      call yfft(c32,wfft1,wfft2,wfft3,wfft4,is)
  !$omp section     
      call yfft(c33,wfft1,wfft2,wfft3,wfft4,is)


  !$omp section     
      call yfft(c11n,wfft1,wfft2,wfft3,wfft4,is)    
  !$omp section     
      call yfft(c12n,wfft1,wfft2,wfft3,wfft4,is)
  !$omp section     
      call yfft(c13n,wfft1,wfft2,wfft3,wfft4,is)
  !$omp section     
      call yfft(c22n,wfft1,wfft2,wfft3,wfft4,is)
  !$omp section     
      call yfft(c23n,wfft1,wfft2,wfft3,wfft4,is)
  !$omp section     
      call yfft(c33n,wfft1,wfft2,wfft3,wfft4,is)

  !$omp section     
      call yfft(str11n,wfft1,wfft2,wfft3,wfft4,is)    
  !$omp section     
      call yfft(str12n,wfft1,wfft2,wfft3,wfft4,is)
  !$omp section     
      call yfft(str13n,wfft1,wfft2,wfft3,wfft4,is)
  !$omp section     
      call yfft(str22n,wfft1,wfft2,wfft3,wfft4,is)
  !$omp section     
      call yfft(str23n,wfft1,wfft2,wfft3,wfft4,is)
  !$omp section     
      call yfft(str33n,wfft1,wfft2,wfft3,wfft4,is)

  !$omp section     
      call yfft(qp11,wfft1,wfft2,wfft3,wfft4,is) 
  !$omp section     
      call yfft(qp12,wfft1,wfft2,wfft3,wfft4,is)
  !$omp section     
      call yfft(qp13,wfft1,wfft2,wfft3,wfft4,is)
  !$omp section     
      call yfft(qp22,wfft1,wfft2,wfft3,wfft4,is)
  !$omp section     
      call yfft(qp23,wfft1,wfft2,wfft3,wfft4,is)
  !$omp section     
      call yfft(qp33,wfft1,wfft2,wfft3,wfft4,is)

! End of extra ffts -- Check to make sure omp sections are ok for this
! ================================================================================ !     
  !$omp end parallel sections
  end if

       call norm(v)
       call norm(omy)
       call norm(fn)
       call norm(gn)
       call norm(omz)

! ------------------------------------------------------------------------------- !
!               More norms - Ryan 2-18-22                                         !
! ------------------------------------------------------------------------------- !

     
       call norm(u)  !added 7.28.97
       call norm(w)  !added 7.28.97 
       call norm(scalar)
       call norm(scn)

       if (ipolyflag .eq. 1 .and. it .ge. (src_start-1)) then
       call norm(c11) !added by rah 11/9/2021
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

       call norm(qp11) !added on 1/15/2020
       call norm(qp12)
       call norm(qp13)
       call norm(qp22)
       call norm(qp23)
       call norm(qp33)
       end if

!  Add polymer forces to the momentum equation
!  First calculate forces in polyforce then 
!  add forces to momentum equation using subforce
!  changed to only have one subforce routine - Ryan 8/11/23
       if (ipolyflag .eq. 1 .and. it .ge. (src_start)) then
            call polyforce(qp11,qp12,qp13,qp22,qp23,qp33,t1,t2,t3)
            call norm(t1)
            call norm(t2)
            call norm(t3)
            call subforce(gn,fn,omz,t1,t2,t3)


         if (it .eq. (src_start-1)) then
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


!        Add Coriolis forces and buoyancy forces

       if (scl_flag .ge. 1) then
         call coriolis1(u,v,w,gn,fn,omz)
         call thermal1(fn,gn,scalar,scn,csource)
       end if

! ================================================================================ !     

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
                scn(k,j,i) = scn(k,j,i) * ysmth 
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
                scn(k,j,i) = scn(k,j,i) * zsmth 
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
! Extra norm - Ryan 2-18-22
       call norm(scn)
  

! Adding a time loop subroutine call to update polymer source location
     if (scl_flag .eq. 2 .and. it .ge. src_start .and. it .le. src_stop) then
       call src_update(psource)

!    Added 8/4/22
       call scram(psource,csource)
       call xyzfft(csource,w1,w2,wfft1,wfft2,wfft3,wfft4,-1)
     end if

! ------------------------------------------------------------------------ !

  !     update the running time
  
       time=time+dt
  write(*,*) '4',it,nsteps
       if(mod(it,iprnfrq) .eq. 0 .or. it .eq. nsteps) then
  
          open(122, file = 'outputs/last-restart', status = 'replace', form = 'unformatted')   ! alex 09/14/2018
          write(122) time
          write(122) v,fn,fnm1,omy,gn,gnm1,u0,h1n,h1nm1,w0,h3n,h3nm1
          write(122) fxintg, fyintg, fzintg
          rewind(122)
          close(122)

! Adding some file write that I'm not sure what is for yet
         if (ipolyflag .eq. 1) then
         open(123, file = 'outputs/c-last-restart', status = 'replace', form = 'unformatted') 
         write(123) time
         write(123) c11,c12,c13,c21,c22,c23,c31,c32,c33,c11n,c12n,c13n,c22n,c23n,c33n,     &
                    c11nm1,c12nm1,c13nm1,c22nm1,c23nm1,c33nm1,str11n,str12n,str13n,str22n, &
                    str23n,str33n,str11nm1,str12nm1,str13nm1,str22nm1,str23nm1,str33nm1
         rewind(123)
         close(123)
         end if

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

    ! Polymer inputs - Ryan 2-24-22
    real omearth,rotang,grav,virtual
    real zlmax,tpoly,qbeta,diff,diffpoly,deltat,c11z,c22z,c33z,alpha_poly,mpoly
    real ampbdy1,ampbdy2,c11amp,tfrac
    real att,btt,abb,bbb ! BCs for scalar ?
    real xshift,yshift,zshift,sigmax,sigmay,sigmaz
    real forbeta,polyrate

    integer src_start, src_stop
    integer ipolyflag,itarget,ipeter,scl_flag
    ! -------------------------------- !
  
    integer irstrt,nsteps,iprnfrq,nt,nb,ixfax(19),izfax(19),iyfax(19),ixfax32(19),izfax32(19)
    integer ib,it,i,j,k,iyt,iyb,is,jj, imatrix(nyp,mz,mx), bfhead, bftail, bfwidth, nyi, nzi, nxi
    integer print3d, crstrt
    common/preforce/ seght
    common/pre2forc/ spread,fxintg
    common/pre4f/ fzintg, fyintg
    common/pre5f/ it
    common/params/ gain, ugain
    common/data1/ pi,dt,theta,wavz,wavx,c,yc
    common/data2/ re,alpha,beta,xl,zl,yl,time,dyde,xstart,uinf
! Adding common blocks - Ryan 2-24-22
    common/thermcor/ omearth,rotang,grav,virtual
    common/polymer1/ zlmax,tpoly,qbeta,diff,diffpoly,deltat,c11z,c22z,c33z,alpha_poly,mpoly
    common/polymer2/ ampbdy1,ampbdy2,c11amp,tfrac,polyrate,src_start,src_stop
    common/polymer3/ ipolyflag,itarget,ipeter,scl_flag
    common/sclgeom/ xshift,yshift,zshift,sigmax,sigmay,sigmaz
    common/data22/ forbeta
! Done adding common blocks
    common/energ/ e1,ckx,tm1
    common/iocntrl/ irstrt,nsteps,iprnfrq,print3d,crstrt
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

    real bdyfx, L, rad, xcenter, ycenter, zcenter 
    common/vortRing/ bdyfx, L, rad, xcenter, ycenter, zcenter 
  
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
    read(110,*) omearth
    read(110,*) rotang
    read(110,*) grav
    read(110,*) virtual

    read(110,*) 
    read(110,*) 
    read(110,*) 
    read(110,*) 
    read(110,*) 

    ! Read in selection flags/switches 
    read(110,*) irstrt    !--- irstrt -- ne 0 means a restart file is used to initialize this run
    read(110,*) crstrt    !--- crstrt -- ne 0 means a restart file is used to initialize this run
    read(110,*) print3d   !--- print3d -- switch how we print the flowfield
    read(110,*) geomtype
    read(110,*) readvdes
    read(110,*) particle_flag
    read(110,*) flow_select
    read(110,*) itarget
    read(110,*) ipeter
    read(110,*) ipolyflag
    read(110,*) scl_flag

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

    gravity = grav ! Ryan 7/14/22 - There's probably a better way to do this, but this was easiest to have the same gravity in two pre-defined common blocks

    read(110,*) 
    read(110,*) 
    read(110,*) 
    read(110,*) 
    read(110,*) 

    ! Read in data for scalar and polymer
    read(110,*) zlmax
    read(110,*) tpoly
    read(110,*) alpha_poly
    read(110,*) diff
    read(110,*) diffpoly
    read(110,*) deltat
    read(110,*) c11z
    read(110,*) c22z
    read(110,*) c33z
    read(110,*) qbeta
    read(110,*) mpoly
    read(110,*) xshift,yshift,zshift
    read(110,*) sigmax,sigmay,sigmaz
    read(110,*) polyrate
    read(110,*) src_start 
    read(110,*) src_stop

    if (src_stop .lt. src_start) then
        src_stop = nsteps
    end if

    read(110,*) 
    read(110,*) 
    read(110,*) 
    read(110,*) 
    read(110,*) 

    ! Read in data for body force terms (used for vortex ring generation
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

        case (0,5,10) ! Vortex Ring or Still fluid
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
            print *,'       0,1,2,3,4,5,11,12'
            print *,''
    end select

    ! For all the cases so far, these numbers have not changed, so I'm defining
    ! them outside the selection loop
 
    p1b = 1.0
    p1t = -1.0
    p2b = -1.0
    p2t = 0.5 

end subroutine setbcs


! Adding new arguments to the end of existing ones (starting with scalar)- Ryan 2-24-22
subroutine initial(u,u0,v,w,w0,omx,omy,omz,fn,fnm1,gn,gnm1,h1n,h1nm1,h3n,h3nm1,  &
     wrk1,wrkc,wfft1,wfft2,wfft3,wfft4,w1,w2,bctop,bcbot,a,t,                    &
     scalar,sclx,scly,sclz,scn,scnm1,T1,T2,T3,c11,c12,c13,c21,c22,c23,c31,c32,c33, &
     u11,u12,u13,u21,u22,u23,u31,u32,u33,dc111,dc112,dc113,dc211,dc212,dc213,    &
     dc311,dc312,dc313,dc121,dc122,dc123,dc221,dc222,dc223,dc321,dc322,dc323,    &
     dc131,dc132,dc133,dc231,dc232,dc233,dc331,dc332,dc333,c11n,c12n,c13n,       &
     c22n,c23n,c33n,c11nm1,c12nm1,c13nm1,c22nm1,c23nm1,c33nm1,str11n,str12n,     &
     str13n,str22n,str23n,str33n,str11nm1,str12nm1,str13nm1,str22nm1,str23nm1,   &
     str33nm1,qp11,qp12,qp13,qp22,qp23,qp33,csource,Lu,Lv,Lw)

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

! ============================================================================== !
!                  Adding variables from polymer code - Ryan 3-2-22              !
! ============================================================================== !


  complex scalar(nyp,nz,nxh),sclz(nyp,nz,nxh)
  complex sclx(nyp,nz,nxh),scly(nyp,nz,nxh) ! These are GF1 and GF2 in GMU code -Ryan 6-15-22
  complex scn(nyp,nz,nxh), scnm1(nyp,nz,nxh)
  complex csource(nyp,nz,nxh), psource(nyp,nz,nx) ! added 8/3/22

  complex c11n(nyp,nz,nxh),c12n(nyp,nz,nxh),c13n(nyp,nz,nxh)
  complex c22n(nyp,nz,nxh),c23n(nyp,nz,nxh),c33n(nyp,nz,nxh) ! nonlinear terms for conformation
  complex c11nm1(nyp,nz,nxh),c12nm1(nyp,nz,nxh),c13nm1(nyp,nz,nxh)
  complex c22nm1(nyp,nz,nxh),c23nm1(nyp,nz,nxh),c33nm1(nyp,nz,nxh) ! nonlinear terms for conformation
  complex str11n(nyp,nz,nxh),str12n(nyp,nz,nxh),str13n(nyp,nz,nxh)
  complex str22n(nyp,nz,nxh),str23n(nyp,nz,nxh),str33n(nyp,nz,nxh) ! nonlinear terms for conformation
  complex str11nm1(nyp,nz,nxh),str12nm1(nyp,nz,nxh),str13nm1(nyp,nz,nxh)
  complex str22nm1(nyp,nz,nxh),str23nm1(nyp,nz,nxh),str33nm1(nyp,nz,nxh) ! nonlinear terms for conformation

!     conformation tensor
  complex c11(nyp,nz,nxh),c12(nyp,nz,nxh),c13(nyp,nz,nxh)
  complex c21(nyp,nz,nxh),c22(nyp,nz,nxh),c23(nyp,nz,nxh)
  complex c31(nyp,nz,nxh),c32(nyp,nz,nxh),c33(nyp,nz,nxh)      
!     velocity gradient tensor 
  complex u11(nyp,nz,nxh),u12(nyp,nz,nxh),u13(nyp,nz,nxh)
  complex u21(nyp,nz,nxh),u22(nyp,nz,nxh),u23(nyp,nz,nxh) 
  complex u31(nyp,nz,nxh),u32(nyp,nz,nxh),u33(nyp,nz,nxh)
  complex,dimension(nyp,nz,nxh) :: Lu,Lv,Lw ! Laplacian terms - Ryan 7/24/23
!     derivatives of conformation tensor
  complex dc111(nyp,nz,nxh),dc112(nyp,nz,nxh),dc113(nyp,nz,nxh) 
  complex dc211(nyp,nz,nxh),dc212(nyp,nz,nxh),dc213(nyp,nz,nxh)
  complex dc311(nyp,nz,nxh),dc312(nyp,nz,nxh),dc313(nyp,nz,nxh)  

  complex dc121(nyp,nz,nxh),dc122(nyp,nz,nxh),dc123(nyp,nz,nxh)
  complex dc221(nyp,nz,nxh),dc222(nyp,nz,nxh),dc223(nyp,nz,nxh)
  complex dc321(nyp,nz,nxh),dc322(nyp,nz,nxh),dc323(nyp,nz,nxh)

  complex dc131(nyp,nz,nxh),dc132(nyp,nz,nxh),dc133(nyp,nz,nxh)  
  complex dc231(nyp,nz,nxh),dc232(nyp,nz,nxh),dc233(nyp,nz,nxh)  
  complex dc331(nyp,nz,nxh),dc332(nyp,nz,nxh),dc333(nyp,nz,nxh)
!     polymer stresses used for targeting
  complex qp11(nyp,nz,nxh),qp12(nyp,nz,nxh),qp13(nyp,nz,nxh)
  complex qp22(nyp,nz,nxh),qp23(nyp,nz,nxh)
  complex qp33(nyp,nz,nxh)

  complex forcebdy1(nyp,nz,nxh),forcebdy2(nyp,nz,nxh)
  complex tempforce(nyp,nz,nxh) ! second force from bforce1

  real wrk11(nyp,nz,nx),wrk12(nyp,nz,nx),wrk13(nyp,nz,nx)
  real wrk21(nyp,nz,nx),wrk22(nyp,nz,nx),wrk23(nyp,nz,nx)
  real wrk31(nyp,nz,nx),wrk32(nyp,nz,nx),wrk33(nyp,nz,nx)

! Done adding variables
! ============================================================================== !

  common/data1/ pi,dt,theta,wavz(nz),wavx(nxh),c(nyp),yc(nyp)
  common/data2/ re,alpha,beta,xl,zl,yl,time,dyde,xstart,uinf
  common/energ/ e1(nz,nxh),ckx(nz,nxh),tm1
  common/grnfct/nt,nb,p1t,p1b,p2t,p2b
  common/iocntrl/irstrt,nsteps,iprnfrq,print3d,crstrt
  common/pre2forc/ spread(mz2,mx2),fxintg(nyp,mz,mx)
  common/pre4f/ fzintg(nyp,mz,mx),fyintg(nyp,mz,mx)
  common/pre5f/ it
  common/params/ gain, ugain

! Polymer common blocks - Ryan 3-2-22
  common/fbdy/ forcebdy1, forcebdy2 
  common/thermcor/ omearth,rotang,grav,virtual
  common/polymer3/ ipolyflag,itarget,scl_flag

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
     
! Adding arguments from polymer code - Ryan 3-2-22 (maybe?)
     call readuv(omx,omy,omz,sclx,psource,wrk11,wrk12,wrk13,wrk21,wrk22,wrk23,wrk31,wrk32,wrk33)


     call scram(omy,v)
     call xyzfft(v,w1,w2,wfft1,wfft2,wfft3,wfft4,-1)

     call scram(omx,u)
     call xyzfft(u,w1,w2,wfft1,wfft2,wfft3,wfft4,-1)

     call scram(omz,w)
     call xyzfft(w,w1,w2,wfft1,wfft2,wfft3,wfft4,-1)

!    Added by Ryan 3-2-22
     call scram(sclx,scalar)
     call xyzfft(scalar,w1,w2,wfft1,wfft2,wfft3,wfft4,-1)

!    Added 8/3/22
     call scram(psource,csource)
     call xyzfft(csource,w1,w2,wfft1,wfft2,wfft3,wfft4,-1)

!    Adding applied body force (?) subroutine calls - Ryan 3-8-22
!----------------------------------------------------------------
! I'm not sure if this should be in the immersed boundary force section for this code
     if (flow_select .lt. 0) then
         call bforce1(omy,tempforce)
         call scram(omy,forcebdy1)
         call xyzfft(forcebdy1,w1,w2,wfft1,wfft2,wfft3,wfft4,-1)
         call scram(tempforce,forcebdy2)
         call xyzfft(forcebdy2,w1,w2,wfft1,wfft2,wfft3,wfft4,-1)
     end if


     call scram(wrk11,c11)
     call xyzfft(c11,w1,w2,wfft1,wfft2,wfft3,wfft4,-1)
     call scram(wrk12,c12)
     call xyzfft(c12,w1,w2,wfft1,wfft2,wfft3,wfft4,-1)
     call scram(wrk13,c13)
     call xyzfft(c13,w1,w2,wfft1,wfft2,wfft3,wfft4,-1)
     call scram(wrk21,c21)
     call xyzfft(c21,w1,w2,wfft1,wfft2,wfft3,wfft4,-1)
     call scram(wrk22,c22)
     call xyzfft(c22,w1,w2,wfft1,wfft2,wfft3,wfft4,-1)
     call scram(wrk23,c23)
     call xyzfft(c23,w1,w2,wfft1,wfft2,wfft3,wfft4,-1)
     call scram(wrk31,c31)
     call xyzfft(c31,w1,w2,wfft1,wfft2,wfft3,wfft4,-1)
     call scram(wrk32,c32)
     call xyzfft(c32,w1,w2,wfft1,wfft2,wfft3,wfft4,-1)
     call scram(wrk33,c33)
     call xyzfft(c33,w1,w2,wfft1,wfft2,wfft3,wfft4,-1)
!--------------------------------------------------------
! Done adding new subroutine calls

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

!    More norms - Ryan 3-8-22

     call norm(scalar)
     call norm(forcebdy1)
     call norm(forcebdy2)
     call norm(csource)

     call norm(c11)
     call norm(c12)
     call norm(c13)
     call norm(c21)
     call norm(c22)
     call norm(c23)
     call norm(c31)
     call norm(c32)
     call norm(c33)
! ------------------------------ !
  
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
  
     call vort(u,v,w,scalar,omx,omz,sclx,scly,sclz,wrkc) ! Extra args - Ryan 6/27/22

     call norm(omx)
     call norm(omz)

! ========================================================== !
!    Adding more norms and conformation tensor derivatives   !
! ========================================================== !
! Ryan 3-8-22

     call norm(sclx)
     call norm(scly)
     call norm(sclz)

! Calculate derivates needed for conformation tensor sources

     call derivscji(u,v,w,c11,c12,c13,c21,c22,c23,c31,c32,c33,   &
          u11,u12,u13,u21,u22,u23,u31,u32,u33,Lu,Lv,Lw,          &
          dc111,dc112,dc113,dc211,dc212,dc213,dc311,dc312,dc313, &
          dc121,dc122,dc123,dc221,dc222,dc223,dc321,dc322,dc323, &
          dc131,dc132,dc133,dc231,dc232,dc233,dc331,dc332,dc333)

! Need to normalize all arrays 

     call norm(c11)
     call norm(c12)
     call norm(c13)
     call norm(c21)
     call norm(c22)
     call norm(c23)
     call norm(c31)
     call norm(c32)
     call norm(c33)

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
! Done adding subroutine calls
! ===================================================== !


!    transform all data into y-physical
  
     is = 1

     call yfft(u,wfft1,wfft2,wfft3,wfft4,is)
     call yfft(v,wfft1,wfft2,wfft3,wfft4,is)
     call yfft(w,wfft1,wfft2,wfft3,wfft4,is)
     call yfft(omx,wfft1,wfft2,wfft3,wfft4,is)
     call yfft(omy,wfft1,wfft2,wfft3,wfft4,is)
     call yfft(omz,wfft1,wfft2,wfft3,wfft4,is)

! ===================================================== !
!           A lot more FFTs - Ryan 3-8-22               !
! ===================================================== !
     call yfft(scalar,wfft1,wfft2,wfft3,wfft4,is) !added scalar 1/15/2020
     call yfft(sclx,wfft1,wfft2,wfft3,wfft4,is)     !  dc/dx
     call yfft(scly,wfft1,wfft2,wfft3,wfft4,is)     !  dc/dy
     call yfft(sclz,wfft1,wfft2,wfft3,wfft4,is)  !  dc/dz
     call yfft(c11,wfft1,wfft2,wfft3,wfft4,is)    
     call yfft(c12,wfft1,wfft2,wfft3,wfft4,is)
     call yfft(c13,wfft1,wfft2,wfft3,wfft4,is)
     call yfft(c21,wfft1,wfft2,wfft3,wfft4,is)
     call yfft(c22,wfft1,wfft2,wfft3,wfft4,is)
     call yfft(c23,wfft1,wfft2,wfft3,wfft4,is)
     call yfft(c31,wfft1,wfft2,wfft3,wfft4,is)
     call yfft(c32,wfft1,wfft2,wfft3,wfft4,is)
     call yfft(c33,wfft1,wfft2,wfft3,wfft4,is)

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

     call yfft(dc111,wfft1,wfft2,wfft3,wfft4,is)    
     call yfft(dc112,wfft1,wfft2,wfft3,wfft4,is)
     call yfft(dc113,wfft1,wfft2,wfft3,wfft4,is)
     call yfft(dc211,wfft1,wfft2,wfft3,wfft4,is)
     call yfft(dc212,wfft1,wfft2,wfft3,wfft4,is)
     call yfft(dc213,wfft1,wfft2,wfft3,wfft4,is)
     call yfft(dc311,wfft1,wfft2,wfft3,wfft4,is)
     call yfft(dc312,wfft1,wfft2,wfft3,wfft4,is)
     call yfft(dc313,wfft1,wfft2,wfft3,wfft4,is)

     call yfft(dc121,wfft1,wfft2,wfft3,wfft4,is)    
     call yfft(dc122,wfft1,wfft2,wfft3,wfft4,is)
     call yfft(dc123,wfft1,wfft2,wfft3,wfft4,is)
     call yfft(dc221,wfft1,wfft2,wfft3,wfft4,is)
     call yfft(dc222,wfft1,wfft2,wfft3,wfft4,is)
     call yfft(dc223,wfft1,wfft2,wfft3,wfft4,is)
     call yfft(dc321,wfft1,wfft2,wfft3,wfft4,is)
     call yfft(dc322,wfft1,wfft2,wfft3,wfft4,is)
     call yfft(dc323,wfft1,wfft2,wfft3,wfft4,is)

     call yfft(dc131,wfft1,wfft2,wfft3,wfft4,is)    
     call yfft(dc132,wfft1,wfft2,wfft3,wfft4,is)
     call yfft(dc133,wfft1,wfft2,wfft3,wfft4,is)
     call yfft(dc231,wfft1,wfft2,wfft3,wfft4,is)
     call yfft(dc232,wfft1,wfft2,wfft3,wfft4,is)
     call yfft(dc233,wfft1,wfft2,wfft3,wfft4,is)
     call yfft(dc331,wfft1,wfft2,wfft3,wfft4,is)
     call yfft(dc332,wfft1,wfft2,wfft3,wfft4,is)
     call yfft(dc333,wfft1,wfft2,wfft3,wfft4,is)
! ===================================================== !


! New vcw3d (includes vcw3dp stuff) - Ryan 6/28/22
      call vcw3d(u,v,w,omx,omy,omz,fn,gn,scalar,sclx,scly,sclz,scn,     &
                 c11,c12,c13,c21,c22,c23,c31,c32,c33,                   &
                 u11,u12,u13,u21,u22,u23,u31,u32,u33,Lu,Lv,Lw,          &
                 dc111,dc112,dc113,dc211,dc212,dc213,dc311,dc312,dc313, &
                 dc121,dc122,dc123,dc221,dc222,dc223,dc321,dc322,dc323, &
                 dc131,dc132,dc133,dc231,dc232,dc233,dc331,dc332,dc333, &
                 c11n,c12n,c13n,c22n,c23n,c33n,                         &
                 str11n,str12n,str13n,str22n,str23n,str33n,             &
                 qp11,qp12,qp13,qp22,qp23,qp33) 

!    transform all data into y-spectral
  
     is = -1

     call yfft(v,wfft1,wfft2,wfft3,wfft4,is)
     call yfft(omy,wfft1,wfft2,wfft3,wfft4,is)
     call yfft(fn,wfft1,wfft2,wfft3,wfft4,is)
     call yfft(gn,wfft1,wfft2,wfft3,wfft4,is)
     call yfft(omz,wfft1,wfft2,wfft3,wfft4,is)

! ===================================================== !
!           A lot more FFTs - Ryan 3-8-22               !
! ===================================================== !
     call yfft(scalar,wfft1,wfft2,wfft3,wfft4,is) !added scalar 1/15/2020
     call yfft(scn,wfft1,wfft2,wfft3,wfft4,is)

     call yfft(c11,wfft1,wfft2,wfft3,wfft4,is)    
     call yfft(c12,wfft1,wfft2,wfft3,wfft4,is)
     call yfft(c13,wfft1,wfft2,wfft3,wfft4,is)
     call yfft(c21,wfft1,wfft2,wfft3,wfft4,is)
     call yfft(c22,wfft1,wfft2,wfft3,wfft4,is)
     call yfft(c23,wfft1,wfft2,wfft3,wfft4,is)
     call yfft(c31,wfft1,wfft2,wfft3,wfft4,is)
     call yfft(c32,wfft1,wfft2,wfft3,wfft4,is)
     call yfft(c33,wfft1,wfft2,wfft3,wfft4,is)

     call yfft(c11n,wfft1,wfft2,wfft3,wfft4,is)    
     call yfft(c12n,wfft1,wfft2,wfft3,wfft4,is)
     call yfft(c13n,wfft1,wfft2,wfft3,wfft4,is)
     call yfft(c22n,wfft1,wfft2,wfft3,wfft4,is)
     call yfft(c23n,wfft1,wfft2,wfft3,wfft4,is)
     call yfft(c33n,wfft1,wfft2,wfft3,wfft4,is)

     call yfft(str11n,wfft1,wfft2,wfft3,wfft4,is)    
     call yfft(str12n,wfft1,wfft2,wfft3,wfft4,is)
     call yfft(str13n,wfft1,wfft2,wfft3,wfft4,is)
     call yfft(str22n,wfft1,wfft2,wfft3,wfft4,is)
     call yfft(str23n,wfft1,wfft2,wfft3,wfft4,is)
     call yfft(str33n,wfft1,wfft2,wfft3,wfft4,is)

     call yfft(qp11,wfft1,wfft2,wfft3,wfft4,is) !added on 1/15/2020
     call yfft(qp12,wfft1,wfft2,wfft3,wfft4,is)
     call yfft(qp13,wfft1,wfft2,wfft3,wfft4,is)
     call yfft(qp22,wfft1,wfft2,wfft3,wfft4,is)
     call yfft(qp23,wfft1,wfft2,wfft3,wfft4,is)
     call yfft(qp33,wfft1,wfft2,wfft3,wfft4,is)

! ===================================================== !


! ===================================================== !
!           A lot more norms - Ryan 3-8-22              !
! ===================================================== !
     call norm(scalar)
     call norm(scn)

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

     call norm(qp11) !added on 1/15/2020
     call norm(qp12)
     call norm(qp13)
     call norm(qp22)
     call norm(qp23)
     call norm(qp33)

! ===================================================== !


     call norm(v)
     call norm(omy)
     call norm(fn)
     call norm(gn)
     call norm(omz)

! ===================================================== !
!           Polymer forces - Ryan 3-8-22                !
! ===================================================== !
! Add polymer forces to the momentum equation
! First calculate forces in polyforce then 
! add forces to momentum equation using subforce2
! changed to only have one subforce routine - Ryan 8/11/23
    
    if (ipolyflag .eq. 1 .and. it .ge. (src_start)) then
         call polyforce(qp11,qp12,qp13,qp22,qp23,qp33,t1,t2,t3)
         call norm(t1)
         call norm(t2)
         call norm(t3)
         call subforce(gn,fn,omz,t1,t2,t3)
    end if


!   add coriolis forces and buoyancy forces

    call coriolis1(u,v,w,gn,fn,omz)
    call thermal1(fn,gn,scalar,scn,csource)
! ===================================================== !

     do i=1,nyp
        h1n(i)=real(gn(i,1,1))
        h3n(i)=real(omz(i,1,1))
     end do

     call nonlin(fn,omz,gn,wrkc,wrk1)

     call norm(fn)
     call norm(gn)
     call norm(scn)

!    for the first time step fnm1=fn & h1nm1=h1n

     do k=1,nxh
        do j=1,nz
           do i=1,nyp
              fnm1(i,j,k)=fn(i,j,k)
              gnm1(i,j,k)=gn(i,j,k)
              scnm1(i,j,k)=scn(i,j,k)
           end do
        end do
     end do
     do i=1,nyp
        h1nm1(i)=h1n(i)
        h3nm1(i)=h3n(i)
     end do

! More stuff - Ryan 3-8-22
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
             enddo
         enddo
     enddo
! Done with more stuff
  else

! Adding body force loops - Ryan 3-8-22
  do k=1,nxh
    do j=1,nz
          do i=1,nyp
              forcebdy1(i,j,k)=(0.0,0.0)  !makes body force zero after restart
              forcebdy2(i,j,k)=(0.0,0.0)
              csource(i,j,k) = (0.0,0.0) ! added 8/3/22
          enddo
      enddo
  enddo

!    read restart file
     open(112, file = 'setup/restart', status = 'old', form = 'unformatted')   ! alex 09/14/2018
     read(112) time
     read(112) v,fn,fnm1,omy,gn,gnm1,u0,h1n,h1nm1,w0,h3n,h3nm1
     read(112) fxintg, fyintg, fzintg
     close(112)                                           ! alex 09/14/2018

!    read polymer restart file
     if (crstrt .eq. 1 .and. ipolyflag .eq. 1) then
         open(113, file = 'setup/c-restart', status = 'replace', form = 'unformatted') 
         read(113) time
         read(113) c11,c12,c13,c21,c22,c23,c31,c32,c33,c11n,c12n,c13n,c22n,c23n,c33n,     &
                    c11nm1,c12nm1,c13nm1,c22nm1,c23nm1,c33nm1,str11n,str12n,str13n,str22n, &
                    str23n,str33n,str11nm1,str12nm1,str13nm1,str22nm1,str23nm1,str33nm1
         close(113)
     end if
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
  real zlmax,tpoly,qbeta,diff,diffpoly,deltat,c11z,c22z,c33z
  common/data1/ pi,dt,theta,wavz(nz),wavx(nxh),c(nyp),yc(nyp)
  common/data2/ re,alpha,beta,xl,zl,yl,time,dyde,xstart,uinf
  common/polymer1/ zlmax,tpoly,qbeta,diff,diffpoly,deltat,c11z,c22z,c33z

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
  real zlmax,tpoly,qbeta,diff,diffpoly,deltat,c11z,c22z,c33z
  common/data1/ pi,dt,theta,wavz(nz),wavx(nxh),c(nyp),yc(nyp)
  common/data2/ re,alpha,beta,xl,zl,yl,time,dyde,xstart,uinf
  common/polymer1/ zlmax,tpoly,qbeta,diff,diffpoly,deltat,c11z,c22z,c33z

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

subroutine vort(u,v,w,scalar,omx,omz,sclx,scly,sclz,wrkc)
  use grid_size
  complex u(nyp,nz,nxh), v(nyp,nz,nxh), w(nyp,nz,nxh)
  complex omx(nyp,nz,nxh), omz(nyp,nz,nxh), wrkc(nyp,nz,nxh)
  complex scalar(nyp,nz,nxh), sclx(nyp,nz,nxh), scly(nyp,nz,nxh), sclz(nyp,nz,nxh)
  complex im
  common/data1/ pi,dt,theta,wavz(nz),wavx(nxh),c(nyp),yc(nyp)

!**********************************************************************
!     omx =  dw/dy-dv/dz                                              *
!     omy =  du/dz-dw/dx   calculated above                           *
!     omz =  dv/dx-du/dy                                              *
!    sclx =  dc/dx                                                    *
!    scly =  dc/dy                                                    *
!  scalar =  dc/dz                                                    *
!     all quantities are spectral                                     *
!     scalar contains passive scalar on input and dc/dz on output     *
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
 
!---------------------------------------------------------------------!
!              Adding in scalar stuff - Ryan 6/24/22                  ! 
!---------------------------------------------------------------------!
  ! dc/dx
  !$omp parallel do
  do k = 1,nxh
    do j = 1,nz
      do i = 1,nyp
        sclx(i,j,k) = im*wavx(k)*scalar(i,j,k)
      end do
    end do
  end do
  !$omp end parallel do

  ! dc/dy
  call cderiv(scalar,scly)

  ! dc/dz
  !$omp parallel do
  do k = 1,nxh
    do j = 1,nz
      do i = 1,nyp
        sclz(i,j,k) = im*wavz(j)*scalar(i,j,k)
      end do
    end do
  end do
  !$omp end parallel do

    ! Note: I don't see why these can't all be lumped into a single nested do
    ! loop, but I'm not going to worry about efficiency until I actually get it
    ! working - Ryan 6/24/22
!---------------------------------------------------------------------!

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

subroutine readuv(u1,u2,u3,sc,psource,c11,c12,c13,c21,c22,c23,c31,c32,c33)
  use grid_size
  implicit none

  integer, parameter :: nzlowres = 4, mzlowres = 3*nzlowres/2

  real,dimension(nyp,nz,nx) :: u1,u2,u3,initu,initv,initw,sc! added sc - Ryan 3-2-22
  real psource(nyp,nz,nx) ! added 8/3/22
  real re,alpha,beta,xl,zl,yl,time,dyde,xstart,uinf
  integer i, j, k, irstrt, nsteps, iprnfrq, numjstrips, jstrip, jcount, n, ns

  common /init/ initu, initv, initw
  common/iocntrl/ irstrt, nsteps, iprnfrq
  common/data2/ re,alpha,beta,xl,zl,yl,time,dyde,xstart,uinf

  real u_lowres(nyp,nzlowres,nx), v_lowres(nyp,nzlowres,nx), w_lowres(nyp,nzlowres,nx)
  real pi, xc1, yc1, zc1, xcor, zcor, xsq, ysq, zsq ! added more variable declarations - Ryan 3-16-22
  real xshift,yshift,zshift,sigmax,sigmay,sigmaz,betax,betay,betaz
  real zlmax,tpoly,qbeta,diff,diffpoly,deltat,c11z,c22z,c33z
  real ampbdy1,ampbdy2,c11amp,tfrac,polyrate
  integer ipolyflag,itarget,ipeter,scl_flag,src_start,src_stop

  common/polymer1/ zlmax,tpoly,qbeta,diff,diffpoly,deltat,c11z,c22z,c33z
  common/polymer2/ ampbdy1,ampbdy2,c11amp,tfrac,polyrate,src_start,src_stop
  common/polymer3/ ipolyflag,itarget,ipeter,scl_flag
  common/sclgeom/ xshift,yshift,zshift,sigmax,sigmay,sigmaz

  ! Conformation tensor - Ryan 3-2-22
  real c11(nyp,nz,nx),c12(nyp,nz,nx),c13(nyp,nz,nx)
  real c21(nyp,nz,nx),c22(nyp,nz,nx),c23(nyp,nz,nx)
  real c31(nyp,nz,nx),c32(nyp,nz,nx),c33(nyp,nz,nx)

  pi = 2.0 * acos(0.0)

  if (npart .ne. 0 .and. scl_flag .ne. 0) then
    ns = npart
  else
    ns = 1
  end if 

  if (irstrt .eq. 0) then
     ! Adding scalar initialization - Ryan 3-16-22
!     call c11init(c11) ! initialize c11 - this is only for the specific body
!     force, so I'm just leaving it as c11z for now. Generalize later if you
!     want to - Ryan 7/14/22

    ! zero out psource first
    do i = 1,nyp
      do j = 1,nz
        do k = 1,nx
          psource(i,j,k) = 0.0
        end do
      end do
    end do

    do n = 1,ns
    ! New source center definitions
    if (scl_flag .eq. 1) then
        yc1 = ycoord(nyh+1) + yshift
        zc1 = zl/2.0 + zshift
        xc1 = xl/2.0 + xshift
    else if (scl_flag .eq. 2) then
        open(95, file = 'setup/particles/particles.dat', status = 'old', action = 'read')
        do j = 1,n
            read(95,*) 
        end do
        read(95,*) xc1, yc1, zc1
        close(95)
    end if

    do  k = 1,nyp
       ysq = (ycoord(k) - yc1)**2
       betay = ysq/(2.*sigmay**2) 
       do  j = 1,nz       
          zcor = zl*(float(j-1)/(float(nz)))
          zsq = (zcor - zc1)**2
          betaz = zsq/(2.*sigmaz**2)
          do  i = 1,nx
             xcor = xl*(float(i-1)/(float(nx)))
             xsq = (xcor - xc1)**2
             betax = xsq/(2.*sigmax**2)

            ! Note (k,j,i) order!
             u1(k,j,i) = initu(k,j,i)
             u2(k,j,i) = initv(k,j,i)
             u3(k,j,i) = initw(k,j,i)

             ! ICs for Cij
             c11(k,j,i) = c11z ! Added by Ryan 7/14/22
             c12(k,j,i) = 0.0
             c13(k,j,i) = 0.0
             c21(k,j,i) = 0.0
             c22(k,j,i) = c22z
             c23(k,j,i) = 0.0
             c31(k,j,i) = 0.0
             c32(k,j,i) = 0.0
             c33(k,j,i) = c33z
             if (scl_flag .eq. 1) then
                 sc(k,j,i) = deltat*exp(-(betax + betay + betaz))
             else if (scl_flag .eq. 2) then
                 psource(k,j,i) = psource(k,j,i) + polyrate*exp(-(betax + betay + betaz))
             end if
          end do
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

! New subroutine for updating source location each time step - added by Ryan
! 8/4/22
subroutine src_update(psource)
  use grid_size
  use omp_lib

  implicit none

  real psource(nyp,nz,nx) ! added 8/3/22
  real re,alpha,beta,xl,zl,yl,time,dyde,xstart,uinf
  integer i, j, k, irstrt, nsteps, iprnfrq,n,ns

  common/iocntrl/ irstrt, nsteps, iprnfrq
  common/data2/ re,alpha,beta,xl,zl,yl,time,dyde,xstart,uinf

  real pi, xc1, yc1, zc1, xcor, zcor, xsq, ysq, zsq ! added more variable declarations - Ryan 3-16-22
  real xshift,yshift,zshift,sigmax,sigmay,sigmaz,betax,betay,betaz
  real ampbdy1,ampbdy2,c11amp,tfrac,polyrate
  integer ipolyflag,itarget,ipeter,scl_flag,src_start,src_stop
  real,dimension(npart) :: xpart,ypart,zpart
  
  common/polymer2/ ampbdy1,ampbdy2,c11amp,tfrac,polyrate,src_start,src_stop
  common/polymer3/ ipolyflag,itarget,ipeter,scl_flag
  common/sclgeom/ xshift,yshift,zshift,sigmax,sigmay,sigmaz
  common/part_traj/ xpart,ypart,zpart


  pi = 2.0 * acos(0.0)

  ! zero out psource first
  !$omp parallel do
  do i = 1,nyp
    do j = 1,nz
      do k = 1,nx
        psource(i,j,k) = 0.0
      end do
    end do
  end do
  !$omp end parallel do

  if (npart .ne. 0 .and. scl_flag .ne. 0) then
    ns = npart
  else
    ns = 1
  end if 

  !$omp parallel do default(shared) private(i,j,k,n,xc1,yc1,zc1,ysq,betay,zcor,zsq,betaz,xcor,xsq,betax)
  do n = 1,ns
    yc1 = ypart(n)
    zc1 = zpart(n)
    xc1 = xpart(n)
  do  k = 1,nyp
     ysq = (ycoord(k) - yc1)**2
     betay = ysq/(2.*sigmay**2) 
     do  j = 1,nz       
        zcor = zl*(float(j-1)/(float(nz)))
        zsq = (zcor - zc1)**2
        betaz = zsq/(2.*sigmaz**2)
        do  i = 1,nx
           xcor = xl*(float(i-1)/(float(nx)))
           xsq = (xcor - xc1)**2
           betax = xsq/(2.*sigmax**2)

           psource(k,j,i) = psource(k,j,i) +  polyrate*exp(-(betax + betay + betaz))
        end do
     end do
  end do
  end do
  !$omp end parallel do


end subroutine src_update

!=======================================================================


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

subroutine output(u,v,w,urms,vrms,wrms,umn,wrk1,wrkc,wfft1,wfft2,wfft3,wfft4,w1,w2,scalar)

  use grid_size
  implicit none 

  complex u(nyp,nz,nxh), v(nyp,nz,nxh), w(nyp,nz,nxh)
  complex wrkc(nyp,nz,nxh)
  complex im
  complex scalar(nyp,nz,nxh)

  real wrk1(nyp,nz,nx), arms(nyp), amn(nyp)
  real wfft1(1),wfft2(1),w1(1), w2(1)
  real wfft3(1),wfft4(1)
  real urms(nyp), vrms(nyp), wrms(nyp), umn(nyp)
  real up(nyp,nz,nx), vp(nyp,nz,nx), wp(nyp,nz,nx), scp(nyp,nz,nx) ! Ryan 6/28/22
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
 
  call ccopy (nz*nxh*nyp,scalar,1,wrkc,1)
  call xyzfft(wrkc,w1,w2,wfft1,wfft2,wfft3,wfft4,1)
  call unscram(wrkc,wrk1)
  write(60,*) (((wrk1(i,j,k),i=1,nyp),j=1,nz),k=1,nx)
  close(60)

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
     

subroutine vcw3d(u,v,w,omx,omy,omz,fn,gn,scalar,sclx,scly,sclz,scn,     &
           c11,c12,c13,c21,c22,c23,c31,c32,c33,                   &
           u11,u12,u13,u21,u22,u23,u31,u32,u33,Lu,Lv,Lw,          &
           dc111,dc112,dc113,dc211,dc212,dc213,dc311,dc312,dc313, &
           dc121,dc122,dc123,dc221,dc222,dc223,dc321,dc322,dc323, &
           dc131,dc132,dc133,dc231,dc232,dc233,dc331,dc332,dc333, &
           c11n,c12n,c13n,c22n,c23n,c33n,                         &
           str11n,str12n,str13n,str22n,str23n,str33n,             &
           qp11,qp12,qp13,qp22,qp23,qp33) 
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
    complex fn(nyp,nz,nxh), gn(nyp,nz,nxh), scn(nyp,nz,nxh) ! Ryan 2-24-22
    complex sclx(nyp,nz,nxh), scly(nyp,nz,nxh), sclz(nyp,nz,nxh) ! Ryan 2-24-22
  
  !    local variables:
    real up(mzp,mxp2), vp(mzp,mxp2), wp(mzp,mxp2)
    real wx(mzp,mxp2), wy(mzp,mxp2), wz(mzp,mxp2)
    real vwx(mzp,mxp2), vwy(mzp,mxp2), vwz(mzp,mxp2)
    real wrk(nwrk)
    real cx(mzp,mxp2), cy(mzp,mxp2), cz(mzp,mxp2) ! Ryan 2-24-22
    real vc(mzp,mxp2) ! Ryan 2-24-22
  
! =================================================================================== !
!                       Adding vcw3dp stuff - Ryan 6/28/22                            !
! =================================================================================== !
    complex scalar(nyp,nz,nxh) ! scalar field added 1/17/2020
    complex c11n(nyp,nz,nxh),c12n(nyp,nz,nxh),c13n(nyp,nz,nxh) ! nonlinear terms 
    complex c22n(nyp,nz,nxh),c23n(nyp,nz,nxh),c33n(nyp,nz,nxh)
    complex str11n(nyp,nz,nxh),str12n(nyp,nz,nxh)
    complex str13n(nyp,nz,nxh),str22n(nyp,nz,nxh)
    complex str23n(nyp,nz,nxh),str33n(nyp,nz,nxh)

    real c11np(mzp,mxp2),c12np(mzp,mxp2),c13np(mzp,mxp2)  ! nonlinear terms
    real c22np(mzp,mxp2),c23np(mzp,mxp2),c33np(mzp,mxp2)
    real str11np(mzp,mxp2),str12np(mzp,mxp2),str13np(mzp,mxp2)
    real str22np(mzp,mxp2),str23np(mzp,mxp2),str33np(mzp,mxp2)

    real qp11np(mzp,mxp2),qp12np(mzp,mxp2),qp13np(mzp,mxp2) ! new targeting arrays 1/16/2020
    real qp22np(mzp,mxp2),qp23np(mzp,mxp2),qp33np(mzp,mxp2)

    real zbeta1

    real scp(mzp,mxp2) ! added by rah 1/16/2020
    real beta_poly(mzp,mxp2) ! added by Ryan 9/23/22

    real c11p(mzp,mxp2),c12p(mzp,mxp2),c13p(mzp,mxp2)
    real c21p(mzp,mxp2),c22p(mzp,mxp2),c23p(mzp,mxp2)
    real c31p(mzp,mxp2),c32p(mzp,mxp2),c33p(mzp,mxp2)      

    real u11p(mzp,mxp2),u12p(mzp,mxp2),u13p(mzp,mxp2)
    real u21p(mzp,mxp2),u22p(mzp,mxp2),u23p(mzp,mxp2) 
    real u31p(mzp,mxp2),u32p(mzp,mxp2),u33p(mzp,mxp2)
    real,dimension(mzp,mxp2) :: Lup,Lvp,Lwp ! Laplacian terms - Ryan 7/24/23
    
    real dc111p(mzp,mxp2),dc112p(mzp,mxp2),dc113p(mzp,mxp2) 
    real dc211p(mzp,mxp2),dc212p(mzp,mxp2),dc213p(mzp,mxp2)
    real dc311p(mzp,mxp2),dc312p(mzp,mxp2),dc313p(mzp,mxp2)  

    real dc121p(mzp,mxp2),dc122p(mzp,mxp2),dc123p(mzp,mxp2)
    real dc221p(mzp,mxp2),dc222p(mzp,mxp2),dc223p(mzp,mxp2)
    real dc321p(mzp,mxp2),dc322p(mzp,mxp2),dc323p(mzp,mxp2)

    real dc131p(mzp,mxp2),dc132p(mzp,mxp2),dc133p(mzp,mxp2)  
    real dc231p(mzp,mxp2),dc232p(mzp,mxp2),dc233p(mzp,mxp2)  
    real dc331p(mzp,mxp2),dc332p(mzp,mxp2),dc333p(mzp,mxp2) 

    real trp(mzp,mxp2)

!      conformation tensor
    complex c11(nyp,nz,nxh),c12(nyp,nz,nxh),c13(nyp,nz,nxh)
    complex c21(nyp,nz,nxh),c22(nyp,nz,nxh),c23(nyp,nz,nxh)
    complex c31(nyp,nz,nxh),c32(nyp,nz,nxh),c33(nyp,nz,nxh)      
!      velocity gradient tensor 
    complex u11(nyp,nz,nxh),u12(nyp,nz,nxh),u13(nyp,nz,nxh)
    complex u21(nyp,nz,nxh),u22(nyp,nz,nxh),u23(nyp,nz,nxh) 
    complex u31(nyp,nz,nxh),u32(nyp,nz,nxh),u33(nyp,nz,nxh)
    complex,dimension(nyp,nz,nxh) :: Lu,Lv,Lw  ! Laplacian terms - Ryan 7/24/23
!      derivatives of conformation tensor
    complex dc111(nyp,nz,nxh),dc112(nyp,nz,nxh),dc113(nyp,nz,nxh) 
    complex dc211(nyp,nz,nxh),dc212(nyp,nz,nxh),dc213(nyp,nz,nxh)
    complex dc311(nyp,nz,nxh),dc312(nyp,nz,nxh),dc313(nyp,nz,nxh)  

    complex dc121(nyp,nz,nxh),dc122(nyp,nz,nxh),dc123(nyp,nz,nxh)
    complex dc221(nyp,nz,nxh),dc222(nyp,nz,nxh),dc223(nyp,nz,nxh)
    complex dc321(nyp,nz,nxh),dc322(nyp,nz,nxh),dc323(nyp,nz,nxh)

    complex dc131(nyp,nz,nxh),dc132(nyp,nz,nxh),dc133(nyp,nz,nxh)  
    complex dc231(nyp,nz,nxh),dc232(nyp,nz,nxh),dc233(nyp,nz,nxh)  
    complex dc331(nyp,nz,nxh),dc332(nyp,nz,nxh),dc333(nyp,nz,nxh)
!      polymer stresses used for targeting
    complex qp11(nyp,nz,nxh),qp12(nyp,nz,nxh),qp13(nyp,nz,nxh)
    complex qp22(nyp,nz,nxh),qp23(nyp,nz,nxh)
    complex qp33(nyp,nz,nxh)

!   xyzfft
    real wfft1(nmax), wfft2(nmax)
    real wfft3(nmax), wfft4(nmax)
    real w1(mx,mz), w2(mx,mz)
    real tmp2(nyp,nz,nx)

!   GMU printing
    real argy, sumens, enstrophy, zmax, locx, locz, sumcirc, dtheta, KE
    real Lx, Ly, Lz, volume, mtotal

    real zlmax,tpoly,qbeta,diff,diffpoly,deltat,c11z,c22z,c33z,alpha_poly,mpoly
    real ampbdy1,ampbdy2,c11amp,tfrac,polyrate
    integer ipolyflag,itarget,ipeter,scl_flag,src_start,src_stop

    common/polymer1/ zlmax,tpoly,qbeta,diff,diffpoly,deltat,c11z,c22z,c33z,alpha_poly,mpoly
    common/polymer2/ ampbdy1,ampbdy2,c11amp,tfrac,polyrate,src_start,src_stop
    common/polymer3/ ipolyflag,itarget,ipeter,scl_flag

! =================================================================================== !

    real bdyfx, L, rad, xcenter, ycenter, zcenter 
    common/vortRing/ bdyfx, L, rad, xcenter, ycenter, zcenter  

    real dragx(mz,mx), dragy(mz,mx), dragz(mz,mx)
    real cfl(nyp)
  
    real cflcheck, cflmax, x, y, z

    ! Extra variables for vortex ring
    real xi, zj, argx, argrad, fx, fr, rsq
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
 
    real :: dudx, dvdx, dwdx, dudy, dvdy, dwdy, dudz, dvdz, dwdz
    real dely
    real swirl
    real, dimension(nyp,mz,mx) :: lamb2_3d, swirl_3d
    integer num_threads 
    integer maxrec

    real :: dPdx, R_tau
    common /pressure/ dPdx, R_tau
   
    real :: xp, yp, zp, u_interp, ureal, perr
  
    integer geomtype, flow_select
    common /geometry/ geomtype, flow_select

    ! TKE stuff
    real,save :: turbKE2
    real :: turbKE, turbKE1, massFlux, C_f, Ubulk
    real, dimension(mz) :: uxmean, vxmean, wxmean
    real,dimension(nyp) :: uzmean, vzmean, wzmean

    ! friction velocity
    real :: dumdy,u_tau

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
 
  ! Adding scalar/C to omp directive - Ryan 3-17-22 
!-------------------------------------------------------------------------------!
  !omp parallel do default(shared) private(k,j,i,jj,i1,i2,up,vp,wp,     &
  !omp   wx,wy,wz,vwx,vwy,vwz,dragx,dragy,dragz,wrk,inc,isign,          &
  !omp   jump,lot,ii,ipii,idif,jdif,segdrag,jpjj,cx,cy,cz,vc,           &
  !omp   xsegdrag,ysegdrag,zsegdrag,cflcheck,y,x,z,iiter,jiter) schedule(dynamic)

  !$omp parallel do default(shared) private(k,j,i,jj,i1,i2,up,vp,wp,scp,beta_poly, &
  !$omp    wrk,inc,isign,jump,lot,ii,wx,wy,wz,vwx,vwy,vwz,dragx,dragy,dragz,       &
  !$omp    xsegdrag,ysegdrag,zsegdrag,cflcheck,y,x,z,iiter,jiter,                  &
  !$omp    ipii,idif,jdif,segdrag,jpjj,cx,cy,cz,vc,Lup,Lvp,Lwp,                    &
  !$omp    c11p,c12p,c13p,c21p,c22p,c23p,c31p,c32p,c33p,u11p,u12p,u13p,u21p,       &
  !$omp    u22p,u23p,u31p,u32p,u33p,dc111p,dc112p,dc113p,dc211p,dc212p,dc213p,     &
  !$omp    dc311p,dc312p,dc313p,dc121p,dc122p,dc123p,dc221p,dc222p,dc223p,         &
  !$omp    dc321p,dc322p,dc323p,dc131p,dc132p,dc133p,dc231p,dc232p,dc233p,         &
  !$omp    dc331p,dc332p,dc333p,trp,c11np,c12np,c13np,c22np,c23np,c33np,           &
  !$omp    str11np,str12np,str13np,str22np,str23np,str33np,                        &
  !$omp    qp11np,qp12np,qp13np,qp22np,qp23np,qp33np,zbeta1) schedule(dynamic)
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
            
! Terms from polymer code - Ryan 2-24-22
             if (scl_flag .ge. 1) then
             scp(j,i) = 0.0
             cx(j,i) = 0.0
             cy(j,i) = 0.0
             cz(j,i) = 0.0
             vc(j,i) = 0.0
             end if
             
             if (ipolyflag .eq. 1 .and. it .ge. src_start-1) then
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
! --------------------------- !
  
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
            
! Terms from polymer code - Ryan 2-24-22 
             if (scl_flag .ge. 1) then
             scp(jj,i1) = real(scalar(k,j,i))
             cx(jj,i1) = real(sclx(k,j,i))
             cy(jj,i1) = real(scly(k,j,i))
             cz(jj,i1) = real(sclz(k,j,i))
             end if
             if (ipolyflag .eq. 1 .and. it .ge. src_start-1) then
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
             end if
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
           
! Terms from polymer code - Ryan 2-24-22 
             if (scl_flag .ge. 1) then
             scp(jj,i2) = aimag(scalar(k,j,i))
             cx(jj,i2) = aimag(sclx(k,j,i))
             cy(jj,i2) = aimag(scly(k,j,i))
             cz(jj,i2) = aimag(sclz(k,j,i))
             end if 
             if (ipolyflag .eq. 1 .and. it .ge. src_start-1) then
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
     
       if (scl_flag .ge. 1) then
       call cfftmlt(scp(1,1),scp(1,2),wrk,trigz32,izfax32,inc,jump,mz,lot,isign)
       call cfftmlt(cx(1,1),cx(1,2),wrk,trigz32,izfax32,inc,jump,mz,lot,isign)
       call cfftmlt(cy(1,1),cy(1,2),wrk,trigz32,izfax32,inc,jump,mz,lot,isign)
       call cfftmlt(cz(1,1),cz(1,2),wrk,trigz32,izfax32,inc,jump,mz,lot,isign)
       end if 
       if (ipolyflag .eq. 1 .and. it .ge. src_start-1) then
       call cfftmlt(c11p(1,1),c11p(1,2),wrk,trigz32,izfax32,inc,jump,mz,lot,isign)
       call cfftmlt(c12p(1,1),c12p(1,2),wrk,trigz32,izfax32,inc,jump,mz,lot,isign)
       call cfftmlt(c13p(1,1),c13p(1,2),wrk,trigz32,izfax32,inc,jump,mz,lot,isign)
       call cfftmlt(c21p(1,1),c21p(1,2),wrk,trigz32,izfax32,inc,jump,mz,lot,isign)
       call cfftmlt(c22p(1,1),c22p(1,2),wrk,trigz32,izfax32,inc,jump,mz,lot,isign)
       call cfftmlt(c23p(1,1),c23p(1,2),wrk,trigz32,izfax32,inc,jump,mz,lot,isign)
       call cfftmlt(c31p(1,1),c31p(1,2),wrk,trigz32,izfax32,inc,jump,mz,lot,isign)
       call cfftmlt(c32p(1,1),c32p(1,2),wrk,trigz32,izfax32,inc,jump,mz,lot,isign)
       call cfftmlt(c33p(1,1),c33p(1,2),wrk,trigz32,izfax32,inc,jump,mz,lot,isign)
     
       call cfftmlt(dc111p(1,1),dc111p(1,2),wrk,trigz32,izfax32,inc,jump,mz,lot,isign)
       call cfftmlt(dc112p(1,1),dc112p(1,2),wrk,trigz32,izfax32,inc,jump,mz,lot,isign)
       call cfftmlt(dc113p(1,1),dc113p(1,2),wrk,trigz32,izfax32,inc,jump,mz,lot,isign)
       call cfftmlt(dc211p(1,1),dc211p(1,2),wrk,trigz32,izfax32,inc,jump,mz,lot,isign)
       call cfftmlt(dc212p(1,1),dc212p(1,2),wrk,trigz32,izfax32,inc,jump,mz,lot,isign)
       call cfftmlt(dc213p(1,1),dc213p(1,2),wrk,trigz32,izfax32,inc,jump,mz,lot,isign)
       call cfftmlt(dc311p(1,1),dc311p(1,2),wrk,trigz32,izfax32,inc,jump,mz,lot,isign)
       call cfftmlt(dc312p(1,1),dc312p(1,2),wrk,trigz32,izfax32,inc,jump,mz,lot,isign)
       call cfftmlt(dc313p(1,1),dc313p(1,2),wrk,trigz32,izfax32,inc,jump,mz,lot,isign)
     
       call cfftmlt(dc121p(1,1),dc121p(1,2),wrk,trigz32,izfax32,inc,jump,mz,lot,isign)
       call cfftmlt(dc122p(1,1),dc122p(1,2),wrk,trigz32,izfax32,inc,jump,mz,lot,isign)
       call cfftmlt(dc123p(1,1),dc123p(1,2),wrk,trigz32,izfax32,inc,jump,mz,lot,isign)
       call cfftmlt(dc221p(1,1),dc221p(1,2),wrk,trigz32,izfax32,inc,jump,mz,lot,isign)
       call cfftmlt(dc222p(1,1),dc222p(1,2),wrk,trigz32,izfax32,inc,jump,mz,lot,isign)
       call cfftmlt(dc223p(1,1),dc223p(1,2),wrk,trigz32,izfax32,inc,jump,mz,lot,isign)
       call cfftmlt(dc321p(1,1),dc321p(1,2),wrk,trigz32,izfax32,inc,jump,mz,lot,isign)
       call cfftmlt(dc322p(1,1),dc322p(1,2),wrk,trigz32,izfax32,inc,jump,mz,lot,isign)
       call cfftmlt(dc323p(1,1),dc323p(1,2),wrk,trigz32,izfax32,inc,jump,mz,lot,isign)
     
       call cfftmlt(dc131p(1,1),dc131p(1,2),wrk,trigz32,izfax32,inc,jump,mz,lot,isign)
       call cfftmlt(dc132p(1,1),dc132p(1,2),wrk,trigz32,izfax32,inc,jump,mz,lot,isign)
       call cfftmlt(dc133p(1,1),dc133p(1,2),wrk,trigz32,izfax32,inc,jump,mz,lot,isign)
       call cfftmlt(dc231p(1,1),dc231p(1,2),wrk,trigz32,izfax32,inc,jump,mz,lot,isign)
       call cfftmlt(dc232p(1,1),dc232p(1,2),wrk,trigz32,izfax32,inc,jump,mz,lot,isign)
       call cfftmlt(dc233p(1,1),dc233p(1,2),wrk,trigz32,izfax32,inc,jump,mz,lot,isign)
       call cfftmlt(dc331p(1,1),dc331p(1,2),wrk,trigz32,izfax32,inc,jump,mz,lot,isign)
       call cfftmlt(dc332p(1,1),dc332p(1,2),wrk,trigz32,izfax32,inc,jump,mz,lot,isign)
       call cfftmlt(dc333p(1,1),dc333p(1,2),wrk,trigz32,izfax32,inc,jump,mz,lot,isign)
       end if
 
       do j = 1, mz
          up(j,nxp) = up(j,2)
          vp(j,nxp) = vp(j,2)
          wp(j,nxp) = wp(j,2)
          wx(j,nxp) = wx(j,2)
          wy(j,nxp) = wy(j,2)
          wz(j,nxp) = wz(j,2)
! New variables - Ryan 2-24-22
          scp(j,nxp)=scp(j,2) ! test
          cx(j,nxp) = cx(j,2)
          cy(j,nxp) = cy(j,2)
          cz(j,nxp) = cz(j,2)
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
! New variables - Ryan 2-24-22
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
     
       if (scl_flag .ge. 1) then
          scp(j,2) = 0.0
          scp(j,nxp2) = 0.0 ! test
          cx(j,2) = 0.0
          cy(j,2) = 0.0
          cz(j,2) = 0.0
       end if 
       if (ipolyflag .eq. 1 .and. it .ge. src_start-1) then
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
       end if
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
! New variables - Ryan 2-24-22
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
     
       if (scl_flag .ge. 1) then
       call rfftmlt(scp,wrk,trigx32,ixfax32,inc,jump,mx,lot,isign)
       call rfftmlt(cx,wrk,trigx32,ixfax32,inc,jump,mx,lot,isign)
       call rfftmlt(cy,wrk,trigx32,ixfax32,inc,jump,mx,lot,isign)
       call rfftmlt(cz,wrk,trigx32,ixfax32,inc,jump,mx,lot,isign)
       end if 
       if (ipolyflag .eq. 1 .and. it .ge. src_start-1) then
       call rfftmlt(c11p,wrk,trigx32,ixfax32,inc,jump,mx,lot,isign)
       call rfftmlt(c12p,wrk,trigx32,ixfax32,inc,jump,mx,lot,isign)
       call rfftmlt(c13p,wrk,trigx32,ixfax32,inc,jump,mx,lot,isign)
       call rfftmlt(c21p,wrk,trigx32,ixfax32,inc,jump,mx,lot,isign)
       call rfftmlt(c22p,wrk,trigx32,ixfax32,inc,jump,mx,lot,isign)
       call rfftmlt(c23p,wrk,trigx32,ixfax32,inc,jump,mx,lot,isign)
       call rfftmlt(c31p,wrk,trigx32,ixfax32,inc,jump,mx,lot,isign)
       call rfftmlt(c32p,wrk,trigx32,ixfax32,inc,jump,mx,lot,isign)
       call rfftmlt(c33p,wrk,trigx32,ixfax32,inc,jump,mx,lot,isign)
     
       call rfftmlt(dc111p,wrk,trigx32,ixfax32,inc,jump,mx,lot,isign)
       call rfftmlt(dc112p,wrk,trigx32,ixfax32,inc,jump,mx,lot,isign)
       call rfftmlt(dc113p,wrk,trigx32,ixfax32,inc,jump,mx,lot,isign)
       call rfftmlt(dc211p,wrk,trigx32,ixfax32,inc,jump,mx,lot,isign)
       call rfftmlt(dc212p,wrk,trigx32,ixfax32,inc,jump,mx,lot,isign)
       call rfftmlt(dc213p,wrk,trigx32,ixfax32,inc,jump,mx,lot,isign)
       call rfftmlt(dc311p,wrk,trigx32,ixfax32,inc,jump,mx,lot,isign)
       call rfftmlt(dc312p,wrk,trigx32,ixfax32,inc,jump,mx,lot,isign)
       call rfftmlt(dc313p,wrk,trigx32,ixfax32,inc,jump,mx,lot,isign)
     
       call rfftmlt(dc121p,wrk,trigx32,ixfax32,inc,jump,mx,lot,isign)
       call rfftmlt(dc122p,wrk,trigx32,ixfax32,inc,jump,mx,lot,isign)
       call rfftmlt(dc123p,wrk,trigx32,ixfax32,inc,jump,mx,lot,isign)
       call rfftmlt(dc221p,wrk,trigx32,ixfax32,inc,jump,mx,lot,isign)
       call rfftmlt(dc222p,wrk,trigx32,ixfax32,inc,jump,mx,lot,isign)
       call rfftmlt(dc223p,wrk,trigx32,ixfax32,inc,jump,mx,lot,isign)
       call rfftmlt(dc321p,wrk,trigx32,ixfax32,inc,jump,mx,lot,isign)
       call rfftmlt(dc322p,wrk,trigx32,ixfax32,inc,jump,mx,lot,isign)
       call rfftmlt(dc323p,wrk,trigx32,ixfax32,inc,jump,mx,lot,isign)
     
       call rfftmlt(dc131p,wrk,trigx32,ixfax32,inc,jump,mx,lot,isign)
       call rfftmlt(dc132p,wrk,trigx32,ixfax32,inc,jump,mx,lot,isign)
       call rfftmlt(dc133p,wrk,trigx32,ixfax32,inc,jump,mx,lot,isign)
       call rfftmlt(dc231p,wrk,trigx32,ixfax32,inc,jump,mx,lot,isign)
       call rfftmlt(dc232p,wrk,trigx32,ixfax32,inc,jump,mx,lot,isign)
       call rfftmlt(dc233p,wrk,trigx32,ixfax32,inc,jump,mx,lot,isign)
       call rfftmlt(dc331p,wrk,trigx32,ixfax32,inc,jump,mx,lot,isign)
       call rfftmlt(dc332p,wrk,trigx32,ixfax32,inc,jump,mx,lot,isign)
       call rfftmlt(dc333p,wrk,trigx32,ixfax32,inc,jump,mx,lot,isign)
       end if  
  !......................................................................
  
  !    compute the cross product of velocity and vorticity.
  
       do j = 1, mz
          do i = 1, mx
             vwx(j,i) =  vp(j,i) * wz(j,i) - wp(j,i) * wy(j,i)
             vwy(j,i) = -up(j,i) * wz(j,i) + wp(j,i) * wx(j,i)
             vwz(j,i) =  up(j,i) * wy(j,i) - vp(j,i) * wx(j,i)
! New variable - Ryan 2-24-22
       if (scl_flag .ge. 1) then
             vc(j,i) = -(up(j,i) * cx(j,i) + vp(j,i) * cy(j,i) + wp(j,i) * cz(j,i))
       end if
       if (ipolyflag .eq. 1 .and. it .ge. src_start-1) then
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
     
             trp(j,i) = c11p(j,i)+c22p(j,i)+c33p(j,i)
     
             if (ipeter .eq. 0) then  ! peterlin function is 1.0
    
             str11np(j,i)= c11p(j,i)
             str12np(j,i)= c12p(j,i)
             str13np(j,i)= c13p(j,i)
             str22np(j,i)= c22p(j,i)
             str23np(j,i)= c23p(j,i)
             str33np(j,i)= c33p(j,i)
         
             else
          
             str11np(j,i)=( (zlmax**2-3.0)/(zlmax**2-trp(j,i)) )*c11p(j,i)
             str12np(j,i)=( (zlmax**2-3.0)/(zlmax**2-trp(j,i)) )*c12p(j,i)
             str13np(j,i)=( (zlmax**2-3.0)/(zlmax**2-trp(j,i)) )*c13p(j,i)
             str22np(j,i)=( (zlmax**2-3.0)/(zlmax**2-trp(j,i)) )*c22p(j,i)
             str23np(j,i)=( (zlmax**2-3.0)/(zlmax**2-trp(j,i)) )*c23p(j,i)
             str33np(j,i)=( (zlmax**2-3.0)/(zlmax**2-trp(j,i)) )*c33p(j,i)
          
             end if
          
             ! Add brownian motion terms
             str11np(j,i)=str11np(j,i)-1.0
             str22np(j,i)=str22np(j,i)-1.0
             str33np(j,i)=str33np(j,i)-1.0

!             beta_poly(j,i) = 0.9

             if (it .lt. src_start) then
                 beta_poly(j,i) = 1.0
             else if (itarget .eq. 0) then ! polymer ocean
                 beta_poly(j,i) = qbeta
             else if (itarget .eq. 1) then
                 ! New nonlinear polymer model 10/3/22 - Ryan
                 ! nu_0 - nu_s = nu_s (exp(alpha*gamma) - 1)
                 ! alpha --> alpha_poly = 2.6E-03 PPM (empirical constant)
                 ! gamma --> scp = scalar concentration (in PPM)

                 beta_poly(j,i) = exp(-alpha_poly*abs(scp(j,i))) ! Define the viscosity ratio based on model
!                beta_poly(j,i) = qbeta ! test case only

             end if

             zbeta1 = (1.0-beta_poly(j,i))/(re*beta_poly(j,i)*tpoly) ! = (nu_0 - nu_s)

             qp11np(j,i)= zbeta1*str11np(j,i)
             qp12np(j,i)= zbeta1*str12np(j,i)
             qp13np(j,i)= zbeta1*str13np(j,i)
             qp22np(j,i)= zbeta1*str22np(j,i)
             qp23np(j,i)= zbeta1*str23np(j,i)
             qp33np(j,i)= zbeta1*str33np(j,i)
 
         end if 

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
       scp3d(k,1:mz,1:mx) = (scp(1:mz,1:mx)) 
       beta3d(k,1:mz,1:mx) = beta_poly(1:mz,1:mx) ! Added by Ryan 9/23/22
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
       if (scl_flag .ge. 1) then
       call rfftmlt(vc ,wrk,trigx32,ixfax32,inc,jump,mx,lot,isign) ! Ryan (sometime in March, idk)
       end if
       if (ipolyflag .eq. 1 .and. it .ge. src_start-1) then
       call rfftmlt(c11np,wrk,trigx32,ixfax32,inc,jump,mx,lot,isign)
       call rfftmlt(c12np,wrk,trigx32,ixfax32,inc,jump,mx,lot,isign)
       call rfftmlt(c13np,wrk,trigx32,ixfax32,inc,jump,mx,lot,isign)
       call rfftmlt(c22np,wrk,trigx32,ixfax32,inc,jump,mx,lot,isign)
       call rfftmlt(c23np,wrk,trigx32,ixfax32,inc,jump,mx,lot,isign)
       call rfftmlt(c33np,wrk,trigx32,ixfax32,inc,jump,mx,lot,isign)
     
       call rfftmlt(str11np,wrk,trigx32,ixfax32,inc,jump,mx,lot,isign)
       call rfftmlt(str12np,wrk,trigx32,ixfax32,inc,jump,mx,lot,isign)
       call rfftmlt(str13np,wrk,trigx32,ixfax32,inc,jump,mx,lot,isign)
       call rfftmlt(str22np,wrk,trigx32,ixfax32,inc,jump,mx,lot,isign)
       call rfftmlt(str23np,wrk,trigx32,ixfax32,inc,jump,mx,lot,isign)
       call rfftmlt(str33np,wrk,trigx32,ixfax32,inc,jump,mx,lot,isign)
     
       call rfftmlt(qp11np,wrk,trigx32,ixfax32,inc,jump,mx,lot,isign)
       call rfftmlt(qp12np,wrk,trigx32,ixfax32,inc,jump,mx,lot,isign)
       call rfftmlt(qp13np,wrk,trigx32,ixfax32,inc,jump,mx,lot,isign)
       call rfftmlt(qp22np,wrk,trigx32,ixfax32,inc,jump,mx,lot,isign)
       call rfftmlt(qp23np,wrk,trigx32,ixfax32,inc,jump,mx,lot,isign)
       call rfftmlt(qp33np,wrk,trigx32,ixfax32,inc,jump,mx,lot,isign)
       end if 
  !now wiping out unwanted last mode.
       do j=1,mz
          vwx(j,2) = 0.0          
          vwy(j,2) = 0.0          
          vwz(j,2) = 0.0
       if (scl_flag .ge. 1) then
          vc(j,2)  = 0.0
       end if
       if (ipolyflag .eq. 1 .and. it .ge. src_start-1) then
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
       end if
       end do
  
       inc = 1
       isign = -1
       jump = 2*mzp
       lot = nx/2
  
       call cfftmlt(vwx(1,1),vwx(1,2),wrk,trigz32,izfax32,inc,jump,mz,lot,isign)
       call cfftmlt(vwy(1,1),vwy(1,2),wrk,trigz32,izfax32,inc,jump,mz,lot,isign)
       call cfftmlt(vwz(1,1),vwz(1,2),wrk,trigz32,izfax32,inc,jump,mz,lot,isign)
       if (scl_flag .ge. 1) then
       call cfftmlt(vc(1,1) ,vc(1,2) ,wrk,trigz32,izfax32,inc,jump,mz,lot,isign)
       end if
       if (ipolyflag .eq. 1 .and. it .ge. src_start-1) then
       call cfftmlt(c11np(1,1),c11np(1,2),wrk,trigz32,izfax32,inc,jump,mz,lot,isign)
       call cfftmlt(c12np(1,1),c12np(1,2),wrk,trigz32,izfax32,inc,jump,mz,lot,isign)
       call cfftmlt(c13np(1,1),c13np(1,2),wrk,trigz32,izfax32,inc,jump,mz,lot,isign)
       call cfftmlt(c22np(1,1),c22np(1,2),wrk,trigz32,izfax32,inc,jump,mz,lot,isign)
       call cfftmlt(c23np(1,1),c23np(1,2),wrk,trigz32,izfax32,inc,jump,mz,lot,isign)
       call cfftmlt(c33np(1,1),c33np(1,2),wrk,trigz32,izfax32,inc,jump,mz,lot,isign)
     
       call cfftmlt(str11np(1,1),str11np(1,2),wrk,trigz32,izfax32,inc,jump,mz,lot,isign)
       call cfftmlt(str12np(1,1),str12np(1,2),wrk,trigz32,izfax32,inc,jump,mz,lot,isign)
       call cfftmlt(str13np(1,1),str13np(1,2),wrk,trigz32,izfax32,inc,jump,mz,lot,isign)
       call cfftmlt(str22np(1,1),str22np(1,2),wrk,trigz32,izfax32,inc,jump,mz,lot,isign)
       call cfftmlt(str23np(1,1),str23np(1,2),wrk,trigz32,izfax32,inc,jump,mz,lot,isign)
       call cfftmlt(str33np(1,1),str33np(1,2),wrk,trigz32,izfax32,inc,jump,mz,lot,isign)
     
       call cfftmlt(qp11np(1,1),qp11np(1,2),wrk,trigz32,izfax32,inc,jump,mz,lot,isign)
       call cfftmlt(qp12np(1,1),qp12np(1,2),wrk,trigz32,izfax32,inc,jump,mz,lot,isign)
       call cfftmlt(qp13np(1,1),qp13np(1,2),wrk,trigz32,izfax32,inc,jump,mz,lot,isign)
       call cfftmlt(qp22np(1,1),qp22np(1,2),wrk,trigz32,izfax32,inc,jump,mz,lot,isign)
       call cfftmlt(qp23np(1,1),qp23np(1,2),wrk,trigz32,izfax32,inc,jump,mz,lot,isign)
       call cfftmlt(qp33np(1,1),qp33np(1,2),wrk,trigz32,izfax32,inc,jump,mz,lot,isign)
       end if 
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
          if (scl_flag .ge. 1) then
             scn(k,j,i) = cmplx(vc(jj,i1) *rmz, vc(jj,i2) *rmz)
          end if
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
                call write_flowfield_ascii(it, up3d, vp3d, wp3d, wx3d, wy3d, wz3d, swirl_3d, real(imatrix), beta3d, xl, zl, print3d) ! Changed scp3d to beta3d - Ryan 9/23/22
            else 
                call write_flowfield_plt(it, up3d, vp3d, wp3d, wx3d, wy3d, wz3d, swirl_3d, real(imatrix), beta3d, xl, zl)
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
          call part_track(up3d,vp3d,wp3d,wx3d,wy3d,wz3d,u_old,v_old,w_old, &
                          u11p3d,u12p3d,u13p3d,u21p3d,u22p3d,u23p3d,u31p3d,u32p3d,u33p3d,Lup3d,Lvp3d,Lwp3d,Lup1,Lvp1,Lwp1)
       else
           call part_track(up1,vp1,wp1,wx1,wy1,wz1,u_old,v_old,w_old, &
                          u11p1,u12p1,u13p1,u21p1,u22p1,u23p1,u31p1,u32p1,u33p1)
       end if
    end if

    if (it .eq. 1 .or. frozen_flag .eq. .false.) then
       u_old = up3d
       v_old = vp3d
       w_old = wp3d
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
                vxmean(j) = sum(vp3d(k,j,:))/mx
                wxmean(j) = sum(wp3d(k,j,:))/mx
            end do
            uzmean(k) = sum(uxmean)/mz
            vzmean(k) = sum(vxmean)/mz
            wzmean(k) = sum(wxmean)/mz

            write(71,"(*(e14.6,1x))") uzmean(k) ! up3d(k,floor(mz/2.0),floor(mx/2.0))
        end do
        close(71)

        ! Computing friction velocity, u_tau = sqrt(nu*d<U>/dy|(y=0))
        ! d<U>/dy
        dumdy = (uzmean(2) - uzmean(1))/(ycoord(2) - ycoord(1)) ! forward diff. derivative
        u_tau = sqrt(dumdy/re)
        print *,'u_tau = ',u_tau

        ! compute global (u(t) - <u>)^2 * dt
        turbKE1 = 0.0
        do i = 1,nyp
          do j = 1,mz
            do k = 1,mx
              turbKE1 = turbKE1 + 0.5*((up3d(i,j,k) - uzmean(i))**2 + (vp3d(i,j,k) - vzmean(i))**2 + (wp3d(i,j,k) - wzmean(i))**2)
            end do
          end do
        end do

        ! Write TKE to file
        if (irstrt .eq. it) then
            open(72,file="outputs/TKE")
            turbKE = 0.0
            turbKE2 = 0.0
        else
            open(72,file="outputs/TKE",position="append")
            ! Add current TKE info to old to integrate
            turbKE2 = turbKE2 + turbKE1*dt

            ! Compute total TKE at this time step
            turbKE = 1.0/(float(it)*dt)*(turbKE2)
        end if

        write(72,*) turbKE
        close(72)

        ! Computing Mass Flux: integrate <U> over y - Ryan 8/2/23
        if (irstrt .eq. it) then
            open(73,file="outputs/mass_flux")
            open(74,file="outputs/C_f")
        else
            open(73,file="outputs/mass_flux",position="append")
            open(74,file="outputs/C_f",position="append")
        end if

        massFlux = 0.0
        Ubulk = 0.0
        do i = 1,ny 
            Ubulk = Ubulk + 1.0/yl*0.5*(uzmean(i+1) + uzmean(i))*(ycoord(i+1) - ycoord(i)) ! bulk velocity
        end do
        C_f = dPdx/(0.5*Ubulk**2)
        massFlux = Ubulk

        write(73,*) massFlux
        write(74,*) C_f
        close(73)
        close(74)

    print *,'Finished writing mean U data'
    end if

    !-----------------------------------------------------------------------------------!

    ! Computing Global KE and Enstrophy for all flow cases
!!!! If I ever want to use this again, I need to fix the Ly calculation (see integration below) !!!!
!    if (print3d .ge. 0) then
!        sumens = 0.0
!        KE = 0.0
!
!        do i = 1,nyp
!          do j = 1,mz
!            do k = 1,mx
!
!              ! Calculate cell volume
!              ! Calc x length
!              if (k .eq. 1 .or. k .eq. mx) then
!                  Lx = delxm/2.0
!              else
!                  Lx = delxm
!              end if
!
!              ! Calc z length
!              if (j .eq. 1 .or. j .eq. mz) then
!                  Lz = delzm/2.0
!              else
!                  Lz = delzm
!              end if
!
!              ! Calc y length
!              if (i .eq. 1) then
!                  Ly = ycoord(2) - ycoord(1)
!              else if (i .eq. nyp) then
!                  Ly = ycoord(nyp) - ycoord(ny)
!              else
!                  Ly = ycoord(i+1) - ycoord(i-1)
!              end if
!
!              ! Calculate cell volume
!              volume = Lx*Ly*Lz
!           
!              ! Calculate global enstrophy (=volume integral of |omega|^2) 
!              sumens = sumens + (wx3d(i,j,k)**2 + wy3d(i,j,k)**2 + wz3d(i,j,k)**2)*volume
!
!              ! Calculate global kinetic energy (=1/2 volume integral of |v|^2)
!              KE = KE + 0.5*(up3d(i,j,k)**2 + vp3d(i,j,k)**2 + wp3d(i,j,i)**2)*volume
!            end do
!          end do
!        end do
!
!        open(45,file='outputs/enstrophy',position="append")
!        write(45,*) sumens
!        close(45)
!
!        open(46,file='outputs/KE',position="append")
!        write(46,*) KE
!        close(46)
!    end if
                     
 
    !-----------------------------------------------------------------------------------!


    ! Integrating scalar to get total mass of polymer added -Ryan 8/23/23
    ! Updating to stop the source once it reaches a given mass - Ryan 8/31/23

!    if (it .lt. src_stop .and. it .ge. src_start + 100) then
    if (it .lt. src_stop .and. it .ge. src_start) then
        ! check
        print *,'integrating mass... '
        mtotal = 0.0
		KE = 0.0
		volume = xl*yl*zl
        !$omp parallel do reduction(+:mtotal,KE) default(shared) private(i,j,k,Lx,Ly,Lz)
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
                    ! Calculate cell volume
                    ! Calc x length
                    if (k .eq. 1 .or. k .eq. mx) then
                        Lx = delxm/2.0
                    else
                        Lx = delxm
                    end if

                    mtotal = mtotal + scp3d(i,j,k)/1000.0*Lx*Ly*Lz
					KE = KE + beta3d(i,j,k)*Lx*Ly*Lz ! Actually avg beta, but using KE since it's already a variable
                end do
            end do
        end do
        !$omp end parallel do

		print *,'avg beta = ',KE/volume
        if (mtotal .ge. mpoly) then
            open(100,file="outputs/total_mass")
            write(100,*) mtotal
            close(100)
            src_stop = it
        end if
    end if
    !-----------------------------------------------------------------------------------!
    !                       Adding GMU Printing Stuff (Ryan 7/15/22)                    !
    !-----------------------------------------------------------------------------------!
    if (print3d .eq. -1) then
        open(unit=34, file="GMU_out/fort.34",form='formatted',status='unknown')
        write(34,*) (((wx3d(i,j,kk),i=1,nyp),j=1,mz),kk=1,mx)
        write(34,*) (((wy3d(i,j,kk),i=1,nyp),j=1,mz),kk=1,mx)
        write(34,*) (((wz3d(i,j,kk),i=1,nyp),j=1,mz),kk=1,mx)
        close(34)
      
10      format(4e12.5)

!       calculate enstrophy

        dtheta = pi/float(ny)
        open(unit=45,file='GMU_out/enstrophy',form='formatted',status='unknown')
        open(unit=50,file='GMU_out/omymax',form='formatted',status='unknown')
        open(unit=51,file='GMU_out/xloc',form='formatted',status='unknown')
        open(unit=52,file='GMU_out/zloc',form='formatted',status='unknown')
        open(unit=53,file='GMU_out/circulation',form='formatted',status='unknown')

        if(mod((it-1),iprnfrq).eq.0) then

         sumens = 0.0
         sumcirc = 0.0

!        compute global average enstrophy
         do i = 1,nyp
         argy = pi*float(i-1)/float(ny)
         do j = 1,mz
         do kk = 1,mx
         sumens= sumens+ (wx3d(i,j,kk)**2 + wy3d(i,j,kk)**2 + wz3d(i,j,kk)**2)*sin(argy)
         enddo
         enddo
         enddo
         enstrophy = ((sumens*dtheta)/2.0)/float(mx*mz)

!        compute max vorty
         zmax = 0.0
         i = ny/2 +1 ! center plane

         do j = 1,mz/2+1
         do kk = 1,mx  ! changed by rah 10/5/21

          if (abs(wy3d(i,j,kk)) .gt. zmax) then
              zmax = abs(wy3d(i,j,kk))
              locx = kk
              locz = j
           endif

         enddo
         enddo

!        compute circulation using center plane i=ny/2+1 which is defined above
         do j = 1,mz
         do kk = 1,mx
         sumcirc = sumcirc + abs(wy3d(i,j,kk))
         enddo
         enddo
         sumcirc = sumcirc/2.0  !divide by 2 because of symmetry of vring

         write(45,1234) time,enstrophy   !enstrophy
         write(50,1234) time,zmax        !omy
         write(51,1235) time,locx        !xloc
         write(52,1235) time,locz        !zloc
         write(53,1234) time,sumcirc
         close(45)
         close(50)
         close(51)
         close(52)
         close(53)
1234        format(e20.13,1x,e20.13)
1235        format(e20.13,1x,i5)
        endif  !outer loop

    end if

                 
    !-----------------------------------------------------------------------------------!


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
subroutine write_flowfield_ascii(it, u, v, w, wx, wy, wz, swirl, ctp, &
                                 scl, xl, zl, print3d)
    use grid_size
    integer :: it, print3d, mz_copy
    integer :: irstrt, nsteps,iprnfrq
    common/iocntrl/ irstrt,nsteps,iprnfrq
    real, dimension(nyp,mz,mx) :: u, v, w, wx, wy, wz, ctp, scl, swirl
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

    delxm = xl /float(mx-1)
    delzm = zl /float(mz-1)

    print *, 'Writing data at time step ', it
    
    write(filename,'("outputs/flowfield/time-",i6.6,".dat")')it
    
    ! Write 3d-data 
    open(6, file = filename)
    ! added scalar - Ryan 6/28/22
    write(6,9711) 'filetype = solution, variables = "swirl", "u", "v", "w", "wx", "wy", "wz", "ctp"'!, "scalar"'
    write(6,9712) 'zone f=point t="Field", solutiontime=', it/iprnfrq,',i=',mx, 'j=',1, 'k=', nyp, new_line('a')
    
    do k = 1,nyp
        y = ycoord(k)
        do j = 1,mz_copy
            z = float(j-1)*delzm 
            do i = 1,mx
                x = float(i-1)*delxm
                ! added scl - Ryan 6/28/22
                write(6,8611) swirl(k,j,i),u(k,j,i),v(k,j,i),w(k,j,i),wx(k,j,i),wy(k,j,i),wz(k,j,i),ctp(k,j,i)!,scl(k,j,i)
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

subroutine write_flowfield_plt(it, u, v, w, wx, wy, wz, swirl, ctp, scl, xl, zl)
! Note that the scl variable is actually beta, but I was too lazy to change the
! name - Ryan 9/23/22
    use iso_c_binding 
    use grid_size
    implicit none

    include 'tecio.f90'

    integer :: it, mz_copy
    integer :: irstrt, nsteps,iprnfrq
    common/iocntrl/ irstrt,nsteps,iprnfrq
    real, dimension(nyp,mz,mx) :: u, v, w, wx, wy, wz, ctp, scl
    real, dimension(nyp,mz,mx) :: swirl
    character*1 :: NULLCHR
    character (35) :: filename
    
    integer ipolyflag,itarget,ipeter,scl_flag
    common/polymer3/ ipolyflag,itarget,ipeter,scl_flag


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
    ! Initialize helper arrays
    numValues = (nyp) * (mz_copy) * (mx)
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
            (mx), &
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
        allocate(coordinates(nyp,mz_copy,mx, 3))
        delxm = xl /float(mx-1)
        delzm = zl /float(mz-1)
        do k = 1,nyp
            y = ycoord(k)
            do j = 1,mz_copy
                z = float(j-1) * delzm 
                do i = 1,mx
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
    
    ! Open files
    if (scl_flag .ge. 1) then
    i = tecini142("Solution", "swirl u v w wx wy wz ctp beta", filename//NULLCHR, ".", &
        fileFormat, outputFileType, 1, 0)
    else
    i = tecini142("Solution", "swirl u v w wx wy wz ctp", filename//NULLCHR, ".", &
        fileFormat, outputFileType, 1, 0)
    end if
    
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
        mx, &
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
    
    ! Write scalar
    if (scl_flag .ge. 1) then
    floatValues = pack(scl(:,1:mz_copy,:), .true.)
    i = tecdat142(numValues, floatValues, isDouble)
    end if

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
