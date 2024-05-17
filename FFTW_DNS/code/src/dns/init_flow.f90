subroutine init_flow(initu, initv, initw)

! This set of subroutines initializes the flow field outside the main dns code

! This is the main subroutine which selects the type of flow. If the flow case
! you want to run is already coded, simply change the flow_type selector and
! adjust any appropriate variables defined here

! NOTE: Be sure the geometry and setup files align with the flow type selected

use grid_size

implicit none

real,dimension(nyp,nz,nx) :: initu, initv, initw

! Flow setup variables
integer :: geomtype,flow_select,perturbtime

common/setup/ geomtype,flow_select,perturbtime

! Flow selection decision
if (flow_select .le. 0 .or. flow_select .eq. 5) then
   call still(initu, initv, initw)
else if (flow_select .eq. 1) then
   call channel(initu, initv, initw)
else if (flow_select .eq. 2) then
   call Couette_2D(initu,initv,initw)
else if (flow_select .eq. 3) then
   call Couette_3D(initu,initv,initw)
else if (flow_select .eq. 4) then
   call BlasiusBL(initu,initv,initw)
else if (flow_select .eq. 10) then
   call vortex_only(initu,initv,initw)
else if (flow_select .eq. 11) then
   call vortex_channel(initu,initv,initw)
else if (flow_select .eq. 12) then
   call vortex_couette(initu,initv,initw)
end if

end subroutine

! ------------------------------------------------------------------------- !
!                                                                           !
!                              Still Fluid                                  !
!                                                                           !
! ------------------------------------------------------------------------- !

subroutine still(initu, initv, initw)

! Still fluid - U = 0 everywhere

! Used for vortex ring generation in addition to separate body forces inside
! dns.f90

  use grid_size

  implicit none

  real,dimension(nyp,nz,nx) :: initu, initv, initw
  integer :: i,j,k

  do k = 1,nyp
    do j = 1,nz
      do i = 1,nx
        initu(k,j,i) = 0.0
        initv(k,j,i) = 0.0
        initw(k,j,i) = 0.0
      end do
    end do
  end do


end subroutine

! ------------------------------------------------------------------------- !
!                                                                           !
!                              Channel Flow                                 !
!                                                                           !
! ------------------------------------------------------------------------- !

subroutine channel(initu, initv, initw)

! Pressure-driven (Poiseuille) flow - set R_tau here and geometry in setup files

    use grid_size
    
    implicit none
    
    ! Passed variables 
    real,dimension(nyp,nz,nx) :: initu, initv, initw
    
    ! Calculation variables
    real :: ReyNo, Uc
    integer :: i, j, k, kwall
    integer, dimension(nyp,mz,mx) :: imatrix
    
    ! Geometry variables
    real :: xl,yl,zl
    
    ! Flow setup variables
    real :: re,Uinf,R_tau,dPdx
    
    
    ! ==================================================================== !
    !                        Define common blocks                          !
    ! ==================================================================== !
    
    common/flow/   re,Uinf,R_tau,dPdx
    common/imat/   imatrix,kwall
    common/domain/ xl,yl,zl

    ! -------------------------------------------------------------------- !

    ! Begin calculations
    
    dPdx = -2/yl*(2*R_tau/(yl*re))**2
    Uc = -((yl/2)**2)*re/2*dPdx 
    ReyNo = Uc*yl/2*re ! Calculate channel half height Reynolds number
    
    write(*,*) "R_tau equals ",R_tau
    write(*,*) "U_c equals ", Uc
    write(*,*) "Reynolds Number equals ",ReyNo
    write(*,*) "dPdx equals ",dPdx
    
    do k = 1,nx
      do j = 1,nz
        do i = 1,nyp
           if (i .le. kwall .or. (nyp-i) .le. kwall) then
              initu(i,j,k) = 0.0  
           else  
              initu(i,j,k) = re*dPdx*(0.5*abs(ycoord(i)-ycoord(nyp))**2 - yl/2*abs(ycoord(i)-ycoord(nyp)))
              initv(i,j,k) = 0.0
              initw(i,j,k) = 0.0
          end if
        end do
      end do
    end do
end subroutine

! ------------------------------------------------------------------------- !
!                                                                           !
!                            2D Couette Flow                                !
!                                                                           !
! ------------------------------------------------------------------------- !

subroutine Couette_2D(initu, initv, initw)

    ! 2D Couette flow (Linear velocity profile in x)

    use grid_size
    
    implicit none
   
    ! Passed variables 
    real,dimension(nyp,nz,nx) :: initu, initv, initw
   
    ! Calculation variables 
    integer :: i,j,k
  
    ! Flow setup variables
	real    :: re,Uinf,R_tau,dPdx

    ! Geometry variables
    integer :: kwall,kmaxsurf
    integer, dimension(nyp,mz,mx) :: imatrix

    ! --------------------------------------------------------------- !
    !                         Common Blocks                           !
    ! --------------------------------------------------------------- !
    common/flow/       re,Uinf,R_tau,dPdx
    common/imat/       imatrix,kwall,kmaxsurf
    ! --------------------------------------------------------------- !
    ! --------------------------------------------------------------- !
    !                        Set Initial Flow                         !
    ! --------------------------------------------------------------- !
        do k = 1,nx
            do j = 1,nz
                do i = 1,nyp
                initu(i,j,k) = Uinf*(abs(ycoord(i)-ycoord(nyp))/abs(ycoord(1)-ycoord(nyp)))
                initv(i,j,k) = 0.0
                initw(i,j,k) = 0.0
                end do
            end do
        end do
end subroutine

! ------------------------------------------------------------------------- !
!                                                                           !
!                            3D Couette Flow                                !
!                                                                           !
! ------------------------------------------------------------------------- !

subroutine Couette_3D(initu, initv, initw)

    ! 3D Couette flow (Linear velocity profile in x and z)

    use grid_size
    
    implicit none
   
    ! Passed variables 
    real,dimension(nyp,nz,nx) :: initu, initv, initw
   
    ! Calculation variables 
    integer :: i,j,k
  
    ! Flow setup variables
	real    :: re,Uinf,R_tau,dPdx

    ! Geometry variables
    integer :: kwall,kmaxsurf
    integer, dimension(nyp,mz,mx) :: imatrix

    ! --------------------------------------------------------------- !
    !                         Common Blocks                           !
    ! --------------------------------------------------------------- !
    common/flow/       re,Uinf,R_tau,dPdx
    common/imat/       imatrix,kwall,kmaxsurf
    ! --------------------------------------------------------------- !
    ! --------------------------------------------------------------- !
    !                        Set Initial Flow                         !
    ! --------------------------------------------------------------- !
    do k = 1,nyp
        do j = 1,nz
            do i = 1,nx
                initu(k,j,i) = Uinf*((ycoord(k)-ycoord(nyp))/(ycoord(nyp)-ycoord(nyp)))
                initv(k,j,i) = 0.0
                initw(k,j,i) = Uinf*((ycoord(k)-ycoord(nyp))/(ycoord(nyp)-ycoord(nyp)))
            end do
        end do
    end do
end subroutine

! ------------------------------------------------------------------------- !
!                                                                           !
!                          Blasius Boundary Layer                           !
!                                                                           !
! ------------------------------------------------------------------------- !

subroutine BlasiusBL(initu, initv, initw)

! Blasius boundary layer profile

! Change the following files:

! 1. ./setup/dns --> BCs, domain, flow speed, re
! 2. ./setup/geometry/Geom.f90 --> suction along top boundary (ctp = 7)

    ! NOTE: I don't think this is working properly yet. Since I haven't used
    ! it, I haven't gone through the trouble of debugging, so good luck

    use omp_lib
    use grid_size
    
    implicit none
    
    real,dimension(nyp,nz,nx) :: initu, initv, initw
    
   
    ! Calculation Variables 
    integer :: i,j,k
    real    :: bconst, eta, U_temp
    
    ! Flow setup variables
	real    :: re,Uinf,R_tau,dPdx

    ! Geometry variables
    integer :: kwall,kmaxsurf
    integer, dimension(nyp,mz,mx) :: imatrix

    ! Buffer variables
    real    :: xstart

    ! --------------------------------------------------------------- !
    !                         Common Blocks                           !
    ! --------------------------------------------------------------- !
    common/flow/       re,Uinf,R_tau,dPdx
    common/imat/       imatrix,kwall,kmaxsurf
    common/buffx/      xstart
    ! --------------------------------------------------------------- !
    ! --------------------------------------------------------------- !
    !                        Set Initial Flow                         !
    ! --------------------------------------------------------------- !
    print *,'Uinf = ',Uinf
   
    ! Initialize
    initu = 0.0
    initv = 0.0
    initw = 0.0

    ! Solve Blasius BL equation 
    if (xstart .le. 0.0001) then
        initu = Uinf ! need to update this -- ??
    else
        bconst = sqrt(re*Uinf/(2*xstart))
        !omp parallel do 
        do i = 1,nyp
            eta = (ycoord(i)-ycoord(nyp))*bconst
            if (eta .le. 8.0) then
                call BlasSolver(eta,U_temp)
            end if
            do j = 1,nz
                do k = 1,nx
                    initu(i,j,k) = U_temp*Uinf
                end do
            end do
        end do
        !omp end parallel do
    end if
!        do k = kwall,(nyp-kwall)
!            ytemp = ycoord(k) - ycoord(kwall)
!            eta = ytemp*bconst
!            if (eta .le. 8.0) then
!                call BlasSolver(eta, U_temp)
!            end if
!            do j = 1,nz
!                do i = 1,nx    
!                    initu(k,j,i) = U_temp * Uinf 
!                end do
!            end do
!        end do
end subroutine

subroutine BlasSolver(eta,u)

    implicit none

    real :: h, y, f, fp, fpp, fppp, u, eta

    h = 1.0e-6
    y = 0.0
    f = 0.0
    fp = 0.0
    fpp = 0.46960
    fppp = (f * fpp)/(-1.0)
    
    do while(y <= eta)
        fpp = fpp + h*fppp
        fp = fp + h*fpp
        f = f + h*fp
        fppp = (f * fpp)/(-1.0)
        y = y + h
    end do
    
    u = fp

end subroutine blassolver 

! ------------------------------------------------------------------------- !
!                                                                           !
!                             Vortex Only                                   !
!                                                                           !
! ------------------------------------------------------------------------- !

subroutine vortex_only(initu,initv,initw)

! This flow is a weak vortex with otherwise still fluid

! vyc and vzc are the center y and z coordinates of the primary vortex,
! respectively. vortY and vortZ are the y and z coordinates of the particular
! vortex being calculated (position changes for images).

! Images are placed above and below the y boundaries for no-flow-through
! conditions to be satisfied, and the primary vortex along with these y images
! are repeated in the positive and negative z-directions (theoretically
! infinitely for the periodic domain, but practically as many as it takes for a
! sufficient solution)

    use grid_size
    
    implicit none

    real,dimension(nyp,nz,nx) :: initu, initv, initw

    ! Vortex variables
    real,dimension(nz,nyp) :: vortVi, vortWi
    real    :: vyc,vzc,vortR2,vSign
    real    :: vortGamma,vortSigma,vortY,vortZ,vortSpace
    integer :: vNum
    
    ! Calculation variables
    real    :: pi, y, z
    integer :: i, j, k
    
    ! Geometry variables
    real :: xl,yl,zl
    
    ! Temporary variables for saving vortY and vortZ (just in case they're used elsewhere)
    real :: tempY, tempZ 
    
    ! ==================================================================== !
    !                        Define common blocks                          !
    ! ==================================================================== !
    
    common/domain/ xl,yl,zl
    common/vort/   vortGamma,vortSigma,vortY,vortZ,vortSpace,vNum

    ! -------------------------------------------------------------------- !

    tempY = vortY
    tempZ = vortZ
    vyc = vortY
    vzc = vortZ
    pi = 2.0*acos(0.0)

    ! initialize vortex velocity arrays
    vortVi = 0.0
    vortWi = 0.0

    ! Add a second vortex if vNum = 2
    do i = 0,vNum-1

        vyc = vyc + i*vortSpace ! Center of vortex
        vSign = float((-1)**i) ! Sign of vortex - default is counter-rotating vortex pair - Ryan 6/8/23

        do j = 1,nz
            z = (j-1)*zl/(nz-1)
            
            do k = 1,nyp
                y = ycoord(k)
    
                ! Primary vortex flowfield
                vortY = y - vyc
                vortZ = z - vzc
                vortR2 = (vortY)**2 + (vortZ)**2
                if (vortR2 .eq. 0.0) then
                    vortR2 = 1.0e-8 ! no NaNs
                end if
                vortVi(j,k) = vortVi(j,k) - vSign*(vortZ)/vortR2*(1 - exp(-vortR2/vortSigma**2))
                vortWi(j,k) = vortWi(j,k) + vSign*(vortY)/vortR2*(1 - exp(-vortR2/vortSigma**2))
   
                ! Testing with 8 image vortices
                ! Right (positive z)
               vortY = y - vyc
               vortZ = z - (vzc + zl) 
               vortR2 = (vortY)**2 + (vortZ)**2
               vortVi(j,k) = vortVi(j,k) - vSign*(vortZ)/vortR2*(1 - exp(-vortR2/vortSigma**2))
               vortWi(j,k) = vortWi(j,k) + vSign*(vortY)/vortR2*(1 - exp(-vortR2/vortSigma**2))
   
               ! Left (negative z)
               vortY = y - vyc
               vortZ = z - (vzc - zl) 
               vortR2 = (vortY)**2 + (vortZ)**2
               vortVi(j,k) = vortVi(j,k) - vSign*(vortZ)/vortR2*(1 - exp(-vortR2/vortSigma**2))
               vortWi(j,k) = vortWi(j,k) + vSign*(vortY)/vortR2*(1 - exp(-vortR2/vortSigma**2))
   
                ! Top (positive y) - flip sign of vortex
                vortY = y - (2*yl - vyc)
                vortZ = z - vzc
                vortR2 = (vortY)**2 + (vortZ)**2
                vortVi(j,k) = vortVi(j,k) + vSign*(vortZ)/vortR2*(1 - exp(-vortR2/vortSigma**2))
                vortWi(j,k) = vortWi(j,k) - vSign*(vortY)/vortR2*(1 - exp(-vortR2/vortSigma**2)) 

                ! Bottom (negative y) - flip sign of vortex
                vortY = y - (-vyc)
                vortZ = z - vzc
                vortR2 = (vortY)**2 + (vortZ)**2
                vortVi(j,k) = vortVi(j,k) + vSign*(vortZ)/vortR2*(1 - exp(-vortR2/vortSigma**2))
                vortWi(j,k) = vortWi(j,k) - vSign*(vortY)/vortR2*(1 - exp(-vortR2/vortSigma**2))
    
                ! Top Right (positive y,z) - flip sign of vortex
                vortY = y - (2*yl - vyc)
                vortZ = z - (vzc + zl) 
                vortR2 = (vortY)**2 + (vortZ)**2
                vortVi(j,k) = vortVi(j,k) + vSign*(vortZ)/vortR2*(1 - exp(-vortR2/vortSigma**2))
                vortWi(j,k) = vortWi(j,k) - vSign*(vortY)/vortR2*(1 - exp(-vortR2/vortSigma**2))
    
                ! Bottom Right (negative y, positive z) - flip sign of vortex
                vortY = y - (-vyc)
                vortZ = z - (vzc + zl) 
                vortR2 = (vortY)**2 + (vortZ)**2
                vortVi(j,k) = vortVi(j,k) + vSign*(vortZ)/vortR2*(1 - exp(-vortR2/vortSigma**2))
                vortWi(j,k) = vortWi(j,k) - vSign*(vortY)/vortR2*(1 - exp(-vortR2/vortSigma**2))
    
                ! Top Left (positive y, negative z) - flip sign of vortex
                vortY = y - (2*yl - vyc)
                vortZ = z - (vzc - zl) 
                vortR2 = (vortY)**2 + (vortZ)**2
                vortVi(j,k) = vortVi(j,k) + vSign*(vortZ)/vortR2*(1 - exp(-vortR2/vortSigma**2))
                vortWi(j,k) = vortWi(j,k) - vSign*(vortY)/vortR2*(1 - exp(-vortR2/vortSigma**2)) 

                ! Bottom Left (negative y,z) - flip sign of vortex
                vortY = y - (-vyc)
                vortZ = z - (vzc - zl) 
                vortR2 = (vortY)**2 + (vortZ)**2
                vortVi(j,k) = vortVi(j,k) + vSign*(vortZ)/vortR2*(1 - exp(-vortR2/vortSigma**2))
                vortWi(j,k) = vortWi(j,k) - vSign*(vortY)/vortR2*(1 - exp(-vortR2/vortSigma**2)) 
            end do
        end do
    end do

    ! Initialize flowfield
    do k = 1, nyp
        do j = 1,nz
            do i = 1,nx
                initu(k,j,i) = 0.0
                initv(k,j,i) = -vortGamma/(2.0*pi)*vortVi(j,k)
                initw(k,j,i) =  vortGamma/(2.0*pi)*vortWi(j,k)
            end do
        end do
    end do

    vorty = tempy
    vortz = tempz

end subroutine

! ------------------------------------------------------------------------- !
!                                                                           !
!                             vortex channel                                !
!                                                                           !
! ------------------------------------------------------------------------- !

subroutine vortex_channel(initu,initv,initw)

! This flow is a steady Channel flow with a weak vortex in the center of the
! channel

! vyc and vzc are the center y and z coordinates of the primary vortex,
! respectively. vortY and vortZ are the y and z coordinates of the particular
! vortex being calculated (position changes for images).

! Images are placed above and below the y boundaries for no-flow-through
! conditions to be satisfied, and the primary vortex along with these y images
! are repeated in the positive and negative z-directions (theoretically
! infinitely for the periodic domain, but practically as many as it takes for a
! sufficient solution)

    use grid_size
    
    implicit none

    real,dimension(nyp,nz,nx) :: initu, initv, initw

    ! Vortex variables
    real,dimension(nz,nyp) :: vortVi, vortWi
    real    :: vyc,vzc,vortR2,vSign
    real    :: vortGamma,vortSigma,vortY,vortZ,vortSpace
    integer :: vNum
    
    ! Calculation variables
    real    :: ReyNo
    real    :: pi, y, z
    integer :: i, j, k
    
    ! Geometry variables
    real :: xl,yl,zl
    
    ! Flow setup variables
    real :: re,Uinf,R_tau,dPdx
   
    ! Temporary variables for saving vortY and vortZ (just in case they're used elsewhere)
    real :: tempY, tempZ 
    
    ! ==================================================================== !
    !                        Define common blocks                          !
    ! ==================================================================== !
    
    common/flow/   re,Uinf,R_tau,dPdx
    common/domain/ xl,yl,zl
    common/vort/   vortGamma,vortSigma,vortY,vortZ,vortSpace,vNum

    ! -------------------------------------------------------------------- !

    tempY = vortY
    tempZ = vortZ
    pi = 2.0*acos(0.0)

    ! Calculate Channel variables
    dPdx = -2/yl*(2*R_tau/(yl*re))**2
    Uinf = -((yl/2)**2)*re/2*dPdx 
    ReyNo = Uinf*yl/2*re ! Calculate channel half height Reynolds number

    write(*,*) "R_tau equals ",R_tau
    write(*,*) "Uinf equals ", Uinf
    write(*,*) "Reynolds Number equals ",ReyNo
    write(*,*) "dPdx equals ",dPdx


    ! initialize vortex velocity arrays
    vortVi = 0.0
    vortWi = 0.0

    ! Add a second vortex if vNum = 2
    do i = 0,vNum-1

        vyc = vyc + i*vortSpace ! Center of vortex
        vSign = float((-1)**i) ! Sign of vortex - default is counter-rotating vortex pair - Ryan 6/8/23

        do j = 1,nz
            z = (j-1)*zl/nz
            
            do k = 1,nyp
                y = ycoord(k)
    
                ! Primary vortex flowfield
                vortY = y - vyc
                vortZ = z - vzc
                vortR2 = (vortY)**2 + (vortZ)**2
                vortVi(j,k) = vortVi(j,k) - vSign*(vortZ)/vortR2*(1 - exp(-vortR2/vortSigma**2))
                vortWi(j,k) = vortWi(j,k) + vSign*(vortY)/vortR2*(1 - exp(-vortR2/vortSigma**2))
    
                ! Testing with 8 image vortices
                ! Right (positive z)
                vortY = y - vyc
                vortZ = z - (vzc + zl) 
                vortR2 = (vortY)**2 + (vortZ)**2
                vortVi(j,k) = vortVi(j,k) - vSign*(vortZ)/vortR2*(1 - exp(-vortR2/vortSigma**2))
                vortWi(j,k) = vortWi(j,k) + vSign*(vortY)/vortR2*(1 - exp(-vortR2/vortSigma**2))
    
                ! Left (negative z)
                vortY = y - vyc
                vortZ = z - (vzc - zl) 
                vortR2 = (vortY)**2 + (vortZ)**2
                vortVi(j,k) = vortVi(j,k) - vSign*(vortZ)/vortR2*(1 - exp(-vortR2/vortSigma**2))
                vortWi(j,k) = vortWi(j,k) + vSign*(vortY)/vortR2*(1 - exp(-vortR2/vortSigma**2))
    
                ! Top (positive y) - flip sign of vortex
                vortY = y - (2*yl - vyc)
                vortZ = z - vzc
                vortR2 = (vortY)**2 + (vortZ)**2
                vortVi(j,k) = vortVi(j,k) + vSign*(vortZ)/vortR2*(1 - exp(-vortR2/vortSigma**2))
                vortWi(j,k) = vortWi(j,k) - vSign*(vortY)/vortR2*(1 - exp(-vortR2/vortSigma**2)) 
                ! Bottom (negative y) - flip sign of vortex
                vortY = y - (-vyc)
                vortZ = z - vzc
                vortR2 = (vortY)**2 + (vortZ)**2
                vortVi(j,k) = vortVi(j,k) + vSign*(vortZ)/vortR2*(1 - exp(-vortR2/vortSigma**2))
                vortWi(j,k) = vortWi(j,k) - vSign*(vortY)/vortR2*(1 - exp(-vortR2/vortSigma**2))
    
                ! Top Right (positive y,z) - flip sign of vortex
                vortY = y - (2*yl - vyc)
                vortZ = z - (vzc + zl) 
                vortR2 = (vortY)**2 + (vortZ)**2
                vortVi(j,k) = vortVi(j,k) + vSign*(vortZ)/vortR2*(1 - exp(-vortR2/vortSigma**2))
                vortWi(j,k) = vortWi(j,k) - vSign*(vortY)/vortR2*(1 - exp(-vortR2/vortSigma**2))
    
                ! Bottom Right (negative y, positive z) - flip sign of vortex
                vortY = y - (-vyc)
                vortZ = z - (vzc + zl) 
                vortR2 = (vortY)**2 + (vortZ)**2
                vortVi(j,k) = vortVi(j,k) + vSign*(vortZ)/vortR2*(1 - exp(-vortR2/vortSigma**2))
                vortWi(j,k) = vortWi(j,k) - vSign*(vortY)/vortR2*(1 - exp(-vortR2/vortSigma**2))
    
                ! Top Left (positive y, negative z) - flip sign of vortex
                vortY = y - (2*yl - vyc)
                vortZ = z - (vzc - zl) 
                vortR2 = (vortY)**2 + (vortZ)**2
                vortVi(j,k) = vortVi(j,k) + vSign*(vortZ)/vortR2*(1 - exp(-vortR2/vortSigma**2))
                vortWi(j,k) = vortWi(j,k) - vSign*(vortY)/vortR2*(1 - exp(-vortR2/vortSigma**2)) 
                ! Bottom Left (negative y,z) - flip sign of vortex
                vortY = y - (-vyc)
                vortZ = z - (vzc - zl) 
                vortR2 = (vortY)**2 + (vortZ)**2
                vortVi(j,k) = vortVi(j,k) + vSign*(vortZ)/vortR2*(1 - exp(-vortR2/vortSigma**2))
                vortWi(j,k) = vortWi(j,k) - vSign*(vortY)/vortR2*(1 - exp(-vortR2/vortSigma**2)) 
            end do
        end do
    end do

    ! Initialize flowfield
    do k = 1, nyp
        do j = 1,nz
            do i = 1,nx
                  initu(k,j,i) = re*dPdx*(0.5*(ycoord(k)-ycoord(nyp))**2 - yl/2*(ycoord(k)-ycoord(nyp)))
                  initv(k,j,i) = -vortGamma/(2.0*pi)*vortVi(j,k)
                  initw(k,j,i) =  vortGamma/(2.0*pi)*vortWi(j,k)
            end do
        end do
    end do

    vortY = tempY
    vortZ = tempZ

end subroutine

! ------------------------------------------------------------------------- !
!                                                                           !
!                             Vortex Couette                                !
!                                                                           !
! ------------------------------------------------------------------------- !

subroutine vortex_Couette(initu,initv,initw)

! This flow is a steady Couette flow with a weak vortex in the center of the
! channel

! vyc and vzc are the center y and z coordinates of the primary vortex,
! respectively. vortY and vortZ are the y and z coordinates of the particular
! vortex being calculated (position changes for images).

! Images are placed above and below the y boundaries for no-flow-through
! conditions to be satisfied, and the primary vortex along with these y images
! are repeated in the positive and negative z-directions (theoretically
! infinitely for the periodic domain, but practically as many as it takes for a
! sufficient solution)

    use grid_size
    
    implicit none

    real,dimension(nyp,nz,nx) :: initu, initv, initw

    ! Vortex variables
    real,dimension(nz,nyp) :: vortVi, vortWi
    real    :: vyc,vzc,vortR2,vSign
    real    :: vortGamma,vortSigma,vortY,vortZ,vortSpace
    integer :: vNum
    
    ! Calculation variables
    real    :: pi, y, z
    integer :: i, j, k
    
    ! Geometry variables
    real :: xl,yl,zl
    
    ! Flow setup variables
    real :: re,Uinf,R_tau,dPdx
   
    ! Temporary variables for saving vortY and vortZ (just in case they're used elsewhere)
    real :: tempY, tempZ 
    
    ! ==================================================================== !
    !                        Define common blocks                          !
    ! ==================================================================== !
    
    common/flow/   re,Uinf,R_tau,dPdx
    common/domain/ xl,yl,zl
    common/vort/   vortGamma,vortSigma,vortY,vortZ,vortSpace,vNum

    ! -------------------------------------------------------------------- !

    tempY = vortY
    tempZ = vortZ
    pi = 2.0*acos(0.0)

    ! initialize vortex velocity arrays
    vortVi = 0.0
    vortWi = 0.0

    ! Add a second vortex if vNum = 2
    do i = 0,vNum-1

        vyc = vyc + i*vortSpace ! Center of vortex
        vSign = float((-1)**i) ! Sign of vortex - default is counter-rotating vortex pair - Ryan 6/8/23

        do j = 1,nz
            z = (j-1)*zl/nz
            
            do k = 1,nyp
                y = ycoord(k)
    
                ! Primary vortex flowfield
                vortY = y - vyc
                vortZ = z - vzc
                vortR2 = (vortY)**2 + (vortZ)**2
                vortVi(j,k) = vortVi(j,k) - vSign*(vortZ)/vortR2*(1 - exp(-vortR2/vortSigma**2))
                vortWi(j,k) = vortWi(j,k) + vSign*(vortY)/vortR2*(1 - exp(-vortR2/vortSigma**2))
    
                ! Testing with 8 image vortices
                ! Right (positive z)
                vortY = y - vyc
                vortZ = z - (vzc + zl) 
                vortR2 = (vortY)**2 + (vortZ)**2
                vortVi(j,k) = vortVi(j,k) - vSign*(vortZ)/vortR2*(1 - exp(-vortR2/vortSigma**2))
                vortWi(j,k) = vortWi(j,k) + vSign*(vortY)/vortR2*(1 - exp(-vortR2/vortSigma**2))
    
                ! Left (negative z)
                vortY = y - vyc
                vortZ = z - (vzc - zl) 
                vortR2 = (vortY)**2 + (vortZ)**2
                vortVi(j,k) = vortVi(j,k) - vSign*(vortZ)/vortR2*(1 - exp(-vortR2/vortSigma**2))
                vortWi(j,k) = vortWi(j,k) + vSign*(vortY)/vortR2*(1 - exp(-vortR2/vortSigma**2))
    
                ! Top (positive y) - flip sign of vortex
                vortY = y - (2*yl - vyc)
                vortZ = z - vzc
                vortR2 = (vortY)**2 + (vortZ)**2
                vortVi(j,k) = vortVi(j,k) + vSign*(vortZ)/vortR2*(1 - exp(-vortR2/vortSigma**2))
                vortWi(j,k) = vortWi(j,k) - vSign*(vortY)/vortR2*(1 - exp(-vortR2/vortSigma**2)) 
                ! Bottom (negative y) - flip sign of vortex
                vortY = y - (-vyc)
                vortZ = z - vzc
                vortR2 = (vortY)**2 + (vortZ)**2
                vortVi(j,k) = vortVi(j,k) + vSign*(vortZ)/vortR2*(1 - exp(-vortR2/vortSigma**2))
                vortWi(j,k) = vortWi(j,k) - vSign*(vortY)/vortR2*(1 - exp(-vortR2/vortSigma**2))
    
                ! Top Right (positive y,z) - flip sign of vortex
                vortY = y - (2*yl - vyc)
                vortZ = z - (vzc + zl) 
                vortR2 = (vortY)**2 + (vortZ)**2
                vortVi(j,k) = vortVi(j,k) + vSign*(vortZ)/vortR2*(1 - exp(-vortR2/vortSigma**2))
                vortWi(j,k) = vortWi(j,k) - vSign*(vortY)/vortR2*(1 - exp(-vortR2/vortSigma**2))
    
                ! Bottom Right (negative y, positive z) - flip sign of vortex
                vortY = y - (-vyc)
                vortZ = z - (vzc + zl) 
                vortR2 = (vortY)**2 + (vortZ)**2
                vortVi(j,k) = vortVi(j,k) + vSign*(vortZ)/vortR2*(1 - exp(-vortR2/vortSigma**2))
                vortWi(j,k) = vortWi(j,k) - vSign*(vortY)/vortR2*(1 - exp(-vortR2/vortSigma**2))
    
                ! Top Left (positive y, negative z) - flip sign of vortex
                vortY = y - (2*yl - vyc)
                vortZ = z - (vzc - zl) 
                vortR2 = (vortY)**2 + (vortZ)**2
                vortVi(j,k) = vortVi(j,k) + vSign*(vortZ)/vortR2*(1 - exp(-vortR2/vortSigma**2))
                vortWi(j,k) = vortWi(j,k) - vSign*(vortY)/vortR2*(1 - exp(-vortR2/vortSigma**2)) 
                ! Bottom Left (negative y,z) - flip sign of vortex
                vortY = y - (-vyc)
                vortZ = z - (vzc - zl) 
                vortR2 = (vortY)**2 + (vortZ)**2
                vortVi(j,k) = vortVi(j,k) + vSign*(vortZ)/vortR2*(1 - exp(-vortR2/vortSigma**2))
                vortWi(j,k) = vortWi(j,k) - vSign*(vortY)/vortR2*(1 - exp(-vortR2/vortSigma**2)) 
            end do
        end do
    end do

    ! Initialize flowfield
    do k = 1, nyp
        do j = 1,nz
            do i = 1,nx
                  initu(k,j,i) = Uinf*(ycoord(k)-ycoord(nyp))/(ycoord(nyp)-ycoord(nyp))
                  initv(k,j,i) = -vortGamma/(2.0*pi)*vortVi(j,k)
                  initw(k,j,i) =  vortGamma/(2.0*pi)*vortWi(j,k)
            end do
        end do
    end do

    vortY = tempY
    vortZ = tempZ

end subroutine
