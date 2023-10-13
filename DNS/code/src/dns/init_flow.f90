subroutine init_flow(initu, initv, initw)

! This set of subroutines initializes the flow field outside the main dns code

! This is the main subroutine which selects the type of flow. If the flow case
! you want to run is already coded, simply change the flow_type selector and
! adjust any appropriate variables defined here

! NOTE: Be sure the geometry and setup files align with the flow type selected

use grid_size

implicit none

real,dimension(nyp,nz,nx) :: initu, initv, initw
integer :: geomtype, flow_select
common/geometry/ geomtype, flow_select


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
  
  real,dimension(nyp,nz,nx) :: initu, initv, initw
  real :: ReyNo, Uinf

  integer :: i, j, k, kwall
  common/kwallpos/ kwall
  
  real :: re, alpha, beta, xl, zl, yl
  common/data2/ re, alpha, beta, xl, zl, yl

  real :: dPdx, R_tau
  common/pressure/ dPdx, R_tau


  dPdx = -2/yl*(2*R_tau/(yl*re))**2
  Uinf = -((yl/2)**2)*re/2*dPdx 
  ReyNo = Uinf*yl/2*re ! Calculate channel half height Reynolds number

  write(*,*) "R_tau equals ",R_tau
  write(*,*) "Uinf equals ", Uinf
  write(*,*) "Reynolds Number equals ",ReyNo
  write(*,*) "dPdx equals ",dPdx

  do k = 1, nyp
    do j = 1,nz
      do i = 1,nx    
         if (k .le. kwall .or. (nyp-k) .le. kwall) then
            initu(k,j,i) = 0.0  
         else   
            initu(k,j,i) = re*dPdx*(0.5*ycoord(k)**2 - yl/2*ycoord(k))
            initv(k,j,i) = 0.0
            initw(k,j,i) = 0.0
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

! Change the following files:

! 1. ./setup/dns --> BCs, domain, flow speed, re
! 2. ./setup/geometry/Geom.f90 --> suction along top boundary (ctp = 7)
  use grid_size
  
  implicit none
  
  real,dimension(nyp,nz,nx) :: initu, initv, initw

  integer :: i,j,k, kwall
  common/kwallpos/ kwall
  
  real :: re, alpha, beta, xl, zl, yl, time, dyde, xstart, Uinf
  common/data2/ re, alpha, beta, xl, zl, yl, time, dyde, xstart, Uinf


  do k = 1, nyp    
    do j = 1,nz
      do i = 1,nx
            initu(k,j,i) = Uinf*(ycoord(k-kwall))/ycoord(nyp - kwall)
            initv(k,j,i) = 0.0
            initw(k,j,i) = 0.0
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

! Change the following files:

! 1. ./setup/dns --> BCs, domain, flow speed, re
! 2. ./setup/geometry/Geom.f90 --> suction along top boundary (ctp = 7)

  use grid_size
  
  implicit none
  
  real,dimension(nyp,nz,nx) :: initu, initv, initw

  integer :: i,j,k, kwall
  common/kwallpos/ kwall
  
  real :: re, alpha, beta, xl, zl, yl, time, dyde, xstart, Uinf
  common/data2/ re, alpha, beta, xl, zl, yl, time, dyde, xstart, Uinf


  do k = 1, nyp  
    do j = 1,nz
      do i = 1,nx  
         if (k .le. kwall) then
            initu(k,j,i) = 0.0  
         else   
            initu(k,j,i) = Uinf*ycoord(k-kwall)/ycoord(nyp-kwall)
            initv(k,j,i) = 0.0
            initw(k,j,i) = Uinf*ycoord(k-kwall)/ycoord(nyp-kwall)
        end if
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


  use omp_lib
  use grid_size
  
  implicit none
  
  real,dimension(nyp,nz,nx) :: initu, initv, initw

  integer :: i,j,k, kwall

  common/kwallpos/ kwall
  
  real :: re, alpha, beta, xl, zl, yl, time, dyde, xstart, Uinf
  common/data2/ re, alpha, beta, xl, zl, yl, time, dyde, xstart, Uinf

  real :: bconst, eta, U_temp, ytemp
 
    print *,'Uinf = ',Uinf
 
  if (xstart .le. 0.0001) then
     initu = Uinf ! need to update this
  else
     bconst = sqrt(re*Uinf/(2*xstart))
     !$omp parallel do default(shared) private(i,j,k,ytemp,eta,U_temp)
     do k = 1, nyp
        print *,'k = ',k
        do j = 1,nz
            do i = 1,nx    
               if (k .le. kwall) then
                  initu(k,j,i) = 0.0  
                  initv(k,j,i) = 0.0
                  initw(k,j,i) = 0.0
               else  
                  ytemp = ycoord(k) - ycoord(kwall)
                  eta = ytemp*bconst
                  if (eta .le. 8.0) then
                     call BlasSolver(eta, U_temp)
                     initu(k,j,i) = U_temp * Uinf 
                     initv(k,j,i) = 0.0
                     initw(k,j,i) = 0.0
                  end if
               end if
            end do
        end do
     end do
     !$omp end parallel do
  end if
end subroutine

subroutine BlasSolver(eta,u)

  implicit none

  real :: h = 1.0/1000000, y, f, fp, fpp, fppp, u, eta

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
!                             Vortex Channel                                !
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

    ! Channel variables
    real :: ReyNo, Uinf

    integer :: i, j, k, kwall
    common/kwallpos/ kwall
    
    real :: re, alpha, beta, xl, zl, yl
    common/data2/ re, alpha, beta, xl, zl, yl

    real :: dPdx, R_tau
    common/pressure/ dPdx, R_tau

    ! Vortex variables
    real :: vyc,vzc
    real :: vortGamma, vortSigma, vortZ, vortY, vortR2
    real :: vortVi(nz,nyp), vortWi(nz,nyp)
    integer :: vNum
    real :: vortSpace

    real :: pi, dij, y, z

    pi = 2.0*acos(0.0)

    ! Calculate Channel variables
    dPdx = -2/yl*(2*R_tau/(yl*re))**2
    Uinf = -((yl/2)**2)*re/2*dPdx 
    ReyNo = Uinf*yl/2*re ! Calculate channel half height Reynolds number

    write(*,*) "R_tau equals ",R_tau
    write(*,*) "Uinf equals ", Uinf
    write(*,*) "Reynolds Number equals ",ReyNo
    write(*,*) "dPdx equals ",dPdx

    ! Read in vortex parameters
    open(999,file='setup/vort.config',status='old')
    do i = 1,14
        read(999,*) 
    end do
    read(999,*) vortGamma, vortSigma
    read(999,*) vyc, vzc
    read(999,*) vNum
    read(999,*) vortSpace
    close(999)



    ! initialize vortex velocity arrays
    do j = 1,nz
        do k = 1,nyp
            vortVi(j,k) = 0.0
            vortWi(j,k) = 0.0
        end do
    end do 

    ! Add a second vortex if vNum = 1
    do i = 0,vNum

    vzc = vzc + i*vortSpace

    do j = 1,nz
        z = (j-1)*zl/nz
        
        do k = 1,nyp
            y = ycoord(k)

            ! Primary vortex flowfield
            vortY = y - vyc
            vortZ = z - vzc
            vortR2 = (vortY)**2 + (vortZ)**2
            vortVi(j,k) = vortVi(j,k) - (vortZ)/vortR2*(1 - exp(-vortR2/vortSigma**2))
            vortWi(j,k) = vortWi(j,k) + (vortY)/vortR2*(1 - exp(-vortR2/vortSigma**2))

            ! Testing with 8 image vortices
            ! Right (positive z)
            vortY = y - vyc
            vortZ = z - (vzc + zl) 
            vortR2 = (vortY)**2 + (vortZ)**2
            vortVi(j,k) = vortVi(j,k) - (vortZ)/vortR2*(1 - exp(-vortR2/vortSigma**2))
            vortWi(j,k) = vortWi(j,k) + (vortY)/vortR2*(1 - exp(-vortR2/vortSigma**2))

            ! Left (negative z)
            vortY = y - vyc
            vortZ = z - (vzc - zl) 
            vortR2 = (vortY)**2 + (vortZ)**2
            vortVi(j,k) = vortVi(j,k) - (vortZ)/vortR2*(1 - exp(-vortR2/vortSigma**2))
            vortWi(j,k) = vortWi(j,k) + (vortY)/vortR2*(1 - exp(-vortR2/vortSigma**2))

            ! Top (positive y) - flip sign of vortex
            vortY = y - (2*yl - vyc)
            vortZ = z - vzc
            vortR2 = (vortY)**2 + (vortZ)**2
            vortVi(j,k) = vortVi(j,k) + (vortZ)/vortR2*(1 - exp(-vortR2/vortSigma**2))
            vortWi(j,k) = vortWi(j,k) - (vortY)/vortR2*(1 - exp(-vortR2/vortSigma**2))

            ! Bottom (negative y) - flip sign of vortex
            vortY = y - (-vyc)
            vortZ = z - vzc
            vortR2 = (vortY)**2 + (vortZ)**2
            vortVi(j,k) = vortVi(j,k) + (vortZ)/vortR2*(1 - exp(-vortR2/vortSigma**2))
            vortWi(j,k) = vortWi(j,k) - (vortY)/vortR2*(1 - exp(-vortR2/vortSigma**2))

            ! Top Right (positive y,z) - flip sign of vortex
            vortY = y - (2*yl - vyc)
            vortZ = z - (vzc + zl) 
            vortR2 = (vortY)**2 + (vortZ)**2
            vortVi(j,k) = vortVi(j,k) + (vortZ)/vortR2*(1 - exp(-vortR2/vortSigma**2))
            vortWi(j,k) = vortWi(j,k) - (vortY)/vortR2*(1 - exp(-vortR2/vortSigma**2))

            ! Bottom Right (negative y, positive z) - flip sign of vortex
            vortY = y - (-vyc)
            vortZ = z - (vzc + zl) 
            vortR2 = (vortY)**2 + (vortZ)**2
            vortVi(j,k) = vortVi(j,k) + (vortZ)/vortR2*(1 - exp(-vortR2/vortSigma**2))
            vortWi(j,k) = vortWi(j,k) - (vortY)/vortR2*(1 - exp(-vortR2/vortSigma**2))

            ! Top Left (positive y, negative z) - flip sign of vortex
            vortY = y - (2*yl - vyc)
            vortZ = z - (vzc - zl) 
            vortR2 = (vortY)**2 + (vortZ)**2
            vortVi(j,k) = vortVi(j,k) + (vortZ)/vortR2*(1 - exp(-vortR2/vortSigma**2))
            vortWi(j,k) = vortWi(j,k) - (vortY)/vortR2*(1 - exp(-vortR2/vortSigma**2))

            ! Bottom Left (negative y,z) - flip sign of vortex
            vortY = y - (-vyc)
            vortZ = z - (vzc - zl) 
            vortR2 = (vortY)**2 + (vortZ)**2
            vortVi(j,k) = vortVi(j,k) + (vortZ)/vortR2*(1 - exp(-vortR2/vortSigma**2))
            vortWi(j,k) = vortWi(j,k) - (vortY)/vortR2*(1 - exp(-vortR2/vortSigma**2))
        end do
    end do
    end do

    ! Initialize flowfield
    do k = 1, nyp
        do j = 1,nz
            do i = 1,nx
                  initu(k,j,i) = re*dPdx*(0.5*ycoord(k)**2 - yl/2*ycoord(k))
!                  initu(k,j,i) = 10.0 ! uniform flow
                  initv(k,j,i) = -vortGamma/(2.0*pi)*vortVi(j,k)
                  initw(k,j,i) =  vortGamma/(2.0*pi)*vortWi(j,k)
            end do
        end do
    end do

end subroutine

! ------------------------------------------------------------------------- !
!                                                                           !
!                             Vortex Couette                                !
!                                                                           !
! ------------------------------------------------------------------------- !

subroutine vortex_couette(initu,initv,initw)

! This is Couette flow with a weak vortex in the center of the
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

    integer :: i, j, k, kwall
    common/kwallpos/ kwall
    
    real :: re, alpha, beta, xl, zl, yl, time, dyde, xstart, Uinf
    common/data2/ re, alpha, beta, xl, zl, yl, time, dyde, xstart, Uinf

    real :: dPdx, R_tau
    common/pressure/ dPdx, R_tau

    ! Vortex variables
    real :: vyc,vzc
    real :: vortGamma, vortSigma, vortZ, vortY, vortR2
    real :: vortVi(nz,nyp), vortWi(nz,nyp)
    integer :: vNum
    real :: vortSpace, vSign

    real :: pi, y, z

    pi = 2.0*acos(0.0)


    ! Read in vortex parameters
    open(999,file='setup/vort.config',status='old')
    do i = 1,14
        read(999,*) 
    end do
    read(999,*) vortGamma, vortSigma
    read(999,*) vyc, vzc
    read(999,*) vNum 
    read(999,*) vortSpace
    close(999)



    ! initialize vortex velocity arrays
    do j = 1,nz
        do k = 1,nyp
            vortVi(j,k) = 0.0
            vortWi(j,k) = 0.0
        end do
    end do 

    ! Add a second vortex if vNum = 1
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
!                  initu(k,j,i) = 10.0 ! uniform flow
                  initu(k,j,i) = Uinf*(ycoord(k-kwall))/ycoord(nyp - kwall)
                  initv(k,j,i) = -vortGamma/(2.0*pi)*vortVi(j,k)
                  initw(k,j,i) =  vortGamma/(2.0*pi)*vortWi(j,k)
            end do
        end do
    end do

end subroutine

! ------------------------------------------------------------------------- !
!                                                                           !
!                             Vortex Only                                   !
!                                                                           !
! ------------------------------------------------------------------------- !

subroutine vortex_only(initu,initv,initw)

! This is a vortex with otherwise still fluid

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

    integer :: i, j, k, kwall
    common/kwallpos/ kwall
    
    real :: re, alpha, beta, xl, zl, yl
    common/data2/ re, alpha, beta, xl, zl, yl

    ! Vortex variables
    real :: vyc,vzc
    real :: vortGamma, vortSigma, vortZ, vortY, vortR2
    real :: vortVi(nz,nyp), vortWi(nz,nyp)
    integer :: vNum
    real :: vortSpace, vSign

    real :: pi, y, z

    pi = 2.0*acos(0.0)


    ! Read in vortex parameters
    open(999,file='setup/vort.config',status='old')
    do i = 1,14
        read(999,*) 
    end do
    read(999,*) vortGamma, vortSigma
    read(999,*) vyc, vzc
    read(999,*) vNum 
    read(999,*) vortSpace
    close(999)



    ! initialize vortex velocity arrays
    do j = 1,nz
        do k = 1,nyp
            vortVi(j,k) = 0.0
            vortWi(j,k) = 0.0
        end do
    end do 

    ! Add a second vortex if vNum = 1
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
                  initu(k,j,i) = 0.0 
                  initv(k,j,i) = -vortGamma/(2.0*pi)*vortVi(j,k)
                  initw(k,j,i) =  vortGamma/(2.0*pi)*vortWi(j,k)
            end do
        end do
    end do

end subroutine


