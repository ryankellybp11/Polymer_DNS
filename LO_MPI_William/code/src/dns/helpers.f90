module helpers
implicit none
! This module contains miscellaneous helper functions used in the dns code
contains
    subroutine yfft(a,wfft1,wfft2,wfft3,wfft4,is)
    !---------------------------------------------------------------------!
    !  Fast Fourier Transform into spectral space:                        !
    !      is = 1  ==> spectral -> physical                               !
    !      is = -1 ==> physical -> spectral                               !
    !---------------------------------------------------------------------!
    
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
    complex, dimension(nyp,nz,nxh) :: a
    real :: wfft1(1),wfft2(1),wfft3(1),wfft4(1)
    integer :: is
    
    ! FFT variables
    real :: trigx(2*nx),trigz(2*nz),trigy(2*ny),sine(ny),cosine(ny)
    integer, dimension(19) :: ixfax,iyfax,izfax,ixfax32,izfax32
    real, dimension(16000) :: trigx32
    real, dimension(4000)  :: trigz32 
    
    ! Calculation variables
    real, dimension(nz) :: sumre,sumim
    real    :: fac
    integer :: i,j,k
    
    !---------------------------------------------------------------------!
    
    ! =================================================================== !
    !                            Common Blocks                            !
    ! =================================================================== !
    common/trig/ trigx,trigy,trigz,sine,cosine,trigx32,trigz32
    common/ifax/ ixfax,iyfax,izfax,ixfax32,izfax32
    !---------------------------------------------------------------------!
    
    ! =================================================================== !
    !                         Begin Calculations                          !
    ! =================================================================== !
    
    ! Physical to spectral
    if (is .eq. -1) then
        fac = 1.0/(2.0*float(ny))
    
        do k = 1,nxh
            do j = 1,nz
                do i = 1,nyp
                    a(i,j,k) = a(i,j,k)*fac
                end do
            end do
        end do
    
        ! y transform
        do k = 1,nxh
            call ccheb(a(1,1,k),wfft1,wfft3,wfft4,wfft2,sumre,sumim,ny,nz,-1,iyfax,trigy,sine,cosine)
        end do
    
    else if (is .eq. 1) then ! Spectral to physical
    
        do k = 1,nxh
            call ccheb(a(1,1,k),wfft1,wfft3,wfft4,wfft2,sumre,sumim,ny,nz,1,iyfax,trigy,sine,cosine)
        end do
    
        fac = 2.0*float(ny)
    
        do k = 1,nxh
            do j = 1,nz
                do i = 1,nyp
                    a(i,j,k) = a(i,j,k)*fac
                end do
            end do
        end do
    
    else ! error check
        write(*,*) ' Error! is = ',is
        stop
    end if
    
    end subroutine yfft
    
    !---------------------------------------------------------------------!
    
    subroutine xyzfft(a,wfft1,wfft2,wfft3,wfft4,is)
    !---------------------------------------------------------------------!
    !  Fast Fourier Transform into spectral space:                        !
    !      is = 1  ==> spectral -> physical                               !
    !      is = -1 ==> physical -> spectral                               !
    !---------------------------------------------------------------------!
    
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
    complex, dimension(nyp,nz,nxh) :: a
    real    :: wfft1(1),wfft2(1),wfft3(1),wfft4(1) ! these were dimension(nmax) in original, but only ever used as a scalar
    integer :: is
    
    ! FFT variables
    real :: trigx(2*nx),trigz(2*nz),trigy(2*ny),sine(ny),cosine(ny)
    integer, dimension(19) :: ixfax,iyfax,izfax,ixfax32,izfax32
    real, dimension(16000) :: trigx32
    real, dimension(4000)  :: trigz32 
    
    ! Calculation variables
    complex, dimension(nxh,nz)     :: c
    complex, dimension(nz,nxh)     :: d
    real, dimension(nz) :: sumre,sumim
    real    :: fac
    integer :: i,j,k
    
    !---------------------------------------------------------------------!
    
    ! =================================================================== !
    !                            Common Blocks                            !
    ! =================================================================== !
    common/trig/       trigx,trigy,trigz,sine,cosine,trigx32,trigz32
    common/ifax/        ixfax,iyfax,izfax,ixfax32,izfax32
    !---------------------------------------------------------------------!
    
    ! =================================================================== !
    !                         Begin Calculations                          !
    ! =================================================================== !
    
    ! Physical to spectral
    if (is .eq. -1) then
    
        ! Do x-z planes, looping on y
        do j = 1,nyp
    
            do k = 1,nz
                do i = 1,nxh
                    c(i,k) = a(j,k,i)
                end do
            end do
    
            ! x transform
            call rcs(c,wfft1,wfft2,nx,nz,ixfax,trigx)
            do i = 1,nxh
                do k = 1,nz
                    d(k,i) = c(i,k)
                end do
            end do
    
            ! z transform
            call ccfft(d,wfft3,wfft4,wfft2,nz,nxh,-1,izfax,trigz)
            fac = 1.0/(2.0*float(ny)*float(nz))
    
    
            do i = 1,nxh
                do k = 1,nz
                    a(j,k,i) = d(k,i)*fac
                end do
            end do
        end do
    
        ! y transform
        do k = 1,nxh
            call ccheb(a(1,1,k),wfft1,wfft2,wfft3,wfft4,sumre,sumim,ny,nz,-1,iyfax,trigy,sine,cosine)
        end do
    
    else if (is .eq. 1) then ! Spectral to physical
    
        do k = 1,nxh
            call ccheb(a(1,1,k),wfft1,wfft2,wfft3,wfft4,sumre,sumim,ny,nz,1,iyfax,trigy,sine,cosine)
        end do
    
        do i = 1,nyp
            fac = 2.0*float(ny)*float(nz)
    
            do k = 1,nxh
                do j = 1,nz
                    d(j,k) = a(i,j,k)*fac
                end do
            end do
    
            call ccfft(d,wfft3,wfft4,wfft2,nz,nxh,1,izfax,trigz)
    
            do k = 1,nxh
                do j = 1,nz
                    c(k,j) = d(j,k)
                end do
            end do
    
            call csr(c,wfft1,wfft2,nx,nz,ixfax,trigx)
    
            do k = 1,nxh
                do j = 1,nz
                    a(i,j,k) = c(k,j)
                end do
            end do
        end do
    
    else ! error check
        write(*,*) ' Error! is = ',is
        stop
    end if
    
    end subroutine xyzfft
    
    !---------------------------------------------------------------------!
    
    subroutine scram(a,b)
    !---------------------------------------------------------------------!
    !  Fills b with data from a                                           !
    !---------------------------------------------------------------------!
    
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
    complex,dimension(nyp,nz,nxh) :: b
    real,   dimension(nyp,nz,nx)  :: a
    
    ! Calculation variables
    integer :: i,j,k
    
    !---------------------------------------------------------------------!
    
    ! =================================================================== !
    !                         Begin Calculations                          !
    ! =================================================================== !
    do i = 1,nyp
        do j = 1,nz
            do k = 1,nxh
                b(i,j,k) = cmplx(a(i,j,k*2-1),a(i,j,k*2))
            end do
        end do
    end do
    
    end subroutine scram
    
    !---------------------------------------------------------------------!
    
    subroutine unscram(b,a)
    !---------------------------------------------------------------------!
    !  Fills a with data from b                                           !
    !---------------------------------------------------------------------!
    
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
    complex,dimension(nyp,nz,nxh) :: b
    real,   dimension(nyp,nz,nx)  :: a
    
    ! Calculation variables
    integer :: i,j,k
    
    !---------------------------------------------------------------------!
    
    ! =================================================================== !
    !                         Begin Calculations                          !
    ! =================================================================== !
    do i = 1,nyp
        do j = 1,nz
            do k = 1,nxh
                a(i,j,2*k-1) = real(b(i,j,k))
                a(i,j,2*k)   = aimag(b(i,j,k))
            end do
        end do
    end do
    
    end subroutine unscram
    
    !---------------------------------------------------------------------!
    
    subroutine norm(a)
    !---------------------------------------------------------------------!
    !  "Normalize" fields, ensure symmetries and reality                  !
    !---------------------------------------------------------------------!
    
    ! =================================================================== !
    !                           Declare Modules                           !
    ! =================================================================== !
    use grid_size
    use omp_lib
    !---------------------------------------------------------------------!
    
    ! =================================================================== !
    !                   Declare all variables explicitly                  !
    ! =================================================================== !
    implicit none
    
    ! Passed variables
    complex, dimension(nyp,nz,nxh) :: a
    
    ! Calculation variables
    integer :: i,j,k,jp
    
    !---------------------------------------------------------------------!
    
    ! =================================================================== !
    !                         Begin Calculations                          !
    ! =================================================================== !
    
    ! k = j = 1 mode
    do i = 1,2
    end do
   
    !$omp parallel do 
    do i = 1,nyp
        a(i,1,1) = real(a(i,1,1))
    end do
    !$omp end parallel do 
    
    ! k = 1, j>1 modes ensure conjugate symmetry
    !$omp parallel do 
    do j = 2,nz/2
        jp = nz + 2 - j
        do i = 1,nyp
            a(i,j,1) = 0.5*(a(i,j,1) + conjg(a(i,jp,1)))
        end do
    end do
    !$omp end parallel do 
    
    !$omp parallel do 
    do j = 2,nz/2
        jp = nz + 2 - j
        do i = 1,nyp
            a(i,jp,1) = conjg(a(i,j,1))
        end do
    end do
    !$omp end parallel do 
    
    ! Zero out highest mode in y-direction, for reality
    !$omp parallel do 
    do k = 1,nxh
        do i = 1,nyp
            a(i,nzhp,k) = 0.0
        end do
    end do
    !$omp end parallel do 

    ! test
    !!$omp parallel do 
    !do k = 1,nxh
    !    do j = 1,nz
    !        do i = 1,nyp
    !            if (abs(real(a(i,j,k))) .lt. 1.0e-16) then
    !                a(i,j,k) = cmplx(0.0,aimag(a(i,j,k)))
    !            end if
    !            if (abs(aimag(a(i,j,k))) .lt. 1.0e-16) then
    !                a(i,j,k) = cmplx(real(a(i,j,k)),0.0)
    !            end if
    !        end do
    !    end do
    !end do
    !!$omp end parallel do 
                
    end subroutine norm
    
    !---------------------------------------------------------------------!
    
    subroutine eig2(u11r,u12r,u21r,u22r,lambda1,lambda2)
    
    ! ==================================================== !
    ! This subroutine is used to get the eigenvalues of a  !
    ! 2x2 matrix, using a pre-determined expression for    !
    ! the values.                                          !
    !                                                      !
    ! The outputs are lambda1, and lambda2, which are the  !
    ! eigenvalues of the matrix U.                         !
    !                                                      !
    ! Created 10/31/2023 by Ryan Kelly                     !
    ! ==================================================== !
    
        implicit none
    
        ! Passed variables  
        real    :: u11r,u12r,u21r,u22r  ! velocity gradients
        complex :: u11,u12,u21,u22  ! velocity gradients
        complex :: lambda1, lambda2,tmp ! Eigenvalues
    

        ! Transfer variables to complex data for complex calculations
        u11 = cmplx(u11r,0.0)
        u12 = cmplx(u12r,0.0)
        u21 = cmplx(u21r,0.0)
        u22 = cmplx(u22r,0.0)
 
        lambda1 = u11/2.0 + u22/2.0 - sqrt(u11**2 - 2*u11*u22 + u22**2 + 4*u12*u21)/2.0
    
        lambda2 = u11/2.0 + u22/2.0 + sqrt(u11**2 - 2*u11*u22 + u22**2 + 4*u12*u21)/2.0
    
    end subroutine eig2
    
    !---------------------------------------------------------------------!
    
    subroutine eig(u11r,u12r,u13r,u21r,u22r,u23r,u31r,u32r,u33r, &
                   lambda1,lambda2,lambda3)
    
    ! ==================================================== !
    ! This subroutine is used to get the eigenvalues of a  !
    ! 3x3 matrix, using a pre-determined expression for    !
    ! the values. 
    !                                                      !
    ! The outputs are lambda1, lambda2, and lambda3, which !
    ! are the eigenvalues of the matrix U.                 !
    !                                                      !
    ! Created 10/26/2022 by Ryan Kelly                     !
    ! ==================================================== !
    
        implicit none
    
        ! Passed variables  
        real    :: u11r,u12r,u13r,u21r,u22r,u23r,u31r,u32r,u33r ! velocity gradients
        complex :: u11,u12,u13,u21,u22,u23,u31,u32,u33 ! velocity gradients
        complex :: lambda1, lambda2, lambda3 ! Eigenvalues
    
        ! Calculation variable
        complex :: im
   
        ! Change to complex data type to avoid NaNs
        u11 = cmplx(u11r,0.0) 
        u12 = cmplx(u12r,0.0) 
        u13 = cmplx(u13r,0.0) 
        u21 = cmplx(u21r,0.0) 
        u22 = cmplx(u22r,0.0) 
        u23 = cmplx(u23r,0.0) 
        u31 = cmplx(u31r,0.0) 
        u32 = cmplx(u32r,0.0) 
        u33 = cmplx(u33r,0.0) 

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
    end subroutine eig
    
    !---------------------------------------------------------------------!
    
    subroutine calcswirl2D(u11,u12,u21,u22,swirl)
    ! This subroutine calculates the swirl strength, defined as the magnitude of the
    ! complex conjugate of the eigenvalues of the velocity gradient tensor
    
    ! Calls subroutine eig(...) to calculate eigenvalues of velocity gradient tensor
    
        use grid_size
    
        implicit none
    
        ! Passed variables  
        real :: u11,u12,u21,u22,swirl
    
        ! Calculation variables
        complex :: l1, l2
        real    :: l1_im, l2_im
   
        print *,'u11 = ',u11 
        print *,'u12 = ',u12 
        print *,'u21 = ',u21 
        print *,'u22 = ',u22
 
        call eig2(u11,u12,u21,u22,l1,l2)

        print *,'l1 = ',l1
        print *,'l2 = ',l2
    
        l1_im = aimag(l1)
        l2_im = aimag(l2)
    
        swirl = sqrt(l1_im**2 + l2_im**2) 
    
        ! Certain velocity gradient tensors produce NaNs using the equations in
        ! eig(...). When this happens, the actual eigenvalues are all 0. This is a
        ! quick check to see if the swirl is a NaN. If it's any real value,
        ! floor(swirl/swirl) will always be 1; otherwise it's a NaN (or Inf, but
        ! that doesn't happen in this case). - Ryan 10/28/22
        if (floor(swirl/swirl) .ne. 1) then
            swirl = 0
        end if
    
    end subroutine calcswirl2D
    
    !---------------------------------------------------------------------!
    
    subroutine calcswirl(u11,u21,u31,u12,u22,u32,u13,u23,u33,swirl)
    ! This subroutine calculates the swirl strength, defined as the magnitude of the
    ! complex conjugate of the eigenvalues of the velocity gradient tensor
    
    ! Calls subroutine eig(...) to calculate eigenvalues of velocity gradient tensor
    
        use grid_size
    
        implicit none
    
        ! Passed variables  
        real :: u11,u12,u13,u21,u22,u23,u31,u32,u33,swirl
    
        ! Calculation variables
        complex :: l1, l2, l3
        real    :: l1_im, l2_im, l3_im
    
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
    
    end subroutine calcswirl
    
    !---------------------------------------------------------------------!
    
    subroutine calc_total_beta(it,delxm,delzm,scp,beta)

        use grid_size
        use omp_lib

        implicit none

        ! Passed variables
        real, dimension(nyp,mz,mx) :: scp,beta
        integer :: it
        real    :: delxm,delzm

        ! Calculation variables
        integer :: i,j,k
        real    :: Lx,Ly,Lz
        real    :: avg_beta,total_scl,volume

        ! Domain variables
        integer :: src_start,src_stop
        real    :: xl,yl,zl

        ! Polymer variables
        integer :: ipolyflag,itarget,ipeter
        real    :: alpha_poly,tpoly,zlmax,diffpoly,qbeta

        !- Common blocks -!
        common/domain/    xl,yl,zl
        common/poly_flgs/ ipolyflag,itarget,ipeter
        common/poly_var/  alpha_poly,tpoly,zlmax,diffpoly,qbeta
        common/src_time/  src_start,src_stop

        !- Begin Calculations -!

        avg_beta = 0.0
        total_scl = 0.0
        volume = xl*yl*zl

        !$omp parallel do reduction(+:avg_beta,total_scl) default(shared) private(i,j,k,Lx,Ly,Lz)
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

                    total_scl = total_scl + scp(i,j,k)*Lx*Ly*Lz
                    avg_beta = avg_beta + beta(i,j,k)*Lx*Ly*Lz
                end do
            end do
        end do
        !$omp end parallel do

        avg_beta = avg_beta/volume
        print *,'avg beta = ',avg_beta
        print *,'total scalar = ',total_scl
        if (avg_beta .le. qbeta) src_stop = it
    end subroutine calc_total_beta
    
    !---------------------------------------------------------------------!

    subroutine polyforce(s11,s12,s13,s22,s23,s33,t1,t2,t3)
    !---------------------------------------------------------------------!
    ! Calculates the divergence of the polymer stress tensor (qp). The    !
    ! resulting force components are stored in t1, t2, and t3.            !
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
    complex, dimension(nyp,nz,nxh) :: s11,s12,s13,s22,s23,s33,t1,t2,t3
    
    ! Calculation variables
    complex, dimension(nyp,nz,nxh) :: wrk1,wrk2,wrk3
    integer :: i,j,k
    complex :: im

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

    ! y-derivatives
    call cderiv(s12,wrk1)
    call cderiv(s22,wrk2)
    call cderiv(s23,wrk3)

    ! Compute div(s)
    do k = 1,nxh
        do j = 1,nz
            do i = 1,nyp
                t1(i,j,k) = im*(wavx(k)*s11(i,j,k) + wavz(j)*s13(i,j,k)) + wrk1(i,j,k) 
                t2(i,j,k) = im*(wavx(k)*s12(i,j,k) + wavz(j)*s23(i,j,k)) + wrk2(i,j,k) 
                t3(i,j,k) = im*(wavx(k)*s13(i,j,k) + wavz(j)*s33(i,j,k)) + wrk3(i,j,k) 
            end do
        end do
    end do

    end subroutine polyforce
    
    !---------------------------------------------------------------------!

    subroutine subforce(gn,fn,omz,t1,t2,t3)
    !---------------------------------------------------------------------!
    !                                                                     !
    !                                                                     !
    !                                                                     !
    !---------------------------------------------------------------------!
    
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
    complex, dimension(nyp,nz,nxh) :: gn,fn,omz,t1,t2,t3
    
    ! Calculation variables
    integer :: i,j,k
    
    !---------------------------------------------------------------------!
    
    ! =================================================================== !
    !                         Begin Calculations                          !
    ! =================================================================== !
    
    do k = 1,nxh
        do j = 1,nz
            do i = 1,nyp
                gn (i,j,k) = gn (i,j,k) + t1(i,j,k)
                fn (i,j,k) = fn (i,j,k) + t2(i,j,k)
                omz(i,j,k) = omz(i,j,k) + t3(i,j,k)
            end do
        end do
    end do

    end subroutine subforce
    
    !---------------------------------------------------------------------!

    subroutine polyNL(c11NL,c12NL,c13NL,c22NL,c23NL,c33NL,         &
                        c11n,c12n,c13n,c22n,c23n,c33n,             &
                        c11nm1,c12nm1,c13nm1,c22nm1,c23nm1,c33nm1, &
                        str11n,str12n,str13n,str22n,str23n,str33n, &
                        str11nm1,str12nm1,str13nm1,str22nm1,str23nm1,str33nm1)
    !---------------------------------------------------------------------!
    ! Calculates the nonlinear terms of the conformation tensor           !
    !---------------------------------------------------------------------!
    
    ! =================================================================== !
    !                           Declare Modules                           !
    ! =================================================================== !
    use grid_size
    use omp_lib
    !---------------------------------------------------------------------!
    
    ! =================================================================== !
    !                   Declare all variables explicitly                  !
    ! =================================================================== !
    implicit none
    
    ! Passed variables
    complex, dimension(nyp,nz,nxh) :: c11n,c12n,c13n,c22n,c23n,c33n
    complex, dimension(nyp,nz,nxh) :: c11nm1,c12nm1,c13nm1,c22nm1,c23nm1,c33nm1
    complex, dimension(nyp,nz,nxh) :: c11NL,c12NL,c13NL,c22NL,c23NL,c33NL
    complex, dimension(nyp,nz,nxh) :: str11n,str12n,str13n,str22n,str23n,str33n
    complex, dimension(nyp,nz,nxh) :: str11nm1,str12nm1,str13nm1,str22nm1,str23nm1,str33nm1

    
    ! Calculation variables
    integer :: i,j,k
    real    :: tpolyinv

    ! Polymer variables
    real    :: alpha_poly,tpoly,zlmax,diffpoly,qbeta

    !---------------------------------------------------------------------!
    
    ! =================================================================== !
    !                            Common Blocks                            !
    ! =================================================================== !
    common/poly_var/  alpha_poly,tpoly,zlmax,diffpoly,qbeta
    !---------------------------------------------------------------------!
    
    ! =================================================================== !
    !                         Begin Calculations                          !
    ! =================================================================== !
    
    tpolyinv = 1.0/tpoly
    do k = 1,nxh
        do j = 1,nz
            do i = 1,nyp
                c11NL(i,j,k) = 3.0*(c11n(i,j,k) - tpolyinv*str11n(i,j,k)) - (c11nm1(i,j,k) - tpolyinv*str11nm1(i,j,k))
                c12NL(i,j,k) = 3.0*(c12n(i,j,k) - tpolyinv*str12n(i,j,k)) - (c12nm1(i,j,k) - tpolyinv*str12nm1(i,j,k))
                c13NL(i,j,k) = 3.0*(c13n(i,j,k) - tpolyinv*str13n(i,j,k)) - (c13nm1(i,j,k) - tpolyinv*str13nm1(i,j,k))
                c22NL(i,j,k) = 3.0*(c22n(i,j,k) - tpolyinv*str22n(i,j,k)) - (c22nm1(i,j,k) - tpolyinv*str22nm1(i,j,k))
                c23NL(i,j,k) = 3.0*(c23n(i,j,k) - tpolyinv*str23n(i,j,k)) - (c23nm1(i,j,k) - tpolyinv*str23nm1(i,j,k))
                c33NL(i,j,k) = 3.0*(c33n(i,j,k) - tpolyinv*str33n(i,j,k)) - (c33nm1(i,j,k) - tpolyinv*str33nm1(i,j,k))
            end do
        end do
    end do

    end subroutine polyNL
    
    !---------------------------------------------------------------------!

    subroutine vortArea(swirl,area)
    !---------------------------------------------------------------------!
    ! Calculate the area inside the vortex core based on swirl > 10       !
    !---------------------------------------------------------------------!
    
    ! =================================================================== !
    !                           Declare Modules                           !
    ! =================================================================== !
    use omp_lib
    use grid_size
    !---------------------------------------------------------------------!
    
    ! =================================================================== !
    !                   Declare all variables explicitly                  !
    ! =================================================================== !
    implicit none
    
    ! Passed variables
    real, dimension(nyp,mz,mx) :: swirl
    real :: area
    
    ! Calculation variables
    integer :: i,j,k
    real    :: delxm, delzm
    real    :: Lx,Ly,Lz

    ! Domain
    real :: xl,yl,zl
    
    common/domain/ xl,yl,zl
    !---------------------------------------------------------------------!
    
    ! =================================================================== !
    !                         Begin Calculations                          !
    ! =================================================================== !
   
    delxm = xl/(mx-1)
    delzm = zl/(mz-1)

    area = 0.0

    !omp parallel do reduction(+:area) default(shared) private(i,j,k,Lx,Ly,Lz) 
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

!                area = area + Ly*Lz*Lx*max((swirl(i,j,k) - 10.0),0.0)/(swirl(i,j,k)-10.0) ! This is a goofy way to avoid an if statement, but it was fun to come up with
                if (swirl(i,j,k) .ge. 10.0) then
                    area = area + Ly*Lz*Lx
                end if
            end do
        end do
    end do
    !omp end parallel do

    if (nx .lt. 10) then
        area = area/xl
    else if (nz .lt. 10) then
        area = area/zl
    end if


    end subroutine vortArea
    
    !---------------------------------------------------------------------!
    
    subroutine calcQ(u11,u21,u31,u12,u22,u32,u13,u23,u33,Q)
    ! This subroutine calculates the Q criterion, which is the local balance
    ! of rotation and strain rate. 

    ! Q > 0 ---> Vortex core (when coupled with local pressure minimum)
    ! Q < 0 ---> Strain-dominant field (presumably what we want to target for polymer drag reduction)
    
        use grid_size
    
        implicit none
    
        ! Passed variables  
        real :: u11,u12,u13,u21,u22,u23,u31,u32,u33,Q
    
        ! Calculation variables
        real :: Rot, Strain
 
        ! Begin calculations
!        Rot = u11**2 + u22**2 + u33**2 ! Rotation tensor magnitude squared
!        Strain = u12*u21 + u13*u31 + u23*u32 ! Strain rate tensor magnitude squared

        Strain = u11**2 + u22**2 + u33**2
        Rot = -2.0*(u12*u21 + u13*u32 + u23*u32)
        Q = 0.5*(Rot - Strain)
         
    end subroutine calcQ
    
    !---------------------------------------------------------------------!

    ! Sorting algorithm taken from Wikipedia
    SUBROUTINE Shell_Sort1(a)
    
      IMPLICIT NONE
      INTEGER :: i, j, increment
      REAL :: temp1,temp2,temp3,temp4
      REAL, INTENT(in out) :: a(:)
        
      increment = SIZE(a) / 2
      DO WHILE (increment > 0)
          DO i = increment+1, SIZE(a)
             j = i
             temp1 = a(i)
             DO WHILE (j >= increment+1 .AND. a(j-increment) > temp1)
                a(j) = a(j-increment)
                j = j - increment
             END DO
             a(j) = temp1
          END DO

          IF (increment == 2) THEN
              increment = 1
          ELSE
              increment = increment * 5 / 11
          END IF      
      END DO
     
    END SUBROUTINE Shell_Sort1
     
    !---------------------------------------------------------------------!

    ! Sorting algorithm taken from Wikipedia - written for multiple variables
    SUBROUTINE Shell_Sort2(a,b)
    
      IMPLICIT NONE
      INTEGER :: i, j, increment
      REAL :: temp1,temp2,temp3,temp4
      integer, INTENT(in out) :: a(:)
      REAL, INTENT(in out) :: b(:)
        
      increment = SIZE(a) / 2
      DO WHILE (increment > 0)
          DO i = increment+1, SIZE(a)
             j = i
             temp1 = a(i)
             temp2 = b(i)
             DO WHILE (j >= increment+1 .AND. a(j-increment) > temp1)
                a(j) = a(j-increment)
                b(j) = b(j-increment)

                j = j - increment
             END DO
             a(j) = temp1
             b(j) = temp2
          END DO

          IF (increment == 2) THEN
              increment = 1
          ELSE
              increment = increment * 5 / 11
          END IF      
      END DO
     
    END SUBROUTINE Shell_Sort2
     
    !---------------------------------------------------------------------!

    ! Sorting algorithm taken from Wikipedia - written for multiple variables
    SUBROUTINE Shell_Sort(a,b,c,d)
    
      IMPLICIT NONE
      INTEGER :: i, j, increment
      REAL :: temp1,temp2,temp3,temp4
      REAL, INTENT(in out) :: a(:)
      INTEGER, INTENT(in out) :: b(:),c(:),d(:)
        
      increment = SIZE(a) / 2
      DO WHILE (increment > 0)
          DO i = increment+1, SIZE(a)
             j = i
             temp1 = a(i)
             temp2 = b(i)
             temp3 = c(i)
             temp4 = d(i)
             DO WHILE (j >= increment+1 .AND. a(j-increment) > temp1)
                a(j) = a(j-increment)
                b(j) = b(j-increment)
                c(j) = c(j-increment)
                d(j) = d(j-increment)

                j = j - increment
             END DO
             a(j) = temp1
             b(j) = temp2
             c(j) = temp3
             d(j) = temp4
          END DO

          IF (increment == 2) THEN
              increment = 1
          ELSE
              increment = increment * 5 / 11
          END IF      
      END DO
     
    END SUBROUTINE Shell_Sort    
     
    !---------------------------------------------------------------------!

    subroutine findmaxQ(u11,u12,u13,u21,u22,u23,u31,u32,u33,scalar,Qmin,Qx,Qy,Qz,beta)
    ! This subroutine is specifically for targeting regions of low Q/high strain.
    ! I suppose it could be easily modified to target other regions without particles,
    ! but I'll worry about that later.

    ! This takes the velocity gradient tensor and scalar field as inputs and outputs
    ! the grid indices of the minimum Q coordinates and the corresponding beta field
    ! for checking local concentration. These quantities must be in 3/2 interpolated 
    ! physical space for whatever reason that I don't really know, but that's why all
    ! these FFTs are necessary.

    ! Transform velocity gradient terms to 3/2 physical space 

    use grid_size
    use omp_lib
    use derivs

    !---------------------------------------------------------------------!
    
    ! =================================================================== !
    !                   Declare all variables explicitly                  !
    ! =================================================================== !

    implicit none

    ! Inputs  
    complex, dimension(nyp,nz,nxh), intent(in) :: u11,u12,u13,u21,u22,u23,u31,u32,u33,scalar

    ! Outputs
    integer, intent(in out)  :: Qx(:),Qy(:),Qz(:)
    real, intent(in out)     :: Qmin(:)
    real, dimension(nyp,mz,mx) :: beta

    ! Intermediate physical y-plane variables
    real, dimension(mzp,mxp2)  :: u11p,u12p,u13p
    real, dimension(mzp,mxp2)  :: u21p,u22p,u23p
    real, dimension(mzp,mxp2)  :: u31p,u32p,u33p
    real, dimension(mzp,mxp2)  :: scp,beta_poly

    ! Intermediate 3D physical variables
    real, dimension(nyp,mz,mx) :: Qcrit
    real, dimension(nyp,mz,mx) :: u11p3d,u12p3d,u13p3d
    real, dimension(nyp,mz,mx) :: u21p3d,u22p3d,u23p3d
    real, dimension(nyp,mz,mx) :: u31p3d,u32p3d,u33p3d

    ! Calculation Variables
    integer :: i,j,k,i1,i2,jj,inc,isgn,jump,lot,cnt
    real    :: QQ
 
    ! FFT variables
    real :: trigx(2*nx),trigz(2*nz),trigy(2*ny),sine(ny),cosine(ny)
    integer, dimension(19) :: ixfax,iyfax,izfax,ixfax32,izfax32
    real, dimension(16000) :: trigx32
    real, dimension(4000)  :: trigz32 
    real,dimension(nwrk)  :: wrk

    ! Polymer variables
    real    :: alpha_poly
    !---------------------------------------------------------------------!

    ! =================================================================== !
    !                            Common Blocks                            !
    ! =================================================================== !
    common/trig/       trigx,trigy,trigz,sine,cosine,trigx32,trigz32
    common/ifax/       ixfax,iyfax,izfax,ixfax32,izfax32
    common/poly_var/   alpha_poly
    !---------------------------------------------------------------------!
    
    ! =================================================================== !
    !                         Begin Calculations                          !
    ! =================================================================== !
    
    !omp parallel do
    do k = 1,nyp
        do j = 1,mzp
            do i = 1,mxp2

                u11p(j,i) = 0.0
                u12p(j,i) = 0.0
                u13p(j,i) = 0.0
                u21p(j,i) = 0.0
                u22p(j,i) = 0.0
                u23p(j,i) = 0.0
                u31p(j,i) = 0.0
                u32p(j,i) = 0.0
                u33p(j,i) = 0.0

                 scp(j,i) = 0.0

            end do ! i
        end do ! j

        do j = 1,nz
            if (j .le. nzh) jj = j
            if (j .gt. nzh) jj = (mz-nz) + j

            do i = 1,nxh
                i1 = 2*(i-1) + 1
                i2 = 2*i

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
                
                scp(jj,i1) =  real(scalar(k,j,i))
                scp(jj,i2) = aimag(scalar(k,j,i))
            end do
        end do
   
    !---------------------------------------------------------------------!
    !        Transform (interpolate) into 3/2 grid physical space         !
    !---------------------------------------------------------------------!
    
        inc  = 1
        isgn = 1
        jump = 2*mzp
        lot  = nx/2
    
        call cfftmlt(u11p(1,1),u11p(1,2),wrk,trigz32,izfax32,inc,jump,mz,lot,isgn)
        call cfftmlt(u12p(1,1),u12p(1,2),wrk,trigz32,izfax32,inc,jump,mz,lot,isgn)
        call cfftmlt(u13p(1,1),u13p(1,2),wrk,trigz32,izfax32,inc,jump,mz,lot,isgn)
        call cfftmlt(u21p(1,1),u21p(1,2),wrk,trigz32,izfax32,inc,jump,mz,lot,isgn)
        call cfftmlt(u22p(1,1),u22p(1,2),wrk,trigz32,izfax32,inc,jump,mz,lot,isgn)
        call cfftmlt(u23p(1,1),u23p(1,2),wrk,trigz32,izfax32,inc,jump,mz,lot,isgn)
        call cfftmlt(u31p(1,1),u31p(1,2),wrk,trigz32,izfax32,inc,jump,mz,lot,isgn)
        call cfftmlt(u32p(1,1),u32p(1,2),wrk,trigz32,izfax32,inc,jump,mz,lot,isgn)
        call cfftmlt(u33p(1,1),u33p(1,2),wrk,trigz32,izfax32,inc,jump,mz,lot,isgn)
    
        call cfftmlt(scp(1,1),scp(1,2),wrk,trigz32,izfax32,inc,jump,mz,lot,isgn)

        do j = 1,mz

            u11p(j,nxp) = u11p(j,2)
            u12p(j,nxp) = u12p(j,2)
            u13p(j,nxp) = u13p(j,2)
            u21p(j,nxp) = u21p(j,2)
            u22p(j,nxp) = u22p(j,2)
            u23p(j,nxp) = u23p(j,2)
            u31p(j,nxp) = u31p(j,2)
            u32p(j,nxp) = u32p(j,2)
            u33p(j,nxp) = u33p(j,2)
    
            scp(j,nxp) = scp(j,2)
            
            u11p(j,2) = 0.0
            u12p(j,2) = 0.0
            u13p(j,2) = 0.0
            u21p(j,2) = 0.0
            u22p(j,2) = 0.0
            u23p(j,2) = 0.0
            u31p(j,2) = 0.0
            u32p(j,2) = 0.0
            u33p(j,2) = 0.0
            
            scp(j,2) = 0.0
     
            u11p(j,nxp2) = u11p(j,2)
            u12p(j,nxp2) = u12p(j,2)
            u13p(j,nxp2) = u13p(j,2)
            u21p(j,nxp2) = u21p(j,2)
            u22p(j,nxp2) = u22p(j,2)
            u23p(j,nxp2) = u23p(j,2)
            u31p(j,nxp2) = u31p(j,2)
            u32p(j,nxp2) = u32p(j,2)
            u33p(j,nxp2) = u33p(j,2)
            
            scp(j,nxp2) = scp(j,2)
        end do

        isgn = 1
        inc  = mzp
        jump = 1
        lot  = mz
    
        call rfftmlt(u11p,wrk,trigx32,ixfax32,inc,jump,mx,lot,isgn)
        call rfftmlt(u12p,wrk,trigx32,ixfax32,inc,jump,mx,lot,isgn)
        call rfftmlt(u13p,wrk,trigx32,ixfax32,inc,jump,mx,lot,isgn)
        call rfftmlt(u21p,wrk,trigx32,ixfax32,inc,jump,mx,lot,isgn)
        call rfftmlt(u22p,wrk,trigx32,ixfax32,inc,jump,mx,lot,isgn)
        call rfftmlt(u23p,wrk,trigx32,ixfax32,inc,jump,mx,lot,isgn)
        call rfftmlt(u31p,wrk,trigx32,ixfax32,inc,jump,mx,lot,isgn)
        call rfftmlt(u32p,wrk,trigx32,ixfax32,inc,jump,mx,lot,isgn)
        call rfftmlt(u33p,wrk,trigx32,ixfax32,inc,jump,mx,lot,isgn)
     
        call rfftmlt(scp,wrk,trigx32,ixfax32,inc,jump,mx,lot,isgn)

    !---------------------------------------------------------------------!
    !       Calculate beta then store it and grad(V) in 3D variables      !
    !---------------------------------------------------------------------!
   
        do j = 1,mz
            do i = 1,mx
                beta_poly(j,i) = exp(-alpha_poly*abs(scp(j,i)))

                beta(k,j,i) = beta_poly(j,i)
                u11p3d(k,j,i) = u11p(j,i) 
                u12p3d(k,j,i) = -u12p(j,i) 
                u13p3d(k,j,i) = u13p(j,i) 
                u21p3d(k,j,i) = -u21p(j,i) 
                u22p3d(k,j,i) = u22p(j,i) 
                u23p3d(k,j,i) = -u23p(j,i) 
                u31p3d(k,j,i) = u31p(j,i) 
                u32p3d(k,j,i) = -u32p(j,i) 
                u33p3d(k,j,i) = u33p(j,i) 
            end do
        end do
    end do ! k
    !omp end parallel do

    !---------------------------------------------------------------------!
    !                        Calculate Q-criterion                        !
    !---------------------------------------------------------------------!

    !$omp parallel do private(i,j,k,QQ)
    do i = 1,nyp
        do j = 1,mz
            do k = 1,mx
                call calcQ(u11p3d(i,j,k),u21p3d(i,j,k),u31p3d(i,j,k),u12p3d(i,j,k),u22p3d(i,j,k), &
                               u32p3d(i,j,k),u13p3d(i,j,k),u23p3d(i,j,k),u33p3d(i,j,k),QQ)
    
                Qcrit(i,j,k) = QQ
            end do
        end do
    end do
    !$omp end parallel do

    ! Check for local (grid) min
    cnt = 1
    do k = 2,mx-1
        do j = 2,mz-1
            do i = 2,ny
                if (Qcrit(i,j,k) .lt. 0 .and. cnt .le. size(Qmin)) then
                    if (Qcrit(i,j,k) .lt. Qcrit(i+1,j,k) .and. Qcrit(i,j,k) .lt. Qcrit(i-1,j,k)) then
                        if (Qcrit(i,j,k) .lt. Qcrit(i,j+1,k) .and. Qcrit(i,j,k) .lt. Qcrit(i,j-1,k)) then
                            if (Qcrit(i,j,k) .lt. Qcrit(i,j,k+1) .and. Qcrit(i,j,k) .lt. Qcrit(i,j,k-1)) then
                                Qmin(cnt) = Qcrit(i,j,k)
                                Qx(cnt) = k
                                Qy(cnt) = i
                                Qz(cnt) = j
                                cnt = cnt + 1
                            end if
                        end if
                    end if
                end if
            end do
        end do
    end do

    call shell_sort(Qmin,Qx,Qy,Qz)

    end subroutine findmaxQ

    !---------------------------------------------------------------------!

    subroutine calcEnergySpectrum(v,smth_flag)
    ! This subroutine takes a y-physical, x-z spectral variable, v, and 
    ! computes the energy spectrum, printed to a file

    use grid_size
    use omp_lib

    !---------------------------------------------------------------------!
    
    ! =================================================================== !
    !                   Declare all variables explicitly                  !
    ! =================================================================== !

    implicit none

    ! Flow variable (input)
    complex, dimension(nyp,nz,nxh), intent(in) :: v

    ! Smoothing flag (output)
    logical, intent(out) :: smth_flag

    ! Wavenumbers
    real,dimension(nxh)        :: wavx
    real,dimension(nz)         :: wavz

    ! Calculation variables
    real,dimension(nxh*nz)    :: vE_sum
    integer,dimension(nxh*nz) :: wn2
    real,dimension(143)       :: tmp ! the size of 143 is dependent on the wavenumbers, and by extension the domain size
    real,dimension(142)       :: tmp2 
    integer :: i,j,k,cnt
    real    :: dx,dy

    !---------------------------------------------------------------------!

    ! =================================================================== !
    !                            Common Blocks                            !
    ! =================================================================== !
    common/waves/     wavx,wavz
    !---------------------------------------------------------------------!
    
    ! =================================================================== !
    !                         Begin Calculations                          !
    ! =================================================================== !
   
    cnt = 1 
    vE_sum = 0.0
    do k = 1,nxh
        do j = 1,nz
            wn2(cnt) = floor(sqrt(wavx(k)**2 + wavz(j)**2))
            do i = 1,nyp
                vE_sum(cnt) = vE_sum(cnt) + (v(i,j,k)*conjg(v(i,j,k)))*0.5
            end do
            cnt = cnt + 1
        end do
    end do

    call Shell_Sort2(wn2,vE_sum)

    ! Initialize arrays
    cnt = 1
    tmp = 0.0

    ! Set first term
    tmp(1) = vE_sum(1)

    ! If duplicate wave number, add to same array element
    do i = 2,nxh*nz
        if (wn2(i) .ne. wn2(i-1)) then
            cnt = cnt + 1
        else
            tmp(cnt) = tmp(cnt) + vE_sum(i)
        end if
    end do

    ! Average values
    do i = 1,142
        tmp2(i) = 0.5*(tmp(i) + tmp(i+1))
    end do    

    open(1,file = 'vE.dat',status = 'replace')
    write(1,*) tmp
    close(1)

    ! Check if we need smoothing
    dy = log(tmp2(100)) - log(tmp2(50))
    dx = log(100.5) - log(50.5)

    print *,'slope = ',dy/dx

    smth_flag = .false.
    if (dy/dx .gt. -3.0) then
        smth_flag = .true.
        print *,'smth_flag = ',smth_flag
    end if


    end subroutine calcEnergySpectrum

    !---------------------------------------------------------------------!

    subroutine calcEnergySpectrum3(u,v,w,smth_flag)
    ! This subroutine takes a y-physical, x-z spectral variable, v, and 
    ! computes the energy spectrum, printed to a file

    use grid_size
    use omp_lib

    !---------------------------------------------------------------------!
    
    ! =================================================================== !
    !                   Declare all variables explicitly                  !
    ! =================================================================== !

    implicit none

    ! Flow variable (input)
    complex, dimension(nyp,nz,nxh), intent(in) :: u,v,w

    ! Wavenumbers
    real,dimension(nxh)        :: wavx
    real,dimension(nz)         :: wavz

    ! Calculation variables
    real,dimension(nxh*nz)    :: vE_sum
    integer,dimension(nxh*nz) :: wn2
    real,dimension(143)       :: tmp ! the size of 143 is dependent on the wavenumbers, and by extension the domain size
    real,dimension(142)       :: tmp2 
    integer :: i,j,k,cnt
    real    :: dx,dy
    logical :: smth_flag
!    real    :: wn2

    !---------------------------------------------------------------------!

    ! =================================================================== !
    !                            Common Blocks                            !
    ! =================================================================== !
    common/waves/     wavx,wavz
    !---------------------------------------------------------------------!
    
    ! =================================================================== !
    !                         Begin Calculations                          !
    ! =================================================================== !
   
    cnt = 1 
    vE_sum = 0.0
    do k = 1,nxh
        do j = 1,nz
            wn2(cnt) = floor(sqrt(wavx(k)**2 + wavz(j)**2))
            do i = 1,nyp
                vE_sum(cnt) = vE_sum(cnt) + 0.5*(u(i,j,k)*conjg(u(i,j,k)) + v(i,j,k)*conjg(v(i,j,k)) + w(i,j,k)*conjg(w(i,j,k)))
            end do
            cnt = cnt + 1
        end do
    end do

    call Shell_Sort2(wn2,vE_sum)

    ! Initialize arrays
    cnt = 1
    tmp = 0.0

    ! Set first term
    tmp(1) = vE_sum(1)

    ! If duplicate wave number, add to same array element
    do i = 2,nxh*nz
        if (wn2(i) .ne. wn2(i-1)) then
            cnt = cnt + 1
        else
            tmp(cnt) = tmp(cnt) + vE_sum(i)
        end if
    end do

    ! Average values
    do i = 1,142
        tmp2(i) = 0.5*(tmp(i) + tmp(i+1))
    end do    

    open(1,file = 'vE.dat',status = 'replace')
    write(1,*) tmp
    close(1)

    ! Check if we need smoothing
    dy = log(tmp2(100)) - log(tmp2(50))
    dx = log(100.5) - log(50.5)

    print *,'slope = ',dy/dx

    smth_flag = .false.
    if (dy/dx .gt. -3.0) then
        smth_flag = .true.
        print *,'smth_flag = ',smth_flag
    end if

    end subroutine calcEnergySpectrum3

    !---------------------------------------------------------------------!

    subroutine calc_enstrophy_terms(omx,torqx,rank)

    use grid_size
    use omp_lib
    use derivs
        
    implicit none

    ! Input variables
    complex, dimension(nyp,nz,nxh) :: omx,torqx,domxdy,domxdz

    ! Calculation variables
    integer :: i,j,k,k1,k2,jj,is,inc,isgn,jump,lot,n
    real    :: poly_ens,flow_ens,Lx,Ly,Lz
    complex :: im = (0.0,1.0)

    ! Physical Space Variables
    real, dimension(mzp,mxp2) :: tx,wx,dwxdy,dwxdz,wrk


    ! FFT Variables
    real, dimension(nmax) :: wfft1,wfft2,wfft3,wfft4 
    real :: trigx(2*nx),trigz(2*nz),trigy(2*ny),sine(ny),cosine(ny)
    integer, dimension(19) :: ixfax,iyfax,izfax,ixfax32,izfax32
    real, dimension(16000) :: trigx32
    real, dimension(4000)  :: trigz32 

    integer :: rank   
    character(10) :: syscommand

    ! Common Blocks
    integer :: irstrt
    real    :: xl,yl,zl 
    real    :: wavx(nxh),wavz(nz) 
    real    :: c(nyp)
    real    :: re
    integer :: it
    real    :: dt

    common/iocontrl/ irstrt
    common/domain/   xl,yl,zl
    common/waves/    wavx,wavz,c
    common/trig/     trigx,trigy,trigz,sine,cosine,trigx32,trigz32
    common/ifax/     ixfax,iyfax,izfax,ixfax32,izfax32
    common/flow/     re
    common/itime/    it,dt

    !! Begin Calculations !!
    ! Compute vorticity derivatives
    call cderiv(omx,domxdy) ! d(wx)/dy

    do k = 1,nxh
        do j = 1,nz
            do i = 1,nyp
                domxdz(i,j,k) = im*wavz(j)*omx(i,j,k) ! d(wx)/dz
            end do
        end do
    end do

    ! Convert variables to physical space for manipulation
    is = 1
  
    ! c2c y-transform 
    call yfft(omx,wfft1,wfft2,wfft3,wfft4,is)
    call yfft(domxdy,wfft1,wfft2,wfft3,wfft4,is)
    call yfft(domxdz,wfft1,wfft2,wfft3,wfft4,is)
    call yfft(torqx,wfft1,wfft2,wfft3,wfft4,is)

    poly_ens = 0.0; flow_ens = 0.0

    !$omp parallel do default(shared) private(i,j,k,k1,k2,jj,inc,isgn,jump,lot,Lx,Ly,Lz, &
    !$omp   wx,dwxdy,dwxdz,tx,wrk) reduction(+:flow_ens,poly_ens) schedule(auto)
    do i = 1,nyp
        do k = 1,mxp2
            do j = 1,mzp
                wx(j,k)    = 0.0
                dwxdy(j,k) = 0.0
                dwxdz(j,k) = 0.0
                tx(j,k)    = 0.0
            end do
        end do

        do j = 1,nz
            if (j .le. nzh) jj = j
            if (j .gt. nzh) jj = (mz-nz) + j

            do k = 1,nxh
                k1 = 2*(k-1) + 1
                k2 = 2*k
                wx(jj,k1) =  real(omx(i,j,k))
                wx(jj,k2) = aimag(omx(i,j,k))
                dwxdy(jj,k1) =  real(domxdy(i,j,k))
                dwxdy(jj,k2) = aimag(domxdy(i,j,k))
                dwxdz(jj,k1) =  real(domxdz(i,j,k))
                dwxdz(jj,k2) = aimag(domxdz(i,j,k))
                tx(jj,k1) =  real(torqx(i,j,k))
                tx(jj,k2) = aimag(torqx(i,j,k))
            end do
        end do

        ! c2r split DFT (z-direction)
        inc  = 1
        isgn = 1
        jump = 2*mzp
        lot  = nx/2
    
        call cfftmlt(wx(1,1),wx(1,2),wrk,trigz32,izfax32,inc,jump,mz,lot,isgn)
        call cfftmlt(dwxdy(1,1),dwxdy(1,2),wrk,trigz32,izfax32,inc,jump,mz,lot,isgn)
        call cfftmlt(dwxdz(1,1),dwxdz(1,2),wrk,trigz32,izfax32,inc,jump,mz,lot,isgn)
        
        call cfftmlt(tx(1,1),tx(1,2),wrk,trigz32,izfax32,inc,jump,mz,lot,isgn)

        do j = 1,mz
            wx(j,nxp) = wx(j,2)        
            dwxdy(j,nxp) = dwxdy(j,2)        
            dwxdz(j,nxp) = dwxdz(j,2)        
            tx(j,nxp) = tx(j,2)        

            wx(j,2) = 0.0
            dwxdy(j,2) = 0.0
            dwxdz(j,2) = 0.0
            tx(j,2) = 0.0

            wx(j,nxp2) = wx(j,2)        
            dwxdy(j,nxp2) = dwxdy(j,2)        
            dwxdz(j,nxp2) = dwxdz(j,2)        
            tx(j,nxp2) = tx(j,2)        
        end do

        ! r2r DFT (x-direction)
        isgn = 1
        inc  = mzp
        jump = 1
        lot  = mz
    
        call rfftmlt(wx,wrk,trigx32,ixfax32,inc,jump,mx,lot,isgn)
        call rfftmlt(dwxdy,wrk,trigx32,ixfax32,inc,jump,mx,lot,isgn)
        call rfftmlt(dwxdz,wrk,trigx32,ixfax32,inc,jump,mx,lot,isgn)
        call rfftmlt(tx,wrk,trigx32,ixfax32,inc,jump,mx,lot,isgn)

        ! Perform physical space calculations
        Lx = 1.0!/mx
        Ly = seght(i)
        Lz = zl/mz
        do k = 1,1!mx
            do j = 1,mz
                poly_ens = poly_ens + 2.0*wx(j,k)*tx(j,k)*Lx*Ly*Lz
                flow_ens = flow_ens - (2.0/re)*(dwxdy(j,k)**2 + dwxdz(j,k)**2)*Lx*Ly*Lz
            end do
        end do

    end do
    !$omp end parallel do

    ! Write outputs to file
    write(syscommand,'("outputs",i2.2,"/")') rank
    if (it .eq. irstrt) then
        open(171,file=syscommand//'dEdt')
    else
        open(171,file=syscommand//'dEdt',position = 'append')
    end if

    write(171,*) flow_ens*dt,poly_ens*dt
    close(171)

    end subroutine calc_enstrophy_terms

end module helpers
