module helpers
implicit none
! This module contains miscellaneous helper functions used in the dns code
contains

    subroutine xyzfft(a,b,is)
    ! Plans and performs a full 3D transform (3-steps) on a either 
    !   is = -1: physical --> spectral
    !   is = +1: spectral --> physical

    use,intrinsic :: iso_c_binding
    use grid_size
    use omp_lib

    implicit none

    include 'fftw3.f03'

    complex(C_DOUBLE_COMPLEX),dimension(nyp,nz,nxh) :: a
    complex(C_DOUBLE_COMPLEX),dimension(nyp,nz,nx) :: aplan
    real(C_DOUBLE),dimension(nyp,nz,nx) :: b,bplan

    integer :: i,j,k,is

    type(C_PTR) :: plan1,plan2,plan3

    real,dimension(nxh) :: wavx
    real,dimension(nz)  :: wavz
    real,dimension(nyp) :: c

    common/waves/ wavx,wavz,c

! ------------------------------------------------------------------------------- !

    ! NOTE: Because plans are created inside the subroutine, it is impossible to execute
    ! this subroutine in parallel because only fftw_execute functions are thread-safe.
    ! This problem can be circumvented by creating plans in the main program and passing them
    ! as arguments to the subroutine, but this works and it's only used inside the initial(...)
    ! subroutine so it has little impact on long runs

        plan1 = fftw_plan_r2r_1d(nyp,bplan,bplan,FFTW_REDFT00,FFTW_ESTIMATE)
        plan2 = fftw_plan_dft_r2c_1d(nx,bplan,aplan,FFTW_ESTIMATE)
        plan3 = fftw_plan_dft_1d(nz,aplan,aplan,is,FFTW_ESTIMATE)

        ! y-transform 
        !$omp parallel do shared(b,plan1) private(j,k) collapse(2)
        do k = 1,nx
            do j = 1,nz
                call fftw_execute_r2r(plan1,b(:,j,k),b(:,j,k))
                b(:,j,k) = b(:,j,k)/float(ny)
                b(1,j,k) = b(1,j,k)/2.0
                b(nyp,j,k) = b(nyp,j,k)/2.0
            end do
        end do
        !$omp end parallel do 
        
        ! x-transform 
        !$omp parallel do shared(a,b,plan2) private(i,j) collapse(2)
        do j = 1,nz
            do i = 1,nyp
                call fftw_execute_dft_r2c(plan2,b(i,j,:),aplan(i,j,:))
                aplan(i,j,:) = aplan(i,j,:)/float(nx)
            end do
        end do
        !$omp end parallel do 
        
        ! z-transform
        !$omp parallel do shared(a,plan3) private(i,k) collapse(2)
        do k = 1,nxh
            do i = 1,nyp
                call fftw_execute_dft(plan3,aplan(i,:,k),aplan(i,:,k))
                aplan(i,:,k) = aplan(i,:,k)/float(nz)
            end do
        end do
        !$omp end parallel do 

        !$omp parallel do shared(a,aplan) private(i,j,k) schedule(dynamic)
        do k = 1,nxh
        do j = 1,nz
        do i = 1,nyp 
        a(i,j,k) = aplan(i,j,k)
        end do
        end do
        end do
        !$omp end parallel do
         
    end subroutine xyzfft
    
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
        complex :: lambda1,lambda2 ! Eigenvalues
    

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
   
        call eig2(u11,u12,u21,u22,l1,l2)

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
        use omp_lib
 
        implicit none
    
        ! Passed variables  
        real, dimension(nyp,mz,mx) :: u11,u12,u13,u21,u22,u23,u31,u32,u33,swirl
    
        ! Calculation variables
        complex :: l1, l2, l3
        real    :: l1_im, l2_im, l3_im
        integer :: i,j,k
   
        !$omp parallel do default(shared) private(i,j,k,l1_im,l2_im,l3_im,l1,l2,l3) collapse(3) schedule(dynamic)
        do k = 1,mx
            do j = 1,mz
                do i = 1,nyp 
                    call eig(u11(i,j,k),u12(i,j,k),u13(i,j,k),u21(i,j,k),u22(i,j,k),u23(i,j,k),u31(i,j,k),u32(i,j,k),u33(i,j,k),l1,l2,l3)
    
                    l1_im = aimag(l1)
                    l2_im = aimag(l2)
                    l3_im = aimag(l3)
    
                    swirl(i,j,k) = sqrt(l1_im**2 + l2_im**2 + l3_im**2) 
    
                    ! Certain velocity gradient tensors produce NaNs using the equations in
                    ! eig(...). When this happens, the actual eigenvalues are all 0. This is a
                    ! quick check to see if the swirl is a NaN. If it's any real value,
                    ! floor(swirl/swirl) will always be 1; otherwise it's a NaN (or Inf, but
                    ! that doesn't happen in this case). - Ryan 10/28/22
                    if (floor(swirl(i,j,k)/swirl(i,j,k)) .ne. 1) then
                        swirl(i,j,k) = 0.0
                    end if
                end do
            end do  
        end do
        !$omp end parallel do 
    
    end subroutine calcswirl
    
    !---------------------------------------------------------------------!
#IFDEF POLYMER    
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
        integer :: irstrt
        integer :: src_start,src_stop
        real    :: xl,yl,zl

        ! Polymer variables
        integer :: ipolyflag,itarget,ipeter
        real    :: alpha_poly,tpoly,zlmax,diffpoly,qbeta,beta_min

        !- Common blocks -!
        common/iocontrl/  irstrt
        common/domain/    xl,yl,zl
        common/poly_flgs/ ipolyflag,itarget,ipeter
        common/poly_var/  alpha_poly,tpoly,zlmax,diffpoly,qbeta,beta_min
        common/src_time/  src_start,src_stop

        !- Begin Calculations -!

        avg_beta = 0.0
        total_scl = 0.0
        volume = xl*yl*zl

        !$omp parallel do reduction(+:avg_beta,total_scl) default(shared) private(i,j,k,Lx,Ly,Lz)
        do i = 1,nyp
            ! Calc y length
            if (i .eq. 1) then
                Ly = abs(ycoord(2) - ycoord(1))/2.0
            else if (i .eq. nyp) then
                Ly = abs(ycoord(nyp) - ycoord(ny))/2.0
            else
                Ly = abs(ycoord(i+1) - ycoord(i))/2.0 + abs(ycoord(i) - ycoord(i-1))/2.0
            end if
            do j = 1,mz
                ! Calc z length
                Lz = delzm
                do k = 1,mx
                    ! Calc x length
                    Lx = delxm

                    total_scl = total_scl + scp(i,j,k)*Lx*Ly*Lz
                    avg_beta = avg_beta + beta(i,j,k)*Lx*Ly*Lz
                end do
            end do
        end do
        !$omp end parallel do

        avg_beta = avg_beta/volume

        if (mod(it,cadence) .eq. 1) then
            print *,'avg beta = ',avg_beta
            print *,'avg scalar = ',total_scl/volume
        end if

!        ! Write data to output files - updated each time step        
!        if (it .eq. irstrt) then
!            open(75,file='outputs/total_scl')
!        else
!            open(75,file='outputs/total_scl',position='append')
!        end if
!        
!        write(75,*) total_scl/volume*1.0e-6 ! PPM --> weight/weight raw #
!        close(75) 

        ! Check stop condition on polymer release
        if (avg_beta .le. beta_min) then
            src_stop = it
            print *,'Final polymer addition @ time step: ',it
        end if
        qbeta = avg_beta
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

    ! y-derivatives
    call cderiv(s12,wrk1)
    call cderiv(s22,wrk2)
    call cderiv(s23,wrk3)

    ! Compute div(s)    
    !$omp parallel do default(shared) private(i,j,k) collapse(3)
    do k = 1,nxh
        do j = 1,nz
            do i = 1,nyp
                t1(i,j,k) = im*(wavx(k)*s11(i,j,k) + wavz(j)*s13(i,j,k)) + wrk1(i,j,k) 
                t2(i,j,k) = im*(wavx(k)*s12(i,j,k) + wavz(j)*s23(i,j,k)) + wrk2(i,j,k) 
                t3(i,j,k) = im*(wavx(k)*s13(i,j,k) + wavz(j)*s33(i,j,k)) + wrk3(i,j,k) 
            end do
        end do
    end do
    !$omp end parallel do

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
   
    !$omp parallel do default(shared) private(i,j,k) schedule(dynamic) 
    do k = 1,nxh
        do j = 1,nz
            do i = 1,nyp
                gn (i,j,k) = gn (i,j,k) + t1(i,j,k)
                fn (i,j,k) = fn (i,j,k) + t2(i,j,k)
                omz(i,j,k) = omz(i,j,k) + t3(i,j,k)
            end do
        end do
    end do
    !$omp end parallel do

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
    !$omp parallel do default(shared) private(i,j,k) schedule(dynamic)
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
    !$omp end parallel do

    end subroutine polyNL

    !---------------------------------------------------------------------!

    subroutine findmaxQ(u11,u12,u13,u21,u22,u23,u31,u32,u33,scalar,Qmin,Qx,Qy,Qz,beta,planY,planZb,planXb)
    ! This subroutine is specifically for targeting regions of low Q/high strain.
    ! I suppose it could be easily modified to target other regions without particles,
    ! but I'll worry about that later.

    ! This takes the velocity gradient tensor and scalar field as inputs and outputs
    ! the grid indices of the minimum Q coordinates and the corresponding beta field
    ! for checking local concentration. These quantities must be in 3/2 interpolated 
    ! physical space for whatever reason that I don't really know, but that's why all
    ! these FFTs are necessary.

    ! Transform velocity gradient terms to 3/2 physical space 

    use,intrinsic :: iso_c_binding
    use grid_size
    use omp_lib
    use derivs

    !---------------------------------------------------------------------!
    
    ! =================================================================== !
    !                   Declare all variables explicitly                  !
    ! =================================================================== !

    implicit none

    include 'fftw3.f03'

    ! Inputs  
    complex, dimension(nyp,nz,nxh), intent(in) :: u11,u12,u13,u21,u22,u23,u31,u32,u33,scalar

    ! Outputs
    integer, dimension(qn), intent(in out) :: Qx(:),Qy(:),Qz(:)
    real,    dimension(qn), intent(in out) :: Qmin(:)
    real,    dimension(nyp,mz,mx)          :: beta

    ! FFTW Variables
    type(C_PTR) :: planY,planZb,planXb

    ! Real/Imaginary parts of spectral variables
    real(C_DOUBLE), dimension(nyp,nz,nxh) :: u11r,u12r,u13r,u21r,u22r,u23r,u31r,u32r,u33r,scr
    real(C_DOUBLE), dimension(nyp,nz,nxh) :: u11i,u12i,u13i,u21i,u22i,u23i,u31i,u32i,u33i,sci
    

    ! Intermediate physical y-plane variables
    complex(C_DOUBLE_COMPLEX), dimension(mz,mx) :: u11s,u12s,u13s
    complex(C_DOUBLE_COMPLEX), dimension(mz,mx) :: u21s,u22s,u23s
    complex(C_DOUBLE_COMPLEX), dimension(mz,mx) :: u31s,u32s,u33s
    complex(C_DOUBLE_COMPLEX), dimension(mz,mx) :: scs

    real(C_DOUBLE), dimension(mz,mx) :: u11p,u12p,u13p
    real(C_DOUBLE), dimension(mz,mx) :: u21p,u22p,u23p
    real(C_DOUBLE), dimension(mz,mx) :: u31p,u32p,u33p
    real(C_DOUBLE), dimension(mz,mx) :: scp,beta_poly

    ! Intermediate 3D physical variables
    real, dimension(nyp,mz,mx) :: Qcrit
    real, dimension(nyp,mz,mx) :: u11p3d,u12p3d,u13p3d
    real, dimension(nyp,mz,mx) :: u21p3d,u22p3d,u23p3d
    real, dimension(nyp,mz,mx) :: u31p3d,u32p3d,u33p3d

    ! Calculation Variables
    integer :: i,j,k,jj,cnt
    real    :: QQ

    real    :: wavx(nxh),wavz(nz) 
    real    :: c(nyp)
    real    :: fac
 
    ! Polymer variables
    real    :: alpha_poly
    !---------------------------------------------------------------------!

    ! =================================================================== !
    !                            Common Blocks                            !
    ! =================================================================== !
    common/waves/      wavx,wavz,c
    common/poly_var/   alpha_poly
    !---------------------------------------------------------------------!
    
    ! =================================================================== !
    !                         Begin Calculations                          !
    ! =================================================================== !

    !---------------------------------------------------------------------!
    !          Transform spectral variables to 3/2 physical space         !
    !---------------------------------------------------------------------!
    ! Convert Chebyshev modes to cosine modes for DCT
    !$omp parallel do default(shared) private(i,j,k,fac) schedule(dynamic)
    do k = 1,nxh
        do j = 1,nz
            do i = 1,nyp

                fac = c(i)/2.0
 
                u11r(i,j,k) = real(u11(i,j,k))*fac
                u12r(i,j,k) = real(u12(i,j,k))*fac
                u13r(i,j,k) = real(u13(i,j,k))*fac
                u21r(i,j,k) = real(u21(i,j,k))*fac
                u22r(i,j,k) = real(u22(i,j,k))*fac
                u23r(i,j,k) = real(u23(i,j,k))*fac
                u31r(i,j,k) = real(u31(i,j,k))*fac
                u32r(i,j,k) = real(u32(i,j,k))*fac
                u33r(i,j,k) = real(u33(i,j,k))*fac
                 scr(i,j,k) = real(scalar(i,j,k))*fac

                u11i(i,j,k) = aimag(u11(i,j,k))*fac
                u12i(i,j,k) = aimag(u12(i,j,k))*fac
                u13i(i,j,k) = aimag(u13(i,j,k))*fac
                u21i(i,j,k) = aimag(u21(i,j,k))*fac
                u22i(i,j,k) = aimag(u22(i,j,k))*fac
                u23i(i,j,k) = aimag(u23(i,j,k))*fac
                u31i(i,j,k) = aimag(u31(i,j,k))*fac
                u32i(i,j,k) = aimag(u32(i,j,k))*fac
                u33i(i,j,k) = aimag(u33(i,j,k))*fac
                 scr(i,j,k) = aimag(scalar(i,j,k))*fac
            end do
        end do
    end do
    !$omp end parallel do

    
    ! Compute Real --> Real DCT-I on real and imaginary parts of spectral variables
    !$omp parallel do default(shared) private(j,k) collapse(2)
    do k = 1,nxh
        do j = 1,nz

            call fftw_execute_r2r(planY,u11r(:,j,k),u11r(:,j,k))
            call fftw_execute_r2r(planY,u12r(:,j,k),u12r(:,j,k))
            call fftw_execute_r2r(planY,u13r(:,j,k),u13r(:,j,k))
            call fftw_execute_r2r(planY,u21r(:,j,k),u21r(:,j,k))
            call fftw_execute_r2r(planY,u22r(:,j,k),u22r(:,j,k))
            call fftw_execute_r2r(planY,u23r(:,j,k),u23r(:,j,k))
            call fftw_execute_r2r(planY,u31r(:,j,k),u31r(:,j,k))
            call fftw_execute_r2r(planY,u32r(:,j,k),u32r(:,j,k))
            call fftw_execute_r2r(planY,u33r(:,j,k),u33r(:,j,k))
            call fftw_execute_r2r(planY, scr(:,j,k), scr(:,j,k))

            call fftw_execute_r2r(planY,u11i(:,j,k),u11i(:,j,k))
            call fftw_execute_r2r(planY,u12i(:,j,k),u12i(:,j,k))
            call fftw_execute_r2r(planY,u13i(:,j,k),u13i(:,j,k))
            call fftw_execute_r2r(planY,u21i(:,j,k),u21i(:,j,k))
            call fftw_execute_r2r(planY,u22i(:,j,k),u22i(:,j,k))
            call fftw_execute_r2r(planY,u23i(:,j,k),u23i(:,j,k))
            call fftw_execute_r2r(planY,u31i(:,j,k),u31i(:,j,k))
            call fftw_execute_r2r(planY,u32i(:,j,k),u32i(:,j,k))
            call fftw_execute_r2r(planY,u33i(:,j,k),u33i(:,j,k))
            call fftw_execute_r2r(planY, sci(:,j,k), sci(:,j,k))
        end do
    end do
    !$omp end parallel do 

    ! Do the rest of the transforms in parallel over y-planes
    !$omp parallel do default(shared) private(           &
    !$omp       u11s,u12s,u13s,u21s,u22s,u23s,u31s,u32s,u33s, &
    !$omp       u11p,u12p,u13p,u21p,u22p,u23p,u31p,u32p,u33p, &
    !$omp       i,j,k,jj,beta_poly,scs,scp) schedule(dynamic)
    do i = 1,nyp    
    !---------------------------------------------------------------------!
    !       Calculate beta then store it and grad(V) in 3D variables      !
    !---------------------------------------------------------------------!
        ! Zero out 2D variables

        u11s = 0.0
        u12s = 0.0
        u13s = 0.0
        u21s = 0.0
        u22s = 0.0
        u23s = 0.0
        u31s = 0.0
        u32s = 0.0
        u33s = 0.0
         scs = 0.0

        ! Copy data to 2D variables
        do k = 1,nxh
            do j = 1,nz
                if (j .le. nzh) jj = j
                if (j .gt. nzh) jj = (mz-nz) + j

                u11s(jj,k) = cmplx(u11r(i,j,k),u11i(i,j,k))
                u12s(jj,k) = cmplx(u12r(i,j,k),u12i(i,j,k))
                u13s(jj,k) = cmplx(u13r(i,j,k),u13i(i,j,k))
                u21s(jj,k) = cmplx(u21r(i,j,k),u21i(i,j,k))
                u22s(jj,k) = cmplx(u22r(i,j,k),u22i(i,j,k))
                u23s(jj,k) = cmplx(u23r(i,j,k),u23i(i,j,k))
                u31s(jj,k) = cmplx(u31r(i,j,k),u31i(i,j,k))
                u32s(jj,k) = cmplx(u32r(i,j,k),u32i(i,j,k))
                u33s(jj,k) = cmplx(u33r(i,j,k),u33i(i,j,k))
                 scs(jj,k) = cmplx( scr(i,j,k), sci(i,j,k))

            end do
        end do
        
        ! Complex --> Complex z-transform
        call fftw_execute_dft(planZb,u11s,u11s)
        call fftw_execute_dft(planZb,u12s,u12s)
        call fftw_execute_dft(planZb,u13s,u13s)
        call fftw_execute_dft(planZb,u21s,u21s)
        call fftw_execute_dft(planZb,u22s,u22s)
        call fftw_execute_dft(planZb,u23s,u23s)
        call fftw_execute_dft(planZb,u31s,u31s)
        call fftw_execute_dft(planZb,u32s,u32s)
        call fftw_execute_dft(planZb,u33s,u33s)
        call fftw_execute_dft(planZb, scs, scs)

        ! Complex --> Real x-transform
        call fftw_execute_dft_c2r(planXb,u11s,u11p)
        call fftw_execute_dft_c2r(planXb,u12s,u12p)
        call fftw_execute_dft_c2r(planXb,u13s,u13p)
        call fftw_execute_dft_c2r(planXb,u21s,u21p)
        call fftw_execute_dft_c2r(planXb,u22s,u22p)
        call fftw_execute_dft_c2r(planXb,u23s,u23p)
        call fftw_execute_dft_c2r(planXb,u31s,u31p)
        call fftw_execute_dft_c2r(planXb,u32s,u32p)
        call fftw_execute_dft_c2r(planXb,u33s,u33p)
        call fftw_execute_dft_c2r(planXb, scs, scp)

 
        do k = 1,mx
            do j = 1,mz
                beta_poly(j,k) = exp(-alpha_poly*abs(scp(j,k)))

                  beta(i,j,k) = beta_poly(j,k)
                u11p3d(i,j,k) = u11p(j,k) 
                u12p3d(i,j,k) = u12p(j,k) 
                u13p3d(i,j,k) = u13p(j,k) 
                u21p3d(i,j,k) = u21p(j,k) 
                u22p3d(i,j,k) = u22p(j,k) 
                u23p3d(i,j,k) = u23p(j,k) 
                u31p3d(i,j,k) = u31p(j,k) 
                u32p3d(i,j,k) = u32p(j,k) 
                u33p3d(i,j,k) = u33p(j,k) 
            end do
        end do
    end do ! k
    !$omp end parallel do

    !---------------------------------------------------------------------!
    !                        Calculate Q-criterion                        !
    !---------------------------------------------------------------------!
    ! If we want to find the maximum values of Q (or strain), we have to flip the signs of the array,
    ! then find local minima and sort smallest to largest. The position vectors are rearranged acccordingly,
    ! so that the indices of the most negative Q-criterion (after being flipped) are first, and so on. The
    ! values inside Qmin are not important after determining the local minima, so they are not flipped back.

    ! Calculate Q-criterion and store in Qcrit
    call calcQ(u11p3d,u21p3d,u31p3d,u12p3d,u22p3d,u32p3d,u13p3d,u23p3d,u33p3d,Qcrit)

    ! Check for local (grid) min | DO NOT PARALLELIZE (unless you use an entirely different algorithm)
    Qcrit = -1.0*Qcrit ! Sort largest -> smallest
    cnt = 1
    do k = 2,mx-1
        do j = 2,mz-1
            do i = 2,ny
                if (Qcrit(i,j,k) .lt. 0 .and. cnt .le. size(Qmin)) then ! First check - Negative requirement and less than the max size of Qmin
                    if (Qcrit(i,j,k) .lt. Qcrit(i+1,j,k) .and. Qcrit(i,j,k) .lt. Qcrit(i-1,j,k)) then ! Check if minimum of y-neighbors
                        if (Qcrit(i,j,k) .lt. Qcrit(i,j+1,k) .and. Qcrit(i,j,k) .lt. Qcrit(i,j-1,k)) then ! Check if minimum of z-neighbors
                            if (Qcrit(i,j,k) .lt. Qcrit(i,j,k+1) .and. Qcrit(i,j,k) .lt. Qcrit(i,j,k-1)) then ! Check if minimum of x-neighbors
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

    ! Sort Qmin from smallest to largest (negative --> positive), and rearrange (x,y,z) vectors accordingly
    call Shell_Sort(Qmin,Qx,Qy,Qz)

    end subroutine findmaxQ

    !---------------------------------------------------------------------!

    subroutine newtarget(u11,u12,u13,u21,u22,u23,u31,u32,u33,scsource,planY,planZb,planXb)
    ! This subroutine is specifically for targeting regions of low Q/high strain.
    ! I suppose it could be easily modified to target other regions without particles,
    ! but I'll worry about that later.

    ! This takes the velocity gradient tensor and scalar field as inputs and outputs
    ! the grid indices of the minimum Q coordinates and the corresponding beta field
    ! for checking local concentration. These quantities must be in 3/2 interpolated 
    ! physical space for whatever reason that I don't really know, but that's why all
    ! these FFTs are necessary.

    ! Transform velocity gradient terms to 3/2 physical space 

    use,intrinsic :: iso_c_binding
    use grid_size
    use omp_lib
    use derivs

    !---------------------------------------------------------------------!
    
    ! =================================================================== !
    !                   Declare all variables explicitly                  !
    ! =================================================================== !

    implicit none

    include 'fftw3.f03'

    ! Inputs  
    complex, dimension(nyp,nz,nxh), intent(in) :: u11,u12,u13,u21,u22,u23,u31,u32,u33

    ! Outputs
    real,    dimension(nyp,mz,mx)          :: scsource

    ! FFTW Variables
    type(C_PTR) :: planY,planZb,planXb

    ! Real/Imaginary parts of spectral variables
    real(C_DOUBLE), dimension(nyp,nz,nxh) :: u11r,u12r,u13r,u21r,u22r,u23r,u31r,u32r,u33r
    real(C_DOUBLE), dimension(nyp,nz,nxh) :: u11i,u12i,u13i,u21i,u22i,u23i,u31i,u32i,u33i
    

    ! Intermediate physical y-plane variables
    complex(C_DOUBLE_COMPLEX), dimension(mz,mx) :: u11s,u12s,u13s
    complex(C_DOUBLE_COMPLEX), dimension(mz,mx) :: u21s,u22s,u23s
    complex(C_DOUBLE_COMPLEX), dimension(mz,mx) :: u31s,u32s,u33s

    real(C_DOUBLE), dimension(mz,mx) :: u11p,u12p,u13p
    real(C_DOUBLE), dimension(mz,mx) :: u21p,u22p,u23p
    real(C_DOUBLE), dimension(mz,mx) :: u31p,u32p,u33p

    ! Intermediate 3D physical variables
    real, dimension(nyp,mz,mx) :: u11p3d,u12p3d,u13p3d
    real, dimension(nyp,mz,mx) :: u21p3d,u22p3d,u23p3d
    real, dimension(nyp,mz,mx) :: u31p3d,u32p3d,u33p3d
    real, dimension(nyp,mz,mx) :: Qcrit

    ! Calculation Variables
    integer :: i,j,k,jj,cnt
    real    :: xl,yl,zl
    real    :: Lx,Ly,Lz,volume,scsum,var1,Qmax

    real    :: wavx(nxh),wavz(nz) 
    real    :: c(nyp)
    real    :: fac
 
    ! Polymer variables
    real    :: alpha_poly
    !---------------------------------------------------------------------!

    ! =================================================================== !
    !                            Common Blocks                            !
    ! =================================================================== !
    common/waves/      wavx,wavz,c
    common/poly_var/   alpha_poly
    common/domain/     xl,yl,zl
    !---------------------------------------------------------------------!
    
    ! =================================================================== !
    !                         Begin Calculations                          !
    ! =================================================================== !

    !---------------------------------------------------------------------!
    !          Transform spectral variables to 3/2 physical space         !
    !---------------------------------------------------------------------!
    ! Convert Chebyshev modes to cosine modes for DCT
    !$omp parallel do default(shared) private(i,j,k,fac) schedule(dynamic)
    do k = 1,nxh
        do j = 1,nz
            do i = 1,nyp

                fac = c(i)/2.0
 
                u11r(i,j,k) = real(u11(i,j,k))*fac
                u12r(i,j,k) = real(u12(i,j,k))*fac
                u13r(i,j,k) = real(u13(i,j,k))*fac
                u21r(i,j,k) = real(u21(i,j,k))*fac
                u22r(i,j,k) = real(u22(i,j,k))*fac
                u23r(i,j,k) = real(u23(i,j,k))*fac
                u31r(i,j,k) = real(u31(i,j,k))*fac
                u32r(i,j,k) = real(u32(i,j,k))*fac
                u33r(i,j,k) = real(u33(i,j,k))*fac

                u11i(i,j,k) = aimag(u11(i,j,k))*fac
                u12i(i,j,k) = aimag(u12(i,j,k))*fac
                u13i(i,j,k) = aimag(u13(i,j,k))*fac
                u21i(i,j,k) = aimag(u21(i,j,k))*fac
                u22i(i,j,k) = aimag(u22(i,j,k))*fac
                u23i(i,j,k) = aimag(u23(i,j,k))*fac
                u31i(i,j,k) = aimag(u31(i,j,k))*fac
                u32i(i,j,k) = aimag(u32(i,j,k))*fac
                u33i(i,j,k) = aimag(u33(i,j,k))*fac
            end do
        end do
    end do
    !$omp end parallel do

    
    ! Compute Real --> Real DCT-I on real and imaginary parts of spectral variables
    !$omp parallel do default(shared) private(j,k) collapse(2)
    do k = 1,nxh
        do j = 1,nz

            call fftw_execute_r2r(planY,u11r(:,j,k),u11r(:,j,k))
            call fftw_execute_r2r(planY,u12r(:,j,k),u12r(:,j,k))
            call fftw_execute_r2r(planY,u13r(:,j,k),u13r(:,j,k))
            call fftw_execute_r2r(planY,u21r(:,j,k),u21r(:,j,k))
            call fftw_execute_r2r(planY,u22r(:,j,k),u22r(:,j,k))
            call fftw_execute_r2r(planY,u23r(:,j,k),u23r(:,j,k))
            call fftw_execute_r2r(planY,u31r(:,j,k),u31r(:,j,k))
            call fftw_execute_r2r(planY,u32r(:,j,k),u32r(:,j,k))
            call fftw_execute_r2r(planY,u33r(:,j,k),u33r(:,j,k))

            call fftw_execute_r2r(planY,u11i(:,j,k),u11i(:,j,k))
            call fftw_execute_r2r(planY,u12i(:,j,k),u12i(:,j,k))
            call fftw_execute_r2r(planY,u13i(:,j,k),u13i(:,j,k))
            call fftw_execute_r2r(planY,u21i(:,j,k),u21i(:,j,k))
            call fftw_execute_r2r(planY,u22i(:,j,k),u22i(:,j,k))
            call fftw_execute_r2r(planY,u23i(:,j,k),u23i(:,j,k))
            call fftw_execute_r2r(planY,u31i(:,j,k),u31i(:,j,k))
            call fftw_execute_r2r(planY,u32i(:,j,k),u32i(:,j,k))
            call fftw_execute_r2r(planY,u33i(:,j,k),u33i(:,j,k))
        end do
    end do
    !$omp end parallel do 

    ! Do the rest of the transforms in parallel over y-planes
    !$omp parallel do default(shared) private(           &
    !$omp       u11s,u12s,u13s,u21s,u22s,u23s,u31s,u32s,u33s, &
    !$omp       u11p,u12p,u13p,u21p,u22p,u23p,u31p,u32p,u33p, &
    !$omp       i,j,k,jj) schedule(dynamic)
    do i = 1,nyp    
    !---------------------------------------------------------------------!
    !       Calculate beta then store it and grad(V) in 3D variables      !
    !---------------------------------------------------------------------!
        ! Zero out 2D variables

        u11s = 0.0
        u12s = 0.0
        u13s = 0.0
        u21s = 0.0
        u22s = 0.0
        u23s = 0.0
        u31s = 0.0
        u32s = 0.0
        u33s = 0.0

        ! Copy data to 2D variables
        do k = 1,nxh
            do j = 1,nz
                if (j .le. nzh) jj = j
                if (j .gt. nzh) jj = (mz-nz) + j

                u11s(jj,k) = cmplx(u11r(i,j,k),u11i(i,j,k))
                u12s(jj,k) = cmplx(u12r(i,j,k),u12i(i,j,k))
                u13s(jj,k) = cmplx(u13r(i,j,k),u13i(i,j,k))
                u21s(jj,k) = cmplx(u21r(i,j,k),u21i(i,j,k))
                u22s(jj,k) = cmplx(u22r(i,j,k),u22i(i,j,k))
                u23s(jj,k) = cmplx(u23r(i,j,k),u23i(i,j,k))
                u31s(jj,k) = cmplx(u31r(i,j,k),u31i(i,j,k))
                u32s(jj,k) = cmplx(u32r(i,j,k),u32i(i,j,k))
                u33s(jj,k) = cmplx(u33r(i,j,k),u33i(i,j,k))

            end do
        end do
        
        ! Complex --> Complex z-transform
        call fftw_execute_dft(planZb,u11s,u11s)
        call fftw_execute_dft(planZb,u12s,u12s)
        call fftw_execute_dft(planZb,u13s,u13s)
        call fftw_execute_dft(planZb,u21s,u21s)
        call fftw_execute_dft(planZb,u22s,u22s)
        call fftw_execute_dft(planZb,u23s,u23s)
        call fftw_execute_dft(planZb,u31s,u31s)
        call fftw_execute_dft(planZb,u32s,u32s)
        call fftw_execute_dft(planZb,u33s,u33s)

        ! Complex --> Real x-transform
        call fftw_execute_dft_c2r(planXb,u11s,u11p)
        call fftw_execute_dft_c2r(planXb,u12s,u12p)
        call fftw_execute_dft_c2r(planXb,u13s,u13p)
        call fftw_execute_dft_c2r(planXb,u21s,u21p)
        call fftw_execute_dft_c2r(planXb,u22s,u22p)
        call fftw_execute_dft_c2r(planXb,u23s,u23p)
        call fftw_execute_dft_c2r(planXb,u31s,u31p)
        call fftw_execute_dft_c2r(planXb,u32s,u32p)
        call fftw_execute_dft_c2r(planXb,u33s,u33p)

 
        do k = 1,mx
            do j = 1,mz
                u11p3d(i,j,k) = u11p(j,k) 
                u12p3d(i,j,k) = u12p(j,k) 
                u13p3d(i,j,k) = u13p(j,k) 
                u21p3d(i,j,k) = u21p(j,k) 
                u22p3d(i,j,k) = u22p(j,k) 
                u23p3d(i,j,k) = u23p(j,k) 
                u31p3d(i,j,k) = u31p(j,k) 
                u32p3d(i,j,k) = u32p(j,k) 
                u33p3d(i,j,k) = u33p(j,k) 
            end do
        end do
    end do ! k
    !$omp end parallel do

    !---------------------------------------------------------------------!
    !                        Calculate Q-criterion                        !
    !---------------------------------------------------------------------!

    ! Calculate Q-criterion and store in Qcrit
    call calcswirl(u11p3d,u21p3d,u31p3d,u12p3d,u22p3d,u32p3d,u13p3d,u23p3d,u33p3d,Qcrit)

    Qmax = maxval(Qcrit)
    Lx = xl/float(mx)
    Lz = zl/float(mz)
    volume = xl*yl*zl
    !$omp parallel do default(shared) private(i,j,k,Ly) reduction(+:scsum)
    do k = 1,mx
        do j = 1,mz
            do i = 1,nyp
                if (i .eq. 1) then
                    Ly = abs(ycoord(2) - ycoord(1))/2.0
                else if (i .eq. nyp) then
                    Ly = abs(ycoord(nyp) - ycoord(ny))/2.0
                else
                    Ly = abs(ycoord(i+1) - ycoord(i))/2.0 + abs(ycoord(i) - ycoord(i-1))/2.0
                end if

                scsource(i,j,k) = Qcrit(i,j,k)/Qmax ! scaled from 0 to 1
                scsum = scsum + scsource(i,j,k)*Lx*Ly*Lz
            end do
        end do  
    end do
    !$omp end parallel do

    ! Scale scalar by total amount to add per time step
    var1 = spts/scsum ! Compare total scalar added to target value (spts)
    do k = 1,mx
        do j = 1,mz
            do i = 1,nyp
                scsource(i,j,k) = scsource(i,j,k)*var1
            end do
        end do
    end do

    end subroutine newtarget
#ENDIF   
 
    !---------------------------------------------------------------------!

    subroutine writeoutputs(u,v,w,wx,wy,wz,u12,umean,swirl &
#IFDEF SCALAR
                            ,beta &
#IFDEF POLYMER
                            ,p12,trC &
#ENDIF
#ENDIF
                           )
    ! Compute output data
    ! If scalar, beta = scp3d

    use grid_size
    use omp_lib

    !---------------------------------------------------------------------!
    
    ! =================================================================== !
    !                   Declare all variables explicitly                  !
    ! =================================================================== !

        implicit none

        ! Input variables
        real, dimension(nyp,mz,mx) :: u,v,w,wx,wy,wz,beta,u12,swirl
        real, dimension(nyp,mz,mx) :: p12,trC

        ! Shear stress (polymer)
        real :: s_tw,s_bw,tpss,bpss,Ckk
        real, dimension(nyp,mz,mx) :: up ! fluctuating streamwise velocity
        real, dimension(nyp)       :: umean,betay,swirly,urms,vrms,wrms,R12,T12,py12

        ! Common block variables
        integer :: irstrt
        integer :: it 
        real    :: xl,yl,zl
        real    :: re

        ! Calculation variables
        ! Shear stress (fluid)
        real :: tau_tw,tau_bw,twss,bwss

        ! Indexing/domain
        integer :: i,j,k
        real    :: Lx,Ly,Lz

        ! Misc
        real    :: volume,massFlux

    !---------------------------------------------------------------------!
   
    ! Common blocks
        common/iocontrl/ irstrt
        common/itime/    it
        common/domain/   xl,yl,zl
        common/flow/     re
     
    ! =================================================================== !
    !                         Begin Calculations                          !
    ! =================================================================== !
 
        Lx = xl/float(mx)
        Lz = zl/float(mz)
        volume = (xl*yl*zl)

        betay = 0.0
        swirly = 0.0
        tau_tw = 0.0
        tau_bw = 0.0
        urms = 0.0
        vrms = 0.0
        wrms = 0.0
        R12 = 0.0
        T12 = 0.0
        py12 = 0.0
        s_tw = 0.0
        s_bw = 0.0
        Ckk = 0.0
        !$omp parallel do reduction(+:tau_tw,tau_bw,s_tw,s_bw,swirly,betay,urms,vrms,wrms,R12,T12,py12,Ckk) default(shared) private(i,j,k,twss,bwss,tpss,bpss,Ly) schedule(dynamic)
        do k = 1,mx
            do j = 1,mz
                do i = 1,nyp
                    ! Calc y length
                    if (i .eq. 1) then
                        Ly = abs(ycoord(2) - ycoord(1))/2.0
                    else if (i .eq. nyp) then
                        Ly = abs(ycoord(nyp) - ycoord(ny))/2.0
                    else
                        Ly = abs(ycoord(i+1) - ycoord(i))/2.0 + abs(ycoord(i) - ycoord(i-1))/2.0
                    end if

                    ! Calculate u'
                    up(i,j,k) = u(i,j,k) - umean(i)

                    ! Calcualte average profiles vs y
                     betay(i) =  betay(i) + beta(i,j,k)*Lx*Lz ! beta vs y
                    swirly(i) = swirly(i) + swirl(i,j,k)*Lx*Lz ! swirl vs y
                      urms(i) =   urms(i) + up(i,j,k)**2.0 ! u'_rms
                      vrms(i) =   vrms(i) +  v(i,j,k)**2.0 ! v'_rms
                      wrms(i) =   wrms(i) +  w(i,j,k)**2.0 ! w'_rms
                       R12(i) =    R12(i) + up(i,j,k)*v(i,j,k) ! RSS
                       T12(i) =    T12(i) + u12(i,j,k)/re ! Mean shear stress
                      py12(i) =   py12(i) + p12(i,j,k) ! polymer shear stress
                       Ckk = Ckk + trC(i,j,k)*Lx*Ly*Lz
                end do

                ! Calc stress on top and bottom walls only
                twss = u12(1,j,k)
                bwss = u12(nyp,j,k)

                tau_tw = tau_tw + twss/re
                tau_bw = tau_bw + bwss/re
                tpss = p12(1,j,k) 
                bpss = p12(nyp,j,k) 

                s_tw = s_tw + tpss
                s_bw = s_bw + bpss
            end do
        end do
        !$omp end parallel do

        massFlux = 0.0
        do i = 1,nyp
             R12(i) = R12(i)/(mx*mz)
             T12(i) = T12(i)/(mx*mz)
            py12(i) = py12(i)/(mx*mz)
            urms(i) = sqrt(urms(i)/(mx*mz))
            vrms(i) = sqrt(vrms(i)/(mx*mz))
            wrms(i) = sqrt(wrms(i)/(mx*mz))

            if (i .lt. nyp) then
            massFlux = massFlux + 1.0/yl*0.5*(umean(i+1) + umean(i))*abs(ycoord(i+1) - ycoord(i)) ! Bulk velocity
            end if
        end do

        ! Finish computing averages and correlation
        betay = betay/(xl*zl)
        swirly = swirly/(xl*zl)

        ! Write data to output files - updated each time step        
        if (it .eq. irstrt) then
            open(10,     file='outputs/Txy0')
            open(11,    file='outputs/urms0')
            open(12,     file='outputs/RSS0')
            open(22,    file='outputs/vrms0')
            open(33,    file='outputs/wrms0')
            open(99,file='outputs/massflux0')
            open(100, file='outputs/mean_u0')
            open(101,  file='outputs/betay0')
            open(102, file='outputs/swirly0')
            open(103,   file='outputs/twss0')
            open(104,   file='outputs/bwss0')
            open(105,   file='outputs/tpss0')
            open(106,   file='outputs/bpss0')
            open(107,    file='outputs/Ckk0')
            open(108,    file='outputs/PSS0')
        else
            open(105,   file='outputs/tpss0',position='append')
            open(106,   file='outputs/bpss0',position='append')
            open(107,    file='outputs/Ckk0',position='append')
            open(108,    file='outputs/PSS0',position='append')
            open(10,     file='outputs/Txy0',position='append')
            open(11,    file='outputs/urms0',position='append')
            open(12,     file='outputs/RSS0',position='append')
            open(22,    file='outputs/vrms0',position='append')
            open(33,    file='outputs/wrms0',position='append')
            open(99,file='outputs/massflux0',position='append')
            open(100, file='outputs/mean_u0',position='append')
            open(101,  file='outputs/betay0',position='append')
            open(102, file='outputs/swirly0',position='append')
            open(103,   file='outputs/twss0',position='append')
            open(104,   file='outputs/bwss0',position='append')
        end if

        ! Write y-profile data in order
        do i = 1,nyp
            write(10,*)  T12(i)
            write(11,*) urms(i)
            write(12,*)  R12(i)
            write(22,*) vrms(i)
            write(33,*) wrms(i)
            write(100,*) umean(i)
            write(101,*) betay(i)
            write(102,*) swirly(i)
            write(108,*) py12(i)
        end do
        write(99,*) massFlux
        write(103,*) tau_tw/(mx*mz)
        write(104,*) tau_bw/(mx*mz)

        close(10)
        close(11)
        close(12)
        close(22)
        close(33)
        close(99)
        close(100)
        close(101)
        close(102)
        close(103)
        close(104)

        write(105,*) s_tw/(mx*mz)
        write(106,*) s_bw/(mx*mz)
        write(107,*) Ckk/volume
        close(105)
        close(106)
        close(107)
        close(108)

    end subroutine writeoutputs

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
   
    delxm = xl/(mx)
    delzm = zl/(mz)

    area = 0.0

    !$omp parallel do reduction(+:area) default(shared) private(i,j,k,Lx,Ly,Lz) 
    do i = 1,nyp
        ! Calc y length
        if (i .eq. 1) then
            Ly = abs(ycoord(2) - ycoord(1))/2.0
        else if (i .eq. nyp) then
            Ly = abs(ycoord(nyp) - ycoord(ny))/2.0
        else
            Ly = abs(ycoord(i+1) - ycoord(i))/2.0 + abs(ycoord(i) - ycoord(i-1))/2.0
        end if
        do j = 1,mz
            ! Calc z length
            Lz = delzm
            
            do k = 1,mx
                ! Calc x length
                Lx = delxm

!                area = area + Ly*Lz*Lx*max((swirl(i,j,k) - 10.0),0.0)/(swirl(i,j,k)-10.0) ! This is a goofy way to avoid an if statement, but it was fun to come up with
                if (swirl(i,j,k) .ge. 10.0) then
                    area = area + Ly*Lz*Lx
                end if
            end do
        end do
    end do
    !$omp end parallel do

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
        use omp_lib
    
        implicit none
    
        ! Passed variables  
        real,dimension(nyp,mz,mx) :: u11,u12,u13,u21,u22,u23,u31,u32,u33,Q
    
        ! Calculation variables
        integer :: i,j,k
        real :: Rot, Strain
 
        ! Begin calculations
        !$omp parallel do default(shared) private(i,j,k,Strain,Rot) schedule(auto)
        do k = 1,mx
            do j = 1,mz
                do i = 1,nyp
                    Strain = u11(i,j,k)**2 + u22(i,j,k)**2 + u33(i,j,k)**2
                    Rot = -2.0*(u12(i,j,k)*u21(i,j,k) + u13(i,j,k)*u31(i,j,k) + u23(i,j,k)*u32(i,j,k))
!                    Strain = u11(i,j,k)**2 + u22(i,j,k)**2 + u33(i,j,k)**2 + 0.5*((u12(i,j,k) + u21(i,j,k))**2 + (u13(i,j,k) + u31(i,j,k))**2 + (u23(i,j,k) + u32(i,j,k))**2)
                    Q(i,j,k) = 0.5*(Rot - Strain)
                end do
            end do
        end do
        !$omp end parallel do
         
    end subroutine calcQ
    
    !---------------------------------------------------------------------!

    ! Sorting algorithm taken from Wikipedia
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

    subroutine correlate_vars(x,y,it)
    ! Calculates the correlation between 2 domain variables - must be in 
    ! physical space and the same size

    use grid_size
    use omp_lib
    !---------------------------------------------------------------------!
    
    ! =================================================================== !
    !                   Declare all variables explicitly                  !
    ! =================================================================== !

        implicit none

        ! Input variables
        real, dimension(nyp,mz,mx) :: x,y
        integer :: it ! Just for file writing ease

        ! Calculation variables
        integer :: i,j,k,N
        real    :: xbar,ybar,xsum,ysum
        real    :: Rnum,Rdenom,Rxy

        ! Simulation control variables
        integer :: irstrt
        common/iocontrl/  irstrt
    !---------------------------------------------------------------------!
    
    ! =================================================================== !
    !                         Begin Calculations                          !
    ! =================================================================== !
   
    N = nyp*mz*mx ! # of elements

    ! Compute mean of x and y
    xbar = sum(x)/N 
    ybar = sum(y)/N

    ! Calculate numerator of correlation coefficient & terms for denominator
    Rnum = 0.0
    xsum = 0.0
    ysum = 0.0
    do k = 1,mx
        do j = 1,mz
            do i = 1,nyp
                Rnum = Rnum + (x(i,j,k) - xbar)*(y(i,j,k) - ybar)
                xsum = xsum + (x(i,j,k) - xbar)**2.0
                ysum = ysum + (y(i,j,k) - ybar)**2.0
            end do
        end do
    end do

    Rdenom = sqrt(xsum*ysum)

    if (Rdenom .eq. 0.0) then
        Rdenom = 1.0
    end if

    Rxy = Rnum/Rdenom

    if (it .eq. irstrt) then
        open(144,file='outputs/correlation',status='new')
    else
        open(144,file='outputs/correlation',status='old',position='append')
    end if

    write(144,*) Rxy
    close(144)

    end subroutine correlate_vars    

    ! ---------------------------------------------------------------------------- !
    
    subroutine write_flowfield_ascii(u, v, w, wx, wy, wz, swirl &
#IFDEF SCALAR
                                     ,scl &
#ENDIF
                                     )
    
        use grid_size
    
        ! Passed variables
        integer :: irstrt, nsteps,iprnfrq,print3d
        real, dimension(nyp,mz,mx) :: u, v, w, wx, wy, wz, swirl
#IFDEF SCALAR
        real, dimension(nyp,mz,mx) :: scl
#ENDIF
        character (33) :: filename
    
        integer :: it
        real :: dt
        real :: x, y, z, delxm, delzm
        integer :: i, j, k, mz_copy
    
        real xl, yl, zl
    
        common/iocontrl/ irstrt,nsteps,iprnfrq,print3d
        common/domain/   xl,yl,zl
        common/itime/    it
        common/dtime/    dt
    
        ! Printing in 2D workaround
        if (nz .gt. 10) then
           mz_copy = mz
        else
            mz_copy = 1
        end if
    
        delxm = xl /float(mx)
        delzm = zl /float(mz)
    
        print *, 'Writing data at time step ', it
        
        write(filename,'("outputs/flowfield/time-",i6.6,".dat")')it
        
        ! Write 3d-data 
        open(6, file = filename)
#IFDEF SCALAR
        write(6,9711) 'filetype = solution, variables = "swirl", "u", "v", "w", "wx", "wy", "wz", "scl"'
#ELSE
        write(6,9711) 'filetype = solution, variables = "swirl", "u", "v", "w", "wx", "wy", "wz"'
#ENDIF
        write(6,9712) 'zone f=point t="Field", solutiontime=', it/iprnfrq,',i=',mx, 'j=',1, 'k=', nyp, new_line('a')
        
        do k = 1,nyp
            y = ycoord(k)
            do j = 1,mz_copy
                z = float(j-1)*delzm 
                do i = 1,mx
                    x = float(i-1)*delxm
#IFDEF SCALAR
                    write(6,8611) swirl(k,j,i),u(k,j,i),v(k,j,i),w(k,j,i),wx(k,j,i),wy(k,j,i),wz(k,j,i),scl(k,j,i)
#ELSE
                    write(6,8611) swirl(k,j,i),u(k,j,i),v(k,j,i),w(k,j,i),wx(k,j,i),wy(k,j,i),wz(k,j,i)
#ENDIF
                 end do
            end do
        end do 
        
        write(*,*) 'Finished writing flowfield'
        close(6)
    
        8611 format(8(e14.6,1x))
        9711 format(a93)
        9712 format(a40,i3,a5,i4,a5,i4,a5,i4,a2)
       
    
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
    
        end if
    end subroutine
    
    !---------------------------------------------------------------------!
    
#ifdef OUTPUTFORM
    
    subroutine write_flowfield_plt(u, v, w, wx, wy, wz, swirl &
#IF DEFINED SCALAR
                                   ,scl &
#ENDIF
                                   )
        use iso_c_binding 
        use grid_size
    
        implicit none
    
        include 'tecio.f90'
    
        integer(c_int64_t) :: mz_copy,mx_copy
        integer :: irstrt, nsteps,iprnfrq,print3d,crstrt
        real, dimension(nyp,mz,mx) :: u, v, w, wx, wy, wz, swirl
#IF DEFINED SCALAR
        real, dimension(nyp,mz,mx) :: scl
#ENDIF
        character*1 :: NULLCHR
        character (35) :: filename
        
    
        real :: x, y, z, delxm, delzm
        integer :: i, j, k
    
        integer :: it
        real    :: dt
        real    :: xl, yl, zl
   
#IFDEF POLYMER
        integer :: ipolyflag
#ENDIF
 
        ! Common blocks
        common/iocontrl/   irstrt,nsteps,iprnfrq,print3d,crstrt
        common/itime/      it
        common/dtime/      dt
        common/domain/     xl,yl,zl
#IFDEF POLYMER
        common/poly_flgs/  ipolyflag
#ENDIF
    
        integer(c_int64_t) :: numValues
        integer(c_int32_t) :: fileFormat = 1 ! 0 = .plt | 1 = .szplt
        integer(c_int32_t) :: gridFileType = 1 ! GRID
        integer(c_int32_t) :: defVarType = 0
        integer(c_int32_t), allocatable :: varTypes(:)
        integer(c_int32_t), allocatable :: shareVarFromZone(:)
        integer(c_int32_t), allocatable :: valueLocation(:)
        integer(c_int32_t), allocatable :: passiveVarList(:)
        integer(c_int32_t) :: shareFaceNeighborsFromZone = 0
        integer(c_int64_t) :: numFaceConnections = 0
        integer(c_int32_t) :: faceNeighborMode = 0
        integer(c_int32_t) :: zone
        real(c_double) :: solutionTime
        real(c_float), allocatable :: floatValues(:)
        real(c_float), allocatable :: coordinates(:,:,:,:)
        type(c_ptr) :: gridFH = C_NULL_PTR
        type(c_ptr) :: FH = C_NULL_PTR
    
        NULLCHR = CHAR(0)
    
        !! Start outputting data
    
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
        numValues = (nyp) * (mz_copy) * (mx_copy )
        allocate(floatValues(numValues))
        allocate(varTypes(numValues))
        allocate(shareVarFromZone(numValues))
        allocate(valueLocation(numValues))
        allocate(passiveVarList(numValues))
    
        varTypes = 1
        shareVarFromZone = 0
        valueLocation = 1
        passiveVarList = 0
    
        i = tecFileWriterOpen("outputs/flowfield/grid.szplt"//NULLCHR, &
                              "Grid", &
                              "x,y,z", &
                              fileFormat, &
                              gridFileType, &
                              defVarType, &
                              C_NULL_PTR, &
                              gridFH)
    
        ! Create zone
        i = tecZoneCreateIJK( gridFH, &
            "Flowfield", &
            nyp, &
            mz_copy, &
            mx_copy, &
            varTypes, &
            shareVarFromZone, &
            valueLocation, &
            passiveVarList, &
            shareFaceNeighborsFromZone, &
            numFaceConnections, &
            faceNeighborMode, &
            zone)
    
        ! Compute Coordinates
        allocate(coordinates(nyp,mz_copy,mx_copy , 3))
        delxm = xl /float(mx)
        delzm = zl /float(mz)
        do k = 1,mx_copy
            x = float(k-1) * delxm
            do j = 1,mz_copy
                z = float(j-1) * delzm 
                do i = 1,nyp
                    y = ycoord(i)
                    coordinates(i,j,k,1) = x
                    coordinates(i,j,k,2) = y
                    coordinates(i,j,k,3) = z
                end do
            end do
        end do 
    
        ! Write x,y,z
        do j = 1,3
            floatValues = pack(coordinates(:,:,:,j), .true.)
            i = tecZoneVarWriteFloatValues(gridFH,zone,j,1,numValues,floatValues)
        end do
    
        deallocate(coordinates)
    
        ! Now write the flowfield
            
        write(filename,'("outputs/flowfield/time-",i6.6,".szplt")')it
      
        ! For whatever reason, opening the grid file as read-only and using it to write the solution file doesn't work
        ! so I have to write it at each print time step and keep it open while I write the solution file. Luckily it 
        ! doesn't take much time so while it's doing more work, it's negligible in wall-clock time.
        ! The grid file is simply overwritten, so this doesn't interfere with the data writing and the solution files
        ! still work perfectly in Tecplot
     
    !    i = tecFileReaderOpen("outputs/flowfield/grid.szplt"//NULLCHR, & 
    !        gridFH)
   
            i = tecFileWriterOpen(filename//NULLCHR, &
                "Solution", &
#IFDEF POLYMER
                "u,v,w,swirl,beta", &
#ELIF DEFINED SCALAR
                "u,v,w,swirl,scalar", &
#ELSE
                "u,v,w,swirl"//NULLCHR, &
#ENDIF
                1, &
                2, &
                defVarType, &
                gridFH, &
                FH)
        
        
        varTypes = 1
        shareVarFromZone = 0
        valueLocation = 1
        passiveVarList = 0
        
        ! Create zone
        i = tecZoneCreateIJK(FH, &
            "Flowfield"//NULLCHR, &
            nyp, &
            mz_copy, &
            mx_copy, &
            varTypes, &
            shareVarFromZone, &
            valueLocation, &
            passiveVarList, &
            shareFaceNeighborsFromZone, &
            numFaceConnections, &
            faceNeighborMode, &
            zone)
       
        solutionTime = float(it)
        i = tecZoneSetUnsteadyOptions(FH, &
            zone, &
            solutionTime, &
            1) ! Strand ID
 
        ! Write u
        floatValues = pack(u(:,1:mz_copy,1:mx_copy), .true.)
        i = tecZoneVarWriteFloatValues(FH,zone,1,1,numValues,floatValues)
        
        ! Write v
        floatValues = pack(v(:,1:mz_copy,1:mx_copy), .true.)
        i = tecZoneVarWriteFloatValues(FH,zone,2,1,numValues,floatValues)
        
        ! Write w
        floatValues = pack(w(:,1:mz_copy,1:mx_copy), .true.)
        i = tecZoneVarWriteFloatValues(FH,zone,3,1,numValues,floatValues)
        
        ! Write swirl
        floatValues = pack(swirl(:,1:mz_copy,1:mx_copy), .true.)
        i = tecZoneVarWriteFloatValues(FH,zone,4,1,numValues,floatValues)
   
        ! We can write vorticity if it's important, but Tecplot can calculate it just fine 
        ! Write wx
!        floatValues = pack(wx(:,1:mz_copy,1:mx_copy), .true.)
!        i = tecZoneVarWriteFloatValues(FH,zone,5,1,numValues,floatValues)
!        
!        ! Write wy
!        floatValues = pack(wy(:,1:mz_copy,1:mx_copy), .true.)
!        i = tecZoneVarWriteFloatValues(FH,zone,6,1,numValues,floatValues)
        
    !    ! Write wz
    !    floatValues = pack(wz(:,1:mz_copy,1:mx_copy), .true.)
    !    i = tecZoneVarWriteFloatValues(FH,zone,7,1,numValues,floatValues)
#IFDEF SCALAR       
        floatValues = pack(scl(:,1:mz_copy,1:mx_copy),.true.)
        i = tecZoneVarWriteFloatValues(FH,zone,5,1,numValues,floatValues)
#ENDIF

        ! Close file
        i = tecFileWriterClose(FH)
        i = tecFileWriterClose(gridFH)
    
        deallocate(floatValues,varTypes,shareVarFromZone,valueLocation,passiveVarList)
    end subroutine
#endif

    subroutine write_FTLE_output(u,v,w)

        use grid_size   

        implicit none

        real,dimension(nyp,mz,mx) :: u,v,w

        integer :: i,j,k
        integer :: it
        real    :: dt

        integer :: irstrt,nsteps,iprnfrq,print3d
        real :: xl,yl,zl
    
        real :: delxm,delzm

        character(17) :: filename

        real,dimension(nyp) :: yvec
        real,dimension(mz)  :: zvec
        real,dimension(mx)  :: xvec


        ! Common blocks
        common/iocontrl/ irstrt,nsteps,iprnfrq,print3d
        common/itime/    it
        common/dtime/    dt
        common/domain/   xl,yl,zl


        ! Calculate domain variables
        delxm = xl/float(mx)
        delzm = zl/float(mz)

        do i = 1,nyp
            yvec(i) = ycoord(i)
        end do
        do j = 1,mz
            zvec(j) = float(j-1)*delzm
        end do
        do k = 1,mx
            xvec(k) = float(k-1)*delxm
        end do


        ! Write grid file at first write
        if (it .eq. 1) then
            ! Write grid file in ASCII format for readability
            open(122,file="outputs/flowfield/FTLE_grid.dat")
            write(122,*) "NYP = ",nyp,", MZ = ",mz,", MX = ",mx
            write(122,*) "LY = ",yl,", LZ = ",zl,", LX = ",xl
            write(122,*) "DT = ",dt,", IPRNFRQ = ",iprnfrq,", NSTEPS = ",nsteps
            write(122,*) xvec
            write(122,*) yvec
            write(122,*) zvec
            close(122)

            ! Write grid file in unformatted binary
            open(123,file="outputs/flowfield/FTLE_grid",status="new",form="unformatted")
            write(123) nyp,mz,mx
            write(123) yl,zl,xl
            write(123) dt,iprnfrq,nsteps
            write(123) xvec
            write(123) yvec
            write(123) zvec
            close(123)
        end if


        ! Write flowfield file every write
        write(filename,'("FTLE_veloc-",i6.6)')it
        open(124,file="outputs/flowfield/"//filename,form="unformatted")
        write(124) u,v,w
        close(124)
    end subroutine write_FTLE_output

    !---------------------------------------------------------------------!
    
    subroutine calc_stress(umean,u,v,w,u12 &
#IFDEF POLYMER
                        ,p12 &
#ENDIF
                        )

        use grid_size
        use omp_lib

        implicit none

        real,dimension(nyp)       :: umean
        real,dimension(nyp,mz,mx) :: u,v,w,u12,up,vp,wp
#IFDEF POLYMER
        real,dimension(nyp,mz,mx) :: p12
        real :: s_tw,s_bw,tpss,bpss
#ENDIF
        real :: tau_tw,tau_bw,twss,bwss

        ! Reynolds stresses (symmetric tensor)
        real,dimension(nyp) :: R11, R12, R13, R22, R23, R33 

        ! Common block variables
        integer :: irstrt
        integer :: it 
        real    :: xl,yl,zl
        real    :: re

        integer :: i,j,k

        ! Common blocks
        common/iocontrl/   irstrt
        common/itime/      it
        common/domain/     xl,yl,zl
        common/flow/       re
     

        ! Begin Calculations
        ! Calculate fluctuating velocities
        !$omp parallel do
        do k = 1,mx
            do j = 1,mz 
                do i = 1,nyp
                    up(i,j,k) = u(i,j,k) - umean(i)
                    vp(i,j,k) = v(i,j,k)
                    wp(i,j,k) = w(i,j,k)
                end do
            end do
        end do
        !$omp end parallel do

        ! Average values to calculate stress terms
        R11 = 0.0; R12 = 0.0; R13 = 0.0
        R22 = 0.0; R23 = 0.0; R33 = 0.0
        !$omp parallel do default(shared) reduction(+:R11,R12,R13,R22,R23,R33) private(i,j,k) collapse(3)
        do k = 1,mx
            do j = 1,mz 
                do i = 1,nyp
                    R11(i) = R11(i) + up(i,j,k)*up(i,j,k)
                    R12(i) = R12(i) + up(i,j,k)*vp(i,j,k)
                    R13(i) = R13(i) + up(i,j,k)*wp(i,j,k)
                    R22(i) = R22(i) + vp(i,j,k)*vp(i,j,k)
                    R23(i) = R23(i) + vp(i,j,k)*wp(i,j,k)
                    R33(i) = R33(i) + wp(i,j,k)*wp(i,j,k)
                end do
            end do
        end do
        !$omp end parallel do

        R11 = R11/(mx*mz)
        R12 = R12/(mx*mz)
        R13 = R13/(mx*mz)
        R22 = R22/(mx*mz)
        R23 = R23/(mx*mz)
        R33 = R33/(mx*mz)

        tau_tw = 0.0
        tau_bw = 0.0
#IFDEF POLYMER
        s_tw = 0.0
        s_bw = 0.0
        !$omp parallel do reduction(+:tau_tw,tau_bw,s_tw,s_bw) default(shared) private(j,k,twss,bwss,tpss,bpss)
#ELSE
        !$omp parallel do reduction(+:tau_tw,tau_bw) default(shared) private(j,k,twss,bwss)
#ENDIF
        do k = 1,mx
            do j = 1,mz
                twss = u12(1,j,k)
                bwss = u12(nyp,j,k)

                tau_tw = tau_tw + twss/re
                tau_bw = tau_bw + bwss/re
#IFDEF POLYMER
                tpss = p12(1,j,k) 
                bpss = p12(nyp,j,k) 

                s_tw = s_tw + tpss
                s_bw = s_bw + bpss
#ENDIF
            end do
        end do
        !$omp end parallel do                

        if (it .eq. irstrt) then
            open(11,file='outputs/Rey_stress1',status='new')
            open(101,file='outputs/twss1')
            open(102,file='outputs/bwss1')
        else
            open(11,file='outputs/Rey_stress1',status='old',position='append')
            open(101,file='outputs/twss1',position='append')
            open(102,file='outputs/bwss1',position='append')
        end if

!        write(11,*) R11,R12,R13,R22,R23,R33
        do i = 1,nyp
            write(11,*) R12(i)
        end do
        write(101,*) tau_tw/(mx*mz)
        write(102,*) tau_bw/(mx*mz)
        close(11)
        close(101)
        close(102)

#IFDEF POLYMER
        if (it .eq. irstrt) then
            open(103,file='outputs/tpss1')
            open(104,file='outputs/bpss1')
        else
            open(103,file='outputs/tpss1',position='append')
            open(104,file='outputs/bpss1',position='append')
        end if

        write(103,*) s_tw/(mx*mz)
        write(104,*) s_bw/(mx*mz)
        close(103)
        close(104)

#ENDIF
    end subroutine calc_stress
end module helpers
