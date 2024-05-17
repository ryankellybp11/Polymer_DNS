program fftw_test

    use,intrinsic :: iso_c_binding
    use grid_size
    implicit none

    include 'fftw3.f03'

!    integer,parameter :: ny = 128, nyp = 129

    real,dimension(nyp) :: y,initu
!    real,dimension(nyp,nz,nx) :: initu
    integer :: i
    real :: pi = 3.14159265358979

    real(C_DOUBLE),dimension(nyp) :: bplan
    type(C_PTR) :: plan

! ----------------------------------------------------- !

    ! Make fftw plan
    plan = fftw_plan_r2r_1d(nyp,bplan,bplan,FFTW_REDFT00,FFTW_ESTIMATE)

    ! Calculate y vector (Chebyshev grid)
    y(1) = 0.0
    do i = 2,nyp
        y(i) = 1-cos(float(i-1)*pi/float(ny)) ! [-1,+1]
    end do
    
    ! Calculate initu
    do i = 1,nyp
        initu(i) = 10.0*abs(y(i)-y(nyp))/abs(y(1)-y(nyp))
    end do

    ! Write data to output file
    open(1)
    write(1,*) initu
    close(1)

    ! Calculate DCT-I of initu
    call xyzfft(initu,plan)

    print *,''
    print *,'Calculated initu transform: '
    print *,initu

!    write(2,*) initu

    ! Read data back into initu variable
    open(1)
    read(1,*) initu
    close(1)

    call xyzfft(initu,plan)

!    write(3,*) initu
    print *,''
    print *,'Read initu transform: '
    print *,initu

    print *,''
    
    call fftw_destroy_plan(plan)
end program

subroutine xyzfft(b,plan)

    use,intrinsic :: iso_c_binding
    use grid_size
    implicit none

    include 'fftw3.f03'

!    integer, parameter :: ny = 16, nyp = 17

    integer :: i

    real(C_DOUBLE),dimension(nyp) :: b,bplan,cheb

    type(C_PTR) :: plan

! ------------------------------------------------------------------------------- !

    do i = 1,nyp
        cheb(i) = 1.0
    end do
    cheb(1) = 2.0
    cheb(nyp) = 2.0

    call fftw_execute_r2r(plan,b,b)
    b = b/float(ny)

    ! Divide first and last terms by 2  
    b(1) = b(1)/2.0
    b(nyp) = b(nyp)/2.0

end subroutine xyzfft
