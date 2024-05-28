program compare
! ==================================================================================== !
! This program compares the DNS FFT routines with the FFTW package. 
! Value verification is done separately.
!
!   Strategy:
!
!       I will read in a sample spectral data file from the DNS code as an input
!       then do repeated transforms on the data to get it to 3D and back again, looped
!       enough that I can get an average computation time for each case.
!
!   Procedure:
!
!       1. Read data from input file
!       2. Initialize DNS FFTs (timed)
!       3. Execute DNS FFTs (timed, separately from initialization)
!       4. Initialize FFTW (timed)
!       5. Execute FFTW FFTs (timed, separately from initialization)
! ==================================================================================== !

    use,intrinsic :: iso_c_binding
    use grid_size
    use mpi
    use omp_lib

    implicit none

    include 'fftw3.f03'

    ! MPI Variables (for timing only right now)
    integer :: ierr

    ! Input data
    complex(C_DOUBLE_COMPLEX),dimension(nyp,nz,nxh) :: us1,us2,us3

    ! DNS Initialization variables
    real :: trigx(2*nx),trigz(2*nz),trigy(2*ny),sine(ny),cosine(ny)
    integer, dimension(19) :: ixfax,iyfax,izfax,ixfax32,izfax32
    real, dimension(16000) :: trigx32
    real, dimension(4000)  :: trigz32 

    ! DNS Transform variables
    integer :: num_ffts = 1

    ! FFTW Initialization variables
    type(C_PTR) :: plan1,plan2,plan3,plan4,plan5
    type(C_PTR) :: plan11,plan12,plan13,plan14,plan15

    ! FFTW Transform variables
    complex(C_DOUBLE_COMPLEX), dimension(nyp,mz,mx) :: uc
    real(C_DOUBLE), dimension(nyp,mz,mx) :: up
    real(C_DOUBLE), dimension(mz,mx,nyp) :: s2
    complex(C_DOUBLE_COMPLEX), dimension(mz,mx) :: uc2d
    real(C_DOUBLE),dimension(mz,mx) :: up2d

    complex(C_DOUBLE_COMPLEX), dimension(mz)  :: cztemp
    complex(C_DOUBLE_COMPLEX), dimension(mx)  :: cxtemp

    real(C_DOUBLE), dimension(nyp) :: rytemp
    real(C_DOUBLE), dimension(mx)  :: rxtemp

    real,dimension(nyp) :: cheb = 1.0

    ! Timing variables
    double precision :: tstart, tend

    ! Loop indices
    integer :: i,j,k,jj

    ! Common block for DNS FFTs
    common/trig/ trigx,trigz,trigy,sine,cosine,trigx32,trigz32
    common/ifax/ ixfax,iyfax,izfax,ixfax32,izfax32

    ! ------------------------------------------------------------------- !
    !                         Begin Calculations                          !
    ! ------------------------------------------------------------------- !
    call MPI_Init(ierr)

    cheb(1) = 2.0
    cheb(nyp) = 2.0

    ! Read input data
    print *,'Read in spectral data'
    open(111,file='outputs/spectral_data',status='old',form='unformatted')
    read(111) us1
    close(111)
   
    us2 = us1
    us3 = us1 

    ! Initialize DNS FFTs
    print *,'   Initializing DNS FFTs... '
    tstart = MPI_Wtime()
    call rcsexp(nx,ixfax,trigx)
    call ccosexp(ny,sine,cosine,iyfax,trigy)
    call ccexp(nz,izfax,trigz)
    call rcsexp(nx32,ixfax32,trigx32)
    call ccexp(nz32,izfax32,trigz32)
    tend = MPI_Wtime()

    print *,''
    print *,'   DNS FFT Initialization: ',tend-tstart,' s'
    print *,'' 

    ! Execute DNS FFTs
    print *,'   Executing DNS FFTs... '
    tstart = MPI_Wtime()
    call DNS_FFTs(us1,num_ffts)
    tend = MPI_Wtime()
    print *,''
    print *,'   DNS FFT total execution time: ',(tend-tstart),' s'
    print *,''
    print *,'   DNS FFT average execution time: ',(tend-tstart)/float(num_ffts),' s'
    print *,''

    write(1,*) us1(:,1,1)

    ! Initialize FFTW FFTs
    print *,''
    print *,' ------------------------------------------------------------- '
    print *,''
    print *,'   Creating FFTW Plans (1D x 3, array slice)'
    print *,''

    tstart = MPI_Wtime()
    ierr = fftw_init_threads() ! Initialize threaded FFTW
    call fftw_plan_with_nthreads(OMP_GET_MAX_THREADS())

    plan1 = fftw_plan_dft_1d(mz,uc,uc,FFTW_BACKWARD,FFTW_MEASURE)
    plan2 = fftw_plan_dft_c2r_1d(mx,uc,up,FFTW_MEASURE)
    plan3 = fftw_plan_r2r_1d(nyp,up,up,FFTW_REDFT00,FFTW_MEASURE)
    plan4 = fftw_plan_dft_r2c_1d(mx,up,uc,FFTW_MEASURE)
    plan5 = fftw_plan_dft_1d(mz,uc,uc,FFTW_FORWARD,FFTW_MEASURE)
    tend = MPI_Wtime()

    print *,'   FFTW Initialization time: ',tend-tstart,' s'
    print *,''

!    write(1,*) us
    ! Execute FFTW FFTs
    print *,'   Executing FFTW FFTs... '
    print *,''
    tstart = MPI_Wtime()
    ! Copy spectral data to bigger array for FFT
    !$omp parallel do shared(uc)
    do k = 1,mx
    do j = 1,mz
    do i = 1,nyp
    uc(i,j,k) = 0.0
    end do
    end do
    end do
    !$omp end parallel do

    !$omp parallel do shared(uc,us2)
    do k = 1,nxh
        do j = 1,nz
            if (j .le. nzh) jj = j
            if (j .gt. nzh) jj = (mz-nz) + j
            do i = 1,nyp
                uc(i,jj,k) = us2(i,j,k)*cheb(i)/2.0
            end do
        end do
    end do
    !$omp end parallel do

!    call FFTW_FFTs(us2,num_ffts,plan1,plan2,plan3,plan4,plan5)
    call perform_FFTs(uc,plan1,plan2,plan3,plan4,plan5,num_ffts)
    ! Fill out spectral variable
    !$omp parallel do shared(us2,uc)
    do k = 1,nxh
        do j = 1,nz
            if (j .le. nzh) jj = j
            if (j .gt. nzh) jj = (mz-nz) + j
            do i = 1,nyp
                us2(i,j,k) = uc(i,jj,k)
            end do
        end do
    end do 
    !$omp end parallel do
    tend = MPI_Wtime()
!    call fftw_destroy_plan(plan1)
!    call fftw_destroy_plan(plan2)
!    call fftw_destroy_plan(plan3)
!    call fftw_destroy_plan(plan4)
!    call fftw_destroy_plan(plan5)

    print *,''
    print *,'   FFTW FFT total execution time: ',(tend-tstart),' s'
    print *,''
    print *,'   FFTW FFT average execution time: ',(tend-tstart)/float(num_ffts),' s'
    print *,''
    write(2,*) us2(:,1,1)

    ! DEBUGGING
!    stop

    ! Initialize FFTW FFts
    print *,''
    print *,' ------------------------------------------------------------- '
    print *,''
    print *,'   Creating FFTW Plans (1D x 3, array --> vector copy)'
    print *,''

    tstart = MPI_Wtime()
!    ierr = fftw_init_threads() ! Initialize threaded FFTW
!    call fftw_plan_with_nthreads(OMP_GET_MAX_THREADS())
    plan11 = fftw_plan_dft_1d(mz,cztemp,cztemp,FFTW_BACKWARD,FFTW_MEASURE)
    plan12 = fftw_plan_dft_c2r_1d(mx,cxtemp,rxtemp,FFTW_MEASURE)
    plan13 = fftw_plan_r2r_1d(nyp,rytemp,rytemp,FFTW_REDFT00,FFTW_MEASURE)
    plan14 = fftw_plan_dft_r2c_1d(mx,rxtemp,cxtemp,FFTW_MEASURE)
    plan15 = fftw_plan_dft_1d(mz,cztemp,cztemp,FFTW_FORWARD,FFTW_MEASURE)
    tend = MPI_Wtime()

    print *,'   FFTW Initialization time: ',tend-tstart,' s'
    print *,''

    ! Execute FFTW FFTs
    print *,'   Executing FFTW FFTs... '
    print *,''
    tstart = MPI_Wtime()
    ! Copy spectral data to bigger array for FFT
    !$omp parallel do shared(uc)
    do k = 1,mx
    do j = 1,mz
    do i = 1,nyp
    uc(i,j,k) = 0.0
    end do
    end do
    end do
    !$omp end parallel do

    call perform_FFTs_2(uc,plan11,plan12,plan13,plan14,plan15,num_ffts)
    ! Fill out spectral variable
    !$omp parallel do shared(uc,us3)
    do k = 1,nxh
        do j = 1,nz
            if (j .le. nzh) jj = j
            if (j .gt. nzh) jj = (mz-nz) + j
            do i = 1,nyp
                us3(i,j,k) = uc(i,jj,k)
            end do
        end do
    end do 
    !$omp end parallel do
    tend = MPI_Wtime()

    print *,''
    print *,'   FFTW FFT total execution time: ',(tend-tstart),' s'
    print *,''
    print *,'   FFTW FFT average execution time: ',(tend-tstart)/float(num_ffts),' s'
    print *,''


    ! Clean up

    call MPI_Finalize(ierr)

end program compare

! ------------------------------------------------------------------------------------ !

subroutine DNS_FFTs(us,num_ffts)

    use grid_size
    use omp_lib

    implicit none

    ! Input data (spectral field)
    complex,dimension(nyp,nz,nxh) :: us
    
    ! Number of FFTs to do
    integer :: num_ffts

    ! Working arrays
    complex,dimension(nyp,nz,nxh) :: w1,w2,w3,w4

    ! FFT variables
    real :: trigx(2*nx),trigz(2*nz),trigy(2*ny),sine(ny),cosine(ny)
    integer, dimension(19) :: ixfax,iyfax,izfax,ixfax32,izfax32
    real, dimension(16000) :: trigx32
    real, dimension(4000)  :: trigz32 
    integer :: is,inc,jump,lot

    ! Transform variables
    real,dimension(mzp,mxp2)  :: up
    real,dimension(nyp,mz,mx) :: up3d
    real,dimension(nwrk)      :: wrk

    ! Loop indices
    integer :: i,j,k,n,jj,k1,k2

    ! Common block for DNS FFTs
    common/trig/ trigx,trigz,trigy,sine,cosine,trigx32,trigz32
    common/ifax/ ixfax,iyfax,izfax,ixfax32,izfax32

    ! ------------------------------------------------------------------- !
    !                         Begin Calculations                          !
    ! ------------------------------------------------------------------- !

    do n = 1,num_ffts
        ! Forward transform (spectral --> physical)
        is = 1

        ! y-transform (Chebyshev coefficients)
        call yfft(us,w1,w2,w3,w4,is)

        ! Loop over y-planes
        !$omp parallel do 
        do i = 1,nyp

            up = 0.0
       
            do j = 1,nz
                if (j.le. nzh) jj = j
                if (j .gt. nzh) jj = (mz-nz) + j

                do k = 1,nxh
                k1 = 2*(k-1) + 1
                k2 = 2*k

                up(jj,k1) = real(us(i,j,k))
                up(jj,k2) = aimag(us(i,j,k))
                end do ! k
            end do ! j

            ! c2c FFT - 1D transform in x-direction
            is = 1
            inc = 1
            jump = 2*mzp
            lot = nxh
 
            call cfftmlt(up(1,1),up(1,2),wrk,trigz32,izfax32,inc,jump,mz,lot,is)
        
            do j = 1,mz
                up(j,2) = 0.0
            end do


            ! c2r FFT - 1D transform in z-direction

            inc = mzp; jump = 1; lot = mz
            call rfftmlt(up,wrk,trigx32,ixfax32,inc,jump,mx,lot,is)

            ! Done! Now you would normally do stuff with the physical variable here, but we're just going to transform it back to
            ! spectral space

            ! Store in 3D variable for FFTW equivalence (and use)
            up3d(i,:,:) = up(1:mz,1:mx)

            ! Physical --> Spectral Transform
            is = -1
            inc = mzp; jump = 1; lot = mz

            ! r2c 1D transform
            call rfftmlt(up,wrk,trigx32,ixfax32,inc,jump,mx,lot,is)

            do j = 1,mz
                up(j,2) = 0.0
            end do

            ! c2c 1D transform
            inc = 1; jump = 2*mzp; lot = nxh
            call cfftmlt(up(1,1),up(1,2),wrk,trigz32,izfax32,inc,jump,mz,lot,is)

            ! Fill out complex spectral variable
            do j = 1,nz
                if (j .le. nzh) jj = j
                if (j .gt. nzh) jj = (mz-nz) + j
                do k = 1,nxh
                k1 = 2*(k-1) + 1
                k2 = 2*k
                us(i,j,k) = cmplx(up(jj,k1)*rmz,up(jj,k2)*rmz)
                end do ! k
            end do !j 
        end do ! i   
        !$omp end parallel do

        ! Do final inverse transform back to y-spectral
        is = -1
        call yfft(us,w1,w2,w3,w4,is) 
    end do ! n


end subroutine DNS_FFts

! ------------------------------------------------------------------------------------ !

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

! ------------------------------------------------------------------------------------ !

subroutine FFTW_FFTs(us,num_ffts,plan1,plan2,plan3,plan4,plan5)

    use,intrinsic :: iso_c_binding
    use grid_size
    use omp_lib

    implicit none

    include 'fftw3.f03'


    ! Input data
    complex(C_DOUBLE_COMPLEX),dimension(nyp,nz,nxh) :: us
    integer :: num_ffts

    ! FFTW Initialization variables
    type(C_PTR) :: plan1,plan2,plan3,plan4,plan5

    ! FFTW Transform variables
    complex(C_DOUBLE_COMPLEX), dimension(nyp,mz,mx) :: uc
    real(C_DOUBLE), dimension(nyp,mz,mx) :: up

    ! Normalization factor
    real,dimension(mz,mx) :: fac

    ! Loop indices
    integer :: i,j,k,jj,n

    ! ------------------------------------------------------------------- !
    !                         Begin Calculations                          !
    ! ------------------------------------------------------------------- !

    do n = 1,num_ffts

    ! Copy spectral data to bigger array for FFT
    !omp parallel do
    do k = 1,nxh
        do j = 1,nz
            if (j .le. nzh) jj = j
            if (j .gt. nzh) jj = (mz-nz) + j
            do i = 1,nyp
                uc(i,jj,k) = us(i,j,k)
            end do
        end do
    end do
    !omp end parallel do

    ! Complex --> Complex z-transform
    !omp parallel do
    do k = 1,nxh
        do i = 1,nyp
            call fftw_execute_dft(plan1,uc(i,:,k),uc(i,:,k))
        end do
    end do
    !omp end parallel do

    ! Complex --> Real x-transform
    !omp parallel do
    do j = 1,mz
        do i = 1,nyp
            call fftw_execute_dft_c2r(plan2,uc(i,j,:),up(i,j,:))
        end do
    end do
    !omp end parallel do

    ! Real --> Real IDCT
    !omp parallel do
    do k = 1,mx
        do j = 1,mz
            call fftw_execute_r2r(plan3,up(:,j,k),up(:,j,k))
        end do
    end do
    !omp end parallel do

    ! Normalize IDCT for real values
    
!    do k = 1,mx
!        do j = 1,mz
!            fac(j,k) = up(1,j,k)
!            do i = 1,nyp
!                up(i,j,k) = (up(i,j,k) - fac(j,k))/2.0
!            end do
!        end do
!    end do
!    print *,'test u = ',up(nyh,1,1)
!    write(1,*) up
    ! Now do everything in reverse!

    ! Un-normalize real values for DCT
    !omp parallel do
    do k = 1,mx
        do j = 1,mz
            do i = 1,nyp
                up(i,j,k) = 2.0*up(i,j,k) + fac(j,k)
            end do
        end do
    end do
    !omp end parallel do

    ! Real --> Real DCT
    !omp parallel do
    do k = 1,mx
        do j = 1,mz
            call fftw_execute_r2r(plan3,up(:,j,k),up(:,j,k))
        end do
    end do
    !omp end parallel do
!    write(1,*) up

    ! Real --> Complex x-transform
    !omp parallel do
    do j = 1,mz
        do i = 1,nyp
            call fftw_execute_dft_r2c(plan4,up(i,j,:),uc(i,j,:))
        end do
    end do
    !omp end parallel do
!    write(2,*) uc*rmz

    ! Complex --> Complex z-transform
    !omp parallel do
    do k = 1,nxh
        do i = 1,nyp
            call fftw_execute_dft(plan5,uc(i,:,k),uc(i,:,k))
        end do
    end do
    !omp end parallel do
!    write(3,*) uc

    ! Fill out spectral variable
    !omp parallel do
    do k = 1,nxh
        do j = 1,nz
            if (j .le. nzh) jj = j
            if (j .gt. nzh) jj = (mz-nz) + j
            do i = 1,nyp
                us(i,j,k) = uc(i,jj,k)
            end do
        end do
    end do 
    !omp end parallel do
    

    end do ! n


end subroutine FFTW_FFTs

subroutine perform_FFTs(uc,plan1,plan2,plan3,plan4,plan5,num_ffts)

    use,intrinsic :: iso_c_binding
    use grid_size
    use omp_lib

    implicit none
    
    include 'fftw3.f03'

    complex(C_DOUBLE_COMPLEX), dimension(nyp,mz,mx) :: uc,vc,wc

    real(C_DOUBLE),dimension(nyp,mz,mx) :: up,vp,wp

    real,dimension(mz,mx) :: fac1,fac2,fac3
    integer :: i,j,k,n,num_ffts

    type(C_PTR) :: plan1,plan2,plan3,plan4,plan5


    do n = 1,num_ffts

    ! Complex --> Complex z-transform
!    print *,'Complex --> Complex Transform'
    !$omp parallel do shared(uc,plan1)
    do k = 1,nxh
        do i = 1,nyp
            call fftw_execute_dft(plan1,uc(i,:,k),uc(i,:,k))
        end do
    end do
    !$omp end parallel do
!    write(200,*) uc

    ! Complex --> Real x-transform
!    print *,'Complex --> Real Transform'
    !$omp parallel do shared(plan2,uc,up)
    do j = 1,mz
    do i = 1,nyp
        call fftw_execute_dft_c2r(plan2,uc(i,j,:),up(i,j,:))
    end do
    end do
    !$omp end parallel do


    ! Real --> Real IDCT
!    print *,'Real --> Real Cosine Transform'
    !$omp parallel do shared(plan3,up)
    do k = 1,mx
        do j = 1,mz
            call fftw_execute_r2r(plan3,up(:,j,k),up(:,j,k))
        end do
    end do
    !$omp end parallel do

    ! Real --> Real DCT
    !$omp parallel do shared(plan3,up)
    do k = 1,mx
        do j = 1,mz
            call fftw_execute_r2r(plan3,up(:,j,k),up(:,j,k))
        end do
    end do
    !$omp end parallel do
    up = up/float(ny)
    up(1,:,:) = up(1,:,:)/2.0
    up(nyp,:,:) = up(nyp,:,:)/2.0
!    write(1,*) up

    ! Real --> Complex x-transform
    !$omp parallel do shared(plan4,up,uc)
    do j = 1,mz
        do i = 1,nyp
            call fftw_execute_dft_r2c(plan4,up(i,j,:),uc(i,j,:))
        end do
    end do
    !$omp end parallel do
    uc = uc/float(mx)
!    write(2,*) uc*rmz

    ! Complex --> Complex z-transform
    !$omp parallel do shared(plan5,uc)
    do k = 1,nxh
        do i = 1,nyp
            call fftw_execute_dft(plan5,uc(i,:,k),uc(i,:,k))
        end do
    end do
    !$omp end parallel do
    uc = uc/float(mz)
!    write(3,*) uc

    end do ! n


end subroutine 

subroutine perform_FFTs_2(uc,plan1,plan2,plan3,plan4,plan5,num_ffts)

    use,intrinsic :: iso_c_binding
    use grid_size
    use omp_lib

    implicit none
    
    include 'fftw3.f03'

    complex(C_DOUBLE_COMPLEX), dimension(nyp,mz,mx) :: uc
    real(C_DOUBLE),dimension(nyp,mz,mx) :: up

    complex(C_DOUBLE_COMPLEX), dimension(mz)  :: cztemp
    complex(C_DOUBLE_COMPLEX), dimension(mx)  :: cxtemp

    real(C_DOUBLE), dimension(nyp) :: rytemp
    real(C_DOUBLE), dimension(mx)  :: rxtemp

    integer :: i,j,k,n,num_ffts

    type(C_PTR) :: plan1,plan2,plan3,plan4,plan5


    do n = 1,num_ffts

    ! Complex --> Complex z-transform
    !$omp parallel do shared(uc,plan1)
    do k = 1,nxh
        do i = 1,nyp
            cztemp = uc(i,:,k)
            call fftw_execute_dft(plan1,cztemp,cztemp)
            uc(i,:,k) = cztemp
        end do
    end do
    !$omp end parallel do

    ! Complex --> Real x-transform
    !$omp parallel do shared(plan2,uc,up)
    do j = 1,mz
    do i = 1,nyp
        cxtemp = uc(i,j,:)
        call fftw_execute_dft_c2r(plan2,cxtemp,rxtemp)
        up(i,j,:) = rxtemp
    end do
    end do
    !$omp end parallel do


    ! Real --> Real IDCT
    !$omp parallel do shared(plan3,up)
    do k = 1,mx
        do j = 1,mz
            rytemp = up(:,j,k)
            call fftw_execute_r2r(plan3,rytemp,rytemp)
            up(:,j,k) = rytemp
        end do
    end do
    !$omp end parallel do

    ! Real --> Real DCT
    !$omp parallel do shared(plan3,up)
    do k = 1,mx
        do j = 1,mz
            rytemp = up(:,j,k)
            call fftw_execute_r2r(plan3,rytemp,rytemp)
            up(:,j,k) = rytemp
        end do
    end do
    !$omp end parallel do
    up = up/float(ny)
    up(1,:,:) = up(1,:,:)/2.0
    up(nyp,:,:) = up(nyp,:,:)/2.0

    ! Real --> Complex x-transform
    !$omp parallel do shared(plan4,up,uc)
    do j = 1,mz
        do i = 1,nyp
            rxtemp = up(i,j,:)
            call fftw_execute_dft_r2c(plan4,rxtemp,cxtemp)
            uc(i,j,:) = cxtemp
        end do
    end do
    !$omp end parallel do
    uc = uc/float(mx)

    ! Complex --> Complex z-transform
    !$omp parallel do shared(plan5,uc)
    do k = 1,nxh
        do i = 1,nyp
            cztemp = uc(i,:,k)
            call fftw_execute_dft(plan5,cztemp,cztemp)
            uc(i,:,k) = cztemp
        end do
    end do
    !$omp end parallel do
    uc = uc/float(mz)

    end do ! n


end subroutine 
