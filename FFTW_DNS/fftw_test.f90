program fftw_test

    use grid_size
    use, intrinsic :: iso_c_binding
    implicit none
    include 'fftw3.f03'

    type(C_PTR) :: plan
    complex(C_DOUBLE_COMPLEX), dimension(nyp,nz,nxh) :: v,v_trans
    real(C_DOUBLE), dimension(nyp,mz,mx) :: vp3d

    print *,'Read in restart file'
    open(113,file='outputs/v_spectral',status='old',form='unformatted')
    read(113) v
    close(113)

!    write(10,*) v
    print *,'Create c2r FFT Plan'
!    plan = fftw_plan_dft_c2r_3d(nxh,nz,nyp,v,v_trans,FFTW_FORWARD,FFTW_ESTIMATE)
    plan = fftw_plan_dft_c2r_3d(nxh,nz,nyp,v,vp3d,FFTW_ESTIMATE)
    print *,'Execute c2r FFT Plan'
!    call fftw_execute_dft(plan,v,v_trans)
    call fftw_execute_dft_c2r(plan,v,vp3d)
    print *,'Success!'
!    write(3,*) v_trans
    write(3,*) vp3d/(nyp*nz*nxh)
    call fftw_destroy_plan(plan)

end program fftw_test
